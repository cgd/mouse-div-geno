#########################################################################
# ReadCel.R
#
# Part of the MouseDivGenotype package
#
#########################################################################

inferGender <- function(
    meanIntensityXPerArray,
    meanIntensityYPerArray,
    meanIntensityPerAutosome)
{
    genderClust <- kmeans(cbind(meanIntensityYPerArray, meanIntensityXPerArray), 2)$cluster
    yMean1 <- mean(meanIntensityYPerArray[genderClust == 1])
    yMean2 <- mean(meanIntensityYPerArray[genderClust == 2])
    
    # we expect that the cluster with a more intense Y chromosome will be males
    if(yMean1 > yMean2)
    {
        maleClust <- 1
        yMeanMale <- yMean1
        yMeanFemale <- yMean2
    }
    else
    {
        maleClust <- 2
        yMeanMale <- yMean2
        yMeanFemale <- yMean1
    }
    
    if(yMeanFemale > min(meanIntensityPerAutosome))
    {
        # if the y intensity of the samples grouped as female is greater than
        # the smallest autosome mean intensity we infer that there are no
        # female samples
        isMale <- rep(TRUE, length(meanIntensityXPerArray))
    }
    else
    {
        isMale <- genderClust == maleClust
    }
    
    isMale
}

genotypeAnyChrChunk <- function(
        chr,
        intensityConts, intensityAvgs,
        hint = NULL,
        isInPAR = NULL,
        transformMethod = c("CCStrans", "MAtrans"),
        isMale = NULL,
        confScoreThreshold = 1e-05,
        logSnp = NULL,
        logfn = NULL)
{
    transformMethod <- match.arg(transformMethod)
    
    if(chr == "X")
    {
        .genotypeXChromosomeChunk(intensityConts, intensityAvgs, hint, isInPAR, transformMethod, isMale, confScoreThreshold, logSnp, logfn)
    }
    else if(chr == "Y")
    {
        .genotypeYChromosomeChunk(intensityConts, intensityAvgs, hint, transformMethod, isMale, confScoreThreshold, logSnp, logfn)
    }
    else if(chr == "M")
    {
        .genotypeHomozygousChunk(intensityConts, intensityAvgs, transformMethod, confScoreThreshold, logSnp, logfn, "M")
    }
    else
    {
        # chr is an autosome
        .genotypeChunk(intensityConts, intensityAvgs, hint, transformMethod, confScoreThreshold, logSnp, logfn, "autosome")
    }
}

# for genotyping probeset "chunks" (where ms and ss rows map to probesets and
# columns map to arrays (CEL files). This function is useful for genotyping
# autosomes and other sections of the genome where it is possible to have 3
# genotyping outcomes. Otherwise use genotypeHomozygousChunk
.genotypeChunk <- function(ms, ss, hint, trans, confScoreThreshold, logSnp=NULL, logfn=NULL, kind = c("autosome", "PAR", "femaleX"))
{
    kind <- match.arg(kind)
    if(!is.null(logfn)) {
        kindDescrip <- switch(kind, autosome="autosome", PAR="pseudoautosomal", "female X chromosome")
        logfn(
            "performing standard genotyping for %i %s SNPs and %i samples",
            nrow(ms),
            kindDescrip,
            ncol(ms))
    }
    
    numArrays <- ncol(ms)
    numProbesets <- nrow(ms)
    
    if(ncol(ss) != numArrays || nrow(ss) != numProbesets)
    {
        stop("Internal error: ms matrix dimensions should match ss matrix dimensions")
    }
    
    if(length(hint) == 0)
    {
        hint <- rep(0, numProbesets)
    }
    
    if(is.null(logfn)) {
        logSnp <- NULL
    }
    
    snpNames <- rownames(ms)
    
    # geno/vinotype each of the probesets in this chunk
    geno <- matrix(0, nrow = numProbesets, ncol = numArrays)
    vino <- matrix(0, nrow = numProbesets, ncol = numArrays)
    conf <- matrix(0, nrow = numProbesets, ncol = numArrays)
    for(probesetIndex in 1 : numProbesets)
    {
        if(!is.null(logSnp) && logSnp[probesetIndex]) {
            if(is.null(snpNames)) {
                logfn("==== GENOTYPING AND VINOTYPING SNP ====")
            } else {
                logfn("==== GENOTYPING AND VINOTYPING %s ====", snpNames[probesetIndex])
            }
        }
        
        currVals <- genotype(
            ms[probesetIndex, ],
            ss[probesetIndex, ],
            hint[probesetIndex],
            trans,
            if(!is.null(logSnp) && logSnp[probesetIndex]) logfn else NULL)
        geno[probesetIndex, ] <- currVals$geno
        vino[probesetIndex, ] <- currVals$vino
        conf[probesetIndex, ] <- currVals$conf

        if(!is.null(logSnp) && logSnp[probesetIndex]) {
            if(is.null(snpNames)) {
                logfn("==== FINISHED GENOTYPING AND VINOTYPING SNP ====")
            } else {
                logfn("==== FINISHED GENOTYPING AND VINOTYPING %s ====", snpNames[probesetIndex])
            }
        }
    }
    chunkResult <- list(geno = geno, vino = vino, conf = conf)
    
    if(confScoreThreshold > 0)
    {
        chunkResult$geno[chunkResult$conf < confScoreThreshold] <- -1
    }
    
    # TODO restore rownames and colnames to chunk result
    chunkResult
}

# for genotyping probeset "chunks" (where ms and ss rows map to probesets and
# columns map to arrays (CEL files). This function is useful for genotyping
# the parts of the genome where you do not expect to observe heterozygous
# alleles
.genotypeHomozygousChunk <- function(ms, ss, trans, confScoreThreshold, logSnp=NULL, logfn=NULL, kind = c("M", "maleX", "Y"))
{
    kind <- match.arg(kind)
    if(!is.null(logfn)) {
        kindDescrip <- switch(kind, M="mitochondrial", maleX="male X chromosome", Y="Y chromosome")
        logfn(
            "performing homozygous genotyping for %i %s SNPs and %i samples",
            nrow(ms),
            kindDescrip,
            ncol(ms))
    }
    
    numArrays <- ncol(ms)
    numProbesets <- nrow(ms)
    
    if(ncol(ss) != numArrays || nrow(ss) != numProbesets)
    {
        stop("Internal error: ms matrix dimensions should match ss matrix dimensions")
    }
    
    # geno/vinotype each of the probesets in this chunk
    geno <- matrix(0, nrow = numProbesets, ncol = numArrays)
    vino <- matrix(0, nrow = numProbesets, ncol = numArrays)
    conf <- matrix(0, nrow = numProbesets, ncol = numArrays)
    for(probesetIndex in 1 : numProbesets)
    {
        currVals <- genotypeHomozygous(
            ms[probesetIndex, ],
            ss[probesetIndex, ],
            trans,
            if(!is.null(logSnp) && logSnp[probesetIndex]) logfn else NULL)
        geno[probesetIndex, ] <- currVals$geno
        vino[probesetIndex, ] <- currVals$vino
        conf[probesetIndex, ] <- currVals$conf
    }
    chunkResult <- list(geno = geno, vino = vino, conf = conf)
    
    if(confScoreThreshold > 0)
    {
        chunkResult$geno[chunkResult$conf < confScoreThreshold] <- -1
    }
    
    # TODO restore rownames to chunk result
    chunkResult
}

.genotypeXChromosomeChunk <- function(ms, ss, hint, isInPAR, trans, isMale, confScoreThreshold, logSnp=NULL, logfn=NULL)
{
    parIndices <- which(as.logical(isInPAR))
    numArrays <- ncol(ms)
    numProbesets <- nrow(ms)
    
    if(ncol(ss) != numArrays || nrow(ss) != numProbesets)
    {
        stop("Internal error: ms matrix dimensions should match ss matrix dimensions")
    }
    
    if(length(hint) == 0)
    {
        hint <- rep(0, numProbesets)
    }
    
    # we'll call the X probesets that aren't in the PAR "normal"
    normalIndices <- 1 : numProbesets
    if(length(parIndices) >= 1)
    {
        normalIndices <- setdiff(normalIndices, parIndices)
    }
    
    maleColumns <- which(isMale)
    femaleColumns <- which(!isMale)
    
    initializeNamesIfMissing <- function(matrixList, itemNames)
    {
        for(itemName in itemNames)
        {
            if(!(itemName %in% names(matrixList)))
            {
                matrixList[[itemName]] <- matrix(nrow = numProbesets, ncol = numArrays)
            }
        }
        
        matrixList
    }
    results <- list()
    
    normalHint <- hint[normalIndices]
    normalLogSnp <- logSnp[normalIndices]
    
    # we can treat the females like autosomes for the purposes of genotyping
    if(length(femaleColumns) >= 2) # TODO is this OK if it is equal to 1?
    {
        normalFemaleMs <- ms[normalIndices, femaleColumns, drop = FALSE]
        normalFemaleSs <- ss[normalIndices, femaleColumns, drop = FALSE]
        
        # we are going to treat the female probesets just like autosomes for
        # the purposes of geno/vinotyping
        normalFemaleResult <- .genotypeChunk(
            normalFemaleMs,
            normalFemaleSs,
            normalHint,
            trans,
            confScoreThreshold,
            normalLogSnp,
            logfn,
            "femaleX")
        results <- initializeNamesIfMissing(results, names(normalFemaleResult))
        for(itemName in names(normalFemaleResult))
        {
            results[[itemName]][normalIndices, femaleColumns] <-
                normalFemaleResult[[itemName]]
        }
    }
    
    # the non-PAR male probesets will get genotyped as homozygous
    if(length(maleColumns) >= 2) # TODO is this OK if it is equal to 1?
    {
        normalMaleMs <- ms[normalIndices, maleColumns, drop = FALSE]
        normalMaleSs <- ss[normalIndices, maleColumns, drop = FALSE]
        
        normalMaleResult <- .genotypeHomozygousChunk(
            normalMaleMs,
            normalMaleSs,
            trans,
            confScoreThreshold,
            normalLogSnp,
            logfn,
            "maleX")
        results <- initializeNamesIfMissing(results, names(normalMaleResult))
        for(itemName in names(normalMaleResult))
        {
            results[[itemName]][normalIndices, maleColumns] <-
                normalMaleResult[[itemName]]
        }
    }
    
    # the PAR probesets will be treated as autosomes whether they're male or female
    if(length(parIndices) >= 1)
    {
        parMs <- ms[parIndices, , drop = FALSE]
        parSs <- ss[parIndices, , drop = FALSE]
        parHint <- hint[parIndices]
        parLogSnp <- logSnp[parIndices]
        
        parResult <- .genotypeChunk(
            parMs,
            parSs,
            parHint,
            trans,
            confScoreThreshold,
            parLogSnp,
            logfn,
            "PAR")
        results <- initializeNamesIfMissing(results, names(parResult))
        for(itemName in names(parResult))
        {
            results[[itemName]][parIndices, ] <- parResult[[itemName]]
        }
    }
    
    results
}

.genotypeYChromosomeChunk <- function(ms, ss, hint, trans, isMale, confScoreThreshold, logSnp=NULL, logfn=NULL)
{
    numArrays <- ncol(ms)
    numProbesets <- nrow(ms)
    
    if(ncol(ss) != numArrays || nrow(ss) != numProbesets)
    {
        stop("Internal error: ms matrix dimensions should match ss matrix dimensions")
    }
    
    maleColumns <- which(isMale)
    
    maleMs <- ms[, maleColumns, drop = FALSE]
    maleSs <- ss[, maleColumns, drop = FALSE]
    maleResult <- .genotypeHomozygousChunk(maleMs, maleSs, trans, confScoreThreshold, logSnp, logfn, "Y")
    
    # fill in the female columns using NA so that the dimensions of the
    # input data match the dimensions of the output data
    fillInNAFemaleCols <- function(maleMatrix)
    {
        allMatrix <- matrix(nrow = numProbesets, ncol = numArrays)
        allMatrix[, maleColumns] <- maleMatrix
        allMatrix
    }
    lapply(maleResult, fillInNAFemaleCols)
}

normalizeCelFile <- function(
        celFile,
        snpProbeInfo,
        referenceDistribution = NULL,
        logFile = NULL) {

    if(!is.null(logFile) && !inherits(logFile, "connection")) {
        stop("the logFile parameter should be a connection")
    }
    
    if(!is.null(logFile)) {
        cat("Reading and normalizing CEL file: ", celFile, "\n", file=logFile, sep="")
    }
    
    celData <- read.celfile(celFile, intensity.means.only = TRUE)
    y <- log2(as.matrix(celData[["INTENSITY"]][["MEAN"]][snpProbeInfo$probeIndex]))
    if (length(snpProbeInfo$correction) > 0)
        # C+G and fragment length correction y
        y <- y + snpProbeInfo$correction
    if (length(referenceDistribution) > 0)
        y <- normalize.quantiles.use.target(y, target = referenceDistribution)
    
    # separate A alleles from B alleles and summarize to the probeset level
    allAint <- y[snpProbeInfo$isAAllele, 1, drop = FALSE]
    allBint <- y[!snpProbeInfo$isAAllele, 1, drop = FALSE]
    allAint <- subColSummarizeMedian(
        matrix(allAint, ncol = 1),
        snpProbeInfo$snpId[snpProbeInfo$isAAllele])
    allBint <- subColSummarizeMedian(
        matrix(allBint, ncol = 1),
        snpProbeInfo$snpId[!snpProbeInfo$isAAllele])
    
    aSnpIds <- rownames(allAint)
    bSnpIds <- rownames(allBint)
    if(!all(aSnpIds == bSnpIds))
    {
        stop("The SNP IDs for the A alleles should match up with the SNP IDs ",
            "for the B alleles but they do not.")
    }
    
    currDataMatrix <- matrix(c(allAint, allBint), ncol = 2)
    rownames(currDataMatrix) <- aSnpIds
    colnames(currDataMatrix) <- c("A", "B")
    
    currDataMatrix
}

.chunkIndices <- function(to, by) {
    chunks <- list()
    chunkNumber <- 0
    for(chunkStart in seq(from = 1, to = to, by = by))
    {
        chunkEnd <- chunkStart + by - 1
        if(chunkEnd > to)
        {
            chunkEnd <- to
        }
        
        chunkNumber <- chunkNumber + 1
        chunks[[chunkNumber]] <- list(start=chunkStart, end=chunkEnd)
    }
    
    chunks
}

# TODO we need more parameters to make sure this is really meaningfully unique
.chunkFileName <- function(baseDir, kind, sampleName, chrName, chunkSize, chunkIndex)
{
    chunkFile <- paste(sampleName, "-", kind, "-", chrName, "-", chunkSize, "-", chunkIndex, ".RData", sep="")
    
    file.path(baseDir, chunkFile)
}

.fileBaseWithoutExtension <- function(celFileName)
{
    sub("\\.[^\\.]*$", "", basename(celFileName))
}

.readSnpIntensitiesFromCEL <- function(
        celFiles,
        snpProbeInfo,
        referenceDistribution = NULL,
        logFile = NULL) {

    readNext <- function(index) {
        headFunc <- function() {
            currDataMatrix <- normalizeCelFile(celFiles[index], snpProbeInfo, referenceDistribution, logFile)
            sampleName <- .fileBaseWithoutExtension(celFiles[index])
            
            list(sampleName = sampleName, sampleData = currDataMatrix)
        }
        
        tailFunc <-
            if(index == length(celFiles)) {
                NULL
            } else {
                function() readNext(index + 1)
            }
        
        list(head = headFunc, tail = tailFunc)
    }
    
    if(1 > length(celFiles)) {
        NULL
    } else {
        function() {readNext(1)}
    }
}

.readSnpIntensitiesFromTab <- function(tabFile, logFile = NULL) {

    if(!is.null(logFile) && !inherits(logFile, "connection")) {
        stop("the logFile parameter should be a connection")
    }
    
    logOn <- !is.null(logFile)
    logf <- function(fmt, ...) {
        if(logOn) {
            cat(sprintf(fmt, ...), file=logFile)
        }
    }
    logfn <- function(fmt, ...) {
        if(logOn) {
            logf(fmt, ...)
            cat("\n", file=logFile)
        }
    }
    
    mustCloseConnection <- FALSE
    if(!inherits(tabFile, "connection")) {
        tabFile <- file(tabFile, "rt")
        mustCloseConnection <- TRUE
    }
    
    # wrap read.table function for convenience
    readSnpIntenTable <- function(file, ...) {
        read.table(
                file,
                header = FALSE,
                sep = "\t",
                row.names = NULL,
                stringsAsFactors = FALSE,
                col.names = c("sampleID", "snpId", "A", "B"),
                colClasses = c("character", "character", "numeric", "numeric"),
                ...)
    }
    
    # for our initial read we need to determine a SNP count
    readFirst <- function() {
        currSampleData <- NULL
        repeat {
            currBlock <- readSnpIntenTable(tabFile, nrows = 10000)
            if(nrow(currBlock) == 0) {
                break
            }
            
            currSampleData <- rbind(currSampleData, currBlock)
            
            if(currSampleData[[1]][1] != currSampleData[[1]][nrow(currSampleData)]) {
                break
            }
        }
        snpCount <- sum(currSampleData[[1]] == currSampleData[[1]][1])
        
        readNext <- function(leftovers) {
            numLeftovers <- nrow(leftovers)
            if(numLeftovers < snpCount + 1) {
                numNewToTake <- 1 + snpCount - numLeftovers
                newData <- readSnpIntenTable(tabFile, nrows = numNewToTake)
                leftovers <- rbind(leftovers, newData)
                numLeftovers <- nrow(leftovers)
            }
            
            if(numLeftovers < snpCount) {
                stop("bad row count")
            }
            
            leftoversToTake <- 1 : snpCount
            currData <- leftovers[leftoversToTake, , drop = FALSE]
            leftovers <- leftovers[-leftoversToTake, , drop = FALSE]
            numLeftovers <- nrow(leftovers)
            
            if(numLeftovers == 0) {
                if(mustCloseConnection) {
                    close(tabFile)
                }
                
                tailFunc <- NULL
            } else {
                tailFunc <- function() {readNext(leftovers)}
            }
            
            logfn("Finished reading data for sample: %s", currData[[1]][1])
            
            currDataMatrix <- matrix(c(currData[[3]], currData[[4]]), ncol = 2)
            rownames(currDataMatrix) <- currData[[2]]
            colnames(currDataMatrix) <- c("A", "B")
            
            headFunc <- function() {
                list(sampleName = currData[[1]][1], sampleData = currDataMatrix)
            }
            
            list(
                head = headFunc,
                tail = tailFunc)
        }
        
        readNext(currSampleData)
    }
    
    readFirst
}

ccsTransform <- function(abMat, k = 4) {
    a <- abMat[ , 1]
    twoPowA <- 2 ^ a
    b <- abMat[ , 2]
    twoPowB <- 2 ^ b
    
    m <- asinh(k * (twoPowA - twoPowB) / (twoPowA + twoPowB)) / asinh(k)
    s <- (a + b) / 2
    
    contAvgMat <- matrix(c(m, s), ncol = 2)
    rownames(contAvgMat) <- rownames(abMat)
    colnames(contAvgMat) <- c("intensityConts", "intensityAvgs")
    
    contAvgMat
}

maTransform <- function(abMat) {
    a <- abMat[ , 1]
    b <- abMat[ , 2]
    
    m <- a - b
    s <- (a + b) / 2
    
    contAvgMat <- matrix(c(m, s), ncol = 2)
    rownames(contAvgMat) <- rownames(abMat)
    colnames(contAvgMat) <- c("intensityConts", "intensityAvgs")
    
    contAvgMat
}

.lazyApply <- function(f, lazyVals) {
    if(is.null(lazyVals)) {
        NULL
    } else {
        newVals <- function() {
            val <- lazyVals()
            list(
                head = function() {f(val$head())},
                tail = .lazyApply(f, val$tail))
        }
        
        newVals
    }
}
