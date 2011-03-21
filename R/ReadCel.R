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
        isMale <- rep(TRUE, length())
    }
    else
    {
        isMale <- genderClust == maleClust
    }
    
    isMale
}

.ccstrans <- function(a, b, k = 4) {
    x <- asinh(k * (a - b)/(a + b))/asinh(k)
    y <- (log2(a) + log2(b))/2
    list(x = x, y = y)
}

genotypeAnyChrChunk <- function(
        chr,
        intensityConts, intensityAvgs,
        hint = NULL,
        parIndices = NULL,
        transformMethod = c("CCStrans", "MAtrans"),
        isMale = NULL,
        confScoreThreshold = 1e-05)
{
    transformMethod <- match.arg(transformMethod)
    
    if(chr == "X")
    {
        .genotypeXChromosomeChunk(intensityConts, intensityAvgs, hint, parIndices, transformMethod, isMale, confScoreThreshold)
    }
    else if(chr == "Y")
    {
        .genotypeYChromosomeChunk(intensityConts, intensityAvgs, hint, transformMethod, isMale, confScoreThreshold)
    }
    else if(chr == "M")
    {
        .genotypeHomozygousChunk(intensityConts, intensityAvgs, transformMethod, confScoreThreshold)
    }
    else
    {
        # chr is an autosome
        .genotypeChunk(intensityConts, intensityAvgs, hint, transformMethod, confScoreThreshold)
    }
}

# for genotyping probeset "chunks" (where ms and ss rows map to probesets and
# columns map to arrays (CEL files). This function is useful for genotyping
# autosomes and other sections of the genome where it is possible to have 3
# genotyping outcomes. Otherwise use genotypeHomozygousChunk
.genotypeChunk <- function(ms, ss, hint, trans, confScoreThreshold)
{
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
    
    # geno/vinotype each of the probesets in this chunk
    geno <- matrix(0, nrow = numProbesets, ncol = numArrays)
    vino <- matrix(0, nrow = numProbesets, ncol = numArrays)
    conf <- matrix(0, nrow = numProbesets, ncol = numArrays)
    for(probesetIndex in 1 : numProbesets)
    {
        # TODO genotype sometimes returns results with colnames and sometimes not.
        currVals <- genotype(
            ms[probesetIndex, ],
            ss[probesetIndex, ],
            hint[probesetIndex],
            trans)
        geno[probesetIndex, ] <- currVals$geno
        vino[probesetIndex, ] <- currVals$vino
        conf[probesetIndex, ] <- currVals$conf
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
.genotypeHomozygousChunk <- function(ms, ss, trans, confScoreThreshold)
{
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
            trans)
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

.genotypeXChromosomeChunk <- function(ms, ss, hint, parIndices, trans, isMale, confScoreThreshold)
{
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
    
    # we can treat the females like autosomes for the purposes of genotyping
    if(length(femaleColumns) >= 2) # TODO is this OK if it is equal to 1?
    {
        normalFemaleMs <- ms[normalIndices, femaleColumns, drop = FALSE]
        normalFemaleSs <- ss[normalIndices, femaleColumns, drop = FALSE]
        normalHint <- hint[normalIndices]
        
        # we are going to treat the female probesets just like autosomes for
        # the purposes of geno/vinotyping
        normalFemaleResult <- .genotypeChunk(
            normalFemaleMs,
            normalFemaleSs,
            normalHint,
            trans,
            confScoreThreshold)
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
            confScoreThreshold)
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
        
        parResult <- .genotypeChunk(
            parMs,
            parSs,
            parHint,
            trans,
            confScoreThreshold)
        results <- initializeNamesIfMissing(results, names(parResult))
        for(itemName in names(parResult))
        {
            results[[itemName]][parIndices, ] <- parResult[[itemName]]
        }
    }
    
    results
}

.genotypeYChromosomeChunk <- function(ms, ss, hint, trans, isMale, confScoreThreshold)
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
    maleResult <- .genotypeHomozygousChunk(maleMs, maleSs, trans, confScoreThreshold)
    
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

# this function will convert SNP A/B intensity pairs into M/S values which are
# suitable as input to the genotypeAnyChrChunk function
convertToContrastAndAverage <- function(
        aIntensities,
        bIntensities,
        transformMethod = c("CCStrans", "MAtrans"),
        intensitiesAreLog2 = FALSE,
        k = 4)
{
    transformMethod <- match.arg(transformMethod)
    
    if (transformMethod == "CCStrans")
    {
        if(intensitiesAreLog2)
        {
            aIntensities <- 2 ^ aIntensities
            bIntensities <- 2 ^ bIntensities
        }
        
        res <- .ccstrans(aIntensities, bIntensities, k)
        M <- res$x
        S <- res$y
    }
    else if (transformMethod == "MAtrans")
    {
        if(!intensitiesAreLog2)
        {
            aIntensities <- log2(aIntensities)
            bIntensities <- log2(bIntensities)
        }
        
        M <- aIntensities - bIntensities
        S <- (aIntensities + bIntensities) / 2
    }
    else
    {
        stop(paste("bad transformation argument:", transformMethod))
    }
    
    list(intensityConts = M, intensityAvgs = S)
}

# reads in the given CEL file, partitions it into groups defined by groupLevels
# then breaks up those groups into chunk sizes no bigger than maxChunkSize
# groups should be factors where
normalizeCelFileByChr <- function(
    celFileName,
    snpProbeInfo,
    snpInfo,
    verbose = FALSE,
    chromosomes = c(as.character(1 : 19), "X", "Y", "M"),
    referenceDistribution = NULL,
    transformMethod = c("CCStrans", "MAtrans", "None"),
    keepSNPIds = FALSE)
{
    if(verbose) cat("Reading and normalizing CEL file: ", celFileName, "\n", sep="")
    
    transformMethod <- match.arg(transformMethod)
    
    celData <- read.celfile(celFileName, intensity.means.only = TRUE)
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
    
    # TODO we could do something here which is more robust
    aSnpIds <- rownames(allAint)
    bSnpIds <- rownames(allBint)
    if(!all(aSnpIds == bSnpIds))
    {
        stop("The SNP IDs for the A alleles should match up with the SNP IDs ",
             "for the B alleles but they do not.")
    }
    
    if (transformMethod == "CCStrans") {
        # fixed K??
        res <- .ccstrans(2^allAint, 2^allBint)
        M <- res$x
        S <- res$y
    }
    else if (transformMethod == "MAtrans") {
        # then prior??
        M <- allAint - allBint
        S <- (allAint + allBint)/2
    }
    
    # divide CEL file up into chromosome pieces
    msList <- list()
    for (chri in chromosomes) {
        #currRows <- which(chrid == chri)
        currChrSnps <- snpInfo$snpId[snpInfo$chrId == chri]
        currRows <- match(currChrSnps, aSnpIds)
        currRows <- currRows[!is.na(currRows)]
        
        numCurrRows <- length(currRows)
        if(numCurrRows == 0)
        {
            stop("Failed to find any probes on chromosome ", chri)
        }
        
        if (transformMethod == "None")
        {
            msList[[chri]]$a <- allAint[currRows]
            msList[[chri]]$b <- allBint[currRows]
        }
        else
        {
            msList[[chri]]$intensityConts <- M[currRows]
            msList[[chri]]$intensityAvgs <- S[currRows]
        }
        
        if(keepSNPIds)
        {
            names(msList[[chri]][[1]]) <- aSnpIds[currRows]
            names(msList[[chri]][[2]]) <- names(msList[[chri]][[1]])
        }
    }
    
    msList
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
        verbose = FALSE)
{
    readNext <- function(index)
    {
        if(index > length(celFiles))
        {
            NULL
        }
        else
        {
            if(verbose) cat("Reading and normalizing CEL file: ", celFiles[index], "\n", sep="")
            
            celData <- read.celfile(celFiles[index], intensity.means.only = TRUE)
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
            
            sampleName <- .fileBaseWithoutExtension(celFiles[index])
            list(
                head = list(sampleName = sampleName, sampleData = currDataMatrix),
                tail = function() readNext(index + 1))
        }
    }
    
    readNext(1)
}

.readSnpIntensitiesFromTab <- function(tabFile, verbose = FALSE)
{
    mustCloseConnection <- FALSE
    if(!inherits(tabFile, "connection"))
    {
        tabFile <- file(tabFile, "rt")
        mustCloseConnection <- TRUE
    }
    
    # wrap read.table function for convenience
    readSnpIntenTable <- function(file, ...)
    {
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
    currSampleData <- NULL
    repeat
    {
        currBlock <- readSnpIntenTable(tabFile, nrows = 10000)
        if(nrow(currBlock) == 0)
        {
            break
        }
        
        currSampleData <- rbind(currSampleData, currBlock)
        
        if(currSampleData[[1]][1] != currSampleData[[1]][nrow(currSampleData)])
        {
            break
        }
    }
    snpCount <- sum(currSampleData[[1]] == currSampleData[[1]][1])
    
    readNext <- function(leftovers)
    {
        numLeftovers <- if(is.null(leftovers)) 0 else nrow(leftovers)
        numLeftoversToTake <- min(snpCount, numLeftovers)
        if(numLeftoversToTake >= 1)
        {
            leftoversToTake <- 1 : numLeftoversToTake
            taken <- leftovers[leftoversToTake, ]
            leftovers <- leftovers[-leftoversToTake, ]
        }
        else
        {
            taken <- NULL
        }
        
        numNewToTake <- snpCount - numLeftoversToTake
        if(numNewToTake >= 1)
        {
            newData <- readSnpIntenTable(tabFile, nrows = numNewToTake)
            currData <- rbind(taken, newData)
        }
        else
        {
            currData <- taken
        }
        
        if(is.null(currData) || nrow(currData) == 0)
        {
            if(mustCloseConnection)
            {
                close(tabFile)
            }
            
            NULL
        }
        else
        {
            if(verbose) cat("Finished reading data for sample: ", currData[[1]][1], "\n", sep="")
            
            currDataMatrix <- matrix(c(currData[[3]], currData[[4]]), ncol = 2)
            rownames(currDataMatrix) <- currData[[2]]
            colnames(currDataMatrix) <- c("A", "B")
            
            list(
                head = list(sampleName = currData[[1]][1], sampleData = currDataMatrix),
                tail = function() readNext(leftovers))
        }
    }
    
    readNext(currSampleData)
}

.ccsTransformSample <- function(sample)
{
    result <- .ccstrans(2 ^ sample$sampleData[ , 1], 2 ^ sample$sampleData[ , 2])
    sampleData <- matrix(c(result$x, result$y), ncol = 2)
    rownames(sampleData) <- rownames(sample$sampleData)
    colnames(sampleData) <- c("intensityConts", "intensityAvgs")
    
    list(sampleName = sample$sampleName, sampleData = sampleData)
}

.maTransformSample <- function(sample)
{
    m <- sample$sampleData[ , 1] - sample$sampleData[ , 2]
    s <- (sample$sampleData[ , 1] + sample$sampleData[ , 2]) / 2
    sampleData <- matrix(c(m, s), ncol = 2)
    rownames(sampleData) <- rownames(sample$sampleData)
    colnames(sampleData) <- c("intensityConts", "intensityAvgs")
    
    list(sampleName = sample$sampleName, sampleData = sampleData)
}

.lazyApply <- function(f, xs)
{
    if(is.null(xs))
    {
        NULL
    }
    else
    {
        list(head = f(xs$head), tail = function() .lazyApply(f, xs$tail()))
    }
}