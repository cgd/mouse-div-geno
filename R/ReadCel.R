#########################################################################
# ReadCel.R
#
# Part of the MouseDivGenotype package
#
#########################################################################
library("cluster")

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

ccstrans = function(a, b, k = 4) {
    x = asinh(k * (a - b)/(a + b))/asinh(k)
    y = (log2(a) + log2(b))/2
    list(x = x, y = y)
}

genotypeAnyChrChunk <- function(chr, ms, ss, hint, trans, doCNV, isMale)
{
    if(chr == "X")
    {
        genotypeXChromosomeChunk(ms, ss, hint, trans, doCNV, isMale)
    }
    else if(chr == "Y")
    {
        genotypeYChromosomeChunk(ms, ss, hint, trans, doCNV, isMale)
    }
    else if(chr == "M")
    {
        genotypeHomozygousChunk(ms, ss, trans, doCNV)
    }
    else
    {
        # chr is an autosome
        genotypeChunk(ms, ss, hint, trans, doCNV)
    }
}

# for genotyping probeset "chunks" (where ms and ss rows map to probesets and
# columns map to arrays (CEL files). This function is useful for genotyping
# autosomes and other sections of the genome where it is possible to have 3
# genotyping outcomes. Otherwise use genotypeHomozygousChunk
genotypeChunk <- function(ms, ss, hint, trans, doCNV)
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
    
    # geno/vinotype each of the probesets in this chunk and
    # append the results to our list of matrices
    chunkResult <- list(geno = NULL, vino = NULL, conf = NULL, baf = NULL, llr = NULL)
    for(probesetIndex in 1 : numProbesets)
    {
        # TODO genotype sometimes returns results with colnames and sometimes
        #      not. try to understand why this is happening
        currVals <- genotype(
            ms[probesetIndex, ],
            ss[probesetIndex, ],
            hint[probesetIndex],
            trans,
            doCNV)
        chunkResult <- mapply(rbind, chunkResult, currVals, SIMPLIFY = FALSE)
    }
    
    # TODO restore rownames and colnames to chunk result
    chunkResult
}

# for genotyping probeset "chunks" (where ms and ss rows map to probesets and
# columns map to arrays (CEL files). This function is useful for genotyping
# the parts of the genome where you do not expect to observe heterozygous
# alleles
genotypeHomozygousChunk <- function(ms, ss, trans, doCNV)
{
    numArrays <- ncol(ms)
    numProbesets <- nrow(ms)
    
    if(ncol(ss) != numArrays || nrow(ss) != numProbesets)
    {
        stop("Internal error: ms matrix dimensions should match ss matrix dimensions")
    }
    
    # geno/vinotype each of the probesets in this chunk and
    # append the results to our list of matrices
    chunkResult <- list(geno = NULL, vino = NULL, conf = NULL, baf = NULL, llr = NULL)
    for(probesetIndex in 1 : numProbesets)
    {
        currVals <- genotypeHomozygous(
            ms[probesetIndex, ],
            ss[probesetIndex, ],
            trans,
            doCNV)
        chunkResult <- mapply(rbind, chunkResult, currVals, SIMPLIFY = FALSE)
    }
    
    # TODO restore rownames to chunk result
    chunkResult
}

genotypeXChromosomeChunk <- function(ms, ss, hint, trans, doCNV, isMale)
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
    
    # TODO need to make pseudoautosomal IDs a parameter somehow. Also think
    #      whether or not rownames is a robust way to do the matching
    
    # the probesets in the pseudo autosomal region should be genotyped as if
    # they come from an autosome
    pseudoautosomalID <- c(
        "JAX00723347", "JAX00723349", "JAX00723350", "JAX00723351", 
        "JAX00187253", "JAX00723353", "JAX00723354", "JAX00723355", "JAX00723356", 
        "JAX00723358", "JAX00723359", "JAX00723360", "JAX00723361", "JAX00723364", 
        "JAX00723365", "JAX00723366", "JAX00723367", "JAX00723368", "JAX00723369", 
        "JAX00723370", "JAX00187255", "JAX00187256", "JAX00725011", "JAX00725012", 
        "JAX00723371", "JAX00723372")
    # TODO it seems that rownames(ms) is NULL. fix that
    parIndices <- which(pseudoautosomalID %in% rownames(ms))
    #parIndices <- parIndices[!is.na(parIndices)]
    
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
    if(length(femaleColumns) >= 2) # TODO what if it is equal to 1? still OK i guess?
    {
        normalFemaleMs <- ms[normalIndices, femaleColumns, drop = FALSE]
        normalFemaleSs <- ss[normalIndices, femaleColumns, drop = FALSE]
        normalHint <- hint[normalIndices]
        
        # we are going to treat the female probesets just like autosomes for
        # the purposes of geno/vinotyping
        normalFemaleResult <- genotypeChunk(
            normalFemaleMs,
            normalFemaleSs,
            normalHint,
            trans,
            doCNV)
        results <- initializeNamesIfMissing(results, names(normalFemaleResult))
        for(itemName in names(normalFemaleResult))
        {
            results[[itemName]][normalIndices, femaleColumns] <-
                normalFemaleResult[[itemName]]
        }
    }
    
    # the non-PAR male probesets will get genotyped as homozygous
    if(length(maleColumns) >= 2) # TODO what if it is equal to 1? what should we do then??
    {
        normalMaleMs <- ms[normalIndices, maleColumns, drop = FALSE]
        normalMaleSs <- ss[normalIndices, maleColumns, drop = FALSE]
        
        normalMaleResult <- genotypeHomozygousChunk(
            normalMaleMs,
            normalMaleSs,
            trans,
            doCNV)
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
        
        parResult <- genotypeChunk(
            parMs,
            parSs,
            parHint,
            trans,
            doCNV)
        results <- initializeNamesIfMissing(results, names(parResult))
        for(itemName in names(parResult))
        {
            results[parIndices, ] <- parResult
        }
    }
    
    results
}

genotypeYChromosomeChunk <- function(ms, ss, hint, trans, doCNV, isMale)
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
    maleResult <- genotypeHomozygousChunk(maleMs, maleSs, trans, doCNV)
    
    # TODO using NA for invalid values here. make sure that's consistent with
    #      the rest of the code
    
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

# reads in the given CEL file, partitions it into groups defined by groupLevels
# then breaks up those groups into chunk sizes no bigger than maxChunkSize
# groups should be factors where
normalizeCelFileByChr <- function(
    celFileName,
    probesetChunkSize,
    verbose,
    allid,
    Aid,
    Bid,
    CGFLcorrection,
    reference,
    SNPname,
    trans,
    mchr1,
    chrid)
{
    if(verbose) cat("Reading and normalizing CEL file: ", celFileName, "\n", sep="")
    
    celData <- read.celfile(celFileName, intensity.means.only = TRUE)
    y <- log2(as.matrix(celData[["INTENSITY"]][["MEAN"]][allid]))
    if (length(CGFLcorrection) > 0)
        # C+G and fragment length correction y
        y = y + CGFLcorrection
    if (length(reference) > 0) 
        y <- normalize.quantiles.use.target(y, target = reference)
    
    # separate A alleles from B alleles
    allAint = y[Aid, 1, drop = FALSE]
    allBint = y[Bid, 1, drop = FALSE]
    allAint <- subColSummarizeMedian(matrix(allAint, ncol = 1), SNPname)
    allBint <- subColSummarizeMedian(matrix(allBint, ncol = 1), SNPname)
    if (trans == "CCStrans") {
        # fixed K??
        res = ccstrans(2^allAint, 2^allBint)
        M = res$x
        S = res$y
    }
    else if (trans == "MAtrans") {
        # then prior??
        M = allAint - allBint
        S = (allAint + allBint)/2
    }
    else {
        stop(paste("bad transformation argument:", trans))
    }
    
    # divide CEL file up into chromosome pieces
    msList <- list()
    for (chri in mchr1) {
        currRows <- which(chrid == chri)
        numCurrRows <- length(currRows)
        if(numCurrRows == 0)
        {
            stop("Failed to find any probes on chromosome ", chri);
        }
        
        msList[[chri]]$M <- M[currRows]
        msList[[chri]]$S <- S[currRows]
    }
    
    msList
}

chunkIndices <- function(to, by) {
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

chunkFileName <- function(baseDir, celFileName, chrName, chunkSize, chunkIndex) {
    fileBase <- sub("\\..*$", "", celFileName)
    chunkFile <- paste(fileBase, "-", chrName, "-", chunkSize, "-", chunkIndex, ".RData", sep="")
    
    file.path(baseDir, chunkFile)
}
