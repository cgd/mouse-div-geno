#########################################################################
#
# MouseDivGenotye.R
#
# Part of the MouseDivGenotype package
#
# This is the main function to genotype the Mouse Diversity Array.
#
#########################################################################

mouseDivGenotype <- function(
        snpProbeInfo, snpInfo, referenceDistribution = NULL,
        transformMethod = c("CCStrans", "MAtrans"),
        celFiles = expandCelFiles(getwd()), isMale = NULL, confScoreThreshold = 1e-05,
        chromosomes = c(1:19, "X", "Y", "M"),
        verbose = FALSE, cluster = NULL,
        probesetChunkSize = 1000, outputDir = NULL, outputFilePrefix = "mouseDivResults_")
{
    transformMethod <- match.arg(transformMethod)
    if(!inherits(snpProbeInfo, "data.frame") ||
       !all(c("probeIndex", "isAAllele", "snpId") %in% names(snpProbeInfo)))
    {
        stop("You must supply a \"snpProbeInfo\" data frame parameter which has ",
                "at a minimum the \"probeIndex\", \"isAAllele\" and \"snpId\" ",
                "components. Please see the help documentation for more details.")
    }
    snpProbeInfo$snpId <- as.factor(snpProbeInfo$snpId)
    
    snpIntensities <- .readSnpIntensitiesFromCEL(
            celFiles                = celFiles,
            snpProbeInfo            = snpProbeInfo,
            referenceDistribution   = referenceDistribution,
            verbose                 = verbose)
    .mouseDivGenotypeInternal(
            snpIntensities      = snpIntensities,
            snpInfo             = snpInfo,
            transformMethod     = transformMethod,
            isMale              = isMale,
            confScoreThreshold  = confScoreThreshold,
            chromosomes         = chromosomes,
            verbose             = verbose,
            cluster             = cluster,
            probesetChunkSize   = probesetChunkSize,
            outputDir           = outputDir,
            outputFilePrefix    = outputFilePrefix)
}

mouseDivGenotypeTab <- function(
        snpIntensityFile, snpInfo,
        transformMethod = c("CCStrans", "MAtrans"),
        isMale = NULL, confScoreThreshold = 1e-05,
        chromosomes = c(1:19, "X", "Y", "M"),
        verbose = FALSE, cluster = NULL,
        probesetChunkSize = 1000, outputDir = NULL, outputFilePrefix = "mouseDivResults_")
{
    transformMethod <- match.arg(transformMethod)
    
    snpIntensities <- .readSnpIntensitiesFromTab(
            tabFile             = snpIntensityFile,
            verbose             = verbose)
    .mouseDivGenotypeInternal(
            snpIntensities      = snpIntensities,
            snpInfo             = snpInfo,
            transformMethod     = transformMethod,
            isMale              = isMale,
            confScoreThreshold  = confScoreThreshold,
            chromosomes         = chromosomes,
            verbose             = verbose,
            cluster             = cluster,
            probesetChunkSize   = probesetChunkSize,
            outputDir           = outputDir,
            outputFilePrefix    = outputFilePrefix)
}

.mouseDivGenotypeInternal <- function(
    snpIntensities, snpInfo,
    transformMethod = c("CCStrans", "MAtrans"),
    isMale = NULL, confScoreThreshold = 1e-05,
    chromosomes = c(1:19, "X", "Y", "M"), cacheDir = tempdir(),
    retainCache = FALSE, verbose = FALSE, cluster = NULL,
    probesetChunkSize = 1000, outputDir = NULL, outputFilePrefix = "mouseDivResults_")
{
    transformMethod <- match.arg(transformMethod)
    if(transformMethod == "CCStrans")
    {
        snpIntensities <- .lazyApply(.ccsTransformSample, snpIntensities)
    }
    else if(transformMethod == "MAtrans")
    {
        snpIntensities <- .lazyApply(.maTransformSample, snpIntensities)
    }
    else
    {
        stop("transformMethod parameter is not valid")
    }
    
    if(!inherits(snpInfo, "data.frame") ||
       !all(c("snpId", "chrId") %in% names(snpInfo)))
    {
        stop("You must supply a \"snpInfo\" data frame parameter which has ",
             "at a minimum the \"snpId\" and \"chrId\" ",
             "components. Please see the help documentation for more details.")
    }
    snpInfo$snpId <- as.factor(snpInfo$snpId)
    
    # it is an error if isInPAR is set to TRUE in a non-X chromosome
    if(!is.null(snpInfo$isInPAR))
    {
        parChrs <- unique(snpInfo$chrId[snpInfo$isInPar])
        if(!all(parChrs == "X"))
        {
            stop("snpInfo$isInPar should only ever be TRUE on the \"X\" ",
                 "chromosome, but TRUE isInPar values were found on chromosomes: ",
                 paste(parChrs, collapse=", "))
        }
    }
    
    # we only need the snow library if the cluster is non-null
    if(!is.null(cluster) && !require("snow"))
    {
        stop("failed to load the snow library")
    }
    
    # make sure that the chromosome vector is not numeric
    # TODO is toupper the right thing to do here?
    chromosomes <- toupper(as.character(chromosomes))
    
    allChr <- as.character(unique(snpInfo$chrId))
    allAutosomes <- setdiff(allChr, c("X", "Y", "M"))
    
    # we must infer gender if we don't yet have it and we need to
    # genotype either sex chromosome
    genderInferenceRequired <- is.null(isMale)
    if(genderInferenceRequired)
    {
        meanIntensityXPerArray <- NULL
        meanIntensityYPerArray <- NULL
        
        meanIntensityPerAutosome <- rep(0, length(allAutosomes))
        names(meanIntensityPerAutosome) <- allAutosomes
    }
    
    # figure out how many chunks there will be per chromosome. The reason we
    # are chunking up the CEL files is for reducing memory requirements and
    # allowing for more fine grained concurrency
    chrChunks <- list()
    for(currChr in allChr)
    {
        chrProbesetCount <- sum(snpInfo$chrId == currChr)
        chrChunks[[currChr]] <- .chunkIndices(chrProbesetCount, probesetChunkSize)
    }
    
    # loop through all the CELL files. We need to read and normalize the CEL
    # file data and save it to file in chunks no bigger than probesetChunkSize
    sampleNames <- character()
    while(!is.null(snpIntensities))
    {
        head <- snpIntensities$head
        sampleNames <- c(sampleNames, head$sampleName)
        
        if(genderInferenceRequired)
        {
            meanIntensityXPerArray[head$sampleName] <- 0
            meanIntensityYPerArray[head$sampleName] <- 0
        }
        
        for(currChr in allChr)
        {
            chrSnpInfo <- snpInfo[snpInfo$chrId == currChr, ]
            chrMatch <- match(chrSnpInfo$snpId, rownames(head$sampleData))
            chrMatch <- chrMatch[!is.na(chrMatch)]
            chrData <- head$sampleData[chrMatch, ]
            for(chunkIndex in 1 : length(chrChunks[[currChr]]))
            {
                chunk <- chrChunks[[currChr]][[chunkIndex]]
                probesetIndices <- chunk$start : chunk$end
                mChunk <- chrData[ , "intensityConts"][probesetIndices]
                sChunk <- chrData[ , "intensityAvgs"][probesetIndices]
                
                if(genderInferenceRequired)
                {
                    if(currChr %in% allAutosomes)
                    {
                        meanIntensityPerAutosome[currChr] <-
                            meanIntensityPerAutosome[currChr] + sum(as.double(sChunk))
                    }
                    else if(currChr == "X")
                    {
                    	meanIntensityXPerArray[head$sampleName] <-
                            meanIntensityXPerArray[head$sampleName] + sum(as.double(sChunk))
                    }
                    else if(currChr == "Y")
                    {
                        meanIntensityYPerArray[head$sampleName] <-
                            meanIntensityYPerArray[head$sampleName] + sum(as.double(sChunk))
                    }
                }
                
                # cache the data
                chunkFile <- .chunkFileName(cacheDir, "snp", head$sampleName, currChr, probesetChunkSize, chunkIndex)
                save(mChunk, sChunk, file = chunkFile)
            }
        }
        
        snpIntensities <- snpIntensities$tail()
    }
    
    nfile <- length(sampleNames)
    if(nfile <= 1)
    {
        stop("Cannot successfully genotype with less than two samples")
    }
    
    if(!is.null(isMale) && length(isMale) != nfile)
    {
        stop("the length of isMale must match the sample count")
    }
    
    # determines which arrays are male and which are female if we need to
    if(genderInferenceRequired)
    {
        # up to this point we've just been summing up intensity values. Now
        # we want to take the mean
        probesetCountPerAutosome <- sapply(
            allAutosomes,
            function(currChr) {sum(snpInfo$chrId == currChr)})
        
        meanIntensityPerAutosome <-
            meanIntensityPerAutosome / probesetCountPerAutosome
        meanIntensityXPerArray <- meanIntensityXPerArray / sum(snpInfo$chrId == "X")
        meanIntensityYPerArray <- meanIntensityYPerArray / sum(snpInfo$chrId == "Y")
        
        isMale <- inferGender(
            meanIntensityXPerArray,
            meanIntensityYPerArray,
            meanIntensityPerAutosome)
    }
    
    results <- list()
    outFileConnections <- NULL
    for (chri in chromosomes)
    {
        chrIndices <- which(snpInfo$chrId == chri)
        argLists <- list()
        
        # pull out the hints for this chromosome (if they exist)
        chrHint <- NULL
        if(!is.null(snpInfo$snpHetHint))
        {
            chrHint <- snpInfo$snpHetHint[chrIndices]
        }
        
        chrPAR <- logical(0)
        if(!is.null(snpInfo$isInPAR))
        {
            chrPAR <- snpInfo$isInPAR[chrIndices]
        }
        
        for(chunkIndex in 1 : length(chrChunks[[chri]]))
        {
            chunk <- chrChunks[[chri]][[chunkIndex]]
            if(verbose)
            {
                cat("geno/vinotyping chromosome ", chri,
                    " from probeset #", chunk$start,
                    " to probeset #", chunk$end, "\n", sep="")
            }
            
            #startTime <- getTime()
        
            # paste the chromosomes together for genotyping
            MM <- NULL
            SS <- NULL
            for (i in 1:nfile)
            {
                chunkFile <- .chunkFileName(cacheDir, "snp", sampleNames[i], chri, probesetChunkSize, chunkIndex)
                load(chunkFile)
                MM <- cbind(MM, mChunk)
                SS <- cbind(SS, sChunk)
                rm(mChunk, sChunk)
            }
            
            colnames(MM) <- sampleNames
            colnames(SS) <- sampleNames
            
            #cat("time it took us to get to genotypethis\n")
            #timeReport(startTime)
            
            #startTime <- getTime()
            chunkRange <- chunk$start : chunk$end
            if(length(cluster) >= 1)
            {
                # the arg list is used to accumulate arguments until we're
                # ready to execute them in parallel
                argLists[[length(argLists) + 1]] <- list(
                    chr = chri,
                    intensityConts = MM,
                    intensityAvgs = SS,
                    hint = chrHint[chunkRange],
                    parIndices = which(chrPAR[chunkRange]),
                    transformMethod = transformMethod,
                    isMale = isMale,
                    confScoreThreshold = confScoreThreshold)
                if(length(argLists) >= length(cluster) || chunkIndex == length(chrChunks[[chri]]))
                {
                    # parallel apply using snow then reset the arg list
                    chunkResultsList <- parLapply(cluster, argLists, .applyGenotypeAnyChrChunk)
                    for(i in 1 : length(chunkResultsList))
                    {
                        chunkResult <- chunkResultsList[[i]]
                        if(is.null(outputDir))
                        {
                            # we only accumulate results if there is no function
                            if(length(results) == 0)
                            {
                                results <- chunkResult
                            }
                            else
                            {
                                results <- mapply(rbind, results, chunkResult, SIMPLIFY = FALSE)
                            }
                        }
                        else
                        {
                            # this is a bit tricky, but we need to grab the chunk object
                            # that matches up correctly with the chunk result that
                            # we are iterating over
                            clustChunk <- chrChunks[[chri]][[chunkIndex - (length(chunkResultsList) - i)]]
                            clustChunkRange <- clustChunk$start : clustChunk$end
                            chunkProbesetInfo <- snpInfo[chrIndices[clustChunkRange], ]
                            
                            for(k in 1 : length(chunkResult))
                            {
                                colnames(chunkResult[[k]]) <- sampleNames
                                rownames(chunkResult[[k]]) <- chunkProbesetInfo$snpId
                            }
                            
                            if(is.null(outFileConnections))
                            {
                                outFileConnections <- .initFlatFileConnections(
                                        outputDir,
                                        chunkResult,
                                        prefix = outputFilePrefix)
                            }
                            .writeResultsToFlatFile(outFileConnections, chunkResult)
                        }
                    }
                    rm(chunkResult)
                    rm(chunkResultsList)
                    
                    argLists <- list()
                }
            }
            else
            {
                chunkResult <- genotypeAnyChrChunk(
                    chr = chri,
                    intensityConts = MM,
                    intensityAvgs = SS,
                    hint = chrHint[chunkRange],
                    parIndices = which(chrPAR[chunkRange]),
                    transformMethod = transformMethod,
                    isMale = isMale,
                    confScoreThreshold = confScoreThreshold)
                if(is.null(outputDir))
                {
                    # we only accumulate results if there is no function
                    if(length(results) == 0)
                    {
                        results <- chunkResult
                    }
                    else
                    {
                        results <- mapply(rbind, results, chunkResult, SIMPLIFY = FALSE)
                    }
                }
                else
                {
                    chunkProbesetInfo <- snpInfo[chrIndices[chunkRange], ]
                    
                    for(k in 1 : length(chunkResult))
                    {
                        colnames(chunkResult[[k]]) <- sampleNames
                        rownames(chunkResult[[k]]) <- chunkProbesetInfo$snpId
                    }
                    
                    if(is.null(outFileConnections))
                    {
                        outFileConnections <- .initFlatFileConnections(
                                outputDir,
                                chunkResult,
                                prefix = outputFilePrefix)
                    }
                    .writeResultsToFlatFile(outFileConnections, chunkResult)
                }
                rm(chunkResult)
            }
            
            #cat("time it took us to genotype the current chunk:\n")
            #timeReport(startTime)
        }
        rm(argLists)
    }
    
    # if the user asked us to use a function to process the results that means
    # we have not accumulated any results to return, otherwise we have a bit
    # of post-processing to do to make sure that the matrices that we return
    # to the user match up nicely with what they passed in
    if(!is.null(outputDir))
    {
        results <- NULL
        for(con in outFileConnections)
        {
            flush(con)
            close(con)
        }
    }
    else
    {
        # using list/unlist because c() will coerce factors into an integer
        snpIdsByChr <- list()
        for(chri in chromosomes)
        {
            snpIdsByChr[[chri]] <- snpInfo$snpId[snpInfo$chrId == chri]
        }
        snpIdsByChr <- unlist(snpIdsByChr, use.names = FALSE)
        
        snpIdsInOrder <- snpInfo$snpId[snpInfo$chrId %in% chromosomes]
        snpOrdering <- match(snpIdsInOrder, snpIdsByChr)
        if(any(is.na(snpOrdering)))
        {
            stop("internal error: failed to match up SNP IDs in results")
        }
        
        # reorder the SNPs so that they match up with the snpInfo frame that
        # was passed in
        for(i in 1 : length(results))
        {
            results[[i]] <- results[[i]][snpOrdering, ]
            rownames(results[[i]]) <- snpIdsInOrder
            colnames(results[[i]]) <- sampleNames
        }
    }
    
    # clean-up unless we were asked to retain the cache
    if(!retainCache)
    {
        if(verbose)
        {
            cat("cleaning up cache files before returning geno/vinotype results")
        }
        
        for(sampleName in sampleNames)
        {
            for(currChr in allChr)
            {
                for(chunkIndex in 1 : length(chrChunks[[currChr]]))
                {
                    chunkFile <- .chunkFileName(cacheDir, "snp", sampleName, currChr, probesetChunkSize, chunkIndex)
                    if(file.exists(chunkFile))
                    {
                        file.remove(chunkFile)
                    }
                }
            }
        }
    }
    
    results$isMale <- isMale
    results
}

# this apply function is defined to take a list in order to allow us to take
# advantage of the snow package's apply functions
.applyGenotypeAnyChrChunk <- function(argList)
{
    genotypeAnyChrChunk(
        argList$chr,
        argList$intensityConts,
        argList$intensityAvgs,
        argList$hint,
        argList$parIndices,
        argList$transformMethod,
        argList$isMale,
        argList$confScoreThreshold)
}
