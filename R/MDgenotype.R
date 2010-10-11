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
    celFiles = expandCelFiles(getwd()), confScoreThreshold = 1e-05,
    chromosomes = c(1:19, "X", "Y", "M"), cacheDir = tempdir(),
    retainCache = FALSE, verbose = FALSE, cluster = NULL,
    probesetChunkSize = 1000, processResultsFunction = NULL)
{
    #library("time")
    
    transformMethod <- match.arg(transformMethod)
    if(!inherits(snpProbeInfo, "data.frame") ||
       !all(c("probeIndex", "isAAllele", "snpId") %in% names(snpProbeInfo)))
    {
        stop("You must supply a \"snpProbeInfo\" data frame parameter which has ",
             "at a minimum the \"probeIndex\", \"isAAllele\" and \"snpId\" ",
             "components. Please see the help documentation for more details.")
    }
    
    if(!inherits(snpInfo, "data.frame") ||
       !all(c("snpId", "chrId") %in% names(snpInfo)))
    {
        stop("You must supply a \"snpInfo\" data frame parameter which has ",
             "at a minimum the \"snpId\" and \"chrId\" ",
             "components. Please see the help documentation for more details.")
    }
    
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
    
    # if we have a list (or data frame) pull out the file names and gender info
    # if it's in there
    isMale <- NULL
    if(is.list(celFiles))
    {
        if(!("fileName" %in% names(celFiles)))
        {
            stop("failed to find \"fileName\" component in the \"celFiles\" list")
        }
        
        if("isMale" %in% names(celFiles))
        {
            isMale <- celFiles$isMale
        }
        celFiles <- celFiles$fileName
    }
    
    nfile <- length(celFiles)
    if(nfile <= 1)
    {
        stop("Cannot successfully genotype with less than two CEL files")
    }
    
    # make sure that the chromosome vector is not numeric
    # TODO is toupper the right thing to do here?
    chromosomes <- toupper(as.character(chromosomes))
    
    allChr <- as.character(unique(snpInfo$chrId))
    allAutosomes <- setdiff(allChr, c("X", "Y", "M"))
    
    # we must infer gender if we don't yet have it and we need to
    # genotype either sex chromosome
    genderInferenceRequired <- is.null(isMale) && any(c("X", "Y") %in% chromosomes)
    if(genderInferenceRequired)
    {
        # in order to be able to infer gender we'll need per-array
        # mean intensity values for the X and Y chromosomes. We'll initialize
        # them as all 0's
        meanIntensityXPerArray <- rep(0, length(celFiles))
        meanIntensityYPerArray <- rep(0, length(celFiles))
        names(meanIntensityXPerArray) <- celFiles
        names(meanIntensityYPerArray) <- celFiles
        
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
    for(celfile in celFiles)
    {
        # wrapping this up as a local function def allows us to avoid doing the
        # normalization work unless it's necessary
        makeMsList <- function()
        {
            .normalizeCelFileByChr(
                celfile,
                verbose,
                snpProbeInfo,
                snpInfo,
                allChr,
                referenceDistribution,
                transformMethod)
        }
        msList <- NULL
        
        for(currChr in allChr)
        {
            for(chunkIndex in 1 : length(chrChunks[[currChr]]))
            {
                chunkFile <- .chunkFileName(cacheDir, "snp", celfile, currChr, probesetChunkSize, chunkIndex)
                chunkFileAlreadyExists <- file.exists(chunkFile)
                if(!chunkFileAlreadyExists)
                {
                    # if we haven't yet normalized the CEL file we'll have to
                    # do that now
                    if(is.null(msList))
                    {
                        msList <- makeMsList()
                    }
                    
                    chunk <- chrChunks[[currChr]][[chunkIndex]]
                    probesetIndices <- chunk$start : chunk$end
                    mChunk <- msList[[currChr]]$M[probesetIndices]
                    sChunk <- msList[[currChr]]$S[probesetIndices]
                    
                    save(mChunk, sChunk, file = chunkFile)
                }
                
                if(genderInferenceRequired)
                {
                    if(chunkFileAlreadyExists)
                    {
                        # we need to bring sChunk into scope
                        # TODO this is loading mChunk too which is unnecessary
                        #      if we save m and s to different files we can
                        #      avoid this
                        load(file = chunkFile)
                    }
                    
                    if(currChr %in% allAutosomes)
                    {
                        meanIntensityPerAutosome[currChr] <-
                            meanIntensityPerAutosome[currChr] + sum(as.double(sChunk))
                    }
                    else if(currChr == "X")
                    {
                    	meanIntensityXPerArray[celfile] <-
                            meanIntensityXPerArray[celfile] + sum(as.double(sChunk))
                    }
                    else if(currChr == "Y")
                    {
                        meanIntensityYPerArray[celfile] <-
                            meanIntensityYPerArray[celfile] + sum(as.double(sChunk))
                    }
                }
            }
        }
        
        rm(msList)
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
        
        isMale <- .inferGender(
            meanIntensityXPerArray,
            meanIntensityYPerArray,
            meanIntensityPerAutosome)
    }
    
    results <- list()
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
                chunkFile <- .chunkFileName(cacheDir, "snp", celFiles[i], chri, probesetChunkSize, chunkIndex)
                load(chunkFile)
                MM <- cbind(MM, mChunk)
                SS <- cbind(SS, sChunk)
                rm(mChunk, sChunk)
            }
            
            colnames(MM) <- celFiles
            colnames(SS) <- celFiles
            
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
                    ms = MM,
                    ss = SS,
                    hint = chrHint[chunkRange],
                    parIndices = which(chrPAR[chunkRange]),
                    trans = transformMethod,
                    isMale = isMale,
                    confScoreThreshold = confScoreThreshold)
                if(length(argLists) >= length(cluster) || chunkIndex == length(chrChunks[[chri]]))
                {
                    # parallel apply using snow then reset the arg list
                    chunkResultsList <- parLapply(cluster, argLists, .applyGenotypeAnyChrChunk)
                    for(i in 1 : length(chunkResultsList))
                    {
                        chunkResult <- chunkResultsList[[i]]
                        if(is.null(processResultsFunction))
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
                                colnames(chunkResult[[k]]) <- .fileBaseWithoutExtension(celFiles)
                            }
                            processResultsFunction(chunkProbesetInfo, chunkResult)
                        }
                    }
                    rm(chunkResult)
                    rm(chunkResultsList)
                    
                    argLists <- list()
                }
            }
            else
            {
                chunkResult <- .genotypeAnyChrChunk(
                    chr = chri,
                    ms = MM,
                    ss = SS,
                    hint = chrHint[chunkRange],
                    parIndices = which(chrPAR[chunkRange]),
                    trans = transformMethod,
                    isMale = isMale,
                    confScoreThreshold = confScoreThreshold)
                if(is.null(processResultsFunction))
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
                        colnames(chunkResult[[k]]) <- .fileBaseWithoutExtension(celFiles)
                    }
                    processResultsFunction(chunkProbesetInfo, chunkResult)
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
    if(!is.null(processResultsFunction))
    {
        results <- NULL
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
            colnames(results[[i]]) <- .fileBaseWithoutExtension(celFiles)
        }
    }
    
    # clean-up unless we were asked to retain the cache
    if(!retainCache)
    {
        if(verbose)
        {
            cat("cleaning up cache files before returning geno/vinotype results")
        }
        
        for(celfile in celFiles)
        {
            for(currChr in allChr)
            {
                for(chunkIndex in 1 : length(chrChunks[[currChr]]))
                {
                    chunkFile <- .chunkFileName(cacheDir, "snp", celfile, currChr, probesetChunkSize, chunkIndex)
                    if(file.exists(chunkFile))
                    {
                        file.remove(chunkFile)
                    }
                }
            }
        }
    }
    
    results
}

# this apply function is defined to take a list in order to allow us to take
# advantage of the snow package's apply functions
.applyGenotypeAnyChrChunk <- function(argList)
{
    .genotypeAnyChrChunk(
        argList$chr,
        argList$ms,
        argList$ss,
        argList$hint,
        argList$parIndices,
        argList$trans,
        argList$isMale,
        argList$confScoreThreshold)
}
