#########################################################################
#
# MouseDivGenotye.R
#
# Part of the MouseDivGenotype package
#
# This is the main function to genotype the Mouse Diversity Array.
#
#########################################################################

MouseDivGenotype = function(
    snpProbeInfo, snpInfo, referenceDistribution = NULL,
    transformMethod = c("CCStrans", "MAtrans"), celFiles = getwd(),
    chromosomes = c(1:19, "X", "Y", "M"), cacheDir = tempdir(),
    verbose = FALSE, cluster = NULL, probesetChunkSize=1000,
    processResultsFunction = NULL) {
    #library("time")
    
    transformMethod = match.arg(transformMethod)
    if(!inherits(snpProbeInfo, "data.frame") ||
       !all(c("probeIndex", "isAAllele", "snpId") %in% names(snpProbeInfo)))
    {
        stop("You must supply a \"snpProbeInfo\" data frame parameter which has ",
             "at a minimum the \"probeIndex\", \"isAAllele\" and \"snpId\"",
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
    
    isMale <- NULL
    filenames <- NULL
    if(is.data.frame(celFiles) || is.matrix(celFiles))
    {
        if(ncol(celFiles) >= 2)
        {
            gender <- toupper(celFiles[, 2])
            isMale <- gender == "MALE" | gender == "M"
            isFemale <- gender == "FEMALE" | gender == "F"

            # TODO in old code there was fallback to using vdist for gender
            #      inference. Ask Hyuna about the vdist fallback
            if(sum(isMale) + sum(isFemale) != nrow(filenames))
            {
                isMale <- NULL
                isFemale <- NULL
            }
        }
        filenames <- celFiles[1, ]
    }
    else if(is.character(celFiles))
    {
        # a simple function that will return all CEL file in a dir if given a
        # dir argument, or a single CEL file if that is what it's given
        celExpander <- function(name)
        {
            retVal <- c()
            if(!file.exists(name))
            {
                stop("failed to find file named \"", name, "\"")
            }
            else if(file.info(name)$isdir)
            {
                # if it's a dir expand it to all CEL file contents
                fileListing <- dir(name, full.names = TRUE)
                retVal <- fileListing[grep("*.CEL$", toupper(fileListing))]
            }
            else if(grep("*.CEL$", toupper(name)))
            {
                retVal <- name
            }
            
            retVal
        }
        
        filenames <- unlist(lapply(celFiles, celExpander))
    }
    else
    {
        stop("expected celFiles to be one of the following: ",
             "a matrix/data frame where the ",
             "first column contains CEL file names and the second column ",
             "contains CEL file gender (\"male\" or \"female\"), the name of ",
             "a directory containing CEL files, or a vector of CEL file names");
    }
    
    nfile <- length(filenames)
    if(nfile <= 1)
    {
        stop("Cannot successfully genotype with less than two CEL files")
    }
    
    # make sure that the chromosome vector is not numeric
    # TODO is toupper the right thing to do here?
    chromosomes <- toupper(as.character(chromosomes))
    
    allChr <- unique(snpInfo$chrId)
    allAutosomes <- setdiff(allChr, c("X", "Y", "M"))
    
    # we must infer gender if we don't yet have it and we need to
    # genotype either sex chromosome
    genderInferenceRequired <- is.null(isMale) && any(c("X", "Y") %in% chromosomes)
    if(genderInferenceRequired)
    {
        # in order to be able to infer gender we'll need per-array
        # mean intensity values for the X and Y chromosomes. We'll initialize
        # them as all 0's
        anyToZero <- function(...) { 0.0 }
        meanIntensityXPerArray <- sapply(filenames, anyToZero)
        meanIntensityYPerArray <- sapply(filenames, anyToZero)
        meanIntensityPerAutosome <- sapply(allAutosomes, anyToZero)
    }
    
    # figure out how many chunks there will be per chromosome. The reason we
    # are chunking up the CEL files is for reducing memory requirements and
    # allowing for more fine grained concurrency
    chrChunks <- list()
    for(currChr in allChr)
    {
        chrProbesetCount <- sum(snpInfo$chrId == currChr)
        chrChunks[[currChr]] <- chunkIndices(chrProbesetCount, probesetChunkSize)
    }
    
    # loop through all the CELL files. We need to read and normalize the CEL
    # file data and save it to file in chunks no bigger than probesetChunkSize
    for(celfile in filenames)
    {
        # wrapping this up as a local function def allows us to avoid doing the
        # normalization work unless it's necessary
        makeMsList <- function()
        {
            normalizeCelFileByChr(
                celfile,
                probesetChunkSize,
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
                chunkFile <- chunkFileName(cacheDir, celfile, currChr, probesetChunkSize, chunkIndex)
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
                        # TODO is as.double needed? (use is.double to check)
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
        
        isMale <- inferGender(
            meanIntensityXPerArray,
            meanIntensityYPerArray,
            meanIntensityPerAutosome)
    }
    
    # TODO storing in a chrResults list may be bad for mem requirements. should
    #      we optionally include a function argument that allows the user to
    #      write the results to file
    chrResults <- list()
    for (chri in chromosomes) {
        chrIndices <- which(snpInfo$chrId == chri)
        argLists <- list()
        chrResults[[chri]] <- list()
        
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
            if(verbose) {
                cat("geno/vinotyping chromosome ", chri,
                    " from probeset #", chunk$start,
                    " to probeset #", chunk$end, "\n", sep="")
            }
            
            #startTime <- getTime()
        
            # paste the chromosomes together for genotyping
            MM <- NULL
            SS <- NULL
            for (i in 1:nfile) {
                chunkFile <- chunkFileName(cacheDir, filenames[i], chri, probesetChunkSize, chunkIndex)
                load(chunkFile)
                MM <- cbind(MM, mChunk)
                SS <- cbind(SS, sChunk)
                rm(mChunk, sChunk)
            }
            
            colnames(MM) <- filenames
            colnames(SS) <- filenames
            
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
                    isMale = isMale)
                if(length(argLists) >= length(cluster) || chunkIndex == length(chrChunks[[chri]]))
                {
                    # parallel apply using snow then reset the arg list
                    chunkResultsList <- parLapply(cluster, argLists, applyGenotypeAnyChrChunk)
                    for(i in 1 : length(chunkResultsList))
                    {
                        chunkResult <- chunkResultsList[[i]]
                        if(is.null(processResultsFunction))
                        {
                            # we only accumulate results if there is no function
                            if(length(chrResults[[chri]]) == 0)
                            {
                                chrResults[[chri]] <- chunkResult
                            }
                            else
                            {
                                chrResults[[chri]] <- mapply(rbind, chrResults[[chri]], chunkResult, SIMPLIFY = FALSE)
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
                                colnames(chunkResult[[k]]) <- fileBaseWithoutExtension(filenames)
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
                chunkResult <- genotypeAnyChrChunk(
                    chr = chri,
                    ms = MM,
                    ss = SS,
                    hint = chrHint[chunkRange],
                    parIndices = which(chrPAR[chunkRange]),
                    trans = transformMethod,
                    isMale = isMale)
                if(is.null(processResultsFunction))
                {
                    # we only accumulate results if there is no function
                    if(length(chrResults[[chri]]) == 0)
                    {
                        chrResults[[chri]] <- chunkResult
                    }
                    else
                    {
                        chrResults[[chri]] <- mapply(rbind, chrResults[[chri]], chunkResult, SIMPLIFY = FALSE)
                    }
                }
                else
                {
                    chunkProbesetInfo <- snpInfo[chrIndices[chunkRange], ]
                    
                    for(k in 1 : length(chunkResult))
                    {
                        colnames(chunkResult[[k]]) <- fileBaseWithoutExtension(filenames)
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
    # we have not accumulated any results to return
    if(!is.null(processResultsFunction))
    {
        chrResults <- NULL
    }
    
    chrResults
}

# this apply function is defined to take a list in order to allow us to take
# advantage of the snow package's apply functions
applyGenotypeAnyChrChunk <- function(argList)
{
    genotypeAnyChrChunk(
        argList$chr,
        argList$ms,
        argList$ss,
        argList$hint,
        argList$parIndices,
        argList$trans,
        argList$isMale)
}
