#########################################################################
#
# MouseDivGenotye.R
#
# Part of the MouseDivGenotype package
#
# This is the main function to genotype the Mouse Diversity Array.
#
#########################################################################

#TODO REMOVE ME KEITH!!!
library("snow")

library("affyio")
library("preprocessCore")
library("cluster")
#library("time")

MouseDivGenotype = function(allid, ABid, chrid, CGFLcorrection = NULL, 
    reference = NULL, hint = NULL, trans = c("CCStrans", "MAtrans"), celnamefile = NULL, 
    mchr = c(1:19, "X", "Y", "M"), celfiledir, outfiledir, subset = FALSE, 
    verbose = FALSE, cluster = NULL, probesetChunkSize=1000) {
    
    trans = match.arg(trans)
    if (missing(celfiledir)) 
        celfiledir = getwd()
    if (missing(outfiledir)) 
        outfiledir = getwd()
    if (missing(allid) | missing(ABid)) 
        stop("No CDF file information")
    SNPname = ABid$SNPname
    Aid = ABid$allAid
    Bid = ABid$allBid
    rm(ABid)
    mpos = chrid$mpos
    
    # TODO this should be handled by which param is passed in
    if (subset) 
        chrid = chrid$schrid
    else
        chrid = chrid$chrid
    
    # TODO we shouldn't be setting the working dir
    setwd(celfiledir)
    
    # TODO rework this section as:
    #       - rename celnamefile to celfiles
    #       - if celfiles is a vector of strings it is treated as a .CEL filename list
    #       - if celfiles is a data.frame it is assumed to have isMale and file
    #         factors
    #       - if celfiles is a single string it is treated as a dir of CEL files
    isMale <- NULL
    if(length(celnamefile) == 0)
    {
        # since no celnamefile is provided assume all CEL files come from the
        # working directory
        tmp = dir()
        isCel = apply(matrix(tmp, ncol = 1), 1, function(x) length(grep(".CEL", x)) > 0)
        filenames = tmp[isCel]
    }
    else
    {
        filenames <- read.delim(celnamefile, header = TRUE)
        if(ncol(filenames) > 1)
        {
            gender <- toupper(filenames[, 2])
            isMale <- gender == "MALE"
            isFemale <- gender == "FEMALE"
            
            # TODO in old code there was fallback to using vdist for gender
            #      inference. Ask Hyuna about the vdist fallback
            if(sum(isMale) + sum(isFemale) != nrow(filenames))
            {
                isMale <- NULL
                isFemale <- NULL
            }
        }
        
        filenames <- filenames[, 1]
    }
    nfile = length(filenames)
    
    # make sure that the chromosome vector is not numeric
    # TODO is toupper the right thing to do here?
    mchr <- toupper(as.character(mchr))
    
    allChr <- unique(chrid)
    allAutosomes <- setdiff(allChr, c("X", "Y", "M"))
    
    # we must infer gender if we don't yet have it and we need to
    # genotype either sex chromosome
    genderInferenceRequired <- is.null(isMale) && any(c("X", "Y") %in% mchr)
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
        chrProbesetCount <- sum(chrid == currChr)
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
                allid,
                Aid,
                Bid,
                CGFLcorrection,
                reference,
                SNPname,
                trans,
                allChr,
                chrid)
        }
        msList <- NULL
        
        for(currChr in allChr)
        {
            for(chunkIndex in 1 : length(chrChunks[[currChr]]))
            {
                chunkFile <- chunkFileName(outfiledir, celfile, currChr, probesetChunkSize, chunkIndex)
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
                        # TODO this is loading mChunk to which is unnecessary
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
            function(currChr) {sum(chrid == currChr)})
        
        meanIntensityPerAutosome <-
            meanIntensityPerAutosome / probesetCountPerAutosome
        meanIntensityXPerArray <- meanIntensityXPerArray / sum(chrid == "X")
        meanIntensityYPerArray <- meanIntensityYPerArray / sum(chrid == "Y")
        
        isMale <- inferGender(
            meanIntensityXPerArray,
            meanIntensityYPerArray,
            meanIntensityPerAutosome)
    }
    
    # TODO storing in a chrResults list may be bad for mem requirements. should
    #      we optionally include a function argument that allows the user to
    #      write the results to file
    chrResults <- list()
    for (chri in mchr) {
        argLists <- list()
        chrResults[[chri]] <- list()
        
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
                chunkFile <- chunkFileName(outfiledir, filenames[i], chri, probesetChunkSize, chunkIndex)
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
            chunkIndices <- chunk$start : chunk$end
            if(length(cluster) >= 1)
            {
                # the arg list is used to accumulate arguments until we're
                # ready to execute them in parallel
                argLists[[length(argLists) + 1]] <- list(
                    chr = chri,
                    ms = MM,
                    ss = SS,
                    hint = hint[[chri]][chunkIndices],
                    trans = trans,
                    isMale = isMale)
                if(length(argLists) >= length(cluster) || chunkIndex == length(chrChunks[[chri]]))
                {
                    # parallel apply using snow then reset the arg list
                    chunkResultsList <- parLapply(cluster, argLists, applyGenotypeAnyChrChunk)
                    for(chunkResult in chunkResultsList)
                    {
                        if(length(chrResults[[chri]]) == 0)
                        {
                            chrResults[[chri]] <- chunkResult
                        }
                        else
                        {
                            chrResults[[chri]] <- mapply(rbind, chrResults[[chri]], chunkResult, SIMPLIFY = FALSE)
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
                    hint = hint[[chri]][chunkIndices],
                    trans = trans,
                    isMale = isMale)
                if(length(chrResults[[chri]]) == 0)
                {
                    chrResults[[chri]] <- chunkResult
                }
                else
                {
                    chrResults[[chri]] <- mapply(rbind, chrResults[[chri]], chunkResult, SIMPLIFY = FALSE)
                }
                rm(chunkResult)
            }
            
            #cat("time it took us to genotype the current chunk:\n")
            #timeReport(startTime)
        }
        rm(argLists)
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
        argList$trans,
        argList$isMale)
}
