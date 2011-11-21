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
        snpIntensities,
        snpProbeInfo, snpInfo, referenceDistribution = NULL,
        transformFunction = ccsTransform,
        isMale = NULL, confScoreThreshold = 1e-05,
        chromosomes = c(1:19, "X", "Y", "M"),
        cacheDir = tempdir(), retainCache = FALSE,
        numCores = NULL,
        probesetChunkSize = 1000, outputDir = NULL, outputFilePrefix = "mouseDivResults_",
        logFile = NULL) {

    if(inherits(logFile, "character")) {
        logFile <- file(logFile, "wt")
    } else if(!is.null(logFile) && !inherits(logFile, "connection")) {
        stop("the logFile parameter should be either NULL, a file name, or a connection")
    }
    
    if(!inherits(snpProbeInfo, "data.frame") ||
        !all(c("probeIndex", "isAAllele", "snpId") %in% names(snpProbeInfo))) {
        stop("You must supply a \"snpProbeInfo\" data frame parameter which has ",
            "at a minimum the \"probeIndex\", \"isAAllele\" and \"snpId\" ",
            "components. Please see the help documentation for more details.")
    }
    snpProbeInfo$snpId <- as.factor(snpProbeInfo$snpId)
    
    if(!inherits(snpIntensities, "snpIntensities")) {
        stop(
            "mouseDivGenotype(...) expects the snpIntensities argument to ",
            "inherit from the snpIntensities class")
    }
    
    # TODO what I'm doing here is very inefficient. It's probably best to have
    # a version of .mouseDivGenotypeInternal which is speciallized just for
    # snpIntensity objects rather than turning them into lazy lists
    colCount <- ncol(snpIntensities)
    lazyReadNext <- function(index) {
        headFunc <- function() {
            snpIntensities[, index, drop=FALSE]
        }
        
        tailFunc <-
            if(index == colCount) {
                NULL
            } else {
                function() {lazyReadNext(index + 1)}
            }
        
        list(head = headFunc, tail = tailFunc)
    }
    
    lazySnpIntensities <-
        if(colCount) {
            function() {lazyReadNext(1)}
        } else {
            NULL
        }
    
    .mouseDivGenotypeInternal(
        snpIntensities      = lazySnpIntensities,
        snpInfo             = snpInfo,
        transformFunction   = transformFunction,
        isMale              = isMale,
        confScoreThreshold  = confScoreThreshold,
        chromosomes         = chromosomes,
        cacheDir            = cacheDir,
        retainCache         = retainCache,
        numCores            = numCores,
        probesetChunkSize   = probesetChunkSize,
        outputDir           = outputDir,
        outputFilePrefix    = outputFilePrefix,
        logFile             = logFile)
}

mouseDivGenotypeCEL <- function(
        celFiles = getwd(),
        snpProbeInfo, snpInfo, referenceDistribution = NULL,
        transformFunction = ccsTransform,
        isMale = NULL, confScoreThreshold = 1e-05,
        chromosomes = c(1:19, "X", "Y", "M"),
        cacheDir = tempdir(), retainCache = FALSE,
        numCores = NULL,
        probesetChunkSize = 1000, outputDir = NULL, outputFilePrefix = "mouseDivResults_",
        logFile = NULL) {

    if(inherits(logFile, "character")) {
        logFile <- file(logFile, "wt")
    } else if(!is.null(logFile) && !inherits(logFile, "connection")) {
        stop("the logFile parameter should be either NULL, a file name, or a connection")
    }
    
    if(!inherits(snpProbeInfo, "data.frame") ||
       !all(c("probeIndex", "isAAllele", "snpId") %in% names(snpProbeInfo))) {
        stop("You must supply a \"snpProbeInfo\" data frame parameter which has ",
                "at a minimum the \"probeIndex\", \"isAAllele\" and \"snpId\" ",
                "components. Please see the help documentation for more details.")
    }
    snpProbeInfo$snpId <- as.factor(snpProbeInfo$snpId)
    
    celFiles <- .expandCelFiles(celFiles)
    snpIntensities <- .readSnpIntensitiesFromCEL(
            celFiles                = celFiles,
            snpProbeInfo            = snpProbeInfo,
            referenceDistribution   = referenceDistribution,
            numCores                = numCores,
            logFile                 = logFile)
    .mouseDivGenotypeInternal(
            snpIntensities      = snpIntensities,
            snpInfo             = snpInfo,
            transformFunction   = transformFunction,
            isMale              = isMale,
            confScoreThreshold  = confScoreThreshold,
            chromosomes         = chromosomes,
            cacheDir            = cacheDir,
            retainCache         = retainCache,
            numCores            = numCores,
            probesetChunkSize   = probesetChunkSize,
            outputDir           = outputDir,
            outputFilePrefix    = outputFilePrefix,
            logFile             = logFile)
}

mouseDivGenotypeTab <- function(
        snpIntensityFile, snpInfo,
        transformFunction = ccsTransform,
        isMale = NULL, confScoreThreshold = 1e-05,
        chromosomes = c(1:19, "X", "Y", "M"),
        cacheDir = tempdir(), retainCache = FALSE,
        numCores = NULL,
        probesetChunkSize = 1000, outputDir = NULL, outputFilePrefix = "mouseDivResults_",
        logFile = NULL) {

    if(inherits(logFile, "character")) {
        logFile <- file(logFile, "wt")
    } else if(!is.null(logFile) && !inherits(logFile, "connection")) {
        stop("the logFile parameter should be either NULL, a file name, ",
             "or a connection")
    }
    
    snpIntensities <- .readSnpIntensitiesFromTab(
            tabFile             = snpIntensityFile,
            logFile             = logFile)
    .mouseDivGenotypeInternal(
            snpIntensities      = snpIntensities,
            snpInfo             = snpInfo,
            transformFunction   = transformFunction,
            isMale              = isMale,
            confScoreThreshold  = confScoreThreshold,
            chromosomes         = chromosomes,
            cacheDir            = cacheDir,
            retainCache         = retainCache,
            numCores            = numCores,
            probesetChunkSize   = probesetChunkSize,
            outputDir           = outputDir,
            outputFilePrefix    = outputFilePrefix,
            logFile             = logFile)
}

.mouseDivGenotypeInternal <- function(
        snpIntensities, snpInfo,
        transformFunction = ccsTransform,
        isMale = NULL, confScoreThreshold = 1e-05,
        chromosomes = c(1:19, "X", "Y", "M"),
        cacheDir = tempdir(), retainCache = FALSE,
        numCores = NULL,
        probesetChunkSize = 1000, outputDir = NULL, outputFilePrefix = "mouseDivResults_",
        logFile = NULL) {

    snpIntensities <- .lazyApply(transformFunction, snpIntensities)
    
    if(!inherits(snpInfo, "data.frame") || !all(c("snpId", "chrId") %in% names(snpInfo))) {
       stop("You must supply a \"snpInfo\" data frame parameter which has ",
             "at a minimum the \"snpId\" and \"chrId\" ",
             "components. Please see the help documentation for more details.")
    }
    snpInfo$snpId <- as.factor(snpInfo$snpId)
    
    if(inherits(logFile, "character")) {
        logFile <- file(logFile, "wt")
    } else if(!is.null(logFile) && !inherits(logFile, "connection")) {
        stop("the logFile parameter should be either NULL, a file name, ",
             "or a connection")
    }
    
    logOn <- !is.null(logFile)
    logfn <- function(fmt, ...) {
        if(logOn) {
            logfnTry <- function() {
                cat(sprintf(fmt, ...), file=logFile)
                cat("\n", file=logFile)
                flush(logFile)
            }
            onError <- function(e) {
                warning("failed to log message")
            }
            tryCatch(logfnTry(), error=onError)
        }
    }
    
    # it is an error if isInPAR is set to TRUE in a non-X chromosome
    if(!is.null(snpInfo$isInPAR)) {
        parChrs <- unique(snpInfo$chrId[snpInfo$isInPAR])
        if(!all(parChrs == "X")) {
            stop("snpInfo$isInPAR should only ever be TRUE on the \"X\" ",
                 "chromosome, but TRUE isInPAR values were found on chromosomes: ",
                 paste(parChrs, collapse=", "))
        }
    }
    
    if(is.null(numCores) || numCores != 1) {
        if(!require("multicore", character.only=T)) {
            numCores <- 1
        } else if(is.null(numCores)) {
            numCores <- multicore:::detectCores()
        }
    }
    
    # make sure that the chromosome vector is not numeric
    # TODO is toupper the right thing to do here?
    chromosomes <- toupper(as.character(chromosomes))
    
    allChr <- as.character(unique(snpInfo$chrId))
    allAutosomes <- setdiff(allChr, c("X", "Y", "M"))
    logfn("genotyping the following chromosomes: %s", paste(chromosomes, collapse = ", "))
    
    # we must infer gender if we don't yet have it and we need to
    # genotype either sex chromosome
    genderInferenceRequired <- is.null(isMale)
    if(genderInferenceRequired) {
        meanIntensityXPerArray <- NULL
        meanIntensityYPerArray <- NULL
        
        meanIntensityPerAutosome <- rep(0, length(allAutosomes))
        names(meanIntensityPerAutosome) <- allAutosomes
    }
    
    # figure out how many chunks there will be per chromosome. The reason we
    # are chunking up the CEL files is for reducing memory requirements and
    # allowing for more fine grained concurrency
    chrChunks <- list()
    for(currChr in allChr) {
        chrProbesetCount <- sum(snpInfo$chrId == currChr)
        chrChunks[[currChr]] <- .chunkIndices(chrProbesetCount, probesetChunkSize)
    }
    
    # loop through all the CELL files. We need to read and normalize the CEL
    # file data and save it to file in chunks no bigger than probesetChunkSize
    sampleNames <- character()
    normSample <- function(sample) {
        logfn("preprocessing intensities from %s ", sample$sampleNames)
        if(ncol(sample) != 1) {
            stop("expected there to be exactly 1 sample but there are ", ncol(sample))
        }
        
        if(genderInferenceRequired) {
            sumIntenX <- 0
            sumIntenY <- 0
            sumIntenPerAuto <- rep(0, length(allAutosomes))
            names(sumIntenPerAuto) <- allAutosomes
        }
        
        for(currChr in allChr) {
            chrSnpInfo <- snpInfo[snpInfo$chrId == currChr, , drop=FALSE]
            chrMatch <- match(chrSnpInfo$snpId, rownames(sample))
            chrMatch <- chrMatch[!is.na(chrMatch)]
            chrData <- sample[chrMatch, , drop=FALSE]
            for(chunkIndex in 1 : length(chrChunks[[currChr]])) {
                chunk <- chrChunks[[currChr]][[chunkIndex]]
                probesetIndices <- chunk$start : chunk$end
                chunkData <- chrData[probesetIndices, , drop=FALSE]
                
                if(genderInferenceRequired) {
                    if(currChr %in% allAutosomes) {
                        sumIntenPerAuto[currChr] <-
                            sumIntenPerAuto[currChr] + sum(as.double(chunkData$snpAverages))
                    } else if(currChr == "X") {
                        sumIntenX <- sumIntenX + sum(as.double(chunkData$snpAverages))
                    } else if(currChr == "Y") {
                        sumIntenY <- sumIntenY + sum(as.double(chunkData$snpAverages))
                    }
                }
                
                # cache the data
                chunkFile <- .chunkFileName(cacheDir, "snp", sample$sampleNames, currChr, probesetChunkSize, chunkIndex)
                save(chunkData, file=chunkFile)
            }
        }
        
        if(genderInferenceRequired) {
            list(
                sampleName=sample$sampleNames,
                sumIntenPerAuto=sumIntenPerAuto,
                sumIntenX=sumIntenX,
                sumIntenY=sumIntenY)
        } else {
            NULL
        }
    }
    
    while(!is.null(snpIntensities)) {
        parIntens <- list()
        for(i in seq_len(numCores)) {
            if(!is.null(snpIntensities)) {
                curr <- snpIntensities()
                currHead <- curr$head()
                parIntens[[currHead$sampleNames]] <- currHead
                snpIntensities <- curr$tail
            }
        }
        
        applyResult <-
            if(numCores >= 2 && length(parIntens) >= 2) {
                mclapply(parIntens, normSample, mc.cores=numCores)
            } else {
                lapply(parIntens, normSample)
            }
        
        for(currRes in applyResult) {
            if(genderInferenceRequired) {
                for(currChr in names(currRes$sumIntenPerAuto)) {
                    meanIntensityPerAutosome[currChr] <-
                        meanIntensityPerAutosome[currChr] + currRes$sumIntenPerAuto[currChr]
                }
                meanIntensityXPerArray[currRes$sampleName] <- currRes$sumIntenX
                meanIntensityYPerArray[currRes$sampleName] <- currRes$sumIntenY
            }
            sampleNames <- c(sampleNames, currRes$sampleName)
        }
    }
    
    nfile <- length(sampleNames)
    if(nfile <= 1) {
        stop("Cannot successfully genotype with less than two samples")
    }
    logfn("finished preprocessing %i samples", nfile)
    
    if(!is.null(isMale) && length(isMale) != nfile) {
        stop("the length of isMale must match the sample count")
    }
    
    # determines which arrays are male and which are female if we need to
    if(genderInferenceRequired) {
        logfn("inferring gender for each sample based on probe intensities")
        
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
        
        logfn("%i samples inferred as male:", sum(isMale))
        for(n in sampleNames[isMale]) {
            logfn("    %s", n)
        }
        logfn("%i samples inferred as female:", sum(!isMale))
        for(n in sampleNames[!isMale]) {
            logfn("    %s", n)
        }
    }
    
    results <- list()
    outFileConnections <- NULL
    argLists <- list()
    for (chri in chromosomes) {
        chrIndices <- which(snpInfo$chrId == chri)
        chrSnpInfo <- snpInfo[chrIndices, , drop=FALSE]
        
        chrLogSNP <- NULL
        if(!is.null(chrSnpInfo$logSNP)) {
            chrLogSNP <- chrSnpInfo$logSNP
            names(chrLogSNP) <- chrSnpInfo$snpId
        }
        
        for(chunkIndex in 1 : length(chrChunks[[chri]])) {
            chunk <- chrChunks[[chri]][[chunkIndex]]
            logfn(
                "geno/vinotyping chromosome %s from probeset #%i to probeset #%i",
                chri,
                chunk$start,
                chunk$end)
            
            #startTime <- getTime()
        
            # paste the samples together for genotyping
            chunkMatrix <- NULL
            for (i in 1:nfile) {
                chunkFile <- .chunkFileName(cacheDir, "snp", sampleNames[i], chri, probesetChunkSize, chunkIndex)
                load(chunkFile)
                chunkMatrix <- cbind(chunkMatrix, chunkData)
                rm(chunkData)
            }
            
            chunkRange <- chunk$start : chunk$end
            
            #cat("time it took us to get to genotypethis\n")
            #timeReport(startTime)
            
            #startTime <- getTime()
            if(numCores >= 2) {
                # the arg list is used to accumulate arguments until we're
                # ready to execute them in parallel
                argLists[[length(argLists) + 1]] <- list(
                    chr = chri,
                    snpContAvg = chunkMatrix,
                    hint = chrSnpInfo$snpHetHint[chunkRange],
                    isInPAR = chrSnpInfo$isInPAR[chunkRange],
                    isMale = isMale,
                    confScoreThreshold = confScoreThreshold,
                    logOn = logOn,
                    logSnp = chrLogSNP[chunkRange])
                
                if(length(argLists) >= numCores ||
                   (chri == chromosomes[length(chromosomes)] &&
                   chunkIndex == length(chrChunks[[chri]]))) {

                    # parallel apply using snow then reset the arg list
                    chunkResultsList <- mclapply(argLists, .applyGenotypeAnyChrChunk, mc.cores=numCores)

                    for(i in 1 : length(chunkResultsList)) {
                        for(logLine in chunkResultsList[[i]]$logLines) {
                            logfn(logLine)
                        }
                        
                        chunkResult <- chunkResultsList[[i]]$genoResult
                        if(is.null(outputDir)) {
                            if(length(results) == 0) {
                                results <- chunkResult
                            } else {
                                results <- mapply(rbind, results, chunkResult, SIMPLIFY = FALSE)
                            }
                        } else {
                            for(k in 1 : length(chunkResult)) {
                                colnames(chunkResult[[k]]) <- sampleNames
                                rownames(chunkResult[[k]]) <- rownames(argLists[[i]][["snpContAvg"]])
                            }
                            
                            if(is.null(outFileConnections)) {
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
            } else {
                chunkResult <- .genotypeAnyChrChunk(
                    chr = chri,
                    snpContAvg = chunkMatrix,
                    hint = chrSnpInfo$snpHetHint[chunkRange],
                    isInPAR = chrSnpInfo$isInPAR[chunkRange],
                    isMale = isMale,
                    confScoreThreshold = confScoreThreshold,
                    logSnp = chrLogSNP[chunkRange],
                    logfn = logfn)
                if(is.null(outputDir)) {
                    # we only accumulate results if there is no function
                    if(length(results) == 0) {
                        results <- chunkResult
                    } else {
                        results <- mapply(rbind, results, chunkResult, SIMPLIFY = FALSE)
                    }
                } else {
                    chunkProbesetInfo <- snpInfo[chrIndices[chunkRange], ]
                    
                    for(k in 1 : length(chunkResult)) {
                        colnames(chunkResult[[k]]) <- sampleNames
                        rownames(chunkResult[[k]]) <- chunkProbesetInfo$snpId
                    }
                    
                    if(is.null(outFileConnections)) {
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
    }
    rm(argLists)
    
    # if the user asked us to use a function to process the results that means
    # we have not accumulated any results to return, otherwise we have a bit
    # of post-processing to do to make sure that the matrices that we return
    # to the user match up nicely with what they passed in
    if(!is.null(outputDir)) {
        results <- NULL
        for(con in outFileConnections) {
            flush(con)
            close(con)
        }
    } else {
        # using list/unlist because c() will coerce factors into an integer
        snpIdsByChr <- list()
        for(chri in chromosomes) {
            snpIdsByChr[[chri]] <- snpInfo$snpId[snpInfo$chrId == chri]
        }
        snpIdsByChr <- unlist(snpIdsByChr, use.names = FALSE)
        
        snpIdsInOrder <- snpInfo$snpId[snpInfo$chrId %in% chromosomes]
        snpOrdering <- match(snpIdsInOrder, snpIdsByChr)
        if(any(is.na(snpOrdering))) {
            stop("internal error: failed to match up SNP IDs in results")
        }
        
        # reorder the SNPs so that they match up with the snpInfo frame that
        # was passed in
        for(i in 1 : length(results)) {
            results[[i]] <- results[[i]][snpOrdering, ]
            rownames(results[[i]]) <- snpIdsInOrder
            colnames(results[[i]]) <- sampleNames
        }
        
        results <- makeProbesetSampleMatrixes(results, snpIdsInOrder, sampleNames)
    }
    
    # clean-up unless we were asked to retain the cache
    if(!retainCache) {
        if(logOn) {
            logfn("cleaning up cache files before returning geno/vinotype results")
        }
        
        for(sampleName in sampleNames) {
            for(currChr in allChr) {
                for(chunkIndex in 1 : length(chrChunks[[currChr]])) {
                    chunkFile <- .chunkFileName(cacheDir, "snp", sampleName, currChr, probesetChunkSize, chunkIndex)
                    if(file.exists(chunkFile)) {
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
.applyGenotypeAnyChrChunk <- function(argList) {
    # To explain the strange logging code below, its purpose is to define a
    # logging function which appends to a list rather than printing to a
    # connection. The reason that we can't just directly print to a connection
    # is that this code could be running on a different machine than the main
    # process (the snow package is used for this), so rather than write
    # to a file we will return the log results along with the genotyping results
    logEnv <- new.env(parent=emptyenv())
    logEnv[["logLines"]] <- list()
    if(argList$logOn) {
        logfn <- function(fmt, ...) {
            logfnTry <- function() {
                currLineCount <- length(logEnv[["logLines"]])
                logEnv[["logLines"]][[currLineCount + 1]] <- sprintf(fmt, ...)
            }
            onError <- function(e) {
                warning("failed to log message")
            }
            tryCatch(logfnTry(), error=onError)
        }
    } else {
        logfn <- NULL
    }

    genoResult <- .genotypeAnyChrChunk(
        argList$chr,
        argList$snpContAvg,
        argList$hint,
        argList$isInPAR,
        argList$isMale,
        argList$confScoreThreshold,
        argList$logSnp,
        logfn)

    list(genoResult=genoResult, logLines=logEnv[["logLines"]])
}

inferGender <- function(
        meanIntensityXPerArray,
        meanIntensityYPerArray,
        meanIntensityPerAutosome) {

    genderClust <- kmeans(cbind(meanIntensityYPerArray, meanIntensityXPerArray), 2)$cluster
    yMean1 <- mean(meanIntensityYPerArray[genderClust == 1])
    yMean2 <- mean(meanIntensityYPerArray[genderClust == 2])
    
    # we expect that the cluster with a more intense Y chromosome will be males
    if(yMean1 > yMean2) {
        maleClust <- 1
        yMeanMale <- yMean1
        yMeanFemale <- yMean2
    } else {
        maleClust <- 2
        yMeanMale <- yMean2
        yMeanFemale <- yMean1
    }
    
    if(yMeanFemale > min(meanIntensityPerAutosome)) {
        # if the y intensity of the samples grouped as female is greater than
        # the smallest autosome mean intensity we infer that there are no
        # female samples
        isMale <- rep(TRUE, length(meanIntensityXPerArray))
    } else {
        isMale <- genderClust == maleClust
    }
    
    isMale
}

inferGenderFromSnpContAvg <- function(snpContAvg, snpInfo) {

    if(!inherits(snpInfo, "data.frame") || !all(c("snpId", "chrId") %in% names(snpInfo))) {
        stop("You must supply a \"snpInfo\" data frame parameter which has ",
            "at a minimum the \"snpId\" and \"chrId\" ",
            "components. Please see the help documentation for more details.")
    }
    
    if(!inherits(snpContAvg, "snpContAvg")) {
        stop("the snpContAvg argument should inherit from the snpContAvg class")
    }
    
    allChr <- as.character(unique(snpInfo$chrId))
    allAutosomes <- setdiff(allChr, c("X", "Y", "M"))
    if(length(allAutosomes) >= 1 && ("X" %in% allChr) && ("Y" %in% allChr)) {
        getChrSnpAvgs <- function(chrName) {
            chrSnpInfo <- snpInfo[snpInfo$chrId == chrName, ]
            chrMatch <- match(chrSnpInfo$snpId, rownames(snpContAvg))
            chrMatch <- chrMatch[!is.na(chrMatch)]
            if(length(chrMatch) == 0) {
                stop("failed to find any data matches for chromosome ", chrName)
            }
            
            snpContAvg$snpAverages[chrMatch, , drop=FALSE]
        }
        
        meanChrIntensityPerArray <- function(chrName) {
            chrSnpAvgs <- getChrSnpAvgs(chrName)
            meanPerArr <- apply(chrSnpAvgs, 2, sum) / nrow(chrSnpAvgs)
            names(meanPerArr) <- colnames(snpContAvg)
            meanPerArr
        }
        
        meanIntensityPerAutosome <- NULL
        for(currChr in allAutosomes) {
            chrSnpAvgs <- getChrSnpAvgs(currChr)
            meanIntensityPerAutosome[currChr] <- sum(chrSnpAvgs) / length(chrSnpAvgs)
        }
        
        inferGender(
            meanIntensityXPerArray = meanChrIntensityPerArray("X"),
            meanIntensityYPerArray = meanChrIntensityPerArray("Y"),
            meanIntensityPerAutosome = meanIntensityPerAutosome)
    } else {
        # if we don't have X, Y and at least one autosome we can't infer gender
        NULL
    }
}

genotypeSnps <- function(
        snpContAvg,
        snpInfo,
        isMale = NULL,
        chromosomes = c(1:19, "X", "Y", "M"),
        confScoreThreshold = 1e-05,
        logFile = NULL) {

    if(!inherits(snpInfo, "data.frame") || !all(c("snpId", "chrId") %in% names(snpInfo))) {
        stop("You must supply a \"snpInfo\" data frame parameter which has ",
             "at a minimum the \"snpId\" and \"chrId\" ",
             "components. Please see the help documentation for more details.")
    }
    snpInfo$snpId <- as.factor(snpInfo$snpId)

    if(inherits(logFile, "character")) {
        logFile <- file(logFile, "wt")
    } else if(!is.null(logFile) && !inherits(logFile, "connection")) {
        stop("the logFile parameter should be either NULL, a file name, ",
             "or a connection")
    }

    logfn <- NULL
    if(!is.null(logFile)) {
        logfn <- function(fmt, ...) {
            logfnTry <- function() {
                cat(sprintf(fmt, ...), file=logFile)
                cat("\n", file=logFile)
                flush(logFile)
            }
            onError <- function(e) {
                warning("failed to log message")
            }
            tryCatch(logfnTry(), error=onError)
        }
    }
    
    if(is.null(isMale)) {
        isMale <- inferGenderFromSnpContAvg(snpContAvg, snpInfo)
    }
    
    results <- NULL
    allChr <- intersect(chromosomes, as.character(unique(snpInfo$chrId)))
    for(currChr in allChr) {
        if(!is.null(logfn)) {
            logfn("genotyping and vinotyping chromosome %s", currChr)
        }
        chrSnpInfo <- snpInfo[snpInfo$chrId == currChr, , drop=FALSE]
        chrMatch <- match(chrSnpInfo$snpId, rownames(snpContAvg))
        chrSnpInfo <- chrSnpInfo[!is.na(chrMatch), , drop=FALSE]
        chrMatch <- chrMatch[!is.na(chrMatch)]
        chrData <- snpContAvg[chrMatch, , drop=FALSE]
        chrResults <- .genotypeAnyChrChunk(
            chr = currChr,
            snpContAvg = chrData,
            hint = chrSnpInfo$snpHetHint,
            isInPAR = chrSnpInfo$isInPAR,
            isMale = isMale,
            confScoreThreshold = confScoreThreshold,
            logSnp = chrSnpInfo$logSNP,
            logfn)
        chrResults <- makeProbesetSampleMatrixes(chrResults, chrSnpInfo$snpId, colnames(snpContAvg))
        results <- rbind(results, chrResults)
    }
    
    results$isMale <- isMale
    
    outOfOrderSnpIds <- rownames(results)
    reorderIndexes <- match(snpInfo$snpId, outOfOrderSnpIds)
    reorderIndexes <- reorderIndexes[!is.na(reorderIndexes)]
    results[reorderIndexes, , drop=FALSE]
}

.genotypeAnyChrChunk <- function(
        chr,
        snpContAvg,
        hint = NULL,
        isInPAR = NULL,
        isMale = NULL,
        confScoreThreshold = 1e-05,
        logSnp = NULL,
        logfn = NULL) {

    if(chr == "X") {
        .genotypeXChromosomeChunk(snpContAvg, hint, isInPAR, isMale, confScoreThreshold, logSnp, logfn)
    } else if(chr == "Y") {
        .genotypeYChromosomeChunk(snpContAvg, hint, isMale, confScoreThreshold, logSnp, logfn)
    } else if(chr == "M") {
        .genotypeHomozygousChunk(snpContAvg, confScoreThreshold, logSnp, logfn, "M")
    } else {
        # chr is an autosome
        .genotypeChunk(snpContAvg, hint, confScoreThreshold, logSnp, logfn, "autosome")
    }
}

# for genotyping probeset "chunks" (where ms and ss rows map to probesets and
# columns map to arrays (CEL files). This function is useful for genotyping
# autosomes and other sections of the genome where it is possible to have 3
# genotyping outcomes. Otherwise use genotypeHomozygousChunk
.genotypeChunk <- function(snpContAvg, hint, confScoreThreshold, logSnp=NULL, logfn=NULL, kind = c("autosome", "PAR", "femaleX")) {
    
    kind <- match.arg(kind)
    if(!is.null(logfn)) {
        kindDescrip <- switch(kind, autosome="autosome", PAR="pseudoautosomal", "female X chromosome")
        logfn(
            "performing standard genotyping for %i %s SNPs and %i samples",
            nrow(snpContAvg),
            kindDescrip,
            ncol(snpContAvg))
    }
    
    numArrays <- ncol(snpContAvg)
    numProbesets <- nrow(snpContAvg)
    
    if(length(hint) == 0) {
        hint <- rep(0, numProbesets)
    }
    
    if(is.null(logfn)) {
        logSnp <- NULL
    }
    
    snpNames <- rownames(snpContAvg)
    ms <- snpContAvg$snpContrasts
    ss <- snpContAvg$snpAverages    
    
    # geno/vinotype each of the probesets in this chunk
    transformMethod <- attr(snpContAvg, "transformMethod", exact=TRUE)
    geno <- matrix(0, nrow = numProbesets, ncol = numArrays)
    vino <- matrix(0, nrow = numProbesets, ncol = numArrays)
    conf <- matrix(0, nrow = numProbesets, ncol = numArrays)
    for(probesetIndex in 1 : numProbesets) {
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
            transformMethod,
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
    
    if(confScoreThreshold > 0) {
        chunkResult$geno[chunkResult$conf < confScoreThreshold] <- -1
    }
    
    # TODO restore rownames and colnames to chunk result
    chunkResult
}

# for genotyping probeset "chunks" (where ms and ss rows map to probesets and
# columns map to arrays (CEL files). This function is useful for genotyping
# the parts of the genome where you do not expect to observe heterozygous
# alleles
.genotypeHomozygousChunk <- function(snpContAvg, confScoreThreshold, logSnp=NULL, logfn=NULL, kind = c("M", "maleX", "Y")) {
    kind <- match.arg(kind)
    if(!is.null(logfn)) {
        kindDescrip <- switch(kind, M="mitochondrial", maleX="male X chromosome", Y="Y chromosome")
        logfn(
            "performing homozygous genotyping for %i %s SNPs and %i samples",
            nrow(snpContAvg),
            kindDescrip,
            ncol(snpContAvg))
    }
    
    numArrays <- ncol(snpContAvg)
    numProbesets <- nrow(snpContAvg)
    
    snpNames <- rownames(snpContAvg)
    ms <- snpContAvg$snpContrasts
    ss <- snpContAvg$snpAverages
    
    # geno/vinotype each of the probesets in this chunk
    transformMethod <- attr(snpContAvg, "transformMethod", exact=TRUE)
    geno <- matrix(0, nrow = numProbesets, ncol = numArrays)
    vino <- matrix(0, nrow = numProbesets, ncol = numArrays)
    conf <- matrix(0, nrow = numProbesets, ncol = numArrays)
    for(probesetIndex in 1 : numProbesets) {
        if(!is.null(logSnp) && logSnp[probesetIndex]) {
            if(is.null(snpNames)) {
                logfn("==== GENOTYPING AND VINOTYPING SNP ====")
            } else {
                logfn("==== GENOTYPING AND VINOTYPING %s ====", snpNames[probesetIndex])
            }
        }
        
        currVals <- genotypeHomozygous(
            ms[probesetIndex, ],
            ss[probesetIndex, ],
            transformMethod,
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
    
    if(confScoreThreshold > 0) {
        chunkResult$geno[chunkResult$conf < confScoreThreshold] <- -1
    }
    
    # TODO restore rownames to chunk result
    chunkResult
}

.genotypeXChromosomeChunk <- function(snpContAvg, hint, isInPAR, isMale, confScoreThreshold, logSnp=NULL, logfn=NULL) {
    parIndices <- which(as.logical(isInPAR))
    numArrays <- ncol(snpContAvg)
    numProbesets <- nrow(snpContAvg)
    
    if(length(hint) == 0) {
        hint <- rep(0, numProbesets)
    }
    
    # we'll call the X probesets that aren't in the PAR "normal"
    normalIndices <- 1 : numProbesets
    if(length(parIndices) >= 1) {
        normalIndices <- setdiff(normalIndices, parIndices)
    }
    
    maleColumns <- which(isMale)
    femaleColumns <- which(!isMale)
    
    initializeNamesIfMissing <- function(matrixList, itemNames) {
        for(itemName in itemNames) {
            if(!(itemName %in% names(matrixList))) {
                matrixList[[itemName]] <- matrix(nrow = numProbesets, ncol = numArrays)
            }
        }
        
        matrixList
    }
    results <- list()
    
    normalHint <- hint[normalIndices]
    normalLogSnp <- logSnp[normalIndices]
    
    # we can treat the females like autosomes for the purposes of genotyping
    if(length(femaleColumns) >= 2) {
        # we are going to treat the female probesets just like autosomes for
        # the purposes of geno/vinotyping
        normalFemaleResult <- .genotypeChunk(
            snpContAvg[normalIndices, femaleColumns, drop = FALSE],
            normalHint,
            confScoreThreshold,
            normalLogSnp,
            logfn,
            "femaleX")
        results <- initializeNamesIfMissing(results, names(normalFemaleResult))
        for(itemName in names(normalFemaleResult)) {
            results[[itemName]][normalIndices, femaleColumns] <- normalFemaleResult[[itemName]]
        }
    }
    
    # the non-PAR male probesets will get genotyped as homozygous
    if(length(maleColumns) >= 2) {
        normalMaleResult <- .genotypeHomozygousChunk(
            snpContAvg[normalIndices, maleColumns, drop = FALSE],
            confScoreThreshold,
            normalLogSnp,
            logfn,
            "maleX")
        results <- initializeNamesIfMissing(results, names(normalMaleResult))
        for(itemName in names(normalMaleResult)) {
            results[[itemName]][normalIndices, maleColumns] <- normalMaleResult[[itemName]]
        }
    }
    
    # the PAR probesets will be treated as autosomes whether they're male or female
    if(length(parIndices) >= 1) {
        parHint <- hint[parIndices]
        parLogSnp <- logSnp[parIndices]
        
        parResult <- .genotypeChunk(
            snpContAvg[parIndices, , drop = FALSE],
            parHint,
            confScoreThreshold,
            parLogSnp,
            logfn,
            "PAR")
        results <- initializeNamesIfMissing(results, names(parResult))
        for(itemName in names(parResult)) {
            results[[itemName]][parIndices, ] <- parResult[[itemName]]
        }
    }
    
    results
}

.genotypeYChromosomeChunk <- function(snpContAvg, hint, isMale, confScoreThreshold, logSnp=NULL, logfn=NULL) {
    numArrays <- ncol(snpContAvg)
    numProbesets <- nrow(snpContAvg)
    
    maleColumns <- which(isMale)
    
    maleResult <- .genotypeHomozygousChunk(
        snpContAvg[, maleColumns, drop = FALSE],
        confScoreThreshold,
        logSnp,
        logfn,
        "Y")
    
    # fill in the female columns using NA so that the dimensions of the
    # input data match the dimensions of the output data
    fillInNAFemaleCols <- function(maleMatrix) {
        allMatrix <- matrix(nrow = numProbesets, ncol = numArrays)
        allMatrix[, maleColumns] <- maleMatrix
        allMatrix
    }
    lapply(maleResult, fillInNAFemaleCols)
}
