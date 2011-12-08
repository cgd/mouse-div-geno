#########################################################################
# ReadCel.R
#
# Part of the MouseDivGenotype package
#
#########################################################################

readCELFiles <- function(
        celFiles,
        snpProbeInfo,
        referenceDistribution = NULL,
        logFile = NULL,
        numCores = NULL) {

    if(!is.null(logFile) && !inherits(logFile, "connection")) {
        stop("the logFile parameter should be a connection")
    }
    
    snpId <- as.factor(snpProbeInfo$snpId)
    
    readCELFile <- function(celFile) {
        if(!is.null(logFile)) {
            cat("Reading and normalizing CEL file: ", celFile, "\n", file=logFile, sep="")
            flush(logFile)
        }
        
        celData <- read.celfile(celFile, intensity.means.only = TRUE)
        y <- log2(as.matrix(celData[["INTENSITY"]][["MEAN"]][snpProbeInfo$probeIndex]))
        
        if (length(snpProbeInfo$correction) > 0) {
            # C+G and fragment length correction y
            y <- y + snpProbeInfo$correction
        }
        
        if (length(referenceDistribution) > 0) {
            y <- normalize.quantiles.use.target(y, target = referenceDistribution)
        }
        
        # separate A alleles from B alleles and summarize to the probeset level
        allAint <- y[snpProbeInfo$isAAllele, 1, drop = FALSE]
        allBint <- y[!snpProbeInfo$isAAllele, 1, drop = FALSE]
        allAint <- subColSummarizeMedian(
            matrix(allAint, ncol = 1),
            snpId[snpProbeInfo$isAAllele])
        allBint <- subColSummarizeMedian(
            matrix(allBint, ncol = 1),
            snpId[!snpProbeInfo$isAAllele])
        
        aSnpIds <- rownames(allAint)
        bSnpIds <- rownames(allBint)
        if(!all(aSnpIds == bSnpIds)) {
            stop("The SNP IDs for the A alleles should match up with the SNP IDs ",
                "for the B alleles but they do not.")
        }
        
        list(allAint, allBint)
    }
    
    celFiles <- .expandCelFiles(celFiles)
    celDataList <- .mylapply(celFiles, readCELFile, numCores=numCores)
    aIntMatrix <- NULL
    bIntMatrix <- NULL
    for(celData in celDataList) {
        aIntMatrix <- cbind(aIntMatrix, celData[[1]])
        bIntMatrix <- cbind(bIntMatrix, celData[[2]])
    }
    
    sampleIds <- .fileBaseWithoutExtension(celFiles)
    colnames(aIntMatrix) <- sampleIds
    colnames(bIntMatrix) <- sampleIds
    makeSnpIntensities(aIntMatrix, bIntMatrix, rownames(aIntMatrix), sampleIds)
}

readCELFilesInvariants <- function(
    celFiles,
    invariantProbeInfo,
    referenceDistribution = NULL,
    logFile = NULL,
    numCores = NULL) {
    
    if(!is.null(logFile) && !inherits(logFile, "connection")) {
        stop("the logFile parameter should be a connection")
    }
    
    trueList <- .isTrueList(invariantProbeInfo)
    invariantProbeInfo <- .listify(invariantProbeInfo)
    if(is.null(referenceDistribution)) {
        referenceDistribution <- replicate(length(invariantProbeInfo), NULL)
    } else {
        referenceDistribution <- .listify(referenceDistribution)
    }
    
    invariantProbeInfo <- .mylapply(
        invariantProbeInfo,
        function(x) {
            x$probesetId <- as.factor(x$probesetId)
            x
        },
        numCores=numCores)
    
    readCELFilesForInvariant <- function(celFile) {
        if(!is.null(logFile)) {
            cat("Reading and normalizing CEL file: ", celFile, "\n", file=logFile, sep="")
            flush(logFile)
        }
        
        celData <- read.celfile(celFile, intensity.means.only = TRUE)
        
        results <- list()
        for(i in seq_along(invariantProbeInfo)) {
            y <- log2(as.matrix(celData[["INTENSITY"]][["MEAN"]][invariantProbeInfo[[i]]$probeIndex]))
            
            if (length(invariantProbeInfo[[i]]$correction) > 0) {
                # C+G and fragment length correction y
                y <- y + invariantProbeInfo[[i]]$correction
            }
            
            if (length(referenceDistribution[[i]]) > 0) {
                y <- normalize.quantiles.use.target(y, target = referenceDistribution[[i]])
            }
            
            results[[i]] <- subColSummarizeMedian(matrix(y, ncol = 1), invariantProbeInfo[[i]]$probesetId)
        }
        results
    }
    
    celFiles <- .expandCelFiles(celFiles)
    celDataList <- .mylapply(celFiles, readCELFilesForInvariant, numCores=numCores)
    invIntMatrixList <- replicate(length(invariantProbeInfo), NULL)
    names(invIntMatrixList) <- names(invariantProbeInfo)
    for(celData in celDataList) {
        for(i in seq_along(celData)) {
            invIntMatrixList[[i]] <- cbind(invIntMatrixList[[i]], celData[[i]])
        }
    }
    for(i in seq_along(invIntMatrixList)) {
        colnames(invIntMatrixList[[i]]) <- .fileBaseWithoutExtension(celFiles)
    }
    
    if(trueList) {
        invIntMatrixList
    } else {
        invIntMatrixList[[1]]
    }
}

.readSnpIntensitiesFromCEL <- function(
        celFiles,
        snpProbeInfo,
        referenceDistribution = NULL,
        numCores = NULL,
        logFile = NULL) {

    if(is.null(numCores)) {
        numCores <- .defaultNumCores
    }
    
    readEmpty <- function(index) {
        endIndex <- min(index + numCores - 1, length(celFiles))
        currCELFiles <- as.list(celFiles[index : endIndex])
        readCELs <- .mylapply(currCELFiles, readCELFiles, snpProbeInfo, referenceDistribution, logFile, 1, numCores=numCores)
        
        readNext <- function(offset) {
            headFunc <- function() {
                readCELs[[offset + 1]]
            }
            
            tailFunc <-
                if(index + offset == length(celFiles)) {
                    NULL
                } else {
                    function() {
                        if(offset + 1 < length(readCELs)) {
                            readNext(offset + 1)
                        } else {
                            readEmpty(index + offset + 1)
                        }
                    }
                }
            
            list(head = headFunc, tail = tailFunc)
        }
        
        readNext(0)
    }
    
    if(length(celFiles)) {
        function() {readEmpty(1)}
    } else {
        NULL
    }
}

.readSnpIntensitiesFromTab <- function(tabFile, logFile = NULL) {

    if(!is.null(logFile) && !inherits(logFile, "connection")) {
        stop("the logFile parameter should be a connection")
    }
    
    logOn <- !is.null(logFile)
    logfn <- function(fmt, ...) {
        if(logOn) {
            cat(sprintf(fmt, ...), file=logFile)
            cat("\n", file=logFile)
            flush(logFile)
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
            
            leftoversToTake <- seq_len(snpCount)
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
            
            headFunc <- function() {
                makeSnpIntensities(currData[[3]], currData[[4]], currData[[2]], currData[[1]][1])
            }
            
            list(head = headFunc, tail = tailFunc)
        }
        
        readNext(currSampleData)
    }
    
    readFirst
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
