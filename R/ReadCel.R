#########################################################################
# ReadCel.R
#
# Part of the MouseDivGenotype package
#
#########################################################################

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
