.initFlatFileConnections <- function(
        dir,
        matrices,
        sep = "\t",
        prefix = "mouseDivResults_",
        suffix = ".txt") {

    connections <- list()
    for(name in names(matrices)) {
        outFileName <- file.path(dir, paste(prefix, name, suffix, sep = ""))
        con <- file(description = outFileName, open = "wt")
        connections[[name]] <- con
        
        # write the header
        write.table(
                matrix(c("ID", colnames(matrices[[name]])), nrow = 1),
                file = con,
                sep = sep,
                row.names = FALSE,
                col.names = FALSE,
                qmethod = "double")
        flush(con)
    }
    
    connections
}

.writeResultsToFlatFile <- function(connections, matrices, sep = "\t") {
    n <- length(connections)
    if(length(matrices) != n) {
        stop("number of file connections must match number of matrices")
    }
    
    for(i in seq_len(n)) {
        write.table(
                matrices[[i]],
                file = connections[[i]],
                sep = sep,
                row.names = TRUE,
                col.names = FALSE,
                qmethod = "double")
        flush(connections[[i]])
    }
}

#### biweight midcovariance from Wilcox(1997)
#### code from Dr.RichHerrington at
#### http://www.unt.edu/benchmarks/archives/2001/december01/rss.htm
.bicov <- function(x, y) {
    mx <- median(x)
    my <- median(y)
    ux <- abs((x - mx)/(9 * qnorm(0.75) * mad(x)))
    uy <- abs((y - my)/(9 * qnorm(0.75) * mad(y)))
    aval <- ifelse(ux <= 1, 1, 0)
    bval <- ifelse(uy <= 1, 1, 0)
    top <- sum(aval * (x - mx) * (1 - ux^2)^2 * bval * (y - my) * (1 - uy^2)^2)
    top <- length(x) * top
    botx <- sum(aval * (1 - ux^2) * (1 - 5 * ux^2))
    boty <- sum(bval * (1 - uy^2) * (1 - 5 * uy^2))
    bi <- top/(botx * boty)
    bi
}

#### Calculate biweight midcorrelation
.bivar <- function(x) {
    mx <- median(x)
    ux <- abs((x - mx)/(9 * qnorm(0.75) * mad(x)))
    aval <- ifelse(ux <= 1, 1, 0)
    top <- sum((aval * (x - mx) * (1 - ux^2)^2)^2)
    top <- length(x) * top
    botx <- (sum(aval * (1 - ux^2) * (1 - 5 * ux^2)))^2
    bi <- top/botx
    bi
}

# a simple function that will return all CEL file in a dir if given a
# dir argument, or a single CEL file if that is what it's given. If the given
# file is not null and does not end in ".CEL" then NULL is returned
.expandCelFiles <- function(filenames) {
    retVal <- NULL
    for(f in filenames) {
        if(!file.exists(f)) {
            stop("failed to find file named \"", f, "\"")
        } else if(file.info(f)$isdir) {
            # if it's a dir expand it to all CEL file contents
            fileListing <- dir(f, full.names=TRUE)
            retVal <- c(retVal, fileListing[grep("*.CEL$", toupper(fileListing))])
        } else {
            retVal <- c(retVal, f)
        }
    }
    retVal
}

.chunkIndices <- function(to, by) {
    chunks <- list()
    chunkNumber <- 0
    for(chunkStart in seq(from = 1, to = to, by = by)) {
        chunkEnd <- chunkStart + by - 1
        if(chunkEnd > to) {
            chunkEnd <- to
        }
        
        chunkNumber <- chunkNumber + 1
        chunks[[chunkNumber]] <- list(start=chunkStart, end=chunkEnd)
    }
    
    chunks
}

# TODO we need more parameters to make sure this is really meaningfully unique
.chunkFileName <- function(baseDir, kind, sampleName, chrName, chunkSize, chunkIndex) {
    chunkFile <- paste(sampleName, "-", kind, "-", chrName, "-", chunkSize, "-", chunkIndex, ".RData", sep="")
    
    file.path(baseDir, chunkFile)
}

.fileBaseWithoutExtension <- function(celFileName) {
    sub("\\.[^\\.]*$", "", basename(celFileName))
}

makeProbesetSampleMatrixes <- function(mats, probesetIds, sampleNames) {
    matNames <- names(mats)
    mats[["probesetIds"]] <- as.character(probesetIds)
    mats[["sampleNames"]] <- as.character(sampleNames)
    class(mats) <- c("probesetSampleMatrixes", class(mats))
    attr(mats, "matrixNames") <- matNames
    
    mats
}

makeSnpIntensities <- function(aIntensities, bIntensities, snpIds, sampleNames) {
    if(!is.matrix(aIntensities)) {
        aIntensities <- matrix(aIntensities, ncol=1)
    }
    if(!is.matrix(bIntensities)) {
        bIntensities <- matrix(bIntensities, ncol=1)
    }
    
    retVal <- list(
        A=aIntensities,
        B=bIntensities,
        probesetIds=as.character(snpIds),
        sampleNames=as.character(sampleNames))
    class(retVal) <- c("snpIntensities", "probesetSampleMatrixes", class(retVal))
    attr(retVal, "matrixNames") <- c("A", "B")
    retVal
}

makeSnpContAvg <- function(snpContrasts, snpAverages, snpIds, sampleNames, transformMethod=NULL) {
    if(!is.matrix(snpContrasts)) {
        snpContrasts <- matrix(snpContrasts, ncol=1)
    }
    if(!is.matrix(snpAverages)) {
        snpAverages <- matrix(snpAverages, ncol=1)
    }
    
    retVal <- list(
        snpContrasts=snpContrasts,
        snpAverages=snpAverages,
        probesetIds=as.character(snpIds),
        sampleNames=as.character(sampleNames))
    class(retVal) <- c("snpContAvg", "probesetSampleMatrixes", class(retVal))
    attr(retVal, "matrixNames") <- c("snpContrasts", "snpAverages")
    if(!is.null(transformMethod)) {
        attr(retVal, "transformMethod") <- transformMethod
    }
    retVal
}

dim.probesetSampleMatrixes <- function(x) {
    c(length(x$probesetIds), length(x$sampleNames))
}

dimnames.probesetSampleMatrixes <- function(x) {
    list(x$probesetIds, x$sampleNames)
}

cbind.probesetSampleMatrixes <- function(x, ...) {
    if(is.null(x)) {
        if(length(list(...))) {
            cbind.probesetSampleMatrixes(...)
        } else {
            x
        }
    } else {
        probesetIds <- x$probesetIds
        numProbesets <- length(probesetIds)
        
        for(curr in list(...)) {
            if(!is.null(curr)) {
                if(!inherits(curr, "probesetSampleMatrixes")) {
                    stop("cbind.probesetSampleMatrixes requires arguments to inherit from the probesetSampleMatrixes class")
                }
                
                if(length(curr$probesetIds) != numProbesets || any(curr$probesetIds != probesetIds)) {
                    stop("in cbind.probesetSampleMatrixes probesetIds do not match up")
                }
                
                for(mName in attr(x, "matrixNames", exact=TRUE)) {
                    x[[mName]] <- cbind(x[[mName]], curr[[mName]])
                }
                
                x$sampleNames <- c(x$sampleNames, curr$sampleNames)
            }
        }
        
        x
    }
}

rbind.probesetSampleMatrixes <- function(x, ...) {
    if(is.null(x)) {
        if(length(list(...))) {
            rbind.probesetSampleMatrixes(...)
        } else {
            x
        }
    } else {
        sampleNames <- x$sampleNames
        numSamples <- length(sampleNames)
        
        for(curr in list(...)) {
            if(!is.null(curr)) {
                if(!inherits(curr, "probesetSampleMatrixes")) {
                    stop("rbind.probesetSampleMatrixes requires arguments to inherit from the probesetSampleMatrixes class")
                }
                
                if(length(curr$sampleNames) != numSamples || any(curr$sampleNames != sampleNames)) {
                    stop("in rbind.probesetSampleMatrixes sampleNames do not match up")
                }
                
                for(mName in attr(x, "matrixNames", exact=TRUE)) {
                    x[[mName]] <- rbind(x[[mName]], curr[[mName]])
                }
                
                x$probesetIds <- c(x$probesetIds, curr$probesetIds)
            }
        }
        
        x
    }
}

`[.probesetSampleMatrixes` <- function(x, i, j, drop=FALSE) {
    if(drop) {
        stop("drop=TRUE is not supported for probesetSampleMatrixes")
    }
    
    hasI <- !missing(i)
    hasJ <- !missing(j)
    
    # convert i and j to integer indexes if necessary
    toIntIndex <- function(index, isRow) {
        if(is.numeric(index)) {
            index
        } else if(is.logical(index)) {
            which(index)
        } else if(is.character(index) || is.factor(index)) {
            match(index, if(isRow) x$probesetIds else x$sampleNames)
        } else {
            stop(paste("bad index type:", class(index), collapse=", "))
        }
    }
    
    if(hasI) {
        i <- toIntIndex(i, TRUE)
        x$probesetIds <- x$probesetIds[i]
    }
    if(hasJ) {
        j <- toIntIndex(j, FALSE)
        x$sampleNames <- x$sampleNames[j]
    }
    
    if(hasI && hasJ) {
        for(mName in attr(x, "matrixNames", exact=TRUE)) {
            x[[mName]] <- x[[mName]][i, j, drop=FALSE]
        }
    } else if(hasI) {
        for(mName in attr(x, "matrixNames", exact=TRUE)) {
            x[[mName]] <- x[[mName]][i, , drop=FALSE]
        }
    } else if(hasJ) {
        for(mName in attr(x, "matrixNames", exact=TRUE)) {
            x[[mName]] <- x[[mName]][, j, drop=FALSE]
        }
    }
    
    x
}

#keepSnpsInChr <- function(mat, chrs, info) {
#    mat[rownames(mat) %in% info$snpId[info$chrId %in% chrs], , drop=FALSE]
#}
#
#dropSnpsInChr <- function(mat, chrs, info) {
#    mat[rownames(mat) %in% info$snpId[!(info$chrId %in% chrs)], , drop=FALSE]
#}

ccsTransform <- function(ab, k = 4) {
    UseMethod("ccsTransform")
}

ccsTransform.matrix <- function(ab, k = 4) {
    a <- ab[ , 1]
    twoPowA <- 2 ^ a
    b <- ab[ , 2]
    twoPowB <- 2 ^ b
    
    m <- asinh(k * (twoPowA - twoPowB) / (twoPowA + twoPowB)) / asinh(k)
    s <- (a + b) / 2
    
    contAvgMat <- matrix(c(m, s), ncol = 2)
    rownames(contAvgMat) <- rownames(ab)
    colnames(contAvgMat) <- c("intensityConts", "intensityAvgs")
    
    contAvgMat
}

ccsTransform.snpIntensities <- function(ab, k = 4) {
    a <- ab$A
    twoPowA <- 2 ^ a
    b <- ab$B
    twoPowB <- 2 ^ b
    
    m <- asinh(k * (twoPowA - twoPowB) / (twoPowA + twoPowB)) / asinh(k)
    s <- (a + b) / 2
    
    makeSnpContAvg(m, s, rownames(ab), colnames(ab), "CCS")
}

maTransform <- function(ab) {
    UseMethod("maTransform")
}

maTransform.matrix <- function(ab) {
    a <- ab[ , 1]
    b <- ab[ , 2]
    
    m <- a - b
    s <- (a + b) / 2
    
    contAvgMat <- matrix(c(m, s), ncol = 2)
    rownames(contAvgMat) <- rownames(ab)
    colnames(contAvgMat) <- c("intensityConts", "intensityAvgs")
    
    contAvgMat
}

maTransform.snpIntensities <- function(ab) {
    a <- ab$A
    b <- ab$B
    
    m <- a - b
    s <- (a + b) / 2
    
    makeSnpContAvg(m, s, rownames(ab), colnames(ab), "MA")
}

.defaultNumCores <-
    if(require("multicore", character.only=T)) {
        multicore:::detectCores()
    } else {
        1
    }

.mylapply <- function(X, FUN, ..., numCores) {
    if(is.null(numCores)) {
        numCores <- .defaultNumCores
    }
    
    if(numCores >= 2 && length(X) >= 2 && require("multicore", character.only=T)) {
        mclapply(X, FUN, ..., mc.cores=numCores)
    } else {
        lapply(X, FUN, ...)
    }
}

# a value is considered a "true" list if it is a list that is not also a
# data frame
.isTrueList <- function(x) {
    is.list(x) && !is.data.frame(x)
}

# a little function to make sure that we turn any item which is not a list into
# a list (which makes some algorithms more general/consistent)
.listify <- function(x) {
    if(is.na(x) || is.null(x)) {
        x <- list()
    } else if(!.isTrueList(x)) {
        x <- list(x)
    }
    
    x
}
