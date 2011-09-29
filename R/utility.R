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
    }
    
    connections
}

.writeResultsToFlatFile <- function(connections, matrices, sep = "\t") {
    n <- length(connections)
    if(length(matrices) != n) {
        stop("number of file connections must match number of matrices")
    }
    
    for(i in 1 : n) {
        write.table(
                matrices[[i]],
                file = connections[[i]],
                sep = sep,
                row.names = TRUE,
                col.names = FALSE,
                qmethod = "double")
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
expandCelFiles <- function(filename) {
    retVal <- NULL
    if(!file.exists(filename)) {
        stop("failed to find file named \"", filename, "\"")
    } else if(file.info(filename)$isdir) {
        # if it's a dir expand it to all CEL file contents
        fileListing <- dir(filename, full.names = TRUE)
        retVal <- fileListing[grep("*.CEL$", toupper(fileListing))]
    } else if(grep("*.CEL$", toupper(filename))) {
        retVal <- filename
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

ccsTransform <- function(abMat, k = 4) {
    a <- abMat[ , 1]
    twoPowA <- 2 ^ a
    b <- abMat[ , 2]
    twoPowB <- 2 ^ b
    
    m <- asinh(k * (twoPowA - twoPowB) / (twoPowA + twoPowB)) / asinh(k)
    s <- (a + b) / 2
    
    contAvgMat <- matrix(c(m, s), ncol = 2)
    rownames(contAvgMat) <- rownames(abMat)
    colnames(contAvgMat) <- c("intensityConts", "intensityAvgs")
    
    contAvgMat
}

maTransform <- function(abMat) {
    a <- abMat[ , 1]
    b <- abMat[ , 2]
    
    m <- a - b
    s <- (a + b) / 2
    
    contAvgMat <- matrix(c(m, s), ncol = 2)
    rownames(contAvgMat) <- rownames(abMat)
    colnames(contAvgMat) <- c("intensityConts", "intensityAvgs")
    
    contAvgMat
}
