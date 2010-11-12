plotSNP <- function(nm, ns, geno, vino, main = "", xlab = "Contrast", ylab = "Average Intensity") {
    plot(nm, ns, xlab = xlab, ylab = ylab, main = main)
    points(nm[geno == 1], ns[geno == 1], col = 2)           # red
    points(nm[geno == 2], ns[geno == 2], col = 3)           # green
    points(nm[geno == 3], ns[geno == 3], col = 4)           # blue
    points(nm[vino == 1], ns[vino == 1], pch = 19, col = 6) # filled point
}

plotMouseDivArrayImage <- function(celFilename)
{
    my <- as.matrix(read.celfile(celFilename, intensity.means.only = TRUE)[["INTENSITY"]][["MEAN"]][1:6892960])
    image(
        1 : 2572,
        1 : 2680,
        matrix(log2(my), byrow = T, ncol = 2680),
        col = c("green", "black", "red"),
        xlab = "",
        ylab = "")
}

mouseDivDensityPlot <- function(celFilenames, snpProbeInfo, type = c("Average", "MatchedSet"))
{
    snpProbeInfo$snpId <- as.factor(snpProbeInfo$snpId)
    
    type <- match.arg(type)
    isMatchedSet <- type == "MatchedSet"
    
    for(i in 1 : length(celFilenames))
    {
        y <- as.matrix(read.celfile(as.character(celFilenames[i]), intensity.means.only = TRUE)[["INTENSITY"]][["MEAN"]][snpProbeInfo$probeIndex])
        y <- log2(y)
        allAint <- y[snpProbeInfo$isAAllele, 1, drop = FALSE]
        allBint <- y[!snpProbeInfo$isAAllele, 1, drop = FALSE]
        allAint <- subColSummarizeMedian(
            matrix(allAint, ncol = 1),
            snpProbeInfo$snpId[snpProbeInfo$isAAllele])
        allBint <- subColSummarizeMedian(
            matrix(allBint, ncol = 1),
            snpProbeInfo$snpId[!snpProbeInfo$isAAllele])
        
        if(isMatchedSet)
        {
            tmp <- allAint
            tmp[allAint < allBint] <- allBint[allAint < allBint]
        }
        else
        {
            tmp <- (allAint + allBint) / 2
        }
        
        if(i == 1)
        {
            if(isMatchedSet)
            {
                xlab <- "Matched sequence intensity"
            }
            else
            {
                xlab <- "Average intensity"
            }
            
            plot(density(tmp), xlab = xlab, ylab = "", main = "")
        }
        else
        {
            lines(density(tmp))
        }
    }
}

.createAppendResultsToCSVFunction <- function(fileName)
{
    # if the file already exists we should delete it (without this we would end
    # up appending to the exising file which we don't want)
    if(file.exists(fileName))
    {
        file.remove(fileName)
    }
    
    function(snpProbesetInfoChunk, genotypingResultsChunk)
    {
        # the genotyping results are a list of matrices right now. We need
        # to bind everything together into a single matrix.
        fileExists <- file.exists(fileName)
        if(!fileExists)
        {
            # start by ensuring that column names are unique and meaningful if
            # we're working on a new file
            for(currName in names(genotypingResultsChunk))
            {
                colnames(genotypingResultsChunk[[currName]]) <-
                    paste(currName, colnames(genotypingResultsChunk[[currName]]), sep = "_")
            }
        }
        
        uberMatrix <- snpProbesetInfoChunk
        for(currMat in genotypingResultsChunk)
        {
        	uberMatrix <- cbind(uberMatrix, currMat)
        }
        
        write.table(
            uberMatrix,
            file = fileName,
            append = fileExists,
            sep = ",",
            row.names = FALSE,
            col.names = !fileExists,
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
expandCelFiles <- function(filename)
{
    retVal <- NULL
    if(!file.exists(filename))
    {
        stop("failed to find file named \"", filename, "\"")
    }
    else if(file.info(filename)$isdir)
    {
        # if it's a dir expand it to all CEL file contents
        fileListing <- dir(filename, full.names = TRUE)
        retVal <- fileListing[grep("*.CEL$", toupper(fileListing))]
    }
    else if(grep("*.CEL$", toupper(filename)))
    {
        retVal <- filename
    }
    
    retVal
}
