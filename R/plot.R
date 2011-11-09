plotSnpContAvg <- function(
        snpContAvg,
        renderMask=matrix(TRUE, nrow=nrow(snpContAvg), ncol=ncol(snpContAvg)),
        xlab="Contrast",
        ylab="Average Intensity",
        main="",
        col=rgb(0, 0, 0, 0.25),
        pch=20,
        xlim=c(-max(snpContAvg$snpContrasts), max(snpContAvg$snpContrasts)),
        ylim=c(0, max(snpContAvg$snpAverages)),
        ...) {

    xVals <- c(snpContAvg$snpContrasts)
    yVals <- c(snpContAvg$snpAverages)
    
    renderMask <- c(renderMask)
    if(!is.null(renderMask)) {
        xVals <- xVals[renderMask]
        yVals <- yVals[renderMask]
    }
    
    plot(
        xVals, yVals,
        xlab=xlab,
        ylab=ylab,
        main=main,
        xlim=xlim,
        ylim=ylim,
        col=col,
        pch=pch,
        ...)
}

pointsSnpContAvg <- function(
        snpContAvg,
        renderMask=matrix(TRUE, nrow=nrow(snpContAvg), ncol=ncol(snpContAvg)),
        col=rgb(0, 0, 1),
        pch=20,
        ...) {
    
    xVals <- c(snpContAvg$snpContrasts)
    yVals <- c(snpContAvg$snpAverages)
    
    renderMask <- c(renderMask)
    if(!is.null(renderMask)) {
        xVals <- xVals[renderMask]
        yVals <- yVals[renderMask]
    }
    
    points(
        xVals, yVals,
        col=col,
        pch=pch,
        ...)
}

boxplot.snpIntensities <- function(x, ...) {
    plotList <- list()
    for(sampleIndex in seq_len(ncol(x))) {
        name <- x$sampleNames[sampleIndex]
        aName <- paste(name, ": A", sep="")
        bName <- paste(name, ": B", sep="")
        plotList[[aName]] <- x$A[ , sampleIndex, drop=TRUE]
        plotList[[bName]] <- x$B[ , sampleIndex, drop=TRUE]
    }
    
    args <- list(...)
    tempPar <- par()
    par(
        las=if("las" %in% names(args)) args$las else 3,
        mar=if("mar" %in% names(args)) args$mar else c(11, 4, 4, 2))
    boxplot(plotList, ...)
    par(las=tempPar$las, mar=tempPar$mar)
}

plotMouseDivArrayImage <- function(celFilename) {
    celData <- as.matrix(read.celfile(celFilename, intensity.means.only = TRUE)[["INTENSITY"]][["MEAN"]][1:6892960])
    image(
        1 : 2572,
        1 : 2680,
        matrix(log2(celData), byrow = T, ncol = 2680),
        col = c("green", "black", "red"),
        xlab = "",
        ylab = "")
}

mouseDivDensityPlot <- function(celFilenames, snpProbeInfo, type = c("Average", "MatchedSet")) {
    snpProbeInfo$snpId <- as.factor(snpProbeInfo$snpId)
    
    type <- match.arg(type)
    isMatchedSet <- type == "MatchedSet"
    
    for(i in 1 : length(celFilenames)) {
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
        
        if(isMatchedSet) {
            tmp <- allAint
            tmp[allAint < allBint] <- allBint[allAint < allBint]
        } else {
            tmp <- (allAint + allBint) / 2
        }
        
        if(i == 1) {
            if(isMatchedSet) {
                xlab <- "Matched sequence intensity"
            } else {
                xlab <- "Average intensity"
            }
            
            plot(density(tmp), xlab = xlab, ylab = "", main = "")
        } else {
            lines(density(tmp))
        }
    }
}
