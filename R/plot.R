plotSNP <- function(nm, ns, geno, vino, main = "", xlab = "Contrast", ylab = "Average Intensity") {
    plot(nm, ns, xlab = xlab, ylab = ylab, main = main)
    points(nm[geno == 1], ns[geno == 1], col = 2)           # red
    points(nm[geno == 2], ns[geno == 2], col = 3)           # green
    points(nm[geno == 3], ns[geno == 3], col = 4)           # blue
    points(nm[vino == 1], ns[vino == 1], pch = 19, col = 6) # filled point
}

plotArrayAvgContDensity <- function(
    arrData,
    main="Mean Intensity vs. Intensity Contrast",
    xlab="Intensity Contrast",
    ylab="Mean Intensity",
    useAbs=TRUE,
    pch=20) {
    
    contFun <-
        if(useAbs) {
            function(x) {abs(x$intensityConts)}
        } else {
            function(x) {x$intensityConts}
        }
    
    allCont <- unlist(lapply(arrData, contFun))
    allAvgs <- unlist(lapply(arrData, function(x) {x$intensityAvgs}))
    
    colors <- densCols(allCont, allAvgs, col=topo.colors)
    plot(allCont, allAvgs, col=colors, main=main, xlab=xlab, ylab=ylab, pch=pch)
}

plotMouseDivArrayImage <- function(celFilename) {
    my <- as.matrix(read.celfile(celFilename, intensity.means.only = TRUE)[["INTENSITY"]][["MEAN"]][1:6892960])
    image(
        1 : 2572,
        1 : 2680,
        matrix(log2(my), byrow = T, ncol = 2680),
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
