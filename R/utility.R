plotSNP = function(nm, ns, geno, vino, main = "", xlab = "Contrast", ylab = "Average Intensity") {
    plot(nm, ns, xlab = xlab, ylab = ylab, main = main)
    points(nm[geno == 1], ns[geno == 1], col = 2)
    points(nm[geno == 2], ns[geno == 2], col = 3)
    points(nm[geno == 3], ns[geno == 3], col = 4)
    points(nm[vino == 1], ns[vino == 1], pch = 19, col = 6)
}

imageplot = function(filename, plotname = "imageplot.jpg") {
    library(affyio)
    my <- as.matrix(read.celfile(filename, intensity.means.only = TRUE)[["INTENSITY"]][["MEAN"]][1:6892960])
    y = matrix(log2(my), byrow = T, ncol = 2680)
    jpeg(plotname, width = 1200, heigh = 1200)
    image(c(1:2572), c(1:2680), y, col = c("green", "black", "red"), xlab = "", ylab = "")
    dev.off()
}

densityplot = function(filenames, plotname = "densityplot.jpg") {
    library(affyio)
    library("preprocessCore")
    
    nfile = length(filenames)
    jpeg(plotname, width = 1200, heigh = 1200)
    y = as.matrix(read.celfile(as.character(filenames[1]), intensity.means.only = TRUE)[["INTENSITY"]][["MEAN"]][fid])
    y = log2(y)
    allAint = y[allAid, 1, drop = FALSE]
    allBint = y[allBid, 1, drop = FALSE]
    allAint <- subColSummarizeMedian(matrix(allAint, ncol = 1), allA)
    allBint <- subColSummarizeMedian(matrix(allBint, ncol = 1), allA)
    plot(density(c(allAint, allBint)), xlab = "", ylab = "")
    for (i in 2:nfile) {
        y = as.matrix(read.celfile(as.character(filenames[i]), intensity.means.only = TRUE)[["INTENSITY"]][["MEAN"]][fid])
        y = log2(y)
        allAint = y[allAid, 1, drop = FALSE]
        allBint = y[allBid, 1, drop = FALSE]
        allAint <- subColSummarizeMedian(matrix(allAint, ncol = 1), allA)
        allBint <- subColSummarizeMedian(matrix(allBint, ncol = 1), allA)
        lines(density(c(allAint, allBint)))
    }
    dev.off()
} 
