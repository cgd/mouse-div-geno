penCNVinput = function(chrid, mpos, exon1info, exon2info, celfiledir, 
    filenames, snpsavefiledir, exonoutfiledir, cnvoutfiledir, mchr) {
    exon1CNV = length(exon1info) > 0
    exon2CNV = length(exon2info) > 0
    if (exon1CNV) 
        igp(celfiledir, exonoutfiledir, outfilename = "exon1", exon1info$exon1id, 
            exon1info$chrid, exon1info$mpos, exon1info$exon1, exon1info$CGFLcorrection, 
            exon1info$reference, filenames, mchr)
    if (exon2CNV) 
        igp(celfiledir, exonoutfiledir, outfilename = "exon2", exon2info$exon2id, 
            exon2info$chrid, exon2info$mpos, exon2info$exon2, exon2info$CGFLcorrection, 
            exon2info$reference, filenames, mchr)
    for (chri in mchr) {
        xname = paste(snpsavefiledir, "/baf", chri, sep = "", collapse = "")
        load(xname)
        file.remove(xname)
        xname = paste(snpsavefiledir, "/llr", chri, sep = "", collapse = "")
        load(xname)
        file.remove(xname)
        pos = mpos[match(rownames(llr), names(mpos))]
        
        nn = rownames(baf)
        n = nrow(baf)
        samplename = colnames(llr)
        for (sampleid in 1:ncol(baf)) {
            allr = llr[, sampleid]
            if (exon1CNV) {
                xname2 = paste(exonoutfiledir, "/", gsub(".CEL", "CHR", filenames[sampleid]), 
                  chri, "exon1", sep = "", collapse = "")
                load(xname2)
                file.remove(xname2)
                allr = c(allr, inte)
                if (sampleid == 1) 
                  nn = c(nn, names(inte))
            }
            if (exon2CNV) {
                xname2 = paste(exonoutfiledir, "/", gsub(".CEL", "CHR", filenames[sampleid]), 
                  chri, "exon2", sep = "", collapse = "")
                load(xname2)
                file.remove(xname2)
                allr = c(allr, inte)
                if (sampleid == 1) 
                  nn = c(nn, names(inte))
            }
            spos = NULL
            if (sampleid == 1) {
                if (exon1CNV) 
                  spos = c(spos, exon1info$mpos[exon1info$chrid == chri])
                if (exon2CNV) 
                  spos = c(spos, exon2info$mpos[exon2info$chrid == chri])
                pos = c(pos, spos)
                exonsize = length(spos)
                n = n + exonsize
            }
            
            tmp = c(baf[, sampleid], rep(2, exonsize))
            s = cbind(nn, allr, tmp)
            colnames(s) = c("Name", paste(gsub("CEL", "Log R Ratio", samplename[sampleid]), 
                sep = "", collapse = ""), paste(gsub("CEL", "B Allele Freq", samplename[sampleid]), 
                sep = "", collapse = ""))
            xname1 = paste(cnvoutfiledir, "/sample", sampleid, "chr", chri, ".txt", 
                sep = "", collapse = "")
            write.table(s, file = xname1, row.names = FALSE, sep = "\t", quote = FALSE)
        }
        
        PFB = c(apply(baf, 1, mean), rep(2, exonsize))
        pfb = cbind(nn, rep(chri, n), pos, PFB)
        colnames(pfb) = c("Name", "Chr", "Position", "PFB")
        xname1 = paste(cnvoutfiledir, "/sampleChr", chri, ".pfb", sep = "", collapse = "")
        write.table(pfb, file = xname1, row.names = FALSE, sep = "\t", quote = FALSE)
    }
}


igp = function(celfiledir, outfiledir, outfilename, id, chrid, mpos, ename, 
    CGFLcorrection, reference, filenames, mchr) {
    library("affyio")
    library("preprocessCore")
    setwd(celfiledir)
    if (missing(id)) 
        stop("No CDF file information")
    nfile = length(filenames)
    for (i in 1:nfile) {
        y = as.matrix(read.celfile(as.character(filenames[i]), intensity.means.only = TRUE)[["INTENSITY"]][["MEAN"]][id])
        y = log2(y)
        y = y + CGFLcorrection
        y <- normalize.quantiles.use.target(y, target = reference)
        y <- subColSummarizeMedian(matrix(y, ncol = 1), ename)
        for (chri in mchr) {
            xname2 = paste(outfiledir, "/", gsub(".CEL", "CHR", filenames[i]), chri, 
                outfilename, sep = "", collapse = "")
            inte = y[chrid == chri, ]
            save(inte, file = xname2)
        }
    }
} 
