pennCNVinput <- function(chrid, mpos, exon1info, exon2info, celfiledir, 
    filenames, snpsavefiledir, exonoutfiledir, cnvoutfiledir, mchr) {
    exon1CNV = length(exon1info) > 0
    exon2CNV = length(exon2info) > 0
    if (exon1CNV) 
        e1 = igp(celfiledir, exonoutfiledir, outfilename = "exon1", exon1info$exon1id, 
            exon1info$chrid, exon1info$mpos, exon1info$exon1, exon1info$CGFLcorrection, 
            exon1info$reference, filenames, mchr)
    if (exon2CNV) 
        e2 = igp(celfiledir, exonoutfiledir, outfilename = "exon2", exon2info$exon2id, 
            exon2info$chrid, exon2info$mpos, exon2info$exon2, exon2info$CGFLcorrection, 
            exon2info$reference, filenames, mchr)
    allbaf = alllrr = allpos = chrid = NULL
    for (chri in mchr) {
        xname = paste(snpsavefiledir, "/baf", chri, sep = "", collapse = "")
        load(xname)
        file.remove(xname)
        xname = paste(snpsavefiledir, "/lrr", chri, sep = "", collapse = "")
        load(xname)
        file.remove(xname)
        pos = mpos[match(rownames(lrr), names(mpos))]
        allpos = c(allpos, pos)
        allbaf = rbind(allbaf, baf)
        alllrr = rbind(alllrr, lrr)
        chrid = c(chrid, rep(chri, nrow(baf)))
    }
    nn = rownames(allbaf)
    exonsize = 0
    
    for (sampleid in 1:ncol(baf)) {
        alrr = alllrr[, sampleid]
        if (exon1CNV) {
            xname2 = paste(exonoutfiledir, "/", gsub(".CEL", "exon1", filenames[sampleid]), 
                sep = "", collapse = "")
            load(xname2)
            file.remove(xname2)
            alrr = c(alrr, log2(y/e1))
            if (sampleid == 1) {
                ny = names(y)
                nn = c(nn, ny)
                g = match(ny, names(exon1info$mpos))
                allpos = c(allpos, exon1info$mpos[g])
                chrid = c(chrid, exon1info$chrid[g])
                exonsize = exonsize + length(y)
            }
        }
        if (exon2CNV) {
            xname2 = paste(exonoutfiledir, "/", gsub(".CEL", "exon2", filenames[sampleid]), 
                sep = "", collapse = "")
            load(xname2)
            file.remove(xname2)
            alrr = c(alrr, log2(y/e2))
            if (sampleid == 1) {
                ny = names(y)
                nn = c(nn, ny)
                g = match(ny, names(exon2info$mpos))
                allpos = c(allpos, exon2info$mpos[g])
                chrid = c(chrid, exon2info$chrid[g])
                exonsize = exonsize + length(y)
            }
        }
        
        tmp = c(allbaf[, sampleid], rep(2, exonsize))
        s = cbind(nn, alrr, tmp)
        colnames(s) = c("Name", paste(gsub("CEL", "Log R Ratio", filenames[sampleid]), 
            sep = "", collapse = ""), paste(gsub("CEL", "B Allele Freq", filenames[sampleid]), 
            sep = "", collapse = ""))
        xname1 = paste(cnvoutfiledir, "/cnv", gsub(".CEL", ".txt", filenames[sampleid]), 
            sep = "", collapse = "")
        write.table(s, file = xname1, row.names = FALSE, sep = "\t", quote = FALSE)
    }
    
    PFB = c(apply(allbaf, 1, mean), rep(2, exonsize))
    pfb = cbind(nn, chrid, allpos, PFB)
    colnames(pfb) = c("Name", "Chr", "Position", "PFB")
    xname1 = paste(cnvoutfiledir, "/cnv.pfb", sep = "", collapse = "")
    write.table(pfb, file = xname1, row.names = FALSE, sep = "\t", quote = FALSE)
}

# calculates LRR and BAF values for a single SNP position
# PARAMETERS:
#   intensityConts:
#       per-sample intensity contrasts (A allele vs B allele)
#   intensityAvgs:
#       per-sample intensity averages (A averaged with B)
#   genos:
#       per-sample genotypes
# RETURNS:
#   A data.frame object with BAF adn LRR components.
#   The row count of this dataframe will equal length(genos)
calcLRRAndBAF <- function(intensityConts, intensityAvgs, genos, mm, ms, ss)
{
    sampleCount <- length(genos)
    uniqueGenos <- sort(unique(genos))
    uniqueGenoCount <- length(uniqueGenos)
    BAF <- rep(0, sampleCount)
    LRR <- rep(0, sampleCount)
    
    medianContPerGeno <- rep(0, 3)
    ms <- rep(0, 3)
    for(genoCode in uniqueGenos)
    {
        genoIndices <- which(genos == genoCode)
        medianContPerGeno[genoCode] <- median(intensityConts[genoIndices])
        ms[genoCode] <- median(intensityAvgs[genoIndices])
    }
    
    if (uniqueGenoCount == 3) {
        c1 = sqrt(ss[[1]][1, 1])
        c2 = sqrt(ss[[2]][1, 1])
        c3 = sqrt(ss[[3]][1, 1])
        s1 = sqrt(ss[[1]][2, 2])
        s2 = sqrt(ss[[2]][2, 2])
        s3 = sqrt(ss[[3]][2, 2])
        tmp = intensityConts >= medianContPerGeno[1]
        if (any(tmp)) {
            BAF[tmp] = 0
            LRR[tmp] = log2(intensityAvgs[tmp]/intensityAvgs[1])
        }
        tmp = intensityConts < medianContPerGeno[1] & intensityConts > medianContPerGeno[2]
        if (any(tmp)) {
            k1 = (medianContPerGeno[1] - intensityConts[tmp])/c1
            k2 = (intensityConts[tmp] - medianContPerGeno[2])/c2
            BAF[tmp] = 0.5 * k1/(k1 + k2)
            LRR[tmp] = log2(intensityAvgs[tmp]/((k2 * intensityAvgs[1] + k1 * intensityAvgs[2])/(k1 + k2)))
        }
        tmp = intensityConts <= medianContPerGeno[2] & intensityConts > medianContPerGeno[3]
        if (any(tmp)) {
            k1 = (medianContPerGeno[2] - intensityConts[tmp])/c2
            k2 = (intensityConts[tmp] - medianContPerGeno[3])/c3
            BAF[tmp] = 0.5 + 0.5 * k1/(k1 + k2)
            LRR[tmp] = log2(intensityAvgs[tmp]/((k2 * intensityAvgs[2] + k1 * intensityAvgs[3])/(k1 + k2)))
        }
        tmp = intensityConts <= medianContPerGeno[3]
        if (any(tmp)) {
            BAF[tmp] = 1
            LRR[tmp] = log2(intensityAvgs[tmp]/intensityAvgs[3])
        }
    }
    else if (uniqueGenoCount == 2) {
        id = sum(uniqueGenos)
        if (id == 4) {
            c1 = sqrt(ss[[1]][1, 1])
            c3 = sqrt(ss[[3]][1, 1])
            s1 = sqrt(ss[[1]][2, 2])
            s3 = sqrt(ss[[3]][2, 2])
            tmp = intensityConts >= medianContPerGeno[1]
            if (any(tmp)) {
                BAF[tmp] = 0
                LRR[tmp] = log2(intensityAvgs[tmp]/intensityAvgs[1])
            }
            tmp = intensityConts < medianContPerGeno[1] & intensityConts > medianContPerGeno[3]
            if (any(tmp)) {
                k1 = (medianContPerGeno[1] - intensityConts[tmp])/c1
                k3 = (intensityConts[tmp] - medianContPerGeno[3])/c3
                BAF[tmp] = k1/(k1 + k3)
                LRR[tmp] = log2(intensityAvgs[tmp]/((k3 * intensityAvgs[1] + k1 * intensityAvgs[3])/(k1 + k3)))
            }
            tmp = intensityConts <= medianContPerGeno[3]
            if (any(tmp)) {
                BAF[tmp] = 1
                LRR[tmp] = log2(intensityAvgs[tmp]/intensityAvgs[3])
            }
        }
        else if (id == 5) {
            c2 = sqrt(ss[[2]][1, 1])
            c3 = sqrt(ss[[3]][1, 1])
            s2 = sqrt(ss[[2]][2, 2])
            s3 = sqrt(ss[[3]][2, 2])
            tmp = intensityConts >= medianContPerGeno[2]
            if (any(tmp)) {
                BAF[tmp] = 0.5
                LRR[tmp] = log2(intensityAvgs[tmp]/intensityAvgs[2])
            }
            tmp = intensityConts < medianContPerGeno[2] & intensityConts > medianContPerGeno[3]
            if (any(tmp)) {
                k2 = (medianContPerGeno[2] - intensityConts[tmp])/c2
                k3 = (intensityConts[tmp] - medianContPerGeno[3])/c3
                BAF[tmp] = 0.5 + 0.5 * k2/(k2 + k3)
                LRR[tmp] = log2(intensityAvgs[tmp]/((k3 * intensityAvgs[2] + k2 * intensityAvgs[3])/(k2 + k3)))
            }
            tmp = intensityConts <= medianContPerGeno[3]
            if (any(tmp)) {
                BAF[tmp] = 1
                LRR[tmp] = log2(intensityAvgs[tmp]/intensityAvgs[3])
            }
        }
        else if (id == 3) {
            c1 = sqrt(ss[[1]][1, 1])
            c2 = sqrt(ss[[2]][1, 1])
            s1 = sqrt(ss[[1]][2, 2])
            s2 = sqrt(ss[[2]][2, 2])
            tmp = intensityConts >= medianContPerGeno[1]
            if (any(tmp)) {
                BAF[tmp] = 0
                LRR[tmp] = log2(intensityAvgs[tmp]/intensityAvgs[1])
            }
            tmp = intensityConts < medianContPerGeno[1] & intensityConts > medianContPerGeno[2]
            if (any(tmp)) {
                k1 = (medianContPerGeno[1] - intensityConts[tmp])/c1
                k2 = (intensityConts[tmp] - medianContPerGeno[2])/c2
                BAF[tmp] = 0.5 * k1/(k1 + k2)
                LRR[tmp] = log2(intensityAvgs[tmp]/((k1 * intensityAvgs[2] + k2 * intensityAvgs[1])/(k1 + k2)))
            }
            tmp = intensityConts <= medianContPerGeno[2]
            if (any(tmp)) {
                BAF[tmp] = 0.5
                LRR[tmp] = log2(intensityAvgs[tmp]/intensityAvgs[2])
            }
        }
    }
    else if (uniqueGenoCount == 1) {
        if (uniqueGenos == 1) {
            BAF = rep(0, sampleCount)
            LRR = log2(intensityAvgs/intensityAvgs[1])
        }
        else if (uniqueGenos == 2) {
            BAF = rep(0.5, sampleCount)
            LRR = log2(intensityAvgs/intensityAvgs[2])
        }
        else if (uniqueGenos == 3) {
            BAF = rep(1, sampleCount)
            LRR = log2(intensityAvgs/intensityAvgs[3])
        }
    }
    
    data.frame(BAF = BAF, LRR = LRR)
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
        g = match(chrid, mchr)
        y = y[match(names(chrid[!is.na(g)]), rownames(y)), ]
        xname2 = paste(outfiledir, "/", gsub(".CEL", outfilename, filenames[i]), 
            sep = "", collapse = "")
        save(y, file = xname2)
        if (i == 1) 
            ty = y
        else ty = ty + y
    }
    ty = ty/nfile
    ty
}

simpleigp = function(celfiledir, id, chrid, mpos, ename, CGFLcorrection, 
    reference, filenames, mchr) {
    library("affyio")
    library("preprocessCore")
    setwd(celfiledir)
    
    y = as.matrix(read.celfile(as.character(filenames), intensity.means.only = TRUE)[["INTENSITY"]][["MEAN"]][id])
    y = log2(y)
    y = y + CGFLcorrection
    y <- normalize.quantiles.use.target(y, target = reference)
    y <- subColSummarizeMedian(matrix(y, ncol = 1), ename)
    out = list()
    for (chri in mchr) {
        pos = mpos[chrid == chri]
        ey = y[chrid == chri]
        out[[chri]] = list(pos = pos, ey = ey)
    }
    out
}

simpleCNV = function(allid, ABid, chrid, CGFLcorrection, reference, exon1info, 
    exon2info, celfiledir, filenames, cnvoutfiledir, mchr = c(1:19), stname, refsample) {
    if (missing(celfiledir)) 
        celfiledir = getwd()
    if (missing(cnvoutfiledir)) 
        outfiledir = getwd()
    if (missing(stname)) 
        stname = filenames
    refid = match(refsample, filenames)
    stname = as.character(stname)
    simpleCNVdata(allid, ABid, chrid, CGFLcorrection, reference, exon1info, exon2info, 
        celfiledir, filenames, cnvoutfiledir, mchr)
    simpleCNVsummary(cnvoutfiledir, filenames, refid, mchr, stname)
}

simpleCNVdata = function(allid, ABid, chrid, CGFLcorrection, reference, 
    exon1info, exon2info, celfiledir, filenames, cnvoutfiledir, mchr) {
    library("affyio")
    library("preprocessCore")
    setwd(celfiledir)
    SNPname = ABid$SNPname
    Aid = ABid$allAid
    Bid = ABid$allBid
    rm(ABid)
    mpos = chrid$mpos
    chrid = chrid$chrid
    nfile = length(filenames)
    for (i in 1:nfile) {
        y = as.matrix(read.celfile(as.character(filenames[i]), intensity.means.only = TRUE)[["INTENSITY"]][["MEAN"]][allid])
        y = log2(y)
        if (length(CGFLcorrection) > 0) 
            y = y + CGFLcorrection
        if (length(reference) > 0) 
            y = normalize.quantiles.use.target(y, target = reference)
        oAint = y[Aid, 1, drop = FALSE]
        oBint = y[Bid, 1, drop = FALSE]
        allAint = subColSummarizeMedian(matrix(oAint, ncol = 1), SNPname)
        allBint = subColSummarizeMedian(matrix(oBint, ncol = 1), SNPname)
        tmp = names(allAint[allAint > allBint, ])
        g = match(SNPname, tmp)
        B = oAint
        B[is.na(g)] = oBint[is.na(g)]
        B = normalize.quantiles.use.target(B, target = exon1info$reference)
        B = subColSummarizeMedian(matrix(B, ncol = 1), SNPname)
        
        exon1 = simpleigp(celfiledir, exon1info$exon1id, exon1info$chrid, exon1info$mpos, 
            exon1info$exon1, exon1info$CGFLcorrection, exon1info$reference, filenames[i], 
            mchr)
        exon2 = simpleigp(celfiledir, exon2info$exon2id, exon2info$chrid, exon2info$mpos, 
            exon2info$exon2, exon2info$CGFLcorrection, exon2info$reference, filenames[i], 
            mchr)
        for (chri in mchr) {
            pos = c(mpos[chrid == chri], exon1[[chri]]$pos, exon2[[chri]]$pos)
            id = c(rep(1, sum(chrid == chri)), rep(2, length(exon1[[chri]]$pos)), 
                rep(3, length(exon2[[chri]]$pos)))
            inte = c(B[chrid == chri], exon1[[chri]]$ey, exon2[[chri]]$ey)
            o = order(pos)
            pos = pos[o]
            id = id[o]
            inte = inte[o]
            xname2 = paste(cnvoutfiledir, "/", gsub(".CEL", "Chr", filenames[i]), 
                chri, sep = "", collapse = "")
            save(inte, file = xname2)
            if (i == 1) {
                xname2 = paste(cnvoutfiledir, "/pos", chri, sep = "", collapse = "")
                save(pos, id, file = xname2)
            }
        }
    }
}

simpleCNVsummary = function(cnvoutfiledir, filenames, refid, mchr, stname) {
    library(HiddenMarkov)
    k = 0.9999
    pi = matrix(c(k, (1 - k)/2, (1 - k)/2, (1 - k)/2, k, (1 - k)/2, (1 - k)/2, (1 - 
        k)/2, k), 3, 3)
    delta = c(0, 1, 0)
    nfile = length(filenames)
    for (chri in mchr) {
        xname2 = paste(cnvoutfiledir, "/pos", chri, sep = "", collapse = "")
        load(xname2)
        xname2 = paste(cnvoutfiledir, "/", gsub(".CEL", "Chr", filenames[refid]), 
            chri, sep = "", collapse = "")
        load(xname2)
        ref = inte
        cnv = cnvtable = NULL
        for (i in 1:nfile) {
            xname2 = paste(cnvoutfiledir, "/", gsub(".CEL", "Chr", filenames[i]), 
                chri, sep = "", collapse = "")
            load(xname2)
            a = inte/ref
            m = mean(a)
            b = dthmm(a, pi, delta, "norm", list(mean = c(m - 0.1, m, m + 0.1), sd = c(0.05, 
                0.05, 0.05)))
            states <- Viterbi(b)
            cnv = cbind(cnv, states)
            
            l = cbind(which(states == 1), pos[states == 1])
            g = cbind(which(states == 3), pos[states == 3])
            ll = length(l)
            lg = length(g)
            if (lg > 2) {
                k = diff(g[, 1])
                ep = which(k > 1)
                sp = c(1, ep + 1)
                ep = c(ep, lg/2)
                s = length(sp)
                if (s == 1) 
                  cnvtable = rbind(cnvtable, c(stname[i], "gain", g[sp, 2], g[ep, 
                    2], g[ep, ] - g[sp, ] + 1, sum(id[g[sp, 1]:g[ep, 1]] == 1), sum(id[g[sp, 
                    1]:g[ep, 1]] == 2), sum(id[g[sp, 1]:g[ep, 1]] == 3), round(mean(ref[g[sp, 
                    1]:g[ep, 1]]), 3), round(mean(inte[g[sp, 1]:g[ep, 1]]), 3)))
                if (s > 1) {
                  k = cbind(g[sp, 2], g[ep, 2], g[ep, ] - g[sp, ] + 1, g[sp, 1], 
                    g[ep, 1])
                  k1 = apply(k, 1, function(x, ref, id, inte) c(sum(id[x[6]:x[5]] == 
                    1), sum(id[x[6]:x[5]] == 2), sum(id[x[6]:x[5]] == 3), mean(ref[x[6]:x[5]]), 
                    mean(inte[x[6]:x[5]])), ref, id, inte)
                  cnvtable = rbind(cnvtable, cbind(rep(stname[i], s), rep("gain", 
                    s), k[, 1:4], t(round(k1, 3))))
                }
            }
            if (ll > 2) {
                k = diff(l[, 1])
                ep = which(k > 1)
                sp = c(1, ep + 1)
                ep = c(ep, ll/2)
                s = length(sp)
                if (s == 1) 
                  cnvtable = rbind(cnvtable, c(stname[i], "loss", l[sp, 2], l[ep, 
                    2], l[ep, ] - l[sp, ] + 1, sum(id[l[sp, 1]:l[ep, 1]] == 1), sum(id[l[sp, 
                    1]:l[ep, 1]] == 2), sum(id[l[sp, 1]:l[ep, 1]] == 3), round(mean(ref[l[sp, 
                    1]:l[ep, 1]]), 3), round(mean(inte[l[sp, 1]:l[ep, 1]]), 3)))
                if (s > 1) {
                  k = cbind(l[sp, 2], l[ep, 2], l[ep, ] - l[sp, ] + 1, l[sp, 1], 
                    l[ep, 1])
                  k1 = apply(k, 1, function(x, ref, id, inte) c(sum(id[x[6]:x[5]] == 
                    1), sum(id[x[6]:x[5]] == 2), sum(id[x[6]:x[5]] == 3), mean(ref[x[6]:x[5]]), 
                    mean(inte[x[6]:x[5]])), ref, id, inte)
                  cnvtable = rbind(cnvtable, cbind(rep(stname[i], s), rep("loss", 
                    s), k[, 1:4], t(round(k1, 3))))
                }
            }
        }
        colnames(cnvtable) = c("Name", "Status", "StartPosition", "EndPosition", 
            "Number of Probe sets", "Size", "Number of SNP probe sets", "Number of exon1 sets", 
            "Number of exon2 sets", "Mean intensity of reference sample", "Mean intensity of sample")
        colnames(cnv) = stname
        xname = paste("~/cnvChr", chri, sep = "", collapse = "")
        save(cnv, cnvtable, file = xname)
    }
} 
