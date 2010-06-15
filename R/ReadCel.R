#########################################################################
# ReadCel.R
#
# Part of the MouseDivGenotype package
#
# This is the function to read .CEL files and normalize/transform the data, and save the normalized intensity by chromosomes.
# It also computes gender, if needed.
# input : dir - directory that .CEL files stored.
#         outfiledir - directory that you want to save normalized log2 intensity. File saves by chromosomes.
#         allid - index file indicating SNP position (download from cgd.jax.org)
#         ABid - A and B allele position and name (download from cgd.jax.org)
#         CGFLcorrection - C+G in 25 mer flanking sequences and restriction enzyme fragment length correction coefficients.
#         reference - targer reference distribution used for quantile normalization.
#         trans - two possible transformation, i.e. CCS and MA transformation
#         celnamefile - name of file containing .CEL files that you want to read. First column should be .CEL name and
#            second column (optional) should be gender information. The file must have a column header.
#
# output : isMale = computed gender
#
#########################################################################
ReadCelFile = function(celfiledir, outfiledir, allid, ABid, chrid, CGFLcorrection = NULL, 
    reference = NULL, trans = c("CCStrans", "MAtrans"), celnamefile = NULL, subset = FALSE) {
    library(affyio)
    library("preprocessCore")
    trans = match.arg(trans)
    if (missing(celfiledir)) 
        celfiledir = getwd()
    if (missing(outfiledir)) 
        outfiledir = getwd()
    if (missing(allid) | missing(ABid)) 
        stop("No CDF file information")
    SNPname = ABid$SNPname
    Aid = ABid$allAid
    Bid = ABid$allBid
    rm(ABid)
    if (subset) 
        chrid = chrid$schrid
    else
        chrid = chrid$chrid
    
    setwd(celfiledir)
    gender = NULL
    if (length(celnamefile) == 0) {
        tmp = dir()
        isCel = apply(matrix(tmp, ncol = 1), 1, function(x) length(grep(".CEL", x)) > 
            0)
        filenames = tmp[isCel]
    }
    else {
        filenames = read.table(celnamefile, header = TRUE)
        if (ncol(filenames) > 1) 
            gender = filenames[, 2]
        filenames = filenames[, 1]
    }
    nfile = length(filenames)
    mchr = c(1:19, "X", "Y", "M")
    
    for (i in 1:nfile) {
        y = as.matrix(read.celfile(as.character(filenames[i]), intensity.means.only = TRUE)[["INTENSITY"]][["MEAN"]][allid])
        y = log2(y)
        if (length(CGFLcorrection) > 0)
            # C+G and fragment length correction y
            y = y + CGFLcorrection
        if (length(reference) > 0) 
            y <- normalize.quantiles.use.target(y, target = reference)
        allAint = y[Aid, 1, drop = FALSE]
        allBint = y[Bid, 1, drop = FALSE]
        allAint <- subColSummarizeMedian(matrix(allAint, ncol = 1), SNPname)
        allBint <- subColSummarizeMedian(matrix(allBint, ncol = 1), SNPname)
        if (is.element(trans, "CCStrans")) {
            # fixed K??
            res = ccstrans(2^allAint, 2^allBint)
            M = res$x
            S = res$y
        }
        else if (is.element(trans, "MAtrans")) {
            # then prior??
            M = allAint - allBint
            S = (allAint + allBint)/2
        }
        xname1 = paste(outfiledir, "/", gsub(".CEL", "intensity", filenames[i]), sep = "", collapse = "")
        save(M, S, file = xname1)
        for (chri in mchr) {
            xname2 = paste(outfiledir, "/", gsub(".CEL", "CHR", filenames[i]), chri, sep = "", collapse = "")
            MM1 = M[chrid == chri, ]
            SS1 = S[chrid == chri, ]
            save(MM1, SS1, file = xname2)
        }
    }
    
    autoint = NULL
    for (chri in mchr) {
        MM = SS = NULL
        for (i in 1:nfile) {
            xname1 = paste(outfiledir, "/", gsub(".CEL", "CHR", filenames[i]), chri, sep = "", collapse = "")
            # origa,origb
            load(xname1)
            MM = cbind(MM, MM1)
            SS = cbind(SS, SS1)
            rm(MM1, SS1)
            file.remove(xname1)
        }
        if (!is.na(match(chri, c(1:19)))) 
            autoint = c(autoint, mean(SS[sample(c(1:nrow(SS)), 100), ]))
        xname = paste(outfiledir, "/rawdataMMchr", chri, sep = "", collapse = "")
        colnames(MM) = colnames(SS) = filenames
        save(MM, file = xname)
        if (chri == "X") 
            intX = apply(SS, 2, mean)
        if (chri == "Y") 
            intY = apply(SS, 2, mean)
        xname = paste(outfiledir, "/rawdataSSchr", chri, sep = "", collapse = "")
        save(SS, file = xname)
    }
    if (length(gender) < nfile) 
        isMale = computegender(intX, intY, autoint)
    else {
        gender1 = rep(-1, nfile)
        gender1[gender == "male"] = 1
        gender1[gender == "female"] = 2
        if (length(gender1 == 1) > 1 & length(gender1 == 2) > 1) {
            gender = vdist(intY, intX, gender1)
            isMale = gender == 1
        }
        else isMale = computegender(intX, intY, autoint)
    }
    isMale
}
computegender = function(intX, intY, autoint) {
    if (length(autoint) < 1) 
        autoint = 7.362069
    n = length(intX)
    gender = kmeans(cbind(intY, intX), 2)$cluster
    y1 = mean(intY[gender == 1])
    y2 = mean(intY[gender == 2])
    if (y1 > y2) {
        if (y2 >= min(autoint)) 
            isMale = rep(TRUE, n)
        else isMale = gender == 1
    }
    else {
        if (y1 >= min(autoint)) 
            isMale = rep(TRUE, n)
        else isMale = gender == 2
    }
    if (abs(y1 - y2) < 0.1) {
        if (y1 > min(autoint)) 
            isMale = rep(TRUE, n)
        else isMale = rep(FALSE, n)
    }
    isMale
}
ccstrans = function(a, b, k = 4) {
    x = asinh(k * (a - b)/(a + b))/asinh(k)
    y = (log2(a) + log2(b))/2
    list(x = x, y = y)
}

#########################################################################
#
# genotypethis.R
#
# Part of the MouseDivGenotype package
#
# This is the function to genotype and vinotype the array
#
#########################################################################
genotypethis = function(savefiledir, MM, SS, hint = NULL, isMale, trans, chr, doCNV = FALSE) {
    library(cluster)
    n = ncol(MM)
    nn = nrow(MM)
    if (length(hint) == 0) 
        hint = rep(0, nn)
    if (missing(savefiledir)) 
        savefiledir = getwd()
    
    autosome = match(chr, c(1:19))
    if (!is.na(autosome)) {
        #cat("in autosome part of genotype this\n")
        
        # process in blocks of 2000 to balance efficiency and memory use
        conf = geno = vino = baf = llr = NULL
        sp = seq(1, nn, 2000)
        lls = length(sp)
        ep = c(sp[2:lls] - 1, nn)
        for (iter in 1:lls) {
            #cat("timing iteration ", iter, "\n", sep="")
            #startTime <- getTime()
            
            tmp = sp[iter]:ep[iter]
            genom = apply(cbind(MM[tmp, ], SS[tmp, ], hint[tmp]), 1, function(x, 
                trans, doCNV, n) {
                b = genotype(x[1:n], x[(n + 1):(2 * n)], x[(2 * n + 1):(2 * n + 1)], 
                  trans, doCNV)
                c(b$geno, b$vino, b$conf, b$baf, b$llr)
            }, trans, doCNV, n)
            genom = t(genom)
            conf = rbind(conf, genom[, c((2 * n + 1):(3 * n))])
            vino = rbind(vino, genom[, c((n + 1):(2 * n))])
            geno = rbind(geno, genom[, c(1:n)])
            if (doCNV) {
                baf = rbind(baf, genom[, c((3 * n + 1):(4 * n))])
                llr = rbind(llr, genom[, c((4 * n + 1):(5 * n))])
            }
            
            #timeReport(startTime)
        }
        rownames(geno) = rownames(conf) = rownames(vino) = rownames(MM)
        colnames(geno) = colnames(conf) = colnames(vino) = colnames(MM)
        xname = paste(savefiledir, "/genodata", chr, sep = "", collapse = "")
        save(geno, file = xname)
        xname = paste(savefiledir, "/vinodata", chr, sep = "", collapse = "")
        save(vino, file = xname)
        xname = paste(savefiledir, "/confdata", chr, sep = "", collapse = "")
        save(conf, file = xname)
        if (doCNV) {
            rownames(baf) = rownames(llr) = rownames(MM)
            colnames(baf) = colnames(llr) = colnames(MM)
            xname = paste(savefiledir, "/baf", chr, sep = "", collapse = "")
            save(baf, file = xname)
            xname = paste(savefiledir, "/llr", chr, sep = "", collapse = "")
            save(llr, file = xname)
        }
    }
    
    chrX = match(chr, "X")
    if (!is.na(chrX)) {
        conf = geno = vino = baf = llr = matrix(0, nn, n)
        pseudoautosomalID = c("JAX00723347", "JAX00723349", "JAX00723350", "JAX00723351", 
            "JAX00187253", "JAX00723353", "JAX00723354", "JAX00723355", "JAX00723356", 
            "JAX00723358", "JAX00723359", "JAX00723360", "JAX00723361", "JAX00723364", 
            "JAX00723365", "JAX00723366", "JAX00723367", "JAX00723368", "JAX00723369", 
            "JAX00723370", "JAX00187255", "JAX00187256", "JAX00725011", "JAX00725012", 
            "JAX00723371", "JAX00723372")
        pid = match(pseudoautosomalID, rownames(MM))
        pid = pid[!is.na(pid)]
        lp = length(pid)
        id = c(1:nn)
        if (lp > 0) 
            id = setdiff(c(1:nn), pid)
        sp = seq(1, length(id), 2000)
        lls = length(sp)
        ep = c(sp[2:lls] - 1, length(id))
        n = sum(!isMale)
        if (n > 1) {
            # for female three groups
            confm = genom = vinom = bafm = llrm = NULL
            sMM = MM[id, !isMale]
            sSS = SS[id, !isMale]
            shint = hint[id]
            for (iter in 1:lls) {
                tmp = sp[iter]:ep[iter]
                genom1 = apply(cbind(sMM[tmp, ], sSS[tmp, ], shint[tmp]), 1, function(x, 
                  trans, doCNV, n) {
                  b = genotype(x[1:n], x[(n + 1):(2 * n)], x[(2 * n + 1):(2 * n + 
                    1)], trans, doCNV)
                  c(b$geno, b$vino, b$conf, b$baf, b$llr)
                }, trans, doCNV, n)
                genom1 = t(genom1)
                confm = rbind(confm, genom1[, c((2 * n + 1):(3 * n))])
                vinom = rbind(vinom, genom1[, c((n + 1):(2 * n))])
                genom = rbind(genom, genom1[, c(1:n)])
                if (doCNV) {
                  bafm = rbind(bafm, genom1[, c((3 * n + 1):(4 * n))])
                  llrm = rbind(llrm, genom1[, c((4 * n + 1):(5 * n))])
                }
            }
            conf[id, !isMale] = confm
            vino[id, !isMale] = vinom
            geno[id, !isMale] = genom
            if (doCNV) {
                baf[id, !isMale] = bafm
                llr[id, !isMale] = llrm
            }
        }
        n = sum(isMale)
        if (n > 1) {
            # for male two groups
            confm = genom = vinom = bafm = llrm = NULL
            sMM = MM[id, isMale]
            sSS = SS[id, isMale]
            for (iter in 1:lls) {
                tmp = sp[iter]:ep[iter]
                genom1 = apply(cbind(sMM[tmp, ], sSS[tmp, ]), 1, function(x, trans, 
                  doCNV, n) {
                  b = genotypeSexchr(x[1:n], x[(n + 1):(2 * n)], trans, doCNV)
                  c(b$geno, b$vino, b$conf, b$baf, b$llr)
                }, trans, doCNV, n)
                genom1 = t(genom1)
                confm = rbind(confm, genom1[, c((2 * n + 1):(3 * n))])
                vinom = rbind(vinom, genom1[, c((n + 1):(2 * n))])
                genom = rbind(genom, genom1[, c(1:n)])
                if (doCNV) {
                  bafm = rbind(bafm, genom1[, c((3 * n + 1):(4 * n))])
                  llrm = rbind(llrm, genom1[, c((4 * n + 1):(5 * n))])
                }
            }
            conf[id, isMale] = confm
            vino[id, isMale] = vinom
            geno[id, isMale] = genom
            if (doCNV) {
                baf[id, isMale] = bafm
                llr[id, isMale] = llrm
            }
        }
        n = ncol(MM)
        if (lp > 1) {
            confm = genom = vinom = bafm = llrm = NULL
            genom = apply(cbind(MM[pid, ], SS[pid, ], hint[pid]), 1, function(x, 
                trans, doCNV, n) {
                b = genotype(x[1:n], x[(n + 1):(2 * n)], x[(2 * n + 1):(2 * n + 1)], 
                  trans, doCNV)
                c(b$geno, b$vino, b$conf, b$baf, b$llr)
            }, trans, doCNV, n)
            genom = t(genom)
            conf[pid, ] = genom[, c((2 * n + 1):(3 * n))]
            vino[pid, ] = genom[, c((n + 1):(2 * n))]
            geno[pid, ] = genom[, c(1:n)]
            if (doCNV) {
                baf[pid, ] = genom[, c((3 * n + 1):(4 * n))]
                llr[pid, ] = genom[, c((4 * n + 1):(5 * n))]
            }
        }
        
        rownames(geno) = rownames(conf) = rownames(vino) = rownames(MM)
        colnames(geno) = colnames(conf) = colnames(vino) = colnames(MM)
        xname = paste(savefiledir, "/genodataX", sep = "", collapse = "")
        save(geno, file = xname)
        xname = paste(savefiledir, "/vinodataX", sep = "", collapse = "")
        save(vino, file = xname)
        xname = paste(savefiledir, "/confdataX", sep = "", collapse = "")
        save(conf, file = xname)
        if (doCNV) {
            rownames(baf) = rownames(llr) = rownames(MM)
            colnames(baf) = colnames(llr) = colnames(MM)
            xname = paste(savefiledir, "/bafX", sep = "", collapse = "")
            save(baf, file = xname)
            xname = paste(savefiledir, "/llrX", sep = "", collapse = "")
            save(llr, file = xname)
        }
    }
    
    chrY = match(chr, "Y")
    if (!is.na(chrY)) {
        conf = geno = vino = baf = llr = matrix(-1, nn, n)
        n = sum(isMale)
        if (n > 1) {
            # for male two groups
            geno1 = apply(cbind(MM[, isMale], SS[, isMale]), 1, function(x, trans, doCNV, n) {
                b = genotypeSexchr(x[1:n], x[(n + 1):(2 * n)], trans, doCNV)
                c(b$geno, b$vino, b$conf, b$baf, b$llr)
            }, trans, doCNV, n)
            geno1 = t(geno1)
            conf[, isMale] = geno1[, c((2 * n + 1):(3 * n))]
            vino[, isMale] = geno1[, c((n + 1):(2 * n))]
            geno[, isMale] = geno1[, c(1:n)]
            if (doCNV) {
                baf[, isMale] = geno1[, c((3 * n + 1):(4 * n))]
                llr[, isMale] = geno1[, c((4 * n + 1):(5 * n))]
            }
        }
        rownames(geno) = rownames(conf) = rownames(vino) = rownames(MM)
        colnames(geno) = colnames(conf) = colnames(vino) = colnames(MM)
        xname = paste(savefiledir, "/genodataY", sep = "", collapse = "")
        save(geno, file = xname)
        xname = paste(savefiledir, "/vinodataY", sep = "", collapse = "")
        save(vino, file = xname)
        xname = paste(savefiledir, "/confdataY", sep = "", collapse = "")
        save(conf, file = xname)
        if (doCNV) {
            rownames(baf) = rownames(llr) = rownames(MM)
            colnames(baf) = colnames(llr) = colnames(MM)
            xname = paste(savefiledir, "/bafY", sep = "", collapse = "")
            save(baf, file = xname)
            xname = paste(savefiledir, "/llrY", sep = "", collapse = "")
            save(llr, file = xname)
        }
    }
    
    chrM = match(chr, "M")
    if (!is.na(chrM)) {
        conf = geno = vino = baf = llr = matrix(0, nn, n)
        geno1 = apply(cbind(MM, SS), 1, function(x, trans, doCNV, n) {
            b = genotypeSexchr(x[1:n], x[(n + 1):(2 * n)], trans, doCNV)
            c(b$geno, b$vino, b$conf, b$baf, b$llr)
        }, trans, doCNV, n)
        geno1 = t(geno1)
        conf = geno1[, c((2 * n + 1):(3 * n))]
        vino = geno1[, c((n + 1):(2 * n))]
        geno = geno1[, c(1:n)]
        
        rownames(geno) = rownames(conf) = rownames(vino) = rownames(MM)
        colnames(geno) = colnames(conf) = colnames(vino) = colnames(MM)
        xname = paste(savefiledir, "/genodataM", sep = "", collapse = "")
        save(geno, file = xname)
        xname = paste(savefiledir, "/vinodataM", sep = "", collapse = "")
        save(vino, file = xname)
        xname = paste(savefiledir, "/confdataM", sep = "", collapse = "")
        save(conf, file = xname)
        if (doCNV) {
            baf = geno1[, c((3 * n + 1):(4 * n))]
            llr = geno1[, c((4 * n + 1):(5 * n))]
            rownames(baf) = rownames(llr) = rownames(MM)
            colnames(baf) = colnames(llr) = colnames(MM)
            xname = paste(savefiledir, "/bafM", sep = "", collapse = "")
            save(baf, file = xname)
            xname = paste(savefiledir, "/llrM", sep = "", collapse = "")
            save(llr, file = xname)
        }
    }
} 
