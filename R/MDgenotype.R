#########################################################################
#
# MouseDivGenotye.R
#
# Part of the MouseDivGenotype package
#
# This is the main function to genotype the Mouse Diversity Array.
#
#########################################################################
MouseDivGenotype = function(allid, ABid, chrid, CGFLcorrection = NULL, 
    reference = NULL, hint = NULL, trans = c("CCStrans", "MAtrans"), celnamefile = NULL, 
    mchr = c(1:19, "X", "Y", "M"), celfiledir, outfiledir, subset = FALSE, doCNV = FALSE, 
    exon1info = NULL, exon2info = NULL, exonoutfiledir, cnvoutfiledir, verbose = FALSE) {
    library("affyio")
    library("preprocessCore")
    library("cluster")
    #library("time")
    
    trans = match.arg(trans)
    if (missing(celfiledir)) 
        celfiledir = getwd()
    if (missing(outfiledir)) 
        outfiledir = getwd()
    if (missing(exonoutfiledir)) 
        exonoutfiledir = getwd()
    if (missing(cnvoutfiledir)) 
        cnvoutfiledir = getwd()
    if (missing(allid) | missing(ABid)) 
        stop("No CDF file information")
    SNPname = ABid$SNPname
    Aid = ABid$allAid
    Bid = ABid$allBid
    rm(ABid)
    mpos = chrid$mpos
    if (subset) 
        chrid = chrid$schrid
    else
        chrid = chrid$chrid
    
    setwd(celfiledir)
    gender = isMale = NULL
    if (length(celnamefile) == 0) {
        tmp = dir()
        isCel = apply(matrix(tmp, ncol = 1), 1, function(x) length(grep(".CEL", x)) > 0)
        filenames = tmp[isCel]
    }
    else {
        filenames = read.table(celnamefile, header = TRUE)
        if (ncol(filenames) > 1) 
            gender = filenames[, 2]
        filenames = filenames[, 1]
    }
    nfile = length(filenames)
    
    # read individual files
    iiy = !is.na(match("Y", mchr))
    iix = !is.na(match("X", mchr))
    mchr1 = mchr
    if (iiy | iix) 
        mchr1 = unique(c(mchr, "X", "Y"))
    
    # iter through CEL files
    for (i in 1:nfile) {
        if(verbose) cat("Reading and normalizing CEL file: ", filenames[i], "\n", sep="")
        
        y = as.matrix(read.celfile(as.character(filenames[i]), intensity.means.only = TRUE)[["INTENSITY"]][["MEAN"]][allid])
        y = log2(y)
        if (length(CGFLcorrection) > 0)
            # C+G and fragment length correction y
            y = y + CGFLcorrection
        if (length(reference) > 0) 
            y <- normalize.quantiles.use.target(y, target = reference)
        
        # separate a alleles from be alleles
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
        if (is.element(trans, "MAtrans")) {
            # then prior??
            M = allAint - allBint
            S = (allAint + allBint)/2
        }
        
        # divide CEL file up into chromosome pieces
        for (chri in mchr1) {
            currRows <- which(chrid == chri)
            if(length(currRows) == 0)
            {
                stop("Failed to find any probes on chromosome ", chri);
            }
            
            xname2 = paste(outfiledir, "/", gsub(".CEL", "CHR", filenames[i]), chri, sep = "", collapse = "")
            if(verbose)
            {
                cat("Saving ", length(currRows),
                    " probes in chromosome ", chri, " to ", xname2, "\n", sep="")
            }
            
            MM1 = M[currRows, ]
            SS1 = S[currRows, ]
            save(MM1, SS1, file = xname2)
        }
    }
    
    # combine all sampels and genotype by chromosomes
    autoint = NULL
    mchr1 = c(1:19, "M")
    ii = match(mchr1, mchr)
    mchr1 = mchr1[!is.na(ii)]
    for (chri in mchr1) {
        if(verbose) cat("geno/vinotyping chromosome ", chri, "\n", sep="")
        
        #startTime <- getTime()
    
        # paste the chromosomes together for genotyping
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
        xname = paste(outfiledir, "/rawdataSSchr", chri, sep = "", collapse = "")
        save(SS, file = xname)
        isMale = rep(TRUE, ncol(MM))
        if (length(hint) == 0) 
            hint1 = NULL
        else
            hint1 = hint[[chri]]
        #cat("time it took us to get to genotypethis\n")
        #timeReport(startTime)
        
        #startTime <- getTime()
        genotypethis(outfiledir, MM, SS, hint1, isMale, trans, chri, doCNV)
        #cat("time it took us to complete genotypethis\n")
        #timeReport(startTime)
    }
    
    if (iix | iiy) {
        MM = SS = NULL
        for (i in 1:nfile) {
            xname1 = paste(outfiledir, "/", gsub(".CEL", "CHR", filenames[i]), "Y", sep = "", collapse = "")
            # origa,origb
            load(xname1)
            MM = cbind(MM, MM1)
            SS = cbind(SS, SS1)
            rm(MM1, SS1)
            file.remove(xname1)
        }
        intY = apply(SS, 2, mean)
        if (iiy) {
            xname = paste(outfiledir, "/rawdataMMchrY", sep = "", collapse = "")
            colnames(MM) = colnames(SS) = filenames
            save(MM, file = xname)
            xname = paste(outfiledir, "/rawdataSSchrY", sep = "", collapse = "")
            save(SS, file = xname)
            MMy = MM
            SSy = SS
        }
        MM = SS = NULL
        for (i in 1:nfile) {
            xname1 = paste(outfiledir, "/", gsub(".CEL", "CHR", filenames[i]), "X", sep = "", collapse = "")
            load(xname1)
            MM = cbind(MM, MM1)
            SS = cbind(SS, SS1)
            rm(MM1, SS1)
            file.remove(xname1)
        }
        intX = apply(SS, 2, mean)
        if (iix) {
            xname = paste(outfiledir, "/rawdataMMchrX", sep = "", collapse = "")
            colnames(MM) = colnames(SS) = filenames
            save(MM, file = xname)
            xname = paste(outfiledir, "/rawdataSSchrX", sep = "", collapse = "")
            save(SS, file = xname)
        }
        # compute the gender
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
        if (iix) {
            if (length(hint) == 0) 
                hint1 = NULL
            else hint1 = hint[[20]]
            genotypethis(outfiledir, MM, SS, hint1, isMale, trans, "X", doCNV)
        }
        if (iiy) 
            genotypethis(outfiledir, MMy, SSy, NULL, isMale, trans, "Y", doCNV)
    }
    if (doCNV) {
        penCNVinput(chrid, mpos, exon1info, exon2info, celfiledir, filenames, outfiledir, exonoutfiledir, cnvoutfiledir, mchr)
    }
}

 
