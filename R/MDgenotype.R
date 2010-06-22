#########################################################################
#
# MouseDivGenotye.R
#
# Part of the MouseDivGenotype package
#
# This is the main function to genotype the Mouse Diversity Array.
#
#########################################################################

#TODO REMOVE ME KEITH!!!
library("snow")

library("affyio")
library("preprocessCore")
library("cluster")
#library("time")

MouseDivGenotype = function(allid, ABid, chrid, CGFLcorrection = NULL, 
    reference = NULL, hint = NULL, trans = c("CCStrans", "MAtrans"), celnamefile = NULL, 
    mchr = c(1:19, "X", "Y", "M"), celfiledir, outfiledir, subset = FALSE, doCNV = FALSE, 
    exon1info = NULL, exon2info = NULL, exonoutfiledir, cnvoutfiledir, verbose = FALSE,
    cluster = NULL, probesetChunkSize=1000) {
    
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
        filenames = read.delim(celnamefile, header = TRUE)
        if (ncol(filenames) > 1) 
            gender = filenames[, 2]
        filenames = filenames[, 1]
    }
    nfile = length(filenames)
    
    # make sure that the chromosome vector is not numeric
    mchr <- as.character(mchr)
    
    # we need both X and Y to determine sex so if the caller asks for either
    # then we need to read both
    iiy = !is.na(match("Y", mchr))
    iix = !is.na(match("X", mchr))
    mchr1 = mchr
    if (iiy | iix) 
        mchr1 = unique(c(mchr, "X", "Y"))
    
    # figure out how many chunks there will be per chromosome
    chrChunks <- list()
    chrCount <- length(mchr1)
    for(i in 1 : chrCount)
    {
        currChr <- mchr1[[i]]
        chrProbesetCount <- sum(chrid == currChr)
        chrChunks[[currChr]] <- chunkIndices(chrProbesetCount, probesetChunkSize)
    }
    
    # loop through all the CELL files. We need to read and normalize the CEL
    # file data and save it to file in chunks no bigger than probesetChunkSize
    for(celfile in filenames)
    {
        # wrapping this up as a local function def
        # allows us to avoid the normalization work unless it's necessary
        makeMsList <- function()
        {
            normalizeCelFileByChr(
                celfile,
                probesetChunkSize,
                verbose,
                allid,
                Aid,
                Bid,
                CGFLcorrection,
                reference,
                SNPname,
                trans,
                mchr1,
                chrid)
        }
        msList <- NULL
        
        for(i in 1 : chrCount)
        {
            currChr <- mchr1[[i]]
            
            for(chunkIndex in 1 : length(chrChunks[[currChr]]))
            {
                chunkFile <- chunkFileName(outfiledir, celfile, currChr, probesetChunkSize, chunkIndex)
                if(!file.exists(chunkFile))
                {
                    # if we haven't yet normalized the CEL file we'll have to
                    # do that now
                    if(is.null(msList))
                    {
                        msList <- makeMsList()
                    }
                    
                    cat("writing to ", chunkFile, "\n", sep="")
                    
                    chunk <- chrChunks[[currChr]][[chunkIndex]]
                    indices <- chunk$start : chunk$end
                    mChunk <- msList[[i]]$M[indices]
                    sChunk <- msList[[i]]$S[indices]
                    
                    save(mChunk, sChunk, file = chunkFile)
                }
            }
        }
        
        msList <- NULL
        
        firstCelFile <- FALSE
    }
    
    # combine all sampels and genotype by chromosomes
    autoint = NULL
    mchr1 = c(1:19, "M")
    ii = match(mchr1, mchr)
    mchr1 = mchr1[!is.na(ii)]
    for (chri in mchr1) {
        #for(chunk in chrChunks[[currChr]])
        for(chunkIndex in 1 : length(chrChunks[[chri]]))
        {
            chunk <- chrChunks[[chri]][[chunkIndex]]
            if(verbose) {
                cat("geno/vinotyping chromosome ", chri,
                    " from probeset #", chunk$start,
                    " to probeset #", chunk$end, "\n", sep="")
            }
            
            #startTime <- getTime()
        
            # paste the chromosomes together for genotyping
            MM = SS = NULL
            for (i in 1:nfile) {
                chunkFile <- chunkFileName(outfiledir, filenames[i], chri, probesetChunkSize, chunkIndex)
                load(chunkFile)
                MM = cbind(MM, mChunk)
                SS = cbind(SS, sChunk)
                rm(mChunk, sChunk)
            }
            
            # TODO reduce the number of samples we're taking based on the
            #      chunk size
            if (!is.na(match(chri, c(1:19)))) 
                autoint = c(autoint, mean(SS[sample(c(1:nrow(SS)), 100), ]))
            colnames(MM) <- filenames
            colnames(SS) <- filenames
            
            isMale = rep(TRUE, ncol(MM))
            
            chunkIndices <- chunk$start : chunk$end
            if (length(hint) == 0)
                hint1 = NULL
            else
                hint1 = hint[[chri]][chunkIndices]
            
            #cat("time it took us to get to genotypethis\n")
            #timeReport(startTime)
            
            #startTime <- getTime()
            genos <- genotypeAutosomeChunk(MM, SS, hint1, isMale, trans, doCNV)
            
            #cat("time it took us to complete genotypethis\n")
            #timeReport(startTime)
        }
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
        pennCNVinput(chrid, mpos, exon1info, exon2info, celfiledir, filenames, outfiledir, exonoutfiledir, cnvoutfiledir, mchr)
    }
}
