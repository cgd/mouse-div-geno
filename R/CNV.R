# TODO make sure it's well documented that the SNPs in "genotypes must match the
#      order of the SNPs in snpProbeInfo (otherwise we should do something
#      allowing for SNP IDs to be used... actually this is what we should do)
buildPennCNVInputFiles <- function(
    outdir = getwd(), allowOverwrite = FALSE,
    genotypes, snpProbeInfo, snpInfo, snpReferenceDistribution = NULL,
    invariantProbeInfo, invariantProbesetInfo, invariantReferenceDistribution = NULL,
    transformMethod = c("CCStrans", "MAtrans"), celFiles = getwd(),
    chromosomes = c(1:19, "X", "Y", "M"), cacheDir = tempdir(),
    retainCache = FALSE, verbose = FALSE,
    probesetChunkSize = 1000)
{
    snpCount <- nrow(snpInfo)
    sampleCount <- ncol(genotypes)
    
    transformMethod <- match.arg(transformMethod)
    
    # validate snp info parameters
    if(!inherits(snpProbeInfo, "data.frame") ||
       !all(c("probeIndex", "isAAllele", "snpId") %in% names(snpProbeInfo)))
    {
        stop("You must supply a \"snpProbeInfo\" data frame parameter which has ",
            "at a minimum the \"probeIndex\", \"isAAllele\" and \"snpId\" ",
            "components. Please see the help documentation for more details.")
    }
    
    if(!inherits(snpInfo, "data.frame") ||
       !all(c("snpId", "chrId", "positionBp") %in% names(snpInfo)))
    {
        stop("You must supply a \"snpInfo\" data frame parameter which has ",
            "at a minimum the \"snpId\", \"chrId\" and \"positionBp\"",
            "components. Please see the help documentation for more details.")
    }
    
    if(!is.null(snpInfo$isInPAR))
    {
        parChrs <- unique(snpInfo$chrId[snpInfo$isInPar])
        if(!all(parChrs == "X"))
        {
            stop("snpInfo$isInPar should only ever be TRUE on the \"X\" ",
                "chromosome, but TRUE isInPar values were found on chromosomes: ",
                paste(parChrs, collapse=", "))
        }
    }
    
    # validate invariant parameters
    if(!inherits(invariantProbeInfo, "data.frame") ||
        !all(c("probeIndex", "probesetId") %in% names(invariantProbeInfo)))
    {
        stop("You must supply a \"invariantProbeInfo\" data frame parameter which has ",
            "at a minimum the \"probeIndex\", and \"probesetId\" ",
            "components. Please see the help documentation for more details.")
    }
    
    if(!inherits(invariantProbesetInfo, "data.frame") ||
        !all(c("probesetId", "chrId", "positionBp") %in% names(invariantProbesetInfo)))
    {
        stop("You must supply a \"invariantProbesetInfo\" data frame parameter which has ",
            "at a minimum the \"probesetId\", \"chrId\" and \"positionBp\"",
            "components. Please see the help documentation for more details.")
    }
    
    # if the user passes us a list (or data frame) rather than a vector of file
    # names we should attempt to pull out a fileName component
    if(is.list(celFiles))
    {
        if(!("fileName" %in% names(celFiles)))
        {
            stop("failed to find \"fileName\" component in the \"celFiles\" list")
        }
        
        celFiles <- celFiles$fileName
    }
    
    if(length(celFiles) != sampleCount)
    {
        stop("the number of CEL files doesn't match the number of columns in ",
            "\"genotypes\" (", sampleCount, "). The expanded CEL files are: ",
            paste(celFiles, sep = ", "))
    }
    
    genoRownames <- rownames(genotypes)
    if(is.null(genoRownames))
    {
        if(nrow(genotypes) != snpCount)
        {
            stop("the number of rows in the \"genotypes\" parameter must match ",
                 "the number of rows in the \"snpInfo\" parameter. Alternatively ",
                 "you can set \"rownames(genotypes)\" to match up with the ",
                 "corresponding \"snpInfo$snpId\" values")
        }
    }
    else
    {
        # make sure that the genotype rows all match up OK
        genoIndexMapping <- match(genoRownames, snpInfo$snpId)
        if(any(is.na(genoIndexMapping)))
        {
            stop("all of the \"rownames(genotypes)\" values must match up with ",
                 "a SNP ID in \"snpInfo$snpId\". Alternatively you can set ",
                 "\"rownames(genotypes)\" to NULL which implies that genotypes ",
                 "will be matched to snpInfo rows on index alone")
        }
        
        # subset the SNP info to just those SNPs that we have genotypes for
        snpInfo <- snpInfo[genoIndexMapping, ]
        rm(genoIndexMapping)
    }
    rm(genoRownames)
    
    # make sure that the chromosome vector is not numeric
    # TODO is toupper the right thing to do here? (I do it in MDgenotype.R too)
    chromosomes <- toupper(as.character(chromosomes))
    
    snpChromosomes <- unique(snpInfo$chrId)
    if(!all(chromosomes %in% snpChromosomes))
    {
        warning(
            "SNP data for the following requested chromosomes are not available: ",
            paste(setdiff(chromosomes, snpChromosomes), collapse = ", "),
            ". These chromosomes will be skipped.")
    }
    snpChromosomes <- intersect(chromosomes, snpChromosomes)
    
    invariantChromosomes <- unique(invariantProbesetInfo$chrId)
    if(!all(chromosomes %in% invariantChromosomes))
    {
        warning(
            "Invariant data for the following requested chromosomes are not available: ",
            paste(setdiff(chromosomes, invariantChromosomes), collapse = ", "),
            ". These chromosomes will be skipped.")
    }
    invariantChromosomes <- unique(chromosomes, invariantChromosomes)
    rm(chromosomes)
    
    if(verbose)
    {
        cat("processing CEL files for SNP probes\n")
    }
    
    # figure out how many chunks there will be per chromosome
    chrChunks <- list()
    for(currChr in snpChromosomes)
    {
        chrProbesetCount <- sum(snpInfo$chrId == currChr)
        chrChunks[[currChr]] <- chunkIndices(chrProbesetCount, probesetChunkSize)
    }
    
    for(celfile in celFiles)
    {
        # wrapping this up as a local function def allows us to avoid doing the
        # normalization work unless it's necessary
        makeMsList <- function()
        {
            normalizeCelFileByChr(
                celfile,
                verbose,
                snpProbeInfo,
                snpInfo,
                snpChromosomes,
                snpReferenceDistribution,
                transformMethod)
        }
        msList <- NULL
        
        for(currChr in snpChromosomes)
        {
            for(chunkIndex in 1 : length(chrChunks[[currChr]]))
            {
                chunkFile <- chunkFileName(cacheDir, "snp", celfile, currChr, probesetChunkSize, chunkIndex)
                chunkFileAlreadyExists <- file.exists(chunkFile)
                
                if(!chunkFileAlreadyExists)
                {
                    # if we haven't yet normalized the CEL file we'll have to
                    # do that now
                    if(is.null(msList))
                    {
                        msList <- makeMsList()
                    }
                    
                    chunk <- chrChunks[[currChr]][[chunkIndex]]
                    probesetIndices <- chunk$start : chunk$end
                    mChunk <- msList[[currChr]]$M[probesetIndices]
                    sChunk <- msList[[currChr]]$S[probesetIndices]
                    
                    save(mChunk, sChunk, file = chunkFile)
                }
            }
        }
        
        rm(msList)
    }
    
    lrrAndBafBaseNames <- paste(fileBaseWithoutExtension(celFiles), ".txt", sep = "")
    lrrAndBafOutputFiles <- file.path(outdir, lrrAndBafBaseNames)
    if(!allowOverwrite && any(file.exists(lrrAndBafOutputFiles)))
    {
        stop("refusing to overwrite the following files: ",
             paste(lrrAndBafOutputFiles[file.exists(lrrAndBafOutputFiles)], collapse = ", "),
             ". In order to procede you can either set allowOverwrite to TRUE or ",
             "you can manually delete the files.")
    }
    
    pfbOutputFile <- file.path(outdir, "pfbdata.txt")
    if(!allowOverwrite && file.exists(pfbOutputFile))
    {
        stop("refusing to overwrite : ", pfbOutputFile,
            ". In order to procede you can either set allowOverwrite to TRUE or ",
            "you can manually delete the files.")
    }
    
    lrrAndBafConnections <- lapply(as.list(lrrAndBafOutputFiles), file, "wt")
    pfbConnection <- file(pfbOutputFile, "wt")
    
    # TODO write the headers!!!!!!!!!!!!!!!!!!!!!!!
    
    if(verbose)
    {
        cat("generating PennCNV input for SNP probes\n")
    }
    
    # append the SNP-based LRR/BAF values
    for(currChr in snpChromosomes)
    {
        chrIndices <- which(snpInfo$chrId == currChr)
        chrGenos <- genotypes[chrIndices, ]
        chrSnpInfo <- snpInfo[chrIndices, ]
        for(chunkIndex in 1 : length(chrChunks[[currChr]]))
        {
            chunk <- chrChunks[[currChr]][[chunkIndex]]
            chunkIndices <- chunk$start : chunk$end
            chunkSize <- 1 + chunk$end - chunk$start
            
            # pre-allocate M and S matrices
            mMatrix <- matrix(0, nrow = chunkSize, ncol = length(celFiles))
            sMatrix <- matrix(0, nrow = chunkSize, ncol = length(celFiles))
            
            # stitch together the M and S chunks
            for(fileIndex in 1 : length(celFiles))
            {
                celfile <- celFiles[fileIndex]
                
                # loads mChunk and sChunk into scope
                chunkFile <- chunkFileName(cacheDir, "snp", celfile, currChr, probesetChunkSize, chunkIndex)
                load(chunkFile)
                
                mMatrix[, fileIndex] <- mChunk
                sMatrix[, fileIndex] <- sChunk
            }
            
            appendToPennCNVForSNPs(
                snpInfo = chrSnpInfo[chunkIndices, ],
                intensityConts = mMatrix,
                intensityAvgs = sMatrix,
                genotypes = chrGenos[chunkIndices, ],
                lrrAndBafConnections = lrrAndBafConnections,
                pfbConnection = pfbConnection)
        }
    }
    
    if(verbose)
    {
        cat("processing CEL files for invariant probes\n")
    }
    
    invariantChrChunks <- list()
    for(currChr in invariantChromosomes)
    {
        chrProbesetCount <- sum(invariantProbesetInfo$chrId == currChr)
        invariantChrChunks[[currChr]] <- chunkIndices(chrProbesetCount, probesetChunkSize)
    }

    # append the invariant LRR/BAF values
    for(celfile in celFiles)
    {
        makeNormInvariantList <- function()
        {
            #TODO remember to add invariantProbesetInfo !!!!!!!!!!!!!!!!!!!!!!!!!!!
            normalizeCelFileForInvariants(
                celfile,
                verbose,
                invariantProbeInfo,
                invariantProbesetInfo,
                invariantChromosomes,
                invariantReferenceDistribution)
        }
        normInvariantList <- NULL
        
        for(currChr in invariantChromosomes)
        {
            for(chunkIndex in 1 : length(invariantChrChunks[[currChr]]))
            {
                chunkFile <- chunkFileName(cacheDir, "invariant", celfile, currChr, probesetChunkSize, chunkIndex)
                chunkFileAlreadyExists <- file.exists(chunkFile)
                
                if(!chunkFileAlreadyExists)
                {
                    # if we haven't yet normalized the CEL file we'll have to
                    # do that now
                    if(is.null(normInvariantList))
                    {
                        normInvariantList <- makeNormInvariantList()
                    }
                    
                    chunk <- invariantChrChunks[[currChr]][[chunkIndex]]
                    probesetIndices <- chunk$start : chunk$end
                    intensityChunk <- normInvariantList[[currChr]][probesetIndices]
                    
                    save(intensityChunk = intensityChunk, file = chunkFile)
                }
            }
        }
        
        rm(normInvariantList)
    }
    
    if(verbose)
    {
        cat("generating PennCNV input for invariant probes\n")
    }
    
    for(currChr in invariantChromosomes)
    {
        chrInvariantProbesetInfo <- invariantProbesetInfo[invariantProbesetInfo$chrId == currChr, ]
        for(chunkIndex in 1 : length(invariantChrChunks[[currChr]]))
        {
            chunk <- invariantChrChunks[[currChr]][[chunkIndex]]
            chunkIndices <- chunk$start : chunk$end
            chunkSize <- 1 + chunk$end - chunk$start
            
            # pre-allocate intensity matrix
            intensityMatrix <- matrix(0, nrow = chunkSize, ncol = length(celFiles))
            
            # stitch together the intensity chunks
            for(fileIndex in 1 : length(celFiles))
            {
                celfile <- celFiles[fileIndex]
                
                # loads intensityChunk into scope
                chunkFile <- chunkFileName(cacheDir, "invariant", celfile, currChr, probesetChunkSize, chunkIndex)
                load(chunkFile)
                
                intensityMatrix[, fileIndex] <- intensityChunk
            }
            
            appendToPennCNVForInvariants(
                probesetInfo = chrInvariantProbesetInfo[chunkIndices, ],
                probesetIntensities = intensityMatrix,
                lrrAndBafConnections = lrrAndBafConnections,
                pfbConnection = pfbConnection)
        }
    }
    
    # close all of the open connections
    for(con in lrrAndBafConnections)
    {
        close(con)
    }
    close(pfbConnection)
}

# append PennCNV data for the given LRR/BAF files and the PFB file using
# invariants (thes files should already have a header in place since this
# function will not create one)
# PARAMETERS:
#   probesetInfo:
#       a data frame with a row per probeset. this parameter should have the
#       following components: probesetId, chrId, positionBp
#   probesetIntensities:
#       the summerized intensities per-probeset
#   lrrAndBafConnection:
#       vector of connections to append to for LRR and BAF data (one connection
#       per sample)
#   pfbConnection:
#       the connection to use for the PBF file.
appendToPennCNVForInvariants <- function(
    probesetInfo,
    probesetIntensities,
    lrrAndBafConnections,
    pfbConnection)
{
    probesetCount <- nrow(probesetInfo)
    sampleCount <- length(lrrAndBafConnections)
    
    if(!all(c("probesetId", "chrId", "positionBp") %in% names(probesetInfo)))
    {
        stop("the \"probesetInfo\" parameter must contain the following components: ",
             "probesetId, chrId, positionBp")
    }
    
    if(!all(sapply(lrrAndBafConnections, inherits, "connection")) ||
       !inherits(pfbConnection, "connection"))
    {
        stop("both lrrAndBafConnections and pfbConnection should inherit from ",
             "the \"connection\" class")
    }
    
    # set all BAFs to 2 indicating that it's not a polymorphic probeset
    bafs <- rep(2, probesetCount)
    lrrs <- log2(probesetIntensities / apply(probesetIntensities, 1, mean))
    for(sampleIndex in 1 : sampleCount)
    {
        write.table(
            data.frame(probesetInfo$probesetId, lrrs[, sampleIndex], bafs),
            file = lrrAndBafConnections[[sampleIndex]],
            sep = "\t",
            row.names = FALSE,
            col.names = FALSE,
            qmethod = "double")
    }
    
    # write the PFB (Population frequency of B allele) file using mean BAF in
    # the PBF column
    write.table(
        data.frame(probesetInfo$probesetId, probesetInfo$chrId, probesetInfo$positionBp, bafs),
        file = pfbConnection,
        sep = "\t",
        row.names = FALSE,
        col.names = FALSE,
        qmethod = "double")
}

# append PennCNV data for the given LRR/BAF files and the PFB file using SNP
# data (thes files should already have a header in place since this function
# will not create one)
# PARAMETERS:
#   snpInfo:
#       a data frame with a row per SNP. this parameter should have the
#       following components: snpId, chrId, positionBp
#   intensityConts:
#       a matrix of intensity contrasts which has a column per sample and a row
#       per SNP
#   intensityAvgs:
#       a matrix of intensity averages which has a column per sample and a row
#       per SNP
#   genotypes:
#       a matrix of genotype codes which has a column per sample and a row
#       per SNP
#   lrrAndBafConnection:
#       vector of connections to append to for LRR and BAF data (one connection
#       per sample)
#   pfbConnection:
#       the connection to use for the PBF file.
appendToPennCNVForSNPs <- function(
    snpInfo,
    intensityConts,
    intensityAvgs,
    genotypes,
    lrrAndBafConnections,
    pfbConnection)
{
    if(!all(sapply(lrrAndBafConnections, inherits, "connection")) ||
       !inherits(pfbConnection, "connection"))
    {
        stop("both lrrAndBafConnections and pfbConnection should inherit from ",
            "the \"connection\" class")
    }
    
    snpCount <- nrow(snpInfo)
    sampleCount <- length(lrrAndBafConnections)
    
    # preallocate BAF and LRR matrices for speed then calculate BAF and LRR for
    # each SNP
    bafs <- matrix(0.0, nrow = snpCount, ncol = sampleCount)
    lrrs <- matrix(0.0, nrow = snpCount, ncol = sampleCount)
    for(snpIndex in 1 : snpCount)
    {
        bafAndLrr <- calcLRRAndBAF(
            intensityConts[snpIndex, ],
            intensityAvgs[snpIndex, ],
            genotypes[snpIndex, ])
        bafs[snpIndex, ] <- bafAndLrr$BAF
        lrrs[snpIndex, ] <- bafAndLrr$LRR
    }
    
    # write the LRR/BAF files per-sample
    for(sampleIndex in 1 : sampleCount)
    {
        write.table(
            data.frame(snpInfo$snpId, lrrs[, sampleIndex], bafs[, sampleIndex]),
            file = lrrAndBafConnections[[sampleIndex]],
            sep = "\t",
            row.names = FALSE,
            col.names = FALSE,
            qmethod = "double")
    }
    
    # write the PFB (Population frequency of B allele) file using mean BAF in
    # the PBF column
    write.table(
        data.frame(snpInfo$snpId, snpInfo$chrId, snpInfo$positionBp, apply(bafs, 1, mean)),
        file = pfbConnection,
        sep = "\t",
        row.names = FALSE,
        col.names = FALSE,
        qmethod = "double")
}

# calculates LRR and BAF values for a single SNP position
# PARAMETERS:
#   intensityConts:
#       per-sample intensity contrasts (A allele vs B allele)
#   intensityAvgs:
#       per-sample intensity averages (A averaged with B)
#   genos:
#       per-sample genotypes 1 = AA, 2 = AB, 3 = BB
# RETURNS:
#   A data.frame object with B-Allele Frequency (BAF) and Log R Ratio (LRR)
#   components. The row count of this dataframe will equal length(genos)
calcLRRAndBAF <- function(intensityConts, intensityAvgs, genos)
{
    sampleCount <- length(genos)
    uniqueGenos <- sort(unique(genos))
    uniqueGenoCount <- length(uniqueGenos)
    BAF <- rep(0, sampleCount)
    LRR <- rep(0, sampleCount)
    
    medianContPerGeno <- numeric(uniqueGenoCount)
    medianAvgsPerGeno <- numeric(uniqueGenoCount)
    contVariancePerGeno <- numeric(uniqueGenoCount)
    for(genoIndex in 1 : uniqueGenoCount)
    {
        genoCode <- uniqueGenos[genoIndex]
        currIndices <- which(genos == genoCode)
        
        currIntensityConts <- intensityConts[currIndices]
        currIntensityAvgs <- intensityAvgs[currIndices]
        
        if(length(currIndices) == 1)
        {
            medianContPerGeno[genoIndex] <- currIntensityConts
            medianAvgsPerGeno[genoIndex] <- currIntensityAvgs
            
            # there's no valid variance for the single-genotype case. So,
            # we'll set it to NA for now and fix it later.
            contVariancePerGeno[genoIndex] <- NA
        }
        else
        {
            medianContPerGeno[genoIndex] <- median(currIntensityConts)
            medianAvgsPerGeno[genoIndex] <- median(currIntensityAvgs)
            contVariancePerGeno[genoIndex] <- bivar(currIntensityConts)
        }
    }
    
    naVars <- is.na(contVariancePerGeno)
    naVarIndices <- which(naVars)
    validVarIndices <- which(!naVars)
    
    # the genotypes with only a single sample will be given the average
    # variance of the other genotypes
    averagedVariances <- sum(contVariancePerGeno[validVarIndices]) / length(validVarIndices)
    contVariancePerGeno[naVarIndices] <- averagedVariances
    
    # shrink the variances toward the averaged variances
    for(genoIndex in validVarIndices)
    {
        contVariancePerGeno[genoIndex] <-
            0.5 * (averagedVariances + contVariancePerGeno[genoIndex])
    }
    contStdDevs <- sqrt(contVariancePerGeno)
    
    # all three genotypes should be represented
    if (uniqueGenoCount == 3) {
        aaIndex <- 1
        abIndex <- 2
        bbIndex <- 3
        
        # for the contrasts that are >= the AA median
        # ----BB----AB----AAxxxx
        tmp <- intensityConts >= medianContPerGeno[aaIndex]
        if (any(tmp)) {
            BAF[tmp] <- 0
            LRR[tmp] <- log2(intensityAvgs[tmp]/medianAvgsPerGeno[aaIndex])
        }
        
        # < AA median and > AB meadian:
        # ----BB----AB-xx-AA----
        tmp <- intensityConts < medianContPerGeno[aaIndex] & intensityConts > medianContPerGeno[abIndex]
        if (any(tmp)) {
            k1 <- (medianContPerGeno[aaIndex] - intensityConts[tmp])/contStdDevs[aaIndex]
            k2 <- (intensityConts[tmp] - medianContPerGeno[abIndex])/contStdDevs[abIndex]
            BAF[tmp] <- 0.5 * k1/(k1 + k2)
            LRR[tmp] <- log2(intensityAvgs[tmp]/((k2 * medianAvgsPerGeno[aaIndex] + k1 * medianAvgsPerGeno[abIndex])/(k1 + k2)))
        }
        
        # <= the AB median and > BB median
        # ----BB-xxxAB----AA----
        tmp <- intensityConts <= medianContPerGeno[abIndex] & intensityConts > medianContPerGeno[bbIndex]
        if (any(tmp)) {
            k1 <- (medianContPerGeno[abIndex] - intensityConts[tmp])/contStdDevs[abIndex]
            k2 <- (intensityConts[tmp] - medianContPerGeno[bbIndex])/contStdDevs[bbIndex]
            BAF[tmp] <- 0.5 + 0.5 * k1/(k1 + k2)
            LRR[tmp] <- log2(intensityAvgs[tmp]/((k2 * medianAvgsPerGeno[abIndex] + k1 * medianAvgsPerGeno[bbIndex])/(k1 + k2)))
        }
        
        # <= BB median
        # xxxxBB----AB----AA----
        tmp <- intensityConts <= medianContPerGeno[bbIndex]
        if (any(tmp)) {
            BAF[tmp] <- 1
            LRR[tmp] <- log2(intensityAvgs[tmp]/medianAvgsPerGeno[bbIndex])
        }
    }
    
    # there are only two genotypes represented in the samples for this SNP
    else if (uniqueGenoCount == 2) {
        id <- sum(uniqueGenos)
        
        # the genotypes are AA and BB: c(1, 3)
        if (id == 4) {
            aaIndex <- 1
            bbIndex <- 2
            
            # >= AA median
            # ----BB----AB----AAxxxx
            tmp <- intensityConts >= medianContPerGeno[aaIndex]
            if (any(tmp)) {
                BAF[tmp] <- 0
                LRR[tmp] <- log2(intensityAvgs[tmp]/medianAvgsPerGeno[aaIndex])
            }
            
            # < AA median and > BB median
            # ----BB-xxxABxxx-AA----
            tmp <- intensityConts < medianContPerGeno[aaIndex] & intensityConts > medianContPerGeno[bbIndex]
            if (any(tmp)) {
                k1 <- (medianContPerGeno[aaIndex] - intensityConts[tmp])/contStdDevs[aaIndex]
                k3 <- (intensityConts[tmp] - medianContPerGeno[bbIndex])/contStdDevs[bbIndex]
                BAF[tmp] <- k1/(k1 + k3)
                LRR[tmp] <- log2(intensityAvgs[tmp]/((k3 * medianAvgsPerGeno[aaIndex] + k1 * medianAvgsPerGeno[bbIndex])/(k1 + k3)))
            }
            
            # <= BB median
            # xxxxBB----AB----AA----
            tmp <- intensityConts <= medianContPerGeno[bbIndex]
            if (any(tmp)) {
                BAF[tmp] <- 1
                LRR[tmp] <- log2(intensityAvgs[tmp]/medianAvgsPerGeno[bbIndex])
            }
        }
        
        # the genotypes are AB and BB: c(2, 3)
        else if (id == 5) {
            abIndex <- 1
            bbIndex <- 2
            
            # >= AB median
            # ----BB----ABxxxxAAxxxx
            tmp <- intensityConts >= medianContPerGeno[abIndex]
            if (any(tmp)) {
                BAF[tmp] <- 0.5
                LRR[tmp] <- log2(intensityAvgs[tmp]/medianAvgsPerGeno[abIndex])
            }
            
            # < AB median and > BB median
            # ----BB-xx-AB----AA----
            tmp <- intensityConts < medianContPerGeno[abIndex] & intensityConts > medianContPerGeno[bbIndex]
            if (any(tmp)) {
                k2 <- (medianContPerGeno[abIndex] - intensityConts[tmp])/contStdDevs[abIndex]
                k3 <- (intensityConts[tmp] - medianContPerGeno[bbIndex])/contStdDevs[bbIndex]
                BAF[tmp] <- 0.5 + 0.5 * k2/(k2 + k3)
                LRR[tmp] <- log2(intensityAvgs[tmp]/((k3 * medianAvgsPerGeno[abIndex] + k2 * medianAvgsPerGeno[bbIndex])/(k2 + k3)))
            }
            
            # <= BB median
            # xxxxBB----AB----AA----
            tmp <- intensityConts <= medianContPerGeno[bbIndex]
            if (any(tmp)) {
                BAF[tmp] <- 1
                LRR[tmp] <- log2(intensityAvgs[tmp]/medianAvgsPerGeno[bbIndex])
            }
        }
        
        # the genotypes are AA and AB: c(1, 2)
        else if (id == 3) {
            aaIndex <- 1
            abIndex <- 2
            
            # >= AA median
            # ----BB----AB----AAxxxx
            tmp <- intensityConts >= medianContPerGeno[aaIndex]
            if (any(tmp)) {
                BAF[tmp] <- 0
                LRR[tmp] <- log2(intensityAvgs[tmp]/medianAvgsPerGeno[aaIndex])
            }
            
            # < AA median and > AB median
            # ----BB----AB-xx-AA----
            tmp <- intensityConts < medianContPerGeno[aaIndex] & intensityConts > medianContPerGeno[abIndex]
            if (any(tmp)) {
                k1 <- (medianContPerGeno[aaIndex] - intensityConts[tmp])/contStdDevs[aaIndex]
                k2 <- (intensityConts[tmp] - medianContPerGeno[abIndex])/contStdDevs[abIndex]
                BAF[tmp] <- 0.5 * k1/(k1 + k2)
                LRR[tmp] <- log2(intensityAvgs[tmp]/((k1 * medianAvgsPerGeno[abIndex] + k2 * medianAvgsPerGeno[aaIndex])/(k1 + k2)))
            }
            
            # <= AB median
            # xxxxBBxxxxAB----AA----
            tmp <- intensityConts <= medianContPerGeno[abIndex]
            if (any(tmp)) {
                BAF[tmp] <- 0.5
                LRR[tmp] <- log2(intensityAvgs[tmp]/medianAvgsPerGeno[abIndex])
            }
        }
    }
    
    # if there is only a single genotype at this SNP for all samples
    else if (uniqueGenoCount == 1) {
        # the genotype for all samples is AA
        if (uniqueGenos == 1) {
            BAF <- rep(0, sampleCount)
            LRR <- log2(intensityAvgs/medianAvgsPerGeno[1])
        }
        
        # the genotype for all samples is AB
        else if (uniqueGenos == 2) {
            BAF <- rep(0.5, sampleCount)
            LRR <- log2(intensityAvgs/medianAvgsPerGeno[2])
        }
        
        # the genotype for all samples is BB
        else if (uniqueGenos == 3) {
            BAF <- rep(1, sampleCount)
            LRR <- log2(intensityAvgs/medianAvgsPerGeno[3])
        }
    }
    
    data.frame(BAF = BAF, LRR = LRR)
}

normalizeCelFileForInvariants <- function(
    celFileName,
    verbose,
    invariantProbeInfo,
    invariantProbesetInfo,
    allChr,
    referenceDistribution)
{
    if(verbose) cat("Reading and normalizing CEL file: ", celFileName, "\n", sep="")
    
    celData <- read.celfile(celFileName, intensity.means.only = TRUE)
    y <- log2(as.matrix(celData[["INTENSITY"]][["MEAN"]][invariantProbeInfo$probeIndex]))
    if (length(invariantProbeInfo$correction) > 0)
        # C+G and fragment length correction y
        y <- y + invariantProbeInfo$correction
    if (length(referenceDistribution) > 0)
        y <- normalize.quantiles.use.target(y, target = referenceDistribution)
    
    y <- subColSummarizeMedian(matrix(y, ncol = 1), invariantProbeInfo$probesetId)
    
    # now that the data is normalized by probeset,
    # organize the data into chromosome groups
    chrIntensityList <- list()
    for(chr in allChr)
    {
        chrIntensityList[[chr]] <- y[invariantProbesetInfo$chrId == chr]
    }
    
    chrIntensityList
}

# this is essentially using the same normalization as normalizeCelFileForInvariants
# except that it keeps track of corresponding "pos" also
simpleigp <- function(
    celfiledir,
    id,
    chrid,
    mpos,
    ename,
    CGFLcorrection,
    reference,
    filenames,
    mchr)
{
    setwd(celfiledir)
    
    y <- as.matrix(read.celfile(as.character(filenames), intensity.means.only = TRUE)[["INTENSITY"]][["MEAN"]][id])
    y <- log2(y)
    y <- y + CGFLcorrection
    y <- normalize.quantiles.use.target(y, target = reference)
    y <- subColSummarizeMedian(matrix(y, ncol = 1), ename)
    out <- list()
    for (chri in mchr) {
        pos <- mpos[chrid == chri]
        ey <- y[chrid == chri]
        out[[chri]] <- list(pos = pos, ey = ey)
    }
    out
}

# top-level function for running the simple CNV analysis
simpleCNV <- function(
    allid,
    ABid,
    chrid,
    CGFLcorrection,
    reference,
    exon1info, 
    exon2info,
    celfiledir,
    filenames,
    cnvoutfiledir,
    mchr = c(1:19),
    stname,
    refsample)
{
    if (missing(celfiledir)) 
        celfiledir <- getwd()
    if (missing(cnvoutfiledir)) 
        outfiledir <- getwd()
    if (missing(stname)) 
        stname <- filenames
    refid <- match(refsample, filenames)
    stname <- as.character(stname)
    simpleCNVdata(
        allid,
        ABid,
        chrid,
        CGFLcorrection,
        reference,
        exon1info,
        exon2info,
        celfiledir,
        filenames,
        cnvoutfiledir,
        mchr)
    simpleCNVsummary(cnvoutfiledir, filenames, refid, mchr, stname)
}

# this function processes the SNP and exon intensity data and then saves the
# results to file
simpleCNVdata <- function(
    allid,
    ABid,
    chrid,
    CGFLcorrection,
    reference,
    exon1info,
    exon2info,
    celfiledir,
    filenames,
    cnvoutfiledir,
    mchr)
{
    # some simple variable initialization
    setwd(celfiledir)
    SNPname <- ABid$SNPname
    Aid <- ABid$allAid
    Bid <- ABid$allBid
    rm(ABid)
    mpos <- chrid$mpos
    chrid <- chrid$chrid
    nfile <- length(filenames)
    
    for (i in 1:nfile) {
        
        # normalize the SNP probesets
        y <- as.matrix(read.celfile(as.character(filenames[i]), intensity.means.only = TRUE)[["INTENSITY"]][["MEAN"]][allid])
        y <- log2(y)
        if (length(CGFLcorrection) > 0) 
            y <- y + CGFLcorrection
        if (length(reference) > 0) 
            y <- normalize.quantiles.use.target(y, target = reference)
        oAint <- y[Aid, 1, drop = FALSE]
        oBint <- y[Bid, 1, drop = FALSE]
        allAint <- subColSummarizeMedian(matrix(oAint, ncol = 1), SNPname)
        allBint <- subColSummarizeMedian(matrix(oBint, ncol = 1), SNPname)
        
        # B will be the higher of the A or B allele (the allele that we think is "on")
        tmp <- names(allAint[allAint > allBint, ])
        g <- match(SNPname, tmp)
        B <- oAint
        B[is.na(g)] <- oBint[is.na(g)]
        B <- normalize.quantiles.use.target(B, target = exon1info$reference)
        B <- subColSummarizeMedian(matrix(B, ncol = 1), SNPname)
        
        # normalize the exon probesets and return the data as per-chromosome
        # lists (each chromosome list item has a pos component and an ey component)
        exon1 <- simpleigp(
            celfiledir,
            exon1info$exon1id,
            exon1info$chrid,
            exon1info$mpos,
            exon1info$exon1,
            exon1info$CGFLcorrection,
            exon1info$reference,
            filenames[i],
            mchr)
        exon2 <- simpleigp(
            celfiledir,
            exon2info$exon2id,
            exon2info$chrid,
            exon2info$mpos,
            exon2info$exon2,
            exon2info$CGFLcorrection,
            exon2info$reference,
            filenames[i],
            mchr)
        
        for (chri in mchr) {
            # for pos, id and inte append snp vals, exon1 vals and exon2 vals
            pos <- c(mpos[chrid == chri], exon1[[chri]]$pos, exon2[[chri]]$pos)
            id <- c(
                rep(1, sum(chrid == chri)),
                rep(2, length(exon1[[chri]]$pos)),
                rep(3, length(exon2[[chri]]$pos)))
            inte <- c(B[chrid == chri], exon1[[chri]]$ey, exon2[[chri]]$ey)
            
            # reorder the values according to thier posisions
            o <- order(pos)
            pos <- pos[o]
            id <- id[o]
            inte <- inte[o]
            
            # save intensities per CEL file, per chromosome
            xname2 <- paste(cnvoutfiledir, "/", gsub(".CEL", "Chr", filenames[i]), chri, sep = "", collapse = "")
            save(inte, file = xname2)
            
            # for the 1st CEL file only, save the position and ID data
            if (i == 1) {
                xname2 <- paste(cnvoutfiledir, "/pos", chri, sep = "", collapse = "")
                save(pos, id, file = xname2)
            }
        }
    }
}

# this function calculates the maximum likelihood CNV states using an HMM,
# then saves those states along with other statistics to file
simpleCNVsummary <- function(cnvoutfiledir, filenames, refid, mchr, stname, th = 0.1, trans = 0.9999) {
    # TODO add HiddenMarkov as a suggested library to NAMESPACE
    library(HiddenMarkov)
    
    # the transition matrix (pi) and initial state probabilities (delta)
    #   State 1: indicates a copy loss w.r.t. reference
    #   State 2: indicates copy count matches reference
    #   State 3: indicates a copy gain w.r.t. reference
    pi <- matrix(
        c(
            trans,
            (1 - trans)/2,
            (1 - trans)/2,
            
            (1 - trans)/2,
            trans,
            (1 - trans)/2,
            
            (1 - trans)/2,
            (1 - trans)/2,
            trans),
        3,
        3)
    delta <- c(0, 1, 0)
    nfile <- length(filenames)
    for (chri in mchr) {
        
        # this load brings pos and id (from simpleCNVdata) into scope
        xname2 <- paste(cnvoutfiledir, "/pos", chri, sep = "", collapse = "")
        load(xname2)
        
        # this load brings inte (from simpleCNVdata) for the reference into scope
        xname2 <- paste(
            cnvoutfiledir,
            "/",
            gsub(".CEL", "Chr", filenames[refid]),
            chri,
            sep = "",
            collapse = "")
        load(xname2)
        
        # rename inte so it doesn't conflict with subsequent loads
        ref <- inte
        
        # cnv stores maximum likelihood CNV states (1=loss, 2=no change, 3=gain)
        # per-sample, per-probeset (probesets are rows and samples are columns)
        cnv <- NULL
        
        # cnvtable is used to hold more detailed info about the gain and loss
        # intervals
        cnvtable <- NULL
        for (i in 1:nfile) {
            # this load brings inte (from simpleCNVdata) for the current CEL file into scope
            xname2 <- paste(
                cnvoutfiledir,
                "/",
                gsub(".CEL", "Chr", filenames[i]),
                chri,
                sep = "",
                collapse = "")
            load(xname2)
            
            # find the maximum likelihood states using viterbi for the current CEL file
            # and append column to the cnv matrix
            a <- inte/ref
            m <- mean(a)
            b <- dthmm(a, pi, delta, "norm", list(mean = c(m - th, m, m + th), sd = c(0.05, 0.05, 0.05)))
            states <- Viterbi(b)
            cnv <- cbind(cnv, states)
            
            l <- cbind(which(states == 1), pos[states == 1])
            g <- cbind(which(states == 3), pos[states == 3])
            ll <- length(l)
            lg <- length(g)
            if (lg > 2) {
                
                # calculate start and end position of gains
                k <- diff(g[, 1])
                ep <- which(k > 1)
                sp <- c(1, ep + 1)
                ep <- c(ep, lg/2)
                s <- length(sp)
                
                if (s == 1) 
                  cnvtable <- rbind(
                      cnvtable,
                      c(
                          stname[i],             # name
                          "gain",                # status
                          g[sp, 2],              # start pos
                          g[ep, 2],              # end pos
                          g[ep, ] - g[sp, ] + 1, # num probesets, size
                          sum(id[g[sp, 1]:g[ep, 1]] == 1),  # num SNPs
                          sum(id[g[sp, 1]:g[ep, 1]] == 2),  # num exon1
                          sum(id[g[sp, 1]:g[ep, 1]] == 3),  # num exon2
                          round(mean(ref[g[sp, 1]:g[ep, 1]]), 3),   # mean ref intensity
                          round(mean(inte[g[sp, 1]:g[ep, 1]]), 3))) # mean sample intensity loss
                if (s > 1) {
                  k <- cbind(
                      g[sp, 2],                # col1: start position
                      g[ep, 2],                # col2: end position
                      g[ep, ] - g[sp, ] + 1,   # col3,4: end index - start index + 1, end pos - start pos + 1
                      g[sp, 1],                # col5: start index
                      g[ep, 1])                # col6: end index
                  k1 <- apply(
                      k,
                      1,
                      function(x, ref, id, inte) {
                          c(
                              sum(id[x[6]:x[5]] == 1), # how many in gain are SNPs
                              sum(id[x[6]:x[5]] == 2), # how many in gain are exon1
                              sum(id[x[6]:x[5]] == 3), # how many in gain are exon2
                              mean(ref[x[6]:x[5]]),    # mean reference intensity of gain
                              mean(inte[x[6]:x[5]]))   # mean sample intensity of gain
                      },
                      ref,
                      id,
                      inte)
                  cnvtable <- rbind(
                      cnvtable,
                      cbind(
                          rep(stname[i], s), # Name
                          rep("gain", s),    # Status
                          k[, 1:4],          # StartPosition, EndPosition, Num Probe sets, Size
                          t(round(k1, 3))))  # #SNPs, #exon1, #exon2, mean ref inten, mean sample inten
                }
            }
            if (ll > 2) {
                
                # calculate the start and end position of losses
                k <- diff(l[, 1])
                ep <- which(k > 1)
                sp <- c(1, ep + 1)
                ep <- c(ep, ll/2)
                s <- length(sp)
                
                if (s == 1) 
                  cnvtable <- rbind(
                      cnvtable,
                      c(
                          stname[i],             # name
                          "loss",                # status
                          l[sp, 2],              # start pos
                          l[ep, 2],              # end pos
                          l[ep, ] - l[sp, ] + 1, # num probesets, size
                          sum(id[l[sp, 1]:l[ep, 1]] == 1),  # num SNPs
                          sum(id[l[sp, 1]:l[ep, 1]] == 2),  # num exon1
                          sum(id[l[sp, 1]:l[ep, 1]] == 3),  # num exon2
                          round(mean(ref[l[sp, 1]:l[ep, 1]]), 3),   # mean ref intensity
                          round(mean(inte[l[sp, 1]:l[ep, 1]]), 3))) # mean sample intensity
                if (s > 1) {
                  k <- cbind(
                      l[sp, 2],                # col1: start position
                      l[ep, 2],                # col2: end position
                      l[ep, ] - l[sp, ] + 1,   # col3,4: end index - start index + 1, end pos - start pos + 1
                      l[sp, 1],                # col5: start index
                      l[ep, 1])                # col6: end index
                  k1 <- apply(
                      k,
                      1,
                      function(x, ref, id, inte) {
                          c(
                              sum(id[x[6]:x[5]] == 1), # how many in loss are SNPs
                              sum(id[x[6]:x[5]] == 2), # how many in loss are exon1
                              sum(id[x[6]:x[5]] == 3), # how many in loss are exon2
                              mean(ref[x[6]:x[5]]),    # mean reference intensity of loss
                              mean(inte[x[6]:x[5]]))   # mean sample intensity of loss
                      },
                      ref,
                      id,
                      inte)
                  cnvtable <- rbind(
                      cnvtable,
                      cbind(
                          rep(stname[i], s), # Name
                          rep("loss", s),    # Status
                          k[, 1:4],          # StartPosition, EndPosition, Num Probe sets, Size
                          t(round(k1, 3))))  # #SNPs, #exon1, #exon2, mean ref inten, mean sample inten
                }
            }
        }
        
        # nothing to do here if we haven't observed any CNVs
        if(length(cnvtable) > 0)
        {
            if(length(cnvtable) == 11)
            {
                cnvtable <- matrix(c(cnvtable[1:2], chri, cnvtable[3:11]), nrow = 1)
            }
            else
            {
                cnvtable <- cbind(cnvtable[, 1:2], rep(chri, nrow(cnvtable)), cnvtable[, 3:11])
            }
            
            colnames(cnvtable) <- c(
                "Name",
                "Status",
                "Chr",
                "StartPosition",
                "EndPosition",
                "Number of Probe sets",
                "Size",
                "Number of SNP probe sets",
                "Number of exon1 sets",
                "Number of exon2 sets",
                "Mean intensity of reference sample",
                "Mean intensity of sample")
            colnames(cnv) <- stname
            rownames(cnvtable) <- NULL
            rownames(cnv) <- NULL
            xname <- paste(cnvoutfiledir, "/cnvChr", chri, sep = "", collapse = "")
            save(cnv, cnvtable, file = xname)
            xname <- paste(cnvoutfiledir, "/cnvChr", chri, ".txt", sep = "", collapse = "")
            write.table(cnvtable, file = xname, quote = FALSE, row.names = FALSE, sep = "\t")
        }
    }
} 
