# TODO make sure it's well documented that the SNPs in "genotypes must match the
#      order of the SNPs in snpProbeInfo (otherwise we should do something
#      allowing for SNP IDs to be used... actually this is what we should do)
buildPennCNVInputFiles <- function(
    outdir = getwd(), allowOverwrite = FALSE,
    genotypes, snpProbeInfo, snpInfo, snpReferenceDistribution = NULL,
    invariantProbeInfo, invariantProbesetInfo, invariantReferenceDistribution = NULL,
    transformMethod = c("CCStrans", "MAtrans"), celFiles = expandCelFiles(getwd()),
    chromosomes = c(1:19, "X", "Y", "M"),
    chromosomeRenameMap = list(X = 20, Y = 21, M = 22),
    cacheDir = tempdir(),
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
        chrChunks[[currChr]] <- .chunkIndices(chrProbesetCount, probesetChunkSize)
    }
    
    for(celfile in celFiles)
    {
        # wrapping this up as a local function def allows us to avoid doing the
        # normalization work unless it's necessary
        makeMsList <- function()
        {
            .normalizeCelFileByChr(
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
                chunkFile <- .chunkFileName(cacheDir, "snp", celfile, currChr, probesetChunkSize, chunkIndex)
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
    
    # initialize the connections for BAF & LRR files and write the header rows
    lrrAndBafBaseNames <- paste(.fileBaseWithoutExtension(celFiles), ".txt", sep = "")
    lrrAndBafOutputFiles <- file.path(outdir, lrrAndBafBaseNames)
    if(!allowOverwrite && any(file.exists(lrrAndBafOutputFiles)))
    {
        stop("refusing to overwrite the following files: ",
             paste(lrrAndBafOutputFiles[file.exists(lrrAndBafOutputFiles)], collapse = ", "),
             ". In order to procede you can either set allowOverwrite to TRUE or ",
             "you can manually delete the files.")
    }
    lrrAndBafOutputFiles <- as.list(lrrAndBafOutputFiles)
    names(lrrAndBafOutputFiles) <- .fileBaseWithoutExtension(celFiles)
    
    pfbOutputFile <- file.path(outdir, "pfbdata.txt")
    if(!allowOverwrite && file.exists(pfbOutputFile))
    {
        stop("refusing to overwrite : ", pfbOutputFile,
            ". In order to procede you can either set allowOverwrite to TRUE or ",
            "you can manually delete the files.")
    }
    
    # write the header row for LRR and BAF files
    for(celName in names(lrrAndBafOutputFiles))
    {
        con <- file(lrrAndBafOutputFiles[[celName]], "wt")
        header <- c(
            "Name",
            paste(celName, "B Allele Freq", sep = "."),
            paste(celName, "Log R Ratio", sep = "."))
        write.table(
            matrix(header, nrow = 1),
            file = con,
            quote = FALSE,
            sep = "\t",
            row.names = FALSE,
            col.names = FALSE)
        close(con)
    }
    
    pfbConnection <- file(pfbOutputFile, "wt")
    write.table(
        matrix(c("Name", "Chr", "Position", "PFB"), nrow = 1),
        file = pfbConnection,
        quote = FALSE,
        sep = "\t",
        row.names = FALSE,
        col.names = FALSE)
    
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
        
        if(currChr %in% names(chromosomeRenameMap))
        {
            chrRename <- chromosomeRenameMap[[currChr]]
            chrSnpInfo$chrId <- rep(as.character(chrRename), length(chrSnpInfo$chrId))
        }
        
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
                chunkFile <- .chunkFileName(cacheDir, "snp", celfile, currChr, probesetChunkSize, chunkIndex)
                load(chunkFile)
                
                mMatrix[, fileIndex] <- mChunk
                sMatrix[, fileIndex] <- sChunk
            }
            
            .appendToPennCNVForSNPs(
                snpInfo = chrSnpInfo[chunkIndices, ],
                intensityConts = mMatrix,
                intensityAvgs = sMatrix,
                genotypes = chrGenos[chunkIndices, ],
                lrrAndBafOutputFiles = lrrAndBafOutputFiles,
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
        invariantChrChunks[[currChr]] <- .chunkIndices(chrProbesetCount, probesetChunkSize)
    }

    # append the invariant LRR/BAF values
    for(celfile in celFiles)
    {
        makeNormInvariantList <- function()
        {
            #TODO remember to add invariantProbesetInfo !!!!!!!!!!!!!!!!!!!!!!!!!!!
            .normalizeCelFileForInvariants(
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
                chunkFile <- .chunkFileName(cacheDir, "invariant", celfile, currChr, probesetChunkSize, chunkIndex)
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
        
        if(currChr %in% names(chromosomeRenameMap))
        {
            chrRename <- chromosomeRenameMap[[currChr]]
            chrInvariantProbesetInfo$chrId <- rep(as.character(chrRename), length(chrInvariantProbesetInfo$chrId))
        }
        
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
                chunkFile <- .chunkFileName(cacheDir, "invariant", celfile, currChr, probesetChunkSize, chunkIndex)
                load(chunkFile)
                
                intensityMatrix[, fileIndex] <- intensityChunk
            }
            
            .appendToPennCNVForInvariants(
                probesetInfo = chrInvariantProbesetInfo[chunkIndices, ],
                probesetIntensities = intensityMatrix,
                lrrAndBafOutputFiles = lrrAndBafOutputFiles,
                pfbConnection = pfbConnection)
        }
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
#   lrrAndBafOutputFiles:
#       vector of file names to append to for LRR and BAF data to (one
#       file name per sample)
#   pfbConnection:
#       the connection to use for the PBF file.
.appendToPennCNVForInvariants <- function(
    probesetInfo,
    probesetIntensities,
    lrrAndBafOutputFiles,
    pfbConnection)
{
    probesetCount <- nrow(probesetInfo)
    sampleCount <- length(lrrAndBafOutputFiles)
    
    if(!all(c("probesetId", "chrId", "positionBp") %in% names(probesetInfo)))
    {
        stop("the \"probesetInfo\" parameter must contain the following components: ",
             "probesetId, chrId, positionBp")
    }
    
    if(!inherits(pfbConnection, "connection"))
    {
        stop("pfbConnection should inherit from the \"connection\" class")
    }
    
    if(!all(sapply(lrrAndBafOutputFiles, inherits, "character")))
    {
        stop("all lrrAndBafOutputFiles should inherit from the \"character\" class")
    }
    
    # set all BAFs to 2 indicating that it's not a polymorphic probeset
    bafs <- rep(2, probesetCount)
    lrrs <- log2(probesetIntensities / apply(probesetIntensities, 1, mean))
    
    for(sampleIndex in 1 : sampleCount)
    {
        con <- file(lrrAndBafOutputFiles[[sampleIndex]], "at")
        write.table(
            data.frame(probesetInfo$probesetId, lrrs[, sampleIndex], bafs),
            file = con,
            quote = FALSE,
            sep = "\t",
            row.names = FALSE,
            col.names = FALSE)
        close(con)
    }
    
    # write the PFB (Population frequency of B allele) file using mean BAF in
    # the PBF column
    write.table(
        data.frame(probesetInfo$probesetId, probesetInfo$chrId, probesetInfo$positionBp, bafs),
        file = pfbConnection,
        quote = FALSE,
        sep = "\t",
        row.names = FALSE,
        col.names = FALSE)
}

# append PennCNV data for the given LRR/BAF files and the PFB file using SNP
# data (thes files should already have a header in place since this function
# will not create one)
#
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
#   lrrAndBafOutputFiles:
#       vector of file names to append to for LRR and BAF data to (one
#       file name per sample)
#   pfbConnection:
#       the connection to use for the PBF file.
.appendToPennCNVForSNPs <- function(
    snpInfo,
    intensityConts,
    intensityAvgs,
    genotypes,
    lrrAndBafOutputFiles,
    pfbConnection)
{
    if(!inherits(pfbConnection, "connection"))
    {
        stop("pfbConnection should inherit from the \"connection\" class")
    }
    
    if(!all(sapply(lrrAndBafOutputFiles, inherits, "character")))
    {
        stop("all lrrAndBafOutputFiles should inherit from the \"character\" class")
    }
    
    snpCount <- nrow(snpInfo)
    sampleCount <- length(lrrAndBafOutputFiles)
    
    # preallocate BAF and LRR matrices for speed then calculate BAF and LRR for
    # each SNP
    bafs <- matrix(0.0, nrow = snpCount, ncol = sampleCount)
    lrrs <- matrix(0.0, nrow = snpCount, ncol = sampleCount)
    for(snpIndex in 1 : snpCount)
    {
        bafAndLrr <- .calcLRRAndBAF(
            intensityConts[snpIndex, ],
            intensityAvgs[snpIndex, ],
            genotypes[snpIndex, ])
        bafs[snpIndex, ] <- bafAndLrr$BAF
        lrrs[snpIndex, ] <- bafAndLrr$LRR
    }
    
    # write the LRR/BAF files per-sample
    for(sampleIndex in 1 : sampleCount)
    {
        con <- file(lrrAndBafOutputFiles[[sampleIndex]], "at")
        write.table(
            data.frame(snpInfo$snpId, lrrs[, sampleIndex], bafs[, sampleIndex]),
            file = con,
            quote = FALSE,
            sep = "\t",
            row.names = FALSE,
            col.names = FALSE)
        close(con)
    }
    
    # write the PFB (Population frequency of B allele) file using mean BAF in
    # the PBF column
    write.table(
        data.frame(snpInfo$snpId, snpInfo$chrId, snpInfo$positionBp, apply(bafs, 1, mean)),
        file = pfbConnection,
        quote = FALSE,
        sep = "\t",
        row.names = FALSE,
        col.names = FALSE)
}

# calculates LRR and BAF values for a single SNP position
#
# PARAMETERS:
#   intensityConts:
#       per-sample intensity contrasts (A allele vs B allele)
#   intensityAvgs:
#       per-sample intensity averages (A averaged with B)
#   genos:
#       per-sample genotypes 1 = AA, 2 = AB, 3 = BB
#
# RETURNS:
#   A data.frame object with B-Allele Frequency (BAF) and Log R Ratio (LRR)
#   components. The row count of this dataframe will equal length(genos)
.calcLRRAndBAF <- function(intensityConts, intensityAvgs, genos)
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
            contVariancePerGeno[genoIndex] <- .bivar(currIntensityConts)
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
            LRR <- log2(intensityAvgs/medianAvgsPerGeno)
        }
        
        # the genotype for all samples is AB
        else if (uniqueGenos == 2) {
            BAF <- rep(0.5, sampleCount)
            LRR <- log2(intensityAvgs/medianAvgsPerGeno)
        }
        
        # the genotype for all samples is BB
        else if (uniqueGenos == 3) {
            BAF <- rep(1, sampleCount)
            LRR <- log2(intensityAvgs/medianAvgsPerGeno)
        }
    }
    
    data.frame(BAF = BAF, LRR = LRR)
}

.normalizeCelFileForInvariants <- function(
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

# a little function to make sure that we turn any item which is not a list into
# a list (which makes some algorithms more general/consistent)
.listify <- function(x)
{
    if(is.na(x) || is.null(x))
    {
        x <- list()
    }
    else if(inherits(x, "data.frame") || !is.list(x))
    {
        x <- list(x)
    }
    
    x
}

simpleCNV <- function(
    snpProbeInfo, snpInfo, snpReferenceDistribution = NULL,
    invariantProbeInfo, invariantProbesetInfo, invariantReferenceDistribution = NULL,
    celFiles = expandCelFiles(getwd()),
    referenceCelFile,
    chromosomes = c(1:19, "X", "Y", "M"),
    verbose = FALSE,
    cluster = NULL)
{
    snpCount <- nrow(snpInfo)
    
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
    
    # make the invariant data all lists for consistency
    invariantProbeInfo <- .listify(invariantProbeInfo)
    invariantProbesetInfo <- .listify(invariantProbesetInfo)
    invariantReferenceDistribution <- .listify(invariantReferenceDistribution)
    if(length(invariantReferenceDistribution) == 0)
    {
        for(i in 1 : length(invariantProbeInfo))
        {
            invariantReferenceDistribution[[i]] <- NULL
        }
    }
    
    invariantGroupCount <- length(invariantProbeInfo)
    if(invariantGroupCount != length(invariantProbesetInfo) ||
       invariantGroupCount != length(invariantReferenceDistribution))
    {
        stop("there is a missmatch between the \"invariantProbeInfo\", ",
            "\"invariantProbesetInfo\", \"invariantReferenceDistribution\"")
    }
    
    for(i in 1 : invariantGroupCount)
    {
        if(!inherits(invariantProbeInfo[[i]], "data.frame") ||
            !all(c("probeIndex", "probesetId") %in% names(invariantProbeInfo[[i]])))
        {
            stop("You must supply a \"invariantProbeInfo\" data frame parameter which has ",
                "at a minimum the \"probeIndex\", and \"probesetId\" ",
                "components. Please see the help documentation for more details.")
        }
        
        if(!inherits(invariantProbesetInfo[[i]], "data.frame") ||
            !all(c("probesetId", "chrId", "positionBp") %in% names(invariantProbesetInfo[[i]])))
        {
            stop("You must supply a \"invariantProbesetInfo\" data frame parameter which has ",
                "at a minimum the \"probesetId\", \"chrId\" and \"positionBp\"",
                "components. Please see the help documentation for more details.")
        }
        
        if(!is.null(invariantReferenceDistribution[[i]]) &&
            !is.numeric(invariantReferenceDistribution[[i]]))
        {
            stop("The \"invariantReferenceDistribution\" should either be ",
                "numeric or NULL")
        }
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
    sampleCount <- length(celFiles)
    
    # make sure that the chromosome vector is not numeric
    # TODO is toupper the right thing to do here? (I do it in MDgenotype.R too)
    # TODO change code to deal correctly with a case where there is SNP
    #      data but no exon data and vice-versa
    chromosomes <- toupper(as.character(chromosomes))
    
    snpChromosomes <- unique(snpInfo$chrId)
    if(!all(chromosomes %in% snpChromosomes))
    {
        warning(
            "SNP data for the following requested chromosomes are not available: ",
            paste(setdiff(chromosomes, snpChromosomes), collapse = ", "),
            ". These chromosomes will be skipped.")
    }
    chromosomes <- intersect(chromosomes, snpChromosomes)
    rm(snpChromosomes)
    
    invariantChromosomes <- as.character(unique(unlist(
            lapply(invariantProbesetInfo, function(x) {x$chrId}),
            use.names=F)))
    if(!all(chromosomes %in% invariantChromosomes))
    {
        warning(
            "Invariant data for the following requested chromosomes are not available: ",
            paste(setdiff(chromosomes, invariantChromosomes), collapse = ", "),
            ". These chromosomes will be skipped.")
    }
    chromosomes <- unique(chromosomes, invariantChromosomes)
    rm(invariantChromosomes)
    
    if(length(chromosomes) == 0)
    {
        stop("Stopping CNV analysis. There are no chromosomes to process")
    }
    
    if(verbose)
    {
        cat("processing CEL files\n")
    }
    
    # the reference intensities
    refIntensities <- normalizeForSimpleCNV(
        referenceCelFile,
        chromosomes,
        snpProbeInfo, snpInfo, snpReferenceDistribution,
        invariantProbeInfo, invariantProbesetInfo, invariantReferenceDistribution,
        verbose)
    normRefPath <- normalizePath(referenceCelFile)
    
    # a list of matrices by chromosome initialized to NULLs
    allCnvsByChr <- list()
    for(currCelFile in celFiles)
    {
        if(normRefPath == normalizePath(currCelFile))
        {
            # the reference cannot have any CNV losses or gains against itself
            noCopyChange <- 2
            cnvs <- list()
            for(currChr in names(refIntensities))
            {
                cnvs[[currChr]] <- rep(noCopyChange, length(refIntensities[[currChr]]))
                names(cnvs[[currChr]]) <- names(refIntensities[[currChr]])
            }
        }
        else
        {
            # normalizes intensities and sorts by position
            currIntensities <- normalizeForSimpleCNV(
                currCelFile,
                chromosomes,
                snpProbeInfo, snpInfo, snpReferenceDistribution,
                invariantProbeInfo, invariantProbesetInfo, invariantReferenceDistribution,
                verbose)
            
            if(length(cluster) == 0)
            {
                cnvs <- list()
                for(chr in names(currIntensities))
                {
                    if(verbose)
                    {
                        cat("inferring CNVs for chromosome ", chr, " of ", currCelFile, "\n")
                    }
                    cnvs[[chr]] <- .inferCNVFromIntensity(
                        currIntensities[[chr]],
                        refIntensities[[chr]])
                }
            }
            else
            {
                cat("Applying cluster resources to infer CNVs for ", currCelFile, "\n")
                
                # we need to do some strange data massaging here so that we
                # can apply cluster resources
                combinedIntList <- mapply(
                    function(intToTest, refInt)
                    {
                        list(testInt = intToTest, refInt = refInt)
                    },
                    currIntensities,
                    refIntensities,
                    SIMPLIFY = FALSE)
                cnvs <- parLapply(cluster, combinedIntList, .applyInferCNVFromIntensity)
            }
        }
        
        for(chr in names(cnvs))
        {
            allCnvsByChr[[chr]] <- cbind(allCnvsByChr[[chr]], cnvs[[chr]])
        }
    }
    
    # give column names based on the CEL file
    for(chr in names(allCnvsByChr))
    {
        colnames(allCnvsByChr[[chr]]) <- .fileBaseWithoutExtension(celFiles)
    }
    
    allCnvsByChr
}

.applyInferCNVFromIntensity <- function(combinedInt)
{
    .inferCNVFromIntensity(combinedInt$testInt, combinedInt$refInt)
}

.inferCNVFromIntensity <- function(
    intensities,
    refIntensities,
    sameStateProb = 0.9999,
    th = NULL,
    stdDev = NULL)
{
    library(HiddenMarkov)
    
    # the transition matrix (pi) and initial state probabilities (delta)
    #   State 1: indicates a copy loss w.r.t. reference
    #   State 2: indicates copy count matches reference
    #   State 3: indicates a copy gain w.r.t. reference
    pi <- matrix(
        c(
            sameStateProb,
            (1 - sameStateProb)/2,
            (1 - sameStateProb)/2,
            
            (1 - sameStateProb)/2,
            sameStateProb,
            (1 - sameStateProb)/2,
            
            (1 - sameStateProb)/2,
            (1 - sameStateProb)/2,
            sameStateProb),
        3,
        3)
    delta <- c(0, 1, 0)
    
    a <- intensities / refIntensities
    m <- mean(a)
    
    if(is.null(stdDev))
    {
        stdDev <- sqrt(var(a))
    }
    
    if(is.null(th))
    {
        th <- 2.58 * stdDev
    }
    
    b <- dthmm(a, pi, delta, "norm", list(mean = c(m - th, m, m + th), sd = rep(stdDev, 3)))
    
    estimatedCNVStates <- Viterbi(b)
    names(estimatedCNVStates) <- names(intensities)
    estimatedCNVStates
}

normalizeForSimpleCNV <- function(
    celFileName,
    chromosomes,
    snpProbeInfo, snpInfo, snpReferenceDistribution,
    invariantProbeInfo, invariantProbesetInfo, invariantReferenceDistribution,
    verbose)
{
    if(verbose) cat("Reading and normalizing CEL file: ", celFileName, "\n", sep="")
    
    invariantGroupCount <- length(invariantProbeInfo)
    
    # extract log2(mean intensity) for each probe
    celData <- read.celfile(celFileName, intensity.means.only = TRUE)
    celData <- log2(as.matrix(celData[["INTENSITY"]][["MEAN"]]))
    
    # normalize invariants
    invY <- list()
    for(i in 1 : invariantGroupCount)
    {
        invY[[i]] <- celData[invariantProbeInfo[[i]]$probeIndex, , drop = FALSE]
        if(length(invariantProbeInfo[[i]]$correction) > 0)
            # C+G and fragment length correction for Y
            invY[[i]] <- invY[[i]] + invariantProbeInfo[[i]]$correction
        if(length(invariantReferenceDistribution[[i]]) > 0)
            invY[[i]] <- normalize.quantiles.use.target(invY[[i]], target = invariantReferenceDistribution[[i]])
        
        invY[[i]] <- subColSummarizeMedian(matrix(invY[[i]], ncol = 1), invariantProbeInfo[[i]]$probesetId)
        
        # "vectorize" invY[[i]]
        invY[[i]] <- invY[[i]][ , 1, drop = TRUE]
    }
    
    # normalize SNPs
    snpY <- celData[snpProbeInfo$probeIndex, , drop = FALSE]
    if(length(snpProbeInfo$correction) > 0)
        # C+G and fragment length correction for Y
        snpY <- snpY + snpProbeInfo$correction
    if(length(snpReferenceDistribution) > 0)
        snpY <- normalize.quantiles.use.target(snpY, target = snpReferenceDistribution)
    
    # TODO pretty much all of the snpId naming used in this function is brittle. Make it robust
    
    # separate A alleles from B alleles and summarize to the probeset level
    aProbeInt <- snpY[snpProbeInfo$isAAllele, 1, drop = FALSE]
    bProbeInt <- snpY[!snpProbeInfo$isAAllele, 1, drop = FALSE]
    aSnpInt <- subColSummarizeMedian(
        matrix(aProbeInt, ncol = 1),
        snpProbeInfo$snpId[snpProbeInfo$isAAllele])
    bSnpInt <- subColSummarizeMedian(
        matrix(bProbeInt, ncol = 1),
        snpProbeInfo$snpId[!snpProbeInfo$isAAllele])
    
    aSnpIds <- rownames(aSnpInt)
    bSnpIds <- rownames(bSnpInt)
    if(!all(aSnpIds == bSnpIds))
    {
        stop("The SNP IDs for the A alleles should match up with the SNP IDs ",
            "for the B alleles but they do not.")
    }
    
    # pull probe intensities from A or B SNPs depending on which one has a
    # higher intensity, then summarize
    aGreaterIds <- aSnpIds[aSnpInt > bSnpInt]
    aGreaterProbeIndices <- which(snpProbeInfo$snpId[snpProbeInfo$isAAllele] %in% aGreaterIds)
    highProbeInt <- bProbeInt
    highProbeInt[aGreaterProbeIndices, ] <- aProbeInt[aGreaterProbeIndices, ]
    #TODO using invariantReferenceDistribution[[1]] here seems questionable
    #     Hyuna's data uses the exact same ref distribution for the two classes
    #     of invariants that we deal with. Investigate whether it's possible
    #     for us to rely on this always being true (in which case
    #     invariantReferenceDistribution could just be a single vector rather
    #     than a list of vectors)
    highProbeInt <- normalize.quantiles.use.target(highProbeInt, target = invariantReferenceDistribution[[1]])
    highSnpInt <- subColSummarizeMedian(
        matrix(highProbeInt, ncol = 1),
        snpProbeInfo$snpId[snpProbeInfo$isAAllele])
    
    # "vectorize" highSnpInt
    highSnpInt <- highSnpInt[ , 1, drop = TRUE]
    
    if(length(snpInfo$snpId) != length(highSnpInt))
    {
        stop("internal error: vector lengths should match but they do not")
    }
    
    snpIdOrdering <- match(snpInfo$snpId, names(highSnpInt))
    
    if(any(is.na(snpIdOrdering)))
        stop("Failed to match up SNP probe IDs")
    highSnpInt <- highSnpInt[snpIdOrdering]
    
    for(i in 1 : invariantGroupCount)
    {
        if(length(invariantProbesetInfo[[i]]$probesetId) != length(invY[[i]]))
        {
            stop("internal error: vector lengths should match but they do not")
        }
        
        invIdOrdering <- match(invariantProbesetInfo[[i]]$probesetId, names(invY[[i]]))
        if(any(is.na(invIdOrdering)))
            stop("Failed to match up invariant probe IDs")
        invY[[i]] <- invY[[i]][invIdOrdering]
    }
    
    combinedIntensities <- list()
    for(chr in chromosomes)
    {
        # combine the SNP and invariant intensities, then sort by position
        chrInvInt <- NULL
        chrInvPos <- NULL
        for(i in 1 : invariantGroupCount)
        {
            chrInvIndices <- which(invariantProbesetInfo[[i]]$chrId == chr)
            chrInvInt <- c(chrInvInt, invY[[i]][chrInvIndices])
            chrInvPos <- c(chrInvPos, invariantProbesetInfo[[i]]$positionBp[chrInvIndices])
        }
        
        chrSnpIndices <- which(snpInfo$chrId == chr)
        chrSnpInt <- highSnpInt[chrSnpIndices]
        chrSnpPos <- snpInfo$positionBp[chrSnpIndices]
        
        combinedInt <- c(chrInvInt, chrSnpInt)
        combinedPos <- c(chrInvPos, chrSnpPos)
        combinedInt <- combinedInt[order(combinedPos)]
        
        combinedIntensities[[chr]] <- combinedInt
    }
    
    combinedIntensities
}
