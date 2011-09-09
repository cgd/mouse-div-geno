#============================================== read the cel files
# Definition : Vinotype the array. Used inside of genotypthis function
#
# input :
# nm : Contrast ( A allele intensity - B allele intensity ) of one probe set.
# ns : Average of one probe set.
# geno : genotype obtained from genotype() function. 1 = AA, 2 = AB, 3 = BB
#
# output :
# vino : 1 indicates vino (all others should be zero)
# <example>
# load('MM'); load('SS');
# nm = MM[1,]; ns = SS[1,]
# geno = genotype(nm, ns, hint = NULL)
# vino = vinotype(nm, ns, geno)
#======================================================================

vinotype <- function(nm, ns, geno, logfn=NULL) {
    doLog <- !is.null(logfn)
    
    # logging code expects nm vector elements to be named
    if(doLog && is.null(names(nm))) {
        names(nm) <- paste("Sample", seq_along(nm), sep="")
    }
    sampleNames <- names(nm)
    
    if(doLog) {
        logfn("performing SNP VInOtyping")
    }
    
    #=========== find centers
    iig <- sort(unique(geno))
    ngeno <- length(iig)
    if(doLog) {
        logfn("vinotyping using observed genotype(s): %s", paste(iig, collapse=", "))
    }
    
    # identifies likely H's
    rmid <- nm >= -0.3 & nm <= 0.3 & ns < quantile(ns, 0.2)
    if(doLog && any(rmid)) {
        logfn(
            "the following %i sample(s) are excluded for the purposes of calculating a VInO threshold:",
            sum(rmid))
        for(n in sampleNames[rmid]) {
            logfn("    %s", n)
        }
    }
    
    # key ideas:
    #   * get some kind of threshold which is mms to decide if there
    #     is no chance of being vino
    #   * if there are many vinos it will bring down the whole threshold
    #   * often vino is near H rather than A or B
    if (ngeno == 1) {
        mms <- median(ns[!rmid])
    }
    else if (ngeno == 3) {
        # for 3 we can use median pair
        mms <- min(tapply(ns[geno == 1 | geno == 3], geno[geno == 1 | geno == 3], median))
    }
    else if (ngeno == 2) {
        mms <- min(tapply(ns[!rmid], geno[!rmid], median))
    }
    mm <- rep(0, 3)
    ms <- rep(0, 3)
    adata <- cbind(nm, ns)
    nsize <- length(geno)
    vino <- rep(-1, nsize)
    
    # 2 indicates it's definitely not a vino
    aboveThresh <- ns > mms
    vino[aboveThresh] <- 2
    vino1 <- vino
    ss <- list()
    m <- matrix(0, 2, 2)
    lm <- 0
    
    if(doLog) {
        if(any(aboveThresh)) {
            logfn(
                "the following %i samples are above the VInO threshold of %f and therefore are not VInOs:",
                sum(aboveThresh),
                mms)
            for(n in sampleNames[aboveThresh]) {
                logfn("    %s", n)
            }
        } else {
            logfn("all samples are below the VInO threshold of %f and will be tested as vinos", mms)
        }
    }
    
    # iterate through the genotypes
    for (ik in iig) {
        # count how many genos match the current genotype
        l <- sum(geno == ik)
        if (l > 1) {
            nns <- ns[geno == ik]
            nnm <- nm[geno == ik]
            onn <- order(nns)
            onns <- nns[onn]
            onnm <- nnm[onn]
            
            # th => which indices might be vinos
            th <- which(onns < mms)
            if (length(th) > 0) {
                th <- max(th)
                if (l > th) 
                  th <- th + 1
                
                # calculates the difference between nearest ns's
                donns <- diff(onns)
                
                # when gap is big that's captured by boxplot's b$out
                # using a big range (5) allows us to detect outliers
                tb <- boxplot(donns, range = 5, plot = FALSE)
                if (sum(tb$out > 0.2) > 0) {
                  # gap is big and > 0.2 (extreme cases)
                  tb$out <- tb$out[tb$out > 0.2]
                  k <- match(tb$out, donns)
                  k <- k[k < th]
                  if (length(k) > 0) {
                    # find them and mark them as vino here
                    vino[geno == 2 & ns <= onns[max(k)]] <- 1
                    if(doLog && any(geno == 2 & ns <= onns[max(k)])) {
                        logfn(
                            "the following hets will be set to vino since we detected an intensity gap and AB averaged intensity <= %f:",
                            onns[max(k)])
                        for(n in sampleNames[geno == 2 & ns <= onns[max(k)]]) {
                            logfn("    %s", n)
                        }
                    }
                  }
                }
            }
            mm[ik] <- median(nnm)
            ms[ik] <- median(nns)
            
            # ss[[ik]] becomes the covariance matrix for the current genotype's
            # nm and ns values
            ssm <- matrix(0, 2, 2)
            ssm[1, 1] <- .bivar(nnm)
            ssm[2, 2] <- .bivar(nns)
            ssm[1, 2] <- .bicov(nnm, nns)
            ssm[2, 1] <- ssm[1, 2]
            ss[[ik]] <- ssm
            
            # sum the covariances as m (we will later make this an average
            m <- m + ss[[ik]]
            lm <- lm + 1
        }
        else {
            mm[ik] <- nm[geno == ik]
            ms[ik] <- ns[geno == ik]
            
            # NA is a placeholder that we will later fill with the average
            # covariance
            ss[[ik]] <- NA
        }
    }
    
    # here we set genotypes where there is only 1 value (the NA's) to the
    # average covariance. For the genotypes that had more than 1 value we
    # give equal weight to the average covariance and the genotype specific
    # covariance
    m <- m/lm
    for (ik in iig) {
        if (is.na(ss[ik])) 
            ss[[ik]] <- m
        else ss[[ik]] <- ss[[ik]] * 0.5 + m * 0.5
    }
    
    #=========== test 2 vs. 3
    if (ngeno == 3) {
        # three groups; # lower BIC better fit.
        l1 <- sum(geno == 1)
        l2 <- sum(geno == 2)
        l3 <- sum(geno == 3)
        mdd <- rep(0, nsize)
        th <- c(0.99, 0.95, 0.99)
        if (l1 > 1 & det(ss[[1]]) > 10^(-10)) {
            tmp <- mahalanobis(adata[geno == 1, ], c(mm[1], ms[1]), ss[[1]])
            
            # dd1 gives (1 - prob) here
            dd1 <- pchisq(tmp, df = 2)
            mdd[geno == 1] <- dd1
            
            # vino1 only for aa group, (=2 means not a vino for sure)
            vino1[geno == 1][dd1 < th[1]] <- 2
            if(doLog && any(dd1 < th[1])) {
                logfn(
                    "excluding the following as candidate VInOs since geno==1 and dd1 < %f:",
                    th[1])
                for(n in sampleNames[geno == 1][dd1 < th[1]]) {
                    logfn("    %s", n)
                }
            }
        }
        if (l2 > 1 & det(ss[[2]]) > 10^(-10)) {
            tmp <- mahalanobis(adata[geno == 2, ], c(mm[2], ms[2]), ss[[2]])
            dd2 <- pchisq(tmp, df = 2)
            mdd[geno == 2] <- dd2
            vino1[geno == 2][dd2 < th[2]] <- 2
            
            if(doLog && any(dd2 < th[2])) {
                logfn(
                    "excluding the following as candidate VInOs since geno==2 and dd2 < %f:",
                    th[2])
                for(n in sampleNames[geno == 2][dd2 < th[2]]) {
                    logfn("    %s", n)
                }
            }
        }
        if (l3 > 1 & det(ss[[3]]) > 10^(-10)) {
            tmp <- mahalanobis(adata[geno == 3, ], c(mm[3], ms[3]), ss[[3]])
            dd3 <- pchisq(tmp, df = 2)
            mdd[geno == 3] <- dd3
            vino1[geno == 3][dd3 < th[3]] <- 2
            
            if(doLog && any(dd3 < th[3])) {
                logfn(
                    "excluding the following as candidate VInOs since geno==3 and dd3 < %f:",
                    th[3])
                for(n in sampleNames[geno == 3][dd3 < th[3]]) {
                    logfn("    %s", n)
                }
            }
        }
    }
    else if (ngeno == 2) {
        l1 <- sum(geno == iig[1])
        l2 <- sum(geno == iig[2])
        mdd <- rep(0, nsize)
        th <- c(0.99, 0.95, 0.99)
        th <- th[iig]
        if (l1 > 1 & det(ss[[iig[1]]]) > 10^(-10)) {
            tmp <- mahalanobis(adata[geno == iig[1], ], c(mm[iig[1]], ms[iig[1]]), 
                ss[[iig[1]]])
            dd1 <- pchisq(tmp, df = 2)
            mdd[geno == iig[1]] <- dd1
            vino1[geno == iig[1]][dd1 < th[1]] <- 2
            
            if(doLog && any(dd1 < th[1])) {
                logfn(
                    "excluding the following as candidate VInOs since geno==%i and dd1 < %f:",
                    iig[1],
                    th[1])
                for(n in sampleNames[geno == iig[1]][dd1 < th[1]]) {
                    logfn("    %s", n)
                }
            }
        }
        if (l2 > 1 & det(ss[[iig[2]]]) > 10^(-10)) {
            tmp <- mahalanobis(adata[geno == iig[2], ], c(mm[iig[2]], ms[iig[2]]), 
                ss[[iig[2]]])
            dd2 <- pchisq(tmp, df = 2)
            mdd[geno == iig[2]] <- dd2
            vino1[geno == iig[2]][dd2 < th[2]] <- 2
            
            if(doLog && any(dd2 < th[2])) {
                logfn(
                    "excluding the following as candidate VInOs since geno==%i and dd2 < %f:",
                    iig[2],
                    th[2])
                for(n in sampleNames[geno == iig[2]][dd2 < th[2]]) {
                    logfn("    %s", n)
                }
            }
        }
    }
    else if (ngeno == 1) {
        tmp <- mahalanobis(adata, c(mm[iig[1]], ms[iig[1]]), ss[[iig[1]]])
        mdd <- pchisq(tmp, df = 2)
        vino1[mdd < 0.99] <- 2
        
        if(doLog && any(mdd < 0.99)) {
            logfn(
                "excluding the following as candidate VInOs since %f < 0.99:",
                mdd)
            for(n in sampleNames[mdd < 0.99]) {
                logfn("    %s", n)
            }
        }
    }
    
    # vino is product of 2 probs. IE: assuming our 2D intensity plot, how far
    # off of center are we in any direction
    # (left, right and up all lead to high mdd)
    # For vinos we are only interested in case where intensity is really low,
    # so we assume intensity follows normal dist using the (==2) vals which we
    # know are not vinos
    alld <- 1 - pnorm(ns, mean = mean(ns[vino1 == 2]), sd = sd(ns[vino1 == 2]))
    
    # vino = outlier with in each group & low intensity
    k <- alld * mdd
    vino[(vino1 != 2 & k > 0.9999)] <- 1
    if(doLog && any(vino1 != 2 & k > 0.9999)) {
        logfn("the following are called as VInO since they haven't yet been excluded and k > 0.9999:")
        for(n in sampleNames[vino1 != 2 & k > 0.9999]) {
            logfn("    %s", n)
        }
    }
    
    ons <- order(ns)
    sn <- ns[ons]
    dsn <- diff(sn)
    
    # if gap is really big give a vino flag
    if (any(dsn > 1.5)) {
        bdsn <- dsn[dsn > 1.5]
        g <- match(bdsn, dsn)
        th <- sum(ns < median(ns))
        g <- g[g < th]
        if (length(g) > 0) {
            vino[ns <= sn[max(g)]] <- 1
            if(doLog && any(ns <= sn[max(g)])) {
                logfn("the following are set to VInO because of a large gap in intensity (despite any previous exclusions):")
                for(n in sampleNames[ns <= sn[max(g)]]) {
                    logfn("    %s", n)
                }
            }
        }
    }
    
    if (any(vino == 1) && any(vino == -1)) {
        if(doLog) {
            naVinos <- which(vino == -1)
        }
        
        # pick up remaining vinos using vdist clustering
        vino <- .vdist(nm/5, ns/(2 * (max(ns) - min(ns))), vino)
        
        if(doLog) {
            newVinos <- c()
            for(i in naVinos) {
                newVinos[length(newVinos) + 1] <- sprintf("%s=%i", sampleNames[i], vino[i])
            }
            
            logfn("the following VInO values were determined using .vdist clustering:")
            for(n in newVinos) {
                logfn("    %s", n)
            }
        }
    }
    
    vino[vino != 1] <- 0
    
    list(vino = vino, conf = 1 - mdd)
}
