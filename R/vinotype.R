#============================================== read the cel files
# Definition : Vinotype the array. Used inside of genotypthis function
#
# input :
# nm : Contrast ( A allele intensity - B allele intensity ) of one probe set.
# ns : Average of one probe set.
# geno : genotype obtained from genotype() function. 1 = AA, 2 = AB, 3 = BB
#
# output :
# vino : 1 = vino 2 and 0 non vino (0 and 2 are space holder, and do not have any meaning)
# <example>
# load('MM'); load('SS');
# nm = MM[1,]; ns = SS[1,]
# geno = genotype(nm, ns, hint = NULL)
# vino = vinotype(nm, ns, geno)
#======================================================================
bicov <- function(x, y) {
    mx <- median(x)
    my <- median(y)
    ux <- abs((x - mx)/(9 * qnorm(0.75) * mad(x)))
    uy <- abs((y - my)/(9 * qnorm(0.75) * mad(y)))
    aval <- ifelse(ux <= 1, 1, 0)
    bval <- ifelse(uy <= 1, 1, 0)
    top <- sum(aval * (x - mx) * (1 - ux^2)^2 * bval * (y - my) * (1 - uy^2)^2)
    top <- length(x) * top
    botx <- sum(aval * (1 - ux^2) * (1 - 5 * ux^2))
    boty <- sum(bval * (1 - uy^2) * (1 - 5 * uy^2))
    bi <- top/(botx * boty)
    bi
}
bivar <- function(x) {
    mx <- median(x)
    ux <- abs((x - mx)/(9 * qnorm(0.75) * mad(x)))
    aval <- ifelse(ux <= 1, 1, 0)
    top <- sum((aval * (x - mx) * (1 - ux^2)^2)^2)
    top <- length(x) * top
    botx <- (sum(aval * (1 - ux^2) * (1 - 5 * ux^2)))^2
    bi <- top/botx
    bi
}
vinotype = function(nm, ns, geno, doCNV) {
    #=========== find centers
    iig = sort(unique(geno))
    ngeno = length(iig)
    
    # key ideas:
    #   * get some kind of threshold which is mms to decide if there
    #     is no chance of being vino
    #   * if there are many vinos it will bring down the whole threshold
    #   * often vino is near H rather than A or B
    if (ngeno == 1) {
        mms = median(ns)
    }
    else if (ngeno == 3) {
        # for 3 we can use median pair
        mms = min(tapply(ns[geno == 1 | geno == 3], geno[geno == 1 | geno == 3], median))
    }
    else if (ngeno == 2) {
        # for 2 group case we want to get rid of H's when calcing the MMS threshold
        k = tapply(nm, geno, median)
        g = match(geno, names(k[k < -0.3 | k > 0.3]))
        if (sum(!is.na(g)) > 1) {
            mms = min(tapply(ns[!is.na(g)], geno[!is.na(g)], median))
        }
        else {
            mms = median(ns)
        }
    }
    mm = ms = rep(0, 3)
    adata = cbind(nm, ns)
    nsize = length(geno)
    vino = rep(-1, nsize)
    
    # 2 indicates it's definitely not a vino
    vino[(ns > mms)] = 2
    vino1 = vino
    ss = list()
    m = matrix(0, 2, 2)
    lm = 0
    
    # iterate through the genotypes
    for (ik in iig) {
        # count how many genos match the current genotype
        l = sum(geno == ik)
        if (l > 1) {
            nns = ns[geno == ik]
            nnm = nm[geno == ik]
            onn = order(nns)
            onns = nns[onn]
            onnm = nnm[onn]
            
            # th => which indices might be vinos
            th = which(onns < mms)
            if (length(th) > 0) {
                th = max(th)
                if (l > th) 
                  th = th + 1
                
                # calculates the difference between nearest ns's
                donns = diff(onns)
                
                # when gap is big that's captured by boxplot's b$out
                # using a big range (5) allows us to detect outliers
                tb = boxplot(donns, range = 5, plot = FALSE)
                if (sum(tb$out > 0.2) > 0) {
                  # gap is big and > 0.2 (extreme cases)
                  tb$out = tb$out[tb$out > 0.2]
                  k = match(tb$out, donns)
                  k = k[k < th]
                  if (length(k) > 0) {
                    # find them and mark them as vino here
                    vino[geno == 2 & ns <= onns[max(k)]] = 1
                  }
                }
            }
            mm[ik] = median(nnm)
            ms[ik] = median(nns)
            
            # ss contains variance
            # TODO more comment needed here to explain what's going on
            ssm = matrix(0, 2, 2)
            ssm[1, 1] = bivar(nnm)
            ssm[2, 2] = bivar(nns)
            ssm[1, 2] = ssm[2, 1] = bicov(nnm, nns)
            ss[[ik]] = ssm
            m = m + ss[[ik]]
            lm = lm + 1
        }
        else {
            mm[ik] = nm[geno == ik]
            ms[ik] = ns[geno == ik]
            ss[[ik]] = NA
        }
    }
    m = m/lm
    for (ik in iig) {
        if (is.na(ss[ik])) 
            ss[[ik]] = m
        else ss[[ik]] = ss[[ik]] * 0.5 + m * 0.5
    }
    BAF = llr = rep(0, nsize)
    if (doCNV) {
        cnvinputs = CNVinput(nm, ns, geno, mm, ms, ss, nsize, ngeno, iig)
        BAF = cnvinputs$BAF
        llr = cnvinputs$llr
    }
    #=========== test 2 vs. 3
    if (ngeno == 3) {
        # three groups; # lower BIC better fit.
        l1 = sum(geno == 1)
        l2 = sum(geno == 2)
        l3 = sum(geno == 3)
        mdd = rep(0, nsize)
        th = c(0.99, 0.95, 0.99)
        if (l1 > 1 & det(ss[[1]]) > 10^(-10)) {
            tmp = mahalanobis(adata[geno == 1, ], c(mm[1], ms[1]), ss[[1]])
            
            # dd1 gives (1 - prob) here
            dd1 = pchisq(tmp, df = 2)
            mdd[geno == 1] = dd1
            
            # vino1 only for aa group, (=2 means not a vino for sure)
            vino1[geno == 1][dd1 < th[1]] = 2
        }
        if (l2 > 1 & det(ss[[2]]) > 10^(-10)) {
            tmp = mahalanobis(adata[geno == 2, ], c(mm[2], ms[2]), ss[[2]])
            dd2 = pchisq(tmp, df = 2)
            mdd[geno == 2] = dd2
            vino1[geno == 2][dd2 < th[2]] = 2
        }
        if (l3 > 1 & det(ss[[3]]) > 10^(-10)) {
            tmp = mahalanobis(adata[geno == 3, ], c(mm[3], ms[3]), ss[[3]])
            dd3 = pchisq(tmp, df = 2)
            mdd[geno == 3] = dd3
            vino1[geno == 3][dd3 < th[3]] = 2
        }
    }
    if (ngeno == 2) {
        l1 = sum(geno == iig[1])
        l2 = sum(geno == iig[2])
        mdd = rep(0, nsize)
        th = c(0.99, 0.95, 0.99)
        th = th[iig]
        if (l1 > 1 & det(ss[[iig[1]]]) > 10^(-10)) {
            tmp = mahalanobis(adata[geno == iig[1], ], c(mm[iig[1]], ms[iig[1]]), 
                ss[[iig[1]]])
            dd1 = pchisq(tmp, df = 2)
            mdd[geno == iig[1]] = dd1
            vino1[geno == iig[1]][dd1 < th[1]] = 2
        }
        if (l2 > 1 & det(ss[[iig[2]]]) > 10^(-10)) {
            tmp = mahalanobis(adata[geno == iig[2], ], c(mm[iig[2]], ms[iig[2]]), 
                ss[[iig[2]]])
            dd2 = pchisq(tmp, df = 2)
            mdd[geno == iig[2]] = dd2
            vino1[geno == iig[2]][dd2 < th[2]] = 2
        }
    }
    if (ngeno == 1) {
        tmp = mahalanobis(adata, c(mm[iig[1]], ms[iig[1]]), ss[[iig[1]]])
        mdd = pchisq(tmp, df = 2)
        vino1[mdd < 0.99] = 2
    }
    
    # vino is product of 2 probs. IE: assuming our 2D intensity plot, how far
    # off of center are we in any direction
    # (left, right and up all lead to high mdd)
    # For vinos we are only interested in case where intensity is really low,
    # so we assume intensity follows normal dist using the (==2) vals which we
    # know are not vinos
    alld = 1 - pnorm(ns, mean = mean(ns[vino1 == 2]), sd = sd(ns[vino1 == 2]))
    
    # vino = outlier with in each group & low intensity
    k = alld * mdd
    vino[(vino1 != 2 & k > 0.9999)] = 1
    ons = order(ns)
    sn = ns[ons]
    dsn = diff(sn)
    
    # if gap is really big give a vino flag
    if (any(dsn > 1.5)) {
        bdsn = dsn[dsn > 1.5]
        g = match(bdsn, dsn)
        th = sum(ns < median(ns))
        g = g[g < th]
        if (length(g) > 0) 
            vino[ns <= sn[max(g)]] = 1
    }
    
    if (any(vino == 1) & any(vino == -1)) {
        # pick up remaining vinos using vdist clustering
        vino = vdist(nm/5, ns/(2 * (max(ns) - min(ns))), vino)
    }
    
    list(vino = vino, conf = mdd, BAF = BAF, llr = llr)
}

CNVinput = function(nm, ns, geno, mm, ms, ss, nsize, ngeno, iig) {
    BAF = llr = rep(0, nsize)
    if (ngeno == 3) {
        c1 = sqrt(ss[[1]][1, 1])
        c2 = sqrt(ss[[2]][1, 1])
        c3 = sqrt(ss[[3]][1, 1])
        s1 = sqrt(ss[[1]][2, 2])
        s2 = sqrt(ss[[2]][2, 2])
        s3 = sqrt(ss[[3]][2, 2])
        tmp = nm >= mm[1]
        if (any(tmp)) {
            BAF[tmp] = 0
            llr[tmp] = log2(ns[tmp]/ns[1])
        }
        tmp = nm < mm[1] & nm > mm[2]
        if (any(tmp)) {
            k1 = (mm[1] - nm[tmp])/c1
            k2 = (nm[tmp] - mm[2])/c2
            BAF[tmp] = 0.5 * k1/(k1 + k2)
            llr[tmp] = log2(ns[tmp]/((k2 * ns[1] + k1 * ns[2])/(k1 + k2)))
        }
        tmp = nm <= mm[2] & nm > mm[3]
        if (any(tmp)) {
            k1 = (mm[2] - nm[tmp])/c2
            k2 = (nm[tmp] - mm[3])/c3
            BAF[tmp] = 0.5 + 0.5 * k1/(k1 + k2)
            llr[tmp] = log2(ns[tmp]/((k2 * ns[2] + k1 * ns[3])/(k1 + k2)))
        }
        tmp = nm <= mm[3]
        if (any(tmp)) {
            BAF[tmp] = 1
            llr[tmp] = log2(ns[tmp]/ns[3])
        }
    }
    if (ngeno == 2) {
        id = sum(iig)
        if (id == 4) {
            # 1,3
            c1 = sqrt(ss[[1]][1, 1])
            c3 = sqrt(ss[[3]][1, 1])
            s1 = sqrt(ss[[1]][2, 2])
            s3 = sqrt(ss[[3]][2, 2])
            tmp = nm >= mm[1]
            if (any(tmp)) {
                BAF[tmp] = 0
                llr[tmp] = log2(ns[tmp]/ns[1])
            }
            tmp = nm < mm[1] & nm > mm[3]
            if (any(tmp)) {
                k1 = (mm[1] - nm[tmp])/c1
                k3 = (nm[tmp] - mm[3])/c3
                BAF[tmp] = k1/(k1 + k3)
                llr[tmp] = log2(ns[tmp]/((k3 * ns[1] + k1 * ns[3])/(k1 + k3)))
            }
            tmp = nm <= mm[3]
            if (any(tmp)) {
                BAF[tmp] = 1
                llr[tmp] = log2(ns[tmp]/ns[3])
            }
        }
        if (id == 5) {
            # 2,3
            c2 = sqrt(ss[[2]][1, 1])
            c3 = sqrt(ss[[3]][1, 1])
            s2 = sqrt(ss[[2]][2, 2])
            s3 = sqrt(ss[[3]][2, 2])
            tmp = nm >= mm[2]
            if (any(tmp)) {
                BAF[tmp] = 0.5
                llr[tmp] = log2(ns[tmp]/ns[2])
            }
            tmp = nm < mm[2] & nm > mm[3]
            if (any(tmp)) {
                k2 = (mm[2] - nm[tmp])/c2
                k3 = (nm[tmp] - mm[3])/c3
                BAF[tmp] = 0.5 + 0.5 * k2/(k2 + k3)
                llr[tmp] = log2(ns[tmp]/((k3 * ns[2] + k2 * ns[3])/(k2 + k3)))
            }
            tmp = nm <= mm[3]
            if (any(tmp)) {
                BAF[tmp] = 1
                llr[tmp] = log2(ns[tmp]/ns[3])
            }
        }
        if (id == 3) {
            # 1,2
            c1 = sqrt(ss[[1]][1, 1])
            c2 = sqrt(ss[[2]][1, 1])
            s1 = sqrt(ss[[1]][2, 2])
            s2 = sqrt(ss[[2]][2, 2])
            tmp = nm >= mm[1]
            if (any(tmp)) {
                BAF[tmp] = 0
                llr[tmp] = log2(ns[tmp]/ns[1])
            }
            tmp = nm < mm[1] & nm > mm[2]
            if (any(tmp)) {
                k1 = (mm[1] - nm[tmp])/c1
                k2 = (nm[tmp] - mm[2])/c2
                BAF[tmp] = 0.5 * k1/(k1 + k2)
                llr[tmp] = log2(ns[tmp]/((k1 * ns[2] + k2 * ns[1])/(k1 + k2)))
            }
            tmp = nm <= mm[2]
            if (any(tmp)) {
                BAF[tmp] = 0.5
                llr[tmp] = log2(ns[tmp]/ns[2])
            }
        }
    }
    if (ngeno == 1) {
        if (iig == 1) {
            BAF = rep(0, nsize)
            llr = log2(ns/ns[1])
        }
        if (iig == 2) {
            BAF = rep(0.5, nsize)
            llr = log2(ns/ns[2])
        }
        if (iig == 3) {
            BAF = rep(1, nsize)
            llr = log2(ns/ns[3])
        }
    }
    list(BAF = BAF, llr = llr)
}
