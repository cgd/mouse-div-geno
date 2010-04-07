#============================================== read the cel files
# Definition : genotype the array. Used inside of genotypthis() function. However, given contrast and average intensity, 
# the function can be used to genotype each SNP. 
#
# input :
# nm : Contrast ( A allele intensity - B allele intensity ) of one probe set.
# ns : Average of one probe set.
# hint : mean and variance of each group provided by MouseDivGeno package. genotype() invoke the hint file, but if this function
#     is used separated, hint can be NULL. 
# trans : transformation method to obtain the contrast intensity
#
# output : 
# genotype : 1 = AA, 2 = AB, 3 = BB
# <example>
# data('MM'); data('SS');
# nm = MM[1,]; ns = SS[1,]
# geno = genotype(nm, ns, hint = NULL) 
#======================================================================
genotype = function(nm, ns, hint1, trans){ # hint : postm$ncenter, postm$nscales, posts$ncenter, posts$nscales 123 456 789 101112
  hint = rep(0,3)
  hint[1] = max(nm); hint[3] = min(nm)
  if( is.na(hint[2]) )  hint[2] = 0 
  else hint[2] = hint1 

  nsize = length(nm)
  thres = c(-0.4, -0.3, -0.15, 0.15, 0.3, 0.4) # for CCStrans 
  MAtransthres = c( -0.8, -0.6, -0.3, 0.3, 0.6, 0.8 ) 
  if( is.element( trans, 'MAtrans' ) ) thres = MAtransthres
  rmid =  nm >= thres[2] & nm <= thres[5] & ns < quantile(ns, .20); lnrmid = sum(!rmid)
  istwo = FALSE; iig = NULL
  #================================================ obvious one group SNPs
  if( all( nm[!rmid] < thres[1] ) ){ geno = rep(3, nsize); v1 = vinotype(nm, ns, geno); vino = v1$vino; conf = v1$conf }
  else if(all(nm[!rmid] >  thres[6] ) ){ geno = rep(1, nsize); v1 = vinotype(nm, ns, geno); vino = v1$vino; conf = v1$conf }
  else if(all(nm[!rmid] >= thres[3] & nm[!rmid] <= thres[4])){ geno = rep(2, nsize); v1 = vinotype(nm, ns, geno); vino = v1$vino; conf = v1$conf }
  else if( all( nm[!rmid] <= thres[2] | nm[!rmid] >= thres[5] ) ){ 
      geno= rep(1, nsize); geno[nm<0]=3; v2 = vinotype(nm, ns, geno); vino = v2$vino; conf = v2$conf } # 1,3
  else{
    ##============================================ obvious two groups SNPs 
    if( all( nm <= thres[5] ) ){ istwo = TRUE; iig = c(2,3) } # 2,3
    if( all( nm >= thres[2] ) ){ istwo = TRUE; iig = c(1,2) } # 1,2
    adata = cbind(nm, ns/(max(ns)-min(ns)))
    nsize = length(nm)
    #======== test 2 ====== based on one dimension at this time rather than two
    otheta1 = theta1 = list( tau = c(.5, .5),
              mu1 = hint[1], mu2 = hint[3], sigma1 = 0.01, sigma2 = 0.01 )
    ok = TRUE; delta = .001; imax = 50; iter = 0
    while(ok & iter <= imax){
      iter = iter + 1
      T <- E.step2(theta1,nm[!rmid]) ; T[is.na(T)] = 0.5
      theta1 <- M.step2(T, nm[!rmid])
      ok1 = TRUE
      if(theta1$tau[1] > 0) ok1 = ok1 * (abs(theta1$mu1[1]-otheta1$mu1[1])<delta)
      if(theta1$tau[2] > 0) ok1 = ok1 * (abs(theta1$mu2[1]-otheta1$mu2[1])<delta)
      if( ok1 | is.na(theta1$sigma1) |  is.na(theta1$sigma2) | is.infinite( theta1$sigma1 ) |  is.infinite(theta1$sigma2) ) ok = FALSE
      if( ok ){ m=(theta1$sigma1+theta1$sigma2)/2; theta1$sigma1=0.5*theta1$sigma1+0.5*m; theta1$sigma2=0.5*theta1$sigma2+0.5*m }   
      otheta1 = theta1
    }
    T = round(T, 4)
    tgeno2 = rep(-1,nsize); geno2 = rep(-1, lnrmid )
    ig = order( theta1$tau, decreasing = TRUE )
    ii = ig[1];  geno2[ T[,ii] >= median( T[T[,ig[2]]<.5,ii] ) ] = ii 
    ii = ig[2];  tt = table(T[ geno2 == -1,ii]); stt = sort(tt, decreasing = TRUE)[1]
    th = sort( as.numeric(names(tt[tt == stt])) , decreasing = TRUE)[1]  
    geno2[ geno2 == -1 & T[,ii] >=  min( max(T[,ii]), max(median( T[T[,ig[1]]<.5,ii] ), th)) ] = ii 
    tgeno2[!rmid] = geno2
    if( length( unique( tgeno2[ tgeno2 != -1 ] ) ) == 1 ){ 
      istwo = TRUE; mm = median(nm)
      if( mm > thres[5] ) geno = rep(1, nsize)
      else if( mm < thres[2] ) geno = rep(3, nsize)
      else geno = rep(2, nsize) 
      v1 = vinotype(nm, ns, geno); vino = v1$vino; conf = v1$conf
    }
    else{
      if( any(tgeno2== -1) ) tgeno2 = vdist(adata[,1], adata[,2]/5, tgeno2)
      tscore2 = silhouette( match( tgeno2[!rmid],unique(tgeno2[!rmid])),  dist(nm[!rmid]) )
      if( istwo ){geno = test12(tscore2, nm, tgeno2, rmid, thres, iig, nsize);  v1 = vinotype(nm, ns, geno); vino = v1$vino; conf = v1$conf}
    }
    if( !istwo ){ #======== test 3
      otheta1 = theta1 = list(tau = c(1/3,1/3,1/3), 
        mu1 = hint[1], mu2 =hint[2], mu3 = hint[3],sigma1 = 0.1 , sigma2 = 0.1, sigma3 = 0.1 )
      ok = TRUE; delta = .001; imax = 50; iter = 0
      while(ok & iter <= imax){
        iter = iter + 1
        T <- E.step3(theta1, nm[!rmid])  ; T[is.na(T)] = 1/3 
        theta1 <- M.step3(T, nm[!rmid]) 
        ok1 = TRUE
        if(theta1$tau[1] > 0 ) ok1 = ok1 * (abs(theta1$mu1[1]-otheta1$mu1[1])<delta)
        if(theta1$tau[2] > 0 ) ok1 = ok1 * (abs(theta1$mu2[1]-otheta1$mu2[1])<delta)
        if(theta1$tau[3] > 0 ) ok1 = ok1 * (abs(theta1$mu3[1]-otheta1$mu3[1])<delta)
        if( ok1 | is.na(theta1$sigma1)  |  is.na(theta1$sigma2) |  is.na(theta1$sigma3) | is.infinite( theta1$sigma1 ) 
             |  is.infinite(theta1$sigma2) |  is.infinite(theta1$sigma3) ) ok = FALSE
        if( ok ){ m=(theta1$sigma1+theta1$sigma2+theta1$sigma3)/3; theta1$sigma1=0.5*theta1$sigma1+0.5*m; theta1$sigma2=0.5*theta1$sigma2+0.5*m;
                 theta1$sigma3=0.5*theta1$sigma3+0.5*m }   
        otheta1 = theta1
      }
      T = round(T, 4)
      tgeno3 = rep(-1, nsize); geno3 = rep(-1, lnrmid )
      ig = order( theta1$tau, decreasing = TRUE )
      mm = matrix( c(3,2,1, 1,3, 2, 1, 2, 3), nrow = 3, byrow = TRUE); mm = mm[ig,]; out1 = apply(T, 2, max)
      out = rep(1,3)
      out[ mm[1,3] ] = median(T[ T[,mm[1,1]] < .5 & T[,mm[1,2]] < .5, mm[1,3] ]); 
      tt = table(T[ geno3 == -1,mm[2,3]]); stt = sort(tt, decreasing = TRUE)[1]
      th = sort( as.numeric(names(tt[tt == stt])) , decreasing = TRUE)[1] 
      out[ mm[2,3] ] = 
       max(median(T[ T[,mm[2,1]] < min(out[mm[2,1]],.5) & T[,mm[2,2]] < min(out[mm[2,2]],.5), mm[2,3] ]),th)
      tt = table(T[ geno3 == -1,mm[3,3]]); stt = sort(tt, decreasing = TRUE)[1]
      th = sort( as.numeric(names(tt[tt == stt])) , decreasing = TRUE)[1] 
      out[ mm[3,3] ] = max(median(T[T[,mm[3,1]] < min(out[mm[3,1]],.5) & T[,mm[3,2]] < min(out[mm[3,2]],.5), mm[3,3] ]), th) 
      if(any(is.na(out))){geno = test12(tscore2, nm, tgeno2, rmid, thres, iig, nsize);  v1 = vinotype(nm, ns, geno); vino = v1$vino; conf = v1$conf}
      else{
        out[ out1 < out ] = out1[ out1 < out]
        geno3[ T[,mm[1,3] ] >= out[mm[1,3]] ] = mm[1,3] ; geno3[ T[,mm[2,3] ] >= out[mm[2,3]] ] = mm[2,3] ;
        geno3[ T[,mm[3,3] ] >= out[mm[3,3]] ] = mm[3,3] ;
        tgeno3[!rmid] = geno3
        # need confidence score here 
        if( length( unique( tgeno3[ tgeno3 != -1 ] ) ) < 3 ){
           geno = test12(tscore2, nm, tgeno2, rmid, thres, iig, nsize);  v1 = vinotype(nm, ns, geno); vino = v1$vino; conf = v1$conf
        }
        else{
          if( any(tgeno3== -1) ) tgeno3 = vdist(adata[,1],adata[,2]/5,tgeno3) 
          tscore3 = silhouette(match(tgeno3[!rmid],unique(tgeno3[!rmid])), dist(nm[!rmid]))
          geno = test123(tscore2, tscore3, nm, tgeno2, tgeno3, rmid, thres, iig, nsize)
          v3 = vinotype(nm, ns, geno); vino = v3$vino; conf = v3$conf
        }
      }
    } # !istwo
  } # all genotype end
  list(geno = geno, vino = vino, conf = conf)
}

genotypeSexchr = function(nm, ns, trans){  # always only two groups
  nsize = length(nm); hint = c(max(nm), min(nm))
  thres = c(-0.4, -0.3, -0.15, 0.15, 0.3, 0.4) # for CCStrans 
  MAtransthres = c( -0.8, -0.6, -0.3, 0.3, 0.6, 0.8 ) 
  if( is.element( trans, 'MAtrans' ) ) thres = MAtransthres
  rmid =  nm >= thres[2] & nm <= thres[5] & ns < quantile(ns, .20); lnrmid = sum(!rmid)

  if( all( nm[!rmid] < thres[1] ) ){ geno = rep(3, nsize); v1 = vinotype(nm, ns, geno); vino = v1$vino; conf = v1$conf }
  else if( all( nm[!rmid] >  thres[6] ) ){ geno = rep(1, nsize); v1 = vinotype(nm, ns, geno); vino = v1$vino; conf = v1$conf }
  else if( all( nm[!rmid] <= thres[3] | nm[!rmid] >= thres[4] ) ){ 
      geno= rep(1, nsize); geno[nm<0]=3; v2 = vinotype(nm, ns, geno); vino = v2$vino; conf = v2$conf } # 1,3
  else{
    iig = c(1,3)
    adata = cbind(nm, ns/(max(ns)-min(ns)))
    #======== test 2
    otheta1 = theta1 = list( tau = c(.5, .5),
              mu1 = hint[1], mu2 = hint[2], sigma1 = 0.01, sigma2 = 0.01 )
    ok = TRUE; delta = .001; imax = 50; iter = 0
    while(ok & iter <= imax){
      iter = iter + 1
      T <- E.step2(theta1,nm[!rmid])  
      theta1 <- M.step2(T, nm[!rmid])
      ok1 = TRUE
      if(theta1$tau[1] > 0) ok1 = ok1 * (abs(theta1$mu1[1]-otheta1$mu1[1])<delta)
      if(theta1$tau[2] > 0) ok1 = ok1 * (abs(theta1$mu2[1]-otheta1$mu2[1])<delta)
      if( ok1 | is.na(theta1$sigma1) |  is.na(theta1$sigma2) | is.infinite( theta1$sigma1 ) |  is.infinite(theta1$sigma2) ) ok = FALSE
      if( ok ){ m=(theta1$sigma1+theta1$sigma2)/2; theta1$sigma1=0.5*theta1$sigma1+0.5*m; theta1$sigma2=0.5*theta1$sigma2+0.5*m }   
      otheta1 = theta1
    }
    T = round(T, 4)
    geno = rep(-1,nsize); geno2 = rep(-1, lnrmid )
    ig = order( theta1$tau, decreasing = TRUE )
    ii = ig[1];  geno2[ T[,ii] >= median( T[T[,ig[2]]<.5,ii] ) ] = ii 
    ii = ig[2];  tt = table(T[ geno2 == -1,ii]); stt = sort(tt, decreasing = TRUE)[1]
    th = sort( as.numeric(names(tt[tt == stt])) , decreasing = TRUE)[1]  
    geno2[ geno2 == -1 & T[,ii] >=  min( max(T[,ii]), max(median( T[T[,ig[1]]<.5,ii] ), th)) ] = ii 
    geno[!rmid] = geno2
    if( length( unique( geno[ geno != -1 ] ) ) == 1 ){ 
      mm = median(nm)
      geno = rep(1, nsize)
      if( mm < 0 ) geno = rep(3, nsize)
      v1 = vinotype(nm, ns, geno); vino = v1$vino; conf = v1$conf
    }
    else{
      if( any(geno== -1) ) geno = vdist(adata[,1], adata[,2]/2, geno)
      geno[ geno == 2 ] = 3
      v1 = vinotype(nm, ns, geno); vino = v1$vino; conf = v1$conf
    }
  } # else
  list(geno = geno, vino = vino, conf = conf)
}


test123 = function(tscore2, tscore3, nm, tgeno2, tgeno3, rmid, thres, iig, nsize){
  mscore2 = mean(tscore2[,3]); mscore3 = mean(tscore3[,3])
  if( length(tscore3) < 3 ) mscore3 = 0 
  if( length(tscore2) < 3 ) mscore2 = 0 
  d1 = quantile( nm[tgeno3==1 & !rmid ], .1) - quantile(nm[tgeno3==2 & !rmid ], .9)  
  d2 = quantile( nm[tgeno3==2 & !rmid ], .1) - quantile(nm[tgeno3==3 & !rmid ], .9)
  if( ( (mscore2 - mscore3 <= 0.1) | mscore3 > 0.8 ) & (d1 > 0.1 & d2 > 0.1 )  ) 
           geno = tgeno3 # three groups 
  else geno = test12( tscore2, nm, tgeno2, rmid, thres, iig, nsize ) 
  geno
}

test12 = function( tscore2, nm, tgeno, rmid, thres, iig, nsize){
  mscore = mean(tscore2[,3])
  if( length(tscore2) < 3 ) mscore = 0 
  d = quantile( nm[ tgeno==1 & !rmid ], .1) - quantile(nm[ tgeno==2 & !rmid ], .9)
  if( mscore > .6 & d > .1 ){ # two group 
    if( length(iig) == 0 ){
      tmp = c(median( nm[tgeno==1]), median(nm[tgeno==2]) )
      if( tmp[1] > thres[5] & tmp[2] < thres[2] ) iig = c(1,3)
      else if( tmp[1] < thres[5] ) iig = c(2,3)
      else iig = c(1,2)
    }
    tgeno = iig[match( tgeno, sort(unique(tgeno)) ) ]
  }
  else{
    mm = median(nm)
    if( mm > thres[5] ) tgeno = rep(1, nsize)
    else if( mm < thres[2] ) tgeno = rep(3, nsize)
    else tgeno = rep(2, nsize)
  }
  tgeno 
}

vdist = function(ta,tb,t2){
  la = which(t2 == -1); l = length(la)  
  dta = as.matrix(dist(cbind(ta, tb))) ; diag(dta) = 1000
  while( l > 0 ){ 
    t3 = t2[ t2 != -1 ]; tmp = dta[la, t2 != -1 ];  
    if( l == 1 ){
      t2[ la ]= t3[ which(min(tmp) == tmp)[1] ]   
    }
    else{
      id = which(tmp == min(tmp))[1]; iid = id %% l
      iid[ iid == 0 ] = l
      t2[ la[iid] ] = t3[ which(tmp[iid,] == min(tmp[iid,]))[1] ]
    }
    la = which(t2 == -1); l = length(la)  
  }
  t2
}
E.step3 <- function(theta,data){
  i1 = i2 = i3 = rep(0, length(data))
  if(theta$tau[1]>0) i1=theta$tau[1]*dnorm(data,mean=theta$mu1,sd=sqrt(theta$sigma1))
  if(theta$tau[2]>0) i2=theta$tau[2]*dnorm(data,mean=theta$mu2,sd=sqrt(theta$sigma2))
  if(theta$tau[3]>0) i3=theta$tau[3]*dnorm(data,mean=theta$mu3,sd=sqrt(theta$sigma3))
  t(apply(cbind(i1,i2,i3),1,function(x) x/sum(x)))
}

M.step3 <- function(T,data){
  mu1 = mu2 = mu3 = sigma1 = sigma2 = sigma3 = NA
  tau= apply(T,2,mean)
  if(tau[1] > 0){
    mu1= weighted.mean(data,T[,1]);
    sigma1= cov.wt(matrix(data,ncol=1),T[,1])$cov
  }
  if(tau[2] > 0){
    mu2= weighted.mean(data,T[,2]);
    sigma2= cov.wt(matrix(data,ncol=1),T[,2])$cov
  }
  if(tau[3] > 0 ){
   mu3= weighted.mean(data,T[,3])
   sigma3= cov.wt(matrix(data,ncol=1),T[,3])$cov
  }
  list(tau = tau, mu1 = mu1, mu2 = mu2, mu3=mu3, sigma1 = sigma1,
       sigma2 = sigma2, sigma3 = sigma3)
}

E.step2 <- function(theta,data){
  i1 = i2 = rep(0, length(data))
  if(theta$tau[1]>0) i1=theta$tau[1]*dnorm(data,mean=theta$mu1,sd=sqrt(theta$sigma1))
  if(theta$tau[2]>0) i2=theta$tau[2]*dnorm(data,mean=theta$mu2,sd=sqrt(theta$sigma2))
  t(apply(cbind(i1,i2),1,function(x) x/sum(x)))
}

M.step2 <- function(T,data){
  mu1 = mu2 = sigma1 = sigma2 = NA
  tau= apply(T,2,mean)
  if(tau[1] > 0){
    mu1= weighted.mean(data,T[,1]);
    sigma1= cov.wt(matrix(data,ncol=1),T[,1])$cov
  }
  if(tau[2] > 0){
    mu2= weighted.mean(data,T[,2]);
    sigma2= cov.wt(matrix(data,ncol=1),T[,2])$cov
  }
  list(tau = tau, mu1 = mu1, mu2 = mu2, sigma1 = sigma1, sigma2 = sigma2)
}
    