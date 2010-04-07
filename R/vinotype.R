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
vinotype = function(nm,ns,geno){
  #=========== find centers
  iig = sort(unique(geno)); ngeno = length(iig);
  if( ngeno == 1) mms = median(ns)
  if( ngeno==3 ) mms = min(tapply( ns[geno==1|geno==3], geno[geno==1|geno==3], median)) 
  if(ngeno==2){ k = tapply( nm, geno, median );  g = match(geno, names(k[ k < -0.3 | k > 0.3 ]) )
        if( sum(!is.na(g)) > 1 )  mms = min( tapply( ns[!is.na(g)], geno[!is.na(g)], median) ) 
        else mms = median(ns)
  }
  mm = ms = NULL; adata = cbind(nm,ns);
  nsize = length(geno)
  vino = rep(-1, nsize); vino[(ns > mms)] = 2; vino1 = vino; ss = list(); m = matrix(0,2,2); lm = 0
  for(ik in iig){
    l = sum( geno == ik )
    if(l > 1){
      nns = ns[geno==ik]; nnm = nm[geno==ik]
      onn = order(nns); onns = nns[onn]; onnm = nnm[onn]; th = which( onns < mms )
      if( length(th) > 0 ){ 
        th = max(th); if( l > th ) th = th+1
        donns = diff(onns); b = boxplot(donns, range = 2, plot = FALSE); b$out = c(b$out, donns[donns>1.5])
        tb = boxplot(donns, range = 5, plot = FALSE)
        if( sum(tb$out > 0.2 )>0 ){
          tb$out = tb$out[ tb$out > 0.2 ]
          k=match(tb$out, donns); k = k[k<th ];
          if( length(k)>0  ){ vino[geno==2&ns <= onns[max(k)]] = 1}
        }  
        if( length(b$out) > 0 ){  k=match(b$out, donns); k = k[k<th ]; 
          if( length(k)>0  ){ k=max(k)+1; nns  = onns[k:l]; nnm  = onnm[k:l] } } 
      }
      mm = c(mm, mean(nnm));  ms = c(ms, mean(nns) )
      if( length(nnm) == 1 ) ss[[ik]] = NA
      else{ ss[[ik]] = cov( cbind(nnm, nns) ); m = m + ss[[ik]]; lm = lm+1 }
    }
    else{  mm = c(mm, nm[geno==ik]); ms = c(ms, ns[geno==ik]); ss[[ik]] = NA }
  }
  m = m/lm
  for(ik in iig){
    if( is.na( ss[ik] ) )  ss[[ik]] = m
    else ss[[ik]] = ss[[ik]]*0.5 + m*0.5
  }   
  #=========== test 2 vs. 3
  if( ngeno == 3 ){ # three groups; # lower BIC better fit.
    l1 = sum(geno==1); l2 = sum(geno==2); l3 = sum(geno==3);  mdd = rep(0,nsize)
    if(l1 > 1 & det(ss[[ 1 ]]) > 10^(-10) ){ 
      th=c(0.99,0.95, 0.99);
      tmp = mahalanobis(adata[geno==1,], c(mm[1], ms[1]), ss[[1]] )
      dd1 = pchisq(tmp, df=2) ; mdd[geno==1] = dd1 ; vino1[geno==1][dd1 < th[1]] = 2 }
    if(l2 > 1 & det(ss[[ 2 ]]) >10^(-10)){ 
      tmp = mahalanobis(adata[geno==2,], c(mm[2], ms[2]), ss[[2]] ); 
      dd2 = pchisq(tmp, df=2); mdd[geno==2] = dd2 ; vino1[geno==2][dd2 < th[2]] = 2}
    if(l3 > 1 & det(ss[[ 3 ]]) > 10^(-10)){ 
      tmp = mahalanobis(adata[geno==3,], c(mm[3], ms[3]), ss[[3]] ); 
      dd3 = pchisq(tmp, df= 2); mdd[geno==3] = dd3; vino1[geno==3][dd3 < th[3]] = 2}
  }
  if( ngeno == 2 ){
    l1 = sum(geno==iig[1]); l2 = sum(geno==iig[2]); mdd = rep(0,nsize); 
    th=c(0.99,0.95, 0.99); th = th[ iig ]
    if(l1 > 1 & det(ss[[ iig[1] ]]) > 10^(-10) ){ 
      tmp = mahalanobis(adata[geno==iig[1],], c(mm[1], ms[1]), ss[[iig[1]]] )
      dd1 = pchisq(tmp, df = 2) ; mdd[geno==iig[1]] = dd1 ; vino1[geno==iig[1]][dd1 < th[1]] = 2}
    if(l2 > 1 & det(ss[[ iig[2] ]]) > 10^(-10) ){ 
      tmp = mahalanobis(adata[geno==iig[2],], c(mm[2], ms[2]), ss[[iig[2]]] )
      dd2 = pchisq(tmp, df = 2) ; mdd[geno==iig[2]] = dd2; vino1[geno==iig[2]][dd2 < th[2]] = 2 }
  }
  if( ngeno == 1 ){
    tmp = mahalanobis( adata, c(mm[1], ms[1]), ss[[iig[1]]] )
    mdd = pchisq(tmp, df = 2); ; vino1[mdd < 0.99] = 2
  }
  alld = 1- pnorm( ns, mean = mean(ns[vino1==2]), sd = sd(ns[vino1==2]) )
  # vino = outlier with in each group & low intensity 
  k = alld*mdd 
  vino[ (vino1!=2 & k > .9999) ] =  1
  ons = order(ns); sn = ns[ons]; dsn = diff(sn); 
  if( any( dsn > 1.5 )){
    bdsn = dsn[dsn > 1.5]; g=match(bdsn, dsn); th = sum(ns < median(ns)); g=g[ g < th]
    vino[ ns <= sn[g] ] = 1 
  }
  if( any(vino==1) & any(vino == -1) ) vino = vdist(nm/5,ns/(2*(max(ns)-min(ns))), vino)
  # update the genotype membership after vino being removed
  #ug = sort( unique( geno[ vino != 1 ] ) )
  #if( length(ug) == 2 ){
  #  tmp = c(median( nm[ geno== ug[1] & vino != 1 ] ), median(nm[ geno == ug[2] & vino != 1]) )
  #  if( tmp[1] > .3 & tmp[2] < -.3 ) iig = c(1,3)
  #  else if( tmp[1] < .3 ) iig = c(2,3)
  #  else if( tmp[2] > -.3 ) iig = c(1,2)
  #  else iig = c(1,2)
  #  geno = iig[ match(geno, ug) ]
  #}
  list( vino = vino, conf = mdd )
}


