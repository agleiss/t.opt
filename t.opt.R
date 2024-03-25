t.opt <- function (                  # PERFORM PERMUTATION T-TEST AND OPT-TEST
  values,                            # data matrix containing log2-expression data for group 0 and group 1
  # each row corresponding to one gene each column to one sample (cell line)
  # rownames of values give gene-identifier
  groups,                            # vector containing 0s and 1s indicating group for each column in values
  type.all      =FALSE,              # max variant: optimize over no trimming or maximum trimming (alternative: "all")
  maxtrim       =2,                  # maximum number of cell lines trimmed in the smaller (!) group
  nperm         =50,                 # number of permuations for FDR estimation and calculation of permutation p-values
  seed          =010772,             # seed for producing permutations of group labels
  fdr.thresh    =0.05,               # p-value threshold for FDR estimation
  fc.thresh     =NA,                 # optional fold change threshold for FDR estimation
  roundcalc     =5,                  # digit after comma at which additional internal information can be stored in expression values
  s0            =c(0,0),             # add constant to denominator in t- and opt-statistic to obtain a SAM-like moderated test statistic
  called.from.s0= FALSE,             # if called from s0.opt, dont echo and dont permute
  trim.method   ="one.sided"
)
{         
  
  if(trim.method=="two.sided")stop("Two-sided trimming and testing currently not supported.")
  echo<-TRUE
  if(called.from.s0) echo<-FALSE
  nrow<-dim(values)[1]
  if(echo) cat("\nNumber of genes selected:  ",nrow)
  
  value0  <-round(values[,(groups==0)],roundcalc)
  value1  <-round(values[,(groups==1)],roundcalc)
  nsamp0  <-dim(value0)[2]
  nsamp1  <-dim(value1)[2]
  row.id  <-rownames(value0)
  
  if (nsamp0<=nsamp1) {                 # maxtrim for smaller group
    if (type.all) trim0=0:maxtrim
    else trim0=c(0,maxtrim)
    trim1=round(trim0*nsamp1/nsamp0)
  }
  if (nsamp0>nsamp1) {
    if (type.all) trim1=0:maxtrim
    else trim1=c(0,maxtrim)
    trim0=round(trim1*nsamp0/nsamp1)
  }
  trims<-rbind(trim0,trim1)
  colnames(trims)<-1:length(trim0)
  rownames(trims)<-c("group 0: ","group 1: ")
  if(echo) cat("\n\nNumbers of trimmed observations:\n")
  if(echo) print(trims)
  
  add<-matrix(NA,nrow=nrow,ncol=nsamp0)
  for (z in 1:nrow) add[z,]<-sample(0:(nsamp0-1))*10^(-roundcalc-2)
  value0.add<-value0+add                  # sorgt für eindeutige, aber zufällige zeilenweise Minima, wird in test.trim berücksichtigt
  value1.add<-value1+add+10^(-roundcalc-1)                  # sorgt für eindeutige zeilenweise Minima
  value.add<-cbind(value0.add,value1.add)
  
  ### original data
  res0<-test.trim(mat0=value0.add,mat1=value1.add,trim0=trim0,trim1=trim1,roundcalc=roundcalc,s0=s0)
  
  ### permutation test and FDR estimation
  
  if (!is.na(seed)) set.seed(seed)
  count.t<-0
  count.topt<-0
  sign.p.t<-NA*1:nperm
  sign.fc.t<-NA*1:nperm
  sign.p.topt<-NA*1:nperm
  sign.fc.topt<-NA*1:nperm
  if (!is.na(fdr.thresh)) {
    pperms.t<-matrix(NA,nrow=nrow,ncol=nperm)
    pperms.topt<-matrix(NA,nrow=nrow,ncol=nperm)
  }
  
  if(!called.from.s0){  
    cat("\nPermutation: ")   
    
    for (p in 1:nperm) {
      if (round(p/(nperm/10))==p/(nperm/10)) cat(p," ")
      rand.col<-sample(1:(nsamp0+nsamp1))
      value0.rand<-value.add[,rand.col[1:nsamp0]]
      value1.rand<-value.add[,rand.col[(nsamp0+1):(nsamp0+nsamp1)]]
      res<-test.trim(mat0=value0.rand,mat1=value1.rand,trim0=trim0,trim1=trim1,roundcalc=roundcalc,s0=s0)
      count.t<-count.t+1*(res$p.t<=res0$p.t)
      count.topt<-count.topt+1*(res$p.topt<=res0$p.topt)
      sign.p.t[p]<-sum(res$p.t<=fdr.thresh)
      sign.fc.t[p]<-sum(res$fc.t>=fc.thresh)
      sign.p.topt[p]<-sum(res$p.topt<=fdr.thresh)
      sign.fc.topt[p]<-sum(res$fc.topt>=fc.thresh)
      pperms.t[,p]<-res$p.t
      pperms.topt[,p]<-res$p.topt
    }
    cat("\nPermutations completed")   
    pperm.t<-(count.t+1)/(nperm+1)
    pperm.topt<-(count.topt+1)/(nperm+1)
    
    cat("\nCalculate pi0")   
    pi0.p.t<-     (sum(res0$p.t>=median(pperms.t)) / (0.5*nrow))
    pi0.p.topt<-  (sum(res0$p.topt>=median(pperms.topt)) / (0.5*nrow))
    cat("\nCalculate FDR")   
    fdr.p.t<-     mean(sign.p.t)    * pi0.p.t     / sum(res0$p.t<=fdr.thresh) 
    fdr.p.topt<-  mean(sign.p.topt) * pi0.p.topt  / sum(res0$p.topt<=fdr.thresh)
    cat("\nCalculate q-values for t-test")   
    f_p0_p.t<- apply(matrix(res0$p.t,length(res0$p.t),1),1,function(x,a)mean(x>=a),pperms.t) 
    f_p_p.t<- apply(matrix(res0$p.t,length(res0$p.t),1),1,function(x,a)mean(x>=a),res0$p.t) 
    q.t <- pi0.p.t* f_p0_p.t / f_p_p.t
    cat("\nCalculate q-values for topt-test")   
    f_p0_p.topt<- apply(matrix(res0$p.topt,length(res0$p.topt),1),1,function(x,a)mean(x>=a),pperms.topt) 
    f_p_p.topt<- apply(matrix(res0$p.topt,length(res0$p.topt),1),1,function(x,a)mean(x>=a),res0$p.topt) 
    q.topt <- pi0.p.topt* f_p0_p.topt / f_p_p.topt
    cat("\nCalculate SEs")   
    se.sign.p.t<-     (sd(sign.p.t)/sqrt(nperm))    * pi0.p.t     / sum(res0$p.t<=fdr.thresh)
    se.sign.p.topt<-  (sd(sign.p.topt)/sqrt(nperm)) * pi0.p.topt  / sum(res0$p.topt<=fdr.thresh)
    
    cat("\nCalculate significances")   
    if (is.na(fc.thresh)) {
      sign.t<-(q.t<=fdr.thresh)
      sign.topt<-(q.topt<=fdr.thresh)
      cat("\n\nFDR threshold of",fdr.thresh,":\n")
      cat("t-test:     ", sum(sign.t),"genes declared significant\n")
      cat("t.opt-test: ", sum(sign.topt),"genes declared significant\n")
    }
    else {
      sign.t<-(q.t<=fdr.thresh & res0$fc.t>=fc.thresh)
      sign.topt<-(q.topt<=fdr.thresh & res0$fc.topt>=fc.thresh)  
      cat("\n\nFDR threshold of",fdr.thresh," and FC threshold of",fc.thresh,":\n")
      cat("t-test:     ", sum(sign.t),"genes declared significant\n")
      cat("t.opt-test: ", sum(sign.topt),"genes declared significant\n")
    }
    
    rownames(res0)<-row.id
    erg<-list(t=res0$t, topt=res0$topt, q.t=q.t, q.topt=q.topt, trim0=res0$trim0.topt, 
              trim1=res0$trim1.topt, pi0.t=pi0.p.t[1], pi0.topt=pi0.p.topt[1], 
              fc.t=res0$fc.t, fc.topt=res0$fc.topt, s0.t=s0[1], s0.topt=s0[2], p.t=res0$p.t, p.topt=res0$p.topt,
              sign.t=sign.t, sign.topt=sign.topt)
  }
  else {
    erg<-list(t=res0$t, topt=res0$topt)
  }
  erg
}  # end t.opt

############################################

min.trim <- function(x, na.rm=F, trim=1) {
  
  if (trim>1) {
    for (i in 1:(trim-1)) {
      x<-x[x>min(x,na.rm=na.rm)]
    }
  }
  min(x, na.rm=na.rm)
}

##################################

trim.matrix <- function(mat, trim=1) {
  # setzt eindeutige (Minimal-)Werte voraus!!
  if (trim>0) {
    matmax<-max(mat,na.rm=T)*10
    mat[is.na(mat)]<-matmax+runif(sum(is.na(mat))) # notwendig ab R-version 2.9.1
    nr<-dim(mat)[1]  
    nc<-dim(mat)[2]
    mat.trim<-apply(mat,1,min.trim,trim=trim,na.rm=T) # NA's werden ignoriert
    mat.minbig<-rep(mat.trim,rep(nc,nr))  # aufgeblasener Vektor mit zeilenweisen Minima
    mat.vec<-as.vector(t(mat))
    isnot.trim<-(mat.vec>mat.minbig)# (gelöst: ACHTUNG: Ties ergeben mehr als ein FALSE in isnot.trim)
    # ACHTUNG: NA ergibt in isnot.trim dann auch NA
    mat.trim1<-mat.vec[isnot.trim] # NA's werden übernommen => passt
    
    if((nc-trim)*nr != length(mat.trim1)) {
      isnot.trim.mat<-1*isnot.trim
      dim(isnot.trim.mat)<-c(nc,nr)
      search<-apply(isnot.trim.mat,2,sum)
      print(search)
      
    }
    
    dim(mat.trim1)<-c((nc-trim),nr)
    mat.trim1<-t(mat.trim1)
    mat.trim1[mat.trim1>=matmax]<-NA
    
    return(mat.trim1)
  }
  else return(mat)
}

###################################

count <- function (x) {
  
  xx <- !is.na(x)
  sum(xx)
}


###################################

test.trim <- function (mat0, mat1, trim0, trim1, roundcalc=100, s0=c(0,0)) {
  
  # only one-sided!!!
  
  nrow<-dim(mat0)[1]
  n.trim<-length(trim0) # =length(trim1)
  mean0<- mean1<- sd0<- sd1<- t.ttest<- p.ttest<-topt.ttest<- popt.ttest<-matrix(NA,nrow=nrow,ncol=n.trim)
  
  for (i in 1:n.trim) {
    mat0.trim<-trim.matrix(mat=mat0,trim=trim0[i])
    mat1.trim<-trim.matrix(mat=mat1,trim=trim1[i])
    mat0.trim<-round(mat0.trim,roundcalc)
    mat1.trim<-round(mat1.trim,roundcalc)
    n0<-apply(mat0.trim,1,count) # berücksichtige NAs und Trimming
    n1<-apply(mat1.trim,1,count)
    mean0[,i]<-apply(mat0.trim,1,sum,na.rm=T)/n0 # schneller als mean!
    mean1[,i]<-apply(mat1.trim,1,sum,na.rm=T)/n1
    sd0[,i]<-sqrt(apply((mat0.trim-mean0[,i])^2,1,sum,na.rm=T)/(n0-1))     # schneller als sd!
    sd1[,i]<-sqrt(apply((mat1.trim-mean1[,i])^2,1,sum,na.rm=T)/(n1-1))     # schneller als sd!
    s<-sqrt(((n1-1)*sd1[,i]^2+(n0-1)*sd0[,i]^2) / (n1+n0-2))
    s.t<-s+s0[1]                       # add constant to denominators
    s.topt<-s+s0[2]                    # add constant to denominators
    t.ttest[,i]<-(mean1[,i]-mean0[,i])/(s.t*sqrt((n1+n0)/(n1*n0)))
    p.ttest[,i]<-1-pt(t.ttest[,i],df=n0+n1-2)
    topt.ttest[,i]<-(mean1[,i]-mean0[,i])/(s.topt*sqrt((n1+n0)/(n1*n0)))
    popt.ttest[,i]<-1-pt(topt.ttest[,i],df=n0+n1-2)
  }
  p.t<-p.ttest[,1]
  which.topt<-apply(popt.ttest,1,which.min) # nimmt bei Gleichheit ersten Wert als Minimum => passt
  trim0.topt<-trim0[which.topt]
  trim1.topt<-trim1[which.topt]
  fc.t<-2^(mean1[,1]-mean0[,1])
  sel<-cbind(1:nrow,which.topt)
  p.topt<-popt.ttest[sel]
  fc.topt<-2^(mean1[sel]-mean0[sel])
  
  data.frame(trim0.topt,trim1.topt,p.t,p.topt,fc.t,fc.topt,t=t.ttest[,1],topt=topt.ttest[sel])
  
}

############################################


read.abi <- function (                      # READING AND PRE-FILTERING ABI DATASETS
  indat0,                            # CSV-file with log2-expression data for group 0
  indat1,                            # CSV-file with log2-expression data for group 1
  select.gene   =0,                  # optional indicator for selecting subgroup of genes
  maxgene       =100000,             # optional restriction to first maxgene genes
  col.value0    =10+7*seq(0,9),      # column numbers corr. to assay normalized signal values in indat0 file (default for ABI)
  col.value1    =10+7*seq(0,9),      # column numbers corr. to assay normalized signal values in indat1 file (default for ABI)
  flag.value0   =13+7*seq(0,9),      # column numbers corr. to quality flags in indat0 file (default for ABI)
  flag.value1   =13+7*seq(0,9),      # column numbers corr. to quality flags in indat1 file (default for ABI)
  sn.value0     =9+7*seq(0,9),       # column numbers corr. to signal-to-noise ratios in indat0 file (default for ABI)
  sn.value1     =9+7*seq(0,9),       # column numbers corr. to signal-to-noise ratios in indat1 file (default for ABI)
  col.gene      =1,                  # column number for unique gene identifier (default for ABI), if 0 then sequential numbers are defined
  flag.thresh   =5000,               # threshold for quality flags (below="good")
  flag.filter   =0.333,              # if non-negative then genes with a proportion of probes with flags<flag.thresh smaller than or equal to flag.filter are deleted from the datase
  flag.na       =TRUE,               # remove expression values corresponding to flag>=flag.thresh
  sn.thresh     =3,                  # threshold for signal-to-noise ratio
  sn.filter     =0.333              # if non-negative then genes with a proportion of probes with |S/N|>=sn.thresh smaller than or equal to sn.filter are deleted from the dataset
)
{          
  ### read data
  
  
  if (select.gene==0) {
    orig.data0<-read.csv2(indat0,nrows=maxgene)
    orig.data1<-read.csv2(indat1,nrows=maxgene)
  }
  else {
    orig.data0<-read.csv2(indat0)[select.gene,]
    orig.data1<-read.csv2(indat1)[select.gene,]
  }
  ngene   <-dim(orig.data0)[1]
  nsamp0  <-length(col.value0)
  nsamp1  <-length(col.value1)
  nrow    <-min(maxgene,ngene)
  value0  <-log2(orig.data0[1:nrow,col.value0])
  value1  <-log2(orig.data1[1:nrow,col.value1])
  flags   <-1*cbind((orig.data0[1:nrow,flag.value0]<flag.thresh),(orig.data1[1:nrow,flag.value1]<flag.thresh))
  sns     <-1*cbind((abs(orig.data0[1:nrow,sn.value0])>=sn.thresh),(abs(orig.data1[1:nrow,sn.value1])>=sn.thresh))
  if (col.gene>0)
    row.id <-orig.data0[1:nrow,col.gene]
  else
    row.id<-1:nrow
  
  ### prefiltering
  
  if (flag.na) {
    sel0<-(orig.data0[1:nrow,flag.value0]>=flag.thresh)
    sel1<-(orig.data1[1:nrow,flag.value1]>=flag.thresh)
    value0[sel0]<-NA
    value1[sel1]<-NA
  }       
  mean.flags<-(apply(flags,1,sum)/(nsamp0+nsamp1)>flag.filter)        # faster than mean
  mean.sns<-(apply(sns,1,sum)/(nsamp0+nsamp1)>sn.filter)
  both<-(mean.flags & mean.sns)
  
  value0.prefilt<-value0[both,]
  value1.prefilt<-value1[both,]
  row.id<-row.id[both]
  
  cat("\nOriginal number of genes selected:  ",nrow)
  cat("\nNumber of genes after pre-filtering:",1*sum(both))
  cat("\n\n")
  
  rownames(value0.prefilt)<-row.id
  rownames(value1.prefilt)<-row.id
  out<-list(values=cbind(value0.prefilt,value1.prefilt), groups=rep(c(0,1),c(nsamp0,nsamp1)))
  
} # end read.abi

#################################


#optimize s0 for topt AND t
s0.opt<-function(values=genexpr, groups=group, bins=NA, ...){
  #  if(is.na(bins)) bins<-floor(nrow(values)/20) # at least 20 genes in each bin
  if(is.na(bins)) bins<-min(floor(nrow(values)/20),100) # at least 20 genes in each bin
  mean0<-apply(values[,groups==0],1,mean)
  mean1<-apply(values[,groups==1],1,mean)
  n1<-sum(groups==0)
  n2<-sum(groups==1)
  sd0<-sqrt(apply((values[,group==0]-mean0)^2,1,sum,na.rm=T)/(n1-1))     # schneller als sd!
  sd1<-sqrt(apply((values[,group==1]-mean1)^2,1,sum,na.rm=T)/(n2-1))     # schneller als sd!
  s<-sqrt(((n2-1)*sd1^2+(n1-1)*sd0^2) / (n1+n2-2))
  hist(s)
  cvmad.i.t<-NULL
  cvmad.i.topt<-NULL
  for(i in seq(-0.025,1,0.025)){ 
    if(i == -0.025) use.s0<-0     else use.s0<-quantile(s,min(i,1))
    stats1<-  t.opt(values=values,groups=groups, s0=c(use.s0,use.s0),called.from.s0=TRUE, ...)
    #print(cvmad.i.t)
    cvmad.i.t<-rbind(cvmad.i.t,cbind(use.s0,cvmad(stats1$t,s,bins=bins) ))
    cvmad.i.topt<-rbind(cvmad.i.topt,cbind(use.s0,cvmad(stats1$topt,s,bins=bins) ))
  }
  print(cvmad.i.t)
  print(cvmad.i.topt)
  mincvmad.t<-min(cvmad.i.t[,2])
  mincvmad.topt<-min(cvmad.i.topt[,2])
  c(cvmad.i.t[cvmad.i.t[,2]==mincvmad.t,1],cvmad.i.topt[cvmad.i.topt[,2]==mincvmad.topt,1])
}

cv<-function(values){
  sd(values)/mean(values)
}

cvmad<-function(t,s, bins=100,plot=FALSE){
  s<-s[rank(s)]
  t<-t[rank(s)]
  
  binsize<-length(t)/100
  mad.bin<-rep(0,bins)
  s.bin<-rep(0,bins)
  for(i in 1:bins){
    indices<-round((i-1)*binsize+1,0):round(i*binsize,0)
    mad.bin[i]<-mad(t[indices])
    s.bin[i]<-median(s[indices])
  }
  if(plot)plot(s.bin,mad.bin,lty=1)
  #print(mad.bin)
  cv(mad.bin)
}

#################################

