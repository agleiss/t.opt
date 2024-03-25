# Adaptive trimmed one-sided t-test (t.opt)

Implements one-sided adaptive trimming for micro-array data

## Usage:

     t.opt(values, groups, type.all=FALSE, maxtrim=2, nperm=50, seed=010772,
       fdr.thresh.p=0.05, fdr.thresh.fc=NA, roundcalc=5, trim.method   ="one.sided")
      

## Arguments:

  values: data matrix containing log2-expression data for group 0 and
          group 1. Each row corresponds to one gene each column to one
          sample (cell line). Rownames of 'values' give
          gene-identifier.

  groups: vector containing 0s and 1s indicating group for each column 
          in 'values'

type.all=FALSE: specifies if the t.opt test optimizes over no trimming
          or maximum trimming given by 'maxtrim' ('type.all=FALSE', the
          default)  or over all variants in between.

maxtrim=2: maximum number of observations trimmed in the smaller group.
          The maximum trimming for the larger group is calculated
          proportional to the sample sizes.

nperm=50: number of permutations for FDR estimation and calculation of 
          permutation p-values.

seed=010772: seed for producing permutations of group labels.

fdr.thresh=0.05: threshold on FDR for declaring significance of a gene.

fc.thresh=NA: fold change threshold for declaring signficance of a
          gene.

roundcalc=5: digit after comma at which additional internal information
           can be stored in expression values.

s0=c(0,0): Stabilization constant for topt/t test statistics.

called.from.s0=FALSE: option to indicate recursive call from within
          s0.opt(). Suppresses output and permutation.

trim.method="one-sided": One or two-sided tests. Currently only
          one-sided tests are supported (Alternative hypothesis: gene
          expression higher in group 1)

## Details:

     Genes with predominantly (but not homogeneously) high expression
     in the  experimental group as compared to a control group exhibit
     a high difference between  the average expression in the control
     group and the average expression among  highly expressed samples
     of the experimental group. To estimate this difference, it is 
     straightforward to compare trimmed means. One way to choose the
     trimming proportion  adaptively is based on a series of
     comparisons of trimmed means. This is done using a  t-statistic
     where trimmed empirical means and trimmed empirical variances of
     the two groups are plugged in. For this test statistic the
     trimming proportion is selected from a pre-determined set, e. g.,
     Gamma = {0.0, 0.1, 0.2, 0.3}.  Raw P-values are computed for each
     trimming proportion by comparing the value of the test statistic
     to the t-distribution with n_x(1-gamma)+ n_y(1-gamma) - 2 degrees
     of freedom. We define that trimming proportion as the 'optimal'
     one which  yields the smallest raw P-value. Of course, the
     optimization of the trimming  proportion has to be accounted for
     when generating the permutation distribution  of the test
     statistic: the same series of tests is performed for each
     permutation  and the minimum of the obtained P-values is selected
     as the P-value corresponding  to that permutation. Obviously,
     there is a trade-off between the number of tests to be corrected
     for and the chance of finding that trimming proportion that best 
     detects PHE. Because of the inherent optimization with respect to
     the trimming proportion, we denote this test as the optimization
     test (Opt-test).

## Value:

     The object returned has the following attributes: 

       t: t statistics for each gene.

    topt: opt-t statistics for each gene.

     q.t: q-value from t-test for each gene.

  q.topt: q-value from Opt-test for each gene.

   trim0: set of number of observations to be trimmed in group 0.

   trim1: set of number of observations to be trimmed in group 1.

   pi0.t: Proportion of genes not related with group estimated from
          t-test.

pi0.topt: Proportion of genes not related with group estimated from
          Opt-test.

    fc.t: fold change calculated as 2^{mean difference}.

 fc.topt: fold change calculated as 2^{difference of trimed means}.

    s0.t: stabilizing constant for t-test (taken from input).

 s0.topt: stabilizing constant for Opt-test (taken from input).

     p.t: approximate p-value for t-test.

  p.topt: approximate minimally selected p-value for Opt-test.

  sign.t: significance indicator for each gene, from t-test.

sign.topt: significance indicator for each gene, from Opt-test.


## Examples:

     # generate data

     p<-500          # total number of genes evaluated
     genes.hhe<-100  # number of genes homogenously higher expressed
     genes.phe<-50  # number of genes predominantly higher expressed
     omega<-0.2     # predominance factor (fraction of samples in experimental group not higher expressed)
     n1<-10         # number of samples in control group
     n2<-10         # number of samples in experimental group
     shift<-6       # shift in gene expression between control and experimental group (for HHE and PHE genes)
     random.shift<-TRUE # whether shift is uniformly distributed between 0 and shift (TRUE) or the same for all HHE+PHE genes (FALSE)

     group<-c(rep(0,n1),rep(1,n2))
     higher<-matrix(rep(rbinom(n2*genes.phe,1,1-omega),genes.phe*n2),genes.phe, n2) # indicator whether tumor sample is not higher expressed
     if(random.shift) shift<-matrix(rep(runif(n2*genes.hhe,0,shift),genes.hhe*n2),genes.hhe, n2)
     genexpr<-matrix(rnorm(p*(n1+n2),5,1),p,n1+n2)
     genexpr[1:(genes.hhe),group==1]<-genexpr[1:(genes.hhe),group==1]+shift

     if(random.shift) shift<-matrix(rep(runif(n2*genes.phe,0,shift),genes.phe*n2),genes.phe, n2)
     genexpr[(genes.hhe+1):(genes.hhe+genes.phe),group==1]<-genexpr[(genes.hhe+1):(genes.hhe+genes.phe),group==1]+shift*higher

     meanintens<-apply(genexpr,1,mean)
     sdfudge<-meanintens
     genexpr<-genexpr/sdfudge

     aa<-t.opt(values=genexpr,groups=group,nperm=200)
     #FDR threshold of 0.05 :
     #t-test:      112 genes declared significant
     #t.opt-test:  109 genes declared significant

     s0<-s0.opt(genexpr,group)

     bb<-t.opt(values=genexpr,groups=group,nperm=200, s0=s0)
     #FDR threshold of 0.05 :
     #t-test:      116 genes declared significant
     #t.opt-test:  113 genes declared significant

     cc<-t.opt(values=genexpr,groups=group,nperm=200, fc.thresh=1.5, s0=s0)
     #FDR threshold of 0.05  and FC threshold of 1.5 :
     #t-test:      12 genes declared significant
     #t.opt-test:  21 genes declared significant


## Reference:

     Gleiss A, Sanchez-Cabo F, Perco P, Tong D, Heinze G (2010).
     Adaptive trimmed t-Tests  for identifying predominantly high
     expression in a microarray experiment. 
     Statist. Med. 2011, 30 52-61. DOI: 10.1002/sim.4093