# VERSION 2.0 

# ANH NGUYEN DUC - OXFORD UNIVERSITY CLINICAL RESEARCH UNIT

# SUPERVISED BY:  MARCEL WOLBERS
# CREATED:        DEC  09  2014
# LAST MODIFIED:  JUN  08  2015

#/////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

# FUNCTIONS INCLUDED:

# covmat.multinom
# w.prop.test
# power.w.prop.test
# n.w.prop.test

#/////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

covmat.multinom <- function(p){
  # covariance matrix of a multinomial distribution
  cov.p <- -matrix(p,ncol=1)%*%matrix(p,nrow=1) 
  diag(cov.p) <- p*(1-p)
  cov.p # not yet consider n. Actually this is correct but the name is wrong. It should be called
  # a categorical distribution, a multinomial counterpart for bernouli and binomial
} # end of covmat.multinom
#-----------------------------------------------------------------------------------------------

w.prop.test <- function(x1,x2,n1,n2,w){
  # Wald test for comparison of weighted proportions
  # x1: observed failure counts for each component in group 1 (vector of length J)
  # x2: observed failure counts for each component in group 2 (vector of length J)
  # n1,n2: sample size of group 1 or 2 (positive integers)
  # w: component weights (vector of length J)
  #-----------------------------------------------------------
  p1 <- x1/n1
  p2 <- x2/n2
  #- weighted estimate
  diff.w <- sum(w*(p1-p2))
  #- standard error of diff.w
  se.diff.w  <- sqrt(matrix(w,nrow=1)%*%
                       (covmat.multinom(p1)/n1+covmat.multinom(p2)/n2)%*%
                       matrix(w,ncol=1))
  # result
  list(diff.w=diff.w,
       lower.95ci=diff.w-qnorm(0.975)*se.diff.w,
       upper.95ci=diff.w+qnorm(0.975)*se.diff.w,
       z=diff.w/se.diff.w,
       p=2*pnorm(-abs(diff.w/se.diff.w)) # if this is one sided then we multiply by 2??
  )  
} # end of w.prop.test
#-----------------------------------------------------------------------------------------------

power.w.prop.test <- function(p1, p2, n, w, sig.level = .05, alternative="two.sided") {
  # Power of Wald test for comparing of weighted proportions
  # (at the one-sided sig.level significance level to detect a decrease or increase OR
  #  at the tow-sided sig.level significance level to detect a difference in either directions)
  # p1,p2: failure proportions for each component in group 1 and 2 (vectors of length J)
  # n: sample size in each group (assumed to be equal, i.e. n=n1=n2)
  #- target treatment effect
  diff.w <- sum(w*(p1-p2)) # weighted expected difference
  #- standard error of diff.w under H0 and HA
  se.diff.w.H0  <- sqrt(matrix(w,nrow=1)%*%
                          (2*covmat.multinom((p1+p2)/2)/n)%*%
                          matrix(w,ncol=1))
  se.diff.w.HA  <- sqrt(matrix(w,nrow=1)%*%
                          ((covmat.multinom(p1)+covmat.multinom(p2))/n)%*%
                          matrix(w,ncol=1))
  #- derive power
  if (alternative=="one.sided") {
    power <- pnorm((diff.w-qnorm(1-sig.level)*se.diff.w.H0)/se.diff.w.HA)
  } else if (alternative=="two.sided") {
    power <- pnorm((diff.w-qnorm(1-sig.level/2)*se.diff.w.H0)/se.diff.w.HA) +  
             pnorm((diff.w+qnorm(1-sig.level/2)*se.diff.w.H0)/se.diff.w.HA, lower.tail=F)
  } else {
    print("alternative must be either one.sided or two.sided!")  
  }# end of if (alternative=="one.sided") else ...
  
  power    
  
} # hopefully now work for p2 >= p1 (i.e. larger sample size larger power) as well as two sided H_A
#-----------------------------------------------------------------------------------------------

n.w.prop.test <- function(p1, p2, alpha=.05, alternative="two.sided", power, w=NULL) {
  if(is.null(w)) w <- rep(1, length.out=length(p1)) / length(p1)
  foo <- function(n) {
    power.w.prop.test(p1,p2,n,w,alpha,alternative) - power
  }
  n <- uniroot(foo, c(2, 1e7))$root
} # end of n.w.prop.test
#-----------------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------------
