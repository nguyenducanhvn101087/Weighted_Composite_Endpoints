# VERSION 2.0 

# ANH NGUYEN DUC - OXFORD UNIVERSITY CLINICAL RESEARCH UNIT

# SUPERVISED BY:  MARCEL WOLBERS
# CREATED:        DEC  09  2014
# LAST MODIFIED:  JUL  22  2015

#/////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

# FUNCTIONS INCLUDED:

# wchibarsqa
# pchibarsqa
# rchibarsqa               
# qchibarsqa


#/////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

require(varComp)
require(Matrix)
require(gdata)
require(quadprog)
#-----------------------------------------------------------------------------------------------

wchibarsqa <- function (V, R1=NULL, R2=NULL, ...) {
# Input:
# - V        Covariance matrix of the statistic
# - R1, R2   The matrices defining the cone, representing the linear inequalities and equalities 
#            respectively. See Shapiro 1988: Towards a Unified Theory of Inequality Constrained 
#            Testing in Multivariate Analysis.
#            In detail: the cone is defined by R1 x w >=0 & R2 x w = 0, where w is a vector in R^m
# - ...      Extra arguments for rankMatrix
  
# Output:
# - wts      the weights for the chi-bar-square distribution
  
  tryCatch({ # To catch unprecedented errors    
    
    if (det(V) == 0) {
      cat("V must be nonsingular!\n")
      return(NULL)
      # stop("V must be nonsingular!")
    } # end of if (det(V) == 0)
    
    
    if (any(eigen(V)$values <= 0)) {
      cat("V must be positive definite!\n")
      return(NULL)
      # stop("V must be positive definite!")
    } # end of if (any(eigen(V)$values <= 0))
    
    
    nr1 <- nrow(R1)
    nc1 <- ncol(R1)
    
    nr2 <- nrow(R2)
    nc2 <- ncol(R2)
    
    
    if (is.null(R1) & is.null(R2)) { # similar to when R1=diag(nr1), i.e nonnegativity
      
      wts <- wchibarsq(solve(V))  
      
    } # end of if (is.null(R1) & is.null(R2))
    
    
    # If the cone is defined by ONLY linear inequality constraints
    if (!is.null(R1) & is.null(R2)) {
      
      if (nr1 == nc1) {
      
        I1 <- diag(nr1)
        
        if (max(abs(R1-I1))==0) { # similar to when both R1 & R2 are null
          
          wts <- wchibarsq(solve(V))
          
        } else {
          
          if (det(R1)!=0) {
            
            wts <- wchibarsq(R1 %*% solve(V) %*% t(R1))
            
          } else {
            cat("R1 must be nonsingular!\n")    
            return(NULL)
            # stop("R1 must be nonsingular!")    
          } # end of if (det(R1)!=0) else 
          
        } # end of if (max(abs(R1-I1))==0) else ...
           
        if(sum(wts)!=1) wts[length(wts)] <- 1-sum(wts[seq_len(length(wts)-1)])
               
      } else if (nr1 < nc1) {
        
        rank1 <- rankMatrix(R1, ...)
        
        if (rank1 == nr1) {
          
          wtmp <- wchibarsq(R1 %*% solve(V) %*% t(R1))
          
          if(sum(wtmp)!=1) wtmp[length(wtmp)] <- 1-sum(wtmp[seq_len(length(wtmp)-1)])
          
          wts  <- rep(0, nrow(V)+1)
          
          wts[(nc1-nr1+1):(nc1+1)] <- wtmp[1:(nr1+1)]
          
        } else {
          cat("R1 must be of full row rank\n")
          return()
        } # end of if (rank1 == nr1) else ...
        
      } else {
        cat("The number of rows of R1 must not be larger than its number of columns\n")
        return(NULL)
        # stop("The number of rows of R1 must not be larger than its number of columns")
      } # end of if (nr1==nc1) else ...
            
    } # end of if (!is.null(R1) & is.null(R2))
    
    
    # If the cone is defined by a number of linear inequality as well as equality constraints
    if (!is.null(R1) & !is.null(R2)) {
      
      k <- nr1 + nr2 # total number of rows in R1 and R2
      
      if (k > nc1 | nc1 != nc2) { # This condition was verified
        cat("R1 and R2 must have the same number of columns and the sum of their row dimensions must not be larger than their common column dimension!\n")
        return(NULL)
        # stop("R1 and R2 must have the same number of columns and the sum of their row dimensions must not be larger than their common column dimension!")
        
      } else {
        
        R   <- rbind(R1, R2) # equivalently t(cbind(t(R1), t(R2)))
        raR <- rankMatrix(R, ...)
          
        if (raR != k) {
          cat("The matrix combined of rows of R1 and R2 must be of full row rank!\n")
          return(NULL)
          # stop("The matrix combined of rows of R1 and R2 must be of full row rank!")
        } else {
          
          Z <- solve(R %*% solve(V) %*% t(R))[1:nr1, 1:nr1]
          
          wtmp <- wchibarsq(solve(Z))
          
          if(sum(wtmp)!=1) wtmp[length(wtmp)] <- 1-sum(wtmp[seq_len(length(wtmp)-1)])
          
          wts  <- rep(0, nrow(V)+1)
          
          wts[(nc1-k+1):(nc1-nr2+1)] <- wtmp[1:(nr1+1)]
          
        }
        
      } # end of if (k > nc1 | nc1 != nc2)
      
    } # end of if (!is.null(R1) & !is.null(R2))
    
    
    if (is.null(R1) & !is.null(R2)) { # Only linear equality
      # The resulting cone is the linear subspace of R^(nc2)
      # The dimension of this subspace is called the nulity of R2, denoted by k and: rank(R2) + k = nc2
      # Ref http://en.wikipedia.org/wiki/Kernel_(linear_algebra)
      
      raR2 <- rankMatrix(R2, ...) 
      k    <- nc2 - raR2
      
      if (k == 0) {
        cat("The cone defined by R2 x w = 0 degenerates into a single vector 0!\n")
        return(NULL)
        # stop("The cone defined by R2 x w = 0 degenerates into a single vector 0!")
      } else {
        
        wts      <- rep(0, nrow(V)+1)
        wts[k+1] <- 1
        
      } # end of if (k == 0) else ...
          
    } # end of if (is.null(R1) & !is.null(R2))
    
    attributes(wts) <- NULL           
    
    wts
    
  }, error=function(ee) {cat("From function: wchibarsqa\n"); ee; browser()})      
} # end of wchibarsqa
#-----------------------------------------------------------------------------------------------

# # test------------ 
# require(clusterGeneration)
# set.seed(69)
# m <- 5
# V <- genPositiveDefMat(m)$Sigma
# 
# R1 <- rbind(c(1,0,0,0,0), c(0,0,0,0,1))
# R2 <- rbind(c(1,0,0,0,0), c(0,0,0,0,1))
# wchibarsqa(V, R1, R2, method="qrLINPACK") # expect error, correct
# 
# R1 <- rbind(c(0,1,0,0,0), c(0,0,0,1,0))
# wchibarsqa(V, R1, R2, method="qrLINPACK", warn.t=F) # correct
# 
# R2 <- matrix(0, nrow=1, ncol=5)
# wchibarsqa(V, NULL, R2, method="qrLINPACK", warn.t=F) # Scheffe, correct
# 
# R2 <- diag(m)
# wchibarsqa(V, NULL, R2, method="qrLINPACK", warn.t=F) # Degenerate -> expect error, correct
#-----------------------------------------------------------------------------------------------

# Anh's pchibarsq allowing for precomputed wts
pchibarsqa <- function (q, V=NULL, wts=NULL, R1=NULL, R2=NULL, lower.tail = TRUE, log.p = FALSE, ...) {
tryCatch({ # To catch unprecedented errors    
  if(is.null(wts)) {
    wts <- wchibarsqa(V, R1, R2, ...)
    if(is.null(wts)) return(NULL)
    
    n <- nrow(V)
  } else{
    n <- length(wts)-1
  }
  ans <- pchisq(q, 0, lower.tail = FALSE) * wts[1L] + pchisq(q, 
                                                            n, lower.tail = FALSE) * wts[n + 1L]
  for (i in seq_len(n - 1)) {
    ans <- ans + pchisq(q, i, lower.tail = FALSE) * wts[i + 1L]
  }
  ans[q <= 0] <- 1
  ans <- if (isTRUE(lower.tail)) 
    1 - ans
  else ans
  if (isTRUE(log.p)) 
    log(ans)
  else ans
}, error=function(ee) {cat("From function: pchibarsqa\n"); ee; browser()})      
} # end of pchibarsqa
#-----------------------------------------------------------------------------------------------

# Anh's rchibarsq allowing for precomputed wts
rchibarsqa <- function(n, wts=NULL, V=NULL, R1=NULL, R2=NULL, ...) {
tryCatch({ # To catch unprecedented errors  
  if(is.null(wts)) {
    wts <- wchibarsqa(V, R1, R2, ...)
    if(is.null(wts)) return(NULL)
    
    nw  <- nrow(V)+1
  } else {
    nw  <- length(wts)
  } # end of if(is.null(wts)) else ...
  
  ids <- sample(x=seq_len(nw), size=n, replace=T, prob=wts)
  
  rchisq(n, df=ids-1)
}, error=function(ee) {cat("From function: rchibarsqa\n"); ee; browser()})    
} # end of rchibarsqa

# # test------------ 
# m   <- 4
# WW  <- crossprod(matrix(rnorm(m^2), m))
# cbs <- sort(rchibarsqa(1000, V=WW))
# pw  <- pchibarsq(cbs, WW)   
# plot(ecdf(cbs))   
# lines(cbs, pw, col=4, lwd=3, lty=3)
# # alles gut!
#-----------------------------------------------------------------------------------------------

# Anh's qchibarsq allowing for precomputed wts based on either resampling or uniroot (exact but longer)
qchibarsqa <- function(p, wts=NULL, V=NULL, R1=NULL, R2=NULL, method="resampling", B=2e4, lower.tail=T, ...) {
# method      either resampling or uniroot
# B           number of resamples, only needed for resampling method  
# ...         futher arguments for quantile or uniroot (causing no issue at all, tested)
  
tryCatch({ # To catch unprecedented errors  
  if (is.null(wts)) {
    wts <- wchibarsqa(V, R1, R2, ...)
    if(is.null(wts)) return(NULL)
    
    nw  <- nrow(V)+1
  } else {
    nw  <- length(wts)
  } # end of if (is.null(wts)) else ...
  
#   args <- list(...)
#   
#   na.rm <- names <- type <- NULL # essential pars for quantiles
#   ss <- NULL # essential pars for uniroot
#   
#   for(i in 1:length(args)) {
#     assign(names(args)[i], args[[i]])
#   }
#   
  if (method=="resampling") {
    x <- sort(rchibarsqa(B, wts=wts))
    q <- quantile(x, probs=p, ...)
  } else if (method=="uniroot") {
    # ix    <- sort(x, index.return=T, decreasing=T) # the lower for uniroot is always 0
    # y     <- x[ix]
    q <- sapply(p, function(x) {
      upper <- qchisq(x, df=nw-1)
      uniroot(f=function(z) {x - pchibarsqa(q=z, wts=wts)}, interval=c(0, upper), ...)$root
    })
    
  } else {
    cat("Unsupported method!\n")
  } # end of if (method=="resampling") else ...
  
  q
}, error=function(ee) {cat("From function: qchibarsqa\n"); ee; browser()})
} # end of qchibarsq
#------------------------------

# # test sampling-based quantile 
# alp <- .05
# m   <- 4
# WW  <- crossprod(matrix(rnorm(m^2), m))
# wts <- wchibarsq(V=WW)
# attributes(wts) <- NULL   
# 
# system.time(qs  <- qchibarsqa(1-alp, wts=wts, method="resampling")) # less precise, faster
# qs
# pchibarsqa(q=qs, wts=wts)
# 
# system.time(qr  <- qchibarsqa(1-alp, wts=wts, method="uniroot")) # more precise, slower
# qr 
# pchibarsqa(q=qr, wts=wts)
#-----------------------------------------------------------------------------------------------
