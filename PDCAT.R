### This R code computes the p-value of PDCAT   
### Four steps in our method
### Firstly, get the true coefficients after projection in these functions pst and pst_dis
### Secondly, get two test statistics by the direction of multiple coefficients in this dst
### Thirdly, get  the paraments of empirical distribution in this chi_paraments
### Finally, combine these p-values of all test statistics in these PTstat and PDCAT

#correlation matrix of genotype
acor <- function(geno_data){
  cor <- cor(as.matrix(geno_data))
  if(is.positive.definite(cor)==FALSE){
    cor <- make.positive.definite(cor, 1.5)
  }
  return(cor)
}

# true coefficient after projection
pst <- function(geno_data, maf, Wald){
  acor <- acor(geno_data)
  s_Z <- sqrt(2*maf*(1-maf))
  # a projection matrix R_star
  R_star <- acor %*% diag(s_Z)
  # transform Wald test statistics into true coefficients after a projection matrix
  RZ <- ginv(R_star) %*% Wald
  return(RZ=RZ)
}

# normal distribution form of true coefficient
pst_dis <- function(geno_data, maf){
  cor <- cor(as.matrix(geno_data))
  acor <- acor(geno_data)
  s_Z <- sqrt(2*maf*(1-maf))
  # a projection matrix R_star
  R_star <- acor %*% diag(s_Z)
  # normal distribution covariance matrix of true coefficient and standard deviation of all elements
  cor_RZ <- t(ginv(R_star)) %*% cor %*% (ginv(R_star))
  s_RZ <- sqrt(diag(cor_RZ))
  return(list(s_RZ=s_RZ,cor_RZ=cor_RZ))
}

# divide multiple true coefficients into two statistics by direction of all coefficients
dst <- function(RZ, s_RZ, lambdaSet) {
  ## divide them by direction of standardized form of multiple coefficients
  U <- RZ / s_RZ
  Uneg <- U[U < 0]
  Upos <- U[U > 0]
  Tneg <- rep(NA, length(lambdaSet))
  Tpos <- rep(NA, length(lambdaSet))
  # consider all situations of sign
  for (k in 1:length(lambdaSet)) {
    if (length(Uneg) > 0 && length(Upos) > 0) { 
      if (lambdaSet[k] == Inf) {
        Tneg[k] <- max(abs(Uneg))
        Tpos[k] <- max(abs(Upos))
      } else {
        Tneg[k] <- sum(Uneg^lambdaSet[k])
        Tpos[k] <- sum(Upos^lambdaSet[k])
      }
    } else {
      if (length(Upos) == 0 ) { 
        if (lambdaSet[k] == Inf) {
          Tneg[k] <- max(abs(Uneg))
        } else {
          Tneg[k] <- sum(Uneg^lambdaSet[k])
        }
      } else{
        if (lambdaSet[k] == Inf) {
          Tpos[k] <- max(abs(Upos))
        } else {
          Tpos[k] <- sum(Upos^lambdaSet[k])
        }
      }  
    }
  }
  statVal <- c(Tneg, Tpos)
  ## return the value of dst for the candidate value of lambda
  return(statVal)
}

# calculate parameter of the empirical distribution
chi_parameter <- function(s_RZ, cor_RZ, lambdaSet, B){
  ### calculate the p-value using generalized chi-square distribution
  lamNum <- length(lambdaSet)
  dst.boot = matrix(NA,nrow = B,ncol = 2*lamNum)
  num_snps <- length(s_RZ)
  RZ.boot <- mvrnorm(B, rep(0, num_snps), cor_RZ)
  for(j in 1:B){
    dst.boot[j,] = dst(RZ=RZ.boot[j, ], s_RZ=s_RZ, lambdaSet=lambdaSet)
  }
  # Convert NA values in dst.boot and dst to zero
  dst.boot[is.na(dst.boot)] <- 0
  dst.boot_square <- (dst.boot)^2
  # use the moment estimate to fit the generalized chi-square distribution
  K <- rep(NA, 2*lamNum)
  L <- rep(NA, 2*lamNum)
  M <- rep(NA, 2*lamNum)
  a <- rep(NA, 2*lamNum)
  b <- rep(NA, 2*lamNum)
  d <- rep(NA, 2*lamNum)
  for (i in (1:(2*lamNum))){
    K[i] <- mean(dst.boot_square[,i])
    L[i] <- (1/B)*sum((dst.boot_square[,i]-K[i])^2)
    M[i] <- (1/B)*sum((dst.boot_square[,i]-K[i])^3)
    a[i] <- (M[i])/(4*L[i])
    b[i] <- (K[i])-2*(L[i])^2/(M[i])
    d[i] <- 8*(L[i])^3/(M[i])^2
  } 
  v1 <- c(paste0('aneg.', lambdaSet), paste0('apos.', lambdaSet))
  names(a) <- v1
  v2 <- c(paste0('bneg.', lambdaSet), paste0('bpos.', lambdaSet))
  names(b) <- v2
  v3 <- c(paste0('dneg.', lambdaSet), paste0('dpos.', lambdaSet))
  names(d) <- v3
  return(c(a,b,d))
}

# calculate p-values of statistics by these parameters
PTstat <- function(RZ, s_RZ, cor_RZ, a, b, d, lambdaSet, B) {
  lamNum <- length(lambdaSet)
  dst_val <- rep(NA, 2 * lamNum)
  dst_val <- dst(RZ = RZ, s_RZ = s_RZ, lambdaSet = lambdaSet)
  dst_val[is.na(dst_val)] <- 0
  dst_square <- dst_val^2
  PTstat <- rep(NA, (2 * lamNum))
  PTstat_pval <- rep(NA, (2 * lamNum))
  PTstat <- (dst_square - b) / a
  PTstat_pval <- pchisq(PTstat, df = d, lower.tail = FALSE)
  #if the dimension of gene_data is extremely small, it may cause all PTstat to be smaller than 0 
  # repeating bootstrap until one of PTstat is not smaller than 0 can solve this phenomenon
  while (all(PTstat < 0)) {
    chi_parameter <- chi_parameter(s_RZ, cor_RZ, lambdaSet, B)
    PTstat <- (dst_square - chi_parameter$b) / chi_parameter$a
    PTstat_pval <- pchisq(PTstat, df = chi_parameter$d, lower.tail = FALSE)
  }
  names(PTstat) <- c(paste0('neg.', lambdaSet), paste0('pos.', lambdaSet))
  names(PTstat_pval) <- c(paste0('neg.', lambdaSet), paste0('pos.', lambdaSet))
  return(c(PTstat = PTstat, PTstat_pval = PTstat_pval))
}

#Combine these p-values use the Cauchy combination test
PDCAT <- function(PTstat, PTstat_pval){
  lamNum <- length(lambdaSet)
  for (i in 1 : (2*lamNum)){
    if((PTstat[i] > 0) & (PTstat_pval[i] == 1)){
      PTstat_pval[i] <-PTstat_pval[i]-0.000001
    }
  }
  if (all(PTstat > 0)){
    t_com <- matrix(NA, nrow =lamNum, ncol = 2)
    for (i in 1:lamNum) {
      t_com[i, ] <- c(PTstat_pval[i], PTstat_pval[i + lamNum])
    }
    t <- rep(NA, lamNum)
    PTstat_lamdma <- rep(NA, lamNum)
    # Aggregated Cauchy association test
    for (i in 1:lamNum) {
      t[i] <- (1/2) * sum(tan((0.5 - t_com[i,]) * pi) )
      PTstat_lamdma[i] <- 0.5 - (atan(t[i])) / pi
    }
    PDCAT_stat <- (1/lamNum) * sum(tan((0.5 - PTstat_lamdma) * pi))
    PDCAT_pval <- 0.5 - (atan(PDCAT_stat)) / pi  
  }else if(sum(PTstat < 0) < 4){
    for (i in 1:lamNum) {
      if (any(PTstat[c(i, i + lamNum)] < 0)) {
        PTstat_pval[i] <- 1
        PTstat_pval[i+lamNum] <- 1
      }
    }
    # Delete elements at PTstat < 0 in PTstat_pval
    t_eta <- rep(NA, sum(PTstat_pval != 1))
    t_eta <- PTstat_pval[PTstat_pval != 1]
    num_t <- (1/2) * length(t_eta)  
    t_com<- matrix(NA, nrow = num_t, ncol = 2)
    for (i in 1:num_t) {
      t_com[i, ] <- c(t_eta[i], t_eta[i + num_t])
    }
    T_eta <- rep(NA, num_t)
    PTstat_lamd <- rep(NA, num_t)
    # Aggregated Cauchy association test
    for (i in 1:num_t) {
      T_eta[i] <- (1/2) * sum(tan((0.5 - t_com[i,]) * pi) )
      PTstat_lamd[i] <- 0.5 - (atan(T_eta[i])) / pi
    }
    PDCAT_stat <- (1/num_t) * sum(tan((0.5 - PTstat_lamd) * pi))
    PDCAT_pval <- 0.5 - (atan(PDCAT_stat)) / pi 
  }else{
    # Delete elements at PTstat < 0 in PTstat_pval
    PTstat_lamd <- rep(NA, sum(PTstat_pval != 1))
    PTstat_lamd <- PTstat_pval[PTstat_pval != 1]
    PDCAT_stat <- (1/length(PTstat_lamd)) * sum(tan((0.5 - PTstat_lamd) * pi))
    PDCAT_pval <- 0.5 - (atan(PDCAT_stat)) / pi 
  }
  names(PDCAT_pval) <- "PDCAT_pval" 
  return(PDCAT_pval) #obtain the final p-value of PDCAT
}
