
### Simulation
### num_pheno : number of phenotype datasets
### num_samples : sample size
### num_snps : number of SNP in a gene
### non_zero_snps : number of effect sizes of SNPs on phenotype is not 0
### r : correlation between SNPs
### b1, b2 : the effect size range (r > 0)
### b3, b4 : the effect size range (r < 0)

# Install packages
my_packages = c("MASS", "COMBAT", "corpcor", "mvtnorm")                             
not_installed = my_packages[!(my_packages %in% installed.packages()[ , "Package"])]          
if(length(not_installed)) install.packages(not_installed)    

library(MASS)
library(COMBAT)
library(corpcor)
library(mvtnorm)
# First generate Wald test statistics

# generate SNP
num_pheno <- 1000
num_samples <- 1000
num_snps <- 50
non_zero_snps <- 0
r=0.5
# if non_zero_snps != 0, we can set effect sizes(b1, b2, b3, b4).
# non_zero_snps <- 25
# b1 <- 0.055
# b2 <- 0.065
# b3 <- 0.05
# b4 <- 0.06 

  maf <- runif(num_snps, 0.05, 0.5)
  geno_data <- matrix(0, nrow = num_samples, ncol = num_snps)
  cov_matrix <- matrix(0, nrow = num_snps, ncol = num_snps)
  for (i in 1:num_snps) {
    for (j in 1:num_snps) {
      cov_matrix[i, j] <- (r) ^ abs(i - j)
    }
  }
  multivar_data <- mvrnorm(num_samples, maf, cov_matrix)
  for (i in 1:num_snps) {
    geno_data[, i] <- rbinom(num_samples, 2, plogis(multivar_data[, i]))
  }
  
# generate phenotype
  y <- vector("list", num_pheno)
  if(r > 0){
    for (i in 1:num_pheno) {
      b <- rep(0, num_snps)
      if (non_zero_snps > 0) {
        non_zero_indices <- sample(1:(num_snps - non_zero_snps + 1), 1)
        if (non_zero_snps %% 2 == 0) {            
          sign <- rep(c(-1, 1), length.out = non_zero_snps)
        } else {
          sign <- c(rep(c(-1, 1), length.out = non_zero_snps - 1), -1)
        }
        b[non_zero_indices:(non_zero_indices + non_zero_snps - 1)] <- runif(non_zero_snps, b1, b2) * sign
      }
      e <- rnorm(num_samples)
      y[[i]] <- as.matrix(geno_data) %*% b + e
    }
  }else{
    for (i in 1:num_pheno) {
      b <- rep(0, num_snps)
      if (non_zero_snps > 0) {
        non_zero_indices <- sample(1:(num_snps - non_zero_snps + 1), 1)
        positive_snps <- floor((4/5) * non_zero_snps)
        negative_snps <- non_zero_snps - positive_snps
        b[non_zero_indices:(non_zero_indices + positive_snps - 1)] <- runif(positive_snps, b3, b4)
        b[(non_zero_indices + positive_snps):(non_zero_indices + non_zero_snps - 1)] <- runif(negative_snps, -b4, -b3)
      }
      e <- rnorm(num_samples)
      y[[i]] <- as.matrix(geno_data) %*% b + e
    }
  }

# GWAS univariate linear regression to get Wald test statistics
z_statistics <- matrix(NA, nrow = num_pheno, ncol = num_snps * 2)  
p <- matrix(NA, nrow = num_pheno, ncol = num_snps) 

for (i in 1:num_pheno) {
  for (j in 1:num_snps){
    lm_results <- summary(lm(y[[i]] ~ geno_data[, j]))
    p[i, j] <- lm_results$coefficients[-1, "Pr(>|t|)"]
    z_statistics[i, (j * 2 - 1)] <- lm_results$coefficients[-1, "Estimate"]
    z_statistics[i, (j * 2 )] <- lm_results$coefficients[-1, "Std. Error"]
  }
}

Z_Wald <- matrix(NA, num_pheno, num_snps)
for (i in 1:num_pheno) {
  for (j in 1:num_snps){
    Z_Wald[i, j] <-  z_statistics[i, (j * 2 - 1)] / z_statistics[i, (j * 2 )]
  }
}



# Next use the function defined in advance to get the needed result
# setwd("C:\\Users\\lying\\Desktop\\article") #path of R document  
# source("PDCAT.R") # functions of our method 


# lambdaSet is a set including all candidate lambda
# B represents times of resampling
lambdaSet <- c(1, 2, 3, Inf)
B <- 1000
# obtain our test statistics by the relationship between Wald test statistics and true coefficients 
pst_results <- lapply(1:nrow(Z_Wald), function(i) pst(geno_data, maf, Z_Wald[i, ]))

# normal distribution form of true coefficient
pst_dis_results <- pst_dis(geno_data, maf)

#divide multiple true coefficients into two statistics by directions of all coefficients
dst_results <- lapply(1:nrow(Z_Wald), function(i) dst(pst_results[[i]], pst_dis_results$s_RZ, lambdaSet))

# calculate parameter of the empirical distribution 

# paraments of the generalized chi-square distribution
paraments <- chi_parameter(pst_dis_results$s_RZ, pst_dis_results$cor_RZ, lambdaSet, B)

# to prevent the appearance of all PTstat < 0 , we repeat the generating paraments procedure until one of PTstat is not smaller than 0
PTstat_results <- matrix(NA, nrow = nrow(Z_Wald), ncol = (4*length(lambdaSet)))
colnames(PTstat_results) <- c("PTstat.neg.1", "PTstat.neg.2", "PTstat.neg.3", "PTstat.neg.Inf", "PTstat.pos.1", "PTstat.pos.2", "PTstat.pos.3", "PTstat.pos.Inf", "PTstat_pval.neg.1", "PTstat_pval.neg.2", "PTstat_pval.neg.3", "PTstat_pval.neg.Inf", "PTstat_pval.pos.1", "PTstat_pval.pos.2", "PTstat_pval.pos.3", "PTstat_pval.pos.Inf")
for (i in 1: nrow(Z_Wald)){
  PTstat_results[i, ] <- PTstat(pst_results[[i]], pst_dis_results$s_RZ, pst_dis_results$cor_RZ, paraments[grep("^a", names(paraments))], paraments[grep("^b", names(paraments))], paraments[grep("^d", names(paraments))], lambdaSet, B)
}
#Combine p-values by the Caubdy Combination test
pval <- rep(NA)
for(i in 1 : (nrow(PTstat_results))){
  pval[i] <- PDCAT(PTstat_results[i, -c(grep("^PTstat_pval", colnames(PTstat_results)))], PTstat_results[i, grep("^PTstat_pval", colnames(PTstat_results))])
}

# If all effect sizes (non_zero_snps) are equal to 0, this below results refer to the type I error rate.
# If any effect sizes are not equal to 0, this below results refer to the empirical power.
sum(pval <= 0.05)/ num_pheno # at the significance level 0.05
sum(pval <= 0.01)/ num_pheno # at the significance level 0.01
# When the number of phenotype datasets (num_pheno) is 1, it refers to the p-value by analyzing the association of multiple SNPs in a gene on a certain phenotype 

#Compare with other methods such as COMBAT, GATES, VEGAS, SimpleM
pval_c <- matrix(NA, num_pheno, 8)
colnames(pval_c) <- c("COMBAT", "GATES", "VEGAS.max", "VEGAS.p0.1", "VEGAS.p0.2", "VEGAS.p0.5", "VEGAS.all", "SimpleM") 
for(i in 1: num_pheno){
  pval_c[i, ] <- COMBAT(p[i, ], geno_data, vegas.pct <- c(0.1,0.2,0.5,1), ncores = 14)
}

# This below has the same meaning as our method
sum(pval_c <= 0.05)/ num_pheno # at the significance level 0.05
sum(pval_c <= 0.01)/ num_pheno # at the significance level 0.01




### Real data analysis

### Input dictionary of relationship between gene and SNP
### Input reference genotype data
### Input summary data of InCHIANTI
### setwd("C:\\Users\\lying\\Desktop\\article") #path of R document 
load('gene_dict_INCHIANTI.Rdata') # dictionary
load('genotype_part1.Rdata')   # reference genotype part 1 data
load('genotype_part2.Rdata')   # reference genotype part 2 data
load('INCHIANTI.Rdata') # summary data
genotype <- rbind(genotype_part1, genotype_part2)

# Calculate p-value of each gene in six phenotypes(AA, ALA, DHA, EDA, EPA, LA)
#AA
results_AA <- data.frame(matrix(NA, nrow = ncol(gene_dict_INCHIANTI), ncol = 6))
colnames(results_AA) <- c("gene", "COMBAT", "GATES", "VEGAS","SimpleM", "PDCAT")

for (i in 1:ncol(gene_dict_INCHIANTI)){
  snp <- c(na.omit(gene_dict_INCHIANTI[,i]))
  geno_data <- t(genotype[c(snp),])
  geno_data <- geno_data[, !apply(geno_data, 2, function(x) all(x == x[1]))] # If two columns of the reference data are identical, delete one of the columns
  maf <- INCHIANTI[c(colnames(geno_data)), "FREQ1"]
  Wald_AA <- INCHIANTI[c(colnames(geno_data)), "Wald_AA"]
  P_AA <- INCHIANTI[c(colnames(geno_data)), "P_AA"]
  pvalue_c <- COMBAT(x=P_AA, snp.ref=geno_data, vegas.pct = c(0.1, 0.2, 0.5, 1))
  pst_AA <- pst(geno_data, maf, Wald_AA)
  pst_dis_AA <- pst_dis(geno_data, maf)
  dst_AA <- dst(pst_AA, pst_dis_AA$s_RZ, lambdaSet)
  paraments <- chi_parameter(pst_dis_AA$s_RZ, pst_dis_AA$cor_RZ, lambdaSet, B)
  PTstat_AA <- PTstat(pst_AA, pst_dis_AA$s_RZ, pst_dis_AA$cor_RZ, paraments[grep("^a", names(paraments))], paraments[grep("^b", names(paraments))], paraments[grep("^d", names(paraments))], lambdaSet, B)
  pvalue <- PDCAT(PTstat_AA[-c(grep("^PTstat_pval", names(PTstat_AA)))], PTstat_AA[grep("^PTstat_pval", names(PTstat_AA))])
  pnum <- length(pvalue_c)
  results_AA[i,  "gene"] <- colnames(gene_dict_INCHIANTI)[i]
  results_AA[i,  2:5] <- c(pvalue_c[1:2], min(pvalue_c[3:(pnum-1)]), pvalue_c[pnum]) #choose the smallest p-value in the VEGAS-set
  results_AA[i,  6] <- pvalue
}

# ALA
results_ALA <- data.frame(matrix(NA,  nrow = ncol(gene_dict_INCHIANTI),  ncol = 6))
colnames(results_ALA) <- c("gene",  "COMBAT",  "GATES",  "VEGAS", "SimpleM",  "PDCAT")

for (i in 1:ncol(gene_dict_INCHIANTI)){
  snp <- c(na.omit(gene_dict_INCHIANTI[,i]))
  geno_data <- t(genotype[c(snp),])
  geno_data <- geno_data[, !apply(geno_data, 2, function(x) all(x == x[1]))] # If two columns of the reference data are identical, delete one of the columns
  maf <- INCHIANTI[c(colnames(geno_data)), "FREQ1"]
  Wald_ALA <- INCHIANTI[c(colnames(geno_data)), "Wald_ALA"]
  P_ALA <- INCHIANTI[c(colnames(geno_data)), "P_ALA"]
  pvalue_c <- COMBAT(x=P_ALA, snp.ref=geno_data, vegas.pct = c(0.1, 0.2, 0.5, 1))
  pst_ALA <- pst(geno_data, maf, Wald_ALA)
  pst_dis_ALA <- pst_dis(geno_data, maf)
  dst_ALA <- dst(pst_ALA, pst_dis_ALA$s_RZ, lambdaSet)
  paraments <- chi_parameter(pst_dis_ALA$s_RZ, pst_dis_ALA$cor_RZ, lambdaSet, B)
  PTstat_ALA <- PTstat(pst_ALA, pst_dis_ALA$s_RZ, pst_dis_ALA$cor_RZ, paraments[grep("^a", names(paraments))], paraments[grep("^b", names(paraments))], paraments[grep("^d", names(paraments))], lambdaSet, B)
  pvalue <- PDCAT(PTstat_ALA[-c(grep("^PTstat_pval", names(PTstat_ALA)))], PTstat_ALA[grep("^PTstat_pval", names(PTstat_ALA))])
  pnum <- length(pvalue_c)
  results_ALA[i,  "gene"] <- colnames(gene_dict_INCHIANTI)[i]
  results_ALA[i,  2:5] <- c(pvalue_c[1:2], min(pvalue_c[3:(pnum-1)]), pvalue_c[pnum]) #choose the smallest p-value in the VEGAS-set
  results_ALA[i,  6] <- pvalue
}

# DHA
results_DHA <- data.frame(matrix(NA,  nrow = ncol(gene_dict_INCHIANTI),  ncol = 6))
colnames(results_DHA) <- c("gene",  "COMBAT",  "GATES",  "VEGAS", "SimpleM",  "PDCAT")

for (i in 1:ncol(gene_dict_INCHIANTI)){
  snp <- c(na.omit(gene_dict_INCHIANTI[,i]))
  geno_data <- t(genotype[c(snp),])
  geno_data <- geno_data[, !apply(geno_data, 2, function(x) all(x == x[1]))] # If two columns of the reference data are identical, delete one of the columns
  maf <- INCHIANTI[c(colnames(geno_data)), "FREQ1"]
  Wald_DHA <- INCHIANTI[c(colnames(geno_data)), "Wald_DHA"]
  P_DHA <- INCHIANTI[c(colnames(geno_data)), "P_DHA"]
  pvalue_c <- COMBAT(x=P_DHA, snp.ref=geno_data, vegas.pct = c(0.1, 0.2, 0.5, 1))
  pst_DHA <- pst(geno_data, maf, Wald_DHA)
  pst_dis_DHA <- pst_dis(geno_data, maf)
  dst_DHA <- dst(pst_DHA, pst_dis_DHA$s_RZ, lambdaSet)
  paraments <- chi_parameter(pst_dis_DHA$s_RZ, pst_dis_DHA$cor_RZ, lambdaSet, B)
  PTstat_DHA <- PTstat(pst_DHA, pst_dis_DHA$s_RZ, pst_dis_DHA$cor_RZ, paraments[grep("^a", names(paraments))], paraments[grep("^b", names(paraments))], paraments[grep("^d", names(paraments))], lambdaSet, B)
  pvalue <- PDCAT(PTstat_DHA[-c(grep("^PTstat_pval", names(PTstat_DHA)))], PTstat_DHA[grep("^PTstat_pval", names(PTstat_DHA))])
  pnum <- length(pvalue_c)
  results_DHA[i,  "gene"] <- colnames(gene_dict_INCHIANTI)[i]
  results_DHA[i,  2:5] <- c(pvalue_c[1:2], min(pvalue_c[3:(pnum-1)]), pvalue_c[pnum]) #choose the smallest p-value in the VEGAS-set
  results_DHA[i,  6] <- pvalue
}

# EDA
results_EDA <- data.frame(matrix(NA,  nrow = ncol(gene_dict_INCHIANTI),  ncol = 6))
colnames(results_EDA) <- c("gene",  "COMBAT",  "GATES",  "VEGAS", "SimpleM",  "PDCAT")

for (i in 1:ncol(gene_dict_INCHIANTI)){
  snp <- c(na.omit(gene_dict_INCHIANTI[,i]))
  geno_data <- t(genotype[c(snp),])
  geno_data <- geno_data[, !apply(geno_data, 2, function(x) all(x == x[1]))] # If two columns of the reference data are identical, delete one of the columns
  maf <- INCHIANTI[c(colnames(geno_data)), "FREQ1"]
  Wald_EDA <- INCHIANTI[c(colnames(geno_data)), "Wald_EDA"]
  P_EDA <- INCHIANTI[c(colnames(geno_data)), "P_EDA"]
  pvalue_c <- COMBAT(x=P_EDA, snp.ref=geno_data, vegas.pct = c(0.1, 0.2, 0.5, 1))
  pst_EDA <- pst(geno_data, maf, Wald_EDA)
  pst_dis_EDA <- pst_dis(geno_data, maf)
  dst_EDA <- dst(pst_EDA, pst_dis_EDA$s_RZ, lambdaSet)
  paraments <- chi_parameter(pst_dis_EDA$s_RZ, pst_dis_EDA$cor_RZ, lambdaSet, B)
  PTstat_EDA <- PTstat(pst_EDA, pst_dis_EDA$s_RZ, pst_dis_EDA$cor_RZ, paraments[grep("^a", names(paraments))], paraments[grep("^b", names(paraments))], paraments[grep("^d", names(paraments))], lambdaSet, B)
  pvalue <- PDCAT(PTstat_EDA[-c(grep("^PTstat_pval", names(PTstat_EDA)))], PTstat_EDA[grep("^PTstat_pval", names(PTstat_EDA))])
  pnum <- length(pvalue_c)
  results_EDA[i,  "gene"] <- colnames(gene_dict_INCHIANTI)[i]
  results_EDA[i,  2:5] <- c(pvalue_c[1:2], min(pvalue_c[3:(pnum-1)]), pvalue_c[pnum]) #choose the smallest p-value in the VEGAS-set
  results_EDA[i,  6] <- pvalue
}

# EPA
results_EPA <- data.frame(matrix(NA,  nrow = ncol(gene_dict_INCHIANTI),  ncol = 6))
colnames(results_EPA) <- c("gene",  "COMBAT",  "GATES",  "VEGAS", "SimpleM",  "PDCAT")

for (i in 1:ncol(gene_dict_INCHIANTI)){
  snp <- c(na.omit(gene_dict_INCHIANTI[,i]))
  geno_data <- t(genotype[c(snp),])
  geno_data <- geno_data[, !apply(geno_data, 2, function(x) all(x == x[1]))] # If two columns of the reference data are identical, delete one of the columns
  maf <- INCHIANTI[c(colnames(geno_data)), "FREQ1"]
  Wald_EPA <- INCHIANTI[c(colnames(geno_data)), "Wald_EPA"]
  P_EPA <- INCHIANTI[c(colnames(geno_data)), "P_EPA"]
  pvalue_c <- COMBAT(x=P_EPA, snp.ref=geno_data, vegas.pct = c(0.1, 0.2, 0.5, 1))
  pst_EPA <- pst(geno_data, maf, Wald_EPA)
  pst_dis_EPA <- pst_dis(geno_data, maf)
  dst_EPA <- dst(pst_EPA, pst_dis_EPA$s_RZ, lambdaSet)
  paraments <- chi_parameter(pst_dis_EPA$s_RZ, pst_dis_EPA$cor_RZ, lambdaSet, B)
  PTstat_EPA <- PTstat(pst_EPA, pst_dis_EPA$s_RZ, pst_dis_EPA$cor_RZ, paraments[grep("^a", names(paraments))], paraments[grep("^b", names(paraments))], paraments[grep("^d", names(paraments))], lambdaSet, B)
  pvalue <- PDCAT(PTstat_EPA[-c(grep("^PTstat_pval", names(PTstat_EPA)))], PTstat_EPA[grep("^PTstat_pval", names(PTstat_EPA))])
  pnum <- length(pvalue_c)
  results_EPA[i,  "gene"] <- colnames(gene_dict_INCHIANTI)[i]
  results_EPA[i,  2:5] <- c(pvalue_c[1:2], min(pvalue_c[3:(pnum-1)]), pvalue_c[pnum]) #choose the smallest p-value in the VEGAS-set
  results_EPA[i,  6] <- pvalue
}

# LA
results_LA <- data.frame(matrix(NA,  nrow = ncol(gene_dict_INCHIANTI),  ncol = 6))
colnames(results_LA) <- c("gene",  "COMBAT",  "GATES",  "VEGAS", "SimpleM",  "PDCAT")

for (i in 1:ncol(gene_dict_INCHIANTI)){
  snp <- c(na.omit(gene_dict_INCHIANTI[,i]))
  geno_data <- t(genotype[c(snp),])
  geno_data <- geno_data[, !apply(geno_data, 2, function(x) all(x == x[1]))] # If two columns of the reference data are identical, delete one of the columns
  maf <- INCHIANTI[c(colnames(geno_data)), "FREQ1"]
  Wald_LA <- INCHIANTI[c(colnames(geno_data)), "Wald_LA"]
  P_LA <- INCHIANTI[c(colnames(geno_data)), "P_LA"]
  pvalue_c <- COMBAT(x=P_LA, snp.ref=geno_data, vegas.pct = c(0.1, 0.2, 0.5, 1))
  pst_LA <- pst(geno_data, maf, Wald_LA)
  pst_dis_LA <- pst_dis(geno_data, maf)
  dst_LA <- dst(pst_LA, pst_dis_LA$s_RZ, lambdaSet)
  paraments <- chi_parameter(pst_dis_LA$s_RZ, pst_dis_LA$cor_RZ, lambdaSet, B)
  PTstat_LA <- PTstat(pst_LA, pst_dis_LA$s_RZ, pst_dis_LA$cor_RZ, paraments[grep("^a", names(paraments))], paraments[grep("^b", names(paraments))], paraments[grep("^d", names(paraments))], lambdaSet, B)
  pvalue <- PDCAT(PTstat_LA[-c(grep("^PTstat_pval", names(PTstat_LA)))], PTstat_LA[grep("^PTstat_pval", names(PTstat_LA))])
  pnum <- length(pvalue_c)
  results_LA[i,  "gene"] <- colnames(gene_dict_INCHIANTI)[i]
  results_LA[i,  2:5] <- c(pvalue_c[1:2], min(pvalue_c[3:(pnum-1)]), pvalue_c[pnum]) #choose the smallest p-value in the VEGAS-set
  results_LA[i,  6] <- pvalue
}
