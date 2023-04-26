rm(list = ls())
library(ggplot2)
library("tidyverse")
library("ks")
library("reldist")
library("boot")
library("MetricsWeighted")
set.seed(9)
setwd("~/Documents/BUaca/research/causal/writing/concep_drafts/22/new_version/rebutal")


######

###
###function to generate randomized trials
###
#total.n: total number of patients
#prob.l: P(L=1)
#prob.d.l0A0: probability of death in patients with L=0 and A=0, i.e., P(D=1|L=0, A=0)
#prob.d.l1A0: probability of death in patients with L=1 and A=0, i.e., P(D=1|L=1, A=0)
#prob.d.l0A1: probability of death in patients with L=0 and A=1, i.e., P(D=1|L=0, A=1)
#prob.d.l1A1: probability of death in patients with L=1 and A=1, i.e., P(D=1|L=1, A=1)
#beta0, beta1, beta2: to generate Y, Y = beta0 + beta1 * A + beta2 * L 
#pats.sd: standard deviation in Y
#d.value: value assigned for those who died.


genrct <- function(total.n, prob.l,
                   prob.d.l0A0, prob.d.l1A0, prob.d.l0A1, prob.d.l1A1,
                   beta0, beta1, beta2,
                   pats.sd, d.value) {
  l <- rbinom(total.n, 1, prob.l)  #generate L=1 according to P(L=1)
  l0.num <- sum(l == 0)      #total length of L = 0
  l1.num <- sum(l == 1)      #total length of L = 1
  
  A <- numeric(total.n)    #a vector of treatment of 0's 
  a1id <- sample(1:total.n, size = (total.n / 2))   ##randomize treatment A, generate ids of A=1
  A[a1id] <- 1    ##  assign A=1                          
    
  rct.pats <- data.frame(l, A)
  l0.A0.num <- length(which(rct.pats$l == 0 & rct.pats$A == 0))   #length of L = 0 and A = 0
  l1.A0.num <- length(which(rct.pats$l == 1 & rct.pats$A == 0))   #length of L = 1 and A = 0
  l0.A1.num <- length(which(rct.pats$l == 0 & rct.pats$A == 1))   #length of L = 0 and A = 1
  l1.A1.num <- length(which(rct.pats$l == 1 & rct.pats$A == 1))   #length of L = 1 and A = 1
  
  D <- numeric(total.n)  #a total vector of death
  
  #generate death vector according to the probability of death in each stratum
  D[which(rct.pats$l == 0 & rct.pats$A == 0)] <- rbinom(l0.A0.num, 1, prob.d.l0A0)     #generate D = 1 in the stratum of L = 0 and A = 0
  D[which(rct.pats$l == 0 & rct.pats$A == 1)] <- rbinom(l0.A1.num, 1, prob.d.l0A1)     #generate D = 1 in the stratum of L = 0 and A = 1
  D[which(rct.pats$l == 1 & rct.pats$A == 0)] <- rbinom(l1.A0.num, 1, prob.d.l1A0)     #generate D = 1 in the stratum of L = 1 and A = 0
  D[which(rct.pats$l == 1 & rct.pats$A == 1)] <- rbinom(l1.A1.num, 1, prob.d.l1A1)     #generate D = 1 in the stratum of L = 1 and A = 1
  
  #generate a whole dataframe of patients
  rct.pats <- data.frame(l, A, D)
  
  #### assign values for Y
  meanY <- beta0 + beta1 * rct.pats$A + beta2 * rct.pats$l  #calculate mean of Y
  Y <- rnorm(total.n, meanY, pats.sd)   #generate all Y values
  rct.pats <- data.frame(l, A, D, Y)      
  rct.pats[which(rct.pats$D == 1), "Y"] <- d.value    ##for those who died, assigned a death value
  
  return(rct.pats)
}




###
### generate randomized trials data for N = 500, 1500, 5000


#generate randomized data for N = 500
pats.rct.500 <- list()
for (d in 1:2000) {
  temp.rct <- genrct(
    total.n = 500, prob.l = 0.6,
    prob.d.l0A0 = 0.20, prob.d.l0A1 = 0.20, prob.d.l1A0 = 0.35, prob.d.l1A1 = 0.15,
    beta0 = 3, beta1 = 0.3, beta2 = -3,
    pats.sd = 1, d.value = -1000
  )
  pats.rct.500[[d]] <- temp.rct
}
save(pats.rct.500, file = "pats.rct.500.rda")


#generate randomized data for N = 1500
pats.rct.1500 <- list()
for (d in 1:2000) {
  temp.rct <- genrct(
    total.n = 1500, prob.l = 0.6,
    prob.d.l0A0 = 0.20, prob.d.l0A1 = 0.20, prob.d.l1A0 = 0.35, prob.d.l1A1 = 0.15,
    beta0 = 3, beta1 = 0.3, beta2 = -3,
    pats.sd = 1, d.value = -1000
  )
  pats.rct.1500[[d]] <- temp.rct
}
save(pats.rct.1500, file = "pats.rct.1500.rda")


#generate randomized data for N = 5000
pats.rct.5000 <- list()
for (d in 1:2000) {
  temp.rct <- genrct(
    total.n = 5000, prob.l = 0.6,
    prob.d.l0A0 = 0.20, prob.d.l0A1 = 0.20, prob.d.l1A0 = 0.35, prob.d.l1A1 = 0.15,
    beta0 = 3, beta1 = 0.3, beta2 = -3,
    pats.sd = 1, d.value = -1000
  )
  pats.rct.5000[[d]] <- temp.rct
}
save(pats.rct.5000, file = "pats.rct.5000.rda")

load("pats.rct.500.rda")
load("pats.rct.1500.rda")
load("pats.rct.5000.rda")

#################################################
################### Medians for RCT
#################################################
##get medians from the simualtion data
##

summary_results <- function(pats_dat_input){
  num_datsets <- length(pats_dat_input)
  
  a0medsurvivors <- numeric(num_datsets)    #to store results for median in survivors for those A = 0
  a1medsurvivors <- numeric(num_datsets)    #to store results for median in survivors for those A = 1
  a0surincmed <- numeric(num_datsets)       #to store results for the survival-incorporated median for those A = 0
  a1surincmed <- numeric(num_datsets)   
  
  for (d in 1:num_datsets){
    temp_pats <- pats_dat_input[[d]]
    
    a0medsurvivors[d] <- median(temp_pats[which(temp_pats$A == 0 & temp_pats$D == 0), ]$Y)  #results for median in survivors for those A = 0
    a1medsurvivors[d] <- median(temp_pats[which(temp_pats$A == 1 & temp_pats$D == 0), ]$Y)
    a0surincmed[d] <- median(temp_pats[which(temp_pats$A == 0), ]$Y)               #results for the survival-incorporated medians for those A = 0
    a1surincmed[d] <- median(temp_pats[which(temp_pats$A == 1), ]$Y)
  }
  
  return(data.frame(a0medsurvivors, a1medsurvivors, a0surincmed, a1surincmed))
}

pats500_results  <- summary_results(pats.rct.500)
pats1500_results <- summary_results(pats.rct.1500)
pats5000_results <- summary_results(pats.rct.5000)



#############
######summary funtion

#function to summarize MSE
#rct: input rct summarized results
#rct.truth: truth of the trials

summary_rMSE <- function(dat_input, rct.truth) {
  a0_medsurvivors <- sqrt(mean((dat_input[, 1] - rct.truth[1])^2))
  a1_medsurvivors <- sqrt(mean((dat_input[, 2] - rct.truth[2])^2))
  
  a0_surincmed <- sqrt(mean((dat_input[, 3] - rct.truth[3])^2))
  a1_surincmed <- sqrt(mean((dat_input[, 4] - rct.truth[4])^2))
  
  dif_medsurvivors <- sqrt(mean((dat_input[, 2] - dat_input[, 1] - rct.truth[5])^2))
  dif_surincmed <- sqrt(mean((dat_input[, 4] - dat_input[, 3] - rct.truth[6])^2))

  
  rMSE <- c(a0_medsurvivors, a1_medsurvivors, a0_surincmed, a1_surincmed, dif_medsurvivors, dif_surincmed)
  names(rMSE) <- c("a0_medsurvivors", "a1_medsurvivors", "a0_surincmed", "a1_surincmed", "dif_medsurvivors", "dif_surincmed")
  return(rMSE)
  #rMSE[1]: rMSE of median in the survivors in A = 0
  #rMSE[2]: rMSE of median in the survivors in A = 1
  #rMSE[3]: rMSE of the survival-incorporated median in A = 0
  #rMSE[4]: rMSE of the survival-incorporated median in A = 1
  #rMSE[5]: rMSE of the difference between the median in the survivors
  #rMSE[6]: rMSE of the difference between the survival-incorporated median
}

#function to summarize bias
#rct: input rct summarized results
#rct.truth: truth of the trials
summary_bias <- function(dat_input, rct.truth) {
  a0_medsurvivors <- mean(dat_input[, 1]) - rct.truth[1]
  a1_medsurvivors <- mean(dat_input[, 2]) - rct.truth[2]
  
  a0_surincmed <- mean(dat_input[, 3]) - rct.truth[3]
  a1_surincmed <- mean(dat_input[, 4]) - rct.truth[4]
  
  dif_medsurvivors <- mean(dat_input[, 2] - dat_input[, 1]) - rct.truth[5]
  dif_surincmed <- mean(dat_input[, 4] - dat_input[, 3]) - rct.truth[6]
  
  
  biases <- c(a0_medsurvivors, a1_medsurvivors, a0_surincmed, a1_surincmed, dif_medsurvivors, dif_surincmed)
  names(biases) <- c("a0_medsurvivors", "a1_medsurvivors", "a0_surincmed", "a1_surincmed", "dif_medsurvivors", "dif_surincmed")
  return(biases)
  #biases[1]: bias of median in the survivors in A = 0
  #biases[2]: bias of median in the survivors in A = 1
  #biases[3]: bias of the survival-incorporated median in A = 0
  #biases[4]: bias of the survival-incorporated median in A = 1
  #biases[5]: bias of the difference between the median in the survivors
  #biases[6]: bias of the difference between the survival-incorporated median
}

#truth
rct.truth <- c(1.18397, 1.155069, 0.09280129, 0.6701869, -0.028901, 0.5773856)


#summarize results for patients with N = 500
colMeans(pats500_results)
round(summary_rMSE(pats500_results, rct.truth), 4)
round(summary_bias(pats500_results, rct.truth), 4)

#summarize results for patients with N = 1500
colMeans(pats1500_results) 
round(summary_rMSE(pats1500_results, rct.truth), 4)
round(summary_bias(pats1500_results, rct.truth), 4)

#summarize results for patients with N = 5000
colMeans(pats5000_results)
round(summary_rMSE(pats5000_results, rct.truth), 4)
round(summary_bias(pats5000_results, rct.truth), 4)
