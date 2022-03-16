rm(list= ls())
library("boot")
set.seed(10)
##read the data
##data are provided by Ding
setwd('/Users/Qingyan/Documents/BUaca/research/causal/application/Ding2011script')
data1=read.table("data1.txt",sep="\t",header=TRUE)  #QOL scores
data2=read.table("data2.txt",sep="\t",header=TRUE)  #other information
GQOL=read.table("GQOL.txt",sep="\t",header=TRUE)    #change of QoL scores ('dscore') and treatment information


##########################################
##data preprocessing: adapted from Ding's codes
j=1
id<-NULL
for(i in 1:length(GQOL[,1]))
{
  repeat
  {
    if(data1$PATNO[j]==GQOL$PATNO[i])
    {
      id=rbind(id,j)
      break
    }
    j=j+1
  }      
}

##combine the data
data615=cbind(GQOL,data1[id,-1])
##delete some columns of the 'dscores' in the data
data615=data615[,-c(6,7,8)]


##Survival status indicator 
S=ifelse(is.na(data615$score12),0,1)
data615=cbind(data615,S)


##score difference, dscore =score12-score0
dscore12=data615$score12-data615$score0
data615=cbind(data615,dscore12)


#################################Oct 30#####################################
##Latent Variable method
##A: the baseline score: continuous

hist(data615$score0)
##discretize A into different categories. For example,  A=1 if score <= 50; A=2 if 50 < score < 75; A = 3 if score >= 75 
N = length(data615[,1])
A=1:N

#we use the 2nd type of decretization
if(T){                 #first type of discretization
for(i in 1:N)
{
  if(data615$score0[i]<=50)    A[i]=1
  else if(data615$score0[i]>=75)   A[i]=3
  else A[i]=2
}
}
if(F){               #second type of discretization; 4 categories, <=25, (25, 50], (50, 75], >75
  for(i in 1:N)
  {
    if(data615$score0[i]<=25)    A[i]=1
    else if(data615$score0[i]>25 & data615$score0[i]<=50)   A[i]=2
    else if(data615$score0[i]>75) A[i]=4 
    else A[i]=3
  } 
}


data615 = cbind(data615,A)
table(data615$A)

input_dat = data615
discret.n = 4
#########################################################
#####The function to estimate the SACE with monotonicity

mu_1ll <- function(input, i, discret.n, only_sace, used_for_boot){
  #input: swog data
  #i: indices for bootstrap function
  #discret.n: level of discretization
  #used_for_boot: if this function is used for bootstrap
  
  #if this function is used for bootstrap, select indices from input data with replacement
  if(used_for_boot == T){     
    input_dat <- input[i, ]
  }
  else{
    input_dat = input
  }
  
  ############
  #######step 1
  pi_ll = sum((input_dat$Z == 0) & (input_dat$S == 1))/sum(input_dat$Z == 0)  #probability of always-survivor
  pi_dd = sum((input_dat$Z == 1) & (input_dat$S == 0))/sum(input_dat$Z == 1)  #probability of never-survivor
  pi_ld = 1 - pi_ll - pi_dd                                                   #probability of protected
  
  #calculate mu_{0,ll}
  mu_0ll = sum(input_dat[which(input_dat$Z == 0 & input_dat$S == 1), "dscore12"])/length(which(input_dat$Z == 0 & input_dat$S == 1))
  
  ############
  #######step 2
  
  #calculate p(a|LL) for each a
  p_a_ll <- numeric(discret.n)    #first initialize an empty vector for p(a|LL)
  for(a in 1:discret.n){          #then calculate this probability
    p_a_ll[a] <- sum((input_dat$A == a & input_dat$Z == 0 & input_dat$S == 1)) / sum(input_dat$Z == 0 & input_dat$S == 1)
  }
  
  #calculate p(a|DD) for each a
  p_a_dd <- numeric(discret.n)
  for(a in 1:discret.n){
    p_a_dd[a] <- sum((input_dat$A == a & input_dat$Z == 1 & input_dat$S == 0)) / sum(input_dat$Z == 1 & input_dat$S == 0)
  }
  
  #calculate marginal p(a) for each a
  p_a <- numeric(discret.n)
  for(a in 1:discret.n){
    p_a[a] <- sum(input_dat$A == a) / nrow(input_dat)
  }
  #calculate p(a|LD) for each a
  p_a_ld <- numeric(discret.n)
  for(a in 1:discret.n){
    p_a_ld[a] <- (p_a[a] - pi_ll * p_a_ll[a] - pi_dd * p_a_dd[a])/pi_ld
  }
  
  ############
  ######Step 3
  
  #calculate beta(a) = P(G=LL|Z=1, S=1, A=a) for each a 
  beta_a <- numeric(discret.n) #initialize beta_a vector with all 0
  for(a in 1:discret.n){
    beta_a[a] <- pi_ll * p_a_ll[a] / (pi_ll * p_a_ll[a] + pi_ld * p_a_ld[a])
  }
 
  ############
  ######Step 4
  
  #calculate beta(a) = P(G=LL|Z=1, S=1, A=a) for each a 
  mu_a <- numeric(discret.n) 
  for(a in 1:discret.n){
    denominator = sum(input_dat$A == a & input_dat$Z == 1 & input_dat$S == 1)
    numerator = sum(input_dat[which(input_dat$A == a & input_dat$Z == 1 & input_dat$S == 1), "dscore12"])
    
    mu_a[a] <- numerator/denominator
  }
  
  #############
  ######step 5
  
  beta_a_c <- 1 - beta_a   #create 1 minus beta(a) for each a
  rg_dat <- data.frame(cbind(mu_a, beta_a, beta_a_c))   #create datasets with mu_a, beta_a, and (1 - beta_a)
  res <- lm(mu_a ~ beta_a + beta_a_c - 1, data = rg_dat) #regression analysis (without intercept)

  
  if(only_sace == T)
    return(res$coefficients[1] - mu_0ll)                                        #return the sace = mu_1_ll - mu_0_ll
  if(only_sace == F)
    return(list(mu_1ll = res$coefficients[1], mu_1ld = res$coefficients[2], mu_0ll = mu_0ll, sace_mon = res$coefficients[1] - mu_0ll))    #return both the effect seperately
}





############
######return the sace and bootstrapped CI
boot.res <- boot(data615,        #data
                 mu_1ll,         #function statistics
                 discret.n=3,    #level of discretization: 3
                 only_sace= T,     #parameter of the function: we estimate the difference of sace = mu_1_ll - mu_0_ll
                 used_for_boot = T,  #parameter of the function:used for bootstrap
                 R=2000)         #2000 bootstrap replicates


print(boot.res)                  #print bootstrap estimate
boot.ci(boot.res, type= "perc")  #print bootstrap CI



############
#####To return the estimates under each treatmetns
mu_1ll(data615, 0, discret.n = 3, only_sace = F, used_for_boot = F)
