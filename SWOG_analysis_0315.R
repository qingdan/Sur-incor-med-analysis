rm(list = ls())
setwd("/Users/Qingyan/Documents/BUaca/research/causal/application/dataset/SWOG")
set.seed(10)

library("survival")
library("survminer")
library("tidyverse")
library("MetricsWeighted")
library("quantreg")
library("grDevices")

tox <- read.csv("tox.csv")
ovrll_qol2 <- read.csv("ovrll_qol2.csv")
psa <- read.csv("psa.csv")
patient <- read.csv("patient.csv")
qlcover2 <- read.csv("qlcover2.csv")

swog_Ding <- read.table("/Users/Qingyan/Documents/BUaca/research/causal/application/dataset/Ding/swogdata.txt")
###########################################################
####################### patients demographics#############################
###########################################################
# total nums
unique(patient$ALTPATID)

# total death
length(which(patient$SURIND == 1)) # SURIND = 1 represent death

# recode trt: 0 = MP; 1 = DE
patient$trt <- abs(as.numeric(as.factor(patient$ARMNAME)) - 2)
table(patient$trt)

# death in each trt
sum(patient$SURIND[which(patient$trt == 1)])
sum(patient$SURIND[which(patient$trt == 0)])


patient$ARMNAME <- as.factor(patient$ARMNAME)
patient$ARMNAME <- relevel(patient$ARMNAME, ref = "Mitoxantrone + Prednisone")

### some information: survival plot
par(mfrow = c(1, 1))

plot(survfit(Surv(SURTIM, SURIND) ~ ARMNAME, data = patient),
  xlab = "days",
  ylab = "survival-prob",
  main = "survival till day 395"
)

fit <- survfit(Surv(SURTIM, SURIND) ~ ARMNAME, data = patient)


pdf("/Users/Qingyan/Documents/BUaca/research/causal/writing/method_draft/figures/surplot.pdf", width = 6, height = 4.5)
ggsurvplot(fit, data = patient, xlab = "Time: days", ylim = c(0.3, 1), xlim = c(0, 380), legend.labs = c("Mitoxanrone + Prednoisone (MP)", "Docetaxel+Estramustine (DE)"), legend.title = "")
dev.off()


ggsurvplot(fit, data = patient, xlab = "Time: days", ylim = c(0.3, 1))
# ggsurvplot(fit, data = patient, xlab = "Time: days", ylim = c(0.3, 1))

#################################################################
####################### data preprocessing########################
#################################################################

# interest of study: (1) 12-month score (2) 12-month change from baseline
interest <- "12month"
# interest = "change"


# step 1
# remove 29 subjects that are missing from baseline
cat("missing at baseline (actually all missing): ", sum(is.na(ovrll_qol2$SCORE0)))
ovrll_qol2 <- ovrll_qol2[!is.na(ovrll_qol2$SCORE0), ]


# step 2
# decide cutoff range: 365 days
patient12 <- patient

patient12$SURIND <- ifelse(patient12$SURTIM > 365, 0, 1)
patient12$SURTIM <- ifelse(patient12$SURTIM > 365, 365, patient12$SURTIM)

die_id <- patient12 %>% filter(SURIND == 1) %>% select(ALTPATID)
sur_id <- patient12 %>% filter(SURIND == 0) %>% select(ALTPATID)

dies <- ovrll_qol2[ovrll_qol2$ALTPATID %in% die_id[[1]], ] # a subset of those who died
surs <- ovrll_qol2[ovrll_qol2$ALTPATID %in% sur_id[[1]], ] # a subset of survivors


# step 3
# impute missing reasons due to illness (assign -5 as their value)
# find the illness id, for month 12 and reason is illness (QNO == 1)
illness_id <- qlcover2 %>%
  filter(QLT9916 == 12 & QNO == 1) %>%
  select(ALTPATID)

# impute 37 survivors with their scores missing due to illness
if (interest == "12month") { # if the interest is the score at month 12
  surs[which(surs$ALTPATID %in% illness_id[[1]] & is.na(surs$SCORE12)), "SCORE12"] <- -5
} else {
  surs[which(surs$ALTPATID %in% illness_id[[1]] & is.na(surs$SCORE12)), "SCORE12"] <- -100 # if the interest is the change from baseline
}


# step 4
# if a subject with a missing score is alive at 365 but dies before 395, they probably has a bad status and we view them as illness
# impute those 12 subjects as illness
d_after365_id <- patient[which(patient$SURTIM > 365 & patient$SURIND == 1), "ALTPATID"]

if (interest == "12month") { # if the interest is the score at month 12
  surs[which(surs$ALTPATID %in% d_after365_id & is.na(surs$SCORE12)), "SCORE12"] <- -5
} else {
  surs[which(surs$ALTPATID %in% d_after365_id & is.na(surs$SCORE12)), "SCORE12"] <- -100 # if the interest is the change from baseline
}


# step 5
# if now survivors are missing at month12 but not missing at cycle 8, we will use their cycle8 score instead
no12but8_id <- is.na(surs$SCORE12) & !is.na(surs$SCORE8)
surs[no12but8_id, "SCORE12"] <- surs[no12but8_id, "SCORE8"]


# step 6
# find subjects with missing values due to insitituition error (MCAR):
# counts:52
isti_error_id <- qlcover2 %>%
  filter(QLT9916 == 12 & QNO == 5) %>%
  select(ALTPATID)
MCAR_patient_id <- surs[which(surs$ALTPATID %in% isti_error_id[[1]] & is.na(surs$SCORE12)), "ALTPATID"]


# step 7
# find subjects with missing values due to MAR
# counts: 50
MAR_patient_id <- surs[which(is.na(surs$SCORE12) & !surs$ALTPATID %in% MCAR_patient_id), "ALTPATID"]


##########################
# other data manipulation
##########################
if (interest == "12month") { # if the interest is the score at month 12
  dies$SCORE12 <- -10
} else {
  dies$SCORE12 <- -500 # if the interest is the change from baseline
}

# combined death and survive all together
clean_qol2 <- rbind(surs, dies)
clean_qol2 <- merge(
  x = clean_qol2, y = patient12[, c("ALTPATID", "ARMNAME", "SURIND", "AGE", "RACECAT", "PS")],
  by = "ALTPATID", all.x = TRUE
) # combine other covariate information

# convert string to numeric
clean_qol2$trt <- abs(as.numeric(as.factor(clean_qol2$ARMNAME)) - 2)

# MCAR indicator
clean_qol2$C_error <- numeric(nrow(clean_qol2))
clean_qol2$C_error[which(clean_qol2$ALTPATID %in% MCAR_patient_id)] <- 1

# MAR indicator
clean_qol2$C_MAR <- numeric(nrow(clean_qol2))
clean_qol2$C_MAR[which(clean_qol2$ALTPATID %in% MAR_patient_id)] <- 1

# calculate p.MCAR
qol_surs <- clean_qol2 %>% filter(clean_qol2$SCORE12 >= 0 | is.na(clean_qol2$SCORE12))
p.c.MCAR <- sum(qol_surs$C_error) / nrow(qol_surs) # global MCAR

p.c.MCAR.0 <- sum(qol_surs[qol_surs$trt == 0, "C_error"]) / nrow(qol_surs[qol_surs$trt == 0, ]) # p.MCAR under A = 0
p.c.MCAR.1 <- sum(qol_surs[qol_surs$trt == 1, "C_error"]) / nrow(qol_surs[qol_surs$trt == 1, ]) # p.MCAR under A = 1


#################################################################
############################analysis#############################
#################################################################
##########################
# some summary
##########################
# total subjects in each trt
sum(clean_qol2$trt == 0)
sum(clean_qol2$trt == 1)


# number of death; prob of death
sum(clean_qol2$trt == 0 & clean_qol2$SURIND == 1)
sum(clean_qol2$trt == 0 & clean_qol2$SURIND == 1) / sum(clean_qol2$trt == 0)
sum(clean_qol2$trt == 1 & clean_qol2$SURIND == 1)
sum(clean_qol2$trt == 1 & clean_qol2$SURIND == 1) / sum(clean_qol2$trt == 1)
hist(clean_qol2[which(clean_qol2$SCORE12 > 0), "SCORE12"], xlab = "QoL scores", main = "histogram of QoL at month 12")


#######################
# First: focus on trt = 0
#########################
### trt0 healthy survivors
qol_trt0_surs <- clean_qol2 %>% filter(clean_qol2$trt == 0 & (clean_qol2$SCORE12 >= 0 | is.na(clean_qol2$SCORE12)))

### trt0 healthy survivors excluding MCAR = 1 (for ps modeling)
trt0_ps_model_dat <- qol_trt0_surs[-which(qol_trt0_surs$C_error == 1), ]

# ps modelling
res <- glm(C_MAR ~ AGE + PS + RACECAT + SCORE0, family = binomial(), data = trt0_ps_model_dat)
trt0_ps_model_dat$p.c.MAR <- round(predict(res, type = "response"), 3)
MAR_0_patient_id <- trt0_ps_model_dat[trt0_ps_model_dat$C_MAR == 0, "ALTPATID"]

# adding weight
qol_trt0_surs$w <- 1
qol_trt0_surs[qol_trt0_surs$ALTPATID %in% MAR_0_patient_id, "w"] <- 1 / ((1 - p.c.MCAR.0) * (1 - trt0_ps_model_dat[trt0_ps_model_dat$C_MAR == 0, "p.c.MAR"]))

## combined it with death under trt = 0
# trt0 all data, includeing survivors and death
trt0_all <- clean_qol2 %>% filter(clean_qol2$trt == 0)
trt0_all$w <- 1 # for dead and ill, weight is 1
trt0_all[match(qol_trt0_surs$ALTPATID, trt0_all$ALTPATID), "w"] <- qol_trt0_surs$w

## sum of weights of those with a score (including imputed)
sum(trt0_all$w[!is.na(trt0_all$SCORE12)])


# result:33.33
if (interest == "12month") {
  weighted_quantile(trt0_all$SCORE12, w = trt0_all$w, probs = 0.5, na.rm = TRUE)
} else {
  trt0_all$delta <- trt0_all$SCORE12 - trt0_all$SCORE0
  weighted_quantile(trt0_all$delta, w = trt0_all$w, probs = 0.5, na.rm = TRUE)
}




#######################
# next: focus on trt = 1
#########################
### trt1 healthy survivors
qol_trt1_surs <- clean_qol2 %>% filter(clean_qol2$trt == 1 & (clean_qol2$SCORE12 >= 0 | is.na(clean_qol2$SCORE12)))

### trt1 healthy survivors excluding MCAR = 1 (for ps modeling)
trt1_ps_model_dat <- qol_trt1_surs[-which(qol_trt1_surs$C_error == 1), ]

# ps modelling
res <- glm(C_MAR ~ AGE + PS + RACECAT + SCORE0, family = binomial(), data = trt1_ps_model_dat)
trt1_ps_model_dat$p.c.MAR <- round(predict(res, type = "response"), 3)
MAR_0_patient_id <- trt1_ps_model_dat[trt1_ps_model_dat$C_MAR == 0, "ALTPATID"]

# adding weight
qol_trt1_surs$w <- 1
qol_trt1_surs[qol_trt1_surs$ALTPATID %in% MAR_0_patient_id, "w"] <- 1 / ((1 - p.c.MCAR.1) * (1 - trt1_ps_model_dat[trt1_ps_model_dat$C_MAR == 0, "p.c.MAR"]))

## combined it with death under trt = 1
# trt1 all data, includeing survivors and death
trt1_all <- clean_qol2 %>% filter(clean_qol2$trt == 1)
trt1_all$w <- 1
trt1_all[match(qol_trt1_surs$ALTPATID, trt1_all$ALTPATID), "w"] <- qol_trt1_surs$w

sum(trt1_all$w[!is.na(trt1_all$SCORE12)])


if (interest == "12month") {
  weighted_quantile(trt1_all$SCORE12, w = trt1_all$w, probs = 0.5, na.rm = TRUE)
} else {
  trt1_all$delta <- trt1_all$SCORE12 - trt1_all$SCORE0
  weighted_quantile(trt1_all$delta, w = trt1_all$w, probs = 0.5, na.rm = TRUE)
}







###############################################################
##### Quantile regression for confidence interval: 12-month score
##################################
########### for survival-incorproated median
### trt=0
all_subs_w <- rbind(trt0_all, trt1_all)
rq.res <- rq(SCORE12 ~ trt, tau = 0.5, data = all_subs_w, weights = w)
rq.sum <- summary(rq.res, se = "iid")

vc.coe <- rq.sum$coefficients
CI <- cbind(vc.coe[, 1] - 1.96 * vc.coe[, 2], vc.coe[, 1] + 1.96 * vc.coe[, 2])

### trt=1
all_subs_w$trt <- abs(all_subs_w$trt - 1)
rq.res <- rq(SCORE12 ~ trt, tau = 0.5, data = all_subs_w, weights = w)
rq.sum <- summary(rq.res, se = "iid")

vc.coe <- rq.sum$coefficients
CI <- cbind(vc.coe[, 1] - 1.96 * vc.coe[, 2], vc.coe[, 1] + 1.96 * vc.coe[, 2])


########################
### for median in the survivors
### trt=0
all_subs_w <- rbind(trt0_all, trt1_all)
all_surs_w <- all_subs_w[all_subs_w$SURIND == 0, ] # get survivors

rq.res <- rq(SCORE12 ~ trt, tau = 0.5, data = all_surs_w, weights = w)
rq.sum <- summary(rq.res, se = "iid")

vc.coe <- rq.sum$coefficients
CI <- cbind(vc.coe[, 1] - 1.96 * vc.coe[, 2], vc.coe[, 1] + 1.96 * vc.coe[, 2])  #calculate conservative CI

### trt=1
all_surs_w$trt <- abs(all_surs_w$trt - 1)

rq.res <- rq(SCORE12 ~ trt, tau = 0.5, data = all_surs_w, weights = w)
rq.sum <- summary(rq.res, se = "iid")

vc.coe <- rq.sum$coefficients
CI <- cbind(vc.coe[, 1] - 1.96 * vc.coe[, 2], vc.coe[, 1] + 1.96 * vc.coe[, 2])


###############
##### Quantile regression for confidence interval: change from baseline
if (interest == "change") {
  ### trt=0
  all_subs_w <- rbind(trt0_all, trt1_all)
  rq.res <- rq(delta ~ trt, tau = 0.5, data = all_subs_w, weights = w)
  rq.sum <- summary(rq.res, se = "iid")
  rq.sum

  vc.coe <- rq.sum$coefficients
  CI <- cbind(vc.coe[, 1] - 1.96 * vc.coe[, 2], vc.coe[, 1] + 1.96 * vc.coe[, 2])

  ### trt=1, which now becomes
  all_subs_w$trt <- abs(all_subs_w$trt - 1)
  rq.res <- rq(delta ~ trt, tau = 0.5, data = all_subs_w, weights = w)
  rq.sum <- summary(rq.res, se = "iid")
  rq.sum

  vc.coe <- rq.sum$coefficients
  CI <- cbind(vc.coe[, 1] - 1.96 * vc.coe[, 2], vc.coe[, 1] + 1.96 * vc.coe[, 2])
}




