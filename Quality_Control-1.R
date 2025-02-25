#################################################################################
#################################################################################
###          R Code for Quality control for `Outlier' test reviewers
###
###                         Yujie Wu and Molin Wang
###
###             For any questions please contact Yujie Wu
###                  Email: yujiewu@hsph.harvard.edu
#################################################################################
#################################################################################

library(readxl)
library(ggplot2)
library(ggthemes)
library(MASS)
library(gee)
################################################################################
################################################################################
########## Part 1. functions needed in the two-stage quality control ###########
################################################################################
################################################################################
#### Power formula used in the two-stage quality control procedure.
f.uniroot.noncentral.chisq <- function(x, other){
  1 - pchisq(qchisq(1-x, df=1), df = 1, ncp = other[1]^2) - other[2]
}


################### fdr.power.chisq() ##########################################
#### Calculates the estimated FDR under a range of powers,
#### where the hypotheses are compared with untruncated mean.
####
#### param@model.fit:   Contains the test reviewers' names (1st column),
####   and the estimated coefficients (2nd column).
####
#### param@Sigma.22:    Contains the estimated variance-covariance matrix of the 
####   coefficient estimates to test reviewers.
####
#### param@alternative: The alternative hypothesis.
####
#### output@results:    The range of power: (0.1,0.95), and corresponding estimated
####   FDR.
#################################################################################
fdr.power.chisq   <- function(model.fit, Sigma.22, alternative){
  num_tester       <- nrow(model.fit)
  power.candidate  <- seq(from = 0.1, to = 0.95, by = 0.01)
  results          <- c()
  for(power.fixed in power.candidate){
    results.temp        <- rep(NA, 2)
    names(results.temp) <- c("FDR.est", "Power")
    #### get the coeffcients for test reviewers
    tester_coef         <- matrix(model.fit[,2], ncol = 1)
    #### Calculate the test-specific significance levels and p-values
    pvalue              <- matrix(NA, ncol = 2, nrow = length(tester_coef))
    colnames(pvalue)    <- c("Significance Level", "P-value")
    for(j in 1:num_tester){
      C                 <- matrix(rep(-1/num_tester, num_tester), nrow = 1)
      C[j]              <- (num_tester-1)/num_tester
      screen.var        <- C %*% Sigma.22 %*% t(C)
      screen.test.stat  <- (C%*%tester_coef)^2/(C %*% Sigma.22 %*% t(C))
      screen.alpha      <- uniroot(f.uniroot.noncentral.chisq,c(0,1),
                                   other=c( alternative/sqrt(screen.var), power.fixed))$root
      pvalue[j,1]       <- screen.alpha
      pvalue[j,2]       <- 1- pchisq(screen.test.stat, df = 1)
    }
    
    results.temp["Power"] <- power.fixed
    
    #### get the test-specifc significance levels
    #### and p-values 
    Alpha.U               <- pvalue[,1]
    p.value               <- pvalue[,2]
    #### Estimate the FDR
    if(sum(p.value <= Alpha.U)!=0){
      Rm.U           <- sum( p.value <= Alpha.U )
      FDR.hat.naive  <- sum(Alpha.U)/Rm.U
      results.temp["FDR.est"] <- ifelse(FDR.hat.naive > 1, 1, FDR.hat.naive)
    } else{
      results.temp["FDR.est"] <- 0
    }
    
    results        <- rbind(results, results.temp)
  }
  return(list(results = results))
}

###################### fdr.power.chisq.trunc.mean() ###################################
#### calculates the estimated FDR under a range of powers,
#### where the hypothesis is compared with \delta*100% truncated mean.
####
#### param@model.fit:   Contains the test reviewers' names (1st column),
####   and the estimated coefficients (2nd column).
####
#### param@Sigma.22:    Contains the estimated variance-covariance matrix of the 
####   coefficient estimates to test reviewers.
####
#### param@alternative: The alternative hypothesis.
####
#### param@trunc.prop:  The truncated proportion for calculating the 
####   \delta*100% truncated mean.
####
####
#### output@results:    The range of power:(0.1,0.95), and corresponding estimated FDR.
#######################################################################################
fdr.power.chisq.trunc.mean   <- function(model.fit, Sigma.22, alternative, trunc.prop){
  num_tester       <- nrow(model.fit)
  power.candidate  <- seq(from = 0.1, to = 0.95, by = 0.01)
  results          <- c()
  for(power.fixed in power.candidate){
    results.temp        <- rep(NA, 2)
    names(results.temp) <- c("FDR.est", "Power")
    #### Fit linear regression directly
    tester_coef         <- matrix(model.fit[,2], ncol = 1)

    num.discard         <- round(num_tester*trunc.prop)
    discard.min         <- order(tester_coef)[1:num.discard]
    discard.max         <- order(tester_coef)[(num_tester-num.discard+1):num_tester]
    discard             <- c(discard.min, discard.max)
    
    pvalue              <- matrix(NA, ncol = 2, nrow = length(tester_coef))
    colnames(pvalue)    <- c("Significance Level", "P-value")
    for(j in 1:num_tester){
      C                 <- matrix(rep(0, num_tester), nrow = 1)
      C[-discard]       <- -1/(num_tester-2*num.discard)
      C[j]              <- C[j] + 1
      screen.var        <- C %*% Sigma.22 %*% t(C)
      screen.test.stat  <- (C%*%tester_coef)^2/(C %*% Sigma.22 %*% t(C))
      screen.alpha      <- uniroot(f.uniroot.noncentral.chisq,c(0,1),
                                   other=c( alternative/sqrt(screen.var), power.fixed))$root
      pvalue[j,1]       <- screen.alpha
      pvalue[j,2]       <- 1- pchisq(screen.test.stat, df = 1)
    }
    
    results.temp["Power"] <- power.fixed

    Alpha.U               <- pvalue[,1] 
    p.value               <- pvalue[,2]
    
    if(sum(p.value <= Alpha.U)!=0){
      Rm.U           <- sum( p.value <= Alpha.U )
      FDR.hat.naive  <- sum( Alpha.U)/Rm.U
      results.temp["FDR.est"] <- ifelse(FDR.hat.naive > 1, 1, FDR.hat.naive)
    } else{
      results.temp["FDR.est"] <- 0
    }
    
    results        <- rbind(results, results.temp)
  }
  return(list(results = results))
}

#################################### detect.real() #########################################
#### Function to detect `outlier' test reviewers, 
####   where the tests are compared with untruncated mean.
#### Both the unadjusted and adjusted procedure are performed.
#### 
#### param@model.fit:   Contains the test reviewers' names (1st column),
####   and the estimated coefficients (2nd column).
####
#### param@Sigma.22:    Contains the estimated variance-covariance matrix of the 
####   coefficient estimates to test reviewers.
####
#### param@power.fixed: The power selected based on the FDR vs. Power decision plot.
####   (It can either be selected directly or through a prespecified FDR)
####
#### param@alternative: The alternative hypothesis.
####
####
#### output@detected.tester: Detected `outlier' test reviewers (without FDR-driven adjustment).
####
#### output@coef.estimate:   Coefficients corresponding to these detected `outlier'
####   test reviewers (without FDR-driven adjustment).
####
#### output@detected.tester.adj: Detected `outlier' test reviewers (with FDR-driven adjustment).
####
#### output@coef.estimate.adj:   Coefficients corresponding to these detected `outlier'
####   test reviewers (with FDR-driven adjustment).
################################################################################################

detect.real    <- function(model.fit, Sigma.22, power.fixed, alternative){
  num_tester        <- nrow(model.fit)
  tester_coef       <- matrix(model.fit[,2], ncol = 1)
  pvalue            <- matrix(NA, ncol = 2, nrow = length(tester_coef))
  colnames(pvalue)  <- c("Significance Level", "P-value")
  
  for(j in 1:num_tester){
    C                <- matrix(rep(-1/num_tester, num_tester), nrow = 1)
    C[j]             <- (num_tester-1)/num_tester
    screen.var       <- C %*% Sigma.22 %*% t(C)
    screen.test.stat <- (C%*%tester_coef)^2/(C %*% Sigma.22 %*% t(C))
    screen.alpha     <- uniroot(f.uniroot.noncentral.chisq,c(0,1),
                                other=c( alternative/sqrt(screen.var), power.fixed))$root
    pvalue[j,1]      <- screen.alpha
    pvalue[j,2]      <- 1- pchisq(screen.test.stat, df = 1)
  }
  
  abnormal.tester    <- which(pvalue[,2] <= pvalue[,1])
  detect.tester      <- model.fit[abnormal.tester,1]
  
  Alpha.U            <- pvalue[,1] 
  ############ FDR Estimation
  p.value            <- pvalue[,2]
  if(sum(p.value <= Alpha.U)!=0){
    Rm.U           <- sum( p.value <= Alpha.U )
    FDR.hat        <- sum(Alpha.U)/Rm.U
    FDR.hat        <- ifelse(FDR.hat > 1, 1, FDR.hat)
  } else{
    FDR.hat        <- 0
  }
  
  num.discard     <- round(FDR.hat * length(abnormal.tester))
  if( num.discard > 0){
    if(length(abnormal.tester) > 1){
      abnormal.pvalue           <- pvalue[abnormal.tester,]
      abnormal.tester.adj       <- abnormal.tester[-order(abnormal.pvalue[,2],
                                                   decreasing = TRUE)[1:num.discard] ]
      detect.tester.adj         <- model.fit[abnormal.tester.adj,1]
    } else {
      abnormal.pvalue           <- pvalue[abnormal.tester,]
      abnormal.tester.adj       <- abnormal.tester[-order(abnormal.pvalue[2],decreasing = TRUE)[1:num.discard] ]
      detect.tester.adj         <- model.fit[abnormal.tester.adj,1]
    }   
  } else{
    abnormal.tester.adj         <- abnormal.tester
    detect.tester.adj           <- detect.tester
  }
  
  return( list(detected.tester      = detect.tester, 
               detected.tester.adj  = detect.tester.adj,
               coef.estimate        = tester_coef[abnormal.tester],
               coef.estimate.adj    = tester_coef[abnormal.tester.adj]))
}

################################ detect.trunc.real() #####################################
#### Function to detect `outlier' test reviewers,
####   where the tests are compared with \deleta*100% truncated mean.
#### Both the unadjusted and adjusted procedure are performed.
#### 
#### param@model.fit:   Contains the test reviewers' names (1st column),
####   and the estimated coefficients (2nd column).
####
#### param@Sigma.22:    Contains the estimated variance-covariance matrix of the 
####   coefficient estimates to test reviewers.
####
#### param@power.fixed: The power selected based on the FDR vs. Power decision plot.
####
#### param@alternative: The alternative hypothesis.
####
#### param#trunc.prop:  The truncated proportion for calculating the 
####   \delta*100% truncated mean.
####
####
#### output$detected.tester: Detected `outlier' test reviewers (without FDR-driven adjustment).
####
####       $coef.estimate:   Coefficients corresponding to these detected `outlier'
####   test reviewers (without FDR-driven adjustment).
####
####       $detected.tester.adj: Detected `outlier' test reviewers (with FDR-driven adjustment).
####
####       $coef.estimate.adj:   Coefficients corresponding to these detected `outlier'
####   test reviewers (with FDR-driven adjustment).
####
####################################################################################

detect.trunc.real <- function(model.fit, Sigma.22, trunc.prop, 
                              power.fixed, alternative){
  
  num_tester        <- nrow(model.fit)
  tester_coef       <- matrix(model.fit[,2], ncol = 1)
  pvalue            <- matrix(NA, ncol = 2, nrow = length(tester_coef))
  colnames(pvalue)  <- c("Significance Level", "P-value")
  
  num.discard       <- round(num_tester*trunc.prop)
  discard.min       <- order(tester_coef)[1:num.discard]
  discard.max       <- order(tester_coef)[(num_tester-num.discard+1):num_tester]
  discard           <- c(discard.min, discard.max)
  
  for(j in 1:num_tester){
    C                 <- matrix(rep(0, num_tester), nrow = 1)
    C[-discard]       <- -1/(num_tester-2*num.discard)
    C[j]              <- C[j] + 1
    screen.var        <- C %*% Sigma.22 %*% t(C)
    screen.test.stat  <- (C%*%tester_coef)^2/(C %*% Sigma.22 %*% t(C))
    screen.alpha      <- uniroot(f.uniroot.noncentral.chisq,c(0,1),
                                 other=c( alternative/sqrt(screen.var), power.fixed))$root
    pvalue[j,1]       <- screen.alpha
    pvalue[j,2]       <- 1- pchisq(screen.test.stat, df = 1)
  }  
  
  abnormal.tester    <- which(pvalue[,2] <= pvalue[,1])
  detect.tester      <- model.fit[abnormal.tester,1]
  
  Alpha.U            <- pvalue[,1] 
  ############ FDR Estimation
  p.value            <- pvalue[,2]
  if(sum(p.value <= Alpha.U)!=0){
    Rm.U           <- sum( p.value <= Alpha.U )
    FDR.hat        <- sum(Alpha.U)/Rm.U
    FDR.hat        <- ifelse(FDR.hat > 1, 1, FDR.hat)
  } else{
    FDR.hat        <- 0
  }
  
  num.discard     <- round(FDR.hat * length(abnormal.tester))
  if( num.discard > 0){
    if(length(abnormal.tester) > 1){
      abnormal.pvalue           <- pvalue[abnormal.tester,]
      abnormal.tester.adj       <- abnormal.tester[-order(abnormal.pvalue[,2],
                                                          decreasing = TRUE)[1:num.discard] ]
      detect.tester.adj         <- model.fit[abnormal.tester.adj,1]
    } else {
      abnormal.pvalue           <- pvalue[abnormal.tester,]
      abnormal.tester.adj       <- abnormal.tester[-order(abnormal.pvalue[2],decreasing = TRUE)[1:num.discard] ]
      detect.tester.adj         <- model.fit[abnormal.tester.adj,1]
    }   
  } else{
    abnormal.tester.adj         <- abnormal.tester
    detect.tester.adj           <- detect.tester
  }
  
  return( list(detected.tester     = detect.tester, 
               detected.tester.adj = detect.tester.adj,
               coef.estimate       = tester_coef[abnormal.tester],
               coef.estimate.adj   = tester_coef[abnormal.tester.adj]))
}


#################################################################################
#################################################################################
########## Part 2. Main function to perform two-stage quality control ###########
#################################################################################
#################################################################################

########################### QC_single_outcome.1() ###############################
####
#### QC_single_outcome.1() performs the first stage regression and output the 
####   FDR vs. Power decision plot
#### The function outputs $fdr_vs_power:      A range of power and correponding FDR.
####                      $fdr_vs_power_plot: FDR vs. Power decision plot.
####                      $coef.est:          Estimated coefficients to test reviewers.
####                      $cov.est:           Estimate variance-covariance matrix.
####                      $truncated, $trunc.prop, $alternative: record input parameters.
####
#### param@my.data:       Data.frame containing the data.
####
#### param@Y.name:        Test measurement on each study participant.
####
#### param@X.name:        Important predictors of the outcome, and confounders between 
####                        the test outcome and test reviewers.
####
#### param@Test_reviewer: The name of test reviewers. 
####
#### param@alternative:   The value of the alternative hypothesis.
####
#### param@truncated:     Either "Y" or "N". If "Y" is chosen, then the hypotheses are 
####                        performed to compare with truncated mean; if "N" is chosen,
####                        hypotheses are compared with untruncated mean.
####
#### param@trunc.prop:    If choose "Y" in truncated option, then you should specify 
####                        the truncation proportion, which should be between (0, 0.5).
####
######################################################################################

QC_single_outcome.1 <- function(my.data, Y.name, X.name, Test_reviewer, alternative, 
                                truncated, trunc.prop = NA){
  #### First stage no intercept linear regression
  regressors      <- c(X.name, Test_reviewer)
  X.length        <- length(X.name)
  Tester.length   <- length(Test_reviewer)
  Full.length     <- length(regressors)
  regressor.part  <- paste0(regressors, collapse = "+")
  formula.full    <- as.formula(paste0(Y.name,"~ -1 + ",regressor.part))
  lm.fit          <- lm(formula = formula.full, data = my.data)
  index.Tester    <- (X.length + 1) : Full.length
  coef.est        <- data.frame(Test_reviewer = Test_reviewer  , 
                                Estimate      = summary(lm.fit)$coefficient[index.Tester,1])
  cov.est         <- vcov(lm.fit)[index.Tester,index.Tester]
  
  if(truncated == "Y"){
    re           <- fdr.power.chisq.trunc.mean(model.fit = coef.est, Sigma.22 = cov.est,
                                               alternative = alternative, trunc.prop = trunc.prop)
    fdr_vs_power <- re$results
  } else{
    re           <- fdr.power.chisq(model.fit = coef.est, Sigma.22 = cov.est,
                                    alternative = alternative)
    fdr_vs_power <- re$results
  }
  fdr_vs_power   <- data.frame(fdr_vs_power)
  rownames(fdr_vs_power) <- 1:nrow(fdr_vs_power)
  
  fdr_vs_power_plot <- ggplot(fdr_vs_power, aes(x = Power, y = FDR.est))+
                         geom_smooth(se = FALSE)+
                         theme_bw()+
                         ggtitle("FDR vs. Power Decision Plot")
  print(fdr_vs_power_plot)
  
  return(list(fdr_vs_power      = fdr_vs_power, 
              fdr_vs_power_plot = fdr_vs_power_plot, 
              coef.est          = coef.est,
              cov.est           = cov.est,
              truncated         = truncated,
              trunc.prop        = trunc.prop,
              alternative       = alternative))
}


################################ QC_multiple_outcome.1() ################################
####
#### QC_multiple_outcome.1() performs the first stage regression using GEE and output the 
####   FDR vs. Power decision plot. 
#### The function outputs $fdr_vs_power:      a range of power and correponding FDR.
####                      $fdr_vs_power_plot: FDR vs. Power decision plot.
####                      $coef.est:          Estimated coefficients to test reviewers.
####                      $cov.est:           Estimate variance-covariance matrix.
####                      $truncated, $trunc.prop, $alternative: record input parameters.
####
#### param@my.data:       Data.frame containing the data.
####
#### param@Y.name:        Test measurement on each study participant.
####
#### param@X.name:        Important predictors of the outcome, and confounders between 
####                        the test outcome and test reviewers.
####
#### param@Test_reviewer: The name of test reviewers.
####
#### param@id:            id of a patient.
####
#### param@corstr:        Working covariance structure. The following are permitted:
####                        "independence", "exchangeable", "ar1", 
####                        "unstructured" and "userdefined".
####
#### param@alternative:   The value of the alternative hypothesis.
####
#### param@truncated:     Either "Y" or "N". If "Y" is chosen, then the hypotheses are 
####                        performed to compare with truncated mean; if "N" is chosen,
####                        hypotheses are compared with untruncated mean.
####
#### param@trunc.prop:    If chosen "Y" in truncated option, then you should specify 
####                        the truncation proportion, which should between (0, 0.5).
####
######################################################################################

QC_multiple_outcome.1 <- function(my.data, Y.name, X.name, Test_reviewer, alternative, 
                                  id, corstr, truncated, trunc.prop = NA){
  #### First stage no intercept linear regression
  regressors      <- c(X.name, Test_reviewer)
  X.length        <- length(X.name)
  Tester.length   <- length(Test_reviewer)
  Full.length     <- length(regressors)
  regressor.part  <- paste0(regressors, collapse = "+")
  formula.full    <- as.formula(paste0(Y.name,"~ -1 + ",regressor.part))
  gee.fit         <- gee(formula = formula.full,id = id, corstr = corstr,
                        data = my.data)
  index.Tester    <- (X.length + 1) : Full.length
  coef.est        <- data.frame( Test_reviewer = Test_reviewer,
                                 Estimate      = gee.fit$coefficients[index.Tester])
  cov.est         <- gee.fit$robust.variance[index.Tester,index.Tester]
  
  if(truncated == "Y"){
    re           <- fdr.power.chisq.trunc.mean(model.fit = coef.est, Sigma.22 = cov.est,
                                               alternative = alternative, trunc.prop = trunc.prop)
    fdr_vs_power <- re$results
  } else{
    re           <- fdr.power.chisq(model.fit = coef.est, Sigma.22 = cov.est,
                                    alternative = alternative)
    fdr_vs_power <- re$results
  }
  fdr_vs_power   <- data.frame(fdr_vs_power)
  rownames(fdr_vs_power) <- 1:nrow(fdr_vs_power)
  
  fdr_vs_power_plot <- ggplot(fdr_vs_power, aes(x = Power, y = FDR.est))+
    geom_smooth(se = FALSE)+
    theme_bw()+
    ggtitle("FDR vs. Power Decision Plot")
  print(fdr_vs_power_plot)
  
  return(list(fdr_vs_power      = fdr_vs_power, 
              fdr_vs_power_plot = fdr_vs_power_plot, 
              coef.est          = coef.est,
              cov.est           = cov.est,
              truncated         = truncated,
              trunc.prop        = trunc.prop,
              alternative       = alternative))
}

################################ QC_stage_2() ################################################
####
#### QC_stage_2() detects `outlier' test reviewers based on the results from
####  QC_single_outcome.1() or QC_multiple_outcome.1().
#### The detected test reviewers are reported from both the unadjusted and 
####  adjusted procedure.
####
#### output$detected.tester: Detected `outlier' test reviewers (without FDR-driven adjustment).
####
####       $coef.estimate:  Coefficients corresponding to these detected `outlier'
####   test reviewers (without FDR-driven adjustment).
####
####       $detected.tester.adj: Detected `outlier' test reviewers (with FDR-driven adjustment).
####
####       $coef.estimate.adj:   Coefficients corresponding to these detected `outlier'
####   test reviewers (with FDR-driven adjustment).
###
#### param@first.stage:  The results returned by QC_single_outcome.1() 
####                       or QC_multiple_outcome.1().
####
#### !!!! For the following two parameters, either specify FDR or Power, NO NEED to specify both !!!!
####
#### param@FDR:          Specify how large the FDR you want to have if you are more interested 
####                       in controling FDR. Notice that the choice of FDR should be based on 
####                       the decision plot.
####
#### param@Power:        Specify how large power you want to have if you are more interested 
####                       in ensuring enough power for the tests. Our function currently can 
####                       allow Power to be (0.1, 0.95).
##############################################################################################

QC_stage_2 <- function(first.stage, FDR = NA, Power = NA){
  
  fdr_vs_power <- first.stage$fdr_vs_power
  coef.est     <- first.stage$coef.est
  cov.est      <- first.stage$cov.est
  truncated    <- first.stage$truncated
  trunc.prop   <- first.stage$trunc.prop
  alternative  <- first.stage$alternative
  
  if(!is.na(FDR)){
    
    loess.fit  <- loess(FDR.est ~ Power, fdr_vs_power)
    find.power <- function(x, loess.fit, FDR){
      predict(loess.fit, data.frame(Power = x)) - FDR
    }
    Power    <-  uniroot(find.power, interval = c(0.1,0.95), loess.fit = loess.fit, FDR = FDR)$root
    
    if(truncated == "Y"){
      second.sta.re <- detect.trunc.real(model.fit   = coef.est, Sigma.22 = cov.est,
                                         trunc.prop  = trunc.prop, power.fixed = Power,
                                         alternative = alternative)
    } else{
      second.sta.re <- detect.real(model.fit   = coef.est, Sigma.22 = cov.est,
                                   power.fixed = Power, alternative = alternative)
    }
    
  } else{
    
    if(truncated == "Y"){
      second.sta.re <- detect.trunc.real(model.fit   = coef.est, Sigma.22 = cov.est,
                                         trunc.prop  = trunc.prop, power.fixed = Power,
                                         alternative = alternative)
    } else{
      second.sta.re <- detect.real(model.fit   = coef.est, Sigma.22 = cov.est,
                                   power.fixed = Power, alternative = alternative)
    }
    
  }
  
  return(second.sta.re)
}

#################################################################################
#################################################################################
################################ Part 3. Example ################################
#################################################################################
#################################################################################

####### Single outcome ######
####### The data for single outcome were generated following the simulation
####### study in the main paper
gene.data <- function(num_audio,          patient, 
                      non_zero_audio,     abn.audio.mean,
                      medium.mean,        medium_number ,
                      nor.audio.mean,     age.mean,       
                      age.std,            goodhear.p,   
                      littletrouble.p,    coef_age, 
                      coef_age.sq,        coef_goodhear, 
                      coef_littletrouble, RMSE){
  
  age            <- round(rnorm(num_audio*patient, age.mean, age.std),0)
  age.sq         <- age^2
  goodhear       <- rbinom(num_audio*patient, 1, goodhear.p)
  littletrouble  <- rbinom(num_audio*patient, 1, littletrouble.p)
  audiologist <- matrix(0, ncol = num_audio, nrow = patient*num_audio)
  colnames(audiologist) <- paste0("audiologist_", 1:num_audio)
  for(i in 1:num_audio){
    index                  <- (patient*(i-1)+1):(patient*i)
    audiologist[,i][index] <- 1
  }
  XX             <- cbind(age, age.sq, goodhear, littletrouble, audiologist)
  
  #### give coefficient values to age, self_report and non-zero audiologists
  
  coef_audio     <- rep(NA, num_audio) 
  index.sig      <- 1:(non_zero_audio + medium_number)
  sig.audio      <- paste0("audiologist_", index.sig)
  nosig.audio    <- paste0("audiologist_", c(1:num_audio)[-index.sig])
  
  abn.audio.coef <- rep(abn.audio.mean, non_zero_audio)
  
  coef_audio[index.sig]   <- c(abn.audio.coef, rep(medium.mean, medium_number))
  coef_audio[-index.sig]  <- rep(nor.audio.mean, (num_audio-non_zero_audio-medium_number))
  
  coef_all   <- c(coef_age, coef_age.sq, coef_goodhear, 
                  coef_littletrouble, coef_audio)
  
  #### generating outcome variable Y, with random noise N(0,1)
  YY         <- XX %*% coef_all + rnorm(num_audio*patient, 0, RMSE)
  df.XY      <- data.frame(XX, Y=YY)
  
  return(list(df.XY = df.XY,
              coef  = coef_all))
}
data <- gene.data(     num_audio          = 100,         patient = 40, 
                       non_zero_audio     = 5,    abn.audio.mean = 75.10, 
                       medium.mean        = 70.10, medium_number = 3,         
                       nor.audio.mean     = 66.95,      age.mean = 56.56,                
                       age.std            = 4.36,     goodhear.p = 0.4358,       
                       littletrouble.p    = 0.2519,     coef_age = -2.73,        
                       coef_age.sq        = 0.03,  coef_goodhear = 0.03, 
                       coef_littletrouble = 3.32,           RMSE = 10)$df.XY

#### Compare with 10% truncated mean.
single.trunc.stage1 <- QC_single_outcome.1(my.data = data, Y.name = "Y", X.name = c("age", "age.sq", "goodhear", "littletrouble"),
                           Test_reviewer = colnames(data)[-c(1:4,105)], alternative = 5,
                           truncated = "Y", trunc.prop = 0.1)
#### We can either choose an FDR or Power to perform 
#### hypothesis to detect `outlier' audiologists.
#### Note that when choosing FDR, it should fall between the range of 
#### the FDR given by the FDR vs. Power decision plot.
QC_stage_2(single.trunc.stage1, FDR = 0.2)
QC_stage_2(single.trunc.stage1, Power = 0.8)


#### Compare with untruncated mean
single.untrunc.stage1 <- QC_single_outcome.1(my.data = data, Y.name = "Y", X.name = c("age", "age.sq", "goodhear", "littletrouble"),
                           Test_reviewer = colnames(data)[-c(1:4,105)], alternative = 5,
                           truncated = "N")
QC_stage_2(single.untrunc.stage1, FDR = 0.2)
QC_stage_2(single.untrunc.stage1, Power = 0.8)



####### Multiple outcomes ######
gene.data.gee <- function(num_audio,          patient, 
                          non_zero_audio,     abn.audio.mean,
                          medium.mean,        medium_number,
                          nor.audio.mean,     age.mean, 
                          age.std,            goodhear.p,     
                          littletrouble.p,    coef_age, 
                          coef_age.sq,        coef_goodhear, 
                          coef_littletrouble, RMSE, rho){
  
  age            <- round(rnorm(num_audio*patient, age.mean, age.std),0)
  age.sq         <- age^2
  goodhear       <- rbinom(num_audio*patient, 1, goodhear.p)
  littletrouble  <- rbinom(num_audio*patient, 1, littletrouble.p)
  audiologist <- matrix(0, ncol = num_audio, nrow = patient*num_audio)
  colnames(audiologist) <- paste0("audiologist_", 1:num_audio)
  for(i in 1:num_audio){
    index                  <- (patient*(i-1)+1):(patient*i)
    audiologist[,i][index] <- 1
  }
  XX             <- cbind(age, age.sq, goodhear, littletrouble, audiologist)
  
  #### give coefficient values to age, self_report and non-zero audiologists
  
  coef_audio     <- rep(NA, num_audio) 
  index.sig      <- 1:(non_zero_audio + medium_number)
  sig.audio      <- paste0("audiologist_", index.sig)
  nosig.audio    <- paste0("audiologist_", c(1:num_audio)[-index.sig])
  
  abn.audio.coef <- rep(abn.audio.mean, non_zero_audio)
  
  coef_audio[index.sig]   <- c(abn.audio.coef, rep(medium.mean, medium_number))
  coef_audio[-index.sig]  <- rep(nor.audio.mean, (num_audio-non_zero_audio-medium_number))
  
  coef_all   <- c(coef_age, coef_age.sq, coef_goodhear, 
                  coef_littletrouble, coef_audio)
  
  #### generating outcome variable Y, with random noise N(0,1)
  mu         <- XX %*% coef_all 
  XX         <- cbind(XX, id = 1:(num_audio*patient))
  mu         <- rep(mu, each = 2)
  Y          <- c()
  for(pair in 1:(length(mu)/2)){
    index      <- (2*(pair-1)+1) : (2*pair)
    YY         <- mvrnorm(n = , mu = mu[index], Sigma = matrix(c(RMSE^2, rho*RMSE^2,
                                                                 rho*RMSE^2, RMSE^2), byrow = TRUE,
                                                                 nrow = 2, ncol = 2))
    id         <- pair
    YY         <- c(YY, pair)
    Y          <- rbind(Y, YY)
  }
  Y.long       <- data.frame(Y  = c(Y[,1], Y[,2]),
                             id = c(Y[,3], Y[,3]))
  Y.long       <- Y.long[order(Y.long$id),]
  
  df.XY        <- merge(Y.long, XX, by.x = "id")

  return(list(df.XY = df.XY,
              coef  = coef_all))
}

data.gee <- gene.data.gee(     num_audio          = 100,         patient = 40, 
                               non_zero_audio     = 5,    abn.audio.mean = 75.10, 
                               medium.mean        = 70.10, medium_number = 3,   
                               nor.audio.mean     = 66.95,      age.mean = 56.56,                
                               age.std            = 4.36,     goodhear.p = 0.4358,       
                               littletrouble.p    = 0.2519,     coef_age = -2.73,        
                               coef_age.sq        = 0.03,  coef_goodhear = 0.03, 
                               coef_littletrouble = 3.32,           RMSE = 10, rho = 0.3)$df.XY
#### First plot the FDR vs. Power decision plot
multiple.trunc.stage1 <- QC_multiple_outcome.1(my.data = data.gee, Y.name = "Y", X.name = c("age", "age.sq", "goodhear", "littletrouble"),
                             Test_reviewer = colnames(data.gee)[-c(1:6)], alternative = 5,
                             id = "id", corstr = "exchangeable", truncated = "Y", trunc.prop = 0.1)
#### Based on the decision plot, we can either choose a power or an FDR to
#### reject the null hypothesis.
QC_stage_2(multiple.trunc.stage1, FDR = 0.2)
QC_stage_2(multiple.trunc.stage1, Power = 0.8)

multiple.untrunc.stage1 <- QC_multiple_outcome.1(my.data = data.gee, Y.name = "Y", X.name = c("age", "age.sq", "goodhear", "littletrouble"),
                                                 Test_reviewer = colnames(data.gee)[-c(1:6)], alternative = 5,
                                                 id = "id", corstr = "exchangeable", truncated = "N")
QC_stage_2(multiple.untrunc.stage1, FDR = 0.2)
QC_stage_2(multiple.untrunc.stage1, Power = 0.8)
