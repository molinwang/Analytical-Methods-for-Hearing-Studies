library(MASS)
library(gee)
library(readr)
library(deployer)

oneres<-read.csv("oneres.csv", header=T)

#########################################################################################################################################################
#Exposure:DASH score quartiles (vq, 1/2/3/4)
#Outcome:hearing threshold 5+ dB decline for low, middle, high PTA
#        (generate chg_r_l_ge5/chg_l_l_ge5, chg_r_m_ge5/chg_l_m_ge5, chg_r_h_ge5/chg_l_h_ge5, 0/1)
#        worse/better/left/right ear outcome are derived from the above variables
#Other covariates:
#1.age(age_3y, yrs) 
#2.race(white/black/multi/other or unknown, 0/1) 
#3.BMI(<25/25-29/30-34/35-39/40+/missing, 0/1)
#4.smoking status(current/past/never/missing, 0/1)
#5.cumulative-avereage energy intake(calor11v, kcal/day)
#6.loud noise exposure(vln, 0/1)
#7.baseline PTA for low, middle and high frequencies(lptar_bl/lptal_bl, mptar_bl/mptal_bl, hptar_bl/hptal_bl, dB)
#  baseline PTA corresponds to the outcome in different models
#8.tinnitus(tinnever, 0/1)
#########################################################################################################################################################

#Generate outcomes(right/left/worse/better for low/middle/high frequency)
oneres$chg_r_l_ge5<-ifelse(oneres$lptar_fu-oneres$lptar_bl>=5, 1, 0)
oneres$chg_l_l_ge5<-ifelse(oneres$lptal_fu-oneres$lptal_bl>=5, 1, 0)
oneres$chg_w_l_ge5<-ifelse(oneres$chg_r_l_ge5==0 & oneres$chg_l_l_ge5==0, 0, 1)
oneres$chg_b_l_ge5<-ifelse(oneres$chg_r_l_ge5==1 & oneres$chg_l_l_ge5==1, 1, 0)
oneres$chg_r_m_ge5<-ifelse(oneres$mptar_fu-oneres$mptar_bl>=5, 1, 0)
oneres$chg_l_m_ge5<-ifelse(oneres$mptal_fu-oneres$mptal_bl>=5, 1, 0)
oneres$chg_w_m_ge5<-ifelse(oneres$chg_r_m_ge5==0 & oneres$chg_l_m_ge5==0, 0, 1)
oneres$chg_b_m_ge5<-ifelse(oneres$chg_r_m_ge5==1 & oneres$chg_l_m_ge5==1, 1, 0)
oneres$chg_r_h_ge5<-ifelse(oneres$hptar_fu-oneres$hptar_bl>=5, 1, 0)
oneres$chg_l_h_ge5<-ifelse(oneres$hptal_fu-oneres$hptal_bl>=5, 1, 0)
oneres$chg_w_h_ge5<-ifelse(oneres$chg_r_h_ge5==0 & oneres$chg_l_h_ge5==0, 0, 1)
oneres$chg_b_h_ge5<-ifelse(oneres$chg_r_h_ge5==1 & oneres$chg_l_h_ge5==1, 1, 0)

######NEW ANALYSIS--WITHOUT INTERACTION######
#generate new DASH exposure using DASH quantile
summary(oneres$dash11v)
vq<-rep(0, nrow(oneres))
for (i in 1:length(vq)) {
  if(oneres$dash11v[i]<=21.5) {vq[i]<-1}
  else if(oneres$dash11v[i]>21.5 & oneres$dash11v[i]<=24.5) {vq[i]<-2}
  else if(oneres$dash11v[i]>24.5 & oneres$dash11v[i]<=27.5) {vq[i]<-3}
  else {vq[i]<-4}
}


#both ear-GEE
#create covariates
id<-rep(1:nrow(oneres), each=2)
age_3y<-rep(oneres$age_3y, each=2)
bmicat<-rep(oneres$bmicat, each=2)
calor11v<-rep(oneres$calor11v, each=2)
vq.gee<-rep(vq, each=2)
race<-rep(oneres$race, each=2)
smk11<-rep(oneres$smk11, each=2)
tinnever<-rep(oneres$tinnever, each=2)
vln<-rep(oneres$vln, each=2)
#create outcomes
both.l<-as.vector(t(cbind(oneres$chg_r_l_ge5, oneres$chg_l_l_ge5)))#low-frequency
both.m<-as.vector(t(cbind(oneres$chg_r_m_ge5, oneres$chg_l_m_ge5)))#middle-frequency
both.h<-as.vector(t(cbind(oneres$chg_r_h_ge5, oneres$chg_l_h_ge5)))#high-frequency
#create baseline
base.l<-as.vector(t(cbind(oneres$lptar_bl, oneres$lptal_bl)))#low-frequency
base.m<-as.vector(t(cbind(oneres$mptar_bl, oneres$mptal_bl)))#middle-frequency
base.h<-as.vector(t(cbind(oneres$hptar_bl, oneres$hptal_bl)))#high-frequency
#run GEE model
gee.pta<-function(outcome, baseline){
  gee.both<-gee(outcome~as.factor(vq.gee)+age_3y+as.factor(bmicat)+calor11v+
                  baseline+as.factor(race)+as.factor(smk11)+as.factor(tinnever)+as.factor(vln),
                family=binomial("logit"), id=id, corstr="exchangeable")
  var.gee<-diag(gee.both$robust.variance)
  con.low<-round(exp(gee.both$coefficients-qnorm(0.975)*sqrt(var.gee)),2)
  con.up<-round(exp(gee.both$coefficients+qnorm(0.975)*sqrt(var.gee)),2)
  p<-round(2*(1-pnorm(abs(coef(summary(gee.both))[,"Robust z"]))),5)
  result<-data.frame(round(exp(gee.both$coefficients),2),con.low, con.up, p)
  colnames(result)<-c("OR", "CI_low", "CI_up", "p-value")
  print(result)
}
gee.pta(both.l,base.l)#low-frequency
gee.pta(both.m,base.m)#middle-frequency
gee.pta(both.h,base.h)#high-frequency
table(rep(vq, each=2), both.l)#low-frequency case number
table(rep(vq, each=2), both.m)#middle-frequency case number
table(rep(vq, each=2), both.h)#high-frequency case number

######NEW ANALYSIS--WITH INTERACTION######
#create covariates
age_3y.l<-rep(oneres$age_3y, each=3)
bmicat.l<-rep(oneres$bmicat, each=3)
calor11v.l<-rep(oneres$calor11v, each=3)
race.l<-rep(oneres$race, each=3)
smk.l<-rep(oneres$smk11, each=3)
tinnever.l<-rep(oneres$tinnever, each=3)
vln.l<-rep(oneres$vln, each=3)
#create exposure
vq.l<-rep(vq, each=3)
vq.l[vq.l==1]<-0
#create frequency indicator--0:low, 1:mid, 3:high
freq.ind<-rep(c(0,1,3),nrow(oneres))
#logistic regression macro

#both ear method--GEE
#create covariates
id.g<-rep(1:nrow(oneres), each=6)
age_3y.g<-rep(oneres$age_3y, each=6)
bmicat.g<-rep(oneres$bmicat, each=6)
calor11v.g<-rep(oneres$calor11v, each=6)
race.g<-rep(oneres$race, each=6)
smk.g<-rep(oneres$smk11, each=6)
tinnever.g<-rep(oneres$tinnever, each=6)
vln.g<-rep(oneres$vln, each=6)
#create exposure
vq.g<-rep(vq, each=6)
vq.g[vq.g==1]<-0
#create outcome
both.g<-as.vector(t(cbind(oneres$chg_r_l_ge5, oneres$chg_l_l_ge5, oneres$chg_r_m_ge5, oneres$chg_l_m_ge5, oneres$chg_r_h_ge5, oneres$chg_l_h_ge5)))
#create baseline
both.b<-as.vector(t(cbind(oneres$lptar_bl, oneres$lptal_bl, oneres$mptar_bl, oneres$mptal_bl, oneres$hptar_bl, oneres$hptal_bl)))
#create frequency indicator--0:low, 1:mid, 3:high
freq.both<-rep(c(0, 0, 1, 1, 3, 3),nrow(oneres))
#GEE analysis
both.gee<-gee(both.g~as.factor(vq.g)+as.factor(freq.both)+as.factor(vq.g*freq.both)+age_3y.g+as.factor(bmicat.g)+calor11v.g+
                both.b+as.factor(race.g)+as.factor(smk.g)+as.factor(tinnever.g)+as.factor(vln.g),
              family=binomial("logit"), id=id.g, corstr="exchangeable")
summary(both.gee)
both.gee$robust.variance[1:12, 1:12]


