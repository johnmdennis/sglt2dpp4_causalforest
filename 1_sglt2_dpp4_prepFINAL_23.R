library(foreign)
library(tidyverse)  # data manipulation and visualization
library(modelr)     # provides easy pipeline modeling functions
library(broom)      # helps to tidy up model outputs
library(stringr)    # work with string variables
library(forcats)    # work with categorical vars
library(lubridate)  # work with dates
library(magrittr)
library(rms)
library(MatchIt)
library(stargazer)
library(dplyr)
library(tidyr)
require(gridExtra)
library(ggExtra)
library(ggpubr)
library(cowplot)
library(corrplot)
library(caret)
require(emmeans)
library(plyr)

######################################
#DPP4 and SGLT2 cohort
######################################

#Outcome:
#posthba1c_final6m - closest HbA1c to 6 month (range 3-15m)

rm(list=ls())
setwd("C:/Users/jmd237/OneDrive - University of Exeter/John/local/CPRD/cprd/CPRD_19/data/update")
load("finalcprd19_analysisreadyv2.Rda") # from: 1_respdatasetprepFINALv2.R
cprd_bckup <- cprd      

#Set posthba1c_final to 6m value
      cprd <- cprd %>% mutate(posthba1c_final=posthba1c_final6m)
      cprd$posthba1c_final6m <- NULL
      
#SGLT2 vs DPP4 > 2013
      cprd <- cprd_bckup %>% filter((drugclass=="DPP4" & yrdrugstart >= 2013 | drugclass=="SGLT2"))

#Are any duplicates i.e. started treatment on same day
      head(cprd)
      cprd$pated <- as.character(cprd$pated)
      px <- data.frame(unique(cprd$pated)) #YES
      
      cprddup <- cprd[cprd$pated %in% cprd$pated[duplicated(cprd$pated)],]
      cprddup$dummy2 <- 1
      cprddup <- cprddup %>% select(patid,dummy2)
      cprddup <- cprddup[!duplicated(cprddup$patid),]
      
      #keep only those starting monotherapy
      cprd <- merge(cprd,cprddup,all.x=TRUE,by="patid")
      table(cprd$dummy2)
      cprd <- cprd %>% filter(is.na(dummy2))
      
#gen clear male sex variable 
      cprd <- cprd %>% mutate(malesex=factor(sex))
      table(cprd$sex)
      
#gen ethnicity variables 
      cprd <- cprd %>% mutate(ethcat=ifelse(cprd$ethnicity5==1,"White",
                                                      ifelse(cprd$ethnicity5==2,"Mixed",
                                                             ifelse(cprd$ethnicity5==3,"Asian",
                                                                    ifelse(cprd$ethnicity5==4,"Black",
                                                                           ifelse(cprd$ethnicity5==5,"Other",
                                                                                  ifelse(cprd$ethnicity5=="6",NA,NA)))))))
      cprd$ethcat <- factor(cprd$ethcat)
      describe(cprd$ethcat)
      
      cprd <- cprd %>% mutate(nonwhite=ifelse(cprd$ethnicity5==1,"White",
                                                        ifelse(cprd$ethnicity5>1 & cprd$ethnicity5 < 6,"Non-White",NA)))
      cprd$nonwhite <- factor(cprd$nonwhite)
      describe(cprd$nonwhite)
      table(cprd$ethnicity5,cprd$nonwhite)
      
#Apply criteria from finalmerge script (Bev)  
    # drop if less than 61 days since started previous line of therapy
      cprd <- cprd %>% mutate(timeprevcom_less61=if_else(timeprevcomb<=61,1,NA_real_))          
      table(cprd$timeprevcom_less61) 
      cprd <- cprd %>% filter(is.na(timeprevcom_less61))
      cprd$timeprevcom_less61 <- NULL
      
#Exclude from this study if first-line treatment
      table(cprd$drugline)
      cprd <- cprd %>% filter(drugline != 1)
      
# drop if treated with INSULIN when starting new drug
      table(cprd$INS)
      cprd <- cprd %>% filter(INS==0)      

#Exclude from this study if egfr <45
      cprd <- cprd %>% mutate(egfr45= if_else(egfr_ckdepi<45,1,0))
      cprd <- cprd %>% mutate(egfr45= if_else(is.na(egfr45),0,egfr45))
      describe(cprd$egfr45)
      table(cprd$egfr45)
      
      cprd.test <- cprd %>% filter(egfr45==1)
      prop.table(table(cprd.test$drugclass))
      cprd <- cprd %>% filter(egfr45==0)

#Require HbA1c between 53 and 120 at baseline (<13%)
      cprd <- cprd %>% mutate(hb_extreme=if_else(prehba1cmmol<53,1,NA))
      describe(cprd$hb_extreme)
      cprd <- cprd %>% mutate(hb_extreme=if_else(prehba1cmmol>120,1,NA))
      describe(cprd$hb_extreme)
      
      cprd <- cprd %>% mutate(hb_extreme=if_else(prehba1cmmol<53|prehba1cmmol>120,1,0))
      describe(cprd$hb_extreme)
      table(cprd$hb_extreme)
      prop.table(table(cprd$hb_extreme)) #5%
      test <- cprd %>% filter(hb_extreme==1)
      describe(test$prehba1cmmol)
      hist(test$prehba1cmmol)
      cprd <- cprd %>% filter(hb_extreme==0)
      
#Require valid baseline and outcome HbA1c      
      describe(cprd$prehba1cmmol)
      describe(cprd$posthba1c_final)
      cprd <- cprd %>% filter(!is.na(prehba1cmmol))
      #cprd <- cprd %>% filter(!is.na(posthba1c_final))

      drop <- cprd %>% filter(is.na(posthba1c_final))
      describe(drop$timetochange)#
      drop$time91 <- ifelse(drop$timetochange<91,1,0)
      table(drop$time91)      
      
      drop91 <- drop %>% filter(time91==1)
      table(drop91$nextdrugchange)
      
      # Tabulation: Freq.   Numeric  Label
      #         1  Add
      #          2  Remove
      #         3  Swap
      #          4  Stop - break
      #          5  Stop - end of px
      
      drop91stop <- drop91 %>% filter(nextdrugchange==2)
      head(drop91stop)
      table(drop91$stopdrug3m_3mFU)
      describe(drop91stop$stopdrug3m_3mFU)
#Code drugline as factor with 5+ group
      table(cprd$drugline)
      cprd$drugline <- ifelse(cprd$drugline>5,5,cprd$drugline)
      cprd$drugline <- factor(cprd$drugline)
      levels(cprd$drugline)

      cprd$druglinena <- ifelse(cprd$druglinena>5,5,cprd$druglinena)
      cprd$druglinena <- factor(cprd$druglinena)
      levels(cprd$druglinena)
      
      describe(cprd$drugline)
      describe(cprd$druglinena)
      
#Define drugclass as 2 level factor
      cprd$drugclass <- as.factor(as.character(cprd$drugclass))
      levels(cprd$drugclass)
      table(cprd$drugclass)

#Define TRG/HDL ratio (both mmol/L originally, work out ratio in mg/dl)    https://www.omnicalculator.com/health/cholesterol-ratio
      cprd <- cprd %>% mutate(trghdlratio=(pretrig*88.5)/(prehdl*38.67))
      
#Patients with valid data on each therapy    
      # sglt2dpp4resp <- cprd[cprd$patid %in% cprd$patid[duplicated(cprd$patid)],]
      # table(sglt2dpp4resp$drugclass)   
      # save(sglt2dpp4resp,file="sglt2dpp4resp.Rda")
      
#Define cohorts based on trials
      # #Trial CANTATA-D2 (Bck MFN + SU)
      # #Inclusin HbA1c 53-91, exclusion eGFr <55, creatinine >124 (m) or 115 (f), FPG >16.7, T1D, CVD, uncontrolled HTN, tx with PPARy agonist, other glucose lowering drug
      # 
      # #Trial CANTATA-D (Bck MFN)
      # #Inclusin age 18-80, HbA1c 53-91, exclusion eGFr <55, creatinine >124 (m) or 115 (f), FPG >15, T1D, CVD (including myocardial infarction,
      # # unstable angina, revascularisation procedure or cerebrovascular
      # # acciden), uncontrolled HTN, tx with PPARy agonist, other glucose lowering drug
      # 
      # #Starting cohort 31215        
      # finalall.trial.cohort <- finalall %>% filter(prehba1cmmol>=53 & prehba1cmmol <= 91) #25761
      # summary(finalall.trial.cohort$pregluc)
      # finalall.trial.cohort <- finalall.trial.cohort %>% filter(pregluc<=15 | is.na(pregluc)) #23495
      # summary(finalall.trial.cohort$egfr)
      # finalall.trial.cohort <- finalall.trial.cohort %>% filter(egfr>=55 | is.na(egfr)) #21201
      # #HTN
      # #CVD
      # #NON MFN SU GLUCOSE LOWERING DRUG
      # table(as.character(finalall.trial.cohort$drugcombo))
      # finalall.trial.cohort <- finalall.trial.cohort %>% 
      #   filter(drugcombo=="MFNDPP4"|drugcombo=="MFNSGLT2"|drugcombo=="MFNSUDPP4"|drugcombo=="MFNSUSGLT2") #N15755
      # final$drugline <- factor(final$drugline)
      # m0 <- ols(iniresp15m~prehba1cmmol+drugclass+drugline,data=finalall.trial.cohort,x=TRUE,y=TRUE)
      # m0

      ####################
      #Centre scale variables
      ####################
      cs. <- function(x) scale(x,center=TRUE,scale=TRUE)
      hist(cprd$t2dmduration)
      cprd$t2dmdurationlog <- log(cprd$t2dmduration)
      hist(cprd$t2dmdurationlog)
      cprd$t2dmdurationcs <- as.numeric(cs.(cprd$t2dmduration))
      cprd$t2dmdurationcslog <- as.numeric(cs.(log(cprd$t2dmduration)))
      
      hist(cprd$agetx)
      cprd$agetxcs <- as.numeric(cs.(cprd$agetx))
      cprd$agetxcslog <- as.numeric(cs.(log(cprd$agetx)))
      
      hist(cprd$precreat)
      cprd$precreatcs <- as.numeric(cs.(cprd$precreat))
      cprd$precreatcslog <- as.numeric(cs.(log(cprd$precreat)))
      
      hist(cprd$egfr_ckdepi)
      cprd$egfr_ckdepics <- as.numeric(cs.(cprd$egfr_ckdepi))
      cprd$egfr_ckdepicslog <- as.numeric(cs.(log(cprd$egfr_ckdepi)))
      
      hist(cprd$prebmi)
      cprd$prebmi <- ifelse(cprd$prebmi<5,NA,cprd$prebmi)
      cprd$prebmics <- as.numeric(cs.(cprd$prebmi))
      cprd$prebmicslog <- as.numeric(cs.(log(cprd$prebmi)))
      hist(cprd$prebmicslog)
      
      hist(cprd$preweight)
      cprd$preweightcs <- as.numeric(cs.(cprd$preweight))
      cprd$preweightcslog <- as.numeric(cs.(log(cprd$preweight)))
      
      hist(cprd$preldl)
      cprd$preldlcs <- as.numeric(cs.(cprd$preldl))
      cprd$preldlcslog <- as.numeric(cs.(log(cprd$preldl)))
      
      hist(cprd$prehdl)
      cprd$prehdlcs <- as.numeric(cs.(cprd$prehdl))
      cprd$prehdlcslog <- as.numeric(cs.(log(cprd$prehdl)))
      
      hist(cprd$pretrig)
      cprd$pretrigcs <- as.numeric(cs.(cprd$pretrig))
      cprd$pretriglogcs <- as.numeric(cs.(log(cprd$pretrig)))
      hist(cprd$pretriglog)
      
      hist(cprd$prealt)
      cprd$prealtcs <- as.numeric(cs.(cprd$prealt))
      cprd$prealtlogcs <- as.numeric(cs.(log(cprd$prealt)))
      hist(cprd$prealtlog)
      
      hist(cprd$preast)
      cprd$preastcs <- as.numeric(cs.(cprd$preast))
      cprd$preastlogcs <- as.numeric(cs.(log(cprd$preast)))
      hist(cprd$preastlog)
      
      hist(cprd$pregluc)
      cprd$pregluccs <- as.numeric(cs.(cprd$pregluc))
      cprd$pregluclogcs <- as.numeric(cs.(log(cprd$pregluc)))
      hist(cprd$pregluclog)
      
      hist(cprd$prealb)
      cprd$prealbcs <- as.numeric(cs.(cprd$prealb))
      cprd$prealblogcs <- as.numeric(cs.(log(cprd$prealb)))
      
      hist(cprd$prebil)
      cprd$prebilcs <- as.numeric(cs.(cprd$prebil))
      cprd$prebillogcs <- as.numeric(cs.(log(cprd$prebil)))
      
      hist(cprd$preplatelets)
      cprd$preplateletscs <- as.numeric(cs.(cprd$preplatelets))
      cprd$preplateletslogcs <- as.numeric(cs.(log(cprd$preplatelets)))
      
      hist(cprd$trghdlratio)
      cprd$trghdlratiocs <- as.numeric(cs.(cprd$trghdlratio))
      cprd$trghdlratiologcs <- as.numeric(cs.(log(cprd$trghdlratio)))
      
######################################
#define hold back cohort
######################################
    #Random sample of 60% of patients on each drug
    # set.seed(8731)
    # cprd.dev <- cprd %>% group_by(drugclass) %>% sample_frac(.6)
    # cprd.dev <- data.frame(cprd.dev)
    # cprd.val <- subset(cprd, !(pateddrug %in% cprd.dev$pateddrug))
      
      
    #Temporal validation - subset starting treatment > 1/1/2018
    # cprd <- cprd %>% mutate(testcohort=if_else(datedrug>=as.Date("2018-01-01"),1,0))
    # table(cprd$testcohort)
    # table(cprd$drugclass,cprd$testcohort)
    # prop.table(table(cprd$drugclass,cprd$testcohort),margin=1)
    # #Jan 18 gives 10% of DPP4 px and 19% of SGLT2 px in test cohort, Jul 18 gives 6% of DPP4 and 11% of SGLT2 px in test cohort
    # cprd.dev <- cprd %>% dplyr::filter(testcohort==0)
    # cprd.val <- cprd %>% dplyr::filter(testcohort==1)
    
####################
#Final cohorts
####################
  # #final dev cohort
  #   final.dev <- cprd.dev
  #   #final.dev <- cprd %>% filter(testcohort == 0)
  #   table(final.dev$drugclass)
  #   table(final.dev$drugline,final.dev$drugclass)
  # 
  # #final validation cohort
  #   final.val <- cprd.val
  #   #final.val <- cprd %>% filter(testcohort == 1)
  #   table(final.val$drugclass)

  #describe val cohort and dev cohort
    final.all <- cprd
    final.all$drug <- as.character(final.all$drugclass)
    final.all$fldrug <- as.factor(as.character(final.all$fldrug))
    
####################
#Final cohort desc
####################
   # #describe val cohort and dev cohort
   #  finalall <- cprd
   #  finalall$drug <- as.character(finalall$drugclass)
   #  finalall$fldrug <- as.factor(as.character(finalall$fldrug))
   # 
   #  ## List variables for the table
   #  Vars <- c("drugline","druglinena","ncurrtx","yrdrugstart","agedx","agetx","gender","malesex","dpp4subtype","sglt2subtype","t2dmduration","ethnicity5","adherence_t")
   #  ## List of variables which are factors
   #  FactorVars <- c("drugline","druglinena","ncurrtx","gender","malesex","dpp4subtype","sglt2subtype","ethnicity5")
   #  ## BY MATCHED / ID Create the table
   #  DemographicsTable <- CreateTableOne(vars=Vars, strata=c("drug"),factorVars = FactorVars, data=finalall)
   #  ## Export as an html document
   #  stargazer(print(DemographicsTable, dropEqual = T)) #, out = "DemographicsTableStrat.html")
   #  
   #  ## List variables for the table
   #  Vars <- c("precreat","egfr","egfr_ckdepi","egfr_cg",
   #            "preldl","prehdl","pretrig","preweight","prebmi",
   #            "preast","prealt","preplatelets","prebil","prealb","api","trghdlratio")    ## List of variables which are factors
   #  ## BY MATCHED / ID Create the table
   #  Table <- CreateTableOne(vars=Vars, strata=c("drug"),factorVars = FactorVars, data=finalall)
   #  ## Export as an html document
   #  stargazer(print(Table, dropEqual = T))
   # 
   #  ## List variables for the table
   #  Vars <- c("prehba1cmmol","pregluc","prehba1cchange","iniresp_final","posthba1c_final") #"prehba1cmmol","posthba1c6mmmol","posthba1c12mmmol","iniresp6m","iniresp12m","iniresp6mmmol","iniresp12mmmol"
   #  ## BY MATCHED / ID Create the table
   #  Table <- CreateTableOne(vars=Vars, strata=c("drug"),factorVars = FactorVars, data=finalall)
   #  ## Export as an html document
   #  stargazer(print(Table, dropEqual = T))
   # 
   #  ## List variables for the table
   #  Vars <- c("fldrug","fltimeon","flt2dmduration","flprehba1cmmol","fliniresp_final","flresid","fladherence_t","timesince1stline")
   #  FactorVars <- c("fldrug")
   #  ## BY MATCHED / ID Create the table
   #  Table <- CreateTableOne(vars=Vars, strata=c("drug"),factorVars = FactorVars, data=finalall)
   #  ## Export as an html document
   #  stargazer(print(Table, dropEqual = T))
  
####################
#Log outcome variable
####################
  # final.dev$posthba1c_finallog <- log(final.dev$posthba1c_final)
  # describe(final.dev$posthba1c_finallog)

####################
#Missing data
####################
#   vars <- final.dev %>% select("drugline","druglinena","ncurrtx","yrdrugstart","agedx","agetx","gender","malesex","dpp4subtype","sglt2subtype","t2dmduration","ethnicity5","adherence_t","hba1cmonth", "nonadh")
#   describe(vars)
#   vars <- final.dev %>% select("precreat","egfr","egfr_ckdepi","egfr_cg",
#                                 "preldl","prehdl","pretrig","preweight","prebmi",
#                                 "preast","prealt","preplatelets","prebil","prealb","api")  
#   describe(vars)
#   vars <- final.dev %>% select("prehba1cmmol","pregluc","prehba1cchange","iniresp_final","posthba1c_final","hba1cmmolm_pre1","hba1ctimediff_pre1","hba1cmmolm_pre2","hba1ctimediff_pre2")
#   describe(vars)
#   vars <- final.dev %>% select("fldrug","flt2dmduration","flprehba1cmmol","fliniresp_final","flresid","fladherence_t","timesince1stline","fldrugstopdate","flstoptime","flstop","flposthba1c_final","flnonadh")
#   describe(vars)
# 
# #################### 
# #Correlated data
# ####################
#   test <- final.dev %>% select("agedx","agetx","precreat","egfr_ckdepi",
#                                "preldl","prehdl","pretrig","preweight","prebmi",
#                                "preast","prealt","preplatelets","prebil","prealb","api","t2dmduration","adherence_t","trghdlratio",
#                                "prehba1cmmol","pregluc","prehba1cchange","fladherence_t","fliniresp_final","flt2dmduration","flprehba1cmmol","timesince1stline","flposthba1c_final")
#   descrCor <- cor(test, use="na.or.complete")
#   print(descrCor)
#   corrplot(descrCor, order = "FPC", method = "color", type = "lower", tl.cex = 0.7, tl.col = rgb(0, 0, 0))
#   highlyCorrelated <- findCorrelation(descrCor, cutoff=0.5)
#   highlyCorCol <- colnames(test)[highlyCorrelated]
#   highlyCorCol
  
####################
#Add additional risk factors
####################
    
###########################      
# Add comorbidity
###########################      
load("C:/Users/jmd237/OneDrive - University of Exeter/John/CPRD/CPRD_19_Comorbidities/cprd.comorbidity.Rda") 

#Merge in comorbidity
cprd.comorbidity$patid <- NULL
final.all <- merge(final.all,cprd.comorbidity,by="pateddrug",all.x=TRUE)

head(final.all)

# Established ASCVD: MI, ischemic stroke, unstable angina with ECG changes, myocardial ischemia on imaging or stress test, revascularisation of coronary, carotid or peripheral arteries.
# High risk ASCVD: Age 55 + Left ventricular hypertrophy or coronary, carotid, lower extremity artery stenosis >50%
# HF especially HFrEF (EF <45%)
# CKD - eGFR 30-60 or UACR>30 (preferably >300) - says 300 in main flow.

final.all <- final.all %>% mutate(ascvd=if_else(predrug.earliest.mi==1|predrug.earliest.stroke==1|predrug.earliest.revasc==1|predrug.earliest.pad==1|predrug.earliest.ihd==1,1,0))
table(final.all$ascvd);prop.table(table(final.all$ascvd))
describe(final.all$ascvd)

final.all <- final.all %>% mutate(ascvd.hf=if_else(ascvd==1|predrug.earliest.heartfailure==1,1,0)) 
table(final.all$ascvd.hf);prop.table(table(final.all$ascvd.hf))
describe(final.all$ascvd.hf)

final.all <- final.all %>% mutate(ckd=if_else(predrug.earliest.ckd5==1|egfr_ckdepi<60,1,0)) 
final.all <- final.all %>% mutate(ckd=if_else(is.na(ckd),0,ckd))
table(final.all$ckd);prop.table(table(final.all$ckd))
describe(final.all$ckd)

final.all <- final.all %>% mutate(ascvd.hf.ckd=if_else(ascvd.hf==1|ckd==1,1,0)) 
table(final.all$ascvd.hf.ckd);prop.table(table(final.all$ascvd.hf.ckd))
describe(final.all$ascvd.hf.ckd)

final.all <- final.all %>% mutate(hypertension=if_else(predrug.earliest.hypertension==1,1,0)) 

final.all <- final.all %>% mutate(neuropathy=if_else(predrug.earliest.neuropathy==1,1,0)) 
final.all <- final.all %>% mutate(nephropathy=if_else(predrug.earliest.nephropathy==1,1,0)) 
final.all <- final.all %>% mutate(retinopathy=if_else(predrug.earliest.retinopathy==1,1,0)) 
final.all <- final.all %>% mutate(hf=if_else(predrug.earliest.heartfailure==1,1,0)) 

##ADD SBP

load("C:/Users/jmd237/OneDrive - University of Exeter/John/local/CPRD/cprd/CPRD_19/data/update/finalcprd19_analysisreadyv3.Rda") 
cprdv3 <- cprd

#Smoking binary
describe(final.all$Category)
final.all <- final.all %>% mutate(smok=as.factor(ifelse(is.na(Category),"Non-smoker",Category)))
describe(final.all$smok)

#chol
cprd.tc.sys <- cprdv3 %>% dplyr::select(pateddrug,pretc,presys)
final.all <- merge(final.all,cprd.tc.sys,all.x=T, by="pateddrug")
describe(final.all$pretc)
final.all$chol <- final.all$pretc

#sbp
describe(final.all$presys)
final.all$sbp <- final.all$presys

#################### 
#save
####################
  # save(final.dev,file="cprd_19_sglt2dpp4_devcohort.Rda")
  # #write.csv(final.dev, file = "cprd_19_sglt2dpp4_devcohort.csv")
  # write.dta(final.dev, "cprd_19_sglt2dpp4_devcohort.dta") 
  # 
  # save(final.val,file="cprd_19_sglt2dpp4_valcohort.Rda")
  # # write.csv(final.val, file = "cprd_19_sglt2dpp4_valcohort.csv")
  # # write.dta(final.val, "cprd_19_sglt2dpp4_valcohort.dta") 
  
  save(final.all,file="cprd_19_sglt2dpp4_allcohort23.Rda")
  # write.csv(finalall, file = "cprd_19_sglt2dpp4_allcohort.csv")
  # write.dta(finalall, "cprd_19_sglt2dpp4_allcohort.dta") 