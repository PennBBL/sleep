library(ggplot2)
library(dplyr)
library(epiDisplay)
library(lme4)
library(reshape2)
library(psych)
library(table1)
library(visreg)
library(reshape2)
library(stringr)
library("lmerTest")

#cd /Users/kruparel/Library/Mobile Documents/com~apple~CloudDocs/BBL/Raquel_Gur/22q_sleep
# read merged data
datas<-read.csv("sleep_raw_scales_demo_cnb_cnbz_dx_sips_clinical_denovo_merged_feb2024.csv")#99

#remove 1 person who does not have a bblid and was not part of bbl
datas<-datas[datas$sleepid!=67,] #99

datas<-datas[which(is.na(datas$psqitot)==FALSE),]#92


#define good and bad sleepers based on psqi known cutoff
datas$psqitot<-as.numeric(datas$psqitot)
datas$sleeper<-NA
datas$sleeper[which(datas$psqitot>5)]<-"Poor Sleeper"
datas$sleeper[datas$psqitot<=5]<-"Good Sleeper"
datas$sleeper<-as.factor(datas$sleeper)
datas$X22q_inherited_denovo<-factor(datas$X22q_inherited_denovo)



#demographics
datas$gender<-NA
datas$gender[datas$sex.x==1]<-"Male"
datas$gender[datas$sex.x==2]<-"Female"
datas$gender<-as.factor(datas$gender)

datas$race[datas$race==1]<-"White"
datas$race[datas$race==2]<-"Black"
datas$race[datas$race==5]<-"Other"
datas$race<-as.factor(datas$race)

#demographics table
table1(~ age + gender+race| sleeper, data=datas)


datas$diagnosis<-NA
datas$diagnosis[datas$dx_pscat=="pro"]<-"PS"
datas$diagnosis[datas$dx_pscat=="psy"]<-"PS"
datas$diagnosis[datas$dx_pscat=="noDSMdx"]<-"TD"
datas$diagnosis[datas$dx_pscat=="other"]<-"Other"


datas$diagnosis<-as.factor(datas$diagnosis)

#table of comorbidity
datas$dx_prodromal<-as.factor(datas$dx_prodromal)
datas$dx_psychosis<-as.factor(datas$dx_psychosis)
datas$dx_scz<-as.factor(datas$dx_scz)
datas$dx_moodnos<-as.factor(datas$dx_moodnos)
datas$dx_mdd<-as.factor(datas$dx_mdd)
datas$dx_bp1<-as.factor(datas$dx_bp1)
datas$dx_adhd<-as.factor(datas$dx_adhd)
datas$dx_anx<-as.factor(datas$dx_anx)
datas$dx_ptsd<-as.factor(datas$dx_ptsd)
datas$dx_other<-as.factor(datas$dx_other)

table1(~ dx_prodromal+dx_psychosis+dx_scz+dx_moodnos+dx_mdd+dx_bp1+dx_adhd+dx_anx+dx_ptsd+dx_other| sleeper, data=datas,overall="Total",
    render.missing=NULL, render.categorical="FREQ (PCTnoNA%)")

#chisquare - diagnosis
chisq.test(table(datas$dx_prodromal,datas$sleeper))

chisq.test(table(datas$dx_psychosis,datas$sleeper))

chisq.test(table(datas$dx_scz,datas$sleeper))
chisq.test(table(datas$dx_moodnos,datas$sleeper))
chisq.test(table(datas$dx_mdd,datas$sleeper))


datas$dxage<-round(abs(as.Date(as.character(datas$dobirth),"%m/%d/%y")-as.Date(as.character(datas$date_diagnosis),"%m/%d/%y"))/365.5)
datas$dxage<-as.integer(datas$dxage)

# summary sips scores
#note those that do not have completed sips will be NA when counting rowSUMS
datas$sum_possips<-rowSums(datas[,c("p1","p2","p3","p4","p5")])
datas$sum_possips2<-rowSums(datas[,c("p1","p2","p3","p4","p5")])
datas$sum_negsips<-rowSums(datas[,c("n1","n2","n3","n4","n5","n6")])
datas$sum_disorgsips<-rowSums(datas[,c("d1","d2","d3","d4")])
datas$sum_gensips<-rowSums(datas[,c("g1","g2","g3","g4")])


#sips table
table1(~ gaf_c+diagnosis+sum_possips+sum_negsips+sum_disorgsips+sum_gensips| sleeper, data=datas,overall="Total",
    render.missing=NULL, render.categorical="FREQ (PCTnoNA%)",extra.col=list(`pval`=pvalue), extra.col.pos=4)

#with missing
table1(~ gaf_c+diagnosis+sum_possips+sum_negsips+sum_disorgsips+sum_gensips| sleeper, data=datas,overall="Total",
    render.categorical="FREQ (PCTnoNA%)",extra.col=list(`pval`=pvalue), extra.col.pos=4)


# t tests two sample- sips
datas %>% t_test(gaf_c ~ sleeper)
datas %>% t_test(sum_possips ~ sleeper)
datas %>% t_test(sum_negsips ~ sleeper)

datas %>% t_test(sum_disorgsips ~ sleeper)
datas %>% t_test(sum_gensips ~ sleeper)




#SIPS models

lm.psips<-lm(sum_possips~sleeper+dxage+gender,data=datas)
    #visreg(lm.psips)
    visreg(lm.psips,"dxage",by="sleeper")

lm.nsips<-lm(sum_negsips~sleeper+dxage+gender,data=datas)
    #visreg(lm.nsips)
    visreg(lm.nsips,"dxage",by="sleeper")

#continuous psqitot

lm2.psips<-lm(sum_possips~+psqitot+dxage+gender,data=datas)
visreg(lm2.psips)

lm2.nsips<-lm(sum_negsips~+psqitot+dxage+gender,data=datas)
visreg(lm2.nsips)

#GAF
lm.gaf<-lm(gaf_c~sleeper+dxage+gender,data=datas)
        visreg(lm.gaf)
lm2.gaf<-lm(gaf_c~+psqitot+dxage+sex,data=datas)
        visreg(lm2.gaf)

#library(tidyr)
#datas_cli2_long2<-datas_cli2_long %>% drop_na(sleep_condition)
#ggplot(datas_cli2_long2, aes(x=sleep_condition, y=psqitot)) + geom_boxplot()+ theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))

#ggplot(data=datas, mapping = aes(x = sleep_problem, y = psqitot))+ geom_boxplot() +
          #theme_bw()


##corrlation plot
#pairs(~ psqitot + soccog_eff + exe_eff+cpxres_eff+memory_eff+delib_speed, data = datas)


##Cognition
ggplot(datas, aes(x=psqitot, y=memory_eff)) + geom_point() + ggtitle("PSQI Total Score vs Memory Efficiency") + geom_smooth(method=lm, se=FALSE)

ggplot(datas, aes(x=psqitot, y=exe_eff)) + geom_point() + ggtitle("PSQI Total Score vs Memory Efficiency") + geom_smooth(method=lm, se=FALSE)

lm1<-lm(exe_eff~+sleeper+test_sessions_v.age+test_sessions_v.gender,data=datas)
lm2<-lm(soccog_eff~+sleeper+test_sessions_v.age+test_sessions_v.gender,data=datas)
lm3<-lm(cpxres_eff~+sleeper+test_sessions_v.age+test_sessions_v.gender,data=datas)
lm4<-lm(memory_speed~+sleeper+test_sessions_v.age+test_sessions_v.gender,data=datas)
lm5<-lm(delib_speed~+sleeper+test_sessions_v.age+test_sessions_v.gender,data=datas)
lm6<-lm(rapid_speed~+sleeper+test_sessions_v.age+test_sessions_v.gender,data=datas)
lm7<-lm(overall_speed~+sleeper+test_sessions_v.age+test_sessions_v.gender,data=datas)

lm8<-lm(gaf_c~+psqitot+test_sessions_v.age+test_sessions_v.gender,data=datas)


datas$age<-round(abs(as.Date(as.character(datas$dobirth),"%m/%d/%y")-as.Date(as.character(datas$date_of_entry),"%m/%d/%y"))/365.5)
datas$age<-as.integer(datas$age)



#tableStack(vars=c(gaf_c,dx_pscat,sum_possips,sum_negsips,sum_disorgsips,sum_gensips), by=sleeper, dataFrame=datas, name.test=FALSE,na.rm=F,total=T,iqr=c(gaf_c,sum_possips,sum_negsips,sum_disorgsips,sum_gensips))




#cognitive

table1(~ ABF_Efficiency+ATT_Efficiency+WM_Efficiency+VMEM_Efficiency+FMEM_Efficiency+SMEM_Efficiency| sleeper, data=datas)

table1(~ LAN_Efficiency+NVR_Efficiency+SPA_Efficiency+EID_Efficiency+EDI_Efficiency+ADI_Efficiency| sleeper, data=datas)


table1(~ memory_acc+soccog_acc+cpxres_acc+memory_eff+soccog_eff+exe_eff+cpxres_eff| sleeper, data=datas)



table1(~ overall_speed+memory_speed+delib_speed+rapid_speed+motor_speed| sleeper, data=datas)



cnb<-datas[,c(1,417,351:376,377:388)]
cnb_long<-melt(cnb, id.vars=c("bblid", "sleeper"))
names(cnb_long)[3:4]<-c("CNB_Test","Score")

source("~/summarySE.R")
smry_cnb<-summarySE(cnb_long, measurevar="Score", groupvars=c("sleeper","CNB_Test"),na.rm=T)

pd <- position_dodge(0.1) # move them .05 to the left and right

smry_cnb$sleeper = factor(smry_cnb$sleeper, levels = c("Poor Sleeper","Good Sleeper"))

smry_cnb_a<-smry_cnb[smry_cnb$CNB_Test

smry_cnb_a<-smry_cnb[grep("_A", smry_cnb$CNB_Test),]
smry_cnb_s<-smry_cnb[grep("_S", smry_cnb$CNB_Test),]
smry_cnb_e<-smry_cnb[grep("_Efficiency", smry_cnb$CNB_Test),]

ggplot(smry_cnb_e, aes(x=CNB_Test, y=Score, colour=sleeper, group=sleeper)) +
    geom_errorbar(aes(ymin=Score-se, ymax=Score+se), colour="black", width=.1, position=pd) +
    geom_line(position=pd) +
    geom_point(position=pd, size=3, shape=21, fill="white") + # 21 is filled circle
    xlab("CNB Tests") +
    ylab("Mean Efficiency Score") +
  scale_x_discrete(guide = guide_axis(angle = 45))


  ggplot(smry_cnb_a, aes(x=CNB_Test, y=Score, colour=sleeper, group=sleeper)) +
      geom_errorbar(aes(ymin=Score-se, ymax=Score+se), colour="black", width=.1, position=pd) +
      geom_line(position=pd) +
      geom_point(position=pd, size=3, shape=21, fill="white") + # 21 is filled circle
      xlab("CNB Tests") +
      ylab("Mean Accuracy Score") +
    scale_x_discrete(guide = guide_axis(angle = 45))


      ggplot(smry_cnb_s, aes(x=CNB_Test, y=Score, colour=sleeper, group=sleeper)) +
          geom_errorbar(aes(ymin=Score-se, ymax=Score+se), colour="black", width=.1, position=pd) +
          geom_line(position=pd) +
          geom_point(position=pd, size=3, shape=21, fill="white") + # 21 is filled circle
          xlab("CNB Tests") +
          ylab("Mean Speed Score") +
        scale_x_discrete(guide = guide_axis(angle = 45))

# ttests added 8/27/2024 on Ruben's request
# Load required R packages
library(tidyverse)
library(rstatix)
library(ggpubr)
stat_test <- cnb_long %>%
  group_by(CNB_Test) %>%
  t_test(Score ~ sleeper, p.adjust.method = "bonferroni")
# Remove unnecessary columns and display the outputs
stat_test %>% select(-.y., -statistic, -df)



##lme model
#eeficieny mixed effect not significant

lmer.cnb<-lmer(Eff_Score ~ CNB_Test * sleeper + (1 | bblid), data=cnb_long)
anova(lmer.cnb)

#accuracy model
cnba<-datas[,c(1,390,351:362)]
cnba_long<-melt(cnba, id.vars=c("bblid", "sleeper"))
names(cnba_long)[3:4]<-c("CNB_Test","Acc_Score")
lmer.cnba<-lmer(Acc_Score ~ CNB_Test * sleeper + (1 | bblid), data=cnba_long)
anova(lmer.cnba)
visreg(lmer.cnba,"CNB_Test",by="sleeper")

visreg(lmer.cnba)
#speed model
cnbs<-datas[,c(1,390,363:376)]
cnbs_long<-melt(cnbs, id.vars=c("bblid", "sleeper"))
names(cnbs_long)[3:4]<-c("CNB_Test","Speed_Score")
lmer.cnbs<-lmer(Speed_Score ~ CNB_Test * sleeper + (1 | bblid), data=cnbs_long)
visreg(lmer.cnbs,"sleeper",by="CNB_Test")

## By ruben's request a mixed model with speed and accuracy

cnb_long2<-cnb_long[-grep("_Efficiency", cnb_long$CNB_Test),]
cnb_long2[c('test', 'measure')] <- str_split_fixed(cnb_long2$CNB_Test, '_', 2)
#remove mot and sm
cnb_long3<-cnb_long2[-which(cnb_long2$test=="MOT" | cnb_long2$test=="SM"),]

#big mixed model
cnb_long3$test<-as.factor(cnb_long3$test)
cnb_long3$measure<-as.factor(cnb_long3$measure)
lmer.cnbas<-lmer(Score ~ measure*test * sleeper + (1 | bblid), data=cnb_long3)

anova(lmer.cnbas)
visreg(lmer.cnbas,"test",by="sleeper")
visreg(lmer.cnbas, 'test', by='sleeper', cond=list(measure="A"), layout=c(2,1))
visreg(lmer.cnbas, 'test', by='sleeper', cond=list(measure="S"), layout=c(2,2))

visreg(lmer.cnbas,"measure",by="sleeper")

##clinical data comparision from M. Sounders
## correlation matrix to compare psqitot to clinical sleep data from  M. Sounders
xcor <- mixedCor(datas_cli)$rho

ggcorrplot(xcor,
  type = "lower",
  insig = "blank",
  lab = TRUE,
  digits = 1,pch = 0.8, pch.col = "black", pch.cex =0.9,
           tl.cex = 8
)


datas_cli2<-datas[,c(2,62,417,394:406,409:410,412:414)]
datas_cli2_long <- melt(datas_cli2,
        # ID variables - all the variables to keep but not split apart on
    id.vars=c("sleepid","psqitot","sleeper"),
        # The source columns
    measure.vars=names(datas_cli2)[4:21],
    variable.name="condition",
    value.name="measurement"
)


#custom  functions

    pvalue <- function(x, ...) {
      x <- x[-length(x)]  # Remove "overall" group
      # Construct vectors of data y, and groups (strata) g
      y <- unlist(x)
      g <- factor(rep(1:length(x), times=sapply(x, length)))
      if (is.numeric(y)) {
        # For numeric variables, perform an ANOVA
        p <- summary(aov(y ~ g))[[1]][["Pr(>F)"]][1]
      } else {
        # For categorical variables, perform a chi-squared test of independence
        p <- chisq.test(table(y, g))$p.value
      }
      # Format the p-value, using an HTML entity for the less-than sign.
      # The initial empty string places the output on the line below the variable label.
      c("", sub("<", "&lt;", format.pval(p, digits=3, eps=0.001)))
    }

#not needed as we do not want rowsums if there are any NA's
    my_rowSums <- function(x) {
      if (is.data.frame(x)) x <- as.matrix(x)
      z <- base::rowSums(x, na.rm = TRUE)
      z[!base::rowSums(!is.na(x))] <- NA
      z
      }

#medications from Oracle
meds<-read.csv("../../DatabaseStuff/psycha1_data_repository/repository_sync/medicineall_v.csv")
meds2<-meds[meds$BBLID %in% datas$bblid,]

meds2$DOCOLLECT<-as.Date(meds2$DOCOLLECT, "%d-%B-%y")
meds2$yearcollect<-format(as.Date(meds2$DOCOLLECT),"%Y")

meds2<-meds2[which(meds2$yearcollect %in% c("2011","2012")),]


##denovo and deletion type
deletion_type<-read.csv("../22q_familial/DeletionType/22q_deletion_type_updated_rarecnv.csv")
datas_del<-merge(datas,deletion_type,by="bblid",all.x=T)
table(datas_del$deletion_type_from_Blaine)