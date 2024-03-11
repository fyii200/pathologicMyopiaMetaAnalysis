# Author: Fabian Yii
# Email: fabian.yii@ed.ac.uk

# For the purpose of open access, the author has applied a creative commons 
# attribution (CC BY) licence to any Author Accepted Manuscript version arising 
# from this work.

################################################################################
#                                 Setting up                                   #
################################################################################
rm(list=ls())
# install.packages("remotes")
# install.packages("devtools")
# library("devtools")
# remotes::install_github("wviechtb/metafor")
# devtools::install_github("NightingaleHealth/ggforestplot")
library(dplyr)
library(metafor)
library(meta)
library(ggforestplot)
library(forestplot)
require(gridExtra)
library(effectsize)
library(ggplot2)
library(forcats)
library(stringr)

dirName <- dirname(getwd())
dataPath <- paste0(dirName, "/cleanedData/Data.csv")
d <- read.csv(dataPath)


################################################################################
#                              Data preparation                                #
################################################################################
# Add a nominal ±0.01 to the 95% CI of OR associated with the prognostic factor
# "age" in Wong et al. (otherwise we would end up with undefined weight); resultant
# 95% CI should therefore be 0.99 to 1.01
d[which(d$adjusted==1 & d$onset==0 & d$study=="wong" & d$factor=="age"),]$OR_lwr <- 0.99
d[which(d$adjusted==1 & d$onset==0 & d$study=="wong" & d$factor=="age"),]$OR_upr <- 1.01

# PM onset: convert Indians vs Malays to Malays vs Indians (reference)
ind <- which(d$factor=="indian" & d$onset==1 & d$adjusted==1)
d[ind,]$OR <- 1/d[ind,]$OR
d[ind,]$OR_lwr <- 1/d[ind,]$OR_upr
d[ind,]$OR_upr <- 1/d[ind,]$OR_lwr
d[ind,]$factor <- "malay"

# PM progression: convert Indians vs Malays to Malays vs Indians (reference)
ind <- which(d$factor=="indian" & d$onset==0 & d$adjusted==1)
d[ind,]$OR <- 1/d[ind,]$OR
d[ind,]$OR_lwr <- 1/d[ind,]$OR_upr
d[ind,]$OR_upr <- 1/d[ind,]$OR_lwr
d[ind,]$factor <- "malay"

# Convert risk ratio to odds ratio, treating incidence rate or progression rate
# as baseline risk, whichever is appropriate.
# Wong et al. (onset); reported 6-year incidence rate = 1.2%
rows <- which(d$note=="rr" & d$adjusted==1 & d$onset==1 & d$study=="wong") 
d[rows,]$OR <- riskratio_to_oddsratio(d[rows,]$OR, 0.012)
d[rows,]$OR_lwr <- riskratio_to_oddsratio(d[rows,]$OR_lwr, 0.012)
d[rows,]$OR_upr <- riskratio_to_oddsratio(d[rows,]$OR_upr, 0.012)
# Wong et al. (progression); reported 6-year progression rate = 17%
rows <- which(d$note=="rr" & d$adjusted==1 & d$onset==0 & d$study=="wong") 
d[rows,]$OR <- riskratio_to_oddsratio(d[rows,]$OR, 0.17)
d[rows,]$OR_lwr <- riskratio_to_oddsratio(d[rows,]$OR_lwr, 0.17)
d[rows,]$OR_upr <- riskratio_to_oddsratio(d[rows,]$OR_upr, 0.17)
# Foo et al. (onset); reported 12-year incidence rate = 10.3%
rows <- which(d$note=="rr" & d$adjusted==1 & d$onset==1 & d$study=="foo") 
d[rows,]$OR <- riskratio_to_oddsratio(d[rows,]$OR, 0.103)
d[rows,]$OR_lwr <- riskratio_to_oddsratio(d[rows,]$OR_lwr, 0.103)
d[rows,]$OR_upr <- riskratio_to_oddsratio(d[rows,]$OR_upr, 0.103)
# Foo et al. (progression); reported 12-year progression rate = 12.0%
rows <- which(d$note=="rr" & d$adjusted==1 & d$onset==0 & d$study=="foo") 
d[rows,]$OR <- riskratio_to_oddsratio(d[rows,]$OR, 0.12)
d[rows,]$OR_lwr <- riskratio_to_oddsratio(d[rows,]$OR_lwr, 0.12)
d[rows,]$OR_upr <- riskratio_to_oddsratio(d[rows,]$OR_upr, 0.12)

# Convert the OR and 95% CI reported by each study into the 
# log OR (saved as "yi") and its variance (saved as "vi")
d <- conv.wald(out=OR, ci.lb=OR_lwr, ci.ub=OR_upr, data=d, n=n, transf=log)
# Capitalise the first letter of the name of each study
d$study <- str_to_title(d$study)
# Capitalise every letter of each factor name
d$factor <- toupper(d$factor)


################################################################################
# Visualising the adjusted OR (displayed on linear scale) associated with each #
#  explored factor using forest plots (with meta-analysis where appropriate)   #
################################################################################

######################## Explored risk factors: PM onset #######################            
png(file=paste0(dirName, "/plots/PMonsetForest.png"), width=8, height=7.5, units="in", res=1500)
d$weight <- "NA"
                                 ## Meta-analyses ##
##  Owing to a limited number of studies (<5) per factor, fixed-effects models  ##
##  are used (https://journals.sagepub.com/doi/full/10.1177/21925682221110527)  ##
# Risk factor 1: baseline age (Fang, Foo & Ueda)
ageOnsetRows <- which(d$factor=="AGE" & d$adjusted==1 & d$onset==1 & d$study!="Foo")
ageOnsetModel <- rma(yi=yi, vi=vi, data=d[ageOnsetRows,], measure="OR", method="EE")
d[ageOnsetRows,]$weight <- round(weights(ageOnsetModel), 0)
# Risk factor 2: baseline AL (Fang, Foo & Ueda)
ALonsetRows <- which(d$factor=="AL" & d$adjusted==1 & d$onset==1 & d$study!="Foo")
ALonsetModel <- rma(yi=yi, vi=vi, data=d[ALonsetRows,], measure="OR", method="EE")
d[ALonsetRows,]$weight <- round(weights(ALonsetModel), 0)
# Risk factor 3: Female (Fang, Foo & Ueda)
femaleOnsetRows <- which(d$factor=="FEMALE" & d$adjusted==1 & d$onset==1 & d$study!="Foo" & d$study!="Fang")
femaleOnsetModel <- rma(yi=yi, vi=vi, data=d[femaleOnsetRows,], measure="OR", method="EE")
d[femaleOnsetRows,]$weight <- round(weights(femaleOnsetModel), 0)
                                 ## Start plotting ##
## set font expansion factor and use a bold font
par(par(cex=0.5, font=2))
## Forest plot
dOnset <- subset(d, adjusted==1 & onset==1)
dOnset <- dOnset[!is.na(dOnset$OR),]
dOnset$factor <- factor(dOnset$factor, levels=c("G/D", "CATARACT", "AL CHANGE", "TESSELLATION", "MALAY", "CHINESE", "HIGHER EDU", "SER", "FEMALE", "AL", "AGE") )
dOnset <- dOnset[order(dOnset$factor), ]
mOnset <- rma(yi, vi, data=dOnset)
forest(mOnset, cex=0.5, psize=1, atransf=exp,
       at=log(c(0.25, 0.5, 1, 5, 10)),
       xlim=c(-4.5, 4.2), ylim=c(-3, 58),
       slab=study, ilab=cbind(n, fu_year, covariates, weight, p), 
       ilab.xpos=c(-4, -3.6, -3.3, 2.45, 2.82),
       ilab.pos=4,
       rows=c(1, 5, 9, 13, 17:19, 23:24, 28:29, 34:36, 41:44, 49:52),
       header=c("Study", "OR [95% CI]   " ),
       xlab="Adjusted Odds Ratio, OR (Log Scale)",
       addfit=FALSE)
## Switch to bold font
par(cex=0.5, font=2)
## Add additional column headings to the plot
text(c(-4.03, -3.62, -2.64, 2.39, 2.93), 57, c("Eyes", "FU", "& Covariates", "W (%)", "P"), pos=4)
text(-3.2, 57, "Variable", pos=4, col="maroon", font=2, cex=1.4)
text(-0.6, 58.5, "PM Onset", pos=4, cex=2)
## Add text for each factor
text(-3.3, c(53.5, 45.5, 37.5, 30.5, 25.5, 20.5, 14.5, 10.5, 6.5, 2.5), pos=4, col='maroon', cex=1.4,
     c("Baseline age", "Baseline AL", "Female", "Baseline SER", "Higher education level", "Ethnicity", "Baseline tessellated fundus", "Change in AL b/w baseline & FU", "Baseline cataract", "G/D"))
## Gridlines
abline(h=c(46.7, 38.7, 31.7, 27, 22, 16, 12, 8, 4, 0), lwd=0.1, lty=2)
## Footnote
text(-3.22, c(38,15), "*", col="maroon", cex=1.7)
text(-2.79, -2, "* p=0.52 (female) and p=0.62 (tessellation) in Fang but OR not reported", col="maroon")
text(-4.42, 19.5, "+", col="darkgreen", cex=1.2)
text(-3.67, -3.5, "+ Chinese vs Indians (reference)", col="darkgreen")
text(-4.42, c(17.5, 18.5), "^", col="darkblue", cex=1.1)
text(-3.7, -5, "^ Malays vs Indians (reference)", col="darkblue")
# Add diamond & text label for each pooled estimate
text(-4.04, c(47.6, 39.6, 32.6), "Pooled w/o Foo", col="red", font=2)
addpoly(ageOnsetModel, atransf=exp, row=47.5, cex=1.1, col="red", border="red", efac=0.4, font=2, mlab="")
text(3.05, 47.6, ifelse(ageOnsetModel$pval<0.001, "<0.001", fmtp(ageOnsetModel$pval)), cex=1.1, col="red", font=2)
addpoly(ALonsetModel, atransf=exp, row=39.5, cex=1.1, col="red", border="red", efac=0.4, fonts="bold", mlab="")
text(3.05, 39.6, ifelse(ALonsetModel$pval<0.001, "<0.001", fmtp(ALonsetModel$pval)), cex=1.1, col="red", font=2)
addpoly(femaleOnsetModel, atransf=exp, row=32.5, cex=1.1, col="red", border="red", efac=0.4, fonts="bold", mlab="")
text(3.05, 32.6, ifelse(femaleOnsetModel$pval<0.001, "<0.001", fmtp(femaleOnsetModel$pval)), cex=1.1, font=2)
## Save plot as png
dev.off()

################### Explored prognostic factors: PM progression ##################
png(file=paste0(dirName, "/plots/PMProgressionForest.png"), width=8, height=7.5, units="in", res=1500)
d$weight <- "NA"
                              ## Meta-analyses ##
##  Owing to a limited number of studies (<5) per factor, fixed-effects models  ##
##  are used (https://journals.sagepub.com/doi/full/10.1177/21925682221110527)  ##
# Prognostic factor 1: baseline age (Fang, Foo, Hopf & Lin)
ageProgRows <- which(d$factor=="AGE" & d$adjusted==1 & d$onset==0 & d$study!="Foo")
ageProgModel <- rma(yi=yi, vi=vi, data=d[ageProgRows,], measure="OR", method="EE")
d[ageProgRows,]$weight <- round(weights(ageProgModel), 0)
# Prognostic factor 2: baseline AL (Fang & Foo)
ALprogRows <- which(d$factor=="AL" & d$adjusted==1 & d$onset==0 & d$study!="Foo")
ALprogModel <- rma(yi=yi, vi=vi, data=d[ALprogRows,], measure="OR", method="EE")
d[ALprogRows,]$weight <- round(weights(ALprogModel), 0)
# Prognostic factor 3: baseline SER (Foo, Hopf & Lin)
SERprogRows <- which(d$factor=="SER" & d$adjusted==1 & d$onset==0 & d$study!="Foo")
SERprogModel <- rma(yi=yi, vi=vi, data=d[SERprogRows,], measure="OR", method="EE")
d[SERprogRows,]$weight <- round(weights(SERprogModel), 0)
# Prognostic factor 4: Female sex (Foo, Hopf & Lin)
femaleProgRows <- which(d$factor=="FEMALE" & d$adjusted==1 & d$onset==0 & d$study!="Foo")
femaleProgModel <- rma(yi=yi, vi=vi, data=d[femaleProgRows,], measure="OR", method="EE")
d[femaleProgRows,]$weight <- round(weights(femaleProgModel), 0)
# Prognostic factor 5: Higher education level (Foo & Lin)
eduProgRows <- which(d$factor=="HIGHER EDU" & d$adjusted==1 & d$onset==0 & d$study!="Foo")
eduProgModel <- rma(yi=yi, vi=vi, data=d[eduProgRows,], measure="OR", method="EE")
d[eduProgRows,]$weight <- round(weights(eduProgModel), 0)
                                   ## Start plotting ##
## set font expansion factor and use a bold font
par(par(cex=0.5, font=2))
## Forest plot
dProg <- subset(d, adjusted==1 & onset==0)
dProg <- dProg[!is.na(dProg$OR),]
dProg$factor <- factor(dProg$factor, levels=c("G/D", "MM 3 OR 4", "HTN", "SER CHANGE", "AL CHANGE", "IOP", "MALAY", "CHINESE", "HIGHER EDU", "FEMALE", "SER", "AL", "AGE") )
dProg <- dProg[order(dProg$factor), ]
mProg <- rma(yi, vi, data=dProg)
forest(mProg, cex=0.5, psize=1, atransf=exp,
       at=log(c(0.25, 0.5, 1, 5, 10)),
       xlim=c(-4.5, 4.2), ylim=c(-1,72),
       slab=study, ilab=cbind(n, fu_year, covariates, weight, p), 
       ilab.xpos=c(-4, -3.6, -3.3, 2.45, 2.83),
       ilab.pos=4,
       rows=c(1, 5, 9, 13, 17, 21, 25:26, 31:33, 38:42, 47:50, 55:57, 62:66),
       header=c("Study", "OR [95% CI]   " ),
       xlab="Adjusted Odds Ratio, OR (Log Scale)",
       addfit=FALSE)
## Switch to bold font
par(cex=0.5, font=2)
## Add additional column headings to the plot
text(c(-4.03,-3.62, -2.66, 2.37, 2.9), 71, c("Eyes", "FU", "& Covariates", "W (%)", "P"), pos=4)
text(-3.2, 71, "Variable", pos=4, col="maroon", cex=1.4)
text(-0.8, 73, "PM Progression", pos=4, cex=2)
## Add text for each factor
text(-3.3, c(67.5, 58.5, 51.5, 43.3, 34.5, 27.3, 22.4, 18.5, 14.5, 10.5, 6.4, 2.3), pos=4, col='maroon', cex=1.4,
     c("Baseline age", "Baseline AL", "Baseline SER", "Female", "Higher education level", "Ethnicity", "Baseline IOP", "Change in AL b/w baseline & FU", "Change in SER b/w baseline & FU", "Baseline hypertension", "Baseline patchy/macular atrophy", "G/D"))
## Gridlines
abline(h=c(59.5, 52.5, 44.5, 35.5, 28.5, 24, 20, 16, 12, 8, 4, -0.2), lwd=0.1, lty=2)
## Footnote
text(-4.42, 26.5, "*", col="darkgreen", cex=1.7)
text(-3.55, -2.2, "* Chinese vs Malays + Indians (reference)", col="darkgreen")
text(-4.42, 25.4, "+", col="darkblue", cex=1.3)
text(-3.78, -3.5, "+ Malays vs Indians (reference)", col="darkblue")
# Add diamond & text label for each pooled estimate
text(-4.06, c(60.6, 53.6, 45.6, 36.6, 29.6), "Pooled w/o Foo", col="red", font=2)
addpoly(ageProgModel, atransf=exp, row=60.5, cex=1.1, col="red", border="red", efac=0.4, font=2, mlab="")
text(3.05, 60.65, ifelse(ageProgModel$pval<0.001, "<0.001", fmtp(ageProgModel$pval)), cex=1.1, col="red", font=2)
addpoly(ALprogModel, atransf=exp, row=53.5, cex=1.1, col="red", border="red", efac=0.4, font=2, mlab="")
text(3.05, 53.6, ifelse(ALprogModel$pval<0.001, "<0.001", fmtp(ALprogModel$pval)), cex=1.1, col="red", font=2)
addpoly(SERprogModel, atransf=exp, row=45.5, cex=1.1, col="red", border="red", efac=0.4, font=2, mlab="")
text(3.05, 45.6, ifelse(SERprogModel$pval<0.001, "<0.001", fmtp(SERprogModel$pval)), cex=1.1, col="red", font=2)
addpoly(femaleProgModel, atransf=exp, row=36.5, cex=1.1, col="red", border="red", efac=0.4, font=2, mlab="")
text(3.05, 36.6, ifelse(femaleProgModel$pval<0.001, "<0.001", fmtp(femaleProgModel$pval)), cex=1.1, col="red", font=2)
addpoly(eduProgModel, atransf=exp, row=29.5, cex=1.1, col="red", border="red", efac=0.4, font=2, mlab="")
text(3.05, 29.6, ifelse(eduProgModel$pval<0.001, "<0.001", fmtp(eduProgModel$pval)), cex=1.1, col="red", font=2)
## Save plot as png
dev.off()

################################################################################
#                             Sensitivity analyses                             #
################################################################################

## Part 1: substitute Wong (6-year FU; SEED cohort) with Foo (12-year FU; SEED cohort)
# PM onset
for(factorName in c("AGE", "AL", "FEMALE")){
  subData <- subset(d, factor==factorName & adjusted==1 & onset==1 & study!="Wong")
  subDataModel <- rma(yi=yi, vi=vi, data=subData, measure="OR", method="EE")
  print(paste("===================", factorName, "==================="))
  print(predict(subDataModel, transf=exp)); print(paste("p =",subDataModel$pval)) }
# PM progression
for(factorName in c("AGE", "AL", "SER", "FEMALE", "HIGHER EDU")){
  subData <- subset(d, factor==factorName & adjusted==1 & onset==0 & study!="Wong")
  subDataModel <- rma(yi=yi, vi=vi, data=subData, measure="OR", method="EE")
  print(paste("===================", factorName, "==================="))
  print(predict(subDataModel, transf=exp)); print(paste("p =",subDataModel$pval)) }

## Part 2: Bias age OR reported by Wong towards lower risk and higher risk 
# Comment out lines 266-272 for increased risk 
# Comment out lines 258-264 for the decrased risk
d <- read.csv(dataPath)

# Bias the estimate towards increased risk with older baseline age
d[which(d$adjusted==1 & d$onset==0 & d$study=="wong" & d$factor=="age"),]$OR <-
  riskratio_to_oddsratio(1.03, 0.17)
d[which(d$adjusted==1 & d$onset==0 & d$study=="wong" & d$factor=="age"),]$OR_lwr <-
  riskratio_to_oddsratio(1.00, 0.17)
d[which(d$adjusted==1 & d$onset==0 & d$study=="wong" & d$factor=="age"),]$OR_upr <-
  riskratio_to_oddsratio(1.04, 0.17)

# Bias the estimate towards decreased risk with older baseline age
d[which(d$adjusted==1 & d$onset==0 & d$study=="wong" & d$factor=="age"),]$OR <-
  riskratio_to_oddsratio(0.97, 0.17)
d[which(d$adjusted==1 & d$onset==0 & d$study=="wong" & d$factor=="age"),]$OR_lwr <-
  riskratio_to_oddsratio(0.96, 0.17)
d[which(d$adjusted==1 & d$onset==0 & d$study=="wong" & d$factor=="age"),]$OR_upr <-
  riskratio_to_oddsratio(1.00, 0.17)

# Convert risk ratio to odds ratio, treating progression rate as baseline risk
# Wong et al. (progression); reported 6-year progression rate = 17%
rows <- which(d$note=="rr" & d$adjusted==1 & d$onset==0 & d$study=="wong") 
d[rows,]$OR <- riskratio_to_oddsratio(d[rows,]$OR, 0.17)
d[rows,]$OR_lwr <- riskratio_to_oddsratio(d[rows,]$OR_lwr, 0.17)
d[rows,]$OR_upr <- riskratio_to_oddsratio(d[rows,]$OR_upr, 0.17)
d <- conv.wald(out=OR, ci.lb=OR_lwr, ci.ub=OR_upr, data=d, n=n, transf=log )

# Meta-analysis
ageProg <- subset(d, factor=="age" & adjusted==1 & onset==0 & study!="foo")
mAge <- rma(yi=yi, vi=vi, data=ageProg, measure="OR", method="EE")
predict(mAge, transf=exp); mAge$pval







  

    




