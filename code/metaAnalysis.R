# Author: Fabian SL Yii
# Email: fabian.yii@ed.ac.uk

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
d <- conv.wald(out=OR, ci.lb=OR_lwr, ci.ub=OR_upr, data=d, n=n, transf=log )
# Capitalise the first letter of the name of each study
d$study <- str_to_title(d$study)
# Capitalise every letter of each factor name
d$factor <- toupper(d$factor)


################################################################################
# Visualising the adjusted OR (displayed on linear scale) associated with each
# explored factor using forest plots (without pooled estimates)
################################################################################
## Explored risk factors: PN onset ##
pdf(file=paste0(dirName, "/plots/PMOnsetAll.pdf"))
dOnset <- subset(d, adjusted==1 & onset==1)
dOnset <- dOnset[!is.na(dOnset$OR),]
mOnset <- rma(yi, vi, data=dOnset)
forest(mOnset, order=factor,  cex=0.5, psize=1, atransf=exp,
       at=log(c(0.25, 0.5, 1, 5)),
       xlim=c(-5,3.5), ylim=c(-1,37),
       slab=study, ilab=cbind(n, fu_year, covariates, p), 
       ilab.xpos=c(-4.4,-3.85, -3.5, 1.7),
       ilab.pos=4,
       rows=c(1:4, 6:9, 11, 13, 15, 17:19, 21, 23:24, 26, 28, 30:31, 33),
       header=c("Study", "OR [95% CI]"),
       xlab="Adjusted Odds Ratio (OR)",
       addfit=FALSE)
# Switch to bold font
par(cex=0.5, font=2)
# Add additional column headings to the plot
text(c(-4.45,-3.85, -3.1, 1.85), 36, 
     c("N eyes", "FU", "Adjusted for", "P"), pos=4)
# Add text for each factor
text(-5.12, c(4.8, 9.8, 11.8, 13.8, 15.8, 19.8, 21.8, 24.8, 26.8, 28.8, 31.8, 33.8), pos=4, 
     c("Age", "AL", "AL change", "Cataract", "Chinese", "Female", "G/D zones", "Higher edu", "Indian", "Malay", "SER", "Tessellation"),
     col="maroon")
dev.off()

## Explored prognostic factors: PN progression ##
pdf(file=paste0(dirName, "/plots/PMProgressionAll.pdf"))
# set font expansion factor (as in forest() above) and use a bold font
op <- par(cex=0.5, font=2)
# set "par" back to the original settings
par(op)
dProg <- subset(d, adjusted==1 & onset==0)
dProg <- dProg[!is.na(dProg$OR),]
mProg <- rma(yi, vi, data=dProg)
forest(mProg, order=factor,  cex=0.5, psize=1, atransf=exp,
       at=log(c(0.25, 0.5, 1, 5)),
       xlim=c(-5,3.5), ylim=c(-1,44),
       slab=study, ilab=cbind(n, fu_year, covariates, p), 
       ilab.xpos=c(-4.4,-3.85, -3.5, 1.9),
       ilab.pos=4,
       rows=c(1:5, 7:9, 11, 13, 15:19, 21, 23:25, 27, 29, 31, 33, 35:38, 40),
       header=c("Study", "OR [95% CI]"),
       xlab="Adjusted Odds Ratio (OR)",
       addfit=FALSE)
# Switch to bold font
par(cex=0.5, font=2)
# Add additional column headings to the plot
text(c(-4.45,-3.85, -3.1, 2), 43, 
     c("N eyes", "FU", "Adjusted for", "P"), pos=4)
# Add text for each factor
text(-5.1, c(5.8, 9.8, 11.8, 13.8, 19.8, 21.8, 25.8, 27.8, 29.8, 31.8, 33.8, 38.8, 40.8), pos=4, col='maroon',
     c("Age", "AL", "AL change", "Chinese", "Female", "G/D", "Higher edu", "HTN", "Indian", "IOP", "MM 3/4", "SER", "SER change"))
dev.off()


################################################################################
#                            Meta-analysis: PM onset                           #
################################################################################
# Owing to a limited number of studies (<5) per factor, fixed-effects models
# are used (https://journals.sagepub.com/doi/full/10.1177/21925682221110527)

## Subset of data containing only risk factors and studies that can be meta-analysed
metaData <- subset(d, adjusted==1 & onset==1 & factor!="TESSELATED FUNDUS" & factor!="CHINESE" & factor!="MALAY" & factor!="CATARACT" & factor!="HIGHER EDU" & factor!="SER" & fu_year!=12 & fu_year!=18)
metaData <- metaData %>% select(study, OR, OR_upr, OR_lwr, n, p, factor)

## Risk factor 1: baseline age (Ueda & Wong)
ageOnset <- subset(d, factor=="AGE" & adjusted==1 & onset==1)[3:4,]
mAge <- rma(yi=yi, vi=vi, data=ageOnset, measure="OR", method="EE")
summary(mAge); predict(mAge, transf=exp)
# Add pooled OR (95% CI), N and p-value associated with baseline age
metaData[nrow(metaData)+1,] <- c("Pooled", 1.0819, 1.1173, 1.0476, 5537, "<0.001", "AGE")

## Risk factor 2: baseline AL (Ueda & Wong)
ALOnset <- subset(d, factor=="AL" & adjusted==1 & onset==1)[3:4,]
mAL <- rma(yi=yi, vi=vi, data=ALOnset, measure="OR", method="EE") 
summary(mAL); predict(mAL, transf=exp)
# Add pooled OR (95% CI), N and p-value associated with baseline AL
metaData[nrow(metaData)+1,] <- c("Pooled", 2.2380, 2.7534, 1.8190, 5537, "<0.001", "AL")

## Risk factor 3: female (Ueda & Wong)
femaleOnset <- subset(d, factor=="FEMALE" & adjusted==1 & onset==1)[3:4,]
mFemale <- rma(yi=yi, vi=vi, data=femaleOnset, measure="OR", method="EE") 
summary(mFemale); predict(mFemale, transf=exp)
# Add pooled OR (95% CI), N and p-value associated with female sex
metaData[nrow(metaData)+1,] <- c("Pooled", 1.1176, 1.9601, 0.6373, 5537, "0.70", "FEMALE")

## Forest plot (with pooled estimate for each risk factor)
# Make sure numeric data are coded as such
metaData[,2:5] <- sapply(metaData[,2:5], as.numeric)
# 95% CI label
metaData$CILabel <- paste0("[", format(round(metaData$OR_lwr,2), nsmall=2), " to ", format(round(metaData$OR_upr,2), nsmall=2), "]")
# Specify label corresponding to the OR associated with each risk factor from each study 
factorLabelLevels <- c("AGE : Pooled", "AGE : Wong", "AGE : Ueda",
                       "AL : Pooled", "AL : Wong", "AL : Ueda",
                       "FEMALE : Pooled", "FEMALE : Wong", "FEMALE : Ueda")
metaData$factorLabel <- factor(paste(metaData$factor, ':', metaData$study), levels=factorLabelLevels)
studyLevels <- c("Wong", "Ueda", "Pooled") # make sure the pooled estimate appears last
metaData$study <- factor(metaData$study, levels=studyLevels )
metaData <- metaData %>% arrange(factorLabel)
# Specify point estimate shape: 10 (circled plus) and 22 (filled square) 
# correspond to pooled estimate and individual estimate, respectively
shapes <- c(10, 22, 22, 10, 22, 22, 10, 22, 22)
# Forest plot
p <- ggplot(metaData, aes(x=OR, y=factorLabel, xmin = OR_lwr, xmax = OR_upr, colour=study, fill=study)) +
  geom_pointrange(shape=shapes, size=0.8) +
  geom_vline(xintercept = 1, linetype = 3) +
  xlab("Adjusted OR (95% CI)") +
  theme_classic() +
  scale_x_log10(limits = c(0.25, 4), 
                breaks = c(0.25, 0.5, 1, 2, 4), 
                labels = c("0.25", "0.5", "1", "2", "4"), expand = c(0,0)) +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.line.y =element_blank(),
        axis.ticks.y=element_blank(),
        plot.margin = margin(30, 0, 20, 0),
        legend.position = "left",
        legend.title=element_blank(),
        title=element_text(size=10),
        strip.text.x = element_text(size = 12),
        strip.background = element_blank()) +
  guides(color = guide_legend(override.aes=list(shape = c(22,22,10)))) +
  facet_wrap(~factor, ncol=1, scales="free_y", strip.position = "top") 
# Table
partitions <- c('—————————————————————————', '—————————————————————————', '—————————————————————————')
names(partitions) <- c("AGE", "AL", "FEMALE")
dataTable <- ggplot(metaData, aes(y = factorLabel)) +
  xlim(c(0,3)) +
  geom_text(aes(x = 0.05, label = n), hjust = 0) +
  geom_text(aes(x = 1, label = round(OR,2) )) +
  geom_text(aes(x = 2.3, label = CILabel), hjust = 1) +
  geom_text(aes(x = 2.7, label = p)) +
  scale_colour_identity() +
  theme_void() + 
  labs(subtitle="  N eyes             OR [95% CI]            P") +
  theme(plot.margin = margin(20, 30, 40, 0)) +
  facet_wrap(~factor, ncol=1, scales="free_y", labeller=labeller(factor=partitions)) 
# Save combined plot (forest plot + table)
metaOnset <- arrangeGrob(p, dataTable, widths=c(2,1.3))
ggsave(file=paste0(dirName, "/plots/metaOnset.pdf"), metaOnset)


################################################################################
#                        Meta-analysis: PM progression                         #
################################################################################
# Owing to a limited number of studies (<5) per prognostic factor, fixed-effects 
# models are used (https://journals.sagepub.com/doi/full/10.1177/21925682221110527)

## Subset of data containing only prognostic factors and studies that can be meta-analysed
metaData <- subset(d, adjusted==1 & onset==0 & factor!="MM 3 OR 4" & factor!="CHINESE" & factor!="SER CHANGE" & factor!="HTN" & factor!="TESSELATED FUNDUS" & factor!="IOP" & factor!="HYPERTENSION" & factor!="AL" & fu_year!=12 & fu_year!=18)
metaData <- metaData %>% select(study, OR, OR_upr, OR_lwr, n, p, factor)

## Prognostic factor 1: baseline age (Hopf & Lin)
ageProg <- subset(d, factor=="AGE" & adjusted==1 & onset==0)[3:5,]
mAge <- rma(yi=yi, vi=vi, data=ageProg, measure="OR", method="EE")
summary(mAge); predict(mAge, transf=exp)
# Add pooled OR (95% CI), N and p-value associated with baseline age
metaData[nrow(metaData)+1,] <- c("Pooled", 0.9967, 1.0084, 0.9851, 373, "0.58", "AGE")

## Prognostic factor 2: female (Hopf, Lin & Wong)
femaleProg <- subset(d, factor=="FEMALE" & adjusted==1 & onset==0)[3:5,]
mFemale <- rma(yi=yi, vi=vi, data=femaleProg, measure="OR", method="EE") 
summary(mFemale); predict(mFemale, transf=exp)
# Add pooled OR (95% CI), N and p-value associated with female sex
metaData[nrow(metaData)+1,] <- c("Pooled", 2.2936, 4.5000, 1.1690, 373, "0.02", "FEMALE")

## Prognostic factor 3: baseline SER (Hopf, Lin & Wong)
SERprog <- subset(d, factor=="SER" & adjusted==1 & onset==0)[2:4,]
mSER <- rma(yi=yi, vi=vi, data=SERprog, measure="OR", method="EE") 
summary(mSER); predict(mSER, transf=exp)
# Add pooled OR (95% CI), N and p-value associated with baseline SER
metaData[nrow(metaData)+1,] <- c("Pooled", 0.8726, 0.9172, 0.8302, 373, "<0.001", "SER")

## Prognostic factor 4: higher education level (Lin & Wong)
eduProg <- subset(d, factor=="HIGHER EDU" & adjusted==1 & onset==0)[2:3,]
mEdu <- rma(yi=yi, vi=vi, data=eduProg, measure="OR", method="EE") 
summary(mEdu); predict(mEdu, transf=exp)
# Add pooled OR (95% CI), N and p-value associated with higher education level
metaData[nrow(metaData)+1,] <- c("Pooled", 3.1653, 7.3468, 1.3638, 339, "0.01", "HIGHER EDU")

## Forest plot (with pooled estimate for each prognostic factor)
# Make sure numeric data are coded as such
metaData[,2:5] <- sapply(metaData[,2:5], as.numeric)
# 95% CI label
metaData$CILabel <- paste0("[", format(round(metaData$OR_lwr,2), nsmall=2), " to ", format(round(metaData$OR_upr,2), nsmall=2), "]")
# Specify label corresponding to the OR associated with each prognostic factor from each study 
factorLabelLevels <- c("AGE : Pooled", "AGE : Hopf", "AGE : Lin", "AGE : Wong",
                       "FEMALE : Pooled", "FEMALE : Hopf","FEMALE : Lin", "FEMALE : Wong",
                       "SER : Pooled", "SER : Hopf", "SER : Lin", "SER : Wong", 
                       "HIGHER EDU : Pooled", "HIGHER EDU : Lin", "HIGHER EDU : Wong")
metaData$factorLabel <- factor(paste(metaData$factor, ":", metaData$study), levels=factorLabelLevels)
studyLevels <- c("Hopf", "Lin", "Wong", "Pooled")
metaData$study <- factor(metaData$study, levels=studyLevels )
metaData <- metaData %>% arrange(factorLabel)
# Clip upper 95% CI at 10
metaData[metaData$OR_upr>10,]$OR_upr <- 10
# Specify point estimate shape: 10 (circled plus) and 22 (filled square) 
# correspond to pooled estimate and individual estimate, respectively
shapes <- c(10,22,22,22,10,22,22,22,10,22,22,22,10,22,22)
# Forest plot
p <- ggplot(metaData, aes(x=OR, y=factorLabel, xmin = OR_lwr, xmax = OR_upr, colour=study, fill=study)) +
  geom_pointrange(shape=shapes, size=0.6) +
  geom_vline(xintercept = 1, linetype = 3) +
  xlab("Adjusted OR (95% CI)") +
  theme_classic() +
  scale_x_log10(limits = c(0.25, 10), 
                breaks = c(0.25, 0.5, 1, 2, 4, 8), 
                labels = c("0.25", "0.5", "1", "2", "4", "8"), expand = c(0,0)) +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.line.y =element_blank(),
        axis.ticks.y=element_blank(),
        plot.margin = margin(10, 0, 2, 0),
        legend.position = "left",
        legend.title=element_blank(),
        title=element_text(size=10),
        strip.text.x = element_text(size = 11),
        strip.background = element_blank()) +
  guides(color = guide_legend(override.aes=list(shape = c(22,22,22,10)))) +
  facet_wrap(~factor, ncol=1, scales="free_y", strip.position="top")
# Table 
partitions <- c("———————————————————————", "———————————————————————", "———————————————————————", "———————————————————————")
names(partitions) <- c("AGE", "FEMALE", "SER", "HIGHER EDU")
dataTable <- ggplot(metaData, aes(y = factorLabel)) +
  xlim(c(0,3)) +
  geom_text(aes(x = 0.1, label = n), hjust = 0, size=4) +
  geom_text(aes(x = 0.8, label = round(OR,2)), size=4) +
  geom_text(aes(x = 2.2, label = CILabel), hjust = 1, size=4) +
  geom_text(aes(x = 2.7, label = p), size=4) +
  scale_colour_identity() +
  theme_void() + 
  labs(subtitle="   N eyes            OR (95% CI)              P") +
  theme(plot.margin = margin(0, 0, 29, 0),
        plot.subtitle=element_text(size=10)) +
  facet_wrap(~factor, ncol=1, scales='free_y', labeller=labeller(factor=partitions)) 
# Save combined plot (forest plot + table)
metaProgression <- arrangeGrob(p, dataTable, widths = c(2,1.3))
ggsave(file=paste0(dirName, "/plots/metaProgression.pdf"), metaProgression)


################################################################################
#                             Sensitivity analyses                             #
################################################################################
# Comment out lines 357 to 363 for the first sensitivity analysis (increased risk
# with older baseline age). Comment out lines 349-355 for the second analysis.

d <- read.csv(dataPath)

# bias the estimate towards increased risk with older baseline age
d[which(d$adjusted==1 & d$onset==0 & d$study=="wong" & d$factor=="age"),]$OR <-
  riskratio_to_oddsratio(1.03, 0.17)
d[which(d$adjusted==1 & d$onset==0 & d$study=="wong" & d$factor=="age"),]$OR_lwr <-
  riskratio_to_oddsratio(1.00, 0.17)
d[which(d$adjusted==1 & d$onset==0 & d$study=="wong" & d$factor=="age"),]$OR_upr <-
  riskratio_to_oddsratio(1.04, 0.17)

# bias the estimate towards decreased risk with older baseline age
d[which(d$adjusted==1 & d$onset==0 & d$study=="wong" & d$factor=="age"),]$OR <-
  riskratio_to_oddsratio(0.97, 0.17)
d[which(d$adjusted==1 & d$onset==0 & d$study=="wong" & d$factor=="age"),]$OR_lwr <-
  riskratio_to_oddsratio(0.96, 0.17)
d[which(d$adjusted==1 & d$onset==0 & d$study=="wong" & d$factor=="age"),]$OR_upr <-
  riskratio_to_oddsratio(1.00, 0.17)

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

d <- conv.wald(out=OR, ci.lb=OR_lwr, ci.ub=OR_upr, data=d, n=n, transf=log )

ageProg <- subset(d, factor=="age" & adjusted==1 & onset==0)[3:5,]
mAge <- rma(yi=yi, vi=vi, data=ageProg, measure="OR", method="EE")
summary(mAge); predict(mAge, transf=exp)







  

    




