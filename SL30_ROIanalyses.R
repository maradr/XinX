#SL30 - ROI analyses
rm(list=ls())
setwd("D:/tegolino-derosa/_Experiments/SL30/ROIanalysis_M/")

library(plyr)
library(ez);# Anova
library(car);# Anova
library(lme4) # (g)lmer
library(multcomp) # glht
library(lattice) # dotplot
library(MASS) # boxcox
library(reshape2) # reshape data format
library(ggplot2)
library(CGPfunctions)
library(afex)
filenames <- list.files()
filenames <- filenames[(grep(".txt$", filenames))]

for(i in 1:length(filenames))
{
  if (i == 1)

  {dataAll <- read.table(filenames[i])
    dataAll$cond <- as.factor(substr(filenames[i], 1,2))
    dataAll$hemisph <- as.factor(toupper(substr(filenames[i], 4,4)))
    dataAll$sbj <- row.names(dataAll)}
  else
  {temp <- read.table(filenames[i])
  temp$cond <- as.factor(substr(filenames[i], 1,2))
  temp$hemisph <- as.factor(toupper(substr(filenames[i], 4,4)))
  temp$sbj <- row.names(temp)
  dataAll <- rbind(dataAll, temp)}
  
}
rm(temp, filenames, i)

colnames(dataAll)[1] <- c("amplitude")





dataAll$combo <- as.factor(paste0(as.character(dataAll$cond), as.character(dataAll$hemisph)))


dataAll <- read.table("SL30_summedBC.txt", header=T)
summary(dataAll)
dataAll$sbj <- as.factor(dataAll$sbj)


ggplot(dataAll, aes(y = amplitude, x = cond, group=combo))+ ylab("Amplitude (mV)") +
  geom_boxplot(aes(fill = hemisph), alpha = .5, outlier.shape = NA) +
  geom_point(aes(col=hemisph), size=2, position=position_jitterdodge(.7)) +
  theme_bw()

dataAll$sbj <- as.factor(dataAll$sbj)
qqnorm(dataAll$amplitude)
qqline(dataAll$amplitude)

hist(dataAll$amplitude)

m1 <- lmer(amplitude ~ cond*hemisph  + (1|sbj), dataAll)
summary(m1)

Anova(m1)

m2 <- lmer(amplitude ~ cond+hemisph  + (1|sbj), dataAll)

anova(m1, m2)


ez1 <- ezANOVA(dataAll, amplitude, sbj, within = .(cond, hemisph), detailed  = T,type=3, return_aov = T)

ez1$ANOVA

x <- Plot2WayANOVA(formula = amplitude ~ hemisph*cond, dataframe = dataAll)

plotting <- data.frame(x$MeansTable)

trad <- aov(amplitude ~ cond*hemisph, data=dataAll)
summary(trad)


#model<- afex::aov_ez(id ="sbj", dv="amplitude", within = c("cond", "hemisph"), data=dataAll)



library(BayesFactor)
bf <- anovaBF(amplitude ~ cond*hemisph, data = dataAll, whichRandom = "sbj",iterations = 1000,
        progress=FALSE)

## ####
#SL30 - ROI analyses
rm(list=ls())
setwd("D:/tegolino-derosa/_Experiments/SL30/ROIanalysis_M/")

library(plyr)
library(ez);# Anova
library(car);# Anova
library(lme4) # (g)lmer
library(multcomp) # glht
library(lattice) # dotplot
library(MASS) # boxcox
library(reshape2) # reshape data format
library(ggplot2)
library(CGPfunctions)
library(tidyr)
library(diplyr)
library(magrittr)
library(rstatix)
library(BayesFactor)


filenames <- list.files()
filenames <- filenames[(grep("_base.txt$", filenames))]

for(i in 1:length(filenames))
{
  if (i == 1)
    
  {dataBase <- read.table(filenames[i])
  dataBase$cond <- as.factor(substr(filenames[i], 1,2))
  dataBase$sbj <- row.names(dataBase)}
  else
  {temp <- read.table(filenames[i])
  temp$cond <- as.factor(substr(filenames[i], 1,2))
  temp$sbj <- row.names(temp)
  dataBase <- rbind(dataBase, temp)}
  
}
rm(temp, filenames, i)

colnames(dataBase)[1] <- c("amplitude")

ggplot(dataBase, aes(y = amplitude, x = cond))+ ylab("Amplitude (mV)") +
  geom_boxplot(aes(fill = cond), alpha = .5, outlier.shape = NA) +
  geom_point(aes(col=cond), size=2, position=position_jitterdodge(.7)) +
  theme_bw()

## data are not normal.
dataBase %>%
  group_by(cond) %>%
  shapiro_test(amplitude)


## Electrodes in occipital ROI: occROI = [9,10,11,12, 13,14,15,16,22,23,24,25,26,27,28,29,40,41,39,38];

qqnorm(dataBase$amplitude)
qqline(dataBase$amplitude)

hist(dataBase$amplitude)
dataBase$sbj <- as.factor(dataBase$sbj)
ez1 <- ezANOVA(dataBase, amplitude, sbj, within = .(cond), detailed  = T,type=3, return_aov = TRUE)

ez1


base.aov <- aov(amplitude~cond, data=dataBase)
summary(base.aov)
# Effect       GGe     p[GG] p[GG]<.05       HFe      p[HF] p[HF]<.05
# 2   cond 0.7062023 0.0482432         * 


pwc <- dataBase %>%
  pairwise_t_test(
    amplitude ~ cond, paired = TRUE,
    p.adjust.method = "bonferroni"
  )


res.aov <- anova_test(data = dataBase, dv = amplitude, wid = sbj, within = cond)
get_anova_table(res.aov)





bf <- anovaBF(amplitude ~ cond, data = dataBase, whichRandom = "sbj",iterations = 1000,
              progress=FALSE)


library(BayesFactor)

m_bayes_int <- lmBF(amplitude ~ cond:hemisph  + 
                      cond + hemisph,
                    data=dataBase, whichRandom = c("sbj"))

m_bayes_noint <- lmBF(amplitude ~ 
                        cond + hemisph,
                      data=dataBase, whichRandom = c("sbj"))


summary(m_bayes)


m_bayes_int/m_bayes_noint

