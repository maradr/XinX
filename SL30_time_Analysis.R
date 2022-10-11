rm(list=ls())

setwd('D:/tegolino-derosa/_Experiments/SL30/TimeAnalysis_M/')

library(plyr)
library(dplyr)
library(lme4)
library(lmerTest)
library(MuMIn)
library(effects)
library(ggplot2)

#### Single trials analysis####

dataTrial <- read.table("SL30_Singletrials_withorder.txt", header=T)
summary(dataTrial)

dataTrial$Cond <- as.factor(dataTrial$Cond)


## tidy this up ####
#dataTrial <- read.table("SL30_singletrials.txt", header=T, sep=",")
#summary(dataTrial)

#dataTrial$Cond <- as.factor(dataTrial$Cond)
#dataTrial$Sbj <- as.factor(paste0("sbj_", substr(as.character(dataTrial$Sbj),8,9)))

#left <- dataTrial[, c(1:8)]
#colnames(left)[4:8] <- c("summed_BC", "first", "second", "third", "fourth")
#left$hemisph <- "Left"
#right <- dataTrial[, c(1:3, 9:13)]
#colnames(right)[4:8] <- c("summed_BC", "first", "second", "third", "fourth")
#right$hemisph <- "Right"

#dataTrial <- rbind(left, right)
#rm(left, right)

#write.table(dataTrial, "tidy_SL30_singletrials.txt", row.names=F)



# Visualizations ####

hist(dataTrial$summed_BC,breaks=100)# About normal, couple of outliers
dataTrial$plotting <- paste0(dataTrial$Trial, dataTrial$hemisph)

dataTrial$Trial <- as.factor(dataTrial$Trial)
ggplot(dataTrial, aes(x=Trial, y=summed_BC)) +  
  geom_boxplot(aes(fill=Trial), alpha=0.5, outlier.shape = NA) +  
  geom_point(aes(col=Trial), size = 2, shape = 21, position = position_jitterdodge()) + 
  theme_bw() + xlab("trial order") + ylab("Amplitude") # + theme(legend.position = "none")


ggplot(subset(dataTrial, Cond==1), aes(x=Trial, y=summed_BC, group=plotting)) +  
  geom_boxplot(aes(fill=hemisph), alpha=0.5, outlier.shape = NA) +  
  geom_point(aes(fill=hemisph), size = 2, shape = 21, position = position_jitterdodge()) + 
  theme_bw()  + labs(title="Pseudofonts") + xlab("trial order") + ylab("Amplitude") # + theme(legend.position = "none")
ggplot(subset(dataTrial, Cond==2), aes(x=Trial, y=summed_BC, group=plotting)) +  
  geom_boxplot(aes(fill=hemisph), alpha=0.5, outlier.shape = NA) +  
  geom_point(aes(fill=hemisph), size = 2, shape = 21, position = position_jitterdodge()) + 
  theme_bw()  + labs(title="Nonwords") + xlab("trial order")  + ylab("Amplitude")
ggplot(subset(dataTrial, Cond==3), aes(x=Trial, y=summed_BC, group=plotting)) +  
  geom_boxplot(aes(fill=hemisph), alpha=0.5, outlier.shape = NA) +  
  geom_point(aes(fill=hemisph), size = 2, shape = 21, position = position_jitterdodge()) + 
  theme_bw()   + labs(title="Pseudowords") + xlab("trial order")  + ylab("Amplitude")
ggplot(subset(dataTrial, Cond==4), aes(x=Trial, y=summed_BC, group=plotting)) +  
  geom_boxplot(aes(fill=hemisph), alpha=0.5, outlier.shape = NA) +  
  geom_point(aes(fill=hemisph), size = 2, shape = 21, position = position_jitterdodge()) + 
  theme_bw() + labs(title="Words") + xlab("trial order") + ylab("Amplitude")






## Analysis ####
dataTrial$Trial <- as.factor(dataTrial$Trial)
contrasts(dataTrial$Trial) <-  contr.helmert(5)#helmertContrast


dataTrial$Cond <- as.factor(dataTrial$Cond)
contrasts(dataTrial$Cond) <- contr.sum(4)

dataTrial$hemisph <- as.factor(dataTrial$hemisph)
contrasts(dataTrial$hemisph) <- contr.sum(2)

m_full <- lmer(summed_BC ~ Trial*Cond*hemisph + order + (1|sbj), data=dataTrial)
m_twoway1 <- lmer(summed_BC ~ Trial*Cond + hemisph + order +  (1|sbj), data=dataTrial)
m_twoway2 <- lmer(summed_BC ~ Trial +Cond*hemisph + order + (1|sbj), data=dataTrial)
m_twoway3 <- lmer(summed_BC ~ Trial*hemisph + Cond+ order  + (1|sbj), data=dataTrial)
m_noway <- lmer(summed_BC ~ Trial + Cond + hemisph + order + (1|sbj), data=dataTrial)

anova(m_full, m_twoway1, m_twoway2, m_twoway3, m_noway)

summary(m_twoway1)
anova(m_twoway1)
plot(allEffects(m_twoway1))
r.squaredGLMM(m_twoway1)




library(BayesFactor)

m1bayes <- lmBF(summed_BC ~ Trial + Cond + hemisph, data=dataTrial, whichRandom = c("Sbj"))
summary(m1bayes)

mb_full <- lmBF(summed_BC ~ Trial + Cond + hemisph +
                  Trial:Cond + Trial:hemisph + Cond:hemisph,
                data=dataTrial, whichRandom = c("Sbj"))

mb_twoway <- lmBF(summed_BC ~ Trial + Cond + hemisph +
                  Trial:Cond,
                data=dataTrial, whichRandom = c("Sbj"))


mb_noway <- lmBF(summed_BC ~ Trial + Cond + hemisph, data=dataTrial, whichRandom = c("Sbj"))

mb_full/mb_twoway



library(emmeans)

post <- emmeans(m_twoway1, ~Trial|Cond)
pairs(post)
confint(post)

library(tidyr); library(magrittr)

data.frame(effect(term=c('Trial*Cond'), m_twoway1)) %>% 
  mutate(Cond = case_when(Cond==1 ~ "Pseudofonts", Cond==2~"Nonwords", Cond==3~"Pseudowords", Cond==4~"Words")) %>% 
ggplot(aes(x = Trial, y = fit, col=Cond, group=Cond)) +
  geom_line(size=1.5) + geom_point() +
geom_ribbon(aes(ymin =lower, ymax =upper, group= Cond, fill = Cond), alpha = 0.2, colour = NA)
#geom_smooth(method = "lm", aes(group=Cond, fill=Cond, col=Cond, size = upper-lower), alpha=.15)


# Ramping Up ####
#filenames <- list.files()
#filenames <- filenames[(grep(".*chunked.*.txt$", filenames))]

#for(i in 1:length(filenames))
#{
  tempChunks <- read.table(filenames[i], sep=',', header = T)
  tempChunks$Sbj <- as.factor(paste0("sbj_", substr(as.character(tempChunks$Sbj),8,9)))
  left_1 <- tempChunks[, c(1,2,3)]
  colnames(left_1)[3] <-"summed_BC"
  left_1$chunk <- "first"
  left_1$hemisph <- "left"

  left_2 <- tempChunks[, c(1,2,5)]
  colnames(left_2)[3] <-"summed_BC"
  left_2$chunk <- "second"
  left_2$hemisph <- "left"
  
  left_3 <- tempChunks[, c(1,2,7)]
  colnames(left_3)[3] <-"summed_BC"
  left_3$chunk <- "third"
  left_3$hemisph <- "left"
  
  r_1 <- tempChunks[, c(1,2,4)]
  colnames(r_1)[3] <-"summed_BC"
  r_1$chunk <- "first"
  r_1$hemisph <- "right"
  
  r_2 <- tempChunks[, c(1,2,6)]
  colnames(r_2)[3] <-"summed_BC"
  r_2$chunk <- "second"
  r_2$hemisph <- "right"
  
  r_3 <- tempChunks[, c(1,2,8)]
  colnames(r_3)[3] <-"summed_BC"
  r_3$chunk <- "third"
  r_3$hemisph <- "right"
  if(i==1)
  {
    dataChunks <- rbind(left_1, left_2, left_3, r_1, r_2, r_3)
  }
    
else
  { dataChunks <- rbind(dataChunks, left_1, left_2, left_3, r_1, r_2, r_3)}
  
#}
#rm(tempChunks, left_1, left_2, left_3, r_1, r_2, r_3, filenames, i)

#write.table(dataChunks, "SL30_Chunked_all.txt", row.names = F)

dataChunks <- read.table("SL30_Chunked_all.txt", header=T)

summary(dataChunks)

# Visualizations ####

hist(dataChunks$summed_BC,breaks=100)# About normal, couple of outliers
dataChunks$plotting <- paste0(dataChunks$chunk, dataChunks$hemisph)

dataChunks$chunk <- as.factor(dataChunks$chunk)
ggplot(dataChunks, aes(x=chunk, y=summed_BC)) +  
  geom_boxplot(aes(fill=chunk), alpha=0.5, outlier.shape = NA) +  
  geom_point(aes(col=chunk), size = 2, shape = 21, position = position_jitterdodge()) + 
  theme_bw() + xlab("trial order") + ylab("Amplitude") # + theme(legend.position = "none")


ggplot(subset(dataChunks, cond==1), aes(x=chunk, y=summed_BC, group=plotting)) +  
  geom_boxplot(aes(fill=hemisph), alpha=0.5, outlier.shape = NA) +  
  geom_point(aes(fill=hemisph), size = 2, shape = 21, position = position_jitterdodge()) + 
  theme_bw()  + labs(title="Pseudofonts") + xlab("trial order") + ylab("Amplitude") # + theme(legend.position = "none")

ggplot(subset(dataChunks, cond==2), aes(x=chunk, y=summed_BC, group=plotting)) +  
  geom_boxplot(aes(fill=hemisph), alpha=0.5, outlier.shape = NA) +  
  geom_point(aes(fill=hemisph), size = 2, shape = 21, position = position_jitterdodge()) + 
  theme_bw()  + labs(title="Nonwords") + xlab("trial order")  + ylab("Amplitude")

ggplot(subset(dataChunks, cond==3), aes(x=chunk, y=summed_BC, group=plotting)) +  
  geom_boxplot(aes(fill=hemisph), alpha=0.5, outlier.shape = NA) +  
  geom_point(aes(fill=hemisph), size = 2, shape = 21, position = position_jitterdodge()) + 
  theme_bw()   + labs(title="Pseudowords") + xlab("trial order")  + ylab("Amplitude")

ggplot(subset(dataChunks, cond==4), aes(x=chunk, y=summed_BC, group=plotting)) +  
  geom_boxplot(aes(fill=hemisph), alpha=0.5, outlier.shape = NA) +  
  geom_point(aes(fill=hemisph), size = 2, shape = 21, position = position_jitterdodge()) + 
  theme_bw() + labs(title="Words") + xlab("trial order") + ylab("Amplitude")

## Analysis ####
#set contrasts
dataChunks$cond <- as.factor(dataChunks$cond)
contrasts(dataChunks$cond) <- contr.sum(4)


dataChunks$chunk <- as.factor(dataChunks$chunk)
contrasts(dataChunks$chunk) <- contr.helmert(3)


dataChunks$hemisph <- as.factor(dataChunks$hemisph)
contrasts(dataChunks$hemisph) <- contr.sum(2)



ramping_full <- lmer(summed_BC ~ cond*chunk*hemisph + (1|sbj), data=dataChunks)
ramping_twoway1 <- lmer(summed_BC ~ cond*chunk + hemisph + (1|sbj), data=dataChunks)
ramping_twoway2 <- lmer(summed_BC ~ chunk +cond*hemisph + (1|sbj), data=dataChunks)
ramping_twoway3 <- lmer(summed_BC ~ chunk*hemisph + cond + (1|sbj), data=dataChunks)
ramping_noway <- lmer(summed_BC ~ chunk + cond + hemisph + (1|sbj), data=dataChunks)

anova(ramping_full, ramping_twoway1, ramping_twoway2, ramping_twoway3, ramping_noway)

anova(ramping_twoway1)
summary(ramping_twoway1)
plot(allEffects(ramping_twoway1))
r.squaredGLMM(ramping_twoway1)

eff <- data.frame(effect(term=c('cond*chunk'), ramping_twoway1))


ggplot(eff, aes(x = chunk, y = fit, group = cond)) +
  theme_bw() + 
  geom_line(aes( color=cond), size = 1.2) +
  geom_point(aes(color=cond), size=3) +
  #geom_errorbar(aes(ymin=lower, ymax=upper, col=cond), width=.1) 
     geom_ribbon(aes(ymin = lower, ymax = upper, color = cond, fill = cond), alpha = 0.06) 
 # facet_wrap(~ cond) 


library(BayesFactor)

ramping_full_bayes <- lmBF(summed_BC ~ cond + chunk + hemisph + 
                             cond:chunk + chunk:hemisph + cond:hemisph,
                             data=dataChunks, whichRandom = c("Sbj"))


ramping_twoway_bayes <- lmBF(summed_BC ~ cond + chunk + hemisph + 
                             cond:chunk, 
                             data=dataChunks, whichRandom = c("Sbj"))


ramping_noway_bayes <- lmBF(summed_BC ~ cond + chunk + hemisph, 
                            data=dataChunks, whichRandom = c("Sbj"))

ramping_full_bayes/ramping_twoway_bayes

ramping_twoway_bayes/ramping_noway_bayes



#
#beh <- read.table('D:/tegolino-derosa/_Experiments/SL30/Behavioral/behavioral_summary.txt', header=T)
#summary(beh)
#colnames(beh)[1:2] <- c("sbj", "cond")

#for (i in 1:30)
#{
  if (i == 1)
  {Sbj <- subset(beh, sbj==i & cond %in% c('C1', 'C2', 'C3', 'C4'))
  Sbj <- Sbj[, c(1, 2)]
  Sbj$order <- 1:20
  for(g in 1:4)
    {
    if(g==1)
    {tempSbj <- subset(Sbj, cond==paste0("C", g))
    tempSbj$orderCond <- 1:5}
  else
  {temptempSbj <- subset(Sbj, cond==paste0("C", g))
  temptempSbj$orderCond <- 1:5
  tempSbj <- rbind(tempSbj, temptempSbj)}}
  Sbj <- tempSbj}
  
  else
  {temp <- subset(beh, sbj==i & cond %in% c('C1', 'C2', 'C3', 'C4'))
  temp <- temp[, c(1, 2)]
  temp$order <- 1:20
  for(g in 1:4)
  {
    if(g==1)
    {tempSbj <- subset(temp, cond==paste0("C", g))
    tempSbj$orderCond <- 1:5}
    else
    {temptempSbj <- subset(temp, cond==paste0("C", g))
    temptempSbj$orderCond <- 1:5
    tempSbj <- rbind(tempSbj, temptempSbj)}}
  temp <- tempSbj
  Sbj <- rbind(Sbj, temp)}
 
#}
#Sbj <- subset(Sbj, sbj !="4")
#summary(Sbj)
#Sbj$tomerge <- paste0(Sbj$cond, Sbj$orderCond)

#Sbj$sbj <- ifelse(Sbj$sbj <9, paste0('sbj_0', Sbj$sbj), paste0("sbj_", Sbj$sbj))

#dataTrial$tomerge <- paste0('C', as.character(dataTrial$Cond), as.character(dataTrial$Trial))
#colnames(dataTrial)[1] <- "sbj"
#trial <- merge(Sbj, dataTrial, by=c("sbj", "tomerge"))
#trial$tomerge <- NULL
#write.table(trial, "SL30_Singletrials_withorder.txt", row.names=F)