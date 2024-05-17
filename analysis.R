########pulling and cleaning the csv file
BIGdata<-read.csv("meta_analysis_final_3.csv")
BIGdata
BIGdata<-BIGdata[-c(seq(422,1023, by = 1)),]
View(BIGdata)
BIGdata_unique <- BIGdata[!duplicated(BIGdata[c("Test_number")]), ] 

BIGdata_unique2<-subset(BIGdata_unique, Statistics!="Regression")
####CV
BIGdata_model1<-lm(omega_squared ~ CV_sd_btw_lvls, data = BIGdata_unique)
summary(BIGdata_model1)
plot(omega_squared ~CV_sd_btw_lvls, data = BIGdata_unique)

#####withoutCV
BIGdata_model1<-lm(omega_squared ~ sd_btw_lvls, data = BIGdata_unique)
summary(BIGdata_model1)
plot(omega_squared ~sd_btw_lvls, data = BIGdata_unique)

library(ggplot2)  # load the package

(prelim_plot <- ggplot(BIGdata_unique, aes(x = CV_sd_btw_lvls, y = omega_squared)) +
    geom_point() +
    geom_smooth(method = "lm"))
plot(BIGdata_model1, which = 1)

plot(BIGdata_model1, which = 2)

boxplot(omega_squared ~ Paper, data = BIGdata_unique)

(colour_plot <- ggplot(BIGdata_unique, aes(x = CV_sd_btw_lvls, y = omega_squared, colour = Paper)) +
    geom_point(size = 2) +
    theme_classic() +
    theme(legend.position = "none"))

paper.lm <- lm(omega_squared ~ CV_sd_btw_lvls + Paper, data = BIGdata_unique)
summary(paper.lm)

library(lme4)

mixed.lmer <- lmer(omega_squared ~ CV_sd_btw_lvls + (1|Paper), data = BIGdata_unique)
summary(mixed.lmer)


mixed.lmer <- lmer(omega_squared ~ CV_sd_btw_lvls + (1|Paper), data = BIGdata_unique)
mixed.lmer <- lmer(omega_squared ~ sd_btw_lvls + (1|Paper), data = BIGdata_unique)

#####withoutregression

BIGdata_model2<-lm(omega_squared ~ CV_sd_btw_lvls, data = BIGdata_unique2)
logomega_squared<-log(BIGdata_unique2$omega_squared +1) ####log data
BIGdata_unique2<- cbind(BIGdata_unique2, logomega_squared) ####log data
BIGdata_modelLOG<-lm(logomega_squared ~ CV_sd_btw_lvls, data = BIGdata_unique2)
summary(BIGdata_model2)
summary(BIGdata_modelLOG)
plot(omega_squared ~CV_sd_btw_lvls, data = BIGdata_unique2)

(prelim_plot2 <- ggplot(BIGdata_unique2, aes(x = CV_sd_btw_lvls, y = omega_squared)) +
        geom_point() +
        geom_smooth(method = "lm"))

(logprelim_plot2 <- ggplot(BIGdata_unique2, aes(x = CV_sd_btw_lvls, y = logomega_squared)) +
        geom_point() +
        geom_smooth(method = "lm"))               #####log plot

#####withoutregression and w/o CV

BIGdata_model2<-lm(omega_squared ~ sd_btw_lvls, data = BIGdata_unique2)
logomega_squared<-log(BIGdata_unique2$omega_squared +1)
View(logomega_squared)
View(logomega_squared2)
logomega_squared2<-log(BIGdata_unique2$sd_btw_lvls +1)####log data
BIGdata_unique2<- cbind(BIGdata_unique2, logomega_squared,logomega_squared2) ####log data
BIGdata_modelLOG<-lm(logomega_squared ~ logomega_squared2, data = BIGdata_unique2)
summary(BIGdata_model2)
summary(BIGdata_modelLOG)
plot(logomega_squared ~logomega_squared2, data = BIGdata_unique2)

(prelim_plot2 <- ggplot(BIGdata_unique2, aes(x = sd_btw_lvls, y = omega_squared)) +
        geom_point() +
        geom_smooth(method = "lm"))

(logprelim_plot2 <- ggplot(BIGdata_unique2, aes(x = logomega_squared2, y = logomega_squared)) +
        geom_point() +
        geom_smooth(method = "lm"))               #####log plot

####normality tests
plot(BIGdata_model2, which = 1)
plot(BIGdata_modelLOG, which =1)
plot(BIGdata_model2, which = 2)
plot(BIGdata_modelLOG, which = 2)

plot(lm(omega_squared~CV_sd_btw_lvls,data=BIGdata_unique2))

shapiro.test(BIGdata_unique2$omega_squared)
shapiro.test(BIGdata_unique2$CV_sd_btw_lvls)

install.packages("lmtest")
library(lmtest)
bptest(BIGdata_model2)
bptest(BIGdata_modelLOG)
boxplot(omega_squared ~ Paper, data = BIGdata_unique2)

###### colour tests and mixed model

(colour_plot <- ggplot(BIGdata_unique2, aes(x = CV_sd_btw_lvls, y = omega_squared, colour = Paper)) +
        geom_point(size = 2) +
        theme_classic() +
        theme(legend.position = "none"))

(colour_plot <- ggplot(BIGdata_unique2, aes(x = CV_sd_btw_lvls, y = logomega_squared, colour = Paper)) +
        geom_point(size = 2) +
        theme_classic() +
        theme(legend.position = "none"))
plot(logomega_squared ~CV_sd_btw_lvls, data = BIGdata_unique2) ###logplot

paper.lm2 <- lm(omega_squared ~ CV_sd_btw_lvls + Paper, data = BIGdata_unique2)
summary(paper.lm2) #no transform
paper.lmLOG <- lm(logomega_squared ~ CV_sd_btw_lvls + Paper, data = BIGdata_unique2)
summary(paper.lmLOG) #log transform

mixed.lmer2 <- lmer(omega_squared ~ CV_sd_btw_lvls + (1|Paper), data = BIGdata_unique2)
summary(mixed.lmer2) #no transform, random factor paper
mixed.lmerlog <- lmer(logomega_squared ~ CV_sd_btw_lvls + (1|Paper), data = BIGdata_unique2)
summary(mixed.lmerlog) #log transform random factor paper

#instead of effect size, analyzing effect of CV on Total_levels_per_factor)
mixed.lmer2F <- lmer(Total_Levels_per_Factor ~ CV_sd_btw_lvls+(1|Paper), data = BIGdata_unique2)
summary(mixed.lmer2F) 
plot(Total_Levels_per_Factor ~CV_sd_btw_lvls, data = BIGdata_unique2)
(prelim_plot3 <- ggplot(BIGdata_unique2, aes(x = CV_sd_btw_lvls, y = Total_Levels_per_Factor)) +
        geom_point() +
        geom_smooth(method = "lm"))

mixed.lmer2U <- lmer(omega_squared~Total_Levels_per_Factor+(1|Paper), data = BIGdata_unique2)
summary(mixed.lmer2U)
plot(omega_squared~Total_Levels_per_Factor, data = BIGdata_unique2)
(prelim_plot4 <- ggplot(BIGdata_unique2, aes(x = Total_Levels_per_Factor, y = omega_squared)) +
        geom_point() +
        geom_smooth(method = "lm"))

#####sample analysis + trial analysis + trial plots
mixed.lmer2i <- lmer(omega_squared ~ CV_sd_btw_lvls+(1|taxa), data = BIGdata_unique2)
summary(mixed.lmer2i)

taxa.lm2 <- lm(omega_squared ~ CV_sd_btw_lvls + taxa, data = BIGdata_unique2)
summary(taxa.lm2)

(colour_plot2 <- ggplot(BIGdata_unique2, aes(x = CV_sd_btw_lvls, y = omega_squared, colour = taxa)) +
        geom_point(size = 2) +
        theme_classic() +
        theme(legend.position = "none"))

(colour_plot <- ggplot(BIGdata_unique2, aes(x = CV_sd_btw_lvls, y = omega_squared, colour = Total_Levels_per_Factor)) +
        geom_point(size = 2) +
        theme_classic() +
        theme(legend.position = "none"))

TLF.lm2 <- lm(omega_squared ~ CV_sd_btw_lvls + as.factor(Total_Levels_per_Factor), data = BIGdata_unique2)
summary(TLF.lm2)

TLF2.lm2 <- lm((CV_sd_btw_lvls) ~ Total_Levels_per_Factor, data = BIGdata_unique2)
summary(TLF2.lm2)

(prelim_plot5 <- ggplot(BIGdata_unique2, aes(x = Total_Levels_per_Factor, y = (CV_sd_btw_lvls))) +
        geom_point() +
        geom_smooth(method = "lm"))

mixed.lmerlog21 <- lmer(logomega_squared ~ Total_Levels_per_Factor+(1|Paper), data = BIGdata_unique2)
summary(mixed.lmerlog21)

library(nlme)

model_1<- lme(logomega_squared ~ CV_sd_btw_lvls*Total_Levels_per_Factor, data = BIGdata_unique2, 
              random = ~1|Paper )
summary(model_1)
anova(model_1)

library(jtools)
summ(model_1)
model_2 <-lmer(logomega_squared~CV_sd_btw_lvls*Total_Levels_per_Factor+(1|Paper), data =BIGdata_unique2)
summ(model_2)               

model_3 <-lmer(logomega_squared~Total_Factors+(1|Paper), data =BIGdata_unique2)
summ(model_3)

model_4 <- lm(omega_squared~CV_sd_btw_lvls*Total_Levels_per_Factor, data =BIGdata_unique2)
summ(model_4)

model_5 <-lmer(logomega_squared~poly(CV_sd_btw_lvls, 2)+(1|Paper), data =BIGdata_unique2)
summ(model_5) 

prelim_plot6 <- ggplot(BIGdata_unique2, aes(x = CV_sd_btw_lvls, y = rank(omega_squared))) +
        geom_point() +
        geom_smooth(method = "lm", formula = y ~ x+I(x^2))

model_6 <- lmer(rank(omega_squared)~(CV_sd_btw_lvls) + (1|Paper), data = BIGdata_unique2)
summ(model_6)

(colour_plot3 <- ggplot(BIGdata_unique2, aes(x = CV_sd_btw_lvls, y = logomega_squared, colour = Paper))) +
        geom_point(size = 2) +
        theme_classic() +xlab("Variance of X")+ ylab("log (??2)") +ggtitle("Effect of Variance of X on Effect Size")
################################end of trials

#####these are  the plots using CV

prelim_plot7 <- ggplot(BIGdata_unique2, aes(x = CV_sd_btw_lvls, y = rank(omega_squared))) +
    geom_point() +
    geom_smooth(method = "lm") + xlab("Coefficient of Variation") + ylab(" ranked (ω^2)") + ylim(0,115)+
    xlim(0,1.5) + theme_classic()  + theme(axis.title = element_text(size = 18))
prelim_plot8 <- ggplot(BIGdata_unique2, aes(x = CV_sd_btw_lvls, y = rank(omega_squared),colour = Paper)) +
    geom_point() +xlab("Coefficient of Variation") + ylab(" ranked (ω^2)") + ylim(0,115)+
    xlim(0,1.5)  + theme_classic()  + theme(axis.title = element_text(size = 18)) +theme(legend.title = element_text(size=18))
   
###############plots without cv and logsd

prelim_plot9 <- ggplot(BIGdata_unique2, aes(x = log(sd_btw_lvls), y = rank(omega_squared))) +
    geom_point() +
    geom_smooth(method = "lm") + xlab("Coefficient of Variation") + ylab(" ranked (ω^2)") + ylim(0,115)+
    xlim(0,1.5) + theme_classic()  + theme(axis.title = element_text(size = 18))
prelim_plot10 <- ggplot(BIGdata_unique2, aes(x = log(sd_btw_lvls), y = rank(omega_squared),colour = Paper)) +
    geom_point() +xlab("Coefficient of Variation") + ylab(" ranked (ω^2)") + ylim(0,115)+
    xlim(0,1.5)  + theme_classic()  + theme(axis.title = element_text(size = 18)) +theme(legend.title = element_text(size=18))
##################with stats

model_7 <- lmer(rank(omega_squared)~(log(sd_btw_lvls)) + (1|Paper), data = BIGdata_unique2)
summ(model_7)

####################rank SD instead of Cv

prelim_plot11 <- ggplot(BIGdata_unique2, aes(x = rank(sd_btw_lvls), y = rank(omega_squared))) +
    geom_point() +
    geom_smooth(method = "lm") + xlab("Ranked Standard Deviation") + ylab(" ranked (ω^2)") + theme_classic()  + theme(axis.title = element_text(size = 18))
prelim_plot12 <- ggplot(BIGdata_unique2, aes(x = log(sd_btw_lvls), y = rank(omega_squared),colour = Paper)) +
    geom_point() +xlab("Coefficient of Variation") + ylab(" ranked (ω^2)") + ylim(0,115)+
    xlim(0,1.5)  + theme_classic()  + theme(axis.title = element_text(size = 18)) +theme(legend.title = element_text(size=18))

model_8 <- lmer(rank(omega_squared)~(rank(sd_btw_lvls)) + (1|Paper), data = BIGdata_unique2)
summ(model_8)
summary(model_8)
