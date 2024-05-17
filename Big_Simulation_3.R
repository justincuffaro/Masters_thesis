# calculate the slope of gain in effect size
# by increasing variance and number of replicates
# according to units of variance.  

generate.matrix.normal <- function(vector.n,vector.mean,vector.sd){
  n.replicates <- length(vector.n)
  matrix.rnd.norm <- list()
  treatments <- list()
  for (i in 1:n.replicates){
    matrix.rnd.norm[[i]] <- rnorm(vector.n[i],vector.mean[i],vector.sd[i])
    treatments[[i]] <- rep(i,vector.n[i])
  }
  result <- list(sim.values=matrix.rnd.norm,treatments=treatments)
  return(result)
}

create.design.linear <- function(slope,var.X,n.treatments,n.replicates){
  # n.replicates is per treatment and treatments are assumed to be balanced, i.e.,
  # same number of replicates for each treatment.
  
  sim.X <- sort(rnorm(n.treatments,1,var.X))
  # create treatments with constant ranges within X
  f <- unique(cut(seq(sim.X[1],sim.X[n.treatments], by = 0.001),breaks=n.treatments))
  labs <- levels(f)[f]
  lower <- as.numeric( sub("\\((.+),.*", "\\1", labs))
  # upper <- as.numeric( sub("[^,]*,([^]]*)\\]", "\\1", labs))
  values.X <- lower
  Y.expected.means <- slope * values.X 
  error <- 1
  sample.size <- n.replicates
  # balanced design (here, for now)
  vector.n <- rep(sample.size,n.treatments)
  vector.sd <- rep(error,n.treatments)
  sim.data <- generate.matrix.normal(vector.n,Y.expected.means,vector.sd)
  treatments <- unlist(sim.data$treatments)
  sim.data <- unlist(sim.data$sim.values)
  # tapply(sim.data,treatments,mean)
  # ANOVA:
  
  F.values <- matrix(0,2)
  p.values <- matrix(0,2)
  names(F.values) <- c("anova","regression")
  names(p.values) <-c("anova","regression")
  
  aov.standard <- anova(lm(sim.data~as.factor(treatments)))
  F.values["anova"] <- aov.standard$`F value`[1]
  p.values["anova"] <- aov.standard$`Pr(>F)`[1]
  # regression:
  X <- rep(values.X,times=vector.n)
  
  aov.regression <- anova(lm(sim.data ~ X))
  F.values["regression"] <- aov.regression$`F value`[1]
  p.values["regression"] <- aov.regression$`Pr(>F)`[1]
  results <- list(F.values=F.values,p.values=p.values)
  return(results)
}

# demonstrating that the function create.design.linear 
# gives the same results as the original loop below
#################################
##remove and replace hashmarks according to parameters of the simulation
################################
n.sims <- 1000
#n.treatments<-10
#n.treatments <-5
n.treatments <-15
###########################################   
############## first result ############### 
###########################################  
# N = 50 and N = 100, i.e., 5 treatments * 10 replicates each and 5*20
sim.slopes <- c(0.05,0.1,0.2)
#n.replicates <- c(3,6)
#n.replicates <- c(6,12)
n.replicates <- c(2,4)
var.X <- c(1,3,5,7,9)

####################################################

total.runs <- n.sims * length(sim.slopes) * length(n.replicates) * length(var.X) 
result.matrix.F <- matrix(0,total.runs,12)
colnames(result.matrix.F) <- c("slopes","n.replicates","var.X","ANOVA","regression","levels_per_treatment",
                               "eta_squaredANOVA","eta_squaredRegression","omega_squaredANOVA","omega_squaredRegression",
                               "R.2ANOVA","R.2Regression")
result.matrix.p <- matrix(0,total.runs,6)
colnames(result.matrix.p) <- c("slopes","n.replicates","var.X","ANOVA","regression","levels_per_treatment")


# we will be running some "waste" here as we won't use all the combinations
# simulated results for graphing, but makes the code more general
#the simulation

count.runs <- 1
for (count.slopes in 1:length(sim.slopes)){
  for (count.n.replicates in 1:length(n.replicates)){
    for (count.sim in 1:n.sims){     
      for (count.var in 1:length(var.X)){
        print(paste0(count.runs, " out of ",total.runs))
        res <- create.design.linear(slope=sim.slopes[count.slopes],var.X=var.X[count.var],
                                    n.treatments=n.treatments,n.replicates=n.replicates[count.n.replicates])
        
        result.matrix.F[count.runs,"slopes"] <- sim.slopes[count.slopes]
        result.matrix.F[count.runs,"n.replicates"] <- n.replicates[count.n.replicates]
        result.matrix.F[count.runs,"var.X"] <- var.X[count.var]
        result.matrix.F[count.runs,"ANOVA"] <- res$F.values["anova"]
        result.matrix.F[count.runs,"regression"] <- res$F.values["regression"]
        result.matrix.F[count.runs,"levels_per_treatment"] <- n.treatments
        result.matrix.F[count.runs,"eta_squaredANOVA"]<- 
          (res$F.values["anova"])*(n.treatments-1)/((res$F.values["anova"])*(n.treatments-1) 
                                +((n.replicates[count.n.replicates])*n.treatments-n.treatments))
        result.matrix.F[count.runs,"omega_squaredANOVA"]<- 
          (res$F.values["anova"]-1)/((res$F.values["anova"]) + 
               (((n.replicates[count.n.replicates])*n.treatments-n.treatments)+1)/(n.treatments -1))
          
        result.matrix.F[count.runs,"eta_squaredRegression"]<-
          (res$F.values["regression"])*(n.treatments-1)/((res$F.values["regression"])*(n.treatments-1) 
                                                    +((n.replicates[count.n.replicates])*n.treatments-n.treatments))
        result.matrix.F[count.runs,"omega_squaredRegression"]<-
          (res$F.values["regression"]-1)/((res$F.values["regression"]) + 
             (((n.replicates[count.n.replicates])*n.treatments-n.treatments)+1)/(n.treatments -1))
        result.matrix.F[count.runs,"R.2ANOVA"]<-
          1 - (1/(1+res$F.values["anova"]*((n.treatments-1)/(
            ((n.replicates[count.n.replicates])*n.treatments)-n.treatments))))
          
        result.matrix.F[count.runs,"R.2Regression"]<-
          1 - (1/(1+res$F.values["regression"]*((n.treatments-1)/(
            ((n.replicates[count.n.replicates])*n.treatments)-n.treatments))))
        
        result.matrix.p[count.runs,"slopes"] <- sim.slopes[count.slopes]
        result.matrix.p[count.runs,"n.replicates"] <- n.replicates[count.n.replicates]
        result.matrix.p[count.runs,"var.X"] <- var.X[count.var]
        result.matrix.p[count.runs,"ANOVA"] <- res$p.values["anova"]
        result.matrix.p[count.runs,"regression"] <- res$p.values["regression"]
        result.matrix.p[count.runs,"levels_per_treatment"] <- n.treatments
        
        count.runs <- count.runs + 1
      }
    }
  }
}

library(dplyr)
library(ggplot2)

#to make bar plots

result.matrix.p <- data.frame(result.matrix.p)
bar.plots <- vector('list',length(sim.slopes)*length(n.replicates))
data.df <- vector('list',length(sim.slopes)*length(n.replicates))
count <- 1
for (count.slopes in 1:length(sim.slopes)){
  for (count.n.replicates in 1:length(n.replicates)){
    # filter doesn't take vector values
    s <- sim.slopes[count.slopes]
    r <- n.replicates[count.n.replicates]
    res.anova <- data.frame(result.matrix.p %>% filter(slopes == s, n.replicates == r) %>% 
                              group_by(var.X) %>% count(ANOVA < 0.05) %>% filter(`ANOVA < 0.05`==TRUE))
    
    res.reg <- data.frame(result.matrix.p %>% filter(slopes == s, n.replicates == r) %>% 
                            group_by(var.X) %>% count(regression < 0.05) %>% filter(`regression < 0.05`==TRUE))
    res.anova <- res.anova[,-2]
    res.anova <- cbind(res.anova,type=rep("anova",nrow(res.anova)))
    res.reg <- res.reg[,-2]
    res.reg <- cbind(res.reg,type=rep("regression",nrow(res.reg)))
    df <- data.frame(rbind(res.anova,res.reg))
    df[,"n"] <- df[,"n"]/n.sims
    data.df[[count]] <- df
    bar.plots[[count]] <- local({
      count <- count
      sp <- ggplot(data=df, aes(fill=type, y=n,x=as.factor(var.X))) +
        geom_bar(aes(fill=type),position="dodge",stat="identity") + 
        xlab("Variance") + ylab("Power") +
        scale_fill_manual("legend", values = c("anova" = "blue", "regression" = "firebrick")) +
        theme_classic() +
        theme(text = element_text(size = 14)) + 
        theme(axis.text = element_text(size = 14)) +
        ggtitle(bquote("n=" ~ .(n.replicates[count.n.replicates]) ~ " slope=" ~ .(sim.slopes[count.slopes]))) +
        theme(plot.title = element_text(hjust = 0.5)) + 
        ylim(0,1) +
        guides(fill=guide_legend(title="  ")) + 
        theme(legend.text = element_text(size=16),
              legend.key.height = unit(0.5, 'cm'),
              legend.key.width = unit(0.5, 'cm'))
              
    })
    count <- count + 1
  }
}

library(patchwork)
library(ggpubr)

###############################################################
###THE BIG PLOT USES THREE DIFFERENT SIMULATION RUNS AS TREATMENT AND REPLICATES VARY
###############################################################
#######archival purposes######use BIGPLOT's 0-> 2 for custom treatment and replicate numbers
BIGPLOT<-ggarrange(bar.plots[[1]]+rremove("xlab"),
                    bar.plots[[2]]+rremove("ylab")+rremove("xlab"),
                    bar.plots[[3]]+rremove("xlab")+rremove("xlab"),
                    bar.plots[[4]]+rremove("ylab")+rremove("xlab"),
                    bar.plots[[5]],
                    bar.plots[[6]]+rremove("ylab"),
                    ncol=2,nrow=3,common.legend = TRUE)

BIGPLOT1<-ggarrange(bar.plots[[1]]+rremove("xlab"),
                    bar.plots[[2]]+rremove("ylab")+rremove("xlab"),
                    bar.plots[[3]]+rremove("xlab")+rremove("xlab"),
                    bar.plots[[4]]+rremove("ylab")+rremove("xlab"),
                    bar.plots[[5]],
                    bar.plots[[6]]+rremove("ylab"),
                    ncol=2,nrow=3,common.legend = TRUE)

BIGPLOT2<-ggarrange(bar.plots[[1]]+rremove("xlab"),
                    bar.plots[[2]]+rremove("ylab")+rremove("xlab"),
                    bar.plots[[3]]+rremove("xlab")+rremove("xlab"),
                    bar.plots[[4]]+rremove("ylab")+rremove("xlab"),
                    bar.plots[[5]],
                    bar.plots[[6]]+rremove("ylab"),
                    ncol=2,nrow=3,common.legend = TRUE)

########archival purposes
####################################################
####################################################
#######Use Big Plot 3 for when treatment is 5, replicates are 6 and 12
#######Use Big Plot 4 for when treatment is 10, replicates 3 and 6
#######Use Big Plot 5 for when treatment is 15 replicates 2 and 4

BIGPLOT3<-ggarrange(bar.plots[[1]]+rremove("xlab"),
                    bar.plots[[2]]+rremove("ylab")+rremove("xlab"),
                    bar.plots[[3]]+rremove("xlab")+rremove("xlab"),
                    bar.plots[[4]]+rremove("ylab")+rremove("xlab"),
                    bar.plots[[5]],
                    bar.plots[[6]]+rremove("ylab"),
                    ncol=2,nrow=3,common.legend = TRUE)

BIGPLOT4<-ggarrange(bar.plots[[1]]+rremove("xlab"),
                    bar.plots[[2]]+rremove("ylab")+rremove("xlab"),
                    bar.plots[[3]]+rremove("xlab")+rremove("xlab"),
                    bar.plots[[4]]+rremove("ylab")+rremove("xlab"),
                    bar.plots[[5]],
                    bar.plots[[6]]+rremove("ylab"),
                    ncol=2,nrow=3,common.legend = TRUE)

BIGPLOT5<-ggarrange(bar.plots[[1]]+rremove("xlab"),
          bar.plots[[2]]+rremove("ylab")+rremove("xlab"),
          bar.plots[[3]]+rremove("xlab")+rremove("xlab"),
          bar.plots[[4]]+rremove("ylab")+rremove("xlab"),
          bar.plots[[5]],
          bar.plots[[6]]+rremove("ylab"),
          ncol=2,nrow=3,common.legend = TRUE)

BIGPLOT_A<-annotate_figure(BIGPLOT,
                top = text_grob("Number of Treatments=3", color = "black", face = "bold", size = 18),
)
BIGPLOT1_A<-annotate_figure(BIGPLOT1,
                           top = text_grob("Number of Treatments=6", color = "black", face = "bold", size = 18),
)
BIGPLOT2_A<-annotate_figure(BIGPLOT2,
                           top = text_grob("Number of Treatments=9", color = "black", face = "bold", size = 18),
)
Largeplot<-BIGPLOT_A+BIGPLOT1_A+BIGPLOT2_A

BIGPLOT3_A<-annotate_figure(BIGPLOT3,
                            top = text_grob("Number of Treatments=5", color = "black", face = "bold", size = 18),
)
BIGPLOT4_A<-annotate_figure(BIGPLOT4,
                            top = text_grob("Number of Treatments=10", color = "black", face = "bold", size = 18),
)
BIGPLOT5_A<-annotate_figure(BIGPLOT5,
                            top = text_grob("Number of Treatments=15", color = "black", face = "bold", size = 18),
)

Largeplot2<-BIGPLOT3_A+BIGPLOT4_A+BIGPLOT5_A

######must run three simulations, each time
#####################################################
#####################################################
######All other graphs use treatment = 10 and replicates = 3,6

#for F values

result.matrix.F <- data.frame(result.matrix.F)
bar.plots <- vector('list',length(sim.slopes)*length(n.replicates))
data.df <- vector('list',length(sim.slopes)*length(n.replicates))
count <- 1
for (count.slopes in 1:length(sim.slopes)){
  for (count.n.replicates in 1:length(n.replicates)){
    # filter doesn't take vector values
    s <- sim.slopes[count.slopes]
    r <- n.replicates[count.n.replicates]
    res.anova <- data.frame(result.matrix.F %>% filter(slopes == s, n.replicates == r) %>% 
                              group_by(var.X))
    
    res.reg <- data.frame(result.matrix.F %>% filter(slopes == s, n.replicates == r) %>% 
                            group_by(var.X))
    res.anova <- res.anova[,-5]
    res.anova <- cbind(res.anova,type=rep("anova",nrow(res.anova)))
    res.reg <- res.reg[,-4]
    res.reg <- cbind(res.reg,type=rep("regression",nrow(res.reg)))
    res.anova.F<-aggregate(ANOVA~var.X,res.anova, FUN = mean)
    colnames(res.anova.F)<- c("var.X","average")
    res.anova.F <- cbind(res.anova.F,type=rep("anova",nrow(res.anova.F)))
    res.reg.F<-aggregate(regression~var.X,res.reg, FUN = mean)
    colnames(res.reg.F)<- c("var.X","average")
    res.reg.F <- cbind(res.reg.F,type=rep("regression",nrow(res.anova.F)))
    df <- data.frame(rbind(res.anova.F,res.reg.F))
    data.df[[count]] <- df
    bar.plots[[count]] <- local({
      count <- count
      sp <- ggplot(data=df, aes(colour=type, y=average,x=as.factor(var.X))) +
        geom_point(aes(colour=type),position="identity",stat="identity") + 
        xlab("variance of treatment") + ylab("average_F-value") +
        scale_colour_manual("legend", values = c("anova" = "blue", "regression" = "firebrick")) +
        theme_classic() +
        theme(text = element_text(size = 10)) + 
        theme(axis.text = element_text(size = 10)) +
        ggtitle(bquote("n=" ~ .(n.replicates[count.n.replicates]) ~ " slope=" ~ .(sim.slopes[count.slopes]))) +
        theme(plot.title = element_text(hjust = 0.5)) + 
        guides(fill=guide_legend(title="  ")) +
        theme(legend.text = element_text(size=12),
              legend.key.height = unit(0.5, 'cm'),
              legend.key.width = unit(0.5, 'cm'))
    })
    count <- count + 1
  }
}

FPLOT1<-ggarrange(bar.plots[[1]]+rremove("xlab"), 
                    bar.plots[[2]]+rremove("ylab")+rremove("xlab"),
                    bar.plots[[3]]+rremove("xlab")+rremove("xlab"),
                    bar.plots[[4]]+rremove("ylab")+rremove("xlab"),
                    bar.plots[[5]],
                    bar.plots[[6]]+rremove("ylab"),
                    ncol=2,nrow=3,common.legend = TRUE)
FPLOT1A<-annotate_figure(FPLOT1,
                        top = text_grob("Number of Treatments=10", color = "black", face = "bold", size = 14))

####eta-squared###########THIS SECTION IS OUT OF USE##########################################
###instead of eta-squared or R^2 I will be using OMEGASQUARED#################################
result.matrix.F <- data.frame(result.matrix.F)
scatter.plots <- vector('list',length(sim.slopes)*length(n.replicates))
data.df <- vector('list',length(sim.slopes)*length(n.replicates))
count <- 1
for (count.slopes in 1:length(sim.slopes)){
  for (count.n.replicates in 1:length(n.replicates)){
    # filter doesn't take vector values
    s <- sim.slopes[count.slopes]
    r <- n.replicates[count.n.replicates]
    res.anova <- data.frame(result.matrix.F %>% filter(slopes == s, n.replicates == r) %>% 
                              group_by(var.X))
    
    res.reg <- data.frame(result.matrix.F %>% filter(slopes == s, n.replicates == r) %>% 
                            group_by(var.X))
    res.anova <- res.anova[,-5]
    res.anova<-res.anova[,-7]
    names(res.anova)[names(res.anova) == 'eta.squaredANOVA'] <- 'eta.squared'
    names(res.anova)[names(res.anova) == 'ANOVA'] <- 'F'
    res.anova.eta<-aggregate(eta.squared~var.X,res.anova, FUN = mean)
    colnames(res.anova.eta)<- c("var.X","averageETA")
    res.anova.eta <- cbind(res.anova.eta,type=rep("anova",nrow(res.anova.eta)))
    res.reg <- res.reg[,-4]
    res.reg <- res.reg[,-6]
    names(res.reg)[names(res.reg) == 'eta.squaredRegression'] <- 'eta.squared'
    names(res.reg)[names(res.reg) == 'regression'] <- 'F'
    res.reg.eta<-aggregate(eta.squared~var.X,res.reg, FUN = mean)
    colnames(res.reg.eta)<- c("var.X","averageETA")
    res.reg.eta <- cbind(res.reg.eta,type=rep("regression",nrow(res.reg.eta)))
    df <- data.frame(rbind(res.anova.eta,res.reg.eta))
    data.df[[count]] <- df
    scatter.plots[[count]] <- local({
      count <- count
      sp <- ggplot(data=df, aes(colour=type, y=averageETA,x=as.factor(var.X))) +
        geom_point(aes(colour=type),position="identity",stat="identity") + 
        xlab("variance of treatment") + ylab("average_eta-squared") +
        scale_colour_manual("legend", values = c("anova" = "blue", "regression" = "firebrick")) +
        theme_classic() +
        theme(text = element_text(size = 10)) + 
        theme(axis.text = element_text(size = 10)) +
        ggtitle(bquote("n=" ~ .(n.replicates[count.n.replicates]) ~ " slope=" ~ .(sim.slopes[count.slopes]))) +
        theme(plot.title = element_text(hjust = 0.5)) + 
        guides(fill=guide_legend(title="  ")) +
        theme(legend.text = element_text(size=12),
              legend.key.height = unit(0.5, 'cm'),
              legend.key.width = unit(0.5, 'cm'))
    })
    count <- count + 1
  }
}

etaPLOT1<-ggarrange(scatter.plots[[1]]+rremove("xlab"), 
                    scatter.plots[[2]]+rremove("ylab")+rremove("xlab"),
                    scatter.plots[[3]]+rremove("xlab")+rremove("xlab"),
                    scatter.plots[[4]]+rremove("ylab")+rremove("xlab"),
                    scatter.plots[[5]],
                    scatter.plots[[6]]+rremove("ylab"),
                  ncol=2,nrow=3,common.legend = TRUE)
etaPLOT1A<-annotate_figure(etaPLOT1,
                         top = text_grob("Number of Treatments=10", color = "black", face = "bold", size = 14))


#####R^2

result.matrix.F <- data.frame(result.matrix.F)
statter.plots <- vector('list',length(sim.slopes)*length(n.replicates))
data.df <- vector('list',length(sim.slopes)*length(n.replicates))
count <- 1
for (count.slopes in 1:length(sim.slopes)){
  for (count.n.replicates in 1:length(n.replicates)){
    # filter doesn't take vector values
    s <- sim.slopes[count.slopes]
    r <- n.replicates[count.n.replicates]
    res.anova <- data.frame(result.matrix.F %>% filter(slopes == s, n.replicates == r) %>% 
                              group_by(var.X))
    res.anova <- res.anova[,-5]
    res.anova<-res.anova[,-7]
    res.anova<-res.anova[,-8]
    
    res.reg <- data.frame(result.matrix.F %>% filter(slopes == s, n.replicates == r) %>% 
                            group_by(var.X))
    res.reg<-res.reg[,-4]
    res.reg<-res.reg[,-6]
    res.reg<-res.reg[,-7]
    names(res.anova)[names(res.anova) == 'R.2ANOVA'] <- 'R.2'
    names(res.anova)[names(res.anova) == 'ANOVA'] <- 'F'
    res.anova.R<-aggregate(R.2~var.X,res.anova, FUN = mean)
    colnames(res.anova.eta)<- c("var.X","averageR.2")
    res.anova.R <- cbind(res.anova.R,type=rep("anova",nrow(res.anova.R)))
    names(res.reg)[names(res.reg) == 'R.2Regression'] <- 'R.2'
    names(res.reg)[names(res.reg) == 'regression'] <- 'F'
    res.reg.R<-aggregate(R.2~var.X,res.reg, FUN = mean)
    colnames(res.reg.R)<- c("var.X","R.2")
    res.reg.R <- cbind(res.reg.R,type=rep("regression",nrow(res.reg.R)))
    df <- data.frame(rbind(res.anova.R,res.reg.R))
    data.df[[count]] <- df
    statter.plots[[count]] <- local({
      count <- count
      sp <- ggplot(data=df, aes(colour=type, y=R.2,x=as.factor(var.X))) +
        geom_point(aes(colour=type),position="identity",stat="identity") + 
        xlab("variance of treatment") + ylab("average_R.2") +
        scale_colour_manual("legend", values = c("anova" = "blue", "regression" = "firebrick")) +
        theme_classic() +
        theme(text = element_text(size = 10)) + 
        theme(axis.text = element_text(size = 10)) +
        ggtitle(bquote("n=" ~ .(n.replicates[count.n.replicates]) ~ " slope=" ~ .(sim.slopes[count.slopes]))) +
        theme(plot.title = element_text(hjust = 0.5)) + 
        guides(fill=guide_legend(title="  ")) +
        theme(legend.text = element_text(size=12),
              legend.key.height = unit(0.5, 'cm'),
              legend.key.width = unit(0.5, 'cm'))
    })
    count <- count + 1
  }
}

RPLOT1<-ggarrange(statter.plots[[1]]+rremove("xlab"), 
                    statter.plots[[2]]+rremove("ylab")+rremove("xlab"),
                    statter.plots[[3]]+rremove("xlab")+rremove("xlab"),
                    statter.plots[[4]]+rremove("ylab")+rremove("xlab"),
                    statter.plots[[5]],
                    statter.plots[[6]]+rremove("ylab"),
                    ncol=2,nrow=3,common.legend = TRUE)
RPLOT1A<-annotate_figure(RPLOT1,
                           top = text_grob("Number of Treatments=10", color = "black", face = "bold", size = 14))

####################################################################################################
####################################################################################################

###omegasquared

result.matrix.F <- data.frame(result.matrix.F)
omega.plots <- vector('list',length(sim.slopes)*length(n.replicates))
data.df <- vector('list',length(sim.slopes)*length(n.replicates))
count <- 1
for (count.slopes in 1:length(sim.slopes)){
  for (count.n.replicates in 1:length(n.replicates)){
    # filter doesn't take vector values
    s <- sim.slopes[count.slopes]
    r <- n.replicates[count.n.replicates]
    res.anova <- data.frame(result.matrix.F %>% filter(slopes == s, n.replicates == r) %>% 
                              group_by(var.X))
    
    res.reg <- data.frame(result.matrix.F %>% filter(slopes == s, n.replicates == r) %>% 
                            group_by(var.X))
    names(res.anova)[names(res.anova) == 'omega_squaredANOVA'] <- 'omega.squared'
    names(res.anova)[names(res.anova) == 'ANOVA'] <- 'F'
    res.anova.omega<-aggregate(omega.squared~var.X,res.anova, FUN = mean)
    colnames(res.anova.omega)<- c("var.X","averageOMEGA")
    res.anova.omega <- cbind(res.anova.omega,type=rep("anova",nrow(res.anova.omega)))
    names(res.reg)[names(res.reg) == 'omega_squaredRegression'] <- 'omega.squared'
    names(res.reg)[names(res.reg) == 'regression'] <- 'F'
    res.reg.omega<-aggregate(omega.squared~var.X,res.reg, FUN = mean)
    colnames(res.reg.omega)<- c("var.X","averageOMEGA")
    res.reg.omega <- cbind(res.reg.omega,type=rep("regression",nrow(res.reg.omega)))
    df <- data.frame(rbind(res.anova.omega,res.reg.omega))
    data.df[[count]] <- df
    omega.plots[[count]] <- local({
      count <- count
      sp <- ggplot(data=df, aes(colour=type, y=averageOMEGA,x=as.factor(var.X))) +
        geom_point(aes(colour=type),position="identity",stat="identity", size = 3) + 
        xlab("Variance") + ylab("Average ω^2") +
        scale_colour_manual("", values = c("anova" = "blue", "regression" = "firebrick")) +
        theme_classic() +
        theme(text = element_text(size = 12)) + 
        theme(axis.text = element_text(size = 12)) +
        ggtitle(bquote("n=" ~ .(n.replicates[count.n.replicates]) ~ " slope=" ~ .(sim.slopes[count.slopes]))) +
        theme(plot.title = element_text(hjust = 0.5)) + 
        guides(fill=guide_legend(title="  ")) +
        ylim(0,1)+
        theme(legend.text = element_text(size=14),
              legend.key.height = unit(0.5, 'cm'),
              legend.key.width = unit(0.5, 'cm'))
    })
    count <- count + 1
  }
}

OPLOT1<-ggarrange(omega.plots[[1]]+rremove("xlab"), 
                  omega.plots[[2]]+rremove("ylab")+rremove("xlab"),
                  omega.plots[[3]]+rremove("xlab")+rremove("xlab"),
                  omega.plots[[4]]+rremove("ylab")+rremove("xlab"),
                  omega.plots[[5]],
                  omega.plots[[6]]+rremove("ylab"),
                  ncol=2,nrow=3,common.legend = TRUE)
OPLOT1A<-annotate_figure(OPLOT1,
                         bottom =  text_grob("Number of Treatments=10", color = "black", face = "bold", size = 18))
View(result.matrix.F)

###########################################################################################
##Packages needed for heatmap
library(lattice)
library(viridisLite)
library(grid)
###############CODE TO PRODUCE HEAT MAP############THIS CODE USES A SIMILIAR, BUT DIFFERENT SIMULATION
##TO PRODUCE THE CODE#################################################################################

Final_results_per_sc<-matrix(NA,10,4,dimnames = list(c(),
                                                     c("p-value_ANOVA","p-value_regression","sd_p_ANOVA","sd_p_Regression")))

slopesloop<-seq(0.01,1,length.out = 10)
VI<-seq(2,10,length.out = 10)

Final_results_sc<-matrix(NA,100,4,dimnames = list(c(),
                                                  c("p-value_ANOVA_percent","p-value_regression_percent","slope","VI")))
for(j in 1:length(VI)){
  focal_VI <-VI[j]
  for(i in 1:length(slopesloop)){
    focal_slope<-slopesloop[i]
    test1<-lapply(1:100,function(x)create.design.linear(focal_slope,focal_VI,5,2))
    results_large<-as.data.frame(t(sapply(test1,function(x)x[[2]])))
    PercentageAnova<-length(which(results_large$anova<=0.05))/nrow(results_large)
    PercentageRegression<-length(which(results_large$regression<=0.05))/nrow(results_large)
    Final_results_sc[i +(j-1)*(length(VI)),c(1,2,3,4)]<-c(PercentageAnova,PercentageRegression,focal_slope,focal_VI)
    # mean_sd<-apply(results_large,2,function(x)c(mean(x),sd(x)))
    # Final_results_per_sc[i,c(1,2,3,4)]<-c(mean_sd[1,1],mean_sd[1,2],mean_sd[2,1],mean_sd[2,2])
    # }
  }}

Heatmap1<-levelplot(pvalue_ANOVA_percent ~ VI*slope, data=Heat_Map_Data,col.regions = inferno(100),
          xlab=list(label ="Variance", cex = 1.4),
          ylab =list(label="Slope", cex = 1.4),
          main="ANOVA",
          scales = list(x=list(cex=1.0),
                        y=list(cex=1.0)),
          colorkey = list(labels=list(cex=1.4),
                          title=list(label ="Power", cex = 1.4, rot = 270)))
Heatmap2<-levelplot(pvalue_Regression_percent ~ VI*slope, data=Heat_Map_Data,col.regions = inferno(100),
                    xlab=list(label ="Variance", cex = 1.4),
                    ylab =list(label="Slope", cex = 1.4),
                    main="Regression",
                    scales = list(x=list(cex=1.0),
                                  y=list(cex=1.0)),
                    colorkey = list(labels=list(cex=1.0),
                                    title=list(label ="Power", cex = 1.4, rot = 270)))

##################################################################################################

