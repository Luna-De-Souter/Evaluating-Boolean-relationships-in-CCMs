#rm(list = ls())

# Load libraries.
library(tidyr)
library(reshape2)
library(ggplot2)
library(viridis)
library(dplyr)

# Set directory and get dataframe (results of simulation experiment
# run in simulation_script.R).
setwd("discr_results")
setwd("C:/Users/lso055/R/CNA1/Contrapositives/discr_results")
df <- read.csv("disc_contr_19jan.csv")


###############################################
### remove results based on few observations ##
###############################################


# Set the number of runs per setting in the simulation experiment.
n_total <- 50000
# Calculate number of outputs containing at least one low contrapositive model.
df$n_low <- n_total - df$few_csfs - df$all_high_contr_AND
df$n_low
# Calculate number of outputs containing at least one high contrapositive model.
df$n_high <- n_total - df$few_csfs - df$all_low_contr_AND
df$n_high

sum(df$n_low < 20)
sum(df$n_high < 20)
# (Some settings have few outputs containing any high contrapositive models.)

# Replace relevant columns of settings for which fewer than 20 outputs contain any 
# high contrapositive models.
df[df$n_high < 20,c(9,11,13,15,27,37)] <- rep(NA,length(c(9,11,13,15,27,37)))
df[df$n_high < 20,]


################################
### standard errors of means ###
################################


## correctness

# For each combination of settings, calculate the standard error of correctness
# of high contrapositive models.
df$SEM_high_corr <- df$high_corr_sd/sqrt(df$n_high)
df$SEM_high_corr
# Get the minimum and maximum standard error.
min(df$SEM_high_corr, na.rm = T)
max(df$SEM_high_corr, na.rm = T)

# for each combination of settings, calculate the standard error of correctness
# of low contrapositive models.
df$SEM_low_corr <- df$low_corr_sd/sqrt(df$n_low)
df$SEM_low_corr
# Get the minimum and maximum standard error.
min(df$SEM_low_corr, na.rm = T)
max(df$SEM_low_corr, na.rm = T)


## degree of completeness

# Calculate the standard errors of degree of completeness means of high
# contrapositive models.
df$SEM_high_comp <- df$high_complx_sd/sqrt(df$n_high_contr_complx)
df$SEM_high_comp
# Get the minimum and maximum standard error.
min(df$SEM_high_comp, na.rm = T)
max(df$SEM_high_comp, na.rm = T)

# Calculate the standard errors of degree of completeness means of low
# contrapositive models.
df$SEM_low_comp <- df$low_complx_sd/sqrt(df$n_low_contr_complx)
df$SEM_low_comp
# Get the minimum and maximum standard error.
min(df$SEM_low_comp)
max(df$SEM_low_comp)


#########################
### correctness plots ###
#########################


selection1 <- df[,c(2,3,6,8,9)]
selection1 <- selection1 %>% group_by(con_threshold, cov_threshold, prevalence,
                                      low_contr_AND_correct,
                                      high_contr_AND_correct
)
selection1

k1 <- melt(selection1, id.vars = c("prevalence", "con_threshold","cov_threshold"))

k1$lab <- rep(c(rep("cov = 0.7",4),rep("cov = 0.75",4),rep("cov = 0.8",4),rep("cov = 0.85",4)), 14)
k1$lab <- factor(k1$lab, levels = c("cov = 0.7", "cov = 0.75", "cov = 0.8", "cov = 0.85"))

k1$variable <- c(rep("low contrapositive", 112), rep("high contrapositive", 112))
k1$variable <- factor(k1$variable, levels = c("low contrapositive", "high contrapositive"))

k1$prevalence <- k1$prevalence*100
k1$prevalence <- as.character(k1$prevalence)
k1$prevalence <- paste0(k1$prevalence, "%")
k1$prevalence <- factor(k1$prevalence, levels = c("60%", "65%", "70%", "75%", "80%", "85%", "90%"))

k1 <- k1 %>% tidyr:::unite("BothLabels", lab , con_threshold , sep = ", con= ", remove = FALSE)
k1$BothLabels <- factor(k1$BothLabels, levels = unique(k1$BothLabels))

plot1 <- ggplot(k1, aes(x = prevalence, y = value, fill = variable)) +
  geom_bar(stat="identity", position = position_dodge(width = 0.5), width = .4) +
  scale_fill_manual(values = c('#303030','#C0C0C0')) +
  theme_bw() + facet_wrap(vars(BothLabels)) +
  theme(plot.title = element_text(size = 3), axis.text.x = element_text(size = 10),
        axis.title = element_text(size = 14),
        legend.spacing.x = unit(0.5, 'cm'),
        legend.title = element_blank() , strip.text = element_text(size = 12),
        legend.text = element_text(size =14,margin = margin(l = -8, r =10, unit = "pt")))+
  scale_x_discrete(name ="prevalence")+ ylim (0, 1) + theme(legend.position ="top") +
  ylab ("correctness")
plot1


####################################
### degree of completeness plots ###
####################################

selection1 <- df[,c(2,3,6,10,11)]
selection1 <- selection1 %>% group_by(con_threshold, cov_threshold, prevalence,
                                      low_contr_AND_complx,
                                      high_contr_AND_complx
)
selection1

k1 <- melt(selection1, id.vars = c("prevalence", "con_threshold","cov_threshold"))

k1$lab <- rep(c(rep("cov = 0.7",4),rep("cov = 0.75",4),rep("cov = 0.8",4),rep("cov = 0.85",4)), 14)
k1$lab <- factor(k1$lab, levels = c("cov = 0.7", "cov = 0.75", "cov = 0.8", "cov = 0.85"))

k1$variable <- c(rep("low contrapositive", 112), rep("high contrapositive", 112))
k1$variable <- factor(k1$variable, levels = c("low contrapositive", "high contrapositive"))

k1$prevalence <- k1$prevalence*100
k1$prevalence <- as.character(k1$prevalence)
k1$prevalence <- paste0(k1$prevalence, "%")
k1$prevalence <- factor(k1$prevalence, levels = c("60%", "65%", "70%", "75%", "80%", "85%", "90%"))

k1 <- k1 %>% tidyr:::unite("BothLabels", lab , con_threshold , sep = ", con= ", remove = FALSE)
k1$BothLabels <- factor(k1$BothLabels, levels = unique(k1$BothLabels))

plot1 <- ggplot(k1, aes(x = prevalence, y = value, fill = variable)) +
  geom_bar(stat="identity", position = position_dodge(width = 0.5), width = .4) +
  scale_fill_manual(values = c('#303030','#C0C0C0')) +
  theme_bw() + facet_wrap(vars(BothLabels)) +
  theme(plot.title = element_text(size = 3), axis.text.x = element_text(size = 10),
        axis.title = element_text(size = 14),
        legend.spacing.x = unit(0.5, 'cm'),
        legend.title = element_blank() , strip.text = element_text(size = 12),
        legend.text = element_text(size =14,margin = margin(l = -8, r =10, unit = "pt")))+
  scale_x_discrete(name ="prevalence")+ ylim (0, 1) + theme(legend.position ="top") +
  ylab ("degree of completeness")
plot1


##############################
### model complexity plots ###
##############################


selection1 <- df[,c(2,3,6,14,15)]
selection1 <- selection1 %>% group_by(con_threshold, cov_threshold, prevalence,
                                      low_contr_AND_nct,
                                      high_contr_AND_nct
)
selection1

k1 <- melt(selection1, id.vars = c("prevalence", "con_threshold","cov_threshold"))

k1$lab <- rep(c(rep("cov = 0.7",4),rep("cov = 0.75",4),rep("cov = 0.8",4),rep("cov = 0.85",4)), 14)
k1$lab <- factor(k1$lab, levels = c("cov = 0.7", "cov = 0.75", "cov = 0.8", "cov = 0.85"))

k1$variable <- c(rep("low contrapositive", 112), rep("high contrapositive", 112))
k1$variable <- factor(k1$variable, levels = c("low contrapositive", "high contrapositive"))

k1$prevalence <- k1$prevalence*100
k1$prevalence <- as.character(k1$prevalence)
k1$prevalence <- paste0(k1$prevalence, "%")
k1$prevalence <- factor(k1$prevalence, levels = c("60%", "65%", "70%", "75%", "80%", "85%", "90%"))

k1 <- k1 %>% tidyr:::unite("BothLabels", lab , con_threshold , sep = ", con= ", remove = FALSE)
k1$BothLabels <- factor(k1$BothLabels, levels = unique(k1$BothLabels))

plot1 <- ggplot(k1, aes(x = prevalence, y = value, fill = variable)) +
  geom_bar(stat="identity", position = position_dodge(width = 0.5), width = .4) +
  scale_fill_manual(values = c('#303030','#C0C0C0')) +
  theme_bw() + facet_wrap(vars(BothLabels)) +
  theme(plot.title = element_text(size = 3), axis.text.x = element_text(size = 10),
        axis.title = element_text(size = 14),
        legend.spacing.x = unit(0.5, 'cm'),
        legend.title = element_blank() , strip.text = element_text(size = 12),
        legend.text = element_text(size =14,margin = margin(l = -8, r =10, unit = "pt")))+
  scale_x_discrete(name ="prevalence")+ ylim (0, 10) + theme(legend.position ="top") +
  ylab ("complexity")
plot1


###############################
### number of outputs plots ###
###############################

selection1 <- df[,c(2,3,6,39,40)]
selection1 <- selection1 %>% group_by(con_threshold, cov_threshold, prevalence,
                                      n_low, n_high)
selection1

k1 <- melt(selection1, id.vars = c("prevalence", "con_threshold","cov_threshold"))

k1$lab <- rep(c(rep("cov = 0.7",4),rep("cov = 0.75",4),rep("cov = 0.8",4),rep("cov = 0.85",4)), 14)
k1$lab <- factor(k1$lab, levels = c("cov = 0.7", "cov = 0.75", "cov = 0.8", "cov = 0.85"))

k1$variable <- c(rep("low contrapositive", 112), rep("high contrapositive", 112))
k1$variable <- factor(k1$variable, levels = c("low contrapositive", "high contrapositive"))

k1$prevalence <- k1$prevalence*100
k1$prevalence <- as.character(k1$prevalence)
k1$prevalence <- paste0(k1$prevalence, "%")
k1$prevalence <- factor(k1$prevalence, levels = c("60%", "65%", "70%", "75%", "80%", "85%", "90%"))

k1 <- k1 %>% tidyr:::unite("BothLabels", lab , con_threshold , sep = ", con= ", remove = FALSE)
k1$BothLabels <- factor(k1$BothLabels, levels = unique(k1$BothLabels))

plot1 <- ggplot(k1, aes(x = prevalence, y = value, fill = variable)) +
  geom_bar(stat="identity", position = position_dodge(width = 0.5), width = .4) +
  scale_fill_manual(values = c('#303030','#C0C0C0')) +
  theme_bw() + facet_wrap(vars(BothLabels)) +
  theme(plot.title = element_text(size = 3), axis.text.x = element_text(size = 10),
        axis.title = element_text(size = 14),
        legend.spacing.x = unit(0.5, 'cm'),
        legend.title = element_blank() , strip.text = element_text(size = 12),
        legend.text = element_text(size =14,margin = margin(l = -8, r =10, unit = "pt")))+
  scale_x_discrete(name ="prevalence")+ ylim (0, 50000) + theme(legend.position ="top") +
  ylab ("number of outputs")
plot1

