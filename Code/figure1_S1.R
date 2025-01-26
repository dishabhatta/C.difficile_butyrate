## Paper: SCFA and C. difficile
## Figure 1: Growth for all scfa prelim data
### 1a: growth curves for CD630
### Supplementary figure S1: growth curves for VPI10463, R20291
### 1b: cd630 colony counts over 6, 12, 18, 24h
### 1c: growth curves for butyrate concentration gradient
### 1d: cd630 colony counts over time vs scfa concentration gradient

rm(list = ls())
setwd("/Data/figure1_S1/")

## Library

library(readxl)
library(tidyverse)
library(tidyr)
library(ggplot2)
library(nlme)
library(stringi)
library(stringr)
library(plyr)
library(dplyr)
library(growthcurver)
library(DescTools)




#### Fig1a: growth curves for CD630
#### This creates the master file
#file_list <- c("SCFACU21", "SCFACU22", "SCFACU29", "SCFACU30", "SCFACU33", "SCFACU34", "SCFACU35")
#filepath <- "~/wherever_the_files_are_stored/"
#list.data <- list()
#for(i in (1:length(file_list))){
  #print(i)
#  list.data[[i]] <- read_tsv(paste0(filepath, file_list[i],'_data.txt'))
  
#}
#list.data.meta <- list()
#for(i in (1:length(file_list))){
  #print(i)
#  list.data.meta[[i]] <- read_tsv(paste0(filepath, file_list[i],'_meta.txt'))
#}

#full.merge.630 <- data.frame()
#time_vec <- list.data[[1]][1]

#for (i in 3:7) {
#  list.data[[i]] <- list.data[[i]][-97,]
#  list.data[[i]][,1] <- time_vec$time
#}

#for (i in 1:7) {
#  list_2 <- list.data[[i]][,-2]
#  colnames(list_2)[1] <- c("time")
#  list_2 <- mutate(list_2, time=as.numeric(gsub("s", "", time)))
#  list_3 <- as.data.frame(t(list_2))
#  list_3$wellID <- rownames(list_3)
#  list_3[1,97] <- "wellID"
#  colnames(list_3) <- list_3[1,]
#  list_4 <- list_3[-1,]
#  merge_630 <- inner_join(list.data.meta[[i]], list_4, by= "wellID")
#  #merge_cdmm <- subset(merge_cdmm, !(merge_cdmm$condition == "Cdiff"| merge_cdmm$condition =="noCdiff" | merge_cdmm$condition == "na"))
#  full.merge.630 <- rbind(full.merge.630, merge_630)
#}

#levels(factor(full.merge.630$expID))
#write_tsv(full.merge.630, file = "/Data/figure1_S1/full.merge.630.tsv")

full.merge.630 <- read_tsv(file = "/Data/figure1_S1/full.merge.630.tsv")
full.merge.630 <- subset(full.merge.630, !(full.merge.630$group =="noCdiff"))
full.merge.630.avg <- full.merge.630 %>% 
  group_by(expID,group) %>% 
  summarise_all(mean) %>% 
  as.data.frame()

#levels(factor(full_merge$group))
#[1] "Ace25" "Ace5"  "BHI"   "But25" "But5"  "Pro25" "Pro5" 
#levels(factor(full_merge$expID))
#[1] "SCFACU21" "SCFACU22" "SCFACU29" "SCFACU30" "SCFACU33" "SCFACU34" "SCFACU35"

## Input the control plate and modify to be able to subtract the controls

controls <- read_tsv("/Data/figure1_S1/controls_new.txt")
controls_meta <- read_tsv("/Data/figure1_S1/controls_new_meta.txt")
controls_GC <- subset(controls_meta, controls_meta$medium == "BHI")

controls$time <- time_vec$time # essentially the time is different every time the machine is run so copy and paste the same time otherwise it wont merge, the difference in time is a second before and after just not exact 
controls_2 <- controls[,-2]
controls_2 <- mutate(controls_2, time = as.numeric(gsub("s", "", time)))
controls_3 <- as.data.frame(t(controls_2))
controls_3$wellID <- rownames(controls_3)
controls_3[1,97] <- "wellID"
colnames(controls_3) <- controls_3[1,]
controls_4 <- controls_3[-1,]

merge_controls <- inner_join(controls_GC, controls_4, by = "wellID")

## Create avg of the controls and subtract from the average of the 3 expts
norm.mat <- merge_controls %>% 
  group_by(expID, group) %>% 
  summarise_all(mean) %>% 
  as.data.frame() %>% 
  select(10:105) %>% 
  as.matrix()

normalized_data <- sweep(as.matrix(full.merge.630.avg[,11:106]), 2, norm.mat)

combined <- select(full.merge.630.avg, 1:2) %>% 
  cbind(., normalized_data)

combined2 <- combined %>% pivot_longer(cols=!(1:2), names_to = "time", values_to = "OD")

final_630 <- combined2 %>% 
  #filter(expID %in% c("SCFACU21", "SCFACU22", "SCFACU29")) %>%
  filter(group %in% c("Ace25", "Ace5",  "BHI", "But25", "But5",  "Pro25", "Pro5")) %>%
  droplevels() %>%
  group_by(time, group) %>%
  dplyr::summarise(mean = mean(OD, na.rm = TRUE), n = length(OD), sd = sd(OD, na.rm = TRUE), se = sd / sqrt(n))

final_630 <- mutate(final_630, time = as.numeric(as.numeric(time)/3600))
final_630$group <- factor(final_630$group, levels = c("BHI", "Ace5", "Ace25", "Pro5", "Pro25", "But5", "But25"))

ggplot(final_630, aes(x=time, y=(mean), colour=group)) +
  geom_smooth(linewidth = 0.5, se = FALSE) + 
  #geom_point(size=.1, alpha = 0.5) +
  #geom_errorbar(aes(ymin=(mean-se), ymax=(mean+se)), width=.01, alpha =0.5) +
  theme(legend.position = "right") +
  scale_colour_manual(values=c("black", "#801100", "#FC6400","#5D78AB","#87A6BD","#145425", "#80987c")) +
  ylab("OD600") +
  xlab("Time(hrs)") +
  theme_minimal() +
  scale_x_continuous(breaks = seq(0, 24, by = 4)) +
  scale_y_log10(breaks=c(0.01,0.03,0.10,0.3)) +
  #scale_y_continuous(trans = log_trans()) +
  theme(legend.key=element_rect(fill='gray99')) +
  theme(axis.title.x=element_text(size=12, face='bold'),
        axis.title.y=element_text(size=12, face='bold'),
        panel.background = element_rect(fill = 'gray99'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(fill=NA, colour = NA))

##STATS

## to calculate area under the curve using growthcurver
#combo <- combined
#combo2 <- combo[,-1]
#combo3 <- as.data.frame(t(combo2))
#combo3$time <- rownames(combo3)
#write_tsv(combo3, "/Data/figure1_S1/combo3_gcvr.txt")
# edits were made to put triplicates together in Excel

combo4_edit_x <- read_tsv(file= "/Data/figure1_S1/combo3_gcvr_edit.txt")
gcvr_out <- SummarizeGrowthByPlate(combo4_edit_x)
gcvr_out$group[gcvr_out$sample %in% c("Ace25...2","Ace25...9","Ace25...16", "Ace25...23", "Ace25...30", "Ace25...37", "Ace25...44")] <- "Ace25"
gcvr_out$group[gcvr_out$sample %in% c("Ace5...3", "Ace5...10", "Ace5...17", "Ace5...24", "Ace5...31", "Ace5...38", "Ace5...45")] <- "Ace5"
gcvr_out$group[gcvr_out$sample %in% c("BHI...4", "BHI...11", "BHI...18", "BHI...25", "BHI...32", "BHI...39", "BHI...46")] <- "BHI"
gcvr_out$group[gcvr_out$sample %in% c("But25...5", "But25...12", "But25...19", "But25...26", "But25...33", "But25...40", "But25...47")] <- "But25"
gcvr_out$group[gcvr_out$sample %in% c("But5...6", "But5...13", "But5...20", "But5...27", "But5...34", "But5...41", "But5...48")] <- "But5"
gcvr_out$group[gcvr_out$sample %in% c("Pro25...7", "Pro25...14", "Pro25...21", "Pro25...28", "Pro25...35", "Pro25...42", "Pro25...49")] <- "Pro25"
gcvr_out$group[gcvr_out$sample %in% c("Pro5...8", "Pro5...15", "Pro5...22", "Pro5...29", "Pro5...36", "Pro5...43", "Pro5...50")] <- "Pro5"

aov_model <- aov(auc_l ~ group, data = gcvr_out)

summary(aov_model)

#           Df   Sum Sq  Mean Sq F value   Pr(>F)    
# group        6 382507395 63751232   30.17 1.01e-13 ***
# Residuals   42  88736832  2112782                      
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

TukeyHSD(aov_model)

#Tukey multiple comparisons of means
#95% family-wise confidence level

#Fit: aov(formula = auc_l ~ group, data = gcvr_out)

#$group
#                    diff         lwr        upr     p adj
#Ace5-Ace25    685.2684 -1719.817  3090.35377 0.9733519
#BHI-Ace25    -342.2521 -2747.337  2062.83327 0.9993787
#But25-Ace25 -6829.8585 -9234.944 -4424.77317 0.0000000
#But5-Ace25  -3193.8543 -5598.940  -788.76892 0.0031519
#Pro25-Ace25 -6308.2745 -8713.360 -3903.18917 0.0000000
#Pro5-Ace25  -2500.8420 -4905.927   -95.75669 0.0368552
#BHI-Ace5    -1027.5205 -3432.606  1377.56483 0.8375676
#But25-Ace5  -7515.1269 -9920.212 -5110.04160 0.0000000
#But5-Ace5   -3879.1227 -6284.208 -1474.03736 0.0002077
#Pro25-Ace5  -6993.5429 -9398.628 -4588.45760 0.0000000
#Pro5-Ace5   -3186.1105 -5591.196  -781.02513 0.0032461
#But25-BHI   -6487.6064 -8892.692 -4082.52110 0.0000000 ***
#But5-BHI    -2851.6022 -5256.688  -446.51686 0.0111465 *
#Pro25-BHI   -5966.0224 -8371.108 -3560.93710 0.0000000 ***
#Pro5-BHI    -2158.5900 -4563.675   246.49537 0.1042441
#But5-But25   3636.0042  1230.919  6041.08957 0.0005574
#Pro25-But25   521.5840 -1883.501  2926.66933 0.9935122
#Pro5-But25   4329.0165  1923.931  6734.10180 0.0000321
#Pro25-But5  -3114.4202 -5519.506  -709.33491 0.0042556
#Pro5-But5     693.0122 -1712.073  3098.09756 0.9718312
#Pro5-Pro25   3807.4325  1402.347  6212.51780 0.0002785


pairwise.t.test(x=gcvr_out$auc_l, g=gcvr_out$group, p.adjust.method = "bonferroni")

#Pairwise comparisons using t tests with pooled SD 

#data:  gcvr_out$auc_l and gcvr_out$group 

#          Ace25  Ace5    BHI    But25  But5   Pro25 
#  Ace5  1.00000 -       -       -       -       -      
#  BHI   1.00000 1.00000 -       -       -       -      
#  But25 9.5e-10 6.3e-11 3.8e-09 -       -       -      
#  But5  0.00376 0.00023 0.01424 0.00063 -       -      
#  Pro25 8.0e-09 4.9e-10 3.3e-08 1.00000 0.00515 -      
#  Pro5  0.05215 0.00388 0.17085 3.4e-05 1.00000 0.00031

#P value adjustment method: bonferroni 

DunnettTest(x=gcvr_out$auc_l, g=gcvr_out$group, control = "BHI")

#Dunnett's test for comparing several treatments with a control :  
#    95% family-wise confidence level

#$BHI
#               diff    lwr.ci     upr.ci    pval    
#Ace25-BHI   342.2521 -1736.682  2421.18644  0.9955    
#Ace5-BHI   1027.5205 -1051.414  3106.45488  0.5996    
#But25-BHI -6487.6064 -8566.541 -4408.67205 1.2e-10 ***
#But5-BHI  -2851.6022 -4930.537  -772.66781  0.0036 ** 
#Pro25-BHI -5966.0224 -8044.957 -3887.08805 4.1e-09 ***
#Pro5-BHI  -2158.5900 -4237.524   -79.65558  0.0389 *  

#---
#Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1


GC_630 <- lme(as.numeric(time)~1+ group, data=full_merge_2, random = ~1|OD, na.action = na.omit)
#GC_630 <- lme(mean~group, data=final, random = ~1|time, na.action = na.omit)
#summary(GC_630)




#####################
## Fig1b: cd630 colony counts over 24h for all scfa

### input

cfu_630_24 <- read_excel(path = "/Data/figure1_S1/new_fig1B_data_630.xlsx")
cfu_630_24 <- cfu_630_24[,-5]
cfu_630_24_5x <- subset(cfu_630_24, !(cfu_630_24$expID %in% c("SCFACU19", "SCFACU22")))
cfu_630_24_3x <- subset(cfu_630_24, !(cfu_630_24$expID %in% c("SCFACU19", "SCFACU22", "SCFACU29")))

cfus_630_24.avg <- cfu_630_24 %>% 
  #filter(expID %in% c("SCFACU21", "SCFACU22", "SCFACU29")) %>%
  filter(group %in% c("Ace25", "Ace5",  "BHI", "But25", "But5",  "Pro25", "Pro5")) %>%
  droplevels() %>%
  group_by(time, group) %>%
  dplyr::summarise(mean = mean(CFU, na.rm = TRUE), n = length(CFU), sd = sd(CFU, na.rm = TRUE), se = sd / sqrt(n))

cfus_630_24.avg$group <- factor(cfus_630_24.avg$group, levels = c("BHI", "Ace5", "Ace25", "Pro5", "Pro25", "But5", "But25"))
cfu_630_24$group <- factor(cfu_630_24$group, levels = c("BHI", "Ace5", "Ace25", "Pro5", "Pro25", "But5", "But25"))

ggplot(cfu_630_24, aes(x=group, y=CFU, color=group)) + geom_point(size =0.5, position = position_jitter(width = 0.2)) +
  #geom_point(size =0.5) + 
  #facet_grid(~time) +
  geom_boxplot(alpha =0.4, size = 0.3, aes(fill = group), outlier.shape = NA) +
  theme_bw() +
  scale_y_continuous(limits= c(0, 10), breaks = seq(0,10,2))+
  #geom_errorbar(aes(ymin=(mean-se), ymax=(mean+se)), width=.01, alpha =0.5) +
  theme(legend.position = "right", panel.grid = element_blank()) +
  scale_color_manual(values=c("black","#801100", "#FC6400", "#5D78AB","#87A6BD","#145425", "#80987c")) +
  scale_fill_manual(values=c("black","#801100", "#FC6400", "#5D78AB","#87A6BD","#145425", "#80987c"))

#STATS

aov_cfu_24 <- aov(CFU~group, data=cfu_630_24)
summary(aov_cfu_24)

#Df Sum Sq Mean Sq F value   Pr(>F)    
#group         6  16.28  2.7135    7.07 1.36e-06 ***
#  Residuals   140  53.73  0.3838                     
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

DunnettTest(x=cfu_630_24$CFU, g=cfu_630_24$group, "BHI")

#Dunnett's test for comparing several treatments with a control :  
#    95% family-wise confidence level

#$BHI
#                 diff      lwr.ci     upr.ci   pval    
#Ace5-BHI   0.28522635 -0.21167321  0.7821259 0.4741    
#Ace25-BHI  0.44952052 -0.04737904  0.9464201 0.0920 .  
#Pro5-BHI   0.02791307 -0.46898649  0.5248126 1.0000    
#Pro25-BHI  0.12676441 -0.37013515  0.6236640 0.9666    
#But5-BHI  -0.29905346 -0.79595302  0.1978461 0.4257    
#But25-BHI -0.62106026 -1.11795982 -0.1241607 0.0078 ** 

#---
#Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1




###########################################

## Fig 1c: growth curves for butyrate concentration gradient


## This will make the final data table for use:
#expt1_dr <- read_tsv("./SCFACU36_data.txt")
#colnames(expt1_dr)[c(1,2)] <- c("time", "temp")
#expt1_meta_dr <- read_tsv("./SCFACU36_meta.txt")
#expt2_dr <- read_tsv("./SCFACU37_data.txt")
#colnames(expt2_dr)[c(1,2)] <- c("time", "temp")
#expt2_meta_dr <- read_tsv("./SCFACU37_meta.txt")
#expt3_dr <- read_tsv("./SCFACU38_data.txt")
#colnames(expt3_dr)[c(1,2)] <- c("time", "temp")
#expt3_meta_dr <- read_tsv("./SCFACU38_meta.txt")

### For experiment1
#expt1_dr_2 <- expt1_dr[,-2]
#expt1_dr_2 <- mutate(expt1_dr_2, time = as.numeric(gsub("s", "", time)))
#expt1_dr_3 <- as.data.frame(t(expt1_dr_2))
#expt1_dr_3$wellID <- rownames(expt1_dr_3)
#expt1_dr_3[1,97] <- "wellID"
#colnames(expt1_dr_3) <- expt1_dr_3[1,]
#expt1_dr_4 <- expt1_dr_3[-1,]
#expt1_dr_4 <- expt1_dr_4[,-97]
#colnames(expt1_dr_4)[97] <- "wellID"
#merge_expt1_dr <- inner_join(expt1_meta_dr, expt1_dr_4, by = "wellID")
#merge_expt1_dr <- merge_expt1_dr[-(49:72),]

### For experiment2
#expt2_dr$time <- expt1_dr$time
#expt2_dr_2 <- expt2_dr[,-2]
#expt2_dr_2 <- mutate(expt2_dr_2, time = as.numeric(gsub("s", "", time)))
#expt2_dr_3 <- as.data.frame(t(expt2_dr_2))
#expt2_dr_3$wellID <- rownames(expt2_dr_3)
#expt2_dr_3[1,97] <- "wellID"
#colnames(expt2_dr_3) <- expt2_dr_3[1,]
#expt2_dr_4 <- expt2_dr_3[-1,]
#expt2_dr_4 <- expt2_dr_4[,-97]
#colnames(expt2_dr_4)[97] <- "wellID"

#merge_expt2_dr <- inner_join(expt2_meta_dr, expt2_dr_4, by = "wellID")
#merge_expt2_dr <- merge_expt2_dr[-(49:72),]

### For experiment3
#expt3_dr$time <- expt1_dr$time
#expt3_dr_2 <- expt3_dr[,-2]
#expt3_dr_2 <- mutate(expt3_dr_2, time = as.numeric(gsub("s", "", time)))
#expt3_dr_3 <- as.data.frame(t(expt3_dr_2))
#expt3_dr_3$wellID <- rownames(expt3_dr_3)
#expt3_dr_3[1,97] <- "wellID"
#colnames(expt3_dr_3) <- expt3_dr_3[1,]
#expt3_dr_4 <- expt3_dr_3[-1,]
#expt3_dr_4 <- expt3_dr_4[,-97]
#colnames(expt3_dr_4)[97] <- "wellID"
#merge_expt3_dr <- inner_join(expt3_meta_dr, expt3_dr_4, by = "wellID")
#merge_expt3_dr <- merge_expt3_dr[-(49:72),]

### All experiments are combined and averaged now

#full_merge_dr <- rbind(merge_expt1_dr, merge_expt2_dr, merge_expt3_dr)
#write_tsv(full_merge_dr, file = "/Data/figure1_S1/full.merge.dr.tsv")
full_merge_dr <- read_tsv(file = "/Data/figure1_S1/full.merge.dr.tsv")
full_merge_dr_2 <- full_merge_dr %>% pivot_longer(cols = !(1:10), names_to = "time", values_to = "OD")

full_merge_avg_dr <- full_merge_dr %>% 
  group_by(group,expID) %>% 
  summarise_all(mean) %>% 
  as.data.frame()

norm.mat_dr <- full_merge_dr %>% 
  group_by(group, expID) %>% 
  summarise_all(mean) %>% 
  as.data.frame() %>% 
  filter(., group == "noCdiff") %>% 
  select(11:106) %>% 
  as.matrix()

## normalize
normalized_data_dr <- sweep(as.matrix(full_merge_avg_dr[,11:106]), 2, norm.mat_dr)

## combine it with metadata
formatted_dr <- select(full_merge_avg_dr, 1:2) %>% 
  cbind(., normalized_data_dr)

## pivot longer
combined_dr <- formatted_dr %>% pivot_longer(cols=!(1:2), names_to = "time", values_to = "OD")

combined_dr <- mutate(combined_dr, time = as.numeric(as.numeric(time)/3600))
final_dr <- combined_dr %>% 
  #filter(expID %in% c("SCFACU21", "SCFACU22", "SCFACU29")) %>%
  filter(group %in% c("BHI","But10","But25","But5","But50")) %>%
  droplevels() %>%
  group_by(time, group) %>%
  dplyr::summarise(mean = mean(OD, na.rm = TRUE), n = length(OD), sd = sd(OD, na.rm = TRUE), se = sd / sqrt(n))

final_dr <- mutate(final_dr, time = as.numeric(as.numeric(time)/3600))

ggplot(final_dr, aes(x=time, y=(mean), colour=group)) +
  geom_smooth(size = 1, se = FALSE) + 
  #geom_point(size=.1, alpha = 0.9) +
  #geom_errorbar(aes(ymin=(mean-se), ymax=(mean+se)), width=.01, alpha =0.5) +
  theme(legend.position = "right") +
  scale_colour_manual(values=c("black", hcl.colors(4, palette = "TealGrn"))) +
  ylab("OD600") +
  xlab("Time(hrs)") +
  theme_minimal() +
  scale_y_log10(breaks = c(0.01,0.03,0.1,0.3), limits = c(0.01, 0.7)) +
  scale_x_continuous(breaks = seq(0, 24, by = 4)) +
  theme(legend.key=element_rect(fill='gray99')) +
  theme(axis.title.x=element_text(size=12, face='bold'),
        axis.title.y=element_text(size=12, face='bold'),
        panel.background = element_rect(fill = 'gray99'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(fill=NA, colour = NA))

##STATS

#combo <- formatted_dr
#combo2 <- combo[,-2]
#combo3 <- as.data.frame(t(combo2))
#combo3$time <- rownames(combo3)
#write_tsv(combo3, "./combo3_gcvr_dr.txt")
# edits were made to put triplicates together in Excel
combo4_edit_x <- read_tsv(file= "/Data/figure1_S1/combo3_gcvr_dr_edit.txt")
gcvr_out_dr <- SummarizeGrowthByPlate(combo4_edit_x)

gcvr_out_dr$group[gcvr_out_dr$sample == "BHI...2" | gcvr_out_dr$sample == "BHI...3" | gcvr_out_dr$sample == "BHI...4"] <- "BHI"
gcvr_out_dr$group[gcvr_out_dr$sample == "But10...5" | gcvr_out_dr$sample == "But10...6" | gcvr_out_dr$sample == "But10...7"] <- "But10"
gcvr_out_dr$group[gcvr_out_dr$sample == "But25...10" | gcvr_out_dr$sample == "But25...9" | gcvr_out_dr$sample == "But25...8"] <- "But25"
gcvr_out_dr$group[gcvr_out_dr$sample == "But5...11" | gcvr_out_dr$sample == "But5...12" | gcvr_out_dr$sample == "But5...13"] <- "But5"
gcvr_out_dr$group[gcvr_out_dr$sample == "But50...14" | gcvr_out_dr$sample == "But50...15" | gcvr_out_dr$sample == "But50...16"] <- "But50"

aov_model <- aov(auc_l ~ group, data = gcvr_out_dr)

summary(aov_model)

#Df   Sum Sq  Mean Sq F value   Pr(>F)    
#group        4 85837113 21459278   13.59 0.000473 ***
#  Residuals   10 15791173  1579117                     
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

TukeyHSD(aov_model)

#Tukey multiple comparisons of means
#95% family-wise confidence level

#Fit: aov(formula = auc_l ~ group, data = gcvr_out_dr)

#$group
#                   diff         lwr        upr     p adj
#But10-BHI   -4461.23034  -7837.9911 -1084.4696 0.0098635
#But25-BHI   -4504.78637  -7881.5471 -1128.0256 0.0092562
#But5-BHI    -2127.51935  -5504.2801  1249.2414 0.3014671
#But50-BHI   -7053.71571 -10430.4765 -3676.9549 0.0003230
#But25-But10   -43.55604  -3420.3168  3333.2047 0.9999991
#But5-But10   2333.71098  -1043.0498  5710.4718 0.2296184
#But50-But10 -2592.48538  -5969.2461   784.2754 0.1600546
#But5-But25   2377.26702   -999.4937  5754.0278 0.2163729
#But50-But25 -2548.92934  -5925.6901   827.8314 0.1702869
#But50-But5  -4926.19636  -8302.9571 -1549.4356 0.0050535

DunnettTest(x=gcvr_out_dr$auc_l, g=gcvr_out_dr$group, "BHI")

#Dunnett's test for comparing several treatments with a control :  
#    95% family-wise confidence level

#$BHI
#               diff     lwr.ci     upr.ci    pval    
#But10-BHI -4461.230  -7430.142 -1492.3190 0.00476 ** 
#But25-BHI -4504.786  -7473.698 -1535.8751 0.00473 ** 
#But5-BHI  -2127.519  -5096.431   841.3919 0.18367    
#But50-BHI -7053.716 -10022.627 -4084.8044 0.00034 ***

#---
#Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1







########################

# fig 1d: cfus different scfa over time

all_cfus <- read_excel(path = "/Data/figure1_S1/full_cfu_3_expt.xlsx")
all_cfus <- all_cfus[,-5]
added_data <- read_xlsx("/Data/figure1_S1/full_fig2_cfus.xlsx")
added_data <- added_data[,-5]

all_cfus2 <- rbind(all_cfus, added_data)
all_cfus3 <- subset(all_cfus2, all_cfus2$type == "unheated")
all_cfus_new <- all_cfus3 %>% 
  filter(group %in% c("Ace25", "Ace5",  "BHI", "But25", "But5",  "Pro25", "Pro5")) %>%
  droplevels()

all_cfus_630 <- all_cfus_new %>% 
  #filter(expID %in% c("SCFACU21", "SCFACU22", "SCFACU29")) %>%
  filter(group %in% c("Ace25", "Ace5",  "BHI", "But25", "But5",  "Pro25", "Pro5")) %>%
  droplevels() %>%
  group_by(time, group) %>%
  dplyr::summarise(mean = mean(CFU, na.rm = TRUE), n = length(CFU), sd = sd(CFU, na.rm = TRUE), se = sd / sqrt(n))

all_cfus_new$group <- factor(all_cfus_new$group, levels = c("BHI", "Ace5", "Ace25", "Pro5", "Pro25", "But5", "But25"))

ggplot(all_cfus_new, aes(x=as.factor(time), y=CFU, color=group)) + geom_point(size =0.5, position = position_jitterdodge(jitter.width = 0.5)) +
  #geom_point(size =0.5) + 
  #facet_grid(~time) +
  geom_boxplot(alpha =0.4, size = 0.3, aes(fill = group), outlier.shape = NA) +
  theme_bw() +
  scale_y_continuous(limits= c(0, 10), breaks = seq(0,10,2))+
  #geom_errorbar(aes(ymin=(mean-se), ymax=(mean+se)), width=.01, alpha =0.5) +
  theme(legend.position = "right") +
  scale_color_manual(values=c("black","#801100", "#FC6400", "#5D78AB","#87A6BD","#145425", "#80987c")) +
  scale_fill_manual(values=c("black","#801100", "#FC6400", "#5D78AB","#87A6BD","#145425", "#80987c"))
#scale_x_continuous(breaks = seq(0, 24, by = 6))

##STATS

all_cfus_6 <- subset(all_cfus_new, all_cfus_new$time == 6)
aov_cfu_6 <- aov(CFU~group, data = all_cfus_6)
summary(aov_cfu_6)

#Df Sum Sq Mean Sq F value  Pr(>F)   
#group         6  35.71   5.951   3.343 0.00439 **
#  Residuals   122 217.17   1.780                   
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

TukeyHSD(aov_cfu_6)
DunnettTest(x=all_cfus_6$CFU, g=all_cfus_6$group, "BHI")

#Dunnett's test for comparing several treatments with a control :  
#    95% family-wise confidence level

#$BHI
#                diff    lwr.ci     upr.ci    pval    
#Ace5-BHI  -0.2088387 -1.329607  0.9119296 0.99361    
#Ace25-BHI -0.6262173 -1.746986  0.4945510 0.51114    
#Pro5-BHI  -0.4432138 -1.563982  0.6775546 0.81285    
#Pro25-BHI -0.5969557 -1.717724  0.5238126 0.56032    
#But5-BHI  -0.3835006 -1.504269  0.7372677 0.89003    
#But25-BHI -1.7616779 -2.882446 -0.6409096 0.00043 ***

#---
#Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

all_cfus_12 <- subset(all_cfus_new, all_cfus_new$time == 12)
aov_cfu_12 <- aov(CFU~group, data = all_cfus_12)
summary(aov_cfu_12)
## no sig

all_cfus_18 <- subset(all_cfus_new, all_cfus_new$time == 18)
aov_cfu_18 <- aov(CFU~group, data = all_cfus_18)
summary(aov_cfu_18)

all_cfus_24 <- subset(all_cfus_new, all_cfus_new$time == 24)
aov_cfu_24 <- aov(CFU~group, data = all_cfus_24)
summary(aov_cfu_24)
#              Df Sum Sq Mean Sq F value   Pr(>F)    
#group         6  21.55   3.592   9.701 1.01e-08 ***
#  Residuals   122  45.17   0.370                     
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


DunnettTest(x=all_cfus_24$CFU, g=all_cfus_24$group, control = "BHI")

#Dunnett's test for comparing several treatments with a control :  
#    95% family-wise confidence level

#$BHI
#                diff     lwr.ci       upr.ci    pval    
#Ace5-BHI  -0.3582234 -0.8693919  0.152945078  0.2811    
#Ace25-BHI -0.3455320 -0.8567004  0.165636501  0.3156    
#Pro5-BHI  -0.5077497 -1.0189181  0.003418808  0.0523 .  
#Pro25-BHI -0.8350152 -1.3461837 -0.323846737  0.0002 ***
#But5-BHI  -0.9635234 -1.4746919 -0.452354972 1.6e-05 ***
#But25-BHI -1.2769698 -1.7881383 -0.765801371 1.0e-08 ***

#---
#Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1





###################################Supplementary Figure S1##################################################################################################

##Fig S1: growth curves with R20291 and VPI10463

## VPI10463
#### I am making the merged meta and expt so that I can average them over 3 biological replicates and then subtract from averaged negative controls

#expt1_VPI <- read_tsv("./SCFACU23_data.txt")
#expt1_meta_VPI <- read_tsv("./SCFACU23_meta.txt")
#expt2_VPI <- read_tsv("./SCFACU24_data.txt")
#expt2_meta_VPI <- read_tsv("./SCFACU24_meta.txt")
#expt3_VPI <- read_tsv("./SCFACU25_data.txt")
#expt3_meta_VPI <- read_tsv("./SCFACU25_meta.txt")

### For experiment1
#expt1_VPI_2 <- expt1_VPI[,-2]
#expt1_VPI_2 <- mutate(expt1_VPI_2, time = as.numeric(gsub("s", "", time)))
#expt1_VPI_3 <- as.data.frame(t(expt1_VPI_2))
#expt1_VPI_3$wellID <- rownames(expt1_VPI_3)
#expt1_VPI_3[1,97] <- "wellID"
#colnames(expt1_VPI_3) <- expt1_VPI_3[1,]
#expt1_VPI_4 <- expt1_VPI_3[-1,]

#merge_expt1_VPI <- inner_join(expt1_meta_VPI, expt1_VPI_4, by = "wellID")
#merge_expt1_VPI <- merge_expt1_VPI[-(85:96),]

### For experiment2
#expt2_VPI$time <- expt1_VPI$time
#expt2_VPI_2 <- expt2_VPI[,-2]
#expt2_VPI_2 <- mutate(expt2_VPI_2, time = as.numeric(gsub("s", "", time)))
#expt2_VPI_3 <- as.data.frame(t(expt2_VPI_2))
#expt2_VPI_3$wellID <- rownames(expt2_VPI_3)
#expt2_VPI_3[1,97] <- "wellID"
#colnames(expt2_VPI_3) <- expt2_VPI_3[1,]
#expt2_VPI_4 <- expt2_VPI_3[-1,]

#merge_expt2_VPI <- inner_join(expt2_meta_VPI, expt2_VPI_4, by = "wellID")
#merge_expt2_VPI <- merge_expt2_VPI[-(85:96),]

### For experiment3
#expt3_VPI$time <- expt1_VPI$time
#expt3_VPI_2 <- expt3_VPI[,-2]
#expt3_VPI_2 <- mutate(expt3_VPI_2, time = as.numeric(gsub("s", "", time)))
#expt3_VPI_3 <- as.data.frame(t(expt3_VPI_2))
#expt3_VPI_3$wellID <- rownames(expt3_VPI_3)
#expt3_VPI_3[1,97] <- "wellID"
#colnames(expt3_VPI_3) <- expt3_VPI_3[1,]
#expt3_VPI_4 <- expt3_VPI_3[-1,]

#merge_expt3_VPI <- inner_join(expt3_meta_VPI, expt3_VPI_4, by = "wellID")
#merge_expt3_VPI <- merge_expt3_VPI[-(85:96),]

### All experiments are combined and averaged now

#full_merge_VPI <- rbind(merge_expt1_VPI, merge_expt2_VPI, merge_expt3_VPI)
#levels(factor(full_merge$group))
#[1] "Ace25" "Ace5"  "BHI"   "But25" "But5"  "Pro25" "Pro5" 
#levels(factor(full_merge$expID))
#[1] "SCFACU21" "SCFACU22" "SCFACU29"
#write_tsv(full_merge_avg_VPI, file = "/Data/figure1_S1/full.merge.VPI.tsv")
full_merge_VPI <- read_tsv(file= "/Data/figure1_S1/full.merge.VPI.tsv")
full_merge_VPI_2 <- full_merge_VPI %>% pivot_longer(cols = !(1:10), names_to = "time", values_to = "OD")
full_merge_avg_VPI <- full_merge_VPI %>% 
  group_by(expID,group) %>% 
  summarise_all(mean) %>% 
  as.data.frame()

## Input the control plate and modify to be able to subtract the controls

controls <- read_tsv("/Data/figure1_S1/controls_new.txt")
controls_meta <- read_tsv("/Data/figure1_S1/controls_new_meta.txt")
controls_GC <- subset(controls_meta, controls_meta$medium == "BHI")

controls$time <- expt1$time
controls_2 <- controls[,-2]
controls_2 <- mutate(controls_2, time = as.numeric(gsub("s", "", time)))
controls_3 <- as.data.frame(t(controls_2))
controls_3$wellID <- rownames(controls_3)
controls_3[1,97] <- "wellID"
colnames(controls_3) <- controls_3[1,]
controls_4 <- controls_3[-1,]

merge_controls <- inner_join(controls_GC, controls_4, by = "wellID")

## Create avg of the controls and subtract from the average of the 3 expts
norm.mat <- merge_controls %>% 
  group_by(expID, group) %>% 
  summarise_all(mean) %>% 
  as.data.frame() %>% 
  select(10:105) %>% 
  as.matrix()

normalized_data_VPI <- sweep(as.matrix(full_merge_avg_VPI[,11:106]), 2, norm.mat)

combined_VPI <- select(full_merge_avg_VPI, 1:2) %>% 
  cbind(., normalized_data_VPI)

combined2_VPI <- combined_VPI %>% pivot_longer(cols=!(1:2), names_to = "time", values_to = "OD")

final_VPI <- combined2_VPI %>% 
  #filter(expID %in% c("SCFACU21", "SCFACU22", "SCFACU29")) %>%
  filter(group %in% c("Ace25", "Ace5",  "BHI", "But25", "But5",  "Pro25", "Pro5")) %>%
  droplevels() %>%
  group_by(time, group) %>%
  dplyr::summarise(mean = mean(OD, na.rm = TRUE), n = length(OD), sd = sd(OD, na.rm = TRUE), se = sd / sqrt(n))

final_VPI <- mutate(final_VPI, time = as.numeric(as.numeric(time)/3600))



## STATS

#combo_VPI <- combined_VPI
#combo2_VPI <- combo_VPI[,-1]
#combo3_VPI <- as.data.frame(t(combo2_VPI))
#combo3_VPI$time <- rownames(combo3_VPI)
#write_tsv(combo3_VPI, "./combo3_VPI10463_gcvr.txt")

##Changes made in excel: put time in front, removed the V1, V2 (row 1)
combo4_VPI <- read_tsv(file= "/Data/figure1_S1/combo3_VPI10463_gcvr_edit.txt")
gcvr_out_VPI <- SummarizeGrowthByPlate(combo4_VPI)
gcvr_out_VPI$group[gcvr_out_VPI$sample == "Ace25...2" | gcvr_out_VPI$sample == "Ace25...9" | gcvr_out_VPI$sample == "Ace25...16"] <- "Ace25"
gcvr_out_VPI$group[gcvr_out_VPI$sample == "Ace5...3" | gcvr_out_VPI$sample == "Ace5...10" | gcvr_out_VPI$sample == "Ace5...17"] <- "Ace5"
gcvr_out_VPI$group[gcvr_out_VPI$sample == "BHI...4" | gcvr_out_VPI$sample == "BHI...11" | gcvr_out_VPI$sample == "BHI...18"] <- "BHI"
gcvr_out_VPI$group[gcvr_out_VPI$sample == "But25...5" | gcvr_out_VPI$sample == "But25...12" | gcvr_out_VPI$sample == "But25...19"] <- "But25"
gcvr_out_VPI$group[gcvr_out_VPI$sample == "But5...6" | gcvr_out_VPI$sample == "But5...13" | gcvr_out_VPI$sample == "But5...20"] <- "But5"
gcvr_out_VPI$group[gcvr_out_VPI$sample == "Pro25...7" | gcvr_out_VPI$sample == "Pro25...14" | gcvr_out_VPI$sample == "Pro25...21"] <- "Pro25"
gcvr_out_VPI$group[gcvr_out_VPI$sample == "Pro5...8" | gcvr_out_VPI$sample == "Pro5...15" | gcvr_out_VPI$sample == "Pro5...22"] <- "Pro5"

aov_model_VPI <- aov(auc_l ~ group, data = gcvr_out_VPI)

summary(aov_model_VPI)
#             Df    Sum Sq  Mean Sq F value Pr(>F)
#group        6  41138321  6856387   0.134   0.99
#Residuals   14 718567141 51326224    

TukeyHSD(aov_model_VPI)

#Tukey multiple comparisons of means
#95% family-wise confidence level

#Fit: aov(formula = auc_e ~ group, data = gcvr_out_VPI)

#$group
#                   diff       lwr      upr     p adj
#Ace5-Ace25  -2212.7856 -22186.65 17761.08 0.9996773
#BHI-Ace25   -1313.6125 -21287.48 18660.25 0.9999846
#But25-Ace25 -4864.2184 -24838.08 15109.65 0.9771394
#But5-Ace25  -1872.8802 -21846.75 18100.99 0.9998770
#Pro25-Ace25 -2933.0836 -22906.95 17040.78 0.9984094
#Pro5-Ace25  -1565.5062 -21539.37 18408.36 0.9999568
#BHI-Ace5      899.1730 -19074.69 20873.04 0.9999984
#But25-Ace5  -2651.4328 -22625.30 17322.43 0.9990962
#But5-Ace5     339.9053 -19633.96 20313.77 1.0000000
#Pro25-Ace5   -720.2980 -20694.16 19253.57 0.9999996
#Pro5-Ace5     647.2794 -19326.59 20621.15 0.9999998
#But25-BHI   -3550.6058 -23524.47 16423.26 0.9954766
#But5-BHI     -559.2677 -20533.13 19414.60 0.9999999
#Pro25-BHI   -1619.4710 -21593.34 18354.40 0.9999474
#Pro5-BHI     -251.8936 -20225.76 19721.97 1.0000000
#But5-But25   2991.3381 -16982.53 22965.20 0.9982259
#Pro25-But25  1931.1348 -18042.73 21905.00 0.9998530
#Pro5-But25   3298.7122 -16675.15 23272.58 0.9969623
#Pro25-But5  -1060.2033 -21034.07 18913.66 0.9999957
#Pro5-But5     307.3740 -19666.49 20281.24 1.0000000
#Pro5-Pro25   1367.5774 -18606.29 21341.44 0.9999805

DunnettTest(x=gcvr_out_R2$auc_l, g=gcvr_out_R2$group, control = "BHI")

##Dunnett's test for comparing several treatments with a control :  
#    95% family-wise confidence level

#$BHI
#               diff    lwr.ci   upr.ci   pval    
#Ace25-BHI  857.5685 -3161.485 4876.622 0.9712    
#Ace5-BHI  1191.2731 -2827.780 5210.326 0.8891    
#But25-BHI -203.9019 -4222.955 3815.151 1.0000    
#But5-BHI  2101.5686 -1917.484 6120.622 0.4814    
#Pro25-BHI -397.9154 -4416.968 3621.138 0.9995    
#Pro5-BHI  1073.7894 -2945.264 5092.843 0.9251    

#---
#Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1


GC_VPI <- lme(as.numeric(time)~1+ group, data=full_merge_VPI_2, random = ~1|OD, na.action = na.omit)
GC_VPI <- lme(mean~group, data=final_VPI, random = ~1|time, na.action = na.omit)
summary(GC_VPI)

## R20291

#### I am making the merged meta and expt so that I can average them over 3 biological replicates and then subtract from averaged negative controls


#expt1_R2 <- read_tsv("./SCFACU26_data.txt")
#expt1_meta_R2 <- read_tsv("./SCFACU26_meta.txt")
#expt2_R2 <- read_tsv("./SCFACU27_data.txt")
#expt2_meta_R2 <- read_tsv("./SCFACU27_meta.txt")
#expt3_R2 <- read_tsv("./SCFACU28_data.txt")
#expt3_meta_R2 <- read_tsv("./SCFACU28_meta.txt")

### For experiment1
#expt1_R2_2 <- expt1_R2[,-2]
#expt1_R2_2 <- mutate(expt1_R2_2, time = as.numeric(gsub("s", "", time)))
#expt1_R2_3 <- as.data.frame(t(expt1_R2_2))
#expt1_R2_3$wellID <- rownames(expt1_R2_3)
#expt1_R2_3[1,97] <- "wellID"
#colnames(expt1_R2_3) <- expt1_R2_3[1,]
#expt1_R2_4 <- expt1_R2_3[-1,]

#merge_expt1_R2 <- inner_join(expt1_meta_R2, expt1_R2_4, by = "wellID")
#merge_expt1_R2 <- merge_expt1_R2[-(85:96),]

### For experiment2
#expt2_R2$time <- expt1_R2$time
#expt2_R2_2 <- expt2_R2[,-2]
#expt2_R2_2 <- mutate(expt2_R2_2, time = as.numeric(gsub("s", "", time)))
#expt2_R2_3 <- as.data.frame(t(expt2_R2_2))
#expt2_R2_3$wellID <- rownames(expt2_R2_3)
#expt2_R2_3[1,97] <- "wellID"
#colnames(expt2_R2_3) <- expt2_R2_3[1,]
#expt2_R2_4 <- expt2_R2_3[-1,]

#merge_expt2_R2 <- inner_join(expt2_meta_R2, expt2_R2_4, by = "wellID")
#merge_expt2_R2 <- merge_expt2_R2[-(85:96),]

### For experiment3
#expt3_R2$time <- expt1_R2$time
#expt3_R2_2 <- expt3_R2[,-2]
#expt3_R2_2 <- mutate(expt3_R2_2, time = as.numeric(gsub("s", "", time)))
#expt3_R2_3 <- as.data.frame(t(expt3_R2_2))
#expt3_R2_3$wellID <- rownames(expt3_R2_3)
#expt3_R2_3[1,97] <- "wellID"
#colnames(expt3_R2_3) <- expt3_R2_3[1,]
#expt3_R2_4 <- expt3_R2_3[-1,]

#merge_expt3_R2 <- inner_join(expt3_meta_R2, expt3_R2_4, by = "wellID")
#merge_expt3_R2 <- merge_expt3_R2[-(85:96),]

### All experiments are combined and averaged now

#full_merge_R2 <- rbind(merge_expt1_R2, merge_expt2_R2, merge_expt3_R2)
#levels(factor(full_merge$group))
#[1] "Ace25" "Ace5"  "BHI"   "But25" "But5"  "Pro25" "Pro5" 
#levels(factor(full_merge$expID))
#[1] "SCFACU21" "SCFACU22" "SCFACU29"

#write_tsv(full_merge_R2, file = "/Data/figure1_S1/full.merge.R2.tsv")
full_merge_R2 <- read_tsv(file = "/Data/figure1_S1/full.merge.R2.tsv")
full_merge_R2_2 <- full_merge_R2 %>% pivot_longer(cols = !(1:10), names_to = "time", values_to = "OD")
full_merge_avg_R2 <- full_merge_R2 %>% 
  group_by(expID,group) %>% 
  summarise_all(mean) %>% 
  as.data.frame()

## Input the control plate and modify to be able to subtract the controls

controls <- read_tsv("/Data/figure1_S1/controls_new.txt")
controls_meta <- read_tsv("/Data/figure1_S1/controls_new_meta.txt")
controls_GC <- subset(controls_meta, controls_meta$medium == "BHI")

controls$time <- expt1_R2$time
controls_2 <- controls[,-2]
controls_2 <- mutate(controls_2, time = as.numeric(gsub("s", "", time)))
controls_3 <- as.data.frame(t(controls_2))
controls_3$wellID <- rownames(controls_3)
controls_3[1,97] <- "wellID"
colnames(controls_3) <- controls_3[1,]
controls_4 <- controls_3[-1,]

merge_controls <- inner_join(controls_GC, controls_4, by = "wellID")

## Create avg of the controls and subtract from the average of the 3 expts
norm.mat <- merge_controls %>% 
  group_by(expID, group) %>% 
  summarise_all(mean) %>% 
  as.data.frame() %>% 
  select(10:105) %>% 
  as.matrix()

normalized_data_R2 <- sweep(as.matrix(full_merge_avg_R2[,11:106]), 2, norm.mat)

combined_R2 <- select(full_merge_avg_R2, 1:2) %>% 
  cbind(., normalized_data_R2)

combined2_R2 <- combined_R2 %>% pivot_longer(cols=!(1:2), names_to = "time", values_to = "OD")

final_R2 <- combined2_R2 %>% 
  #filter(expID %in% c("SCFACU21", "SCFACU22", "SCFACU29")) %>%
  filter(group %in% c("Ace25", "Ace5",  "BHI", "But25", "But5",  "Pro25", "Pro5")) %>%
  droplevels() %>%
  group_by(time, group) %>%
  dplyr::summarise(mean = mean(OD, na.rm = TRUE), n = length(OD), sd = sd(OD, na.rm = TRUE), se = sd / sqrt(n))

final_R2 <- mutate(final_R2, time = as.numeric(as.numeric(time)/3600))

## STATS

#combo_R2 <- combined_R2
#combo2_R2 <- combo_R2[,-1]
#combo3_R2 <- as.data.frame(t(combo2_R2))
#combo3_R2$time <- rownames(combo3_R2)
#write_tsv(combo3_R2, "./combo3_R20291_gcvr.txt")

##Changes made in excel: put time in front, removed the V1, V2 (row 1)
combo4_R2 <- read_tsv(file= "/Data/figure1_S1/combo3_R20291_gcvr_edit.txt")
gcvr_out_R2 <- SummarizeGrowthByPlate(combo4_R2)
gcvr_out_R2$group[gcvr_out_R2$sample == "Ace25...2" | gcvr_out_R2$sample == "Ace25...9" | gcvr_out_R2$sample == "Ace25...16"] <- "Ace25"
gcvr_out_R2$group[gcvr_out_R2$sample == "Ace5...3" | gcvr_out_R2$sample == "Ace5...10" | gcvr_out_R2$sample == "Ace5...17"] <- "Ace5"
gcvr_out_R2$group[gcvr_out_R2$sample == "BHI...4" | gcvr_out_R2$sample == "BHI...11" | gcvr_out_R2$sample == "BHI...18"] <- "BHI"
gcvr_out_R2$group[gcvr_out_R2$sample == "But25...5" | gcvr_out_R2$sample == "But25...12" | gcvr_out_R2$sample == "But25...19"] <- "But25"
gcvr_out_R2$group[gcvr_out_R2$sample == "But5...6" | gcvr_out_R2$sample == "But5...13" | gcvr_out_R2$sample == "But5...20"] <- "But5"
gcvr_out_R2$group[gcvr_out_R2$sample == "Pro25...7" | gcvr_out_R2$sample == "Pro25...14" | gcvr_out_R2$sample == "Pro25...21"] <- "Pro25"
gcvr_out_R2$group[gcvr_out_R2$sample == "Pro5...8" | gcvr_out_R2$sample == "Pro5...15" | gcvr_out_R2$sample == "Pro5...22"] <- "Pro5"

aov_model_R2 <- aov(auc_l ~ group, data = gcvr_out_R2)

summary(aov_model_R2)

#             Df   Sum Sq Mean Sq F value Pr(>F)
#group        6 14615210 2435868   0.855   0.55
#Residuals   14 39901027 2850073              

TukeyHSD(aov_model_R2)

#Tukey multiple comparisons of means
#95% family-wise confidence level

#Fit: aov(formula = auc_l ~ group, data = gcvr_out_R2)

#$group
#diff       lwr      upr     p adj
#Ace5-Ace25    333.7045 -4373.039 5040.448 0.9999761
#BHI-Ace25    -857.5685 -5564.312 3849.175 0.9948382
#But25-Ace25 -1061.4704 -5768.214 3645.273 0.9843775
#But5-Ace25   1244.0001 -3462.743 5950.743 0.9661392
#Pro25-Ace25 -1255.4839 -5962.227 3451.259 0.9646441
#Pro5-Ace25    216.2209 -4490.522 4922.964 0.9999982
#BHI-Ace5    -1191.2731 -5898.016 3515.470 0.9724461
#But25-Ace5  -1395.1749 -6101.918 3311.568 0.9427568
#But5-Ace5     910.2956 -3796.448 5617.039 0.9929175
#Pro25-Ace5  -1589.1884 -6295.932 3117.555 0.9003398
#Pro5-Ace5    -117.4836 -4824.227 4589.260 1.0000000
#But25-BHI    -203.9019 -4910.645 4502.841 0.9999987
#But5-BHI     2101.5686 -2605.175 6808.312 0.7272226
#Pro25-BHI    -397.9154 -5104.658 4308.828 0.9999327
#Pro5-BHI     1073.7894 -3632.954 5780.533 0.9834467
#But5-But25   2305.4705 -2401.273 7012.214 0.6427666
#Pro25-But25  -194.0135 -4900.757 4512.730 0.9999990
#Pro5-But25   1277.6913 -3429.052 5984.434 0.9616260
#Pro25-But5  -2499.4840 -7206.227 2207.259 0.5606055
#Pro5-But5   -1027.7792 -5734.522 3678.964 0.9867270
#Pro5-Pro25   1471.7048 -3235.038 6178.448 0.9277265

DunnettTest(x=gcvr_out_R2$auc_l, g=gcvr_out_R2$group, control = "BHI")

#Dunnett's test for comparing several treatments with a control :  
#    95% family-wise confidence level

#$BHI
#               diff    lwr.ci   upr.ci   pval    
#Ace25-BHI  857.5685 -3161.485 4876.622 0.9712    
#Ace5-BHI  1191.2731 -2827.780 5210.326 0.8890    
#But25-BHI -203.9019 -4222.955 3815.151 1.0000    
#But5-BHI  2101.5686 -1917.484 6120.622 0.4814    
#Pro25-BHI -397.9154 -4416.968 3621.138 0.9995    
#Pro5-BHI  1073.7894 -2945.264 5092.843 0.9251      

#---
#Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

GC_R2 <- lme(as.numeric(time)~1+ group, data=full_merge_R2_2, random = ~1|OD, na.action = na.omit)
GC_R2 <- lme(mean~group, data=final_R2, random = ~1|time, na.action = na.omit)
summary(GC_R2)

## Final plotting

## Plotting VPI and R2 on the same y axis


final_figS1 <- final_VPI
final_figS1$spec <- "VPI"
final_figS1 <- rbind(final_figS1, final_R2)
final_figS1$spec[is.na(final_figS1$spec)] <- "R2"

ggplot(final_figS1, aes(x=time, y=(mean), colour=group)) +
  geom_smooth(size = 0.5, se = FALSE) + 
  #geom_point(size=.1, alpha = 0.5) +
  #geom_errorbar(aes(ymin=(mean-se), ymax=(mean+se)), width=.01, alpha =0.5) +
  theme(legend.position = "right") +
  scale_colour_manual(values=c("black","#801100", "#FC6400","#5D78AB","#87A6BD","#145425", "#80987c")) +
  ylab("log10(OD600)") +
  xlab("Time(hrs)") +
  theme_minimal() +
  facet_wrap(~spec) +
  scale_y_log10() +
  scale_x_continuous(breaks = seq(0, 24, by = 4)) +
  theme(legend.key=element_rect(fill='gray99')) +
  theme(axis.title.x=element_text(size=12, face='bold'),
        axis.title.y=element_text(size=12, face='bold'),
        panel.background = element_rect(fill = 'gray99'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(fill=NA, colour = NA))
