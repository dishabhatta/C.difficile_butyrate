## Paper: SCFA and C. difficile
## Figure 2: Growth CD630 in CDMM + butyrate + sugar
### 2: growth curves for CD630 + CDMM + butyrate +sugar
### Supplementary figure S2: ph6.2 7.2 and 8 for CD630

rm(list = ls())
setwd("/Data/figure2_S2/")

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

## Fig 2: 
## Making the combined table for full information
#file_list <- c("SCFACU56", "SCFACU57", "SCFACU58")
#filepath <- "/"
#list.data.cdmm <- list()
#for(i in (1:length(file_list))){
  #print(i)
#  list.data.cdmm[[i]] <- read_tsv(paste0(filepath, file_list[i],'_data.txt'))
  
#}
#list.data.cdmm.meta <- list()
#for(i in (1:length(file_list))){
  #print(i)
#  list.data.cdmm.meta[[i]] <- read_tsv(paste0(filepath, file_list[i],'_meta.txt'))
#}

#full.merge.cdmm <- data.frame()
#time_vec <- list.data.cdmm[[1]][1]
#for (i in 1:3) {
#  list_2 <- list.data.cdmm[[i]][,-2]
#  list_2[,1] <- time_vec$...1
#  colnames(list_2)[1] <- c("time")
#  list_2 <- mutate(list_2, time=as.numeric(gsub("s", "", time)))
#  list_3 <- as.data.frame(t(list_2))
#  list_3$wellID <- rownames(list_3)
#  list_3 <- list_3[,-97]
#  list_3[1,97] <- "wellID"
#  colnames(list_3) <- list_3[1,]
#  list_4 <- list_3[-1,]
#  merge_cdmm <- inner_join(list.data.cdmm.meta[[i]], list_4, by= "wellID")
  #merge_cdmm <- subset(merge_cdmm, !(merge_cdmm$condition == "Cdiff"| merge_cdmm$condition =="noCdiff" | merge_cdmm$condition == "na"))
#  full.merge.cdmm <- rbind(full.merge.cdmm, merge_cdmm)
#}

#write_tsv(full.merge.cdmm, file = "./full.merge.cdmm.tsv")

full.merge.cdmm <- read_tsv(file = "/Data/figure2_S2/full.merge.cdmm.tsv")
full.merge.cdmm <- subset(full.merge.cdmm, !(full.merge.cdmm$condition == "Cdiff"| full.merge.cdmm$condition =="noCdiff" | full.merge.cdmm$condition == "na"))
full.merge.cdmm.avg <- full.merge.cdmm %>% 
  group_by(expID,group) %>% 
  summarise_all(mean) %>% 
  as.data.frame()


## controls
controls_cdmm <- read_tsv("/Data/figure2_S2/controls_new.txt")
controls_meta_cdmm <- read_tsv("/Data/figure2_S2/controls_new_meta.txt")
controls_meta_cdmm <- subset(controls_meta_cdmm, !(controls_meta_cdmm$medium == "BHI" | controls_meta_cdmm$medium == "none" | controls_meta_cdmm$group == "noCdiff"))

#controls$time <- expt1$time
controls_cdmm_2 <- controls_cdmm[,-2]
time_vec2 <- time_vec[-97,]
controls_cdmm_2$time <- time_vec2$...1
controls_cdmm_2 <- mutate(controls_cdmm_2, time = as.numeric(gsub("s", "", time)))
controls_cdmm_3 <- as.data.frame(t(controls_cdmm_2))
controls_cdmm_3$wellID <- rownames(controls_cdmm_3)
controls_cdmm_3[1,97] <- "wellID"
colnames(controls_cdmm_3) <- controls_cdmm_3[1,]
controls_cdmm_4 <- controls_cdmm_3[-1,]

merge_controls_cdmm <- inner_join(controls_meta_cdmm, controls_cdmm_4, by = "wellID")


## Create avg of the controls and subtract from the average of the 3 expts
norm.mat.cdmm <- merge_controls_cdmm %>% 
  group_by(group) %>% 
  summarise_all(mean) %>% 
  as.data.frame() %>% 
  select(10:105) %>% 
  as.matrix()

normalized_data.cdmm <- sweep(as.matrix(full.merge.cdmm.avg[,11:106]), 2, norm.mat.cdmm)

combined.cdmm <- select(full.merge.cdmm.avg, 1:2) %>% 
  cbind(., normalized_data.cdmm)

combined.cdmm2 <- combined.cdmm %>% pivot_longer(cols=!(1:2), names_to = "time", values_to = "OD")
group_lev <- levels(factor(combined.cdmm2$group))
final.cdmm <- combined.cdmm2 %>% 
  #filter(expID %in% c("SCFACU21", "SCFACU22", "SCFACU29")) %>%
  filter(group %in% group_lev) %>%
  droplevels() %>%
  group_by(time, group) %>%
  dplyr::summarise(mean = mean(OD, na.rm = TRUE), n = length(OD), sd = sd(OD, na.rm = TRUE), se = sd / sqrt(n))

final.cdmm <- mutate(final.cdmm, time = as.numeric(as.numeric(time)/3600))

final.cdmm$group_init[final.cdmm$group == "Cellobiose" | final.cdmm$group == "Cellobiose+But"] <- "Cel"
final.cdmm$group_init[final.cdmm$group == "Fructose" | final.cdmm$group == "Fructose+But"] <- "Fru"
final.cdmm$group_init[final.cdmm$group == "Glucose" | final.cdmm$group == "Glucose+But"] <- "Glu"
final.cdmm$group_init[final.cdmm$group == "Lactose" | final.cdmm$group == "Lactose+But"] <- "Lac"
final.cdmm$group_init[final.cdmm$group == "Maltose" | final.cdmm$group == "Maltose+But"] <- "Mal"
final.cdmm$group_init[final.cdmm$group == "Mannitol" | final.cdmm$group == "Mannitol+But"] <- "Mnl"
final.cdmm$group_init[final.cdmm$group == "Mannose" | final.cdmm$group == "Mannose+But"] <- "Man"
final.cdmm$group_init[final.cdmm$group == "Raffinose" | final.cdmm$group == "Raffinose+But"] <- "Raf"
final.cdmm$group_init[final.cdmm$group == "Sucrose" | final.cdmm$group == "Sucrose+But"] <- "Suc"
final.cdmm$group_init[final.cdmm$group == "Trehalose" | final.cdmm$group == "Trehalose+But"] <- "Tre"

sug <- c("Cellobiose","Fructose","Glucose","Lactose","Maltose","Mannitol","Mannose","Raffinose", "Sucrose","Trehalose")
final.cdmm$but[final.cdmm$group %in% sug] <- "Y"
final.cdmm$but[is.na(final.cdmm$but)] <- "N"

ggplot(final.cdmm, aes(x=time, y=(mean), colour=but)) +
  geom_smooth(size = 0.5, se = FALSE) + 
  #geom_point(size=.1, alpha = 0.5) +
  #geom_errorbar(aes(ymin=(mean-se), ymax=(mean+se)), width=.01, alpha =0.5) +
  theme(legend.position = "right") +
  facet_wrap(~group_init)+
  scale_colour_manual(values=c("black", "indianred")) +
  ylab("OD600") +
  xlab("Time(hrs)") +
  theme_minimal() +
  scale_y_log10() +
  scale_x_continuous(breaks = seq(0, 24, by = 4)) +
  theme(legend.key=element_rect(fill='gray99')) +
  theme(axis.title.x=element_text(size=12, face='bold'),
        axis.title.y=element_text(size=12, face='bold'),
        panel.background = element_rect(fill = 'gray99'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(fill=NA, colour = NA))

##STATS

### Separating out sugars
mannose2 <- subset(combined.cdmm2, combined.cdmm2$group == "Mannose" | combined.cdmm2$group == "Mannose+But")
aov_man2 <- aov(OD~group, data = mannose2)
summary(aov_man2)

#Df Sum Sq Mean Sq F value  Pr(>F)    
#group         1   2.06  2.0634   22.55 2.6e-06 ***
#  Residuals   574  52.54  0.0915                    
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

TukeyHSD(aov_man2)
#Tukey multiple comparisons of means
#95% family-wise confidence level

#Fit: aov(formula = OD ~ group, data = mannose2)

#$group
#diff        lwr         upr   p adj
#Mannose+But-Mannose -0.1197049 -0.1692221 -0.07018766 2.6e-06
pairwise.t.test(x=mannose2$OD, g=mannose2$group, p.adjust.method = "bonferroni")
#Pairwise comparisons using t tests with pooled SD 

#data:  mannose2$OD and mannose2$group 

#Mannose
#Mannose+But 2.6e-06

#P value adjustment method: bonferroni 

cellobiose2 <- subset(combined.cdmm2, combined.cdmm2$group == "Cellobiose" | combined.cdmm2$group == "Cellobiose+But")
aov_cel2 <- aov(OD~group, data = cellobiose2)
summary(aov_cel2)

#Df Sum Sq  Mean Sq F value Pr(>F)
#group         1  0.006 0.005688   0.535  0.465
#Residuals   574  6.101 0.010628  

TukeyHSD(aov_cel2)

#Tukey multiple comparisons of means
#95% family-wise confidence level

#Fit: aov(formula = OD ~ group, data = cellobiose2)

#$group
#diff         lwr        upr     p adj
#Cellobiose+But-Cellobiose 0.006284722 -0.01058902 0.02315847 0.4647465

fructose <- subset(combined.cdmm2, combined.cdmm2$group == "Fructose" | combined.cdmm2$group == "Fructose+But")
aov_fru <- aov(OD~group, data = fructose)
summary(aov_fru)

#Df Sum Sq Mean Sq F value  Pr(>F)   
#group         1   0.79  0.7858   9.594 0.00205 **
#  Residuals   574  47.01  0.0819                   
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

TukeyHSD(aov_fru)

#Tukey multiple comparisons of means
#95% family-wise confidence level

#Fit: aov(formula = OD ~ group, data = fructose)

#$group
#diff        lwr         upr     p adj
#Fructose+But-Fructose -0.07387153 -0.1207137 -0.02702935 0.0020473

glucose <- subset(combined.cdmm2, combined.cdmm2$group == "Glucose" | combined.cdmm2$group == "Glucose+But")
aov_glu <- aov(OD~group, data = glucose)
summary(aov_glu)

#Df Sum Sq Mean Sq F value Pr(>F)  
#group         1   0.26 0.25926   3.535 0.0606 .
#Residuals   574  42.10 0.07334                 
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

lactose <- subset(combined.cdmm2, combined.cdmm2$group == "Lactose" | combined.cdmm2$group == "Lactose+But")
aov_lac <- aov(OD~group, data = lactose)
summary(aov_lac)

#Df Sum Sq Mean Sq F value Pr(>F)  
#group         1 0.0167 0.01668   4.633 0.0318 *
#  Residuals   574 2.0664 0.00360                 
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

lactose <- subset(combined.cdmm2, combined.cdmm2$group == "Lactose" | combined.cdmm2$group == "Lactose+But")
aov_lac <- aov(OD~group, data = lactose)
summary(aov_lac)

#Df Sum Sq Mean Sq F value Pr(>F)  
#group         1 0.0167 0.01668   4.633 0.0318 *
#  Residuals   574 2.0664 0.00360                 
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

TukeyHSD(aov_lac)

#Tukey multiple comparisons of means
#95% family-wise confidence level

#Fit: aov(formula = OD ~ group, data = lactose)

#$group
#diff          lwr        upr     p adj
#Lactose+But-Lactose 0.01076273 0.0009421321 0.02058333 0.0317717

maltose <- subset(combined.cdmm2, combined.cdmm2$group == "Maltose" | combined.cdmm2$group == "Maltose+But")
aov_mal <- aov(OD~group, data = maltose)
summary(aov_mal)

#Df Sum Sq  Mean Sq F value Pr(>F)
#group         1 0.0063 0.006340   1.524  0.218
#Residuals   574 2.3886 0.004161        

mannose <- subset(combined.cdmm2, combined.cdmm2$group == "Mannose" | combined.cdmm2$group == "Mannose+But")
aov_man <- aov(OD~group, data = mannose)
summary(aov_man)

#Df Sum Sq Mean Sq F value  Pr(>F)    
#group         1   2.06  2.0634   22.55 2.6e-06 ***
#  Residuals   574  52.54  0.0915                    
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

TukeyHSD(aov_man)

#Tukey multiple comparisons of means
#95% family-wise confidence level

#Fit: aov(formula = OD ~ group, data = mannose)

#$group
#diff        lwr         upr   p adj
#Mannose+But-Mannose -0.1197049 -0.1692221 -0.07018766 2.6e-06

mannitol <- subset(combined.cdmm2, combined.cdmm2$group == "Mannitol" | combined.cdmm2$group == "Mannitol+But")
aov_mnl <- aov(OD~group, data = mannitol)
summary(aov_mnl)

#Df Sum Sq Mean Sq F value Pr(>F)  
#group         1  0.143 0.14251   6.002 0.0146 *
#  Residuals   574 13.628 0.02374                 
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

TukeyHSD(aov_mnl)

#Tukey multiple comparisons of means
#95% family-wise confidence level

#Fit: aov(formula = OD ~ group, data = mannitol)

#$group
#diff         lwr          upr     p adj
#Mannitol+But-Mannitol -0.03145833 -0.05667802 -0.006238651 0.0145844

raffinose <- subset(combined.cdmm2, combined.cdmm2$group == "Raffinose" | combined.cdmm2$group == "Raffinose+But")
aov_raf <- aov(OD~group, data = raffinose)
summary(aov_raf)

#Df Sum Sq Mean Sq F value Pr(>F)  
#group         1  0.022 0.02197   5.705 0.0172 *
#  Residuals   574  2.210 0.00385                 
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

TukeyHSD(aov_raf)

#Tukey multiple comparisons of means
#95% family-wise confidence level

#Fit: aov(formula = OD ~ group, data = raffinose)

#$group
#diff         lwr        upr     p adj
#Raffinose+But-Raffinose 0.01235069 0.002194764 0.02250662 0.0172366

sucrose <- subset(combined.cdmm2, combined.cdmm2$group == "Sucrose" | combined.cdmm2$group == "Sucrose+But")
aov_suc <- aov(OD~group, data = sucrose)
summary(aov_suc)

#Df Sum Sq Mean Sq F value Pr(>F)  
#group         1  0.106 0.10634   3.014 0.0831 .
#Residuals   574 20.250 0.03528                 
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

TukeyHSD(aov_suc)

#Tukey multiple comparisons of means
#95% family-wise confidence level

#Fit: aov(formula = OD ~ group, data = sucrose)

#$group
#diff         lwr         upr     p adj
#Sucrose+But-Sucrose -0.02717477 -0.05791703 0.003567492 0.0830685

trehalose <- subset(combined.cdmm2, combined.cdmm2$group == "Trehalose" | combined.cdmm2$group == "Trehalose+But")
aov_tre <- aov(OD~group, data = trehalose)
summary(aov_tre)

#Df Sum Sq Mean Sq F value Pr(>F)
#group         1   0.02 0.01543   0.176  0.675
#Residuals   574  50.42 0.08783

###############################################################Supplementary Figure S2###############################################

# figS2: pH growth curves: 

#file_list.ph <- c("SCFACU47", "SCFACU50", "SCFACU48", "SCFACU46", "SCFACU49", "SCFACU51")
#filepath.ph <- ""
#list.data.ph <- list()
#for(i in (1:length(file_list.ph))){
  #print(i)
#  list.data.ph[[i]] <- read_tsv(paste0(filepath.ph, file_list.ph[i],'_data.txt'))
  
#}
#list.data.ph.meta <- list()
#for(i in (1:length(file_list.ph))){
  #print(i)
#  list.data.ph.meta[[i]] <- read_tsv(paste0(filepath.ph, file_list.ph[i],'_meta.txt'))
#}

#full.merge.ph <- data.frame()
#time_vec.ph <- list.data.ph[[1]][1]

#for (i in 1:6) {
#  list_2 <- list.data.ph[[i]][,-2]
#  list_2[,1] <- time_vec.ph$...1
#  colnames(list_2)[1] <- c("time")
#  list_2 <- mutate(list_2, time=as.numeric(gsub("s", "", time)))
#  list_3 <- as.data.frame(t(list_2))
#  list_3$wellID <- rownames(list_3)
#  list_3 <- list_3[,-97]
#  list_3[1,97] <- "wellID"
#  colnames(list_3) <- list_3[1,]
#  list_4 <- list_3[-1,]
#  merge_ph <- inner_join(list.data.ph.meta[[i]], list_4, by= "wellID")
#  full.merge.ph <- rbind(full.merge.ph, merge_ph)
#}

#write_tsv(full.merge.ph, file = "/Data/figure2_S2/full.merge.ph.tsv")
full.merge.ph <- read_tsv(file = "/Data/figure2_S2/full.merge.ph.tsv")

full.merge.ph <- subset(full.merge.ph, !(full.merge.ph$group =="noCdiff"))
full.merge.ph.avg <- full.merge.ph %>% 
  group_by(expID,group,ph) %>% 
  summarise_all(mean) %>% 
  as.data.frame()

###controls

controls <- read_tsv("/Data/figure2_S2/controls_new.txt")
controls_meta <- read_tsv("/Data/figure2_S2/controls_new_meta.txt")
controls_GC <- subset(controls_meta, controls_meta$medium == "BHI")
time_vec.ph.2 <- time_vec.ph[-97,]
controls$time <- time_vec.ph.2$...1
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

normalized_data_ph <- sweep(as.matrix(full.merge.ph.avg[,12:107]), 2, norm.mat)

combined_ph <- select(full.merge.ph.avg, 1:3) %>% 
  cbind(., normalized_data_ph)

combined2_ph <- combined_ph %>% pivot_longer(cols=!(1:3), names_to = "time", values_to = "OD")

final_ph <- combined2_ph %>% 
  filter(ph %in% c(6.2, 7.2, 8)) %>%
  filter(group %in% c("Ace25", "Ace5",  "BHI", "But25", "But5",  "Pro25", "Pro5")) %>%
  droplevels() %>%
  group_by(time, ph, group) %>%
  dplyr::summarise(mean = mean(OD, na.rm = TRUE), n = length(OD), sd = sd(OD, na.rm = TRUE), se = sd / sqrt(n))

final_ph <- mutate(final_ph, time = as.numeric(as.numeric(time)/3600))


ggplot(final_ph, aes(x=time, y=(mean), colour=group)) +
  geom_smooth(size = 0.5, se = FALSE) + 
  facet_wrap(~as.factor(ph)) +
  #geom_point(size=.1, alpha = 0.5) +
  #geom_errorbar(aes(ymin=(mean-se), ymax=(mean+se)), width=.01, alpha =0.5) +
  #theme(legend.position = "right") +
  scale_colour_manual(values=c("black","#801100", "#FC6400", "#5D78AB","#87A6BD","#145425", "#80987c")) +
  #scale_color_manual(values = c("gold", "black", "darkorange")) +
  ylab("OD600") +
  xlab("Time(hrs)") +
  theme_minimal() +
  scale_y_log10(breaks=c(0.01,0.03,0.10,0.3)) +
  scale_x_continuous(breaks = seq(0, 24, by = 4)) +
  theme(legend.key=element_rect(fill='gray99')) +
  theme(axis.title.x=element_text(size=12, face='bold'),
        axis.title.y=element_text(size=12, face='bold'),
        panel.background = element_rect(fill = 'gray99'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(fill=NA, colour = NA))

##STATS

ph6.2 <- subset(combined2_ph, combined2_ph$ph == 6.2)
aov_ph6.2 <- aov(OD~group, data = ph6.2)
summary(aov_ph6.2)
pairwise.t.test(x=ph6.2$OD, g=ph6.2$group, p.adjust.method = "bonferroni")
DunnettTest(x=ph6.2$OD, g=ph6.2$group, "BHI")



ph7.2 <- subset(final_ph, final_ph$ph == 7.2)
aov_ph7.2 <- aov(mean~group, data = ph7.2)
summary(aov_ph7.2)

#Df Sum Sq Mean Sq F value  Pr(>F)   
#group         6  0.739 0.12325    3.33 0.00307 **
#  Residuals   665 24.609 0.03701                   
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

TukeyHSD(aov_ph7.2)

#Tukey multiple comparisons of means
#95% family-wise confidence level

#Fit: aov(formula = mean ~ group, data = ph7.2)

#$group
#diff         lwr           upr     p adj
#Ace5-Ace25   0.0373884549 -0.04472719  0.1195041024 0.8295898
#BHI-Ace25    0.0414296875 -0.04068596  0.1235453350 0.7496148
#But25-Ace25 -0.0517530382 -0.13386869  0.0303626093 0.5052192
#But5-Ace25   0.0119908854 -0.07012476  0.0941065329 0.9995048
#Pro25-Ace25 -0.0403619792 -0.12247763  0.0417536684 0.7720897
#Pro5-Ace25   0.0124613715 -0.06965428  0.0945770190 0.9993820
#BHI-Ace5     0.0040412326 -0.07807441  0.0861568802 0.9999992
#But25-Ace5  -0.0891414931 -0.17125714 -0.0070258455 0.0233903 *
#But5-Ace5   -0.0253975694 -0.10751322  0.0567180781 0.9702273
#Pro25-Ace5  -0.0777504340 -0.15986608  0.0043652135 0.0770363
#Pro5-Ace5   -0.0249270833 -0.10704273  0.0571885642 0.9728823
#But25-BHI   -0.0931827257 -0.17529837 -0.0110670782 0.0145993 *
#But5-BHI    -0.0294388021 -0.11155445  0.0526768454 0.9393118
#Pro25-BHI   -0.0817916667 -0.16390731  0.0003239809 0.0516852
#Pro5-BHI    -0.0289683160 -0.11108396  0.0531473315 0.9437145
#But5-But25   0.0637439236 -0.01837172  0.1458595711 0.2475064
#Pro25-But25  0.0113910590 -0.07072459  0.0935067065 0.9996317
#Pro5-But25   0.0642144097 -0.01790124  0.1463300572 0.2393481
#Pro25-But5  -0.0523528646 -0.13446851  0.0297627829 0.4906967
#Pro5-But5    0.0004704861 -0.08164516  0.0825861336 1.0000000
#Pro5-Pro25   0.0528233507 -0.02929230  0.1349389982 0.4793659

DunnettTest(x=ph7.2$mean, g=ph7.2$group, "BHI")

#Dunnett's test for comparing several treatments with a control :  
#    95% family-wise confidence level

#$BHI
#                  diff      lwr.ci      upr.ci   pval    
#Ace25-BHI -0.041429687 -0.11289537  0.03003599 0.4713    
#Ace5-BHI  -0.004041233 -0.07550691  0.06742445 1.0000    
#But25-BHI -0.093182726 -0.16464841 -0.02171705 0.0047 ** 
#But5-BHI  -0.029438802 -0.10090448  0.04202688 0.7828    
#Pro25-BHI -0.081791667 -0.15325735 -0.01032599 0.0174 *  
#Pro5-BHI  -0.028968316 -0.10043400  0.04249736 0.7940    

#---
#Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

pairwise.t.test(x=ph7.2$mean, g=ph7.2$group, p.adjust.method = "bonferroni")

#Pairwise comparisons using t tests with pooled SD 

#data:  ph7.2$mean and ph7.2$group 

#Ace25 Ace5  BHI   But25 But5  Pro25
#Ace5  1.000 -     -     -     -     -    
#  BHI   1.000 1.000 -     -     -     -    
#  But25 1.000 0.029 0.018 -     -     -    
#  But5  1.000 1.000 1.000 0.462 -     -    
#  Pro25 1.000 0.110 0.070 1.000 1.000 -    
#  Pro5  1.000 1.000 1.000 0.442 1.000 1.000

#P value adjustment method: bonferroni 

ph8 <- subset(final_ph, final_ph$ph == 8)
aov_ph8 <- aov(mean~group, data = ph8)
summary(aov_ph8)
pairwise.t.test(x=ph8$mean, g=ph8$group, p.adjust.method = "bonferroni")
DunnettTest(x=ph8$mean, g=ph8$group, "BHI")








