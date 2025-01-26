## Paper: SCFA and C. difficile
## Figure 3: toxin production in presence of scfa
### 3a: toxin production CD630 in different scfa
### Supplementary figure S3: toxin production in VPI10463, R20291; and toxin production over time 630
### 3b: toxin production in presence of butyrate concentration gradient
### 3d: tcdC, tcdR qPCR results


rm(list = ls())
setwd("/Data/figure3_S3/")

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

# Figure 3a: effect of butyrate on toxin production

tox_630 <- read_tsv(file = "/Data/figure3_S3/ToxinData_ForGraph630.txt")

cfu_630_24 <- read_excel(path = "/Data/figure3_S3/new_fig1B_data_630.xlsx")
#cfu_630_24 <- cfu_630_24[,-5]
#cfu_630_24_5x <- subset(cfu_630_24, !(cfu_630_24$expID %in% c("SCFACU19", "SCFACU22")))
#cfu_630_24_3x <- subset(cfu_630_24, !(cfu_630_24$expID %in% c("SCFACU19", "SCFACU22", "SCFACU29")))

cfus_630_24.avg <- cfu_630_24 %>% 
  #filter(expID %in% c("SCFACU21", "SCFACU22", "SCFACU29")) %>%
  filter(group %in% c("Ace25", "Ace5",  "BHI", "But25", "But5",  "Pro25", "Pro5")) %>%
  droplevels() %>%
  group_by(sampleID, group) %>%
  dplyr::summarise(mean = mean(Colonization, na.rm = TRUE))

cfus_630_24.avg$group <- factor(cfus_630_24.avg$group, levels = c("BHI", "Ace5", "Ace25", "Pro5", "Pro25", "But5", "But25"))
#cfu_630_24$group <- factor(cfu_630_24$group, levels = c("BHI", "Ace5", "Ace25", "Pro5", "Pro25", "But5", "But25"))
cfus_630_24.avg$CFU_avg <- log10(cfus_630_24.avg$mean)
new.tox.avgcfu_24 <- inner_join(tox_630, cfus_630_24.avg, by = c("sampleID", "group"))
new.tox.avgcfu_24$avgtox_cfu <- new.tox.avgcfu_24$Toxin/new.tox.avgcfu_24$CFU_avg

new.tox.avgcfu_24$group <- factor(new.tox.avgcfu_24$group, levels = c("BHI", "Ace5", "Ace25", "Pro5", "Pro25", "But5", "But25"))

ggplot(new.tox.avgcfu_24, aes(x=group, y=avgtox_cfu, color=group)) + geom_point(size =2, position = position_jitter(width = 0.2)) +
  #geom_point(size =0.5) + 
  #facet_grid(~time) +
  geom_boxplot(alpha =0.4, size = 0.3, aes(fill = group), outlier.shape = NA) +
  theme_bw() +
  scale_y_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1), limits = c(0,1)) +
  #ylim(0,8.5) +
  #geom_errorbar(aes(ymin=(mean-se), ymax=(mean+se)), width=.01, alpha =0.5) +
  theme(legend.position = "right", panel.grid = element_blank()) +
  scale_color_manual(values=c("black","#801100", "#FC6400", "#5D78AB","#87A6BD","#145425", "#80987c")) +
  scale_fill_manual(values=c("black","#801100", "#FC6400", "#5D78AB","#87A6BD","#145425", "#80987c"))

##STATS

aov_new.tox.24.all <- aov(avgtox_cfu~group, data = new.tox.avgcfu_24)
summary(aov_new.tox.24.all)

#Df Sum Sq Mean Sq F value  Pr(>F)    
#group        6  1.076 0.17937   10.62 2.5e-08 ***
#  Residuals   69  1.165 0.01689                    
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

DunnettTest(x=new.tox.avgcfu_24$avgtox_cfu, g=new.tox.avgcfu_24$group, "BHI")

# Dunnett's test for comparing several treatments with a control :  
#95% family-wise confidence level

#$BHI
#diff       lwr.ci    upr.ci    pval    
#Ace5-BHI  0.1441605  0.004591483 0.2837296 0.03986 *  
#  Ace25-BHI 0.1226203 -0.016948753 0.2621893 0.10857    
#Pro5-BHI  0.2066103  0.067041294 0.3461794 0.00110 ** 
#  Pro25-BHI 0.2237659  0.084196843 0.3633349 0.00034 ***
#  But5-BHI  0.2114287  0.071859671 0.3509977 0.00082 ***
#  But25-BHI 0.3992535  0.259684447 0.5388225 3.7e-10 ***

#  ---
#  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1


############################################################################


## Figure 3b: butyrate concentration gradient on C.diff tox


all_cfus <- read_excel(path = "/Data/figure3_S3/full_cfu_3_expt.xlsx")
all_cfus <- all_cfus[,-5]
added_data <- read_xlsx("/Data/figure3_S3/full_fig2_cfus.xlsx")
added_data <- added_data[,-5]

all_cfus2 <- rbind(all_cfus, added_data)
all_cfus3 <- subset(all_cfus2, all_cfus2$type == "unheated")
all_cfus_new <- all_cfus3 %>% 
  filter(group %in% c("Ace25", "Ace5",  "BHI", "But25", "But5",  "Pro25", "Pro5")) %>%
  droplevels()

tox_630_time <- read_tsv(file = "/Data/figure3_S3/toxin_time_allSCFA_Julian.txt")


avg.cfu <- all_cfus_new %>% 
  group_by(sampleID, group, time) %>% 
  dplyr::summarise(mean_Colonisation = mean(Colonization, na.rm = TRUE))

avg.cfu$CFU_avg <- log10(avg.cfu$mean_Colonisation)
new.tox.avgcfu <- inner_join(tox_630_time, avg.cfu, by = c("sampleID", "group", "time"))
new.tox.avgcfu$avgtox_cfu <- new.tox.avgcfu$Toxin/new.tox.avgcfu$CFU_avg

ggplot(new.tox.avgcfu, aes(x=as.factor(time), y=avgtox_cfu, color=group)) + geom_point(size =2, position = position_jitterdodge(jitter.width = 0.5)) +
  #geom_point(size =0.5) + 
  #facet_grid(~time) +
  geom_boxplot(alpha =0.4, size = 0.3, aes(fill = group), outlier.shape = NA) +
  theme_bw() +
  #geom_errorbar(aes(ymin=(mean-se), ymax=(mean+se)), width=.01, alpha =0.5) +
  theme(legend.position = "right", panel.grid = element_blank()) +
  scale_color_manual(values=c("black","#801100", "#FC6400", "#5D78AB","#87A6BD","#145425", "#80987c")) +
  scale_fill_manual(values=c("black","#801100", "#FC6400", "#5D78AB","#87A6BD","#145425", "#80987c")) + 
  scale_y_continuous(breaks = seq(0, 1.2, by = 0.2), limits = c(0,1.2))



##STAT

newtox_630_time6 <- subset(new.tox.avgcfu, new.tox.avgcfu$time == 6)
aov_newtox_630_time6 <- aov(avgtox_cfu~group, data = newtox_630_time6)
summary(aov_newtox_630_time6)

#Df  Sum Sq Mean Sq F value   Pr(>F)    
#group        6 0.13665 0.02278   17.79 8.21e-09 ***
#  Residuals   31 0.03969 0.00128                     
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

DunnettTest(x=newtox_630_time6$avgtox_cfu, g=newtox_630_time6$group, "BHI")

#Dunnett's test for comparing several treatments with a control :  
#    95% family-wise confidence level

#$BHI
#               diff    lwr.ci    upr.ci    pval    
#Ace5-BHI  0.2186444 0.1432473 0.2940416 4.1e-08 ***
#Ace25-BHI 0.2643788 0.1889817 0.3397760 5.5e-11 ***
#Pro5-BHI  0.2801255 0.2047284 0.3555227 1.7e-11 ***
#Pro25-BHI 0.2511355 0.1757384 0.3265326 8.8e-10 ***
#But5-BHI  0.2592702 0.1838731 0.3346673 1.8e-10 ***
#But25-BHI 0.2613088 0.1859116 0.3367059 2.3e-10 ***

#---
#Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

newtox_630_time12 <- subset(new.tox.avgcfu, new.tox.avgcfu$time == 12)
aov_newtox_630_time12 <- aov(avgtox_cfu~group, data = newtox_630_time12)
summary(aov_newtox_630_time12)

#Df  Sum Sq  Mean Sq F value  Pr(>F)    
#group        6 0.10732 0.017887   7.202 7.1e-05 ***
#  Residuals   31 0.07699 0.002484                    
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

DunnettTest(x=newtox_630_time12$avgtox_cfu, g=newtox_630_time12$group, "BHI")

#Dunnett's test for comparing several treatments with a control :  
#    95% family-wise confidence level

#$BHI
#               diff     lwr.ci    upr.ci    pval    
#Ace5-BHI  0.1447319 0.03972238 0.2497415 0.00474 ** 
#Ace25-BHI 0.2207522 0.11574266 0.3257618 1.8e-05 ***
#Pro5-BHI  0.1790520 0.07404247 0.2840616 0.00038 ***
#Pro25-BHI 0.2183049 0.11329532 0.3233144 2.3e-05 ***
#But5-BHI  0.1193038 0.01429428 0.2243134 0.02268 *  
#But25-BHI 0.1794620 0.07445250 0.2844716 0.00060 ***

#---
#Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

newtox_630_time18 <- subset(new.tox.avgcfu, new.tox.avgcfu$time == 18)
aov_newtox_630_time18 <- aov(avgtox_cfu~group, data = newtox_630_time18)
summary(aov_newtox_630_time18)

#Df  Sum Sq  Mean Sq F value Pr(>F)
#group        6 0.02763 0.004605   0.795  0.583
#Residuals   25 0.14477 0.005791             

DunnettTest(x=newtox_630_time18$avgtox_cfu, g=newtox_630_time18$group, "BHI")

#no sig

newtox_630_time24 <- subset(new.tox.avgcfu, new.tox.avgcfu$time == 24)
aov_newtox_630_time24 <- aov(avgtox_cfu~group, data = newtox_630_time24)
summary(aov_newtox_630_time24)

#Df Sum Sq Mean Sq F value Pr(>F)  
#group        6 0.7023 0.11705   2.937 0.0137 *
#  Residuals   63 2.5106 0.03985                 
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

DunnettTest(x=newtox_630_time24$avgtox_cfu, g=newtox_630_time24$group, "BHI")

#Dunnett's test for comparing several treatments with a control :  
#    95% family-wise confidence level

#$BHI
#                diff       lwr.ci    upr.ci   pval    
#Ace5-BHI  0.25147586  0.003747443 0.4992043 0.0453 *  
#Ace25-BHI 0.02447727 -0.213899511 0.2628540 0.9996    
#Pro5-BHI  0.07058279 -0.183188667 0.3243543 0.9379    
#Pro25-BHI 0.05731089 -0.185361313 0.2999831 0.9699    
#But5-BHI  0.20630013 -0.036372076 0.4489723 0.1215    
#But25-BHI 0.25260525 -0.001166208 0.5063767 0.0515 .  

#---
#Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

##################################################################################

## Figure 3C: 24h dose response butyrate

tox_but_630 <- read_xlsx(path = "/Data/figure3_S3/toxin_Julian_multibut_time.xlsx")
tox_but_630$group <- factor(tox_but_630$group, levels = c("BHI","But5", "But10", "But25", "But50"))

tox_but_630_24 <- subset(tox_but_630, tox_but_630$time==24)

but_24_cfu <- read_xlsx(path = "/Data/figure3_S3/but_dr_24_cfus.xlsx")


avg.cfu.dr <- but_24_cfu %>% 
  group_by(sampleID, group) %>% 
  dplyr::summarise(mean_Colonisation = mean(Colonization, na.rm = TRUE))

avg.cfu.dr$CFU_avg <- log10(avg.cfu.dr$mean_Colonisation)
new.tox.avgcfu <- inner_join(tox_but_630_24, avg.cfu.dr, by = c("sampleID", "group"))
new.tox.avgcfu$avgtox_cfu <- new.tox.avgcfu$Toxin/new.tox.avgcfu$CFU_avg

new.tox.avgcfu$group <- factor(new.tox.avgcfu$group, levels = c("BHI","But5", "But10", "But25", "But50"))


ggplot(new.tox.avgcfu, aes(x=group, y=avgtox_cfu, color=group)) + geom_point(size =2, position = position_jitterdodge(jitter.width = 0.5)) +
  #geom_line() + 
  #facet_grid(~time) +
  geom_boxplot(alpha=0.4, size=0.3, aes(fill = group), outlier.shape = NA) +
  theme_bw() +
  #ylim(0,6.5)+
  #geom_errorbar(aes(ymin=(mean-se), ymax=(mean+se)), width=.01, alpha =0.5) +
  theme(legend.position = "right", panel.grid = element_blank()) +
  scale_color_manual(values=c("black", hcl.colors(4, palette = "TealGrn")))+
  scale_fill_manual(values=c("black",hcl.colors(4, palette = "TealGrn")))+
  scale_y_continuous(breaks = c(0,0.2,0.4,0.6,0.8,1.0), limits = c(0,1))


# STATS

aov_new.tox <- aov(avgtox_cfu~group, data = new.tox.avgcfu)
summary(aov_new.tox)

#Df Sum Sq Mean Sq F value   Pr(>F)    
#group        4 0.2481 0.06202   12.25 1.03e-07 ***
#  Residuals   75 0.3798 0.00506                     
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

DunnettTest(x=new.tox.avgcfu$avgtox_cfu, g=new.tox.avgcfu$group, "BHI")

#Dunnett's test for comparing several treatments with a control :  
#    95% family-wise confidence level

#$BHI
#               diff     lwr.ci    upr.ci    pval    
#But5-BHI  0.1221473 0.05937495 0.1849196 3.3e-05 ***
#But10-BHI 0.1187892 0.05601691 0.1815615 3.5e-05 ***
#But25-BHI 0.1680618 0.10528952 0.2308341 9.9e-09 ***
#But50-BHI 0.1101487 0.04737641 0.1729210 0.00014 ***

#---
#Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1





###################################################################################


## Figure 3D: qPCR of tcdC and tcdR

qpcr <- read_tsv(file = "/Data/figure3_S3/qPCR_full.txt")
qpcr <- subset(qpcr, !is.na(qpcr$time))

qpcr.bhi.log <- qpcr %>% 
  group_by(gene, time) %>% 
  dplyr::summarise(mean = mean(logfc), n = length(logfc), sd = sd(logfc), se = sd / sqrt(n))

qpcr.bhi.log.tox <- subset(qpcr.bhi.log, qpcr.bhi.log$gene == "tcdC" | qpcr.bhi.log$gene == "tcdR")

qpcr.tcd <- subset(qpcr, qpcr$gene == "tcdC" | qpcr$gene == "tcdR")

ggplot(qpcr.bhi.log.tox, aes(x=gene, y=mean, fill=time)) + geom_col(position = position_dodge(), alpha=0.5) + 
  #geom_boxplot(alpha =0.4, size = 0.3, aes(fill = time), outlier.shape = NA) +
  theme_bw() +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.1, position = position_dodge(0.9)) +
  geom_point(aes(x=gene, y=logfc, color=time), data= qpcr.tcd, position = position_dodge(width = 0.8, preserve = "total"), size=2) +
  scale_fill_manual(values = c("orangered3", "dodgerblue")) +
  scale_color_manual(values = c("orangered3", "dodgerblue")) +
  ylab("mean logfc") +
  theme(axis.title.x=element_text(size=12, face='bold'),
        axis.title.y=element_text(size=12, face='bold'),
        panel.background = element_rect(fill = 'gray99'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

## STATS

#### nothing was found significant

qpcr.tcdC <- subset(qpcr.tcd, qpcr.tcd$gene=="tcdC")
aov.tcdC <- aov(logfc~time, data = qpcr.tcdC)
summary(aov.tcdC)

# no sig

qpcr.tcdR <- subset(qpcr.tcd, qpcr.tcd$gene=="tcdR")
aov.tcdR <- aov(logfc~time, data = qpcr.tcdR)
summary(aov.tcdR)

#no sig

kruskal.test(logfc~time, data = qpcr.tcdC)
wilcox.test(logfc~time, data = qpcr.tcdC)
# no sig

kruskal.test(logfc~time, data = qpcr.tcdR)
wilcox.test(logfc~time, data = qpcr.tcdR)
# no sig

tcdC_qpcr.logfc <- subset(qpcr, qpcr$gene == "tcdC")
aov_tox_qpcr.logfc <- aov(logfc~time, data = tcdC_qpcr.logfc)
summary(aov_tox_qpcr.logfc)
DunnettTest(x=tcdC_qpcr.logfc$logfc, g=tcdC_qpcr.logfc$time)
kruskal.test(logfc~time, data = tcdC_qpcr.logfc)
wilcox.test(logfc~time, data = tcdC_qpcr.logfc)

tcdR_qpcr.logfc <- subset(qpcr, qpcr$gene == "tcdR")
aov_tcdR_qpcr.logfc <- aov(logfc~time, data = tcdR_qpcr.logfc)
summary(aov_tcdR_qpcr.logfc)
kruskal.test(logfc~time, data = tcdR_qpcr.logfc)
wilcox.test(logfc~time, data = tcdR_qpcr.logfc)

qpcr.bhi <- qpcr %>% 
  group_by(gene, time) %>% 
  dplyr::summarise(mean = mean(BHI_delCt), n = length(BHI_delCt), sd = sd(BHI_delCt), se = sd / sqrt(n))

qpcr.bhi$type <- "BHI"

qpcr.but <- qpcr %>% 
  group_by(gene, time) %>% 
  dplyr::summarise(mean = mean(Butyrate_delCt), n = length(Butyrate_delCt), sd = sd(Butyrate_delCt), se = sd / sqrt(n))

qpcr.but$type <- "Butyrate"

qpcr.avg <- rbind(qpcr.bhi, qpcr.but)
qpcr.avg.tcdC <- subset(qpcr.avg, qpcr.avg$gene == "tcdC")
qpcr.avg.tcdC.aov <- aov(mean ~ type, data = qpcr.avg.tcdC)
summary(qpcr.avg.tcdC.aov)
qpcr.early.tcdC <- subset(qpcr, qpcr$time== "early" & qpcr$gene == "tcdC")
kruskal.test(BHI_delCt ~ Butyrate_delCt, data = qpcr.early.tcdC)

qpcr.late.tcdC <- subset(qpcr, qpcr$time== "early" & qpcr$gene == "tcdC")
kruskal.test(BHI_delCt ~ Butyrate_delCt, data = qpcr.early.tcdC)

qpcr.avg.tcdR <- subset(qpcr.avg, qpcr.avg$gene == "tcdR")
qpcr.avg.tcdR.aov <- aov(mean ~ type, data = qpcr.avg.tcdR)
summary(qpcr.avg.tcdR.aov)
qpcr.avg.early.tcdR <- subset(qpcr.avg.tcdR, qpcr.avg.tcdR$time== "early")
qpcr.early.tcdR <- subset(qpcr, qpcr$time== "early" & qpcr$gene == "tcdR")
kruskal.test(BHI_delCt ~ Butyrate_delCt, data = qpcr.early.tcdR)

qpcr.late.tcdR <- subset(qpcr, qpcr$time== "late" & qpcr$gene == "tcdR")
kruskal.test(BHI_delCt ~ Butyrate_delCt, data = qpcr.late.tcdR)


#######################Supplementary Figure S3##################################

## Figure S3A: R20291 toxin

tox_r2 <- read_tsv(file = "/Data/figure3_S3/ToxinData_ForGraphR20291.txt")

r2_cfus <- read_xlsx(path = "/Data/figure3_S3/r20291_24_cfus.xlsx")

avg.cfu.r2 <- r2_cfus %>% 
  group_by(sampleID, group) %>% 
  dplyr::summarise(mean_Colonisation = mean(Colonization, na.rm = TRUE))

avg.cfu.r2$CFU_avg <- log10(avg.cfu.r2$mean_Colonisation)
new.tox.avgcfu.r2 <- inner_join(tox_r2, avg.cfu.r2, by = c("sampleID", "group"))
new.tox.avgcfu.r2$avgtox_cfu <- new.tox.avgcfu.r2$Toxin/new.tox.avgcfu.r2$CFU_avg

new.tox.avgcfu.r2$group <- factor(new.tox.avgcfu.r2$group, levels = c("BHI", "Ace5", "Ace25", "Pro5", "Pro25", "But5", "But25"))

ggplot(new.tox.avgcfu.r2, aes(x=group, y=avgtox_cfu, color=group)) + geom_point(size =2, position = position_jitterdodge(jitter.width = 2)) +
  #geom_point(size =0.5) + 
  #facet_grid(~time) +
  geom_boxplot(alpha =0.4, size = 0.3, aes(fill = group), outlier.shape = NA) +
  theme_bw() +
  #geom_errorbar(aes(ymin=(mean-se), ymax=(mean+se)), width=.01, alpha =0.5) +
  theme(legend.position = "right", panel.grid = element_blank()) +
  scale_color_manual(values=c("black","#801100", "#FC6400", "#5D78AB","#87A6BD","#145425", "#80987c")) +
  scale_fill_manual(values=c("black","#801100", "#FC6400", "#5D78AB","#87A6BD","#145425", "#80987c"))
#scale_x_continuous(breaks = seq(0, 24, by = 6))

##STATS
aov.tox.r2 <- aov(avgtox_cfu~group, data = new.tox.avgcfu.r2)
summary(aov.tox.r2)

#Df Sum Sq Mean Sq F value   Pr(>F)    
#group        6 0.7104 0.11839   8.789 3.58e-06 ***
#  Residuals   41 0.5523 0.01347                     
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

DunnettTest(x=new.tox.avgcfu.r2$avgtox_cfu, g=new.tox.avgcfu.r2$group, "BHI")

#Dunnett's test for comparing several treatments with a control :  
#    95% family-wise confidence level

#$BHI
#               diff    lwr.ci    upr.ci    pval    
#Ace5-BHI  0.3911216 0.1951613 0.5870819 5.2e-05 ***
#Ace25-BHI 0.4024002 0.2164959 0.5883044 2.9e-05 ***
#Pro5-BHI  0.5052564 0.3092961 0.7012167 9.0e-07 ***
#Pro25-BHI 0.4307986 0.2448944 0.6167028 2.2e-06 ***
#But5-BHI  0.4136605 0.2277562 0.5995647 2.0e-06 ***
#But25-BHI 0.3571693 0.1712650 0.5430735 6.0e-05 ***

#---
#Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1


########################################################

## Figure S3B: VPI10463 toxin

tox_vpi <- read_tsv(file = "/Data/figure3_S3/ToxinData_ForGraphVPI.txt")
vpi_cfu <- read_xlsx(path = "/Data/figure3_S3/vpi_24_cfus.xlsx")

avg.cfu.vpi <- vpi_cfu %>% 
  group_by(sampleID, group) %>% 
  dplyr::summarise(mean_Colonisation = mean(Colonization, na.rm = TRUE))

avg.cfu.vpi$CFU_avg <- log10(avg.cfu.vpi$mean_Colonisation)
new.tox.avgcfu.vpi <- inner_join(tox_vpi, avg.cfu.vpi, by = c("sampleID", "group"))
new.tox.avgcfu.vpi$avgtox_cfu <- new.tox.avgcfu.vpi$Toxin/new.tox.avgcfu.vpi$CFU_avg

new.tox.avgcfu.vpi$group <- factor(new.tox.avgcfu.vpi$group, levels = c("BHI", "Ace5", "Ace25", "Pro5", "Pro25", "But5", "But25"))

ggplot(new.tox.avgcfu.vpi, aes(x=group, y=avgtox_cfu, color=group)) + geom_point(size =2, position = position_jitterdodge(jitter.width = 2)) +
  #geom_point(size =0.5) + 
  #facet_grid(~time) +
  geom_boxplot(alpha =0.4, size = 0.3, aes(fill = group), outlier.shape = NA) +
  theme_bw() +
  #geom_errorbar(aes(ymin=(mean-se), ymax=(mean+se)), width=.01, alpha =0.5) +
  theme(legend.position = "right", panel.grid = element_blank()) +
  scale_color_manual(values=c("black","#801100", "#FC6400", "#5D78AB","#87A6BD","#145425", "#80987c")) +
  scale_fill_manual(values=c("black","#801100", "#FC6400", "#5D78AB","#87A6BD","#145425", "#80987c"))




## STATS

aov.vpi.tox <- aov(avgtox_cfu~group, data = new.tox.avgcfu.vpi)
summary(aov.vpi.tox)

#Df Sum Sq Mean Sq F value Pr(>F)  
#group        6  0.713 0.11882   2.252 0.0469 *
#  Residuals   77  4.062 0.05275                 
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

DunnettTest(x=new.tox.avgcfu.vpi$avgtox_cfu, g=new.tox.avgcfu.vpi$group, "BHI")

#Dunnett's test for comparing several treatments with a control :  
#    95% family-wise confidence level

#$BHI
#               diff       lwr.ci    upr.ci   pval    
#Ace5-BHI  0.1733475 -0.108099699 0.4547947 0.3694    
#Ace25-BHI 0.2108361 -0.052134975 0.4738072 0.1559    
#Pro5-BHI  0.2908873  0.009440102 0.5723345 0.0401 *  
#Pro25-BHI 0.2992978  0.036326725 0.5622689 0.0199 *  
#But5-BHI  0.3333864  0.070415280 0.5963575 0.0078 ** 
#But25-BHI 0.2402284 -0.022742716 0.5031995 0.0841 .  

#---
#Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1




