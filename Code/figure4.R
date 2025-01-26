## Paper: SCFA and C. difficile
## Figure 4: spore production
### 4a: spore production CD630 in different scfa
### 4b: spore production in presence of scfa over time
### 4c: spore production in butyrate concentration gradient
### 4d: phase contrast of spores in BHI
### 4e: phase contrast of spores in BHI+butyrate
### 4f: spore efficiency calculation


rm(list = ls())
setwd("/Data/figure4/")

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
library(ggpubr)

## Figure 4: Spore production

## Figure 4a: 24h heated scfa panel - 29,32,33

scfa_24 <- read_xlsx(path = "./new_fig4A_heatscfa_24.xlsx")

scfa_24$group <- factor(scfa_24$group, levels = c("BHI", "Ace5", "Ace25", "Pro5", "Pro25", "But5", "But25"))

ggplot(scfa_24, aes(x=group, y=CFU, color=group)) + geom_point(size =2, position = position_jitter(width = 0.4)) +
  #geom_point(size =0.5) + 
  #facet_grid(~time) +
  geom_boxplot(alpha =0.4, size = 0.3, aes(fill = group), outlier.shape = NA) +
  theme_bw() +
  scale_y_continuous(limits = c(0,8), breaks = c(0,2,4,6,8)) +
  #geom_errorbar(aes(ymin=(mean-se), ymax=(mean+se)), width=.01, alpha =0.5) +
  theme(legend.position = "right", panel.grid = element_blank()) +
  scale_color_manual(values=c("black","#801100", "#FC6400", "#5D78AB","#87A6BD","#145425", "#80987c")) +
  scale_fill_manual(values=c("black","#801100", "#FC6400", "#5D78AB","#87A6BD","#145425", "#80987c"))

## STAT

aov_scfa_24 <- aov(CFU~group, data = scfa_24)
summary(aov_scfa_24)

#Df Sum Sq Mean Sq F value   Pr(>F)    
#group        6  8.904  1.4840    6.59 2.62e-05 ***
#  Residuals   56 12.610  0.2252                     
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

DunnettTest(x=scfa_24$CFU, g=scfa_24$group, "BHI")

#Dunnett's test for comparing several treatments with a control :  
#    95% family-wise confidence level

#$BHI
#                 diff      lwr.ci    upr.ci   pval    
#Ace5-BHI  -0.24216344 -0.83449746 0.3501706 0.7678    
#Ace25-BHI -0.20636130 -0.79869532 0.3859727 0.8649    
#Pro5-BHI   0.04246961 -0.54986440 0.6348036 1.0000    
#Pro25-BHI  0.75752290  0.16518889 1.3498569 0.0069 ** 
#But5-BHI  -0.22259775 -0.81493176 0.3697363 0.8235    
#But25-BHI  0.57819334 -0.01414068 1.1705274 0.0579 .  

#---
#Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1



## Figure 4b:scfa over time spore production

full_cfus <- read_xlsx("./full_fig2_cfus.xlsx")

cfus_630_heat <- subset(full_cfus, full_cfus$type == "heat")

cfus_630_heat$group <- factor(cfus_630_heat$group, levels = c("BHI", "Ace5", "Ace25", "Pro5", "Pro25", "But5", "But25"))


ggplot(cfus_630_heat, aes(x=as.factor(time), y=CFU, color=group)) + geom_point(size =2, position = position_jitterdodge(jitter.width = 0.5)) +
  #geom_point(size =0.5) + 
  #facet_grid(~time) +
  geom_boxplot(alpha =0.4, size = 0.3, aes(fill = group), outlier.shape = NA) +
  theme_bw() +
  scale_y_continuous(limits = c(0,8), breaks = c(0,2,4,6,8)) +
  #geom_errorbar(aes(ymin=(mean-se), ymax=(mean+se)), width=.01, alpha =0.5) +
  theme(legend.position = "right", panel.grid = element_blank()) +
  scale_color_manual(values=c("black","#801100", "#FC6400", "#5D78AB","#87A6BD","#145425", "#80987c")) +
  scale_fill_manual(values=c("black","#801100", "#FC6400", "#5D78AB","#87A6BD","#145425", "#80987c"))


### STATS

cfus_630_heat <- subset(full_cfus, full_cfus$type == "heat")
cfus_6_heat <- subset(cfus_630_heat, cfus_630_heat$time == 6)
aov_cfus_6 <- aov(CFU~group, data = cfus_6_heat)
summary(aov_cfus_6)
DunnettTest(x=cfus_6_heat$CFU, g=cfus_6_heat$group, "BHI")

#Dunnett's test for comparing several treatments with a control :  
#    95% family-wise confidence level

#$BHI
#                 diff     lwr.ci      upr.ci   pval    
#Ace5-BHI   0.03913139 -0.9189772  0.99723994 1.0000    
#Ace25-BHI  0.57355185 -0.3845567  1.53166039 0.4183    
#Pro5-BHI   0.09313985 -0.8649687  1.05124839 0.9998    
#Pro25-BHI  0.02153556 -0.9365730  0.97964410 1.0000    
#But5-BHI  -1.02950805 -1.9876166 -0.07139951 0.0306 *  
#But25-BHI -0.19560250 -1.1537110  0.76250605 0.9872    

#---
#Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1


cfus_12_heat <- subset(cfus_630_heat, cfus_630_heat$time == 12)
aov_cfus_12 <- aov(CFU~group, data = cfus_12_heat)
summary(aov_cfus_12)
DunnettTest(x=cfus_12_heat$CFU, g=cfus_12_heat$group, "BHI")

# No sig

TukeyHSD(aov_cfus_12)
pairwise.t.test(x=cfus_12_heat$CFU, g=cfus_12_heat$group, "bonferroni")

cfus_18_heat <- subset(cfus_630_heat, cfus_630_heat$time == 18)
aov_cfus_18 <- aov(CFU~group, data = cfus_18_heat)
summary(aov_cfus_18)
#Df Sum Sq Mean Sq F value  Pr(>F)   
#group        6  28.10   4.684   3.622 0.00417 **
#  Residuals   56  72.41   1.293                   
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

DunnettTest(x=cfus_18_heat$CFU, g=cfus_18_heat$group, "BHI")

#Dunnett's test for comparing several treatments with a control :  
#   95% family-wise confidence level

#$BHI
#                diff     lwr.ci   upr.ci   pval    
#Ace5-BHI  -0.1975724 -1.6170000 1.221855 0.9983    
#Ace25-BHI  0.8432135 -0.5762142 2.262641 0.4258    
#Pro5-BHI   0.1836903 -1.2357374 1.603118 0.9989    
#Pro25-BHI  1.0427679 -0.3766598 2.462196 0.2273    
#But5-BHI   0.1023132 -1.3171145 1.521741 1.0000    
#But25-BHI  1.8142061  0.3947784 3.233634 0.0069 ** 

#---
#Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1


TukeyHSD(aov_cfus_18)
#Tukey multiple comparisons of means
#95% family-wise confidence level

#Fit: aov(formula = CFU ~ group, data = cfus_18_heat)

#$group
#diff         lwr       upr     p adj
#Ace5-BHI    -0.4197946 -3.21689482 2.3773056 0.9992371
#Ace25-BHI    1.5098801 -1.28722009 4.3069804 0.6506784
#Pro5-BHI     0.1836903 -2.61340995 2.9807905 0.9999940
#Pro25-BHI    1.8205456 -0.97655458 4.6176459 0.4321170
#But5-BHI    -0.1199090 -2.91700926 2.6771912 0.9999995
#But25-BHI    2.9253172  0.12821693 5.7224174 0.0347697



pairwise.t.test(x=cfus_18_heat$CFU, g=cfus_18_heat$group, "bonferroni")

#Pairwise comparisons using t tests with pooled SD 

#data:  cfus_18_heat$CFU and cfus_18_heat$group 

#         BHI   Ace5  Ace25 Pro5  Pro25 But5 
#  Ace5  1.000 -     -     -     -     -    
#  Ace25 1.000 0.827 -     -     -     -    
#  Pro5  1.000 1.000 1.000 -     -     -    
#  Pro25 1.000 0.367 1.000 1.000 -     -    
#  But5  1.000 1.000 1.000 1.000 0.805 -    
#  But25 0.048 0.012 1.000 0.085 1.000 0.032

#P value adjustment method: bonferroni 

cfus_24_heat <- subset(cfus_630_heat, cfus_630_heat$time == 24)
aov_cfus_24 <- aov(CFU~group, data = cfus_24_heat)
summary(aov_cfus_24)
DunnettTest(x=cfus_24_heat$CFU, g=cfus_24_heat$group, "BHI")

## No sig



## Figure 4c: butyrate concentration gradient

cfus_dr <- read_xlsx(path = "~/Box Sync/Disha/Projects/michelle_raw_data/figure_1/dose_but_3_expt.xlsx")
library(ggpubr)
cfus_dr$group <- factor(cfus_dr$group, levels = c("BHI","But5", "But10", "But25", "But50"))
cfus_dr <- subset(cfus_dr, cfus_dr$type == "heat")

ggplot(cfus_dr, aes(x=group, y=(CFU), color=group)) + geom_point(size =0.5, position = position_jitterdodge(jitter.width = 0.75)) +
  #geom_line() + 
  #facet_grid(~time) +
  geom_boxplot(alpha =0.2, size = 0.3, aes(fill = group), outlier.shape = NA) +
  theme_bw() + #stat_compare_means(aes(label = ..p.signif..), method = "kruskal") +
  ylim(0,6) +
  #geom_errorbar(aes(ymin=(mean-se), ymax=(mean+se)), width=.01, alpha =0.5) +
  theme(legend.position = "right", panel.grid = element_blank()) +
  scale_colour_manual(values=c("black", hcl.colors(4, palette = "TealGrn")))+
  scale_fill_manual(values=c("black", hcl.colors(4, palette = "TealGrn")))

#scale_color_manual(values=c("black", "seagreen")) #+
#scale_x_continuous(breaks = seq(0, 24, by = 6))


## STAT
cfus_dr_heat <- subset(cfus_dr, cfus_dr$type == "heat")
DunnettTest(x=cfus_dr$CFU, g=cfus_dr$group, "BHI")

#Dunnett's test for comparing several treatments with a control :  
#    95% family-wise confidence level

#$BHI
#               diff     lwr.ci   upr.ci    pval    
#But5-BHI  0.2667127 -0.6362040 1.169629  0.8660    
#But10-BHI 0.1478881 -0.7550286 1.050805  0.9812    
#But25-BHI 1.8232092  0.9202925 2.726126 2.9e-05 ***
#But50-BHI 2.2740819  1.3711652 3.176999 3.1e-07 ***

#---
#Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1





## figure 4f: SPORE EFFCIENCY

spore_eff <- read.delim(file = "./spore.txt")
spore_eff_2 <- spore_eff %>% pivot_longer(cols = !c(Date., expt), names_to = "media", values_to = "spores")
spore_eff_3 <- spore_eff_2 %>% 
  group_by(expt, media) %>% 
  dplyr::summarise(sum = sum(spores), percent = (sum(spores)/1000) *100)

spore_eff_4 <- spore_eff_3 %>% 
  group_by(media) %>% 
  dplyr::summarise(mean = mean(percent), n = length(percent), sd = sd(percent), se = sd / sqrt(n))

ggplot(spore_eff_4, aes(x=media, y=mean)) + geom_col(aes(fill=media), alpha=0.5) +
  geom_point(aes(x=media, y=percent, color=media), data=spore_eff_3, size =2, position = position_jitterdodge(jitter.width = 0.5)) +
  theme_bw() +
  ylim(0, 15) +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.1) +
  scale_fill_manual(values = c("goldenrod2", "aquamarine4")) +
  scale_color_manual(values = c("goldenrod2", "aquamarine4")) +
  ylab("Sporulation efficiency %") +
  theme(axis.title.x=element_text(size=12, face='bold'),
        axis.title.y=element_text(size=12, face='bold'),
        panel.background = element_rect(fill = 'gray99'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

##STAT

aov_sporeff <- aov(spores~media, spore_eff_2)
summary(aov_sporeff)
#TukeyHSD(aov_sporeff)
t.test(spores~media, data = spore_eff_2)

#Welch Two Sample t-test

#data:  spores by media
#t = -10.586, df = 48.613, p-value = 3.245e-14
#alternative hypothesis: true difference in means between group BHI and group Butyrate is not equal to 0
#95 percent confidence interval:
#  -13.573372  -9.241443
#sample estimates:
#  mean in group BHI mean in group Butyrate 
#3.333333              14.740741 

pairwise.t.test(x=spore_eff_2$spores, g=spore_eff_2$media, "bonferroni")

#Pairwise comparisons using t tests with pooled SD 

#data:  spore_eff_2$spores and spore_eff_2$media 

#           BHI    
#Butyrate 1.4e-14

#P value adjustment method: bonferroni 


