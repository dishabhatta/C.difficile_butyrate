## Paper: SCFA and C. difficile
## Figure 5, 6, S4, S5: RNA seq data analysis
### 5a: NMDS; 5b: early log NES gsea; 5C: late log NES gsea 
### 6a: volcano plot early log; 6b: volcano plot late log; 6c: heatmap relative abundances genes early; 6d: heatmap relative abundances genes late
### S4: heatmap of log2fc
### S5: relative abundances between BHI and butyrate

setwd("/Data/figure5_6_S4_S5/")
#libraries:
library(stringi)
library(stringr)
library(tidyverse)
library(plyr)
library(dplyr)
library(ggplot2)
library(DESeq2)
library(ggrepel)
library(EnhancedVolcano)
library(vegan)
library(clusterProfiler)
library(edgeR)
library(DescTools)
library(egg)
library(cowplot)

# Figure 6a:

## gene counts modified from All.gene.counts.final.txt: remove everything except the abundances and headers and gene names
gca9205_cd630 <- read_tsv(file = "./GCA9205_featcounts/genecounts_GCA9205.txt")

geneid_name_gca9205 <- read_tsv(file = "./GCA9205_featcounts/GCA_000009205.2_ASM920v2_feature_table.txt")
gene_meta <- read_tsv(file = "./gene_meta.txt")

geneid_name_gca9205_2 <- subset(geneid_name_gca9205, !is.na(geneid_name_gca9205$product_accession))

gca9205_2 <- gca9205_cd630
gca9205_2$Geneid <- substr(gca9205_2$Geneid, 5, nchar(gca9205_2$Geneid))

## for early log
early <- subset(gene_meta, gene_meta$time == "early")
coldata_early <- early[,2:5]
rownames(coldata_early) <- early$sample

early_gene_9205 <- gca9205_2[,c(2:4, 8:10)]
rownames(early_gene_9205) <- gca9205_2$Geneid

### To see if Rownames of conditions need to match column names of early_gene, so verify this with the statement below
all(rownames(coldata_early) == colnames(early_gene_9205))

###Set up input matrix for DESeq2. This is where the experimental design is specified
dds_early_gca9205 <- DESeqDataSetFromMatrix(countData = early_gene_9205, colData = coldata_early, design = ~ condition)

###Run analysis
dds_early_gca92051 <- DESeq(dds_early_gca9205)

###Create and view table of results. Adjusted p-value will be the things to look at.
res_early_gca9205 <- results(dds_early_gca92051,tidy = T)
#View(res)
#write.csv(res_early_gca9205_2, file="./GCA9205_Early_DESeq_output_table.csv")

#res_early_gca9205_2 <- read_csv(file = "./GCA9205_Early_DESeq_output_table.csv")
### inner join with prokka tsv to identify COG assignments
res_early_gca9205_2 <- inner_join(res_early_gca9205, geneid_name_gca9205_2, by= c("row" = "product_accession") )

## plot with enhanced volcano
EnhancedVolcano(res_early_gca9205_2,
                lab = res_early_gca9205_2$name,
                pCutoff = 0.000010, pointSize = 2, title= "Early phase GCA9205", titleLabSize = 6,
                subtitleLabSize = 6, legendLabSize = 6,
                labSize = 3, col = c('#509BCD', '#A9A0E9', '#4E3994', '#000861'),
                x = 'log2FoldChange',
                y = 'padj', axisLabSize = 6)

## plot with ggplot
res_early_gca9205_2$diffexpressed <- "NO"
res_early_gca9205_2$diffexpressed[res_early_gca9205_2$log2FoldChange > 1.0 & res_early_gca9205_2$padj < 0.0000010] <- "UP"
res_early_gca9205_2$diffexpressed[res_early_gca9205_2$log2FoldChange < -1.0 & res_early_gca9205_2$padj < 0.0000010] <- "DOWN"
res_early_gca9205_2$label[res_early_gca9205_2$diffexpressed == "UP"] <- res_early_gca9205_2$symbol[res_early_gca9205_2$diffexpressed == "UP"]
res_early_gca9205_2$label[res_early_gca9205_2$diffexpressed == "DOWN"] <- res_early_gca9205_2$symbol[res_early_gca9205_2$diffexpressed == "DOWN"]

up_early <- levels(factor(res_early_gca9205_2$name[res_early_gca9205_2$diffexpressed == "UP"]))
down_early <- levels(factor(res_early_gca9205_2$name[res_early_gca9205_2$diffexpressed == "DOWN"]))
ggplot(res_early_gca9205_2, aes(x=log2FoldChange, y=-log10(padj), color=diffexpressed, alpha=0.5)) + geom_point() + #geom_text(aes(label=label), size=2) +
  theme_bw() +
  theme(text = element_text(size=10), legend.position = "none") +
  ylim(-5,35) + xlim(-10,15) +
  geom_hline(yintercept=-log10(0.0000010), col="brown") +
  geom_vline(xintercept = c(-1,1), col="navy") +
  scale_color_manual(values = c("#509BCD", "#A9A0E9","#000861"))



# Figure 6B:

# for Late log

late <- subset(gene_meta, gene_meta$time == "late")
coldata_late <- late[,2:5]
rownames(coldata_late) <- late$sample

late_gene_gca9205 <- gca9205_2[,c(5:7, 11:13)]
rownames(late_gene_gca9205) <- gca9205_2$Geneid

# To see if Rownames of conditions need to match column names of late_gene, so verify this with the statement below
all(rownames(coldata_late) == colnames(late_gene_gca9205))

##Set up input matrix for DESeq2. This is where the experimental design is specified
dds_lab_late_gca9205 <- DESeqDataSetFromMatrix(countData = late_gene_gca9205, colData = coldata_late, design = ~ condition)

## Run analysis
dds_lab_late_gca92051 <- DESeq(dds_lab_late_gca9205)

##Create and view table of results. Adjusted p-value will be the things to look at.
res_lab_late_gca9205 <- results(dds_lab_late_gca92051,tidy = T)

res_lab_late_gca9205_2 <- inner_join(res_lab_late_gca9205, geneid_name_gca9205_2, by= c("row" = "product_accession") )
#write.csv(res_lab_late_gca9205_2, file="./GCA9205_Late_DESeq_output_table.csv")

#res_lab_late_gca9205_2 <- read_csv(file = "./GCA9205_Late_DESeq_output_table.csv")
## plot with enhanced volcano

EnhancedVolcano(res_lab_late_gca9205_2,
                lab = res_lab_late_gca9205_2$name,
                pCutoff = 0.01, pointSize = 2, title= "Late phase GCA9205", titleLabSize = 6,
                labSize = 3, subtitleLabSize = 6, legendLabSize = 6, col = c('#509BCD', '#A9A0E9', '#4E3994', '#000861'),
                x = 'log2FoldChange', axisLabSize = 6,
                y = 'padj')

res_lab_late_gca9205_2$diffexpressed <- "NO"
res_lab_late_gca9205_2$diffexpressed[res_lab_late_gca9205_2$log2FoldChange > 1.0 & res_lab_late_gca9205_2$padj < 0.0000010] <- "UP"
res_lab_late_gca9205_2$diffexpressed[res_lab_late_gca9205_2$log2FoldChange < -1.0 & res_lab_late_gca9205_2$padj < 0.0000010] <- "DOWN"
res_lab_late_gca9205_2$label[res_lab_late_gca9205_2$diffexpressed == "UP"] <- res_lab_late_gca9205_2$symbol[res_lab_late_gca9205_2$diffexpressed == "UP"]
res_lab_late_gca9205_2$label[res_lab_late_gca9205_2$diffexpressed == "DOWN"] <- res_lab_late_gca9205_2$symbol[res_lab_late_gca9205_2$diffexpressed == "DOWN"]

up_late <- levels(factor(res_lab_late_gca9205_2$name[res_lab_late_gca9205_2$diffexpressed == "UP"]))
down_late <- levels(factor(res_lab_late_gca9205_2$name[res_lab_late_gca9205_2$diffexpressed == "DOWN"]))

intersect(up_early, up_late)
setdiff(up_early, up_late)
ggplot(res_lab_late_gca9205_2, aes(x=log2FoldChange, y=-log10(padj), color=diffexpressed, alpha=0.5)) + geom_point() + #geom_text(aes(label=label), size=2) +
  theme_bw() +
  theme(text = element_text(size=10), legend.position = "none", panel.background = element_rect()) +
  ylim(-5,35) + xlim(-10,15) +
  geom_hline(yintercept=-log10(0.0000010), col="brown") +
  geom_vline(xintercept = c(-1,1), col="navy") +
  scale_color_manual(values = c("#509BCD", "#A9A0E9","#000861"))

# Figure 6c

list_siggenes_early_gca9205 <- subset(res_early_gca9205_2, res_early_gca9205_2$log2FoldChange > 1.0 & res_early_gca9205_2$padj < 0.0000010 | res_early_gca9205_2$log2FoldChange < -1.0 & res_early_gca9205_2$padj < 0.0000010)
list_siggenes_early_gca9205_2 <- c((list_siggenes_early_gca9205$row))
sig_early_up <- c(list_siggenes_early_gca9205$row[list_siggenes_early_gca9205$diffexpressed == "UP"])
sig_early_down <- c(list_siggenes_early_gca9205$row[list_siggenes_early_gca9205$diffexpressed == "DOWN"])
ntd_early_gca9205 <- normTransform(dds_early_gca92051)
ntd_assay_early_gca9205 <- assay(ntd_early_gca9205)
ntd_assay_early_gca9205_sig <- as.data.frame(ntd_assay_early_gca9205)
ntd_assay_early_gca9205_sig$geneid <- rownames(ntd_assay_early_gca9205)

ntd_early_up <- subset(ntd_assay_early_gca9205_sig, ntd_assay_early_gca9205_sig$geneid %in% sig_early_up)
ntd_meta_early_up <- inner_join(ntd_early_up, geneid_name_gca9205_2, by = c("geneid" = "product_accession"))
ntd_meta_early_up2 <- inner_join(ntd_meta_early_up, egg_short, by = c("geneid" = "query"))
cog.def <- read_tsv("./cogcatdef.txt")
ntd_meta_early_up3 <- inner_join(ntd_meta_early_up2, cog.def, by = c("COG_category" = "cogcategory"))
meta_heat_early_up <- ntd_meta_early_up2 %>% pivot_longer(cols = !c(7:27), names_to = "condition", values_to = "relab") %>% arrange(relab)

gg1 <- ggplot(meta_heat_early_up, aes(x=condition, y=name, fill=relab)) +geom_tile(width = 1, height = 1) +
  theme_bw() +
  theme(text = element_text(size = 4), axis.text.x = element_text(angle = 45, hjust=1)) +
  coord_fixed() +
  scale_fill_gradient(low = "white", high = "#0c2759") +
  scale_y_discrete(limits = unique(meta_heat_early_up$name))

ntd_early_down <- subset(ntd_assay_early_gca9205_sig, ntd_assay_early_gca9205_sig$geneid %in% sig_early_down)
ntd_meta_early_down <- inner_join(ntd_early_down, geneid_name_gca9205_2, by = c("geneid" = "product_accession"))
meta_heat_early_down <- ntd_meta_early_down %>% pivot_longer(cols = !c(7:26), names_to = "condition", values_to = "relab") %>% arrange(relab)

gg2 <- ggplot(meta_heat_early_down, aes(x=condition, y=name, fill=relab)) +geom_tile(width = 1, height=1) +
  theme_bw() +
  theme(text = element_text(size = 4), axis.text.x = element_text(angle = 45, hjust=1)) +
  coord_fixed() +
  scale_fill_gradient(low = "white", high = "#0c2759") +
  scale_y_discrete(limits = unique(meta_heat_early_down$name))

p1 <- egg::set_panel_size( gg1, height=unit(5, "cm"),
                           width=unit(10, "cm") )
p2 <- egg::set_panel_size( gg2, height=unit(1.5, "cm"),
                           width=unit(1.1299, "cm") )
cowplot::plot_grid( gg1, p2, ncol = 1 )


# Figure 6D

list_siggenes_late_gca9205 <- subset(res_lab_late_gca9205_2, res_lab_late_gca9205_2$log2FoldChange > 1.0 & res_lab_late_gca9205_2$padj < 0.0000010| res_lab_late_gca9205_2$log2FoldChange < -1.0 & res_lab_late_gca9205_2$padj < 0.0000010)
list_siggenes_late_gca9205_2 <- c((list_siggenes_late_gca9205$row))
sig_late_up <- c(list_siggenes_late_gca9205$row[list_siggenes_late_gca9205$diffexpressed == "UP"])
sig_late_down <- c(list_siggenes_late_gca9205$row[list_siggenes_late_gca9205$diffexpressed == "DOWN"])
ntd_late_gca9205 <- normTransform(dds_lab_late_gca92051)
ntd_assay_late_gca9205 <- assay(ntd_late_gca9205)
ntd_assay_late_gca9205_sig <- as.data.frame(ntd_assay_late_gca9205)
ntd_assay_late_gca9205_sig$geneid <- rownames(ntd_assay_late_gca9205)
#ntd_assay_late_gca9205_sig <- subset(ntd_assay_late_gca9205_sig, ntd_assay_late_gca9205_sig$geneid %in% list_siggenes_late_gca9205_2)

meta_ntd_late_gca9205_sig <- inner_join(ntd_assay_late_gca9205_sig, geneid_name_gca9205_2, by = c("geneid" = "product_accession"))
meta_heat_late_gca9205 <- meta_ntd_late_gca9205_sig %>% pivot_longer(cols = !c(7:26), names_to = "condition", values_to = "relab") %>% arrange((relab)) 

#meta_heat_late_100 <- subset(meta_heat_late, meta_heat_late$relab >10)

ggplot(meta_heat_late_gca9205, aes(x=condition, y=name, fill=relab)) +geom_tile() +
  theme_bw() +
  theme(text = element_text(size = 4), axis.text.x = element_text(angle = 45, hjust=1)) +
  coord_fixed() +
  scale_fill_gradient(low = "white", high = "#0c2759") +
  scale_y_discrete(limits = unique(meta_heat_late_gca9205$name))

ntd_late_up <- subset(ntd_assay_late_gca9205_sig, ntd_assay_late_gca9205_sig$geneid %in% sig_late_up)
top <- ntd_late_up[order(rowSums(ntd_late_up[,1:6]), decreasing = TRUE)[1:50], ]
ntd_meta_late_up <- inner_join(top, geneid_name_gca9205_2, by = c("geneid" = "product_accession"))
meta_heat_late_up <- ntd_meta_late_up %>% pivot_longer(cols = !c(7:26), names_to = "condition", values_to = "relab") %>% arrange(relab)

gg3 <- ggplot(meta_heat_late_up, aes(x=condition, y=name, fill=relab)) +geom_tile() +
  theme_bw() +
  theme(text = element_text(size = 4), axis.text.x = element_text(angle = 45, hjust=1)) +
  coord_fixed() +
  scale_fill_gradient(low = "white", high = "#0c2759") +
  scale_y_discrete(limits = unique(meta_heat_late_up$name))

ntd_late_down <- subset(ntd_assay_late_gca9205_sig, ntd_assay_late_gca9205_sig$geneid %in% sig_late_down)
ntd_meta_late_down <- inner_join(ntd_late_down, geneid_name_gca9205_2, by = c("geneid" = "product_accession"))
meta_heat_late_down <- ntd_meta_late_down %>% pivot_longer(cols = !c(7:26), names_to = "condition", values_to = "relab") %>% arrange(relab)

gg4 <- ggplot(meta_heat_late_down, aes(x=condition, y=name, fill=relab)) +geom_tile() +
  theme_bw() +
  theme(text = element_text(size = 4), axis.text.x = element_text(angle = 45, hjust=1)) +
  coord_fixed() +
  scale_fill_gradient(low = "white", high = "#0c2759") +
  scale_y_discrete(limits = unique(meta_heat_late_down$name))

p3 <- egg::set_panel_size( gg3, height=unit(5, "cm"),
                           width=unit(10, "cm") )
p4 <- egg::set_panel_size( gg4, height=unit(1.5, "cm"),
                           width=unit(1.5, "cm") )
cowplot::plot_grid( gg3, p4, ncol = 1 )


#######################

# Figure 5a


full_lab <- gca9205_2
rownames(full_lab) <- gca9205_2$Geneid
full_lab <- full_lab[,-1]
rownames(full_lab) <- gca9205_2$Geneid
full_lab_1 <- as.matrix(full_lab)
full_lab_2 <- full_lab_1[rowSums(full_lab_1)>0, ]
full_lab_3 <- as.matrix(full_lab_2)
full_lab_4 <- t(full_lab_3)
dist_mat_gca9205 <- vegdist(full_lab_4, method = "bray")
gca9205_early_mds <- metaMDS(dist_mat_gca9205, distance = "bray")
gca9205_early_mds$stress
#[1] 0.03581383 --- this is the final stress overall, < 0.05 so very good
mds_scores_gca9205 <- as.data.frame(scores(gca9205_early_mds))
mds_scores_gca9205$sample <- rownames(mds_scores_gca9205)
mds_score_gca9205_merge <- inner_join(mds_scores_gca9205, gene_meta, by = "sample")

ggplot(mds_score_gca9205_merge, aes(x=NMDS1, y=NMDS2, color=condition, shape=time)) +geom_point(size = 3) + #geom_text(aes(label = sample), size=2, nudge_x = 0.3) +
  theme_bw() +
  theme(panel.grid = element_blank(), legend.position = "right") +
  #stat_ellipse() +
  scale_color_manual(values = c("#dc8c24", "aquamarine4"))


# Figure 5b:

## Over-representation of KEGG and gene set enrichment analysis of KEGG/KO values for supplementary

## the cd630 faa was annotated with EGGNOG-MAPPER for ko and KEGG
## additional annotation can be done using Blast2GO and blastp and interproscan

egg_gca9205 <- read_tsv(file = "./cd630_gca9205_emapper.txt")
egg_gca9205_2 <- egg_gca9205[,c(1,5)]
egg_gca9205_3 <- egg_gca9205_2 %>% dplyr::mutate(KEGG_ko = strsplit(KEGG_ko, ',')) %>% tidyr::unnest(cols = c(KEGG_ko))

## gsea early

kegg_res_early_gca9205 <- inner_join(res_early_gca9205_2, egg_gca9205_3, by= c("row" = "query"))
kegg_res_early_gca9205_2 <- subset(kegg_res_early_gca9205, !(kegg_res_early_gca9205$KEGG_ko=="-"))
kegg_res_early_gca9205_3 <- mutate(kegg_res_early_gca9205_2, KEGG = str_sub(KEGG_ko, 4, -1)) %>% arrange(desc(stat))
geneList_early_gca9205 <- kegg_res_early_gca9205_3$log2FoldChange
names(geneList_early_gca9205) <- kegg_res_early_gca9205_3$KEGG
geneList_early_gca9205_2 <- sort(geneList_early_gca9205, decreasing = TRUE)

kegg_gsea_early_gca9205 <- gseKEGG(geneList = geneList_early_gca9205_2, organism = "ko", keyType = "kegg", pvalueCutoff = 0.05, pAdjustMethod = "fdr")
results_gsea_early_gca9205 <- kegg_gsea_early_gca9205@result %>% arrange(NES)

ggplot(results_gsea_early_gca9205, aes(x=NES , y=Description, fill=NES)) + geom_col() +
  theme_bw() +
  theme(panel.grid = element_blank(), legend.position = "none", text = element_text(size=8)) +
  scale_y_discrete(limits = results_gsea_early_gca9205$Description) +
  scale_fill_gradient(low = "#000861", high = "#D73502")


# Figure 5C:

# gsea late

kegg_res_late_gca9205 <- inner_join(res_lab_late_gca9205_2, egg_gca9205_3, by= c("row" = "query"))
kegg_res_late_gca9205_2 <- subset(kegg_res_late_gca9205, !(kegg_res_late_gca9205$KEGG_ko=="-"))
kegg_res_late_gca9205_3 <- mutate(kegg_res_late_gca9205_2, KEGG = str_sub(KEGG_ko, 4, -1))
#kegg_res_late3$padj_log10 <- -log10(kegg_res_late3$padj)
geneList_late_gca9205 <- kegg_res_late_gca9205_3$log2FoldChange
names(geneList_late_gca9205) <- kegg_res_late_gca9205_3$KEGG
geneList_late_gca9205_2 <- sort(geneList_late_gca9205, decreasing = TRUE)

kegg_gsea_late_gca9205 <- gseKEGG(geneList = geneList_late_gca9205_2, organism = "ko", keyType = "kegg",pvalueCutoff = 0.05, pAdjustMethod = "fdr")
results_gsea_late_gca9205 <- kegg_gsea_late_gca9205@result %>% arrange(NES)

ggplot(results_gsea_late_gca9205, aes(x=NES , y=Description, fill=NES)) + geom_col() +
  theme_bw() +
  theme(panel.grid = element_blank(), legend.position = "none", text = element_text(size=8)) +
  scale_y_discrete(limits = results_gsea_late_gca9205$Description) +
  scale_fill_gradient(low = "#000861", high = "#D73502")


###########################Supplementary Figures##############################

## Supplementary figure 4: boxplots with specific genes

res_early3 <- res_early_gca9205_2
res_early3$time <- "early"

res_late3 <- res_lab_late_gca9205_2
res_late3$time <- "late"

full_res <- rbind(res_early3, res_late3)
pts_genes <- full_res %>% filter(str_detect(name, "PTS"))
pts_genes2 <- subset(pts_genes, pts_genes$padj < 0.05 & pts_genes$log2FoldChange < -1 | pts_genes$padj <0.05 & pts_genes$log2FoldChange > 1)
pts_genes2$gene_type <- "pts"
ggplot(pts_genes2, aes(x= name, y=log2FoldChange, color=time)) + 
  geom_boxplot(outlier.shape = NA) + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust=1, size = 5)) +
  scale_color_manual(values = c("orangered3", "dodgerblue"))

prd_genes <- full_res %>% filter(str_detect(symbol, "prd"))

prd_genes2 <- subset(prd_genes, prd_genes$padj < 0.05 & prd_genes$log2FoldChange < -1 | prd_genes$padj <0.05 & prd_genes$log2FoldChange > 1)
prd_genes2$gene_type <- "prd"
ggplot(prd_genes2, aes(x= name, y=log2FoldChange, color=time)) + 
  geom_boxplot(outlier.shape = NA) + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust=1, size = 5)) +
  scale_color_manual(values = c("orangered3", "dodgerblue"))

rex_genes <- full_res %>% filter(str_detect(name, "Rex")) #not significant
rex_genes$gene_type <- "rex"
cody_genes <- full_res %>% filter(str_detect(name, "cod")) #not significant
cody_genes$gene_type <- "cody"
ccpA_genes <- full_res %>% filter(str_detect(symbol, "ccp")) #this one is significant
ccpA_genes$gene_type <- "ccpA"


spo0A_genes <- full_res %>% filter(str_detect(name, "spo0A") | str_detect(symbol, "spo0A"))
spo0A_genes2 <- subset(spo0A_genes, spo0A_genes$padj < 0.05 & spo0A_genes$log2FoldChange < -1 | spo0A_genes$padj <0.05 & spo0A_genes$log2FoldChange > 1)
spo0A_genes2$gene_type <- "spo0A"
ggplot(spo0A_genes2, aes(x= name, y=log2FoldChange, color=time)) + 
  geom_boxplot(outlier.shape = NA) + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust=1, size = 5)) +
  scale_color_manual(values = c("orangered3", "dodgerblue"))

grd_genes <- full_res %>% filter(str_detect(name, "lycine re") | str_detect(symbol, "grd"))
grd_genes2 <- subset(grd_genes, grd_genes$padj < 0.05 & grd_genes$log2FoldChange < -1 | grd_genes$padj <0.05 & grd_genes$log2FoldChange > 1)
grd_genes2$gene_type <- "grd"
ggplot(grd_genes2, aes(x= name, y=log2FoldChange, color=time)) + 
  geom_boxplot(outlier.shape = NA) + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust=1, size = 5)) +
  scale_color_manual(values = c("orangered3", "dodgerblue"))

tcd_genes <- full_res %>% filter(str_detect(name, "toxin") | str_detect(symbol, "tcd"))
tcd_genes2 <- subset(tcd_genes, tcd_genes$padj < 0.05 & tcd_genes$log2FoldChange < -1 | tcd_genes$padj <0.05 & tcd_genes$log2FoldChange > 1)
tcd_genes2$gene_type <- "tcd"
ggplot(tcd_genes2, aes(x= name, y=log2FoldChange, color=time)) + 
  geom_boxplot(outlier.shape = NA) + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust=1, size = 5)) +
  scale_color_manual(values = c("orangered3", "dodgerblue"))

cst_genes <- full_res %>% filter(str_detect(name, "carbon starv") | str_detect(symbol, "cst"))
cst_genes2 <- subset(cst_genes, cst_genes$padj < 0.05 & cst_genes$log2FoldChange < -1 | cst_genes$padj <0.05 & cst_genes$log2FoldChange > 1)
cst_genes2$gene_type <- "cst"
ggplot(cst_genes2, aes(x= name, y=log2FoldChange, color=time)) + 
  geom_boxplot(outlier.shape = NA) + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust=1, size = 5)) +
  scale_color_manual(values = c("orangered3", "dodgerblue"))


fullgenelist <- rbind(pts_genes2, prd_genes2, spo0A_genes2, grd_genes2, tcd_genes2, cst_genes2, cody_genes, rex_genes, ccpA_genes)
fullgenelist$color[fullgenelist$gene_type == "pts"] <- "blue"
fullgenelist$color[fullgenelist$gene_type == "prd"] <- "mediumorchid4"  
fullgenelist$color[fullgenelist$gene_type == "spo0A"] <- "forestgreen"
fullgenelist$color[fullgenelist$gene_type == "grd"] <- "lightgoldenrod"
fullgenelist$color[fullgenelist$gene_type == "tcd"] <- "firebrick1"
fullgenelist$color[fullgenelist$gene_type == "cst"] <- "peru"
fullgenelist$color[fullgenelist$gene_type == "cody"] <- "yellowgreen"
fullgenelist$color[fullgenelist$gene_type == "rex"] <- "magenta"  
fullgenelist$color[fullgenelist$gene_type == "ccpA"] <- "skyblue"
                  
fullgenelist_color <- fullgenelist[,c(20,29)] %>% distinct()
                
ggplot(fullgenelist, aes(x= time, y=name, fill=log2FoldChange)) + 
                  geom_tile() + 
                  coord_fixed() +
                  theme_bw() + 
                  #facet_grid(~gene_type, scales = "free_x") +
                  theme(axis.text.x = element_text(angle = 45, hjust=1, size = 5), axis.text.y = element_text(size=3, color = fullgenelist_color$color)) +
                  scale_fill_gradient2(low = "white", mid = "yellowgreen", high = "darkgreen", midpoint = 0) +
                  scale_y_discrete(limits = unique(fullgenelist_color$name))


## Supplementary Figure S5:


ntd_assay_early_gca9205_sig$time <- "early"
ntd_assay_late_gca9205_sig$time <- "late"
#full_ntd <- rbind(ntd_assay_early_gca9205_sig, ntd_assay_late_gca9205_sig)
#fullgenelist2 <- inner_join(fullgenelist, ntd_assay_early_gca9205_sig, by = c("row" = "geneid", "time" = "time"))
#fullgenelist2 <- fullgenelist2 %>% pivot_longer(cols = !c(1:28), names_to = "species", values_to = "relab")
#fullgenelist3 <- inner_join(fullgenelist, ntd_assay_late_gca9205_sig, by = c("row" = "geneid", "time" = "time"))
#fullgenelist3 <- fullgenelist3 %>% pivot_longer(cols = !c(1:28), names_to = "species", values_to = "relab")
#fullgenelist_ntd <- rbind(fullgenelist2, fullgenelist3)

#ggplot(fullgenelist_ntd, aes(x= name, y=relab, color=time)) + 
#  geom_boxplot(alpha =0.4, size = 0.3, aes(fill = time), outlier.shape = NA) + geom_point(position = position_jitterdodge(jitter.width = 0.3)) +
#  theme_bw() + 
#  facet_grid(~gene_type, scales = "free_x") +
#  theme(axis.text.x = element_text(angle = 45, hjust=1, size = 1), axis.text.y = element_text(size=5)) #+
#scale_color_manual(values = c("orangered3", "dodgerblue"))

ntd_assay_early_gca9205_sig2 <- ntd_assay_early_gca9205_sig %>% pivot_longer(cols = !c("geneid", "time"), values_to = "relab", names_to = "samples")
ntd_assay_late_gca9205_sig2 <- ntd_assay_late_gca9205_sig %>% pivot_longer(cols = !c("geneid", "time"), values_to = "relab", names_to = "samples")
ntd_early <- ntd_assay_early_gca9205_sig2[ntd_assay_early_gca9205_sig2$geneid %in% fullgenelist$row, ]
ntd_late <- ntd_assay_late_gca9205_sig2[ntd_assay_late_gca9205_sig2$geneid %in% fullgenelist$row, ]
fullntd <- rbind(ntd_early, ntd_late)

fullntd2 <- inner_join(fullntd, fullgenelist, by = c("geneid"="row", "time"))
fullntd2$type[(str_detect(fullntd2$samples,"BHI"))] <- "BHI"
fullntd2$type[(str_detect(fullntd2$samples,"utyrate"))] <- "Butyrate"

fullntd3 <- subset(fullntd2, !(fullntd2$gene_type == "pts"))
fullntd4 <- rbind(pts_genes2, cody_genes, rex_genes)
ntd.pts <- ntd_assay_early_gca9205_sig2[ntd_assay_early_gca9205_sig2$geneid %in% fullntd4$row, ]
ntd.pts.late <- ntd_assay_late_gca9205_sig2[ntd_assay_late_gca9205_sig2$geneid %in% fullntd4$row, ]
fullntd5 <- rbind(ntd.pts, ntd.pts.late)
fullntd6 <- inner_join(fullntd5, fullgenelist, by = c("geneid"="row", "time"))
fullntd6$type[(str_detect(fullntd6$samples,"BHI"))] <- "BHI"
fullntd6$type[(str_detect(fullntd6$samples,"utyrate"))] <- "Butyrate"

ggplot(fullntd3, aes(x= name, y=relab, color=type)) + 
  geom_boxplot(alpha =0.4, size = 0.3, aes(fill = type), outlier.shape = NA) + 
  geom_point(aes(shape = time), position = position_jitterdodge(jitter.width = 0.5)) +
  theme_bw() + 
  facet_grid(~gene_type, scales = "free_x") +
  theme(axis.text.x = element_text(angle = 45, hjust=1, size = 1), axis.text.y = element_text(size=5)) +
  scale_color_manual(values = c("goldenrod2", "aquamarine4"))

fullntd.pts <- subset(fullntd6, fullntd6$gene_type == "pts")
aov.pts <- aov(relab ~ name + type, data = fullntd.pts)
summary(aov.pts)
pts <- TukeyHSD(aov.pts)
#DunnettTest(x=fullntd.pts$relab, g=fullntd.pts$type, "")
list_pts <- split(fullntd.pts, f=fullntd.pts$name)
aov.pts1 <- aov(relab ~ type, data = list_pts[[1]])
summary(aov.pts1)
#Df Sum Sq Mean Sq F value Pr(>F)  
#type         1  29.42   29.42   8.103 0.0466 *
#  Residuals    4  14.52    3.63                 
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
aov.pts2 <- aov(relab ~ type, data = list_pts[[2]])
summary(aov.pts2)
#Df Sum Sq Mean Sq F value Pr(>F)  
#type         1  6.315   6.315   14.05   0.02 *
#  Residuals    4  1.798   0.449                 
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
aov.pts3 <- aov(relab ~ type, data = list_pts[[3]])
summary(aov.pts3) # sig *
aov.pts4 <- aov(relab ~ type, data = list_pts[[4]])
summary(aov.pts4)# sig*
aov.pts5 <- aov(relab ~ type, data = list_pts[[5]])
summary(aov.pts5)# sig *
aov.pts6 <- aov(relab ~ type, data = list_pts[[6]])
summary(aov.pts6)# not sig
aov.pts7 <- aov(relab ~ type, data = list_pts[[7]])
summary(aov.pts7)# sig *
aov.pts8 <- aov(relab ~ type, data = list_pts[[8]])
summary(aov.pts8) # no sig
aov.pts9 <- aov(relab ~ type, data = list_pts[[9]])
summary(aov.pts9)# no sig
aov.pts10 <- aov(relab ~ type, data = list_pts[[10]])
summary(aov.pts10) #no sig
aov.pts11 <- aov(relab ~ type, data = list_pts[[11]])
summary(aov.pts11) # sig *
aov.pts12 <- aov(relab ~ type, data = list_pts[[12]])
summary(aov.pts12)# **
aov.pts13 <- aov(relab ~ type, data = list_pts[[13]])
summary(aov.pts13)# not sig





ggplot(fullntd.pts, aes(x= name, y=relab, color=type)) + 
  geom_boxplot(alpha =0.4, size = 0.3, aes(fill = type), outlier.shape = NA) + 
  geom_point(aes(shape = time), position = position_jitterdodge(jitter.width = 0.5)) +
  theme_bw() + 
  facet_grid(~gene_type) +
  theme(axis.text.x = element_text(angle = 45, hjust=1, size = 1), axis.text.y = element_text(size=5)) +
  scale_color_manual(values =c("goldenrod2", "aquamarine4"))

ggplot(fullntd6, aes(x= name, y=relab, color=type)) + 
  geom_boxplot(alpha =0.4, size = 0.3, aes(fill = type), outlier.shape = NA) + 
  geom_point(aes(shape = time), position = position_jitterdodge(jitter.width = 0.5)) +
  theme_bw() + 
  facet_grid(~gene_type, scales = "free_x") +
  theme(axis.text.x = element_text(angle = 45, hjust=1, size = 1), axis.text.y = element_text(size=5)) +
  scale_color_manual(values =c("goldenrod2", "aquamarine4"))


fullntd.cody <- subset(fullntd6, fullntd6$gene_type == "cody")
aov.cody <- aov(relab ~ type, data = fullntd.cody)
summary(aov.cody)
cody <- TukeyHSD(aov.cody)

fullntd.rex <- subset(fullntd6, fullntd6$gene_type == "rex")
aov.rex <- aov(relab ~ type, data = fullntd.rex)
summary(aov.rex)
cody <- TukeyHSD(aov.rex)

fullntd.spo0A <- subset(fullntd3, fullntd3$gene_type == "spo0A")
aov.spo0A <- aov(relab~type, data=fullntd.spo0A)
summary(aov.spo0A)

fullntd.ccpA <- subset(fullntd6, fullntd6$gene_type == "ccpA")
aov.ccpA <- aov(relab ~ type, data = fullntd.ccpA)
summary(aov.ccpA)
#Df Sum Sq Mean Sq F value  Pr(>F)   
#type         1  4.782   4.782    18.5 0.00156 **
#  Residuals   10  2.585   0.259                   
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


fullntd.prd <- subset(fullntd3, fullntd3$gene_type == "prd")
aov.prd <- aov(relab ~ name + type, data = fullntd.prd)
summary(aov.prd)
prd <- TukeyHSD(aov.prd)
pairwise.t.test(x=fullntd.prd$relab, g=fullntd.prd$type, "bonferroni")
list_prd <- split(fullntd.prd, f=fullntd.prd$name)
aov.prd1 <- aov(relab ~ type, data = list_prd[[1]])
summary(aov.prd1)
aov.prd2 <- aov(relab ~ type, data = list_prd[[2]])
summary(aov.prd2)
aov.prd3 <- aov(relab ~ type, data = list_prd[[3]])
summary(aov.prd3)
aov.prd4 <- aov(relab ~ type, data = list_prd[[4]])
summary(aov.prd4)
aov.prd5 <- aov(relab ~ type, data = list_prd[[5]])
summary(aov.prd5)
#Df Sum Sq Mean Sq F value Pr(>F)   
#type         1 2.3467  2.3467   38.66 0.0034 **
#  Residuals    4 0.2428  0.0607                  
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

fullntd.grd <- subset(fullntd3, fullntd3$gene_type == "grd")
list_grd <- split(fullntd.grd, f=fullntd.grd$name)
aov.grd1 <- aov(relab ~ type, data = list_grd[[1]])# sig
summary(aov.grd1)
aov.grd2 <- aov(relab ~ type, data = list_grd[[2]])
summary(aov.grd2)
aov.grd3 <- aov(relab ~ type, data = list_grd[[3]])
summary(aov.grd3)
aov.grd4 <- aov(relab ~ type, data = list_grd[[4]])
summary(aov.grd4)
aov.grd5 <- aov(relab ~ type, data = list_grd[[5]])# sig
summary(aov.grd5)

fullntd.tcd <- subset(fullntd3, fullntd3$gene_type == "tcd")
list_tcd <- split(fullntd.tcd, f=fullntd.tcd$name)
aov.tcd1 <- aov(relab ~ type, data = list_tcd[[1]])
summary(aov.tcd1)
#Df Sum Sq Mean Sq F value   Pr(>F)    
#type         1  34.12   34.12   23.18 0.000708 ***
#  Residuals   10  14.72    1.47                     
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
aov.tcd2 <- aov(relab ~ type, data = list_tcd[[2]])
summary(aov.tcd2)
#Df Sum Sq Mean Sq F value Pr(>F)  
#type         1 2.0079  2.0079   19.07  0.012 *
#  Residuals    4 0.4212  0.1053                 
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

fullntd.cst <- subset(fullntd3, fullntd3$gene_type == "cst")
list_cst <- split(fullntd.cst, f=fullntd.cst$name)
aov.cst1 <- aov(relab ~ type, data = list_cst[[1]])
summary(aov.cst1)
#Df Sum Sq Mean Sq F value Pr(>F)  
#type         1 12.641  12.641   13.64  0.021 *
#  Residuals    4  3.707   0.927                 
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
aov.cst2 <- aov(relab ~ type, data = list_cst[[2]])
summary(aov.cst2)
#Df Sum Sq Mean Sq F value  Pr(>F)   
#type         1  13.29  13.289   11.22 0.00737 **
#  Residuals   10  11.84   1.184                   
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


