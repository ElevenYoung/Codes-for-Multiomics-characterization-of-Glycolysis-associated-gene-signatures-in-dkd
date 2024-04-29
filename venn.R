#### 提取差异表达基因，WGCNA特征基因以及与糖酵解基因的交集 ####

rm(list = ls())
setwd("E:/本地生信项目/DKD/RNAseq")
options(stringsAsFactors = FALSE)
suppressMessages(library(readxl))
suppressMessages(library(dplyr))
suppressMessages(library(data.table))
suppressMessages(library(tidyverse))
suppressMessages(library(stringr))


#### 1. 读取数据 ####

##糖酵解基因集
gly_geneset <- read.csv("E:/本地生信项目/糖酵解基因集/gly_geneset.csv", row.names = 1)
# 4200


##差异基因集
#logFC0.5
deg_filtered_0.5 <- read.csv("E:/本地生信项目/DKD/RNAseq/Train_Cohort/deg/deg_logFC0.5_filtered.csv",row.names = 1)
table(deg_filtered_0.5$deg)
#提取上调或下调基因
filtered_df_0.5 <- deg_filtered_0.5[deg_filtered_0.5$deg %in% c("Up", "Down"), ]
# 296


##WGCNA的turquoise模块的基因
WGCNA <- read.csv("E:/本地生信项目/DKD/RNAseq/WGCNA/turquoise_gene_filtered.csv", row.names = 1)
# 370


#### 2. 提取相同的基因 ####

## 先看logFC0.5的情况
sameSample=intersect(rownames(filtered_df_0.5), gly_geneset$gene_symbol)
sameSample=as.data.frame(sameSample) # 90
## 然后看WGCNA的情况
sameSample=intersect(sameSample$sameSample, rownames(WGCNA))
sameSample=as.data.frame(sameSample) # 45

## 更改列名
colnames(sameSample) <- c("gene_symbol")

## 提取表达数据
#表达矩阵
expr <- read.csv("E:/本地生信项目/DKD/RNAseq/Train_Cohort/exprmatrix_remove_batch_remove.csv",row.names = 1)
#提取仅含有候选基因的表达矩阵
expr_filtered <- expr[rownames(expr) %in% sameSample$gene_symbol,]


## 保存数据
write.csv(sameSample,"shared_gene_symbol_45.csv")
write.csv(expr_filtered,"shared_gene_expr_45.csv")


#### venn图可视化 ####

library(ggvenn)
library(tidyverse)
library(ggtext)
library(magrittr)
library(ggpubr)
library(cowplot)

##糖酵解基因集
gly_geneset <- read.csv("E:/本地生信项目/糖酵解基因集/gly_geneset.csv", row.names = 1)
gly_geneset$X1 <- gly_geneset$gene_symbol
# 4200


##差异基因集
#logFC0.5
DEG <- read.delim("E:/本地生信项目/DKD/RNAseq/Train_Cohort/deg/deg_all.txt",row.names = 1)
## 设置pvalue和logFC的阈值
#阈值确定：
Pvalue = 0.05
log2FC = 0.5
# 根据阈值分别为上调基因设置‘up’，下调基因设置‘Down’，无差异设置‘Stable’，保存到change列
# 这里的change列用来设置火山图点的颜色
DEG$deg = ifelse(DEG$P.Value < Pvalue & 
                   abs(DEG$logFC) >= log2FC, 
                 ifelse(DEG$logFC> log2FC ,'Up','Down'),'Stable')

head(DEG)
table(DEG$deg)
up <- filter(DEG, deg == "Up")
down <- filter(DEG, deg == "Down")
deg <- rbind(up,down)
deg$X1 <- rownames(deg)
# 301

##WGCNA的turquoise模块的特征基因
WGCNA <- read.csv("E:/本地生信项目/DKD/RNAseq/WGCNA/turquoise_gene_filtered.csv", row.names = 1)
WGCNA$X1 <- rownames(WGCNA)
# 370

p1 <- list(Glycolysis=gly_geneset$X1,DEGs=deg$X1,WGCNA=WGCNA$X1) %>% 
  ggvenn(show_percentage = T,show_elements = F,label_sep = ",",
         digits = 1,stroke_color = "white",
         fill_color = c("#FA7F6F","#82B0D2","#FFBE7A"),
         set_name_color = c("#FA7F6F","#82B0D2","#FFBE7A"))
p1
