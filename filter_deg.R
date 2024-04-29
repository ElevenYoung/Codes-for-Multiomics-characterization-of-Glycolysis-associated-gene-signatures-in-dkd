#### 筛选差异基因并提取相应的表达矩阵 ####


#### 1. 提取差异基因 ####

#我们把logFC 0.5~1的都筛选了一遍

rm(list = ls())
setwd("E:/本地生信项目/DKD/RNAseq/Train_Cohort/deg")
options(stringsAsFactors = FALSE)
suppressMessages(library(readxl))
suppressMessages(library(dplyr))
suppressMessages(library(data.table))
suppressMessages(library(GEOquery))
suppressMessages(library(limma)) 
suppressMessages(library(affy))
suppressMessages(library(stringr))

##读取差异表达分析后结果
deg <- read.delim("deg_all.txt", row.names = 1)


## 设置筛选阈值

## logFC0.5
logFC_filter = 0.5
adj.P.Val_filter = 0.05

## 给差异基因打标签
# logFC > 0.5且 adj.P.Val < 0.05认为是上调基因
# logFC < -0.5且 adj.P.Val < 0.05认为是下调基因
# 其它为非差异基因
nrDEG_limma <- deg %>% 
  mutate(deg = case_when(logFC > logFC_filter & adj.P.Val < adj.P.Val_filter ~ "Up",
                         abs(logFC) < logFC_filter | adj.P.Val > adj.P.Val_filter ~ "None",
                         logFC < -logFC_filter & adj.P.Val < adj.P.Val_filter ~ "Down"))

## 查看差异基因情况
table(nrDEG_limma$deg)
# 101下调
# 195上调

## 保存数据
write.csv(nrDEG_limma,"deg_logFC0.5_filtered.csv")


## logFC1
logFC_filter = 1
adj.P.Val_filter = 0.05

## 给差异基因打标签
# logFC > 1且 adj.P.Val < 0.05认为是上调基因
# logFC < -1且 adj.P.Val < 0.05认为是下调基因
# 其它为非差异基因
nrDEG_limma <- deg %>% 
  mutate(deg = case_when(logFC > logFC_filter & adj.P.Val < adj.P.Val_filter ~ "Up",
                         abs(logFC) < logFC_filter | adj.P.Val > adj.P.Val_filter ~ "None",
                         logFC < -logFC_filter & adj.P.Val < adj.P.Val_filter ~ "Down"))

## 查看差异基因情况
table(nrDEG_limma$deg)
# 20下调
# 8上调

## 保存数据
write.csv(nrDEG_limma,"deg_logFC1_filtered.csv")


#### 2. 提取表达矩阵 ####

# 这里我们只提取logFC0.5的差异基因表达矩阵
#提取上调以及下调基因
up_0.5 <- nrDEG_limma[nrDEG_limma$deg %in% c("Up"), ]
down_0.5 <- nrDEG_limma[nrDEG_limma$deg %in% c("Down"), ]

#表达矩阵
expr <- read.csv("E:/本地生信项目/DKD/RNAseq/Train_Cohort/exprmatrix_remove_batch_remove.csv",row.names = 1)

#提取仅含有候选基因的表达矩阵
expr_filtered_up_0.5 <- expr[rownames(expr) %in% rownames(up_0.5),]
expr_filtered_down_0.5 <- expr[rownames(expr) %in% rownames(down_0.5),]

## 保存数据
write.csv(expr_filtered_up_0.5,"expr_filtered_up_0.5.csv")
write.csv(expr_filtered_down_0.5,"expr_filtered_down_0.5.csv")