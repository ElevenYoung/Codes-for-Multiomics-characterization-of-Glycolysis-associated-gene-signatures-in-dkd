#### 基因相关性分析与可视化 ####

# 我们做一个DKD和Normal组分别的差异基因相关性分析

rm(list = ls())
setwd("E:/本地生信项目/DKD/RNAseq/gene_cor")
options(stringsAsFactors = FALSE)
suppressMessages(library(readxl))
suppressMessages(library(dplyr))
suppressMessages(library(data.table))
suppressMessages(library(GEOquery))
suppressMessages(library(limma)) 
suppressMessages(library(affy))
suppressMessages(library(stringr))


#### 1.读取数据 ####

##候选基因
#读取差异表达分析后的所有结果
DEG <- read.delim("E:/本地生信项目/DKD/RNAseq/Train_Cohort/deg/deg_all.txt", row.names = 1)
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
geneset <- rbind(up,down)


##表达矩阵
expr <- read.csv("E:/本地生信项目/DKD/RNAseq/Train_Cohort/exprmatrix_remove_batch_remove.csv",row.names = 1)


##分组信息
groupinfo <- read.csv("E:/本地生信项目/DKD/RNAseq/Train_Cohort/groupinfo.csv",row.names = 1)


#### 2. 提取表达数据 ####

##首先提取仅含有差异基因的表达矩阵
expr_filtered <- expr[rownames(expr) %in% rownames(geneset),]
write.csv(expr_filtered, "expr_filtered_deg.csv")

##然后将DKD与Normal的表达数据分开
DKDSample <- subset(groupinfo, group_list == "DKD")
expr_DKD <- expr_filtered[,colnames(expr_filtered) %in% rownames(DKDSample)]
write.csv(expr_DKD, "expr_DKD.csv")

NormalSample <- subset(groupinfo, group_list == "Normal")
expr_Normal <- expr_filtered[,colnames(expr_filtered) %in% rownames(NormalSample)]
write.csv(expr_Normal, "expr_Normal.csv")


#### 3. 批量计算基因相关性 ####

##首先创建空向量
gene_name1 <- c()
gene_name2 <- c()
cor_r <- c()
pvalue <- c()


##然后批量计算基因相关性

#先是DKD数据
for (i in 1:nrow(expr_DKD)){
  print(i)
  for (r in i:nrow(expr_DKD)){
    print(r)
    g1 = rownames(expr_DKD)[i]
    g2 = rownames(expr_DKD)[r]
    c_r = cor(as.numeric(expr_DKD[i,]),as.numeric(expr_DKD[r,]),method="pearson")
    p = cor.test(as.numeric(expr_DKD[i,]),as.numeric(expr_DKD[r,]),method ="pearson")[[3]]
    ##保存每一步的数据，而不可直接以空向量作为每一步运行的结果
    gene_name1 = c(gene_name1,g1)
    gene_name2 = c(gene_name2,g2)
    cor_r = c(cor_r,c_r)
    pvalue = c(pvalue,p)
  }
}

data_cor_dkd <- data.frame(gene_name1,gene_name2,cor_r,pvalue)

#然后是Normal数据
for (i in 1:nrow(expr_Normal)){
  print(i)
  for (r in i:nrow(expr_Normal)){
    print(r)
    g1 = rownames(expr_Normal)[i]
    g2 = rownames(expr_Normal)[r]
    c_r = cor(as.numeric(expr_Normal[i,]),as.numeric(expr_Normal[r,]),method="pearson")
    p = cor.test(as.numeric(expr_Normal[i,]),as.numeric(expr_Normal[r,]),method ="pearson")[[3]]
    ##保存每一步的数据，而不可直接以空向量作为每一步运行的结果
    gene_name1 = c(gene_name1,g1)
    gene_name2 = c(gene_name2,g2)
    cor_r = c(cor_r,c_r)
    pvalue = c(pvalue,p)
  }
}

data_cor_normal <- data.frame(gene_name1,gene_name2,cor_r,pvalue)


#### 4. 提取信息 ####

##这里我们要查看每个基因与其他基因正相关/负相关的数量

#首先我们先去除掉p值>0.05的数据
data_cor_normal <- data_cor_normal %>% filter(pvalue < 0.05)
data_cor_dkd <- data_cor_dkd %>% filter(pvalue < 0.05)


#dkd
dkd_result <- data_cor_dkd %>%
  group_by(gene_name1) %>%
  summarise(
    positive_count = sum(cor_r > 0, na.rm = TRUE),
    negative_count = sum(cor_r < 0, na.rm = TRUE),
    zero_count = sum(cor_r == 0, na.rm = TRUE)
)

#normal
normal_result <- data_cor_normal %>%
  group_by(gene_name1) %>%
  summarise(
    positive_count = sum(cor_r > 0, na.rm = TRUE),
    negative_count = sum(cor_r < 0, na.rm = TRUE),
    zero_count = sum(cor_r == 0, na.rm = TRUE)
)



#### 5. 保存数据 ###

#dkd
write.csv(dkd_result, "dkd_cor_count_result.csv")

#normal
write.csv(normal_result, "normal_cor_count_result.csv")


#### 6. 可视化 ####
library(corrplot)
library(ggplot2)
library(ggpubr)

#读取数据
expr_Normal <- read.csv("E:/本地生信项目/DKD/RNAseq/gene_cor/expr_Normal.csv",row.names = 1)
expr_DKD <- read.csv("E:/本地生信项目/DKD/RNAseq/gene_cor/expr_DKD.csv",row.names = 1)

#转置数据
expr_Normal <- t(expr_Normal)
expr_DKD <- t(expr_DKD)

#分别计算相关性
data_Normal <- cor(expr_Normal, method="pearson")
data_DKD <- cor(expr_DKD, method="pearson")


### 基因相关性热图可视化 

#我们先看看最基本的情况
corrplot(data_DKD)
corrplot(data_Normal)


library(RColorBrewer)
library(patchwork)

#没有显著性
p1 <- corrplot(data_DKD, method = "square", type = "upper",
               col =  colorRampPalette(c("#1D75B5","white","#ED7C72"))(100),
               order = "hclust",hclust.method = "ward.D2",
               tl.col = "black", tl.cex = 0.8, tl.srt = 90)

p2 <- corrplot(data_Normal, method = "square", type = "lower",
               col =  colorRampPalette(c("#1D75B5","white","#ED7C72"))(100), 
               order = "hclust",hclust.method = "ward.D2",
               tl.col = "black", tl.cex = 0.8, tl.srt = 90)

#有显著性
testRes_Normal = cor.mtest(data_Normal, method="pearson",conf.level = 0.95)
testRes_dkd = cor.mtest(data_DKD, method="pearson",conf.level = 0.95)

p3 <- corrplot(data_DKD, method = "color", type = "upper",
               col =  colorRampPalette(c("#1D75B5","white","#ED7C72"))(100),
               order = "hclust",hclust.method = "ward.D2",
               #p.mat = testRes_dkd$p,sig.level = c(0.001, 0.01, 0.05),pch.cex = 0.3,insig = 'label_sig', pch.col = 'grey20', # 不显示p值
               #tl.col = "black", tl.cex = 0.5, tl.srt = 90
               tl.pos = 'n')#不显示基因名

p4 <- corrplot(data_Normal, method = "color", type = "lower",
               col =  colorRampPalette(c("#1D75B5","white","#ED7C72"))(100), 
               order = "hclust",hclust.method = "ward.D2",
               #p.mat = testRes_Normal$p,sig.level = c(0.001, 0.01, 0.05),pch.cex = 0.3,insig = 'label_sig', pch.col = 'grey20', # 不显示p值
               #tl.col = "black", tl.cex = 0.5, tl.srt = 90
               tl.pos = 'n')#不显示基因名



