#### 基因相关性分析与可视化 ####

# 我们先做一个DKD和Normal组分别的特征基因的相关性分析

rm(list = ls())
setwd("E:/本地生信项目/DKD/RNAseq/ML_combination")
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
#读取特征基因
fea_df <- read.delim("./Results/fea_df.txt",sep = "\t")
fea_df <- fea_df[fea_df$algorithm == "SVM+Enet[alpha=0.6]",]
fea_df

##表达矩阵
expr <- read.csv("E:/本地生信项目/DKD/RNAseq/ML_combination/MetaGSE_feadf_expr.csv",row.names = 1)


##分组信息
groupinfo <- read.csv("E:/本地生信项目/DKD/RNAseq/ML_combination/MetaGSE_feadf_groupinfo.csv",row.names = 1)


#### 2. 提取表达数据 ####

##将DKD与Normal的表达数据分开
DKDSample <- subset(groupinfo, Type == "DKD")
expr_DKD <- expr[,colnames(expr) %in% rownames(DKDSample)]


NormalSample <- subset(groupinfo, Type == "Normal")
expr_Normal <- expr[,colnames(expr) %in% rownames(NormalSample)]


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
# write.csv(data_cor_dkd,"./data_cor_dkd.csv")


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
# write.csv(data_cor_normal,"./data_cor_normal.csv")


#### 4. 可视化 ####
library(corrplot)
library(ggplot2)
library(ggpubr)

#读取数据
# expr_Normal <- read.csv("E:/本地生信项目/DKD/RNAseq/gene_cor/expr_Normal.csv",row.names = 1)
# expr_DKD <- read.csv("E:/本地生信项目/DKD/RNAseq/gene_cor/expr_DKD.csv",row.names = 1)

#转置数据
expr_Normal <- t(expr_Normal)
expr_DKD <- t(expr_DKD)

#分别计算相关性
data_Normal <- cor(expr_Normal, method="pearson")
data_DKD <- cor(expr_DKD, method="pearson")


### 基因相关性热图可视化 

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
               p.mat = testRes_dkd$p,sig.level = c(0.001, 0.01, 0.05),pch.cex = 1.5,insig = 'label_sig', pch.col = 'grey20', # 显示p值
               tl.col = "black", tl.cex = 0.5, tl.srt = 90
               )

p4 <- corrplot(data_Normal, method = "color", type = "lower",
               col =  colorRampPalette(c("#1D75B5","white","#ED7C72"))(100), 
               order = "hclust",hclust.method = "ward.D2",
               p.mat = testRes_Normal$p,sig.level = c(0.001, 0.01, 0.05),pch.cex = 1.5,insig = 'label_sig', pch.col = 'grey20', # 显示p值
               tl.col = "black", tl.cex = 0.5, tl.srt = 90
               )


## 我们批量做四个data_cor_dkd中最相关的两个基因的折线图
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(cowplot)
library(ggExtra)
# 读取表达量矩阵

expst <- expr_DKD
expst <- expr_Normal
expst <- expr

new_exp <- expst

# 转置
t_exp <- t(new_exp)
# 转化为data.frame
t_exp <- as.data.frame(t_exp)

# 选择目的基因

#DKD
target_gene <- 'ELF3'
target_gene <- 'CYBB'
target_gene <- 'PROM1'

#Normal
target_gene <- 'FCN1'
target_gene <- 'ELF3'
target_gene <- 'LCN2'
target_gene <- 'CD163'

# 计算每个基因与目的基因的相关性系数及p值
cor_list <- list()
for (i in colnames(t_exp)) {
  # 提取目的基因表达值
  tar <- t_exp[,target_gene]
  # 相关性检验
  cor_res <- cor(x = tar,y = t_exp[,i],method = 'pearson')
  # pvalue
  cor_pval <- cor.test(x = tar,y = t_exp[,i])$p.value
  # 合并结果
  final_res <- data.frame(tar_genename = target_gene,
                          gene_name = i,
                          cor_results = cor_res,
                          cor_pvalue = cor_pval)
  # 储存
  cor_list[[i]] <- final_res
}

# 合并结果
gene_corres <- do.call('rbind',cor_list)
# 查看结果
head(gene_corres,4)

# 筛选相关性系数> 0.7 或< -0.7 且p值小于0.05的基因
high_cor <- gene_corres %>% filter(abs(cor_results) >= 0.3 ,cor_pvalue < 0.05 )
# 保存基因名
high_corgene <- high_cor$gene_name
# 查看符合条件的基因数量


# 批量绘图
plotlst <- list()
for (i in high_corgene) {
  # 绘图
  p <- ggplot(t_exp,aes_string(x = target_gene,y = i)) +
    geom_point(size = 2,color = '#F19765') +
    theme_bw() +
    # 主题细节调整
    theme(axis.title = element_text(size = 16),
          axis.text = element_text(size = 14),
          axis.ticks.length = unit(0.25,'cm'),
          axis.ticks = element_line(size = 1),
          panel.border = element_rect(size = 1.5),
          panel.grid = element_blank()
    ) +
    # 添加回归线
    geom_smooth(method = 'lm',se = T,color = 'red',size = 1.5,fill = 'grey') +
    # 添加相关性系数及p值
    stat_cor(method = "pearson",digits = 3,size=6)
  
  # 添加边际柱形密度图
  p1 <- ggMarginal(p,type = "densigram",
                   xparams = list(binwidth = 0.1, fill = "#B3E283",size = .7),
                   yparams = list(binwidth = 0.1, fill = "#8AB6D6",size = .7))
  
  # 保存在list里
  plotlst[[i]] <- p1
}

# 拼图，4个一行
allplot <- plot_grid(plotlist = plotlst,ncol = 4 ,align = "hv")
allplot





# 最后，我们做一个所有样本的特征基因的相关性分析
# 发现结果没什么用，就不保存了

rm(list = ls())
setwd("E:/本地生信项目/DKD/RNAseq/ML_combination")
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
fea_df <- read.delim("./Results/fea_df.txt",sep = "\t")
fea_df <- fea_df[fea_df$algorithm == "SVM+Enet[alpha=0.6]",]
fea_df

##表达矩阵
expr <- read.csv("E:/本地生信项目/DKD/RNAseq/ML_combination/MetaGSE_feadf_expr.csv",row.names = 1)


##分组信息
groupinfo <- read.csv("E:/本地生信项目/DKD/RNAseq/ML_combination/MetaGSE_feadf_groupinfo.csv",row.names = 1)


#### 2. 计算基因相关性 ####

##首先创建空向量
gene_name1 <- c()
gene_name2 <- c()
cor_r <- c()
pvalue <- c()


for (i in 1:nrow(expr)){
  print(i)
  for (r in i:nrow(expr)){
    print(r)
    g1 = rownames(expr)[i]
    g2 = rownames(expr)[r]
    c_r = cor(as.numeric(expr[i,]),as.numeric(expr[r,]),method="pearson")
    p = cor.test(as.numeric(expr[i,]),as.numeric(expr[r,]),method ="pearson")[[3]]
    ##保存每一步的数据，而不可直接以空向量作为每一步运行的结果
    gene_name1 = c(gene_name1,g1)
    gene_name2 = c(gene_name2,g2)
    cor_r = c(cor_r,c_r)
    pvalue = c(pvalue,p)
  }
}

data_cor<- data.frame(gene_name1,gene_name2,cor_r,pvalue)


#### 3. 可视化 ####
library(corrplot)
library(ggplot2)
library(ggpubr)

#转置数据
expr <- t(expr)


#计算相关性
data_all <- cor(expr, method="pearson")



### 基因相关性热图可视化 

library(RColorBrewer)
library(patchwork)


#有显著性
testRes_all = cor.mtest(data_all, method="pearson",conf.level = 0.95)

p5 <- corrplot(data_all, method = "color", type = "upper",
               col =  colorRampPalette(c("#1D75B5","white","#ED7C72"))(100),
               order = "hclust",hclust.method = "ward.D2",
               p.mat = testRes_all$p,sig.level = c(0.001, 0.01, 0.05),pch.cex = 1.5,insig = 'label_sig', pch.col = 'grey20', # 不显示p值
               tl.col = "black", tl.cex = 0.5, tl.srt = 90
)





