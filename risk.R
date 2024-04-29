#### MetaGSE计算GScore ####


# 清空当前环境
rm(list = ls())
getwd()
# 设置当前工作路径
setwd("E:/本地生信项目/DKD/RNAseq/Risk/")
options(stringsAsFactors = FALSE)
suppressMessages(library(GSVA))
suppressMessages(library(data.table))
suppressMessages(library(stringr))
suppressMessages(library(dplyr))
suppressMessages(library(tibble))
suppressMessages(library(openxlsx))


#### 1. 读取数据 ####

##读取免疫细胞Marker基因文件
geneset <- read.xlsx("E:/本地生信项目/DKD/RNAseq/Risk/GScore.xlsx")

##表达矩阵
expr <- read.csv("E:/本地生信项目/DKD/RNAseq/Train_Cohort/exprmatrix_remove_batch_remove.csv",row.names = 1)

##分组信息
groupinfo <- read.csv("E:/本地生信项目/DKD/RNAseq/Train_Cohort/groupinfo.csv",row.names = 1)


# 构建免疫细胞Marker背景基因集(需要组合成list格式)：
geneset2 <- split(geneset$Metagene,geneset$Type)


#### 2. 通过ssgsea计算每个样本的GScore ####
res <- gsva(
  as.matrix(expr), #表达量矩阵，可以是log-CPMs、log-RPKMs或log-TPMs
  geneset2, #Marker基因list
  method = "ssgsea", #单样本基因富集分数估算方法,可选"gsva", "ssgsea", "zscore"和"plage"
  kcdf = "Gaussian",
  mx.diff = F,
  verbose   = F #是否提供每个步骤详细信息
)


resm <- t(res) %>% as.data.frame()

# 合并groupinfo
groupinfo_gscore <- cbind(groupinfo,resm)

# 判断高低
groupinfo_gscore$risk = as.vector(ifelse(groupinfo_gscore$GScore>median(groupinfo_gscore$GScore),"high","low"))
rownames(groupinfo_gscore)
# 保存数据
write.csv(groupinfo_gscore,"./groupinfo_gscore.csv")


rt = groupinfo_gscore
rt=rt[order(rt$GScore),]


"#FFBE7A"
"#FA7F6F"
"#82B0D2"


## High与Low的GScore散点图 ####

riskClass=rt[,"risk"]
lowLength=length(riskClass[riskClass=="low"])
highLength=length(riskClass[riskClass=="high"])
line=rt[,"GScore"]
line[line>10]=10


plot(line, type="p", pch=20, cex = 2.5,
     xlab="Patients (increasing GScore)", ylab="GScore",
     col=c(rep("#82B0D2",lowLength),rep("#FA7F6F",highLength)) )
abline(h=median(rt$riskscore),v=lowLength,lty=2,lwd=2.5)
abline(h=median(rt$GScore),lty=2,lwd=2.5)
legend("topleft", c("High", "Low"),bty="n",pch=19,col=c("#FA7F6F","#82B0D2"),cex=1.2)


## DKD与Normal的GScore散点图 ####

riskClass=rt[,"group_list"]
lowLength=length(riskClass[riskClass=="Normal"])
highLength=length(riskClass[riskClass=="DKD"])
line=rt[,"GScore"]
line[line>10]=10

plot(line, type="p", pch=20, cex = 2.5,
     xlab="Patients (increasing GScore)", ylab="GScore",
     col=c(rep("#82B0D2",lowLength),rep("#FFBE7A",highLength)) )
abline(h=median(rt$riskscore),v=lowLength,lty=2,lwd=2.5)
abline(h=median(rt$GScore),lty=2,lwd=2.5)
legend("topleft", c("DKD", "Normal"),bty="n",pch=19,col=c("#FFBE7A","#82B0D2"),cex=1.2)


## 热图 ####

library(ComplexHeatmap)
setwd("E:/本地生信项目/DKD/RNAseq/Risk/")

## 表达矩阵
expr <- read.csv("E:/本地生信项目/DKD/RNAseq/Train_Cohort/exprmatrix_remove_batch_remove.csv",row.names = 1)

## 分组信息
groupinfo <- groupinfo_gscore

## 读取特征基因
fea_df <- read.delim("E:/本地生信项目/DKD/RNAseq/ML_combination/Results/fea_df.txt",sep = "\t")
fea_df <- fea_df[fea_df$algorithm == "SVM+Enet[alpha=0.6]",]
fea_df

#提取仅含有候选基因的表达矩阵
expr <- expr[rownames(expr) %in% fea_df$features,]

## 进行zscore标准化
n=t(scale(t(expr)))
n[n>1]=1 #限定上限，使表达量大于1的等于1
n[n< -1]= -1 #限定下限，使表达量小于-1的等于-1


## 增加标签的图（这里用的pheatmap）
# 列注释
annotation_col = data.frame(Sample = factor(groupinfo$group_list),
                            Tissue = factor(groupinfo$tissue), 
                            Dataset = factor(groupinfo$dataset),
                            Type = factor(groupinfo$risk))
row.names(annotation_col) = rownames(groupinfo) #这一行必须有，否则会报错：Error in check.length("fill") :  'gpar' element 'fill' must not be length 0

# 下面需要完全一致
row.names(annotation_col)
colnames(expr)

table(annotation_col$Sample)
table(annotation_col$Tissue)
table(annotation_col$Dataset)
table(annotation_col$Type)

# 设置颜色
ann_colors <- list(
  Sample = c(DKD = "#2878B5", Normal = "#9AC9DB"),
  Tissue = c(glomeruli = "#54B345", tubuli = "#32B897"),
  Dataset = c(GSE104948 = "#8ECFC9", GSE104954 = "#FFBE7A", GSE47183 = "#FA7F6F", GSE47184 = "#82B0D2"),
  Type = c(high = "#C24976",low = "#469393")
)

pheatmap(n,
         scale = "none",
         name = "z-score",
         clustering_distance_rows = "correlation", #表示行聚类使用皮尔森相关系数聚类
         treeheight_row = 30, #设置行聚类树的高度
         treeheight_col = 30, #设置列聚类树的高度
         #cutree_cols = 3, #根据样品列聚类情况将热图的行方向隔开为3份
         #cutree_rows = 3, #根据样品行聚类情况将热图的行方向隔开为3份
         display_numbers = F, # 热图上显示数值
         show_colnames = FALSE, #设置列标签的显示
         show_rownames = TRUE, #设置行标签的显示
         #border="white", # 设置边框为白色
         annotation_col = annotation_col, #显示样品列的分组信息及图例
         annotation_colors = ann_colors, #使用annotation_colors参数设定样品列分组的颜色
         col = colorRampPalette(c("#82B0D2","white","#FA7F6F"))(100) #设置颜色
)


