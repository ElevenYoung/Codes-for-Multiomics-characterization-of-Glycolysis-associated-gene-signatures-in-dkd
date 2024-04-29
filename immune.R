#### ssgsea免疫浸润分析 ####



# 清空当前环境
rm(list = ls())
getwd()
# 设置当前工作路径
setwd("E:/本地生信项目/DKD/RNAseq/IMMUNE/")
options(stringsAsFactors = FALSE)
suppressMessages(library(GSVA))
suppressMessages(library(data.table))
suppressMessages(library(stringr))
suppressMessages(library(dplyr))
suppressMessages(library(tibble))
suppressMessages(library(openxlsx))


#### 1. 读取数据 ####

##读取免疫细胞Marker基因文件
geneset <- read.xlsx("E:/本地生信教程/数据分析与可视化（谁删了自己去死）/RNA测序分析/免疫浸润/ssGSEA/mmc3.xlsx")

##表达矩阵
expr <- read.csv("E:/本地生信项目/DKD/RNAseq/Train_Cohort/exprmatrix_remove_batch_remove.csv",row.names = 1)

##分组信息
groupinfo <- read.csv("E:/本地生信项目/DKD/RNAseq/Train_Cohort/groupinfo.csv",row.names = 1)


# 构建免疫细胞Marker背景基因集(需要组合成list格式)：
geneset2 <- split(geneset$Metagene,geneset$Cell.type)


#ssGSEA分析/计算免疫细胞在不同样本的富集分数：
res <- gsva(
  as.matrix(expr), #表达量矩阵，可以是log-CPMs、log-RPKMs或log-TPMs
  geneset2, #Marker基因list
  method = "ssgsea", #单样本基因富集分数估算方法,可选"gsva", "ssgsea", "zscore"和"plage"
  kcdf = "Gaussian",
  mx.diff = F,
  verbose   = F #是否提供每个步骤详细信息
)


#Min-Max标准化：
#将每个样本中不同的免疫细胞富集分数标准化到0-1间；
resm <- res
for (i in colnames(res)) {
  resm[,i] <- (res[,i] -min(res[,i]))/(max(res[,i] )-min(res[,i] ))
}

#聚类热图展示ssGSEA结果：
library(pheatmap)
pheatmap(resm, show_colnames = F)


## 增加标签的图（这里用的pheatmap）
# 列注释
annotation_col = data.frame(Sample = factor(groupinfo$group_list),Tissue = factor(groupinfo$tissue), Dataset = factor(groupinfo$dataset))
row.names(annotation_col) = rownames(groupinfo) #这一行必须有，否则会报错：Error in check.length("fill") :  'gpar' element 'fill' must not be length 0

# 下面需要完全一致
row.names(annotation_col)
colnames(expr)

table(annotation_col$Sample)
table(annotation_col$Tissue)
table(annotation_col$Dataset)

# 设置颜色
ann_colors <- list(
  Sample = c(DKD = "#2878B5", Normal = "#9AC9DB"),
  Tissue = c(glomeruli = "#54B345", tubuli = "#32B897"),
  Dataset = c(GSE104948 = "#8ECFC9", GSE104954 = "#FFBE7A", GSE47183 = "#FA7F6F", GSE47184 = "#82B0D2")
)

library(ComplexHeatmap)
ComplexHeatmap::pheatmap(resm,
                         name = "z-score",
                         show_colnames = FALSE, #设置列标签的显示
                         annotation_col = annotation_col, #显示样品列的分组信息及图例
                         annotation_colors = ann_colors, #使用annotation_colors参数设定样品列分组的颜色
                         col = colorRampPalette(c("#2878B5","white","#C82423"))(100)
)


##计算相关性系数并绘图：
rescor <- cor(t(res))

#添加显著性标记:
resorp <- cor.mtest(rescor, conf.level = .95) #使用cor.mtest做显著性检验;

#提取p值矩阵；
p.mat <- resorp$p

library(corrplot)
corrplot(rescor,
         method = "square",
         order = "hclust",
         type = "lower",
         p.mat = resorp$p, sig.level = c(.001, .01, .05),pch.cex = 0.5,insig = 'label_sig', pch.col = 'grey20',
         tl.cex = 0.6,
         tl.srt = 90,
         col = colorRampPalette(c("#2878B5","white","#C82423"))(100),
         tl.col = "black")



## 绘制分组箱线图
View(resm)
a <- resm %>% t() %>% as.data.frame()

#分组信息
groupinfo <- read.csv("E:/本地生信项目/DKD/RNAseq/ML_combination/MetaGSE_feadf_groupinfo.csv",row.names = 1)

#查看是否一致
identical(rownames(a),rownames(groupinfo))

a$Type <- groupinfo$Type

a <- a %>% rownames_to_column("sample")

library(ggsci)
library(tidyr)
library(ggpubr)

# 将宽数据变为长数据
b <- gather(a, key = cell_type, value = Expression, -c(Type,sample))

#将value根据样本转换为百分比形式(新增一列)：
dt <- b %>%
  group_by(sample) %>%
  mutate(proportion = round(Expression/sum(Expression),3))

head(dt)

colnames(dt) <- c("sample","group","cell_type","Expression","proportion")

#使用t test或wilcox test进行两两比较(T检验为例)：
t <- t_test(group_by(dt, cell_type), proportion ~ group)
tj <- adjust_pvalue(t, method = 'fdr') #p值矫正；
tj

#根据p.adj添加显著性标记符号；
tj <- add_significance(tj, 'p.adj')
tj

#在图表中添加 p 值或者显著性标记；
lab <- add_xy_position(tj, x = 'cell_type', dodge = 0.65)

#ggpubr绘图：
p1 <- ggboxplot(dt, x = "cell_type", y = "proportion",
                fill = "group", color = "black",
                width=0.7, alpha = 0.6,
                outlier.shape = 21, outlier.size = 1.2) +
  scale_fill_manual(values = c("#2878B5","#C82423")) +
  labs(x = "cell type", y = "proportion") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90)) +
    theme(legend.key = element_rect(fill = 'transparent'), 
        legend.background = element_rect(fill = 'transparent'), 
        legend.position = "top",
        legend.title = element_blank(),
        legend.text = element_text(size = 10,margin = margin(t = 6)), 
        axis.text.x = element_text(hjust = 0.5,size = 15), 
        axis.text.y = element_text(size = 10), 
        axis.title.x = element_text(size = 10), 
        axis.title.y = element_text(size = 10), 
        axis.line = element_line(linewidth = 1),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  stat_pvalue_manual(lab, label = 'p.adj.signif', label.size=3, bracket.size=0.5, tip.length = 0.02)
p1


## 特征基因与免疫细胞相关性棒棒糖图

#首先我们计算特征基因与免疫细胞之间的相关性
#表达矩阵
expr <- read.csv("E:/本地生信项目/DKD/RNAseq/Train_Cohort/exprmatrix_remove_batch_remove.csv",row.names = 1)

#分组信息
groupinfo <- read.csv("E:/本地生信项目/DKD/RNAseq/Train_Cohort/groupinfo.csv",row.names = 1)

#读取差异表达分析后的所有结果
fea_df <- read.delim("E:/本地生信项目/DKD/RNAseq/ML_combination/Results/fea_df.txt",sep = "\t")
fea_df <- fea_df[fea_df$algorithm == "SVM+Enet[alpha=0.6]",]
fea_df

#依次提取fea_df的表达矩阵
expr_filtered1 <- expr[rownames(expr) %in% fea_df$features[1],] %>% t() %>% as.data.frame()
expr_filtered2 <- expr[rownames(expr) %in% fea_df$features[2],] %>% t() %>% as.data.frame()
expr_filtered3 <- expr[rownames(expr) %in% fea_df$features[3],] %>% t() %>% as.data.frame()
expr_filtered4 <- expr[rownames(expr) %in% fea_df$features[4],] %>% t() %>% as.data.frame()
expr_filtered5 <- expr[rownames(expr) %in% fea_df$features[5],] %>% t() %>% as.data.frame()
expr_filtered6 <- expr[rownames(expr) %in% fea_df$features[6],] %>% t() %>% as.data.frame()
expr_filtered7 <- expr[rownames(expr) %in% fea_df$features[7],] %>% t() %>% as.data.frame()
expr_filtered8 <- expr[rownames(expr) %in% fea_df$features[8],] %>% t() %>% as.data.frame()
expr_filtered9 <- expr[rownames(expr) %in% fea_df$features[9],] %>% t() %>% as.data.frame()
expr_filtered10 <- expr[rownames(expr) %in% fea_df$features[10],] %>% t() %>% as.data.frame()
expr_filtered11 <- expr[rownames(expr) %in% fea_df$features[11],] %>% t() %>% as.data.frame()
expr_filtered12 <- expr[rownames(expr) %in% fea_df$features[12],] %>% t() %>% as.data.frame()

#提取基因
gene1=colnames(expr_filtered1)[1]
gene2=colnames(expr_filtered2)[1]
gene3=colnames(expr_filtered3)[1]
gene4=colnames(expr_filtered4)[1]
gene5=colnames(expr_filtered5)[1]
gene6=colnames(expr_filtered6)[1]
gene7=colnames(expr_filtered7)[1]
gene8=colnames(expr_filtered8)[1]
gene9=colnames(expr_filtered9)[1]
gene10=colnames(expr_filtered10)[1]
gene11=colnames(expr_filtered11)[1]
gene12=colnames(expr_filtered12)[1]


#取平均值
library(limma)
library(reshape2)
library(ggpubr)
library(vioplot)
library(ggExtra)
data1=avereps(expr_filtered1)
data2=avereps(expr_filtered2)
data3=avereps(expr_filtered3)
data4=avereps(expr_filtered4)
data5=avereps(expr_filtered5)
data6=avereps(expr_filtered6)
data7=avereps(expr_filtered7)
data8=avereps(expr_filtered8)
data9=avereps(expr_filtered9)
data10=avereps(expr_filtered10)
data11=avereps(expr_filtered11)
data12=avereps(expr_filtered12)

#根据目标基因表达量对样品进行分组
data1=as.data.frame(data1)
data1$gene=ifelse(data1[,gene1]>median(data1[,gene1]), "High", "Low")

data2=as.data.frame(data2)
data2$gene=ifelse(data2[,gene2]>median(data2[,gene2]), "High", "Low")

data3=as.data.frame(data3)
data3$gene=ifelse(data3[,gene3]>median(data3[,gene3]), "High", "Low")

data4=as.data.frame(data4)
data4$gene=ifelse(data4[,gene4]>median(data4[,gene4]), "High", "Low")

data5=as.data.frame(data5)
data5$gene=ifelse(data5[,gene5]>median(data5[,gene5]), "High", "Low")

data6=as.data.frame(data6)
data6$gene=ifelse(data6[,gene6]>median(data6[,gene6]), "High", "Low")

data7=as.data.frame(data7)
data7$gene=ifelse(data7[,gene7]>median(data7[,gene7]), "High", "Low")

data8=as.data.frame(data8)
data8$gene=ifelse(data8[,gene8]>median(data8[,gene8]), "High", "Low")

data9=as.data.frame(data9)
data9$gene=ifelse(data9[,gene9]>median(data9[,gene9]), "High", "Low")

data10=as.data.frame(data10)
data10$gene=ifelse(data10[,gene10]>median(data10[,gene10]), "High", "Low")

data11=as.data.frame(data11)
data11$gene=ifelse(data11[,gene11]>median(data11[,gene11]), "High", "Low")

data12=as.data.frame(data12)
data12$gene=ifelse(data12[,gene12]>median(data12[,gene12]), "High", "Low")


#免疫浸润结果
rownames(resm)

#我们只要dkd与normal组间非ns的免疫细胞
resm1 <- resm[-c(4,9,12,13,15,20,21,22,25,28),]
#取平均值
resm1=avereps(resm1)
resm1 <- t(resm1)


#数据合并
sameSample=intersect(row.names(resm1), row.names(data1))
rt1=cbind(resm1[sameSample,,drop=F], data1[sameSample,,drop=F])

sameSample=intersect(row.names(resm1), row.names(data2))
rt2=cbind(resm1[sameSample,,drop=F], data2[sameSample,,drop=F])

sameSample=intersect(row.names(resm1), row.names(data3))
rt3=cbind(resm1[sameSample,,drop=F], data3[sameSample,,drop=F])

sameSample=intersect(row.names(resm1), row.names(data4))
rt4=cbind(resm1[sameSample,,drop=F], data4[sameSample,,drop=F])

sameSample=intersect(row.names(resm1), row.names(data5))
rt5=cbind(resm1[sameSample,,drop=F], data5[sameSample,,drop=F])

sameSample=intersect(row.names(resm1), row.names(data6))
rt6=cbind(resm1[sameSample,,drop=F], data6[sameSample,,drop=F])

sameSample=intersect(row.names(resm1), row.names(data7))
rt7=cbind(resm1[sameSample,,drop=F], data7[sameSample,,drop=F])

sameSample=intersect(row.names(resm1), row.names(data8))
rt8=cbind(resm1[sameSample,,drop=F], data8[sameSample,,drop=F])

sameSample=intersect(row.names(resm1), row.names(data9))
rt9=cbind(resm1[sameSample,,drop=F], data9[sameSample,,drop=F])

sameSample=intersect(row.names(resm1), row.names(data10))
rt10=cbind(resm1[sameSample,,drop=F], data10[sameSample,,drop=F])

sameSample=intersect(row.names(resm1), row.names(data11))
rt11=cbind(resm1[sameSample,,drop=F], data11[sameSample,,drop=F])

sameSample=intersect(row.names(resm1), row.names(data12))
rt12=cbind(resm1[sameSample,,drop=F], data12[sameSample,,drop=F])


#检查是否有基因遗漏


#把数据转换成ggplot2输入文件
data1=rt1[,-(ncol(rt1)-1)]
data1=melt(data1,id.vars=c("gene"))
colnames(data1)=c("gene", "Immune", "Expression")

data2=rt2[,-(ncol(rt2)-1)]
data2=melt(data2,id.vars=c("gene"))
colnames(data2)=c("gene", "Immune", "Expression")

data3=rt3[,-(ncol(rt3)-1)]
data3=melt(data3,id.vars=c("gene"))
colnames(data3)=c("gene", "Immune", "Expression")

data4=rt4[,-(ncol(rt4)-1)]
data4=melt(data4,id.vars=c("gene"))
colnames(data4)=c("gene", "Immune", "Expression")

data5=rt5[,-(ncol(rt5)-1)]
data5=melt(data5,id.vars=c("gene"))
colnames(data5)=c("gene", "Immune", "Expression")

data6=rt6[,-(ncol(rt6)-1)]
data6=melt(data6,id.vars=c("gene"))
colnames(data6)=c("gene", "Immune", "Expression")

data7=rt7[,-(ncol(rt7)-1)]
data7=melt(data7,id.vars=c("gene"))
colnames(data7)=c("gene", "Immune", "Expression")

data8=rt8[,-(ncol(rt8)-1)]
data8=melt(data8,id.vars=c("gene"))
colnames(data8)=c("gene", "Immune", "Expression")

data9=rt9[,-(ncol(rt9)-1)]
data9=melt(data9,id.vars=c("gene"))
colnames(data9)=c("gene", "Immune", "Expression")

data10=rt10[,-(ncol(rt10)-1)]
data10=melt(data10,id.vars=c("gene"))
colnames(data10)=c("gene", "Immune", "Expression")

data11=rt11[,-(ncol(rt11)-1)]
data11=melt(data11,id.vars=c("gene"))
colnames(data11)=c("gene", "Immune", "Expression")

data12=rt12[,-(ncol(rt12)-1)]
data12=melt(data12,id.vars=c("gene"))
colnames(data12)=c("gene", "Immune", "Expression")


##########计算每个基因与免疫细胞的相关性##########

fea_df

outTab=data.frame()

for(i in colnames(rt12)[1:(ncol(rt12)-2)]){
  x=as.numeric(rt12[,gene12])
  y=as.numeric(rt12[,i])
  if(sd(y)==0){y[1]=0.00001}
  cor=cor.test(x, y, method="spearman")
  outVector=cbind(Cell=i, cor=cor$estimate, pvalue=cor$p.value)
  outTab=rbind(outTab,outVector)

}

#输出相关性的结果文件
write.table(outTab,file="TNFAIP8_cor.result.txt",sep="\t",row.names=F,quote=F)

fea_df


##棒棒糖图可视化

#这里只能一个一个基因做
data <- read.delim("TNFAIP8_cor.result.txt",sep="\t")

#定义圆圈颜色的函数
p.col = c('gold','pink','orange','LimeGreen','darkgreen')
fcolor = function(x,p.col){
  color = ifelse(x>0.8,p.col[1],ifelse(x>0.6,p.col[2],ifelse(x>0.4,p.col[3],
                                                             ifelse(x>0.2,p.col[4], p.col[5])
  )))
  return(color)
}



#定义设置圆圈大小的函数
p.cex = seq(2.5, 5.5, length=5)
fcex = function(x){
  x=abs(x)
  cex = ifelse(x<0.1,p.cex[1],ifelse(x<0.2,p.cex[2],ifelse(x<0.3,p.cex[3],
                                                           ifelse(x<0.4,p.cex[4],p.cex[5]))))
  return(cex)
}


#根据pvalue定义圆圈的颜色
points.color = fcolor(x=data$pvalue,p.col=p.col)
data$points.color = points.color

#根据相关系数定义圆圈的大小
points.cex = fcex(x=data$cor)
data$points.cex = points.cex
data=data[order(data$cor),]

data %>% 
  ggplot(aes(cor,forcats::fct_reorder(Cell,cor)))+
  geom_segment(aes(xend=0,yend=Cell))+
  geom_point(aes(col=pvalue,size=abs(cor)))+
  scale_colour_gradientn(colours=c("#C82423","#2878B5"))+
  scale_size_continuous(range =c(2,4))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90)) +
  geom_hline(yintercept = 0, linetype = "solid", color = "black", size = 1) +
  coord_flip()





