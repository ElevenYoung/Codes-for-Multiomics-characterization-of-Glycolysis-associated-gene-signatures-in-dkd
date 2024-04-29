#### WGCNA筛选基因以及特征Module基因的富集分析与网络可视化 ####

#### 0. 环境准备 ####

# 清空当前环境
rm(list = ls())
getwd()
# 设置当前工作路径
setwd("E:/本地生信项目/DKD/RNAseq/WGCNA")
suppressMessages(library(WGCNA))
suppressMessages(library(ggplot2))
suppressMessages(library(ggrepel))
enableWGCNAThreads(8) 


#### 1. 读入表达数据和表型数据  ####

##读入表达矩阵
bulk_data <- read.csv("E:/本地生信项目/DKD/RNAseq/Train_Cohort/exprmatrix_remove_batch_remove.csv",row.names = 1)

##读入需要匹配的临床信息
datTraits <- read.csv("E:/本地生信项目/DKD/RNAseq/Train_Cohort/groupinfo.csv",row.names = 1)


#### 2. 数据预处理 ####

## 分析数据需要转置为行为样本，列为基因
bulk_data <- t(bulk_data)

# 查看数据维度
dim(bulk_data)

#92个样本，10484个基因

# 转换为数据框
bulk_data <- as.data.frame(bulk_data)

#查看行名
rownames(bulk_data)


#查看临床信息行名
rownames(datTraits)
colnames(datTraits)
table(datTraits$group_list)


#### 3. 数据检查 ####

# 我们使用WGCNA提供的goodSamplesGenes检查数据
# 检查下数据，检查数据中的缺失、离群
gsg = goodSamplesGenes(bulk_data, verbose = 3)
gsg$allOK
gsg$goodSamples

rm(gsg)
# 本示例数据数据没问题



## 做一下聚类图看看样本情况，是否有离群值
# 这里我们取了group_list一列
datTraits1 <- datTraits[,1]
datTraits1 <- as.data.frame(datTraits1)
rownames(datTraits1) <- rownames(datTraits)

## 将临床信息进行二项值转换
pheno <- binarizeCategoricalColumns(datTraits1, 
                                    dropFirstLevelVsAll = F, # 一定要选F，不然可能把第一列整没了
                                    minCount=0)
colnames(pheno)
rownames(pheno)
rownames(datTraits1)

#重新排序
#pheno的列顺序没问题，我们重命名即可
colnames(pheno) <- c("DKD","Normal")
colnames(pheno)


#### 4. 样本聚类 ####

# hclust样本聚类
sampleTree = hclust(dist(bulk_data), method = "average")
# 设置颜色
traitColors = numbers2colors(pheno,signed=T)

# 做一个聚类图，下面显示样本
pdf("1.traits_samples.pdf", width = 25, height = 10)
plotDendroAndColors(sampleTree, traitColors,
                    groupLabels = names(pheno),
                    main = "Sample dendrogram and trait heatmap",
                    addTextGuide=T,cex.colorLabels=1,
                    cex.dendroLabels = 1)
dev.off()


# 这里我们发现有离群样本:
# GSM1146337,GSM2810903,GSM1146204,GSM2811043,GSM1146027,GSM1146273,GSM2810774,GSM1146274,GSM2810775,GSM1146205,GSM1146272,GSM2810773
# 但为了数据完整性，我们还是不剔除


#### 5. 确定软阈值 ####


# 选择1-10个power值
powers = c(c(1:10), seq(from = 11, to=20, by=1))
sft = pickSoftThreshold(bulk_data,# 转置的表达矩阵
                        powerVector = powers,
                        verbose = 5)

# 作图
#pdf("2.softpower.pdf", width = 10, height = 5)
par(mfrow = c(1,2))
# 横轴是Soft threshold (power)，纵轴是无标度网络的评估参数，数值越高，网络越符合无标度特征 (non-scale)
plot(sft$fitIndices[,1], 
     -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",
     ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));

text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], labels=powers,cex=0.9,col="red");
abline(h=0.9, col="#FF3300")

# 筛选标准：R-square=0.9

plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=0.9,col="#FF3300")
dev.off()


#我们用ggplot2修饰下
data_sft = data.frame(x=sft$fitIndices[,1],
                      y=-sign(sft$fitIndices[,3])*sft$fitIndices[,2])

p1 <- ggplot(data_sft, aes(x,y))+
  geom_point(color="#1D75B5", size=5)+
  labs(x='Soft Threshold (power)',y='Scale Free Topology Model Fit,signed R^2',
       title = 'Scale independence')+
  theme_bw()+
  theme(axis.text = element_text(colour = 'black', size = 10),
        legend.position = "none")+
  geom_text(aes(label=x, color="#ED7C72"), size=3)+
  geom_hline(yintercept = 0.9, linewidth=0.5)

data_mean <- data.frame(x=sft$fitIndices[,1], y=sft$fitIndices[,5])
p2 <- ggplot(data_mean, aes(x,y))+
  geom_point(color="#1D75B5", size=5)+
  labs(x='Soft Threshold (power)',y='Mean Connectivity',
       title = 'Mean connectivity')+
  theme_bw()+
  theme(axis.text = element_text(colour = 'black', size = 10),
        legend.position = "none")+
  geom_text(aes(label=x, color="#ED7C72"), size=3)

library(patchwork)
p1 + p2

#ggsave("2.softpower_ggplot2.pdf", width = 8, height = 5)

sft$powerEstimate # 返回的是算法推荐的最佳的阈值，通常在五满足条件的power时选用


#### 6. 构建基因共表达网络 ####

##我们一步一步做

# power选择6以构建邻接矩阵
softPower = 6

adjacency = adjacency(bulk_data, power = softPower,type="signed")
# signed表示将不连接强负相关的基因表达谱
# 这里的type可以根据自己的数据进行选择，我们推荐选择signed或unsigned

# 邻接矩阵转为拓扑重叠
TOM = TOMsimilarity(adjacency, TOMType = "signed")
dissTOM = 1-TOM

# 对基因继续进行聚类
geneTree = hclust(as.dist(dissTOM), method = "average")

# 绘制基因聚类树图
#pdf("3.TomGeneCluteringTree.pdf", width = 8, height = 5)
sizeGrWindow(12,9)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04)
dev.off()

# 使用动态混合树切割算法来切割层次聚类树

# 最小模块大小
minModuleSize = 50;

# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, 
                            distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize, 
                            method = 'tree')

# 为每个module赋予颜色
dynamicColors = labels2colors(dynamicMods)

# 查看基因分成了多少module以及每个module的基因有多少
table(dynamicColors)

# Plot the dendrogram and colors underneath
sizeGrWindow(8,6)
#pdf("4.DynamicTreeCut.pdf", width = 8, height = 5)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")

dev.off()

# 合并模块（对相关性比较高的模块）
# 计算特征基因eigengenes
MEList = moduleEigengenes(bulk_data, colors = dynamicColors)
MEs = MEList$eigengenes

# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);

# 对模块eigengenes进行聚类
METree = hclust(as.dist(MEDiss), method = "average");

# Plot the result
sizeGrWindow(7, 6)
#pdf("5.ModulesEigengenesClustering.pdf", width = 10, height = 5)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")

MEDissThres = 0.3 # cut为0.3，对应相关>=0.7及以上的module合并,自行设置
# abline(h=MEDissThres, col = "red")
dev.off()

# Call an automatic merging function
merge = mergeCloseModules(bulk_data, dynamicColors, cutHeight = MEDissThres, verbose = 3)
mergedColors = merge$colors
mergedMEs = merge$newMEs
table(mergedColors)

# 作图--看一下合并前后的对比
#pdf("6.DynamicTreeCutAndMergedDynamicTreeCutContrast.pdf", width = 8, height = 5)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

# 接下来可视化合并的模块，并将每个样本也添加在上面，看看效果
# 单独看merge module
#pdf("7.module_merge.pdf", width = 6, height = 5)
plotDendroAndColors(geneTree, mergedColors,c("colors"),
                    dendroLabels = FALSE, hang = 0.01,
                    addGuide = TRUE, guideHang = 0.01)
dev.off()


####  7. 样本与Module相关性可视化 ####

unique(datTraits1$group_list)

colnames(pheno)

# 将数据变成numeric
DKD = as.numeric(pheno$DKD)
Normal = as.numeric(pheno$Normal)

# 生成一个空的matrix用于存放相关系数
geneCor=matrix(NA,nrow=2,ncol=ncol(bulk_data))

# 计算相关性
for(i in 1:ncol(geneCor)) {
  
  # 提取第i列
  expr=as.numeric(bulk_data[,i])
  
  DKD_r=bicor(expr,DKD,use="pairwise.complete.obs")
  Normal_r=bicor(expr, Normal,use="pairwise.complete.obs")
  
  geneCor[,i]=c(DKD_r, Normal_r)
  cat('Done for gene...',i,'\n')
}

range(geneCor)


# 赋予颜色
library(RColorBrewer)
for(i in 1:2){
  geneCor[i,] =numbers2colors(as.numeric(geneCor[i,]),
                              signed=TRUE,centered=TRUE,
                              colorRampPalette(c("#1D75B5","#d8d8d8","#ED7C72"))(50),
                              lim=c(-0.8,0.8))
}

rownames(geneCor) = c("DKD","Normal")

# 作图
#pdf("8.module_cor.pdf", width = 6, height = 5)
plotDendroAndColors(geneTree, cbind(mergedColors,
                                    geneCor[1,], 
                                    geneCor[2,]),
                    groupLabels=c("colors","DKD","Normal"),
                    dendroLabels = FALSE, hang = 0.01,
                    addGuide = TRUE, guideHang = 0.01)
dev.off()

#### 8. Module与表型的相关性 ####

# 这就是很常规的了，大多数文献常见的内容，看一下module与性状的相关性
# 使用merge后的每个module的eigengene与表型数据进行相关分析

# 提取merge后的特征基因
ME_merge <- merge$newMEs
moduleTraitCor <- cor(ME_merge, pheno, use = 'p')
range(moduleTraitCor)

# 相关p值、看显著性
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nrow(ME_merge))
range(moduleTraitPvalue)

# 添加每个module的基因数（如果你不需要要添加，跳过这一步即可）
ME_gene = rep(NA,length(names(ME_merge)))
for(i in 1:length(names(ME_merge))){
  modules <- names(ME_merge)[i]
  modules <- gsub("ME", "", modules)
  sites <- which(names(table(mergedColors))==modules)
  ME_gene[i] = paste0(names(ME_merge)[i],"(",table(mergedColors)[[sites]],")",sep="")
}


# 替换moduleTraitCor和moduleTraitPvalue行名
rownames(moduleTraitCor) <- ME_gene
rownames(moduleTraitPvalue) <- ME_gene

# 将显著性p值替换为*号
textMatrix <-  signif(moduleTraitPvalue,2)#保留小数点后两位
textM = ifelse(textMatrix>=0.05, "Not sig", 
               ifelse(textMatrix <0.05 & textMatrix >0.01,"*",
                      ifelse(textMatrix <=0.01 & textMatrix >0.001, "**",
                             ifelse(textMatrix <=0.001 & textMatrix >0.0001, "***", "****"))))


# 相关系数-这里我选择展示所有相关系数
range(moduleTraitCor)
textCor <- signif(moduleTraitCor,2)#保留小数点后2位
range(textCor)

textMatrix1 = paste(textCor,"\n",'(',textM ,')', sep = '')
textMatrix1= matrix(textMatrix1,
                    ncol=ncol( moduleTraitPvalue),
                    nrow=nrow(moduleTraitPvalue))

textMatrix1[textMatrix1=="\n()"] <- "" 


# 有了这个数据，我们就可以做热图，WGCNA自带的labeledHeatmap函数就可以完成
# 这个函数出来的图也就是我们平时大多数文章中见到的图,这里我们稍微修饰了一下
par(mar = c(5, 10, 5, 2))
labeledHeatmap( Matrix = moduleTraitCor,
                xLabels = c("DKD","Normal"),
                yLabels = names(ME_merge),
                ySymbols = rownames(moduleTraitPvalue),
                colorLabels = FALSE,
                colors = blueWhiteRed(50),
                textMatrix = textMatrix1,
                setStdMargins = FALSE,
                cex.text = 0.5,
                zlim = c(-1, 1),
                cex.lab.x = 0.8,
                cex.lab.y = 0.8,
                xLabelsAdj=1,
                xLabelsAngle=0,
                main = paste("Module-trait relationships"))
#pdf("9.module_trait_relationships.pdf", width = 10, height = 10)
dev.off()


# 有了矩阵，至于热图的修饰，你可以选择Complexheatmap或者ggplot2做各种个性化形式展现
# 例如使用ggplot2做热图、这里首先把宽数据转化为长数据(这里的演示只是提供一种优化的思路)
suppressMessages(library(tidyr))
datPlot = as.data.frame(moduleTraitCor)
datPlot$module <- rownames(datPlot)
datPlot <-gather(datPlot, sample, Cor, 1:2)

datPvalue = as.data.frame(textMatrix)
datPvalue$module <- rownames(datPvalue)
datPvalue <-gather(datPvalue, sample, p, 1:2)

datPlot$p <- datPvalue$p
datPlot$log10p = -log10(as.numeric(datPlot$p) + 1e-5)
datPlot$star = ifelse(datPlot$p >=0.05, "", 
                      ifelse(datPlot$p <0.05 & datPlot$p >0.01,"*",
                             ifelse(datPlot$p <=0.01 & datPlot$p >0.001, "**",
                                    ifelse(datPlot$p <=0.001 & datPlot$p >0.0001, "***", "****"))))


datPlot$module = factor(datPlot$module,levels=rev(unique(datPlot$module)))
datPlot$sample = factor(datPlot$sample, levels = c("DKD","Normal"))

anno_module  = unique(datPlot$module)
anno_color = c()
for (i in 1:length(anno_module)) {
  
  a = anno_module[i]
  a <- gsub("ME|(\\d+)", "", a)
  a <- gsub("\\()", "", a)
  
  anno_color <- append(anno_color, a)
}

anno_data <- data.frame(Module = anno_module,
                        color = anno_color,
                        x = rep(1, length(anno_module)))

colnames(anno_data) <- c("Module","color","x")
p1 = ggplot(anno_data, aes(x=1,y=Module, fill=Module))+
  geom_tile()+ 
  theme_minimal()+xlab(NULL) + ylab(NULL) +
  theme(axis.text.x = element_blank(),
        legend.position = 'none',
        axis.ticks = element_blank(),
        panel.grid = element_blank())+
  scale_fill_manual(values = rev(anno_data$color))



p2 = ggplot(datPlot, aes(sample, module,fill=log10p)) + 
  geom_tile(aes(fill = log10p),colour = "grey50") + 
  scale_fill_gradient2(low = "#0D8CFF",high = "#FF3300") +
  geom_text(aes(label = star), color = "black", size = 5,vjust=0.75) +
  xlab("") +
  ylab("") +
  labs(fill="log10(FDR)") +
  theme_bw() +
  theme(legend.title = element_text(size = 8),
        legend.text = element_text(size = 8),
        axis.text.x = element_text(size = 8,angle=0,hjust=1),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.direction = "vertical",
        panel.grid = element_blank())


suppressMessages(library(patchwork))
plot =p1+p2+plot_layout(ncol = 2, widths  = c(0.2,4),guides = 'collect') # 左右组合，guides = 'collect'表示三个图的legend组合在一起，widths表示从左到右图的宽度
plot & theme(plot.margin = margin(0.5,0.5,0.5,0.5)) # 排列紧密
ggsave("9.module_trait_relationships_ggplot2.pdf", width = 10, height = 10)


#### 9. 提取模块基因 ####

# 这里我们要提取turquoise模块的基因
# 要提取的是merge之后的

# 提取module基因
moduleColors = mergedColors
moduleColors

colorOrder = c("grey", unique(c(standardColors(50), unique(mergedColors))))

moduleLabels = match(moduleColors, colorOrder)-1

# Modules
module_dataframe <- data.frame(gene_id=colnames(bulk_data), 
                               module_name=paste0('module_', moduleLabels), 
                               module_color=moduleColors)

# 输出所有的基因
write.csv(module_dataframe, file = "module_dataframe_gene_all.csv")

# 只输出MEmagenta以及MEroyalblue这两个模块的基因
table(module_dataframe$module_color)
library(dplyr)
module_dataframe_filtered <- subset(module_dataframe, module_color %in% c("turquoise"))
write.csv(module_dataframe, file = "module_dataframe_gene_filtered.csv")

# 2586个特征基因


#### 10. Define variable weight containing all column of datTraits ####

modNames = substring(names(mergedMEs), 3)

#首先获取modulememership
geneModuleMembership = as.data.frame(cor(bulk_data, mergedMEs, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nrow(bulk_data)))

names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="")

#names of those trait
traitNames=c("DKD","Normal")
traitNames

#这里我们读取预先准备好的datTraits
#1代表是，0代表否
datTraits2 <- read.csv("E:/本地生信项目/DKD/RNAseq/WGCNA/datTraits.csv",row.names = 1)

#然后获得GS
geneTraitSignificance = as.data.frame(cor(bulk_data, datTraits2, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nrow(bulk_data)))

names(geneTraitSignificance) = paste("GS.", traitNames, sep="")
names(GSPvalue) = paste("p.GS.", traitNames, sep="")


####绘制 MM VS GS图

#这里我们绘制所有module的 MM VS GS图
traitNames
modNames

# 默认可视化
for (trait in traitNames){
  traitColumn=match(trait,traitNames)
  
  for (module in modNames){
    column = match(module, modNames)
    moduleGenes = moduleColors==module
    
    if (nrow(geneModuleMembership[moduleGenes,]) > 1){
      
      pdf(file=paste("./10_", trait, "_", module,"_Module membership vs gene significance.pdf"),width = 7,height = 7)
      verboseScatterplot(x = abs(geneModuleMembership[moduleGenes, column]),
                         y = abs(geneTraitSignificance[moduleGenes, traitColumn]),
                         xlab = paste("Module Membership in", module, "module"),
                         ylab = paste("Gene significance for ",trait),
                         main = paste("Module membership vs. gene significance\n"),
                         cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module,
                         abline = T,abline.color = "#FF3300")
      
      dev.off()
    }
  }
}

names(bulk_data)
probes = names(bulk_data)


##导出 MM VS GS 的所有模块基因
geneInfo0 = data.frame(probes= probes,
                       moduleColor = moduleColors)

for (Tra in 1:ncol(geneTraitSignificance))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneTraitSignificance[,Tra],
                         GSPvalue[, Tra])
  names(geneInfo0) = c(oldNames,names(geneTraitSignificance)[Tra],
                       names(GSPvalue)[Tra])
}

for (mod in 1:ncol(geneModuleMembership))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneModuleMembership[,mod],
                         MMPvalue[, mod])
  names(geneInfo0) = c(oldNames,names(geneModuleMembership)[mod],
                       names(MMPvalue)[mod])
}

geneOrder =order(geneInfo0$moduleColor)
geneInfo = geneInfo0[geneOrder, ]

write.csv(geneInfo, file = paste0("./11_GS_and_MM.csv"),row.names=F)


##输出达到过滤阈值的基因

#我们只要turquoise中GS >0.25 和 MM >0.7的基因
module = "turquoise"
column = match(module, modNames)
column

# 查看module基因
moduleGenes = moduleColors == module
table(moduleGenes)

# 提取基因名
turquoise_module <- as.data.frame(dimnames(data.frame(bulk_data))[[2]][moduleGenes]) 
names(turquoise_module) = "genename"

MM <- abs(geneModuleMembership[moduleGenes,column])
GS <- abs(geneTraitSignificance[moduleGenes, 1])
c <- as.data.frame(cbind(MM,GS))
c
rownames(c) = turquoise_module$genename
head(c)

#我们只要turquoise中GS >0.25 和 MM >0.7的基因
turquoise_hub <- abs(c$MM)>0.7 & abs(c$GS)>0.25
table(turquoise_hub)
# 370个基因

write.csv(turquoise_hub, "hubgene_MMGS_turquoise.csv")


####重新使用ggplot2绘制图

turquoise_hub <- read.csv("hubgene_MMGS_turquoise.csv")

turquoise_hub <- as.data.frame(turquoise_hub)

head(turquoise_hub)
table(turquoise_hub$x)
colnames(turquoise_hub) <- c("gene","group")

#将匹配信息添加到散点图矩阵最后一列
c <- cbind(c,turquoise_hub$group)
head(c)
colnames(c) <- c("MM","GS","group")

#提取符合阈值的基因
table(c$group)
turquoisegene <- c %>% filter(group == TRUE)

#提取基因
write.csv(turquoisegene, file = paste0("./turquoise_gene_filtered.csv"))

#提取相应基因的表达量
expr_filtered <- bulk_data[,colnames(bulk_data) %in% rownames(turquoisegene)]
write.csv(expr_filtered, file = paste0("./turquoise_gene_filtered_expr.csv"))

#利用ggplot2绘图
library(ggplot2)
pdf("MM vs. GS_turquoise_DKD_Normal.pdf",width = 7,height = 7)
ggplot(data=c, aes(x=MM, y=GS, color=group))+
  geom_point(size=1.5)+
  scale_colour_manual(values=c("grey60", "turquoise"))+ 
  theme_bw()+  
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+  
  labs(x="Module Membership in turquoise module", y="Gene significance for TL",title = "Module membership vs. gene significance ")+
  theme(axis.title.x =element_text(size=14), axis.title.y=element_text(size=14),axis.text = element_text(size = 12),axis.text.x = element_text(colour = "black"),axis.text.y = element_text(colour = "black"),plot.title = element_text(hjust = 0.5,size = 16,face = "bold"),plot.margin = unit(rep(2,4),'lines')) +
  theme(legend.position = 'none')+geom_hline(aes(yintercept=0.25),colour="#FF3300",lwd=1,linetype=5)+
  geom_vline(aes(xintercept=0.7),colour="#FF3300",lwd=1,linetype=5)
dev.off()


