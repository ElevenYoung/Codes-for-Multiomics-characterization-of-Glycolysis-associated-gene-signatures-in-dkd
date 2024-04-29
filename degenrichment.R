#### deg的GO/KEGG/REACTOME/GSEA/GSVA富集分析 ####

#-----------------------------------------------------------------------------------------------------------------




#### 1. deg的GO富集分析 ####

# 清空当前环境
rm(list = ls())
getwd()
options(stringsAsFactors = FALSE)
suppressMessages(library(ggpubr)) # 绘图使用
suppressMessages(library(ggplot2)) # 绘图使用
suppressMessages(library(clusterProfiler)) # 富集分析使用
suppressMessages(library(org.Hs.eg.db)) # 转换基因ID使用
suppressMessages(library(dplyr)) # 数据处理使用


#### 1.2 数据读取 ####

# 这里我们用limma差异分析后的数据
marker <- read.csv("E:/本地生信项目/DKD/RNAseq/Train_Cohort/enrichment/deg_Enrichment/deg_logFC0.5_filtered.csv", row.names = 1)

table(marker$deg)

up <- filter(marker,deg == "Up")
down <- filter(marker,deg == "Down")

#### 1.3 数据预处理 ####

# Gene Symbol转化为ENTREZID
gene.df1 <- bitr(rownames(up), # 基因列表(GENE SYMBOL)
                fromType = "SYMBOL", # fromType是指你的数据ID类型是属于哪一类的
                toType = c("ENTREZID", "SYMBOL"), # toType是指你要转换成哪种ID类型，可以写多种，也可以只写一种
                OrgDb = org.Hs.eg.db # Orgdb是指对应的注释包是哪个
)

gene.df2 <- bitr(rownames(down), # 基因列表(GENE SYMBOL)
                fromType = "SYMBOL", # fromType是指你的数据ID类型是属于哪一类的
                toType = c("ENTREZID", "SYMBOL"), # toType是指你要转换成哪种ID类型，可以写多种，也可以只写一种
                OrgDb = org.Hs.eg.db # Orgdb是指对应的注释包是哪个
)

# 命名第一列为gene_name
colnames(gene.df1)[1] <- "gene_name"
colnames(gene.df2)[1] <- "gene_name"



#### 1.3 GO富集分析 ####

# 进行GO富集分析
GO_all1 <- enrichGO(gene = gene.df1$ENTREZID,  # 基因列表(ENTREZID)
                   keyType = "ENTREZID",  # 指定的基因ID类型，默认为ENTREZID
                   OrgDb = org.Hs.eg.db,  # 物种对应的org包
                   ont = "ALL",   # CC细胞组件，MF分子功能，BF生物学过程，ALL以上三个
                   pvalueCutoff = 0.05,  # p值阈值
                   pAdjustMethod = "fdr",  # 多重假设检验校正方式
                   qvalueCutoff = 0.05,  # q值阈值
                   readable = TRUE # 是否基因ID转换为基因名
)

# 进行GO富集分析
GO_all2 <- enrichGO(gene = gene.df2$ENTREZID,  # 基因列表(ENTREZID)
                   keyType = "ENTREZID",  # 指定的基因ID类型，默认为ENTREZID
                   OrgDb = org.Hs.eg.db,  # 物种对应的org包
                   ont = "ALL",   # CC细胞组件，MF分子功能，BF生物学过程，ALL以上三个
                   pvalueCutoff = 0.05,  # p值阈值
                   pAdjustMethod = "fdr",  # 多重假设检验校正方式
                   qvalueCutoff = 0.05,  # q值阈值
                   readable = TRUE # 是否基因ID转换为基因名
)

# 将结果转换为dataframe
GO_result1 <- data.frame(GO_all1@result)
GO_result2 <- data.frame(GO_all2@result)

#### 1.4 保存结果 ####

# 保存为csv
write.csv(GO_result1, file = 'E:/本地生信项目/DKD/RNAseq/Train_Cohort/enrichment/deg_Enrichment/deg_Enrichment/GO_result_up.csv')
write.csv(GO_result2, file = 'E:/本地生信项目/DKD/RNAseq/Train_Cohort/enrichment/deg_Enrichment/GO_result_down.csv')


#-----------------------------------------------------------------------------------------------------------------



#### 2. deg的KEGG富集分析 ####


# 进行KEGG富集分析
kegg1 <- enrichKEGG(gene = gene.df1$ENTREZID,  # 基因列表(ENTREZID)
                   keyType = "kegg",  # 指定的基因ID类型，默认为kegg
                   organism = 'hsa',  # 对应的物种
                   pAdjustMethod = "fdr",  # 多重假设检验校正方式
)

kegg2 <- enrichKEGG(gene = gene.df2$ENTREZID,  # 基因列表(ENTREZID)
                   keyType = "kegg",  # 指定的基因ID类型，默认为kegg
                   organism = 'hsa',  # 对应的物种
                   pAdjustMethod = "fdr",  # 多重假设检验校正方式
)


# 将结果转换为dataframe
kegg_result1 <- data.frame(kegg1@result)
kegg_result2 <- data.frame(kegg2@result)

#### 2.4 保存结果 ####

# 保存为csv
write.csv(kegg_result1, file = 'E:/本地生信项目/DKD/RNAseq/Train_Cohort/enrichment/deg_Enrichment/kegg_up.csv')
write.csv(kegg_result2, file = 'E:/本地生信项目/DKD/RNAseq/Train_Cohort/enrichment/deg_Enrichment/kegg_down.csv')


#-----------------------------------------------------------------------------------------------------------------


#### 3. deg的REACTOME富集分析 ####

suppressMessages(library(ReactomePA))

t1 <- as.vector(gene.df1$ENTREZID)
t2 <- as.vector(gene.df2$ENTREZID)


# 进行REACTOME富集分析
reactome1 <- enrichPathway(gene = t1,  # 基因列表(ENTREZID)
                           organism = 'human',  # 对应的物种
                           pAdjustMethod = "fdr",  # 多重假设检验校正方式
                           readable = T
)

reactome2 <- enrichPathway(gene = t2,  # 基因列表(ENTREZID)
                           organism = 'human',  # 对应的物种
                           pAdjustMethod = "fdr",  # 多重假设检验校正方式
                           readable = T
)


# 将结果转换为dataframe
reactome_result1 <- data.frame(reactome1@result)
reactome_result2 <- data.frame(reactome2@result)


#### 保存结果

# 保存为csv
write.csv(reactome_result1, file = 'E:/本地生信项目/DKD/RNAseq/Train_Cohort/enrichment/deg_Enrichment/reactome_up.csv')
write.csv(reactome_result2, file = 'E:/本地生信项目/DKD/RNAseq/Train_Cohort/enrichment/deg_Enrichment/reactome_down.csv')


#-----------------------------------------------------------------------------------------------------------------


#### 4. deg的GSEA富集分析 ####


# 清空当前环境
rm(list = ls())
getwd()
options(stringsAsFactors = FALSE)
suppressMessages(library(ggpubr)) # 绘图使用
suppressMessages(library(ggplot2)) # 绘图使用
suppressMessages(library(clusterProfiler)) # 富集分析使用
suppressMessages(library(org.Hs.eg.db)) # 转换基因ID使用
suppressMessages(library(dplyr)) # 数据处理使用
suppressMessages(library(msigdbr)) # 数据处理使用

#### 4.2 数据读取 ####

# 这里我们用limma差异分析后的数据
marker <- read.csv("E:/本地生信项目/DKD/RNAseq/Train_Cohort/enrichment/deg_Enrichment/deg_logFC0.5_filtered.csv", row.names = 1)

table(marker$deg)

up <- filter(marker,deg == "Up")
up$genesymbol <- rownames(up)

down <- filter(marker,deg == "Down")
down$genesymbol <- rownames(down)

marker$genesymbol <- rownames(marker)


#### 4.3 数据预处理 ####

# Gene Symbol转化为ENTREZID
gene.df1 <- bitr(rownames(up), # 基因列表(GENE SYMBOL)
                 fromType = "SYMBOL", # fromType是指你的数据ID类型是属于哪一类的
                 toType = c("ENTREZID", "SYMBOL"), # toType是指你要转换成哪种ID类型，可以写多种，也可以只写一种
                 OrgDb = org.Hs.eg.db # Orgdb是指对应的注释包是哪个
)

gene.df2 <- bitr(rownames(down), # 基因列表(GENE SYMBOL)
                 fromType = "SYMBOL", # fromType是指你的数据ID类型是属于哪一类的
                 toType = c("ENTREZID", "SYMBOL"), # toType是指你要转换成哪种ID类型，可以写多种，也可以只写一种
                 OrgDb = org.Hs.eg.db # Orgdb是指对应的注释包是哪个
)

gene.df3 <- bitr(rownames(marker), # 基因列表(GENE SYMBOL)
                 fromType = "SYMBOL", # fromType是指你的数据ID类型是属于哪一类的
                 toType = c("ENTREZID", "SYMBOL"), # toType是指你要转换成哪种ID类型，可以写多种，也可以只写一种
                 OrgDb = org.Hs.eg.db # Orgdb是指对应的注释包是哪个
)

# 命名第一列为gene_name
colnames(gene.df1)[1] <- "SYMBOL"
colnames(gene.df2)[1] <- "SYMBOL"
colnames(gene.df3)[1] <- "SYMBOL"

# up
gene_entrezid <- merge(gene.df1,up,by.x = "SYMBOL", by.y = "genesymbol")
genelist <- gene_entrezid$logFC
names(genelist) <- gene_entrezid$ENTREZID
genelist1 <- sort(genelist,decreasing = T)

# down
gene_entrezid <- merge(gene.df2,down,by.x = "SYMBOL", by.y = "genesymbol")
genelist <- gene_entrezid$logFC
names(genelist) <- gene_entrezid$ENTREZID
genelist2 <- sort(genelist,decreasing = T)

# all
gene_entrezid <- merge(gene.df3,marker,by.x = "SYMBOL", by.y = "genesymbol")
genelist <- gene_entrezid$logFC
names(genelist) <- gene_entrezid$ENTREZID
genelist3 <- sort(genelist,decreasing = T)

# 做GSEA分析对数据格式有要求，需要是一个有序的数值型向量，其名字是基因的ID


#### 4.4 GSEA分析 ####

library(msigdbr)
m_t2g <- msigdbr(species = "Homo sapiens", category = "H") %>% 
  dplyr::select(gs_name, entrez_gene)
head(m_t2g)

table(m_t2g$gs_name)

# up
gsea_res1 <- GSEA(genelist1, 
                 TERM2GENE = m_t2g,
                 pvalueCutoff = 0.05,
                 pAdjustMethod = "BH"
)


# down
gsea_res2 <- GSEA(genelist2, 
                  TERM2GENE = m_t2g,
                  pvalueCutoff = 0.05,
                  pAdjustMethod = "BH"
)

# all
gsea_res3 <- GSEA(genelist3, 
                  TERM2GENE = m_t2g,
                  pvalueCutoff = 0.05,
                  pAdjustMethod = "BH"
)


## 保存结果
df <- as.data.frame(gsea_res3)
write.csv(df,"E:/本地生信项目/DKD/RNAseq/Train_Cohort/enrichment/deg_Enrichment/gsea_res.csv")


#### 4.5 gsea结果可视化 ####

library(GseaVis)

##上调的通路
# bacth plot
terms <- c('HALLMARK_XENOBIOTIC_METABOLISM',
           'HALLMARK_KRAS_SIGNALING_DN',
           'HALLMARK_OXIDATIVE_PHOSPHORYLATION',
           'HALLMARK_FATTY_ACID_METABOLISM',
           'HALLMARK_BILE_ACID_METABOLISM',
           'HALLMARK_PEROXISOME')

# plot
lapply(terms, function(x){
  gseaNb(object = gsea_res3,
         geneSetID = x,
         addPval = T,
         pvalX = 0.75,pvalY = 0.75,
         pCol = 'black',
         pHjust = 0)
}) -> gseaList

# combine
cowplot::plot_grid(plotlist = gseaList,ncol = 3,align = 'hv')


##下调的通路
# bacth plot
terms <- c('HALLMARK_ALLOGRAFT_REJECTION',
           'HALLMARK_INTERFERON_GAMMA_RESPONSE',
           'HALLMARK_INTERFERON_ALPHA_RESPONSE',
           'HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION',
           'HALLMARK_IL6_JAK_STAT3_SIGNALING',
           'HALLMARK_INFLAMMATORY_RESPONSE')

# plot
lapply(terms, function(x){
  gseaNb(object = gsea_res3,
         geneSetID = x,
         addPval = T,
         pvalX = 0.75,pvalY = 0.75,
         pCol = 'black',
         pHjust = 0)
}) -> gseaList

# combine
cowplot::plot_grid(plotlist = gseaList,ncol = 3,align = 'hv')
