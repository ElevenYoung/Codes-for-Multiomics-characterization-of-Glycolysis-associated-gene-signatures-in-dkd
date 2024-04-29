#### 提取糖酵解基因 ####

rm(list = ls())
setwd("E:/本地生信项目/糖酵解基因集")
options(stringsAsFactors = FALSE)
suppressMessages(library(readxl))
suppressMessages(library(dplyr))
suppressMessages(library(data.table))
suppressMessages(library(tidyverse))
suppressMessages(library(stringr))


#### 1. 读取基因集 ####

##genecards
genecards <- read.csv("E:/本地生信项目/糖酵解基因集/rawdata/GeneCards/GeneCards-SearchResults.csv")

##msigdb
msigdb <- clusterProfiler::read.gmt("E:/本地生信项目/糖酵解基因集/rawdata/msigdb/genesets.v2023.2.Hs.gmt")


#### 2. 整理基因集 ####

# 我们要分别把基因集整理成第一列为gene symbol，第二列为source的格式

##genecards
genecards <- data.frame(gene_symbol = genecards$Gene.Symbol, source = "Genecards")

#去除重复基因
genecards <- genecards[!duplicated(genecards$gene_symbol),]

##msigdb
table(msigdb$term)

#更改名称
msigdb$term <- gsub("HALLMARK_GLYCOLYSIS", "Msigdb_HALLMARK_GLYCOLYSIS", msigdb$term)
msigdb$term <- gsub("REACTOME_GLYCOLYSIS", "Msigdb_REACTOME_GLYCOLYSIS", msigdb$term)

table(msigdb$term)

#去除重复基因
msigdb <- msigdb[!duplicated(msigdb$gene),]

#整理成dataframe
msigdb <- data.frame(gene_symbol = msigdb$gene, source = msigdb$term)


##合并基因集
gly_geneset <- rbind(genecards,msigdb)

##去除重复行
gly_geneset <- gly_geneset[!duplicated(gly_geneset$gene_symbol), ]

#我们发现，REACTOME的基因全被删了
#秀

##保存数据
write.csv(gly_geneset,"gly_geneset.csv")