#### ML 相关信息 可视化 ###

# 上接script.R的内容

# 我们选择SVM+Enet[alpha=0.6]选择的特征基因

#### 1. 模型可视化 ####

modeldf <- model$`SVM+Enet[alpha=0.6]`
modeldf


# 第二次我们就只运行SVM+Enet[alpha=0.6]
methods <- "SVM+Enet[alpha=0.6]"
methods

## Train the model --------------------------------------------------------
## 设置所要预测的变量名（仅支持[0,1]二元变量格式）
classVar = "outcome"
## 设置模型最少纳入的变量数
min.selected.var = 3

Train_set <- Train_expr
Test_set <- Test_expr

## Pre-training --------------------------------------------------------
## 提取基因
Variable = colnames(Train_set)

## 检视所有方法，分析各方法是否需要进行变量预筛选(pre-training)
preTrain.method =  strsplit(methods, "\\+")
preTrain.method

## 删除各方法用于构建分类模型的算法，保留用于变量筛选的算法
preTrain.method = lapply(preTrain.method, function(x) rev(x)[-1])
preTrain.method

## 汇总所有变量筛选算法，去除重复计算
preTrain.method = unique(unlist(preTrain.method))
preTrain.method

## 用于保存各算法筛选的变量
preTrain.var <- list()

for (method in preTrain.method){
  preTrain.var[[method]] = RunML(method = method, # 变量筛选所需要的机器学习方法
                                 Train_set = Train_set, # 训练集有潜在预测价值的变量
                                 Train_label = Train_class, # 训练集分类标签
                                 mode = "Variable",       # 运行模式，Variable(筛选变量)和Model(获取模型)
                                 classVar = classVar) # 用于训练的分类变量，必须出现在Train_class中
}

# 记录未经筛选的变量集（以便后续代码撰写），可视为使用simple方法（无筛选功能）的变量筛选结果
preTrain.var[["simple"]] <- colnames(Train_set)
preTrain.var


## Model training --------------------------------------------------------
## 创建空list用于保存各模型的所有信息
model <- list()

## 构建基于SVM+Enet[alpha=0.6]的ML模型
  
  method <- methods
  # 循环每一种方法组合
  
  # 输出当前方法
  cat(match(method, methods), ":", method, "\n")
  
  # 本轮算法名称
  method_name = method
  print(method_name)
  
  # 各步骤算法名称
  method <- strsplit(method, "\\+")[[1]]
  print(method)
  
  # 如果本方法没有预筛选变量，则认为本方法使用simple方法进行了变量筛选
  if (length(method) == 1) method <- c("simple", method)   
  
  Variable = preTrain.var$SVM # 根据方法名称的第一个值，调用先前变量筛选的结果
  Train_set = Train_set_bk[, Variable]   # 对训练集取子集，因为有一个算法原作者写的有点问题，无法正常传参
  Train_label = Train_class            # 所以此处需要修改变量名称，以免函数错误调用对象
  method = method[2]        # 根据方法名称第二个值，调用构建的函数分类模型
  Train_set = Train_set     # 训练集有潜在预测价值的变量
  Train_label = Train_label # 训练集分类标签
  mode = "Model"            # 运行模式，Variable(筛选变量)和Model(获取模型)
  classVar = classVar       # 用于训练的分类变量，必须出现在Train_class中


  # 使用训练数据（Train_set）和训练集标签（Train_label）来进行交叉验证，并选择最优的lambda值
  cv.fit = cv.glmnet(x = Train_set,
                     y = Train_label[[classVar]],
                     family = "binomial", alpha = 0.6, nfolds = 10)
  
  plot(cv.fit,xvar = "lambda",label = T)
  
  # 使用选定的最优lambda值来拟合最终的Logistic模型
  lasso = glmnet(x = Train_set,
               y = Train_label[[classVar]],
               family = "binomial", alpha = 0.6)
  
  plot(lasso,xvar = "lambda",label = T,lwd=2)

## 下面使用ggplot2对cv.fit进行美化
  library(ggplot2)
  library(glmnet)
  library(survival)
  library(tidyr)
  
  
  # 提取cv.glmnet数据
  cv.fit
  
  cv.df <- data.frame(lambda =cv.fit$lambda,#交叉验证中的lambda
                      mse =cv.fit$cvm,#mse
                      sd =cv.fit$cvsd)#sd
  
  # 绘制10折交叉验证过程图
  p1<-ggplot(cv.df, aes(log(lambda), mse)) +
    geom_errorbar(aes(ymin = mse - sd, ymax = mse + sd), width = 0.1, size = 0.7,color = "grey") +#添加误差棒
    geom_point(color = "red", size = 2) +#点图
    scale_x_continuous(name = "Log lambda") +
    scale_y_continuous(name = "Binomial Deviance") +
    ggtitle("10-fold Cross-validation using ElasticNet Regression")+
    geom_vline(xintercept = log(cv.fit$lambda.min), linetype = "dashed", color = "#ED7C72",size=1) +
    geom_vline(xintercept = log(cv.fit$lambda.1se), linetype = "dashed", color = "#1D75B5",size=1) +
    annotate(geom = "text",label=("lambda.min: 16 non-zero variables"),x=-7.5,y=12,color = "#ED7C72",size=5)+
    annotate(geom = "text",label=("lambda.lse: 1 non-zero variables"),x=-7.5,y=11,color = "#1D75B5",size=5)+
    theme_bw()+
    theme(axis.text = element_text(size = 10),
          axis.title = element_text(size = 12),
          )
  p1
  
  
  ## 查看每个基因的coef
  coef.min <- coef(cv.fit,s="lambda.min")
  coef.min
  coef.min <- as.matrix(coef.min)
  coef.min <- as.data.frame(coef.min[-1,])
  colnames(coef.min) <- c("coef")
  coef.min <- coef.min %>%
    tibble::rownames_to_column(var = "gene_name")
  coef.min <- coef.min[coef.min$coef != 0, ]#保留有coef的基因
  # 保存数据
  #write.csv(coef.min,"./coef.csv")

  ##棒棒糖图可视化
  dat <- mutate(coef.min,group = ifelse(coef.min$coef<0,"n","p")) %>%  
    arrange(coef)
  dat$gene = factor(dat$gene_name,levels = dat$gene_name)
  head(dat)
  
  ggplot(dat, aes(x=gene, y=coef)) +
    geom_segment( aes(x=gene, xend=gene, y=0, yend=coef), color="grey",linewidth = 2) +
    geom_point( aes(color = group), size=4) +
    scale_color_manual(values = c("#1D75B5", "#ED7C72"))+ 
    theme_bw() +
    theme(legend.position = "none") +
    coord_flip()
    
  
  ## 训练集和验证集的汇总ROC曲线
  library(pROC)
  
  ## 对"SVM+Enet[alpha=0.6]"模型计算AUC值
  table(Test_class$outcome)
  table(Train_class$outcome)
  unique(Test_class$outcome) %in% unique(Train_class$outcome)
  
  fit = modeldf    # 分类预测模型
  Test_set = Test_set      # 测试集预测变量，应当包含训练集中所有的变量，否则会报错
  Test_label = Test_class   # 训练集分类数据，应当包含训练集中所有的变量，否则会报错
  Train_set = Train_set    # 若需要同时评估训练集，则给出训练集表达谱，否则置NULL
  Train_label = Train_class # 若需要同时评估训练集，则给出训练集分类数据，否则置NULL
  Train_name = "MetaGSE"       # 若需要同时评估训练集，可给出训练集的标签，否则按“Training”处理
  cohortVar = "Cohort"      # 重要：用于指定队列的变量，该列必须存在且指定[默认为“Cohort”]，否则会报错
  classVar = classVar       # 用于评估的二元分类变量，必须出现在Test_class中
 

    if(!is.element(cohortVar, colnames(Test_label))) {
      stop(paste0("There is no [", cohortVar, "] indicator, please fill in one more column!"))
    } 
    
    if((!is.null(Train_set)) & (!is.null(Train_label))) {
      new_data <- rbind.data.frame(Train_set[, fit$subFeature],
                                   Test_set[, fit$subFeature])
      
      if(!is.null(Train_name)) {
        Train_label$Cohort <- Train_name
      } else {
        Train_label$Cohort <- "Training"
      }
      colnames(Train_label)[ncol(Train_label)] <- cohortVar
      Test_label <- rbind.data.frame(Train_label[,c(cohortVar, classVar)],
                                     Test_label[,c(cohortVar, classVar)])
      Test_label[,1] <- factor(Test_label[,1], 
                               levels = c(unique(Train_label[,cohortVar]), setdiff(unique(Test_label[,cohortVar]),unique(Train_label[,cohortVar]))))
    } else {
      new_data <- Test_set[, fit$subFeature]
    }
    
    RS <- suppressWarnings(CalPredictScore(fit = fit, new_data = new_data))
    
    Predict.out <- Test_label
    Predict.out$RS <- as.vector(RS)
    Predict.out <- split(x = Predict.out, f = Predict.out[,cohortVar])
    
    #下面我们一个一个计算每个集的auc
    MetaGSE <- Predict.out[[1]]
    GSE96804 <- Predict.out[[2]]
    GSE30528 <- Predict.out[[3]]
    GSE30529 <- Predict.out[[4]]
    
    #MetaGSE
    roc1 <- roc(MetaGSE[[classVar]], MetaGSE$RS)
    ci1=ci.auc(roc1, method="bootstrap")
    plot(roc1, print.auc = TRUE, auc.polygon = TRUE, grid = TRUE, col = "black", auc.polygon.col="#ED7C72")
    text(0.39, 0.39, paste0("95% CI: ",sprintf("%.03f",ci1[1]),"-",sprintf("%.03f",ci1[3])), col="black")
    
    roc2 <- roc(GSE96804[[classVar]], GSE96804$RS)
    ci2=ci.auc(roc2, method="bootstrap")
    plot(roc2, print.auc = TRUE, auc.polygon = TRUE, grid = TRUE, col = "black", auc.polygon.col="#ED7C72")
    text(0.39, 0.39, paste0("95% CI: ",sprintf("%.03f",ci2[1]),"-",sprintf("%.03f",ci2[3])), col="black")
    
    roc3 <- roc(GSE30528[[classVar]], GSE30528$RS)
    ci3=ci.auc(roc3, method="bootstrap")
    plot(roc3, print.auc = TRUE, auc.polygon = TRUE, grid = TRUE, col = "black", auc.polygon.col="#ED7C72")
    text(0.39, 0.39, paste0("95% CI: ",sprintf("%.03f",ci3[1]),"-",sprintf("%.03f",ci3[3])), col="black")
    
    roc4 <- roc(GSE30529[[classVar]], GSE30529$RS)
    ci4=ci.auc(roc4, method="bootstrap")
    plot(roc1, print.auc = TRUE, auc.polygon = TRUE, grid = TRUE, col = "black", auc.polygon.col="#ED7C72")
    text(0.39, 0.39, paste0("95% CI: ",sprintf("%.03f",ci4[1]),"-",sprintf("%.03f",ci4[3])), col="black")
    
    
    ##标志基因的nomogram可视化（这里为MetaGSE）
    #读取表达文件
    rt = read.csv("E:/本地生信项目/DKD/RNAseq/Train_Cohort/exprmatrix_remove_batch_remove.csv", row.names = 1)
    rt = as.matrix(rt)

    #读取分组信息
    sample = read.csv("E:/本地生信项目/DKD/RNAseq/Train_Cohort/groupinfo.csv", row.names = 1)
    sample$ID <- rownames(sample)
    colnames(sample) = c("Type","Tissue","Dataset","ID")
    data = rt[,sample$ID]
    
    data = t(data)
    #行为样本，列为基因
    
    #提取特征基因的表达矩阵
    fea_df <- fea_df[fea_df$algorithm == "SVM+Enet[alpha=0.6]",]

    data <- data[,colnames(data) %in% fea_df$features]

    aSAH = cbind(sample,data)
    
    dd <- datadist(aSAH)
    options(datadist = "dd")
    fit <- lrm(formula = Type ~ CD163+CYBB+ELF3+FCN1+GPR65+LCN2+LTF+PROM1+S100A4+SOX4+TGFBI+TNFAIP8, data = aSAH)
    print(fit)
    coef = as.data.frame(fit$coefficients)[-1,,drop = F]
    coefout = cbind(ID = rownames(coef),coef)
    write.table(coefout,file = "E:/本地生信项目/DKD/RNAseq/ML_combination/nomograme_coefficients.txt",sep = "\t",quote = F,row.names = F)
    
    #绘图
    pdf(file = "E:/本地生信项目/DKD/RNAseq/ML_combination/nomogram.pdf", width = 9, height = 7.5)
    plot(nomogram(fit,fun.at = seq(0.05,0.95,0.05)),funlabel = "nomogram model")
    dev.off()
    
    plot(regplot(fit,plots = c("density","boxes"), observation = T, title = "Prediction Nomogram", clickable = T, points = TRUE, droplines = TRUE))
    
    nomoscore = predict(fit, data = t(aSAH))
    aSAH$nomoscore = nomoscore
    write.table(aSAH,file = "E:/本地生信项目/DKD/RNAseq/ML_combination/nomoscore.txt", sep = "\t", quote = F, row.names = F)

    
    ## 十二个特征基因在训练集与验证集的ROC曲线汇总    
  
    #下面我们一个一个计算每个基因在每个集的auc
    #查看特征基因
    fea_df
    
    
    ## MetaGSE
    #读取表达数据
    rt = read.csv("E:/本地生信项目/DKD/RNAseq/Train_Cohort/exprmatrix_remove_batch_remove.csv", row.names = 1)
    rt = as.matrix(rt)
    #读取分组信息
    sample = read.csv("E:/本地生信项目/DKD/RNAseq/Train_Cohort/groupinfo.csv", row.names = 1)
    sample$ID <- rownames(sample)
    colnames(sample) = c("Type","Tissue","Dataset","ID")
    data = rt[,sample$ID]
    #数据转置
    data = t(data)#行为样本，列为基因
    #提取特征基因的表达矩阵
    data <- data[,colnames(data) %in% fea_df$features]
    MetaGSE = cbind(sample,data)
    
    
    ## GSE96804
    #读取表达数据
    rt = read.csv("E:/本地生信项目/DKD/RNAseq/Validation_Cohort/GSE96804_expr_normed.csv", row.names = 1)
    rt = as.matrix(rt)
    #读取分组信息
    sample = read.csv("E:/本地生信项目/DKD/RNAseq/Validation_Cohort/GSE96804_groupinfo.csv", row.names = 1)
    sample$ID <- rownames(sample)
    data = rt[,sample$ID]
    #数据转置
    data = t(data)#行为样本，列为基因
    #提取特征基因的表达矩阵
    data <- data[,colnames(data) %in% fea_df$features]
    GSE96804 = cbind(sample,data)
    
    
    ## GSE30528
    #读取表达数据
    rt = read.csv("E:/本地生信项目/DKD/RNAseq/Validation_Cohort/GSE30528_expr_normed.csv", row.names = 1)
    rt = as.matrix(rt)
    #读取分组信息
    sample = read.csv("E:/本地生信项目/DKD/RNAseq/Validation_Cohort/GSE30528_groupinfo.csv", row.names = 1)
    sample$ID <- rownames(sample)
    data = rt[,sample$ID]
    #数据转置
    data = t(data)#行为样本，列为基因
    #提取特征基因的表达矩阵
    data <- data[,colnames(data) %in% fea_df$features]
    GSE30528 = cbind(sample,data)
    
    
    ## GSE30529
    #读取表达数据
    rt = read.csv("E:/本地生信项目/DKD/RNAseq/Validation_Cohort/GSE30529_expr_normed.csv", row.names = 1)
    rt = as.matrix(rt)
    #读取分组信息
    sample = read.csv("E:/本地生信项目/DKD/RNAseq/Validation_Cohort/GSE30529_groupinfo.csv", row.names = 1)
    sample$ID <- rownames(sample)
    data = rt[,sample$ID]
    #数据转置
    data = t(data)#行为样本，列为基因
    #提取特征基因的表达矩阵
    data <- data[,colnames(data) %in% fea_df$features]
    GSE30529 = cbind(sample,data)
    
    fea_df
   
    ## CD163
    aflist1 = roc(Type~CD163 , data = MetaGSE)
    plot(aflist1, print.auc = F, auc.polygon = F, grid = TRUE, col = "#A6CEE3")
    aflist2 = roc(outcome~CD163 , data = GSE96804)
    plot(aflist2, add = T, print.auc = F, auc.polygon = F, grid = TRUE, col = "#1F78B4")
    aflist3 = roc(outcome~CD163 , data = GSE30528)
    plot(aflist3, add = T, print.auc = F, auc.polygon = F, grid = TRUE, col = "#B2DF8A")
    aflist4 = roc(outcome~CD163 , data = GSE30529)
    plot(aflist4, add = T, print.auc = F, auc.polygon = F, grid = TRUE, col = "#33A02C")
    
    round(auc(aflist1),3)##AUC

    round(auc(aflist2),3)##AUC

    round(auc(aflist3),3)##AUC

    round(auc(aflist4),3)##AUC

    text(0.5, 0.25, paste0("MetaGSE-AUC: ",sprintf("%.03f",round(auc(aflist1),3))), col = "#A6CEE3")
    text(0.5, 0.18, paste0("GSE96804-AUC: ",sprintf("%.03f",round(auc(aflist1),3))), col = "#1F78B4")
    text(0.5, 0.11, paste0("GSE30528-AUC: ",sprintf("%.03f",round(auc(aflist1),3))), col = "#B2DF8A")
    text(0.5, 0.03, paste0("GSE30529-AUC: ",sprintf("%.03f",round(auc(aflist1),3))), col = "#33A02C")
    
    

    ## CYBB
    aflist1 = roc(Type~CYBB , data = MetaGSE)
    plot(aflist1, print.auc = F, auc.polygon = F, grid = TRUE, col = "#A6CEE3")
    aflist2 = roc(outcome~CYBB , data = GSE96804)
    plot(aflist2, add = T, print.auc = F, auc.polygon = F, grid = TRUE, col = "#1F78B4")
    aflist3 = roc(outcome~CYBB , data = GSE30528)
    plot(aflist3, add = T, print.auc = F, auc.polygon = F, grid = TRUE, col = "#B2DF8A")
    aflist4 = roc(outcome~CYBB , data = GSE30529)
    plot(aflist4, add = T, print.auc = F, auc.polygon = F, grid = TRUE, col = "#33A02C")
    
    round(auc(aflist1),3)##AUC
    
    round(auc(aflist2),3)##AUC
    
    round(auc(aflist3),3)##AUC
    
    round(auc(aflist4),3)##AUC
    
    text(0.5, 0.25, paste0("MetaGSE-AUC: ",sprintf("%.03f",round(auc(aflist1),3))), col = "#A6CEE3")
    text(0.5, 0.18, paste0("GSE96804-AUC: ",sprintf("%.03f",round(auc(aflist1),3))), col = "#1F78B4")
    text(0.5, 0.11, paste0("GSE30528-AUC: ",sprintf("%.03f",round(auc(aflist1),3))), col = "#B2DF8A")
    text(0.5, 0.03, paste0("GSE30529-AUC: ",sprintf("%.03f",round(auc(aflist1),3))), col = "#33A02C")
    
    
    ## ELF3
    aflist1 = roc(Type~ELF3 , data = MetaGSE)
    plot(aflist1, print.auc = F, auc.polygon = F, grid = TRUE, col = "#A6CEE3")
    aflist2 = roc(outcome~ELF3 , data = GSE96804)
    plot(aflist2, add = T, print.auc = F, auc.polygon = F, grid = TRUE, col = "#1F78B4")
    aflist3 = roc(outcome~ELF3 , data = GSE30528)
    plot(aflist3, add = T, print.auc = F, auc.polygon = F, grid = TRUE, col = "#B2DF8A")
    aflist4 = roc(outcome~ELF3 , data = GSE30529)
    plot(aflist4, add = T, print.auc = F, auc.polygon = F, grid = TRUE, col = "#33A02C")
    
    round(auc(aflist1),3)##AUC
    
    round(auc(aflist2),3)##AUC
    
    round(auc(aflist3),3)##AUC
    
    round(auc(aflist4),3)##AUC
    
    text(0.5, 0.25, paste0("MetaGSE-AUC: ",sprintf("%.03f",round(auc(aflist1),3))), col = "#A6CEE3")
    text(0.5, 0.18, paste0("GSE96804-AUC: ",sprintf("%.03f",round(auc(aflist1),3))), col = "#1F78B4")
    text(0.5, 0.11, paste0("GSE30528-AUC: ",sprintf("%.03f",round(auc(aflist1),3))), col = "#B2DF8A")
    text(0.5, 0.03, paste0("GSE30529-AUC: ",sprintf("%.03f",round(auc(aflist1),3))), col = "#33A02C")
    
    
    ## FCN1
    aflist1 = roc(Type~FCN1 , data = MetaGSE)
    plot(aflist1, print.auc = F, auc.polygon = F, grid = TRUE, col = "#A6CEE3")
    aflist2 = roc(outcome~FCN1 , data = GSE96804)
    plot(aflist2, add = T, print.auc = F, auc.polygon = F, grid = TRUE, col = "#1F78B4")
    aflist3 = roc(outcome~FCN1 , data = GSE30528)
    plot(aflist3, add = T, print.auc = F, auc.polygon = F, grid = TRUE, col = "#B2DF8A")
    aflist4 = roc(outcome~FCN1 , data = GSE30529)
    plot(aflist4, add = T, print.auc = F, auc.polygon = F, grid = TRUE, col = "#33A02C")
    
    round(auc(aflist1),3)##AUC
    
    round(auc(aflist2),3)##AUC
    
    round(auc(aflist3),3)##AUC
    
    round(auc(aflist4),3)##AUC
    
    text(0.5, 0.25, paste0("MetaGSE-AUC: ",sprintf("%.03f",round(auc(aflist1),3))), col = "#A6CEE3")
    text(0.5, 0.18, paste0("GSE96804-AUC: ",sprintf("%.03f",round(auc(aflist1),3))), col = "#1F78B4")
    text(0.5, 0.11, paste0("GSE30528-AUC: ",sprintf("%.03f",round(auc(aflist1),3))), col = "#B2DF8A")
    text(0.5, 0.03, paste0("GSE30529-AUC: ",sprintf("%.03f",round(auc(aflist1),3))), col = "#33A02C")
    
    
    ## GPR65
    aflist1 = roc(Type~GPR65 , data = MetaGSE)
    plot(aflist1, print.auc = F, auc.polygon = F, grid = TRUE, col = "#A6CEE3")
    aflist2 = roc(outcome~GPR65 , data = GSE96804)
    plot(aflist2, add = T, print.auc = F, auc.polygon = F, grid = TRUE, col = "#1F78B4")
    aflist3 = roc(outcome~GPR65 , data = GSE30528)
    plot(aflist3, add = T, print.auc = F, auc.polygon = F, grid = TRUE, col = "#B2DF8A")
    aflist4 = roc(outcome~GPR65 , data = GSE30529)
    plot(aflist4, add = T, print.auc = F, auc.polygon = F, grid = TRUE, col = "#33A02C")
    
    round(auc(aflist1),3)##AUC
    
    round(auc(aflist2),3)##AUC
    
    round(auc(aflist3),3)##AUC
    
    round(auc(aflist4),3)##AUC
    
    text(0.5, 0.25, paste0("MetaGSE-AUC: ",sprintf("%.03f",round(auc(aflist1),3))), col = "#A6CEE3")
    text(0.5, 0.18, paste0("GSE96804-AUC: ",sprintf("%.03f",round(auc(aflist1),3))), col = "#1F78B4")
    text(0.5, 0.11, paste0("GSE30528-AUC: ",sprintf("%.03f",round(auc(aflist1),3))), col = "#B2DF8A")
    text(0.5, 0.03, paste0("GSE30529-AUC: ",sprintf("%.03f",round(auc(aflist1),3))), col = "#33A02C")
    
    
    ## LCN2
    aflist1 = roc(Type~LCN2 , data = MetaGSE)
    plot(aflist1, print.auc = F, auc.polygon = F, grid = TRUE, col = "#A6CEE3")
    aflist2 = roc(outcome~LCN2 , data = GSE96804)
    plot(aflist2, add = T, print.auc = F, auc.polygon = F, grid = TRUE, col = "#1F78B4")
    aflist3 = roc(outcome~LCN2 , data = GSE30528)
    plot(aflist3, add = T, print.auc = F, auc.polygon = F, grid = TRUE, col = "#B2DF8A")
    aflist4 = roc(outcome~LCN2 , data = GSE30529)
    plot(aflist4, add = T, print.auc = F, auc.polygon = F, grid = TRUE, col = "#33A02C")
    
    round(auc(aflist1),3)##AUC
    
    round(auc(aflist2),3)##AUC
    
    round(auc(aflist3),3)##AUC
    
    round(auc(aflist4),3)##AUC
    
    text(0.5, 0.25, paste0("MetaGSE-AUC: ",sprintf("%.03f",round(auc(aflist1),3))), col = "#A6CEE3")
    text(0.5, 0.18, paste0("GSE96804-AUC: ",sprintf("%.03f",round(auc(aflist1),3))), col = "#1F78B4")
    text(0.5, 0.11, paste0("GSE30528-AUC: ",sprintf("%.03f",round(auc(aflist1),3))), col = "#B2DF8A")
    text(0.5, 0.03, paste0("GSE30529-AUC: ",sprintf("%.03f",round(auc(aflist1),3))), col = "#33A02C")
    
    
    ## LTF
    aflist1 = roc(Type~LTF , data = MetaGSE)
    plot(aflist1, print.auc = F, auc.polygon = F, grid = TRUE, col = "#A6CEE3")
    aflist2 = roc(outcome~LTF , data = GSE96804)
    plot(aflist2, add = T, print.auc = F, auc.polygon = F, grid = TRUE, col = "#1F78B4")
    aflist3 = roc(outcome~LTF , data = GSE30528)
    plot(aflist3, add = T, print.auc = F, auc.polygon = F, grid = TRUE, col = "#B2DF8A")
    aflist4 = roc(outcome~LTF , data = GSE30529)
    plot(aflist4, add = T, print.auc = F, auc.polygon = F, grid = TRUE, col = "#33A02C")
    
    round(auc(aflist1),3)##AUC
    
    round(auc(aflist2),3)##AUC
    
    round(auc(aflist3),3)##AUC
    
    round(auc(aflist4),3)##AUC
    
    text(0.5, 0.25, paste0("MetaGSE-AUC: ",sprintf("%.03f",round(auc(aflist1),3))), col = "#A6CEE3")
    text(0.5, 0.18, paste0("GSE96804-AUC: ",sprintf("%.03f",round(auc(aflist1),3))), col = "#1F78B4")
    text(0.5, 0.11, paste0("GSE30528-AUC: ",sprintf("%.03f",round(auc(aflist1),3))), col = "#B2DF8A")
    text(0.5, 0.03, paste0("GSE30529-AUC: ",sprintf("%.03f",round(auc(aflist1),3))), col = "#33A02C")
    
    
    ## PROM1
    aflist1 = roc(Type~PROM1 , data = MetaGSE)
    plot(aflist1, print.auc = F, auc.polygon = F, grid = TRUE, col = "#A6CEE3")
    aflist2 = roc(outcome~PROM1 , data = GSE96804)
    plot(aflist2, add = T, print.auc = F, auc.polygon = F, grid = TRUE, col = "#1F78B4")
    aflist3 = roc(outcome~PROM1 , data = GSE30528)
    plot(aflist3, add = T, print.auc = F, auc.polygon = F, grid = TRUE, col = "#B2DF8A")
    aflist4 = roc(outcome~PROM1 , data = GSE30529)
    plot(aflist4, add = T, print.auc = F, auc.polygon = F, grid = TRUE, col = "#33A02C")
    
    round(auc(aflist1),3)##AUC
    
    round(auc(aflist2),3)##AUC
    
    round(auc(aflist3),3)##AUC
    
    round(auc(aflist4),3)##AUC
    
    text(0.5, 0.25, paste0("MetaGSE-AUC: ",sprintf("%.03f",round(auc(aflist1),3))), col = "#A6CEE3")
    text(0.5, 0.18, paste0("GSE96804-AUC: ",sprintf("%.03f",round(auc(aflist1),3))), col = "#1F78B4")
    text(0.5, 0.11, paste0("GSE30528-AUC: ",sprintf("%.03f",round(auc(aflist1),3))), col = "#B2DF8A")
    text(0.5, 0.03, paste0("GSE30529-AUC: ",sprintf("%.03f",round(auc(aflist1),3))), col = "#33A02C")
    
    
    ## S100A4
    aflist1 = roc(Type~S100A4 , data = MetaGSE)
    plot(aflist1, print.auc = F, auc.polygon = F, grid = TRUE, col = "#A6CEE3")
    aflist2 = roc(outcome~S100A4 , data = GSE96804)
    plot(aflist2, add = T, print.auc = F, auc.polygon = F, grid = TRUE, col = "#1F78B4")
    aflist3 = roc(outcome~S100A4 , data = GSE30528)
    plot(aflist3, add = T, print.auc = F, auc.polygon = F, grid = TRUE, col = "#B2DF8A")
    aflist4 = roc(outcome~S100A4 , data = GSE30529)
    plot(aflist4, add = T, print.auc = F, auc.polygon = F, grid = TRUE, col = "#33A02C")
    
    round(auc(aflist1),3)##AUC
    
    round(auc(aflist2),3)##AUC
    
    round(auc(aflist3),3)##AUC
    
    round(auc(aflist4),3)##AUC
    
    text(0.5, 0.25, paste0("MetaGSE-AUC: ",sprintf("%.03f",round(auc(aflist1),3))), col = "#A6CEE3")
    text(0.5, 0.18, paste0("GSE96804-AUC: ",sprintf("%.03f",round(auc(aflist1),3))), col = "#1F78B4")
    text(0.5, 0.11, paste0("GSE30528-AUC: ",sprintf("%.03f",round(auc(aflist1),3))), col = "#B2DF8A")
    text(0.5, 0.03, paste0("GSE30529-AUC: ",sprintf("%.03f",round(auc(aflist1),3))), col = "#33A02C")
    
    
    ## SOX4
    aflist1 = roc(Type~SOX4 , data = MetaGSE)
    plot(aflist1, print.auc = F, auc.polygon = F, grid = TRUE, col = "#A6CEE3")
    aflist2 = roc(outcome~SOX4 , data = GSE96804)
    plot(aflist2, add = T, print.auc = F, auc.polygon = F, grid = TRUE, col = "#1F78B4")
    aflist3 = roc(outcome~SOX4 , data = GSE30528)
    plot(aflist3, add = T, print.auc = F, auc.polygon = F, grid = TRUE, col = "#B2DF8A")
    aflist4 = roc(outcome~SOX4 , data = GSE30529)
    plot(aflist4, add = T, print.auc = F, auc.polygon = F, grid = TRUE, col = "#33A02C")
    
    round(auc(aflist1),3)##AUC
    
    round(auc(aflist2),3)##AUC
    
    round(auc(aflist3),3)##AUC
    
    round(auc(aflist4),3)##AUC
    
    text(0.5, 0.25, paste0("MetaGSE-AUC: ",sprintf("%.03f",round(auc(aflist1),3))), col = "#A6CEE3")
    text(0.5, 0.18, paste0("GSE96804-AUC: ",sprintf("%.03f",round(auc(aflist1),3))), col = "#1F78B4")
    text(0.5, 0.11, paste0("GSE30528-AUC: ",sprintf("%.03f",round(auc(aflist1),3))), col = "#B2DF8A")
    text(0.5, 0.03, paste0("GSE30529-AUC: ",sprintf("%.03f",round(auc(aflist1),3))), col = "#33A02C")
    
    
    
    ## TGFBI
    aflist1 = roc(Type~TGFBI , data = MetaGSE)
    plot(aflist1, print.auc = F, auc.polygon = F, grid = TRUE, col = "#A6CEE3")
    aflist2 = roc(outcome~TGFBI , data = GSE96804)
    plot(aflist2, add = T, print.auc = F, auc.polygon = F, grid = TRUE, col = "#1F78B4")
    aflist3 = roc(outcome~TGFBI , data = GSE30528)
    plot(aflist3, add = T, print.auc = F, auc.polygon = F, grid = TRUE, col = "#B2DF8A")
    aflist4 = roc(outcome~TGFBI , data = GSE30529)
    plot(aflist4, add = T, print.auc = F, auc.polygon = F, grid = TRUE, col = "#33A02C")
    
    round(auc(aflist1),3)##AUC
    
    round(auc(aflist2),3)##AUC
    
    round(auc(aflist3),3)##AUC
    
    round(auc(aflist4),3)##AUC
    
    text(0.5, 0.25, paste0("MetaGSE-AUC: ",sprintf("%.03f",round(auc(aflist1),3))), col = "#A6CEE3")
    text(0.5, 0.18, paste0("GSE96804-AUC: ",sprintf("%.03f",round(auc(aflist1),3))), col = "#1F78B4")
    text(0.5, 0.11, paste0("GSE30528-AUC: ",sprintf("%.03f",round(auc(aflist1),3))), col = "#B2DF8A")
    text(0.5, 0.03, paste0("GSE30529-AUC: ",sprintf("%.03f",round(auc(aflist1),3))), col = "#33A02C")
    
    
    ## TNFAIP8
    aflist1 = roc(Type~TNFAIP8 , data = MetaGSE)
    plot(aflist1, print.auc = F, auc.polygon = F, grid = TRUE, col = "#A6CEE3")
    aflist2 = roc(outcome~TNFAIP8 , data = GSE96804)
    plot(aflist2, add = T, print.auc = F, auc.polygon = F, grid = TRUE, col = "#1F78B4")
    aflist3 = roc(outcome~TNFAIP8 , data = GSE30528)
    plot(aflist3, add = T, print.auc = F, auc.polygon = F, grid = TRUE, col = "#B2DF8A")
    aflist4 = roc(outcome~TNFAIP8 , data = GSE30529)
    plot(aflist4, add = T, print.auc = F, auc.polygon = F, grid = TRUE, col = "#33A02C")
    
    round(auc(aflist1),3)##AUC
    
    round(auc(aflist2),3)##AUC
    
    round(auc(aflist3),3)##AUC
    
    round(auc(aflist4),3)##AUC
    
    text(0.5, 0.25, paste0("MetaGSE-AUC: ",sprintf("%.03f",round(auc(aflist1),3))), col = "#A6CEE3")
    text(0.5, 0.18, paste0("GSE96804-AUC: ",sprintf("%.03f",round(auc(aflist1),3))), col = "#1F78B4")
    text(0.5, 0.11, paste0("GSE30528-AUC: ",sprintf("%.03f",round(auc(aflist1),3))), col = "#B2DF8A")
    text(0.5, 0.03, paste0("GSE30529-AUC: ",sprintf("%.03f",round(auc(aflist1),3))), col = "#33A02C")
    
    
    ## calibration curve绘制
    library(rms)
    table(MetaGSE$Type)
    
    fea_df
    
    # 将Type一列的DKD变为1，其余为0
    MetaGSE$Type <- ifelse(MetaGSE$Type == "DKD", 1,0)
    table(MetaGSE$Type)
    dd <- datadist(MetaGSE)
    options(datadist="dd")
    
    fea_df
    
    fit2 <- lrm(Type ~ CD163+CYBB+ELF3+FCN1+GPR65+LCN2+LTF+PROM1+S100A4+SOX4+TGFBI+TNFAIP8,
                data = MetaGSE,x=T,y=T)
    cal2 <- calibrate(fit2, method='boot', B=1000)
    
    # 可视化
    plot(cal2,
         xlim = c(0,1),
         ylim = c(0,1),
         xlab = "Prediced Probability",
         ylab = "Observed Probability",
         cex.lab=1.2, cex.axis=1, cex.main=1.2, cex.sub=0.8,
         #subtitles = FALSE,
         legend = FALSE
    ) 
    
    lines(cal2[,c("predy","calibrated.corrected")], 
          type = 'l', #连线的类型，可以是"p","b","o"
          lwd = 3, #连线的粗细
          pch = 16, #点的形状，可以是0-20
          col = "#1D75B5") #连线的颜色
    
    lines(cal2[,c("predy","calibrated.orig")],type="l",pch=16,lwd=3,col="#ED7C72")
    
    abline(0,1,
           lty = 2, #对角线为虚线
           lwd = 2, #对角线的粗细
           col = "#224444") #对角线的颜色
    legend(0.6,0.2,
           c("Apparent","Bias-corrected","Ideal"), 
           lty = c(2,1,1), 
           lwd = c(2,3,3), 
           col = c("black","#1D75B5","#ED7C72"), 
           bty = "n")

    ## 决策曲线DCAplot绘制
    #这里我们一共绘制十二个基因的汇总DCAPlot
    library(rmda)
    
    # 将Type一列的DKD变为1，其余为0
    MetaGSE$Type <- ifelse(MetaGSE$Type == "DKD", 1,0)
    Meta <- str(MetaGSE)
    
    # CD163
    fit1 <- decision_curve(Type ~ CD163, # R语言里常见的公式类型
                           data = MetaGSE, 
                           bootstraps = 50 # 重抽样次数
    )
    
    # CYBB
    fit2 <- decision_curve(Type ~ CYBB, # R语言里常见的公式类型
                           data = MetaGSE, 
                           bootstraps = 50 # 重抽样次数
    )
    
    # ELF3
    fit3 <- decision_curve(Type ~ ELF3, # R语言里常见的公式类型
                           data = MetaGSE, 
                           bootstraps = 50 # 重抽样次数
    )
    
    # FCN1
    fit4 <- decision_curve(Type ~ FCN1, # R语言里常见的公式类型
                           data = MetaGSE, 
                           bootstraps = 50 # 重抽样次数
    )
    
    # GPR65
    fit5 <- decision_curve(Type ~ GPR65, # R语言里常见的公式类型
                           data = MetaGSE, 
                           bootstraps = 50 # 重抽样次数
    )
    
    # LCN2
    fit6 <- decision_curve(Type ~ LCN2, # R语言里常见的公式类型
                           data = MetaGSE, 
                           bootstraps = 50 # 重抽样次数
    )
    
    # LTF
    fit7 <- decision_curve(Type ~ LTF, # R语言里常见的公式类型
                           data = MetaGSE, 
                           bootstraps = 50 # 重抽样次数
    )
    
    # PROM1
    fit8 <- decision_curve(Type ~ PROM1, # R语言里常见的公式类型
                           data = MetaGSE, 
                           bootstraps = 50 # 重抽样次数
    )
    
    # S100A4
    fit9 <- decision_curve(Type ~ S100A4, # R语言里常见的公式类型
                           data = MetaGSE, 
                           bootstraps = 50 # 重抽样次数
    )
    
    # SOX4
    fit10 <- decision_curve(Type ~ SOX4, # R语言里常见的公式类型
                           data = MetaGSE, 
                           bootstraps = 50 # 重抽样次数
    )
    
    # TGFBI
    fit11 <- decision_curve(Type ~ TGFBI, # R语言里常见的公式类型
                           data = MetaGSE, 
                           bootstraps = 50 # 重抽样次数
    )
    
    # TNFAIP8
    fit12 <- decision_curve(Type ~ TNFAIP8, # R语言里常见的公式类型
                           data = MetaGSE, 
                           bootstraps = 50 # 重抽样次数
    )
    
    # 总和
    fit13 <- decision_curve(Type ~ CD163+CYBB+ELF3+FCN1+GPR65+LCN2+LTF+PROM1+S100A4+SOX4+TGFBI+TNFAIP8, # R语言里常见的公式类型
                            data = MetaGSE, 
                            bootstraps = 50 # 重抽样次数
    )
    
    # 绘制DCA曲线
    plot_decision_curve( list(fit1, fit2,fit3,fit4,fit5,fit6,fit7,fit8,fit9,fit10,fit11,fit12,fit13),
                         curve.names = c("CD163", "CYBB","ELF3","FCN1","GPR65","LCN2","LTF","PROM1","S100A4","SOX4","TGFBI","TNFAIP8","Normogram"),
                         col = c("#A6CEE3",
                                 "#1F78B4",
                                 "#B2DF8A",
                                 "#33A02C",
                                 
                                 "#8ECFC9",
                                 "#FFBE7A",
                                 "#FA7F6F",
                                 "#82B0D2",
                                 
                                 "#2878B5",
                                 "#9AC9DB",
                                 "#C82423",
                                 "#FF8884",
                                 "#4F70AE",
                                 
                                 "grey",
                                 "black"),
                         confidence.intervals = FALSE,  #remove confidence intervals
                         cost.benefit.axis = FALSE, #remove cost benefit axis
                         legend.position = "topright") #add the legend
    
    
    ### 特征基因在各集合的DKD与Normal中的表达量
 
    ### MetaGSE
    #读取表达数据
    rt = read.csv("E:/本地生信项目/DKD/RNAseq/Train_Cohort/exprmatrix_remove_batch_remove.csv", row.names = 1)
    #提取仅含有差异基因的表达矩阵
    rt <- rt[rownames(rt) %in% fea_df$features,]
    # write.csv(rt,"./MetaGSE_feadf_expr.csv") # 保存数据
    bar_mat <- t(rt)
    #读取分组信息
    sample = read.csv("E:/本地生信项目/DKD/RNAseq/Train_Cohort/groupinfo.csv", row.names = 1)
    sample$Sample <- rownames(sample)
    colnames(sample) <- c("Type","Tissue","Dataset","Sample")
    # write.csv(sample,"./MetaGSE_feadf_groupinfo.csv") # 保存数据
    sample$type2 <- sample$Type
    
    anno <- sample
    
    anno <- anno[rownames(bar_mat),]
    bar_mat <- bar_mat[rownames(anno),]
    bar_mat <- as.data.frame(bar_mat)
    bar_mat$sam = anno$Type
    
    #绘制
    library(RColorBrewer)
    library(ggpubr)
    library(ggplot2)
    
    table(bar_mat$sam)
    
    # 因子水平
    bar_mat$sam <- factor(bar_mat$sam,levels=c("DKD","Normal"))
    # 颜色、分组比较设置
    color <-c("#1D75B5","#ED7C72")
    my_comparisons <- list(c("DKD","Normal"))
    
    # 提取需要循环绘制的基因名
    gc <- colnames(bar_mat)
    #开始批量绘制
    plist<-list()
    for (i in 1:length(gc)){
      bar_tmp<-bar_mat[,c(gc[i],"sam")]
      colnames(bar_tmp)<-c("Expression","sam")
      pb1<-ggboxplot(bar_tmp,
                     x="sam",
                     y="Expression",
                     color="sam",
                     fill=NULL,
                     add = "jitter",
                     bxp.errorbar.width = 0.6,
                     width = 0.4,
                     size=0.01,
                     font.label = list(size=30), 
                     palette = color)+
        theme(panel.background =element_blank())
      
      pb1 <- pb1 + theme(axis.line=element_line(colour="black"))+theme(axis.title.x = element_blank())
      pb1 <- pb1 + theme(axis.title.y = element_blank())+theme(axis.text.x = element_text(size = 15,angle = 0,vjust = 1,hjust = 1))
      pb1 <- pb1 + theme(axis.text.y = element_text(size = 15))+ggtitle(gc[i])+theme(plot.title = element_text(hjust = 0.5,size=15,face="bold"))
      pb1 <- pb1 + theme(legend.position = "NA")
      pb1 <- pb1 + stat_compare_means(method="t.test",hide.ns = F,comparisons =my_comparisons,label="p.signif")+
        stat_boxplot(geom = "errorbar",
                     width=0.3, color = color)
      plist[[i]] <- pb1
    } 
    
    library(cowplot)
    pall <- plot_grid(plist[[1]],plist[[2]],plist[[3]],
                      plist[[4]],plist[[5]],plist[[6]],
                      plist[[7]],plist[[8]],plist[[9]],
                      plist[[10]],plist[[11]],plist[[12]],
                      ncol=3)
    pall
    
    

    ### GSE96804
    #读取表达数据
    rt = read.csv("E:/本地生信项目/DKD/RNAseq/Validation_Cohort/GSE96804_expr_normed.csv", row.names = 1)
    #提取仅含有差异基因的表达矩阵
    rt <- rt[rownames(rt) %in% fea_df$features,]
    # write.csv(rt,"./GSE96804_feadf_expr.csv") # 保存数据
    bar_mat <- t(rt)
    #读取分组信息
    sample = read.csv("E:/本地生信项目/DKD/RNAseq/Validation_Cohort/GSE96804_groupinfo.csv", row.names = 1)
    sample$Sample <- rownames(sample)
    # 将Outcome一列的1变为DKD，其余为Normal
    sample$outcome <- ifelse(sample$outcome == 1, "DKD","Normal")
    colnames(sample) <- c("Dataset","Type","Sample")
    # write.csv(sample,"./GSE96804_feadf_groupinfo.csv") # 保存数据
    sample$type2 <- sample$Type
    
    anno <- sample
    
    anno <- anno[rownames(bar_mat),]
    bar_mat <- bar_mat[rownames(anno),]
    bar_mat <- as.data.frame(bar_mat)
    bar_mat$sam = anno$Type
    
    #绘制
    library(RColorBrewer)
    library(ggpubr)
    library(ggplot2)
    
    table(bar_mat$sam)
    
    # 因子水平
    bar_mat$sam <- factor(bar_mat$sam,levels=c("DKD","Normal"))
    # 颜色、分组比较设置
    color <-c("#1D75B5","#ED7C72")
    my_comparisons <- list(c("DKD","Normal"))
    
    # 提取需要循环绘制的基因名
    gc <- colnames(bar_mat)
    #开始批量绘制
    plist<-list()
    for (i in 1:length(gc)){
      bar_tmp<-bar_mat[,c(gc[i],"sam")]
      colnames(bar_tmp)<-c("Expression","sam")
      pb1<-ggboxplot(bar_tmp,
                     x="sam",
                     y="Expression",
                     color="sam",
                     fill=NULL,
                     add = "jitter",
                     bxp.errorbar.width = 0.6,
                     width = 0.4,
                     size=0.01,
                     font.label = list(size=30), 
                     palette = color)+
        theme(panel.background =element_blank())
      
      pb1 <- pb1 + theme(axis.line=element_line(colour="black"))+theme(axis.title.x = element_blank())
      pb1 <- pb1 + theme(axis.title.y = element_blank())+theme(axis.text.x = element_text(size = 15,angle = 0,vjust = 1,hjust = 1))
      pb1 <- pb1 + theme(axis.text.y = element_text(size = 15))+ggtitle(gc[i])+theme(plot.title = element_text(hjust = 0.5,size=15,face="bold"))
      pb1 <- pb1 + theme(legend.position = "NA")
      pb1 <- pb1 + stat_compare_means(method="t.test",hide.ns = F,comparisons =my_comparisons,label="p.signif")+
        stat_boxplot(geom = "errorbar",
                     width=0.3, color = color)
      plist[[i]] <- pb1
    } 
    
    library(cowplot)
    pall <- plot_grid(plist[[1]],plist[[2]],plist[[3]],
                      plist[[4]],plist[[5]],plist[[6]],
                      plist[[7]],plist[[8]],plist[[9]],
                      plist[[10]],plist[[11]],plist[[12]],
                      ncol=3)
    pall
    
    
    
    ### GSE30528
    #读取表达数据
    rt = read.csv("E:/本地生信项目/DKD/RNAseq/Validation_Cohort/GSE30528_expr_normed.csv", row.names = 1)
    #提取仅含有差异基因的表达矩阵
    rt <- rt[rownames(rt) %in% fea_df$features,]
    # write.csv(rt,"./GSE30528_feadf_expr.csv") # 保存数据
    bar_mat <- t(rt)
    #读取分组信息
    sample = read.csv("E:/本地生信项目/DKD/RNAseq/Validation_Cohort/GSE30528_groupinfo.csv", row.names = 1)
    sample$Sample <- rownames(sample)
    # 将Outcome一列的1变为DKD，其余为Normal
    sample$outcome <- ifelse(sample$outcome == 1, "DKD","Normal")
    colnames(sample) <- c("Dataset","Type","Sample")
    # write.csv(sample,"./GSE30528_feadf_groupinfo.csv") # 保存数据
    sample$type2 <- sample$Type
    
    anno <- sample
    
    anno <- anno[rownames(bar_mat),]
    bar_mat <- bar_mat[rownames(anno),]
    bar_mat <- as.data.frame(bar_mat)
    bar_mat$sam = anno$Type
    
    #绘制
    library(RColorBrewer)
    library(ggpubr)
    library(ggplot2)
    
    table(bar_mat$sam)
    
    # 因子水平
    bar_mat$sam <- factor(bar_mat$sam,levels=c("DKD","Normal"))
    # 颜色、分组比较设置
    color <-c("#1D75B5","#ED7C72")
    my_comparisons <- list(c("DKD","Normal"))
    
    # 提取需要循环绘制的基因名
    gc <- colnames(bar_mat)
    #开始批量绘制
    plist<-list()
    for (i in 1:length(gc)){
      bar_tmp<-bar_mat[,c(gc[i],"sam")]
      colnames(bar_tmp)<-c("Expression","sam")
      pb1<-ggboxplot(bar_tmp,
                     x="sam",
                     y="Expression",
                     color="sam",
                     fill=NULL,
                     add = "jitter",
                     bxp.errorbar.width = 0.6,
                     width = 0.4,
                     size=0.01,
                     font.label = list(size=30), 
                     palette = color)+
        theme(panel.background =element_blank())
      
      pb1 <- pb1 + theme(axis.line=element_line(colour="black"))+theme(axis.title.x = element_blank())
      pb1 <- pb1 + theme(axis.title.y = element_blank())+theme(axis.text.x = element_text(size = 15,angle = 0,vjust = 1,hjust = 1))
      pb1 <- pb1 + theme(axis.text.y = element_text(size = 15))+ggtitle(gc[i])+theme(plot.title = element_text(hjust = 0.5,size=15,face="bold"))
      pb1 <- pb1 + theme(legend.position = "NA")
      pb1 <- pb1 + stat_compare_means(method="t.test",hide.ns = F,comparisons =my_comparisons,label="p.signif")+
        stat_boxplot(geom = "errorbar",
                     width=0.3, color = color)
      plist[[i]] <- pb1
    } 
    
    library(cowplot)
    pall <- plot_grid(plist[[1]],plist[[2]],plist[[3]],
                      plist[[4]],plist[[5]],plist[[6]],
                      plist[[7]],plist[[8]],plist[[9]],
                      plist[[10]],plist[[11]],plist[[12]],
                      ncol=3)
    pall
    
    
    
    ### GSE30529
    #读取表达数据
    rt = read.csv("E:/本地生信项目/DKD/RNAseq/Validation_Cohort/GSE30529_expr_normed.csv", row.names = 1)
    #提取仅含有差异基因的表达矩阵
    rt <- rt[rownames(rt) %in% fea_df$features,]
    # write.csv(rt,"./GSE30529_feadf_expr.csv") # 保存数据
    bar_mat <- t(rt)
    #读取分组信息
    sample = read.csv("E:/本地生信项目/DKD/RNAseq/Validation_Cohort/GSE30529_groupinfo.csv", row.names = 1)
    sample$Sample <- rownames(sample)
    # 将Outcome一列的1变为DKD，其余为Normal
    sample$outcome <- ifelse(sample$outcome == 1, "DKD","Normal")
    colnames(sample) <- c("Dataset","Type","Sample")
    # write.csv(sample,"./GSE30529_feadf_groupinfo.csv") # 保存数据
    sample$type2 <- sample$Type
    
    anno <- sample
    
    anno <- anno[rownames(bar_mat),]
    bar_mat <- bar_mat[rownames(anno),]
    bar_mat <- as.data.frame(bar_mat)
    bar_mat$sam = anno$Type
    
    #绘制
    library(RColorBrewer)
    library(ggpubr)
    library(ggplot2)
    
    table(bar_mat$sam)
    
    # 因子水平
    bar_mat$sam <- factor(bar_mat$sam,levels=c("DKD","Normal"))
    # 颜色、分组比较设置
    color <-c("#1D75B5","#ED7C72")
    my_comparisons <- list(c("DKD","Normal"))
    
    # 提取需要循环绘制的基因名
    gc <- colnames(bar_mat)
    #开始批量绘制
    plist<-list()
    for (i in 1:length(gc)){
      bar_tmp<-bar_mat[,c(gc[i],"sam")]
      colnames(bar_tmp)<-c("Expression","sam")
      pb1<-ggboxplot(bar_tmp,
                     x="sam",
                     y="Expression",
                     color="sam",
                     fill=NULL,
                     add = "jitter",
                     bxp.errorbar.width = 0.6,
                     width = 0.4,
                     size=0.01,
                     font.label = list(size=30), 
                     palette = color)+
        theme(panel.background =element_blank())
      
      pb1 <- pb1 + theme(axis.line=element_line(colour="black"))+theme(axis.title.x = element_blank())
      pb1 <- pb1 + theme(axis.title.y = element_blank())+theme(axis.text.x = element_text(size = 15,angle = 0,vjust = 1,hjust = 1))
      pb1 <- pb1 + theme(axis.text.y = element_text(size = 15))+ggtitle(gc[i])+theme(plot.title = element_text(hjust = 0.5,size=15,face="bold"))
      pb1 <- pb1 + theme(legend.position = "NA")
      pb1 <- pb1 + stat_compare_means(method="t.test",hide.ns = F,comparisons =my_comparisons,label="p.signif")+
        stat_boxplot(geom = "errorbar",
                     width=0.3, color = color)
      plist[[i]] <- pb1
    } 
    
    library(cowplot)
    pall <- plot_grid(plist[[1]],plist[[2]],plist[[3]],
                      plist[[4]],plist[[5]],plist[[6]],
                      plist[[7]],plist[[8]],plist[[9]],
                      plist[[10]],plist[[11]],plist[[12]],
                      ncol=3)
    pall
    
    
    