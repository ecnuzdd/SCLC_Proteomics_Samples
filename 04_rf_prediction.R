if(TRUE){
  ### work dictory ###
  BASE_DIR = "/home/ddzhan/SCLC/20210726-SX"
  BASE_DIR = normalizePath(BASE_DIR)
  setwd(BASE_DIR)
  seed = seed0 = seed1 = 129
  ### library ###
  library('readxl')
  library('xlsx')
  library('stringr')
  library("ctv")
  library("kernlab")
  library("survival")
  library("survminer")
  library("pracma")
  library('Seurat')
  library('ggplot2')
  library('sctransform')
  library('CancerSubtypes')
  library('ggsci')
  library('pheatmap')
  
  ### custom library ###
  
  custom_library_path = normalizePath(
    file.path(BASE_DIR, '..', 'library_directory')
  )
  source(
    normalizePath(
      file.path(custom_library_path, 'whole_gc_analysis_library.R')
    ), 
    encoding = 'UTF-8', 
    echo=TRUE
  )
  source(
    normalizePath(
      file.path(custom_library_path, 'km_survival_analysis_library.R')
    ), 
    encoding = 'UTF-8', 
    echo=TRUE
  )
  source(
    normalizePath(
      file.path(custom_library_path, 'nmf_cc_library.R')
    ), 
    encoding = 'UTF-8', 
    echo=TRUE
  )
  source(
    normalizePath(
      file.path(custom_library_path, 'nmf_cc_library2.R')
    ), 
    encoding = 'UTF-8', 
    echo=TRUE
  )
  source(
    normalizePath(
      file.path(custom_library_path, 'generate_feature_matrix_library.R')
    ), 
    encoding = 'UTF-8', 
    echo=TRUE
  )
  source(
    normalizePath(
      file.path(custom_library_path, 'compute_5year_survival_rate.R')
    ), 
    encoding = 'UTF-8', 
    echo=TRUE
  )
  
  #### QCDATA_DIR ####
  QCDATA_DIR = normalizePath(
    file.path(
      BASE_DIR, 
      'input'
    ), mustWork = F
  )
  clinical_data_path = normalizePath(
    file.path(
      QCDATA_DIR,
      'clinical_data_discovery.csv'
    )
  )
  clinical_data_all = read.csv(clinical_data_path) %>% as.data.frame()
  clinical_data_all = clinical_data_all[!is.na(clinical_data_all$Proteome_fno_T),]
  dim_ifot_data_path = normalizePath(
    file.path(QCDATA_DIR, 'dset_ifot_data_0.csv')
    # file.path(QCDATA_DIR, 'sclc_ifot_aft_qa.csv')
  )
  dim_ifot_data = read.csv(dim_ifot_data_path)
  
  redundant_gene_path = normalizePath(
    file.path(QCDATA_DIR, 'redundant_genes.csv')
  )
  redundant_genes = read.csv(redundant_gene_path)
  redundant_genes = as.vector(redundant_genes$redundant_genes)
  
  #### function ####
  get_deps_info_df <- function(edata, cdata){
    # edata = wh_edata1
    # cdata = wh_cdata
    
    group = cdata$Cluster
    all_data_mat = as.matrix(edata[,-1])
    g = factor(group, levels = sort(unique(group)))
    g.id = as.vector(table(g))
    gi_df = data.frame()
    gi_total = nrow(all_data_mat)
    for (gi in 1:gi_total) {
      x = all_data_mat[gi,]
      x.id = tapply(x, g, function(x){length(which(x>0))})
      x.mean = tapply(x, g, mean)
      
      x.mean.desc_order = order(x.mean, decreasing = T)
      
      g.id_1 = g.id[x.mean.desc_order]
      x.id_1 = x.id[x.mean.desc_order]
      x.id_1_idp = x.id_1/g.id_1
      x.mean_1 = x.mean[x.mean.desc_order]
      x.mean_1_label = names(x.mean_1)
      
      gi.idp = x.id_1_idp[1]
      gi.ratio = x.mean_1[1]/x.mean_1[2]
      x1 = x[g == x.mean_1_label[1]]
      x2 = x[g == x.mean_1_label[2]]
      # x12.test = t.test(log2(x1+1), log2(x2+1), paired = F)
      x12.test = wilcox.test((x1), (x2), paired = F)
      gi.pvalue = x12.test$p.value
      gi.symbol = edata$Symbol[gi]
      gi.slable = x.mean_1_label[1]
      
      coxph_data = data.frame(
        OS_time = cdata$OS_time,
        OS_status = cdata$OS_status,
        Expr = log1p(x)
      )
      cox_model = coxph(
        Surv(OS_time, OS_status) ~ Expr, 
        data = coxph_data
      )
      cox_model.summary = summary(cox_model)
      # plot(
      #   ggforest(cox_model, data = coxph_data, main = gi.symbol)
      # )
      gi.cox_pvalue = cox_model.summary$coefficients[5]
      gi.hr = cox_model.summary$coefficients[2]
      
      
      gi_vector = c(gi.symbol, gi.slable, gi.idp, gi.ratio, gi.pvalue, gi.cox_pvalue, gi.hr)
      
      gi_df = rbind(gi_df, gi_vector)
      
      if(gi%%500==0 | gi == gi_total){
        print(paste0(gi, ' / ', gi_total))
      }
    }
    gi_df$Mean = rowMeans(all_data_mat)
    colnames(gi_df) = c('Symbol', 'SLabel', 'SIDP', 'FC', 'Pvalue', 'Mean', 'CoxP', 'CoxHR')
    
    gi_df$SIDP = as.numeric(gi_df$SIDP)
    gi_df$FC = as.numeric(gi_df$FC)
    gi_df$Pvalue = as.numeric(gi_df$Pvalue)
    gi_df$Mean = as.numeric(gi_df$Mean)
    gi_df$CoxP = as.numeric(gi_df$CoxP)
    gi_df$CoxHR = as.numeric(gi_df$CoxHR)
    gi_df$adj.Pvalue = p.adjust(gi_df$Pvalue, method = 'fdr')
    return(gi_df)
  }
}

topN = 500

#### wh data #### 
wh_cdata =  read_excel(
  normalizePath(file.path(
    BASE_DIR, 'input', "SCLC_clinical_data_S75_C3_F445_plus_batch2_20210118_v1.3_jq.xlsx"
  )),
  sheet = 'Batch1'
)
wh_edata = read.csv(
  normalizePath(
    file.path(QCDATA_DIR, 'dset_ifot_data_0.csv')
  )
)


wh_edata.prot_count = apply(wh_edata[,-1], 2, function(x){length(which(x>0))})
print(quantile(wh_edata.prot_count, seq(0, 1, 0.1)))

wh_edata1 = wh_edata[-match(
  redundant_genes, wh_edata$Symbol
),]
wh_edata1_topN = fot_df_correction(get_topx_df_1(wh_edata1, topN))

# topN
wh_edata2 = data.frame(
  Symbol = wh_edata1_topN$Symbol,
  as.matrix(wh_edata1_topN[,-1])
)

# KM-Plot: Global
library("RColorBrewer")
col.palette = brewer.pal(n = 12, name = "Paired")
if(T){
  # *** KM curve for OS *** #
  survival_data = data.frame(
    # OS_time = survival_data_v$OS_time,
    OS_time = wh_cdata$OS_time,
    OS_status = wh_cdata$OS_status,
    Cluster = wh_cdata$Cluster,
    TNM_stage = wh_cdata$TNM_stage
  )
  survival_data = survival_data
  os_fit <- survfit(
    Surv(as.numeric(survival_data$OS_time), as.numeric(survival_data$OS_status)) ~ Cluster, data = survival_data
  )
  main_title = paste('Survival curve | N=', nrow(survival_data), sep = '')
  cluster_num = 3 
  os_fit_ggsurvplot <- ggsurvplot(
    os_fit, data = survival_data,
    title = main_title,
    surv.scale = 'percent',
    xlab = 'Months',
    legend = 'right', 
    legend.title = paste("NMFCluster", '', sep = ''),
    legend.labs = c(
      # 'NoChemo',
      paste('NMF', sort(unique(survival_data$Cluster)), sep = '')
    ),
    pval = TRUE,
    # palette = c('grey', col.palette, '#000000')[-1],
    palette = c('#157E3A', '#11499D', '#E62017'),
    # risk.table.col = "Group",
    risk.table = TRUE, 
    risk.table.y.text.col = TRUE,
    surv.median.line = 'hv',
    # break.time.by = 12,
    # xlim = c(0,12)
  )
  
  print(os_fit_ggsurvplot)
  
  os_pvalue = surv_pvalue(os_fit)$pval
  print(paste('os_p_value: ', os_pvalue, sep = ''))
  
  # table_5ysr = get_n_year_surv_rate_using_clinical_data(
  #   na.omit(survival_data), n = 5, cluster_label = names(table(survival_data$Cluster))
  # )
  # print(table_5ysr)
}

#### sx data #### 
sx_cdata = read_excel(
  normalizePath(file.path(
    BASE_DIR, 'input', "SX-SCLC-24-ClinicalData-20210722.xlsx"
  )),
  sheet = 1
)
sx_edata = read.csv(
  normalizePath(file.path(
    BASE_DIR, 'input', "SX-SCLC-24-Profiling-v1.csv"
  ))
)
colnames(sx_edata) = c(
  'Symbol', apply(data.frame(colnames(sx_edata)), 1, function(x){
    strsplit(x, split = '.', fixed = T)[[1]][1]
  })[-1]
)

sx_edata.prot_count = apply(sx_edata[,-1], 2, function(x){length(which(x>0))})
print(quantile(sx_edata.prot_count, seq(0, 1, 0.1)))

sx_edata1 = sx_edata[-na.omit(match(
  redundant_genes, sx_edata$Symbol
)),]
sx_edata1_topN = fot_df_correction(get_topx_df_1(sx_edata1, topN))

# topN
sx_edata2 = data.frame(
  Symbol = sx_edata1_topN$Symbol,
  as.matrix(sx_edata1_topN[,-1])
)


#### hn data #### 
hn_cdata = read_excel(
  normalizePath(file.path(
    BASE_DIR, 'input', "HN-SCLC-52-ClinicalData-20210916-v1.xlsx"
  )),
  sheet = 1
)
hn_edata = read.csv(
  normalizePath(file.path(
    BASE_DIR, 'input', "HN-SCLC-52-Profiling.csv"
  ))
)
colnames(hn_edata) = c(
  'Symbol', apply(data.frame(colnames(hn_edata)), 1, function(x){
    strsplit(x, split = '.', fixed = T)[[1]][1]
  })[-1]
)

hn_edata.prot_count = apply(hn_edata[,-1], 2, function(x){length(which(x>0))})
print(quantile(hn_edata.prot_count, seq(0, 1, 0.1)))

hn_edata1 = hn_edata[-na.omit(match(
  redundant_genes, hn_edata$Symbol
)),]
hn_edata1_topN = fot_df_correction(get_topx_df_1(hn_edata1, topN))

# topN
hn_edata2 = data.frame(
  Symbol = hn_edata1_topN$Symbol,
  as.matrix(hn_edata1_topN[,-1])
)

# prot_count distribution
df = data.frame(
  ProtCount = c(wh_edata.prot_count, sx_edata.prot_count, hn_edata.prot_count),
  Site = c(
    rep('WH', length(wh_edata.prot_count)),
    rep('SX', length(sx_edata.prot_count)),
    rep('HN', length(hn_edata.prot_count))
  )
)
df$Site = factor(df$Site, levels = c('WH', 'SX', 'HN'))
p <- ggplot(df, mapping = aes(x=Site, y=ProtCount, fill = Site)) # fill = Site
# p <- p + geom_violin(trim=FALSE) 
p <- p + geom_boxplot(width=0.5)
# p <- p + geom_dotplot(binaxis='y', stackdir='center', position=position_dodge(1))
p <- p + scale_fill_brewer(palette="RdBu")
p <- p + theme_classic2()
print(p)



# deps
# group = wh_cdata$Cluster
# edata = wh_edata1
wh_edata_model = wh_edata2
sx_edata_model = sx_edata2
hn_edata_model = hn_edata2


wh_deps_info_df = get_deps_info_df(wh_edata_model, wh_cdata)

# features
# wh_deps_index = which(
#   wh_deps_info_df$SIDP > 0.25
#   & wh_deps_info_df$FC > 1.5
#   & wh_deps_info_df$adj.Pvalue < 0.25
#   # & wh_deps_info_df$Pvalue < 0.1
#   # & wh_deps_info_df$CoxP < 0.2
# )
wh_deps_index = which(
  wh_deps_info_df$SIDP > 0.25
  & wh_deps_info_df$FC > 1.5
  & wh_deps_info_df$adj.Pvalue < 0.25
  # & wh_deps_info_df$Pvalue < 0.1
  # & wh_deps_info_df$CoxP < 0.2
)

wh_deps_info_df_subset = wh_deps_info_df[wh_deps_index,]
write.csv(
  wh_deps_info_df_subset,
  'SCLC-75-Deps-58-Classifier-20220222.csv',
  row.names = F
)


print(length(wh_deps_index))
wh_features = wh_deps_info_df$Symbol[wh_deps_index]
wh_edata_model.dep_mat = as.matrix(wh_edata_model[na.omit(match(wh_features, wh_edata_model$Symbol)),-1])


ph_cluster = wh_cdata$Cluster[order(wh_cdata$Cluster)]
ph_mat = wh_edata_model.dep_mat[,order(wh_cdata$Cluster)]
annotation_col = data.frame(
  Cluster = ph_cluster
)
rownames(annotation_col) = colnames(ph_mat)
pheatmap(
  ph_mat,
  scale = 'row',
  cluster_cols = F,
  annotation_col = annotation_col,
  col = colorRampPalette(c( "green2","green1","green","black","red","red1","red2"))(100),
  # clustering_distance_rows = 'euclidean',
  main = paste0('Features = ', nrow(ph_mat)),
  clustering_method = 'ward.D2',
  border_color=NA,
  show_rownames = F,
  show_colnames = F
)

# build model
features = wh_features
wh_train_mat  = t(as.matrix(wh_edata_model[na.omit(match(features, wh_edata_model$Symbol)),-1]))
wh_train_mat = log1p(wh_train_mat)

wh_train_mat.mean = apply(wh_train_mat, 2, mean)
wh_train_mat.sd = apply(wh_train_mat, 2, sd)

# wh_train_mat = t(scale(t(wh_train_mat)))

colnames(wh_train_mat) = features
range(wh_train_mat)
dim(wh_train_mat)


wh_train_label = as.factor(wh_cdata$Cluster)
library('caret')
set.seed(11)
model_1 <- train(
  x = wh_train_mat, # Survived is a function of the variables we decided to include
  y = wh_train_label,
  # data = train, # Use the train data frame as the training data
  method = 'rf',# Use the 'random forest' algorithm
  trControl = trainControl(
    method = 'cv', # Use cross-validation
    number = 10
  ) # Use 10 folds for cross-validation
)
print(model_1)





# sx prediction
sx_pred_mat = as.matrix(sx_edata_model[(match(features, sx_edata_model$Symbol)),-1])
print(length(which(is.na(match(features, sx_edata_model$Symbol)))))
sx_pred_mat[is.na(sx_pred_mat)] = 0
rownames(sx_pred_mat) = features
sx_pred_mat = t(log1p(sx_pred_mat))

# for (i in 1:ncol(sx_pred_mat)) {
#   sx_pred_mat[,i] = (sx_pred_mat[,i] - wh_train_mat.mean[i])/wh_train_mat.sd[i]
# }


# sx_pred_mat = t(scale(t(sx_pred_mat)))

sx_pred_crs <- predict(model_1, newdata = sx_pred_mat)
print(table(sx_pred_crs))

# KM-Plot: Global
library("RColorBrewer")
col.palette = brewer.pal(n = 12, name = "Paired")
if(T){
  # *** KM curve for OS *** #
  survival_data = data.frame(
    # OS_time = survival_data_v$OS_time,
    OS_time = sx_cdata$OS_time_day/30,
    OS_status = sx_cdata$OS_status,
    Cluster = sx_pred_crs,
    TNM_stage = sx_cdata$TNM_stage
  )
  survival_data = survival_data
  os_fit <- survfit(
    Surv(as.numeric(survival_data$OS_time), as.numeric(survival_data$OS_status)) ~ Cluster, data = survival_data
  )
  main_title = paste('Survival curve | N=', nrow(survival_data), ' | SX', sep = '')
  cluster_num = 3 
  os_fit_ggsurvplot <- ggsurvplot(
    os_fit, data = survival_data,
    title = main_title,
    surv.scale = 'percent',
    xlab = 'Months',
    legend = 'right', 
    legend.title = paste("Prediction", '', sep = ''),
    legend.labs = c(
      # 'NoChemo',
      paste('Pred', sort(unique(survival_data$Cluster)), sep = '')
    ),
    pval = TRUE,
    # palette = c('grey', col.palette, '#000000')[-1],
    palette = c('#157E3A', '#11499D', '#E62017'),
    # risk.table.col = "Group",
    risk.table = TRUE, 
    risk.table.y.text.col = TRUE,
    surv.median.line = 'hv',
    break.time.by = 12,
    # xlim = c(0,12)
  )
  
  print(os_fit_ggsurvplot)
  
  os_pvalue = surv_pvalue(os_fit)$pval
  print(paste('os_p_value: ', os_pvalue, sep = ''))
  
  # table_5ysr = get_n_year_surv_rate_using_clinical_data(
  #   na.omit(survival_data), n = 5, cluster_label = names(table(survival_data$Cluster))
  # )
  # print(table_5ysr)
}

# KM-Plot: Subgroups
if(F){
  yl_tmc_data = data.frame(
    OS_time = sx_cdata$OS_time_day/30,
    OS_status = sx_cdata$OS_status,
    Cluster = sx_pred_crs,
    Chemotherapy_status = sx_cdata$Chemo_status,
    TNM_stage = sx_cdata$TNM_stage
  )
  yl_tms_cluster = sort(unique(yl_tmc_data$Cluster))
  yl_tmc_count = length(yl_tms_cluster)
  for (j in 1:yl_tmc_count) {
    
    j.clinical_data = yl_tmc_data[yl_tmc_data$Cluster==yl_tms_cluster[j],]
    
    
    j.index = which(j.clinical_data$TNM_stage>1 & j.clinical_data$TNM_stage<4)
    clinical_data.subset = j.clinical_data[j.index,]
    
    
    if(j>0){
      # *** KM curve for OS *** #
      survival_data = data.frame(
        OS_time = clinical_data.subset$OS_time,
        OS_status = clinical_data.subset$OS_status,
        Cluster= clinical_data.subset$Chemotherapy_status
      )
      os_fit <- survfit(
        Surv(as.numeric(survival_data$OS_time), as.numeric(survival_data$OS_status)) ~ Cluster, data = survival_data
      )
      main_title = paste('Survival curve | TNM=II/III | Pred', yl_tms_cluster[j], sep = '')
      os_fit_ggsurvplot <- ggsurvplot(
        os_fit, data = survival_data,
        title = main_title,
        surv.scale = 'percent',
        xlab = 'Months',
        legend = 'right', 
        legend.title = "Chemo status",
        legend.labs = c('NoChemo', 'Chemo'),
        pval = TRUE,
        palette = c('grey', 'blue'),
        # risk.table.col = "Group",
        risk.table = TRUE, 
        risk.table.y.text.col = TRUE,
        surv.median.line = 'hv',
        break.time.by = 12,
        # xlim = c(0,12)
      )
      print(os_fit_ggsurvplot)
      os_pvalue = surv_pvalue(os_fit)$pval
      os_pvalue = round(os_pvalue, 3)
      print(paste('os_p_value: ', os_pvalue, sep = ''))
      
      
    }
    
    
  }
  
  
  
  
}






# hn prediction
hn_pred_mat = as.matrix(hn_edata_model[(match(features, hn_edata_model$Symbol)),-1])
print(length(which(is.na(match(features, hn_edata_model$Symbol)))))
hn_pred_mat[is.na(hn_pred_mat)] = 0
rownames(hn_pred_mat) = features
hn_pred_mat = t(log1p(hn_pred_mat))

# hn_pred_mat = t(scale(t(hn_pred_mat)))

hn_pred_crs <- predict(model_1, newdata = hn_pred_mat)
print(table(hn_pred_crs))

# KM-Plot: Global
library("RColorBrewer")
col.palette = brewer.pal(n = 12, name = "Paired")
if(T){
  # *** KM curve for OS *** #
  survival_data = data.frame(
    # OS_time = survival_data_v$OS_time,
    OS_time = hn_cdata$OS_time_day/30,
    OS_status = hn_cdata$OS_status,
    Cluster = hn_pred_crs,
    TNM_stage = 0
  )
  survival_data = survival_data# [survival_data$Cluster!='C1',]
  os_fit <- survfit(
    Surv(as.numeric(survival_data$OS_time), as.numeric(survival_data$OS_status)) ~ Cluster, data = survival_data
  )
  main_title = paste('Survival curve | N=', nrow(survival_data), ' | HN', sep = '')
  cluster_num = 3 
  os_fit_ggsurvplot <- ggsurvplot(
    os_fit, data = survival_data,
    title = main_title,
    surv.scale = 'percent',
    xlab = 'Months',
    legend = 'right', 
    legend.title = paste("Prediction", '', sep = ''),
    legend.labs = c(
      # 'NoChemo',
      paste('Pred', sort(unique(survival_data$Cluster)), sep = '')
    ),
    pval = TRUE,
    # palette = c('grey', col.palette, '#000000')[-1],
    palette = c('#157E3A', '#11499D', '#E62017'),
    # risk.table.col = "Group",
    risk.table = TRUE, 
    risk.table.y.text.col = TRUE,
    surv.median.line = 'hv',
    break.time.by = 12,
    # xlim = c(0,12)
  )
  
  print(os_fit_ggsurvplot)
  
  os_pvalue = surv_pvalue(os_fit)$pval
  print(paste('os_p_value: ', os_pvalue, sep = ''))
  
  # table_5ysr = get_n_year_surv_rate_using_clinical_data(
  #   na.omit(survival_data), n = 5, cluster_label = names(table(survival_data$Cluster))
  # )
  # print(table_5ysr)
}

# 

#### wh immune data #### 
wh_im_cdata = read_excel(
  normalizePath(file.path(
    BASE_DIR, 'input', "WH-SCLC-PD-L1-19-ClinicalData-20210722.xlsx"
  )),
  sheet = 1
)
wh_im_edata = read.csv(
  normalizePath(file.path(
    BASE_DIR, 'input', "WH-SCLC-PD-L1-19-Profiling-v1.csv"
  ))
)
colnames(wh_im_edata) = c(
  'Symbol', apply(data.frame(colnames(wh_im_edata)), 1, function(x){
    strsplit(x, split = '.', fixed = T)[[1]][1]
  })[-1]
)

wh_im_edata.prot_count = apply(wh_im_edata[,-1], 2, function(x){length(which(x>0))})
print(quantile(wh_im_edata.prot_count, seq(0, 1, 0.1)))

wh_im_edata1 = wh_im_edata[-na.omit(match(
  redundant_genes, wh_im_edata$Symbol
)),]
wh_im_edata1_topN = fot_df_correction(get_topx_df_1(wh_im_edata1, topN))

# topN
wh_im_edata2 = data.frame(
  Symbol = wh_im_edata1_topN$Symbol,
  as.matrix(wh_im_edata1_topN[,-1])
)

# wh im prediction
wh_im_edata_model = wh_im_edata2
wh_im_pred_mat = as.matrix(wh_im_edata_model[(match(features, wh_im_edata_model$Symbol)),-1])
print(length(which(is.na(match(features, wh_im_edata_model$Symbol)))))
wh_im_pred_mat[is.na(wh_im_pred_mat)] = 0
rownames(wh_im_pred_mat) = features
wh_im_pred_mat = t(log1p(wh_im_pred_mat))


wh_im_pred_crs <- predict(model_1, newdata = wh_im_pred_mat)
print(table(Pred = wh_im_pred_crs, IM = wh_im_cdata$RECIST_Best))

# KM-Plot: Global
library("RColorBrewer")
col.palette = brewer.pal(n = 12, name = "Paired")
if(T){
  # *** KM curve for OS *** #
  survival_data = data.frame(
    # OS_time = survival_data_v$OS_time,
    OS_time = wh_im_cdata$OS_time,
    OS_status = wh_im_cdata$OS_status,
    Cluster = wh_im_pred_crs,
    TNM_stage = 0
  )
  survival_data = survival_data# [survival_data$Cluster!='C1',]
  os_fit <- survfit(
    Surv(as.numeric(survival_data$OS_time), as.numeric(survival_data$OS_status)) ~ Cluster, data = survival_data
  )
  main_title = paste('Survival curve | N=', nrow(survival_data), ' | WH-IM', sep = '')
  cluster_num = 3 
  os_fit_ggsurvplot <- ggsurvplot(
    os_fit, data = survival_data,
    title = main_title,
    surv.scale = 'percent',
    xlab = 'Months',
    legend = 'right', 
    legend.title = paste("Prediction", '', sep = ''),
    legend.labs = c(
      # 'NoChemo',
      paste('Pred', sort(unique(survival_data$Cluster)), sep = '')
    ),
    pval = TRUE,
    # palette = c('grey', col.palette, '#000000')[-1],
    palette = c('#157E3A', '#11499D', '#E62017'),
    # risk.table.col = "Group",
    risk.table = TRUE, 
    risk.table.y.text.col = TRUE,
    surv.median.line = 'hv',
    # break.time.by = 12,
    # xlim = c(0,12)
  )
  
  print(os_fit_ggsurvplot)
  
  os_pvalue = surv_pvalue(os_fit)$pval
  print(paste('os_p_value: ', os_pvalue, sep = ''))
  
  # table_5ysr = get_n_year_surv_rate_using_clinical_data(
  #   na.omit(survival_data), n = 5, cluster_label = names(table(survival_data$Cluster))
  # )
  # print(table_5ysr)
}


# wh_cdata
sx_cdata$PredCluster = sx_pred_crs
hn_cdata$PredCluster = hn_pred_crs

write.xlsx2(
  wh_cdata, 'SCLC_Clinical_Data_With_Subtypes_Label_20211009.xlsx', sheetName = 'WH', append = T
)

write.xlsx2(
  sx_cdata, 'SCLC_Clinical_Data_With_Subtypes_Label_20211009.xlsx', sheetName = 'SX', append = T
)

write.xlsx2(
  hn_cdata, 'SCLC_Clinical_Data_With_Subtypes_Label_20211009.xlsx', sheetName = 'HN', append = T
)




if(T){
  # *** KM curve for OS *** #
  survival_data = data.frame(
    # OS_time = survival_data_v$OS_time,
    OS_time = wh_cdata$OS_time,
    OS_status = wh_cdata$OS_status,
    Cluster = wh_cdata$Cluster,
    TNM_stage = wh_cdata$TNM_stage
  )
  survival_data = survival_data
  os_fit <- survfit(
    Surv(as.numeric(survival_data$OS_time), as.numeric(survival_data$OS_status)) ~ Cluster, data = survival_data
  )
  main_title = paste('Survival curve (D-set)',  sep = '')
  cluster_num = 3 
  os_fit_ggsurvplot <- ggsurvplot(
    os_fit, data = survival_data,
    title = main_title,
    surv.scale = 'percent',
    xlab = 'Months',
    legend = 'right', 
    legend.title = paste("Subtypes", '', sep = ''),
    legend.labs = c(
      'S-I', 'S-II', 'S-III'
    ),
    pval = TRUE,
    # palette = c('grey', col.palette, '#000000')[-1],
    palette = c('#157E3A', '#11499D', '#E62017'),
    # risk.table.col = "Group",
    risk.table = TRUE, 
    risk.table.y.text.col = TRUE,
    surv.median.line = 'hv',
    # break.time.by = 12,
    # xlim = c(0,12)
  )
  
  print(os_fit_ggsurvplot)
  
  os_pvalue = surv_pvalue(os_fit)$pval
  print(paste('os_p_value: ', os_pvalue, sep = ''))
  
 
}


if(T){
  # *** KM curve for OS *** #
  survival_data = data.frame(
    # OS_time = survival_data_v$OS_time,
    OS_time = hn_cdata$OS_time_day/30,
    OS_status = hn_cdata$OS_status,
    Cluster = hn_pred_crs,
    TNM_stage = 0
  )
  survival_data = survival_data# [survival_data$Cluster!='C1',]
  os_fit <- survfit(
    Surv(as.numeric(survival_data$OS_time), as.numeric(survival_data$OS_status)) ~ Cluster, data = survival_data
  )
  main_title = paste('Survival curve (V-set)', sep = '')
  cluster_num = 3 
  os_fit_ggsurvplot <- ggsurvplot(
    os_fit, data = survival_data,
    title = main_title,
    surv.scale = 'percent',
    xlab = 'Months',
    legend = 'right', 
    legend.title = paste("Subtypes", '', sep = ''),
    legend.labs = c(
      'S-I', 'S-II', 'S-III'
    ),
    pval = TRUE,
    # palette = c('grey', col.palette, '#000000')[-1],
    palette = c('#157E3A', '#11499D', '#E62017'),
    # risk.table.col = "Group",
    risk.table = TRUE, 
    risk.table.y.text.col = TRUE,
    surv.median.line = 'hv',
    break.time.by = 12,
    # xlim = c(0,12)
  )
  
  print(os_fit_ggsurvplot)
  
  os_pvalue = surv_pvalue(os_fit)$pval
  print(paste('os_p_value: ', os_pvalue, sep = ''))
  
  # table_5ysr = get_n_year_surv_rate_using_clinical_data(
  #   na.omit(survival_data), n = 5, cluster_label = names(table(survival_data$Cluster))
  # )
  # print(table_5ysr)
}


if(T){
  #### Cox ####
  if(T){
    
    clinical_data.temp1 = wh_cdata
    
    coxph_data = data.frame(
      ExpCode = clinical_data.temp1$Proteome_fno_T,
      OS_time = as.numeric(clinical_data.temp1$OS_time),
      OS_status = as.numeric(clinical_data.temp1$OS_status),
      Cluster = (clinical_data.temp1$Cluster),
      Age = ifelse(clinical_data.temp1$Age>=60, '>=60', '<60'),
      Gender = as.factor(clinical_data.temp1$Gender),
      Radiotherapy = ifelse(clinical_data.temp1$Radiotherapy==0, 'No', 'Yes'),
      Chemotherapy = ifelse(clinical_data.temp1$Chemotherapy_status==0, 'No', 'Yes'),
      LN_metastasis = ifelse(clinical_data.temp1$LN_metastasis==0, 'No', 'Yes'),
      Smoke_history = ifelse(clinical_data.temp1$Smoke_history==0, 'No', 'Yes'),
      CVD = ifelse(clinical_data.temp1$CVD==0, 'No', 'Yes'),
      TNM_stage = as.factor(clinical_data.temp1$TNM_stage),
      VALG_stage = factor(clinical_data.temp1$VALG_stage, levels = c('ED-SCLC', 'LD-SCLC'))
    )
    c1.index = which(coxph_data$Cluster == 'C1')
    c2.index = which(coxph_data$Cluster == 'C2')
    c3.index = which(coxph_data$Cluster == 'C3')
    coxph_data$Cluster[c1.index] = 'S-I'
    coxph_data$Cluster[c2.index] = 'S-II'
    coxph_data$Cluster[c3.index] = 'S-III'
    
    
    cox_model = coxph(
      Surv(OS_time, OS_status) ~ Cluster 
      + VALG_stage
      + TNM_stage
      + Chemotherapy
      + Radiotherapy
      + Smoke_history
      + LN_metastasis
      + CVD 
      + Gender
      + Age
      
      
      
      , 
      data = coxph_data
    )
    cox_model.summary = summary(cox_model)
    # plot(ggforest(cox_model, data = coxph_data, main = 'Multivariate Cox analysis (TNM=II/III)'))
    plot(ggforest(cox_model, data = coxph_data, main = 'Multivariate Cox analysis'))
    coxpvalue = c(coxpvalue, cox_model.summary$coefficients[1,5])
  }
}

write.csv(wh_edata2, 'WH-SCLC-75-Profiling-Top500.csv', row.names = F)
write.csv(hn_edata2, 'HN-SCLC-52-Profiling-Top500.csv', row.names = F)
