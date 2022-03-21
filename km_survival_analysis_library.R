# survival_data
# colors = c(
#   4, #1
#   3, #2
#   "#00FF40FF", #3
#   "#00FFFFFF", #4
#   6, #5
#   2, #6
#   'purple', #7
#   1  #8
# )
# dataset_type = 'Global'

km_survival_analysis <- function(survival_data, dataset_type = 'Global', colors){
  # *** KM curve for OS *** #
  cat("                                                     \n")
  cat("*****************************************************\n")
  cat(paste("KM curve for OS", "  "))
  os_fit <- survfit(
    Surv(as.numeric(survival_data$OS_time), as.numeric(survival_data$OS_status)) ~ Cluster, data = survival_data
  )
  main_title = paste('Survival curve (OS) | ', dataset_type, ' (N=', nrow(survival_data), ')',sep = '')
  os_fit_ggsurvplot <- ggsurvplot(
    os_fit, data = survival_data,
    title = main_title,
    surv.scale = 'percent',
    xlab = 'Years',
    legend = 'right', 
    legend.title = "Subtypes",
    legend.labs = paste('GCP', sort(unique(survival_data$Cluster)), sep = ''),
    pval = TRUE,
    palette = colors,
    # risk.table.col = "Group",
    risk.table = TRUE, 
    risk.table.y.text.col = TRUE,
    surv.median.line = 'hv',
    break.time.by = 2,
    xlim = c(0,13)
  )
  print(os_fit_ggsurvplot)
  os_pvalue = surv_pvalue(os_fit)$pval
  print(paste(dataset_type, '_os_p_value: ', os_pvalue, sep = ''))
  
  
  
  # *** KM curve for chemotherapy *** #
  cat("                                                     \n")
  cat("*****************************************************\n")
  cat(paste("KM curve for chemotherapy", "  \n"))
  OS_time_cutoff = 0
  for(st in 1:length(unique(survival_data$Cluster))){
    na_index = which(survival_data$Chemotherapy_status == 'NA' | survival_data$Chemotherapy_status == 'unknown')
    if(length(na_index)>0){
      survival_data$Chemotherapy_status[na_index] = NA
    }
    
    if(T){
      cluster_index = which(
        survival_data$Cluster == st 
        & survival_data$TNM_stage < 4 
        & (survival_data$TNM_stage > 1 
           & (survival_data$OS_time > OS_time_cutoff)) # 
      )
    }
    
    if(length(cluster_index)>2){
      # cluster_index = which(survival_data$TNM_stage>1 & survival_data$Hospital == '301')
      cluster_data = survival_data[cluster_index,]
      os_fit <- survfit(
        # Surv(as.numeric(round(OS_time,0)), as.numeric(OS_status)) ~ as.numeric(Group) + as.numeric(Age) + as.numeric(TNM_stage) + as.numeric(Gender), data = survival_data
        Surv(as.numeric(OS_time), as.numeric(OS_status)) ~ Chemotherapy_status, data = cluster_data
      )
      
      # Change the legend title and labels
      os_fit_ggsurvplot <- ggsurvplot(
        os_fit, data = cluster_data,
        title = paste('Survival curve (OS) | ', 'Cluster = GCP', st, ' (', dataset_type, ')', sep = ''),
        surv.scale = 'percent',
        xlab = 'Years',
        legend = 'right', 
        legend.title = "Chemotherapy",
        legend.labs = c('No chemo', 'Chemo'),
        palette = c('Blue', 'red'),
        pval = TRUE,
        surv.median.line = 'hv',
        risk.table = TRUE, 
        risk.table.y.text.col = TRUE,
        break.time.by = 2,
        xlim = c(0,13)
      )
      print(os_fit_ggsurvplot)
    }
  }
  
}




km_survival_analysis_ACT <- function(
  survival_data, dataset_type = 'Global', colors, act_name
){
  # *** KM curve for OS *** #
  cat("                                                     \n")
  cat("*****************************************************\n")
  cat(paste("KM curve for OS", "  "))
  act_index = match(act_name, colnames(survival_data))
  survival_data$ACT = survival_data[,act_index]
  if(act_name == 'ACT1_LS'){
    legend_labs = c('No chemo', 'FU+Pt', 'Others')
    palettes = c('Blue', '#FF00FF', 'green')
  }
  if(act_name == 'ACT1_LS_JQ'){
    legend_labs = c('No chemo', 'FOLFOX', 'XELOX', 'SOX', 'Others')
    palettes = c('Blue', '#FF00FF', 'green', 'orange', 'purple')
  }
  
  main_title = paste('Survival curve (OS) | ', dataset_type, ' (N=', nrow(survival_data), ')',sep = '')
  
  
  
  # *** KM curve for ACT *** #
  cat("                                                     \n")
  cat("*****************************************************\n")
  cat(paste("KM curve for ACT", "  "))
  OS_time_cutoff = 0
  for(st in 1:length(unique(survival_data$Cluster))){
    na_index = which(survival_data$Chemotherapy_status == 'NA' | survival_data$Chemotherapy_status == 'unknown')
    if(length(na_index)>0){
      survival_data$Chemotherapy_status[na_index] = NA
    }
    
    if(T){
      cluster_index = which(
        survival_data$Cluster == st 
        & survival_data$TNM_stage < 4 
        & (survival_data$TNM_stage > 1 
           & (survival_data$OS_time > OS_time_cutoff)) # 
      )
    }
    
    if(length(cluster_index)>2){
      # cluster_index = which(survival_data$TNM_stage>1 & survival_data$Hospital == '301')
      cluster_data = survival_data[cluster_index,]
      os_fit <- survfit(
        # Surv(as.numeric(round(OS_time,0)), as.numeric(OS_status)) ~ as.numeric(Group) + as.numeric(Age) + as.numeric(TNM_stage) + as.numeric(Gender), data = survival_data
        Surv(as.numeric(OS_time), as.numeric(OS_status)) ~ ACT, data = cluster_data
      )
      
      
      # Change the legend title and labels
      os_fit_ggsurvplot <- ggsurvplot(
        os_fit, data = cluster_data,
        title = paste('Survival curve (OS) | ', 'Cluster = GCP', st, ' (', dataset_type, ')', sep = ''),
        surv.scale = 'percent',
        xlab = 'Years',
        legend = 'right', 
        legend.title = "Chemotherapy",
        legend.labs = legend_labs,
        palette = palettes,
        pval = TRUE,
        surv.median.line = 'hv',
        risk.table = TRUE, 
        risk.table.y.text.col = TRUE,
        break.time.by = 2,
        xlim = c(0,13)
      )
      print(os_fit_ggsurvplot)
    }
  }
  
}

cox_survival_analysis <- function(dim_summary_data){
  if(T){
    library('readxl')
    library('xlsx')
    library('stringr')
    library("ctv")
    library("kernlab")
    library("survival")
    library("survminer")
    library("pracma")
    
    # dim_summary_data = survival_data_global
    dim_summary_data$Cluster = factor(dim_summary_data$Cluster, levels = seq(1,length(unique(dim_summary_data$Cluster))))
    dim_summary_data$MCluster = rep(1, nrow(dim_summary_data))
    if(T){
      
      # Age
      dim_summary_data$Age_label = rep(0, nrow(dim_summary_data))
      dim_summary_data$Age_label[dim_summary_data$Age>=60] = 1
      
      # TNM_stage
      # dim_summary_data$TNM_stage = dim_summary_data$pStageSimple
      dim_summary_data$TNM_stage_new = dim_summary_data$TNM_stage
      # dim_summary_data$TNM_stage_new[which(dim_summary_data$TNM_stage_new<=2)] = 1
      # dim_summary_data$TNM_stage_new[which(dim_summary_data$TNM_stage_new>2)] = 2
      
      
      dim_summary_data$Chemotherapy_status[which(dim_summary_data$Chemotherapy_status == 'unknown')] = NA
    }
    
    
    if(T){
      # classification
      coxph_data = dim_summary_data
      coxph_ggforest(coxph_data, class_flag = 'Cluster', data_flag = 1, main = 'Hazzard ratio of Cluster6')
    }
    
    if(T){
      # p-trend
      coxph_ggforest(coxph_data = dim_summary_data, class_flag = 'Cluster', data_flag = 2, main = 'Hazzard ratio of Cluster6 (P-trend)')
    }
    
    chemo_k = seq(1, length(unique((dim_summary_data$Cluster))))
    for (k in chemo_k) {
      if(T){
        dim_summary_data_1_index = which(
          dim_summary_data$Cluster==k
          & dim_summary_data$TNM_stage > 1
          & dim_summary_data$TNM_stage < 4
        )
        dim_summary_data_1 = dim_summary_data[dim_summary_data_1_index,]
        coxph_ggforest(
          coxph_data = dim_summary_data_1, 
          class_flag = 'Chemo', 
          data_flag = 1,
          main = paste('Hazzard ratio of Cluster6 (chemo_k = ', k, ')', sep = '')
        )
        print(k)
      }
    }
  }
}





km_survival_analysis_total <- function(survival_data, dataset_type = 'Global', colors){
  # *** KM curve for OS *** #
  cat("                                                     \n")
  cat("*****************************************************\n")
  cat(paste("KM curve for OS", "  "))
  os_fit <- survfit(
    Surv(as.numeric(survival_data$OS_time), as.numeric(survival_data$OS_status)) ~ Cluster, 
    data = survival_data
  )
  main_title = paste('Survival curve (OS) | ', dataset_type, ' (N=', nrow(survival_data), ')',sep = '')
  os_fit_ggsurvplot <- ggsurvplot(
    os_fit, data = survival_data,
    title = main_title,
    surv.scale = 'percent',
    xlab = 'Years',
    legend = 'right', 
    legend.title = "Subtypes",
    legend.labs = paste('GCP', sort(unique(survival_data$Cluster)), sep = ''),
    pval = TRUE,
    palette = colors,
    # risk.table.col = "Group",
    risk.table = TRUE, 
    risk.table.y.text.col = TRUE,
    surv.median.line = 'hv',
    break.time.by = 2,
    xlim = c(0,13)
  )
  print(os_fit_ggsurvplot)
  os_pvalue = surv_pvalue(os_fit)$pval
  print(paste(dataset_type, '_os_p_value: ', os_pvalue, sep = ''))
  
}


coxph_ggforest_act <- function(coxph_data, class_flag, main){
  library(survival)
  library(survminer)
  # coxph_data = dim_summary_data
  main = paste(main, '\n', class_flag, sep = '')
  coxph_data$OS_time = as.numeric(coxph_data$OS_time)
  coxph_data$OS_status = as.numeric(coxph_data$OS_status)
  coxph_data$Age_label = as.factor(coxph_data$Age_label)
  coxph_data$Gender = as.factor(coxph_data$Gender)
  coxph_data$TNM_stage_new = as.factor(coxph_data$TNM_stage_new)
  coxph_data$Chemotherapy_status = as.factor(coxph_data$Chemotherapy_status)
  coxph_data$Cluster = as.factor(coxph_data$Cluster)
  cluster_num = length(unique(coxph_data$Cluster))
  act_index = which(colnames(coxph_data) == class_flag)
  act_flag = coxph_data[,act_index]
  act_name = rep(NA, length(act_flag))
  if(class_flag == 'ACT1_LS'){
    legend_labs = c('No chemo', 'FU+Pt', 'Others')
    palettes = c('Blue', '#FF00FF', 'green')
  }
  if(class_flag == 'ACT1_LS_JQ'){
    legend_labs = c('No chemo', 'FOLFOX', 'XELOX', 'SOX', 'Others')
    palettes = c('Blue', '#FF00FF', 'green', 'orange', 'purple')
  }
  act_flag_unique = na.omit(unique(act_flag))
  for (act_i in 1:length(act_flag_unique)) {
    act_i_index = which(act_flag == act_flag_unique[act_i])
    act_name[act_i_index] = legend_labs[act_flag_unique[act_i]+1]
  }
  act_name = factor(act_name, levels = legend_labs)
  coxph_data$ACT = act_name
  
  for (ci in 1:cluster_num) {
    ci_coxph_data = coxph_data[which(coxph_data$Cluster==ci),]
    ci_cox_model = coxph(
      # Surv(OS_time, OS_status) ~ Age_label + Gender + TNM_stage_new + Chemotherapy_status + ACT, 
      Surv(OS_time, OS_status) ~ ACT, 
      data = ci_coxph_data
    )
    plot(ggforest(ci_cox_model, data = ci_coxph_data, main = paste(main, '_GCP', ci, sep = '')))
  }
  
}


cox_survival_analysis_for_act <- function(dim_summary_data, class_flag = 'ACT1_LS', main = 'Hazzard ratio of Cluster6'){
  if(T){
    library('readxl')
    library('xlsx')
    library('stringr')
    library("ctv")
    library("kernlab")
    library("survival")
    library("survminer")
    library("pracma")
    
    # dim_summary_data = survival_data_global
    dim_summary_data$Cluster = factor(dim_summary_data$Cluster, levels = seq(1,length(unique(dim_summary_data$Cluster))))
    dim_summary_data$MCluster = rep(1, nrow(dim_summary_data))
    if(T){
      
      # Age
      dim_summary_data$Age_label = rep(0, nrow(dim_summary_data))
      dim_summary_data$Age_label[dim_summary_data$Age>=60] = 1
      
      # TNM_stage
      # dim_summary_data$TNM_stage = dim_summary_data$pStageSimple
      dim_summary_data$TNM_stage_new = dim_summary_data$TNM_stage
      # dim_summary_data$TNM_stage_new[which(dim_summary_data$TNM_stage_new<=2)] = 1
      # dim_summary_data$TNM_stage_new[which(dim_summary_data$TNM_stage_new>2)] = 2
      
      
      dim_summary_data$Chemotherapy_status[which(dim_summary_data$Chemotherapy_status == 'unknown')] = NA
    }
    
    if(T){
      # classification
      coxph_data = dim_summary_data[
        which(dim_summary_data$TNM_stage > 1 & dim_summary_data$TNM_stage < 4),
        ]
      coxph_ggforest_act(coxph_data, class_flag = class_flag, main = main)
    }
    
  }
}



