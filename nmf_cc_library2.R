# *** Get 5-year survival rate *** #
# get_n_year_surv_rate <- function(surv_model, surv_df, n = 5, cluster_num){
#   res_sum <- surv_summary(surv_model, data = surv_df)
#   res_sum_subset = res_sum[res_sum$n.event>0,]
#   res_sum_subset[,c(1,5:8)] = round(res_sum_subset[,c(1,5:8)], 3)
#   n_year_surv_rate = NULL
#   strata_prefix = strsplit(names(surv_model$strata[1]), split = '=')[[1]][1]
#   for (i in 1:cluster_num) {
#     strata_i = paste(strata_prefix, '=', i, sep = '')
#     strata_i_index = which(res_sum_subset$strata == strata_i)
#     strata_i_year = res_sum_subset$time[strata_i_index]
#     strata_i_n_year_surv_rate = max(strata_i_year[strata_i_year < n])
#     n_year_surv_rate = c(n_year_surv_rate, strata_i_n_year_surv_rate)
#   }
#   n_year_surv_rate_index = match(n_year_surv_rate, res_sum_subset$time)
#   n_year_surv_rate_df = res_sum_subset[n_year_surv_rate_index, c(1,5,7,8)]
#   n_year_surv_rate_df$group = seq(1, cluster_num)
#   order_index_desc = order(n_year_surv_rate_df$surv, decreasing = T)
#   n_year_surv_rate_df = n_year_surv_rate_df[order_index_desc,]
#   n_year_surv_rate_df$cluster = as.factor(seq(1,cluster_num))
#   return(n_year_surv_rate_df)
# }


# *** Get 5-year survival rate *** #
get_n_year_surv_rate <- function(surv_model, surv_df, n = 5){
  res_sum <- surv_summary(surv_model, data = surv_df)
  res_sum_subset = res_sum[res_sum$n.event>0,]
  res_sum_subset[,c(1,5:8)] = round(res_sum_subset[,c(1,5:8)], 3)
  n_year = NULL
  n_year_surv_rate = NULL
  n_upper = NULL
  n_lower = NULL
  strata_prefix = strsplit(names(surv_model$strata[1]), split = '=')[[1]][1]
  cluster_num =  str_split_fixed(names(surv_model$strata), '=', 2)[, 2]
  for (i in cluster_num) {
    strata_i = paste(strata_prefix, '=', i, sep = '')
    strata_i_index = which(res_sum_subset$strata == strata_i)
    strata_i_year = res_sum_subset$time[strata_i_index]
    strata_i_year_surv_rate = res_sum_subset$surv[strata_i_index]
    strata_i_upper = res_sum_subset$upper[strata_i_index]
    strata_i_lower = res_sum_subset$lower[strata_i_index]
    strata_i_n_year = max(strata_i_year[which(strata_i_year<n)])
    strata_i_n_year_index = which(strata_i_year == strata_i_n_year)
    strata_i_n_year_surv_rate = strata_i_year_surv_rate[strata_i_n_year_index]
    strata_i_n_upper = strata_i_upper[strata_i_n_year_index]
    strata_i_n_lower = strata_i_lower[strata_i_n_year_index]
    n_year = c(n_year, strata_i_n_year)
    n_year_surv_rate = c(n_year_surv_rate, strata_i_n_year_surv_rate)
    n_upper = c(n_upper, strata_i_n_upper)
    n_lower = c(n_lower, strata_i_n_lower)
  }
  n_year_surv_rate_df = data.frame(
    time = n_year,
    surv = n_year_surv_rate,
    upper = n_upper,
    lower = n_lower
  )
  n_year_surv_rate_df$group = cluster_num
  order_index_desc = order(n_year_surv_rate_df$surv, decreasing = T)
  n_year_surv_rate_df = n_year_surv_rate_df[order_index_desc,]
  n_year_surv_rate_df$cluster = 1:length(cluster_num)
  # n_year_surv_rate_df$count = as.vector(table(n_year_surv_rate_df$cluster))
  return(n_year_surv_rate_df)
}



# *** Input *** #
# main_title = paste("Global set (N=681)", sep = '')
# seed = 129
# nmf_data = data_all # feature matrix
# cluster_num = 7
# nrun=50
# time = os_time
# status = os_status
# ylab = "Probability of Overall Survival"
# cluster_colors = NULL
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

get_nmf_result_list <- function(
  main_title,
  seed = 129,
  nmf_data,
  cluster_num = 6,
  nrun = 50,
  time = os_time,
  status = os_status,
  n_year = 5,
  ylab = "Probability of Overall Survival",
  cluster_colors = NULL
){
  require("survival")
  require("survminer")
  require('CancerSubtypes')
  
  
  # *** Execute NMF *** #
  cat("                                                     \n")
  cat("*****************************************************\n")
  cat(paste("Execute NMF: ", "Cluster=", cluster_num, "  "))
  set.seed(seed)
  system.time({
    nmf_result = ExecuteCNMF(nmf_data, clusterNum = cluster_num, nrun = nrun)
  })
  nmf_group = nmf_result$group
  nmf_distanceMatrix = nmf_result$distanceMatrix
  
  
  # *** Survival analysis - group *** #
  surv_df = data.frame(
    time = time, 
    status = status, 
    group = nmf_group
  )
  surv_model <- survfit(
    Surv(as.numeric(time), as.numeric(status)) ~ as.factor(group), data = surv_df
  )
  
  # *** 5-year survival rate *** #
  n_year_surv_rate_df = get_n_year_surv_rate(
    surv_model, surv_df, n = n_year #, cluster_num
  )
  rownames(n_year_surv_rate_df) = n_year_surv_rate_df$cluster
  cat("                                                     \n")
  cat("*****************************************************\n")
  cat(paste(n_year, "-year survival rate table", "  ", sep = ''))
  print(n_year_surv_rate_df)
  
  
  # *** Adjust cluster *** #
  surv_df$cluster = rep(NA, nrow(surv_df))
  for (i in 1:cluster_num) {
    group_i = n_year_surv_rate_df$group[i]
    cluster_i = as.numeric(as.vector(n_year_surv_rate_df$cluster[i]))
    index_i = which(surv_df$group == group_i)
    surv_df$cluster[index_i] = cluster_i
  }
  cat("                                                     \n")
  cat("*****************************************************\n")
  cat(paste(n_year, "-Cluster table", "  ", sep = ''))
  print(table(surv_df$cluster))
  
  
  # *** Survival analysis - cluster *** #
  surv_model <- survfit(
    Surv(as.numeric(time), as.numeric(status)) ~ as.factor(cluster), data = surv_df
  )
  sm_pvalue = surv_pvalue(surv_model, data = surv_df)$pval
  
  
  # *** Layout *** #
  cat("                                                     \n")
  cat("*****************************************************\n")
  cat(paste("Layout", "  ", sep = ''))
  if(cluster_num > 1) {
    cat("                                                     \n")
    cat("*****************************************************\n")
    cat(paste(main_title, "Cluster=", cluster_num, "  "))
    p_value = signif(sm_pvalue, 3)
  }else {
    cat("There is only one cluster in the group")
    p_value = 1
  }
  if(!is.null(nmf_distanceMatrix[1, 1])) {
    layout(
      matrix(c(1, 2, 3, 3), 2, 2, byrow = FALSE), 
      widths = c(2.2, 2), 
      heights = c(2, 2)
    )
  }
  
  # *** KM curve *** #
  cat("                                                     \n")
  cat("*****************************************************\n")
  cat(paste("KM curve", "  ", sep = ''))
  title = paste(main_title, " | Cluster =", cluster_num)
  if(is.null(cluster_colors)){
    if(cluster_num <= 7){
      cluster_colors = seq(1, cluster_num)
    }else{
      cluster_colors = c(seq(1, cluster_num), seq(1, cluster_num-7))
    }
    
  }
  # cluster_colors = seq(1,6)
  plot(
    surv_model, 
    lty = 1, 
    col = cluster_colors, 
    lwd = 2,
    xscale = 1,
    xlab = "Years", 
    ylab = ylab, 
    main = title, 
    font.main = 4, 
    cex.main = 0.9
  )
  legend(
    'topright',
    paste("Subtpye", 1:cluster_num), 
    lty = 1, 
    lwd = 3, 
    cex = 0.8, 
    text.font = 4, 
    text.col = cluster_colors, 
    bty = "n", 
    col = cluster_colors, 
    seg.len = 0.3
  )
  digit = ceiling(-log10(p_value) + 2)
  text(
    x = par("usr")[2] * 0.25, 
    y = par("usr")[4] * 0.1, 
    paste("Log-rank P-value = ", p_value, sep = ''), 
    col = "red", 
    font = 3, 
    cex = 1
  )
  
  # *** Heat map *** #
  cat("                                                     \n")
  cat("*****************************************************\n")
  cat(paste("Heat map", "  ", sep = ''))
  if(!is.null(nmf_distanceMatrix[1, 1])){
    nmf_cluster = surv_df$cluster
    if (class(nmf_distanceMatrix) == "Similarity") {
      si = silhouette_SimilarityMatrix(nmf_cluster, nmf_distanceMatrix)
    }else {
      si = silhouette(nmf_cluster, nmf_distanceMatrix)
    }
    attr(nmf_distanceMatrix, "class") = NULL
    ind = order(nmf_cluster, -si[, "sil_width"])
    num = length(unique(nmf_cluster))
    annotation = data.frame(Cluster = as.factor(nmf_cluster))
    Var1 = cluster_colors
    names(Var1) = sort(unique(nmf_cluster))
    ann_colors = list(Cluster = Var1)
    consensusmap(
      nmf_distanceMatrix, Rowv = ind, Colv = ind, 
      main = "Clustering display", annCol = annotation, 
      annColors = ann_colors, 
      # labRow = "Sample", labCol = "Sample", 
      scale = "none"
    )
  }else{
    print(paste("nmf_distanceMatrix is NULL", "  ", sep = ''))
  }
  
  # *** silhouette *** #
  cat("                                                     \n")
  cat("*****************************************************\n")
  cat(paste("Silhouette", "  ", sep = ''))

  si_col = rep(
    cluster_colors[1:cluster_num],
    as.vector(table(nmf_cluster))
  )
  plot(si, col = si_col)
  par(mfrow = c(1, 1))
  
  # *** Return NMF result *** #
  cat("                                                     \n")
  cat("*****************************************************\n")
  # cat(paste("Cluster labels", "  ", sep = ''))
  # cat(paste("Log-rank Pvalue", "  ", sep = ''))
  # cat(paste("N-year survival rate", "  ", sep = ''))
  nmf_result_list = list(
    cluster = nmf_cluster,
    Pvalue = sm_pvalue,
    n_year_surv_rate = n_year_surv_rate_df
  )
  
  return(nmf_result_list)
  
}










get_nmf_result_list_2 <- function(
  main_title,
  seed = 129,
  nmf_data,
  cluster_num = 6,
  nrun = 50,
  time = os_time,
  status = os_status,
  nmf_clinical_data,
  n_year = 5,
  ylab = "Probability of Overall Survival",
  cluster_colors = NULL
){
  require("survival")
  require("survminer")
  require('CancerSubtypes')
  
  
  # *** Execute NMF *** #
  cat("                                                     \n")
  cat("*****************************************************\n")
  cat(paste("Execute NMF: ", "Cluster=", cluster_num, "  "))
  set.seed(seed)
  system.time({
    nmf_result = ExecuteCNMF(nmf_data, clusterNum = cluster_num, nrun = nrun)
  })
  nmf_group = nmf_result$group
  nmf_distanceMatrix = nmf_result$distanceMatrix
  # nmf_clinical_data = clinical_data_discovery
  nmf_clinical_data$Cluster = nmf_group
  
  # *** Survival analysis - cox - chemo*** #
  coxph_chemo_index_df = get_coxph_index_df(nmf_clinical_data, cluster_num)
  cat("                                                     \n")
  cat("*****************************************************\n")
  cat(paste('Cox analysis', "-chemotheray", " \n ", sep = ''))
  print(coxph_chemo_index_df)
  
  # *** Survival analysis *** #
  surv_df = data.frame(
    time = time, 
    status = status, 
    group = nmf_group
  )
  
  # *** Adjust cluster *** #
  surv_df$cluster = rep(NA, nrow(surv_df))
  for (i in 1:cluster_num) {
    group_i = coxph_chemo_index_df$group[i]
    cluster_i = as.numeric(as.vector(coxph_chemo_index_df$cluster[i]))
    index_i = which(surv_df$group == group_i)
    surv_df$cluster[index_i] = cluster_i
  }
  cat("                                                     \n")
  cat("*****************************************************\n")
  cat(paste('Cox analysis', "-Cluster table", " \n ", sep = ''))
  print(table(surv_df$cluster))
  
  # *** Survival analysis - KM *** #
  surv_model <- survfit(
    Surv(as.numeric(time), as.numeric(status)) ~ as.factor(cluster), data = surv_df
  )
  
  # *** 5-year survival rate *** #
  n_year_surv_rate_df = get_n_year_surv_rate(
    surv_model, surv_df, n = n_year
  )
  rownames(n_year_surv_rate_df) = n_year_surv_rate_df$cluster
  cat("                                                     \n")
  cat("*****************************************************\n")
  cat(paste(n_year, "-year survival rate table", " \n  ", sep = ''))
  print(n_year_surv_rate_df)
  
  
 
  
  
  
  
  
  # *** Survival analysis - cluster *** #
  surv_model <- survfit(
    Surv(as.numeric(time), as.numeric(status)) ~ as.factor(cluster), data = surv_df
  )
  sm_pvalue = surv_pvalue(surv_model, data = surv_df)$pval
  
  
  # *** Layout *** #
  cat("                                                     \n")
  cat("*****************************************************\n")
  cat(paste("Layout", "  ", sep = ''))
  if(cluster_num > 1) {
    cat("                                                     \n")
    cat("*****************************************************\n")
    cat(paste(main_title, "Cluster=", cluster_num, "  "))
    p_value = signif(sm_pvalue, 3)
  }else {
    cat("There is only one cluster in the group")
    p_value = 1
  }
  if(!is.null(nmf_distanceMatrix[1, 1])) {
    layout(
      matrix(c(1, 2, 3, 3), 2, 2, byrow = FALSE), 
      widths = c(2.2, 2), 
      heights = c(2, 2)
    )
  }
  
  # *** KM curve *** #
  cat("                                                     \n")
  cat("*****************************************************\n")
  cat(paste("KM curve", "  ", sep = ''))
  title = paste(main_title, " | Cluster =", cluster_num)
  if(is.null(cluster_colors)){
    if(cluster_num <= 7){
      cluster_colors = seq(1, cluster_num)
    }else{
      cluster_colors = c(seq(1, cluster_num), seq(1, cluster_num-7))
    }
    
  }
  # cluster_colors = seq(1,6)
  plot(
    surv_model, 
    lty = 1, 
    col = cluster_colors, 
    lwd = 2,
    xscale = 1,
    xlab = "Years", 
    ylab = ylab, 
    main = title, 
    font.main = 4, 
    cex.main = 0.9
  )
  legend(
    'topright',
    paste("Subtpye", 1:cluster_num), 
    lty = 1, 
    lwd = 3, 
    cex = 0.8, 
    text.font = 4, 
    text.col = cluster_colors, 
    bty = "n", 
    col = cluster_colors, 
    seg.len = 0.3
  )
  digit = ceiling(-log10(p_value) + 2)
  text(
    x = par("usr")[2] * 0.25, 
    y = par("usr")[4] * 0.1, 
    paste("Log-rank P-value = ", p_value, sep = ''), 
    col = "red", 
    font = 3, 
    cex = 1
  )
  
  # *** Heat map *** #
  cat("                                                     \n")
  cat("*****************************************************\n")
  cat(paste("Heat map", "  ", sep = ''))
  if(!is.null(nmf_distanceMatrix[1, 1])){
    nmf_cluster = surv_df$cluster
    if (class(nmf_distanceMatrix) == "Similarity") {
      si = silhouette_SimilarityMatrix(nmf_cluster, nmf_distanceMatrix)
    }else {
      si = silhouette(nmf_cluster, nmf_distanceMatrix)
    }
    attr(nmf_distanceMatrix, "class") = NULL
    ind = order(nmf_cluster, -si[, "sil_width"])
    num = length(unique(nmf_cluster))
    annotation = data.frame(Cluster = as.factor(nmf_cluster))
    Var1 = cluster_colors
    names(Var1) = sort(unique(nmf_cluster))
    ann_colors = list(Cluster = Var1)
    consensusmap(
      nmf_distanceMatrix, Rowv = ind, Colv = ind, 
      main = "Clustering display", annCol = annotation, 
      annColors = ann_colors, 
      # labRow = "Sample", labCol = "Sample", 
      scale = "none"
    )
  }else{
    print(paste("nmf_distanceMatrix is NULL", "  ", sep = ''))
  }
  
  # *** silhouette *** #
  cat("                                                     \n")
  cat("*****************************************************\n")
  cat(paste("Silhouette", "  ", sep = ''))
  
  si_col = rep(
    cluster_colors[1:cluster_num],
    as.vector(table(nmf_cluster))
  )
  plot(si, col = si_col)
  par(mfrow = c(1, 1))
  
  # *** Return NMF result *** #
  cat("                                                     \n")
  cat("*****************************************************\n")
  # cat(paste("Cluster labels", "  ", sep = ''))
  # cat(paste("Log-rank Pvalue", "  ", sep = ''))
  # cat(paste("N-year survival rate", "  ", sep = ''))
  nmf_result_list = list(
    cluster = nmf_cluster,
    Pvalue = sm_pvalue,
    n_year_surv_rate = n_year_surv_rate_df,
    coxph_chemo_index_df = coxph_chemo_index_df
  )
  
  return(nmf_result_list)
  
}




 
get_coxph_index_df <- function(dim_summary_data, cluster_num){
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
      dim_summary_data$TNM_stage_new = dim_summary_data$TNM_stage
      
      # Chemotherapy_status
      dim_summary_data$Chemotherapy_status[which(dim_summary_data$Chemotherapy_status == 'unknown')] = NA
    }
    coxph_index_df = NULL
    for (k in 1:cluster_num) {
      if(T){
        dim_summary_data_1_index = which(
          dim_summary_data$Cluster==k
          & dim_summary_data$TNM_stage > 1
          & dim_summary_data$TNM_stage < 4
        )
        dim_summary_data_1 = dim_summary_data[dim_summary_data_1_index,]
        coxph_data = dim_summary_data_1
        coxph_index = get_coxph_index(coxph_data)
        coxph_index_df  = rbind(coxph_index_df, coxph_index)
        # print(k)
      }
    }
    colnames(coxph_index_df) = c('HR', 'Pvalue')
    coxph_index_df = as.data.frame(coxph_index_df)
    
    coxph_index_df$group = seq(1, cluster_num)
    order_index_asc = order(coxph_index_df$Pvalue, decreasing = F)
    coxph_index_df = coxph_index_df[order_index_asc,]
    coxph_index_df$cluster = as.factor(seq(1, cluster_num))
    rownames(coxph_index_df) = seq(1, cluster_num)
    
    return(coxph_index_df)
  }
}


get_coxph_index <- function(coxph_data){
  library(survival)
  library(survminer)
  # coxph_data = dim_summary_data
  coxph_data$OS_time = as.numeric(coxph_data$OS_time)
  coxph_data$OS_status = as.numeric(coxph_data$OS_status)
  coxph_data$Age_label = as.factor(coxph_data$Age_label)
  coxph_data$Gender = as.factor(coxph_data$Gender)
  coxph_data$TNM_stage_new = as.factor(coxph_data$TNM_stage_new)
  coxph_data$Chemotherapy_status = as.factor(coxph_data$Chemotherapy_status)
  cox_model = coxph(
    Surv(OS_time, OS_status) ~ Age_label + Gender + TNM_stage_new + Chemotherapy_status, data = coxph_data
  )
  cox_model_summary = summary(cox_model)
  cox_model_summary_coefficients = cox_model_summary$coefficients
  hr = cox_model_summary_coefficients['Chemotherapy_status1', 'exp(coef)']
  pvalue = cox_model_summary_coefficients['Chemotherapy_status1', 'Pr(>|z|)']
  cox_reuslt = c(hr, pvalue)
  return(cox_reuslt)
}

# *** Get 5-year survival rate *** #
get_n_year_surv_rate_cox <- function(surv_model, surv_df, n = 5, cluster_num){
  res_sum <- surv_summary(surv_model, data = surv_df)
  res_sum_subset = res_sum[res_sum$n.event>0,]
  res_sum_subset[,c(1,5:8)] = round(res_sum_subset[,c(1,5:8)], 3)
  n_year = NULL
  n_year_surv_rate = NULL
  n_upper = NULL
  n_lower = NULL
  strata_prefix = strsplit(names(surv_model$strata[1]), split = '=')[[1]][1]
  for (i in cluster_num) {
    strata_i = paste(strata_prefix, '=', i, sep = '')
    strata_i_index = which(res_sum_subset$strata == strata_i)
    strata_i_year = res_sum_subset$time[strata_i_index]
    strata_i_year_surv_rate = res_sum_subset$surv[strata_i_index]
    strata_i_upper = res_sum_subset$upper[strata_i_index]
    strata_i_lower = res_sum_subset$lower[strata_i_index]
    strata_i_n_year = max(strata_i_year[which(strata_i_year<n)])
    strata_i_n_year_index = which(strata_i_year == strata_i_n_year)
    strata_i_n_year_surv_rate = strata_i_year_surv_rate[strata_i_n_year_index]
    strata_i_n_upper = strata_i_upper[strata_i_n_year_index]
    strata_i_n_lower = strata_i_lower[strata_i_n_year_index]
    n_year = c(n_year, strata_i_n_year)
    n_year_surv_rate = c(n_year_surv_rate, strata_i_n_year_surv_rate)
    n_upper = c(n_upper, strata_i_n_upper)
    n_lower = c(n_lower, strata_i_n_lower)
  }
  n_year_surv_rate_df = data.frame(
    time = n_year,
    surv = n_year_surv_rate,
    upper = n_upper,
    lower = n_lower
  )
  n_year_surv_rate_df$cluster = cluster_num
  return(n_year_surv_rate_df)
}
