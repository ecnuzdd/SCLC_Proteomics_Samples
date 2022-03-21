
library("RColorBrewer")
col.palette = brewer.pal(n = 12, name = "Paired")



if(T){
  cluster_num = 3
  
  print(length(features))
  match_index = match(as.vector(features), all_data$Symbol)
  seed = seed0 = seed1 = 129
  
  discovery_data = all_data[na.omit(match_index),]
  # *** log1p *** #
  discovery_data_mat = log1p(as.matrix(discovery_data[,-1]))
  rownames(discovery_data_mat) = discovery_data[, 1]
  discovery_data_mat_col_id = apply(discovery_data_mat, 2, function(x){
    length(which(x > 0))
  })
  which(discovery_data_mat_col_id  == 0)
  discovery_data_mat_col_id_probs = discovery_data_mat_col_id/length(features)
  print(quantile(discovery_data_mat_col_id_probs))
  
  features_count = length(features)
  features_kept_index = which(discovery_data_mat_col_id > 0)
  print(length(features_kept_index))
  
  discovery_time = clinical_data_discovery$OS_time[features_kept_index]
  discovery_status = clinical_data_discovery$OS_status[features_kept_index]
  print(paste('discovery set:'))
  print(dim(discovery_data_mat))
  print(range(discovery_data_mat))
  
  colors = brewer.pal(n = 9, name = "Set1")

  # for(cluster_num in 5:9){
  nmf_mat = discovery_data_mat[, features_kept_index]
  plot_flag = T
  
  if(plot_flag){
    pdf_dir = normalizePath(
      file.path(BASE_DIR, 
                'Dset_NMF_1026'),
      mustWork = F)
    if(!dir.exists(pdf_dir)){
      dir.create(pdf_dir)
    }
    pdf(
      normalizePath(file.path(
        pdf_dir,
        paste0("Dset_NMF_Cluster=", cluster_num, "_features=", nrow(nmf_mat), ".pdf")
      ), mustWork = F)
    )
  }
  
  nmf_result_list_discovery <- get_nmf_result_list(
    main_title = paste(
      "Discovery set (N=", ncol(nmf_mat), 
      ") | seed=", seed, '\n',
      " | Features=", nrow(nmf_mat), 
      # " | cv=", cv_cutpoint, 
      sep = ''),
    seed = seed1,
    nmf_mat,
    cluster_num = cluster_num,
    nrun = 50,
    time = discovery_time[features_kept_index],
    status = discovery_status[features_kept_index],
    n_year = 3,
    ylab = "Probability of Overall Survival",
    cluster_colors = colors
  )
  if(plot_flag){
    dev.off()
  }
  
  if(F){

    write.xlsx(
      survival_data,
      normalizePath(file.path(
        pdf_dir,
        paste0("Dset_NMF_Cluster=", cluster_num, "_features=", nrow(nmf_mat), ".xlsx")
      ), mustWork = F),
      row.names = F
    )
  }
  
  
}

if(T){
  clinical_data.d = clinical_data_discovery
  clinical_data.d$Cluster = paste('C', nmf_result_list_discovery$cluster, sep = '')
  write.csv(
    clinical_data.d, 
    normalizePath(file.path(pdf_dir, 'clinical_data_simple_1010.csv')),
    row.names = F
  )
  write.xlsx(
    clinical_data.d, 
    normalizePath(file.path(pdf_dir, 'clinical_data_simple_1010.xlsx')),
    row.names = F
  )
  
  TMC_DIR = pdf_dir
  output_path = normalizePath(
    file.path(TMC_DIR, paste('output', '_', nrow(clinical_data.d), '_5ysr.csv', sep = '')),
    mustWork = F
  )
  
  
  yl_tmc_data = clinical_data.d
  yl_tms_cluster = sort(unique(yl_tmc_data$Cluster))
  yl_tms_cluster_count = length(yl_tms_cluster)
  cluster_num = yl_tms_cluster_count
  
  os_pvalue_list = NULL
  chemo_5ysr_list = NULL
  nochemo_5ysr_list = NULL
  total_count = NULL
  chemo_count = NULL
  nochemo_count = NULL
  pdf_path = normalizePath(
    file.path(TMC_DIR, paste('output', '_', nrow(clinical_data.d), '_5.pdf', sep = '')),
    mustWork = F
  )
  pdf(pdf_path)
  for (j in 1:yl_tms_cluster_count) {
    
    j.clinical_data = yl_tmc_data[yl_tmc_data$Cluster==yl_tms_cluster[j],]
    chemo_count = c(chemo_count, length(which(j.clinical_data$Chemotherapy_status==1)))
    nochemo_count = c(nochemo_count, length(which(j.clinical_data$Chemotherapy_status==0)))
    total_count = c(total_count, nrow(j.clinical_data))
    
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
      main_title = paste('Survival curve (OS) | TNM = II-III | ', yl_tms_cluster[j], sep = '')
      os_fit_ggsurvplot <- ggsurvplot(
        os_fit, data = survival_data,
        title = main_title,
        surv.scale = 'percent',
        xlab = 'Years',
        legend = 'right', 
        legend.title = "Chemo status",
        legend.labs = c('NoChemo', 'Chemo'),
        pval = TRUE,
        palette = c('grey', 'blue'),
        # risk.table.col = "Group",
        risk.table = TRUE, 
        risk.table.y.text.col = TRUE,
        surv.median.line = 'hv',
        break.time.by = 1,
        xlim = c(0,8)
      )
      print(os_fit_ggsurvplot)
      os_pvalue = surv_pvalue(os_fit)$pval
      os_pvalue = round(os_pvalue, 3)
      print(paste('os_p_value: ', os_pvalue, sep = ''))
      table_5ysr_each = get_n_year_surv_rate_using_clinical_data(
        na.omit(survival_data), n = 5, cluster_label = names(table(survival_data$Cluster))
      )
      print(table_5ysr_each)
      chemo_5ysr = table_5ysr_each$surv[table_5ysr_each$group==2]
      nochemo_5ysr = table_5ysr_each$surv[table_5ysr_each$group==1]
      
    }else{
      chemo_5ysr = NA
      nochemo_5ysr = NA
      os_pvalue = NA
    }
    
    
    
    chemo_5ysr_list = c(chemo_5ysr_list, chemo_5ysr)
    nochemo_5ysr_list = c(nochemo_5ysr_list, nochemo_5ysr)
    os_pvalue_list = c(os_pvalue_list, os_pvalue)
    
    
  }
  dev.off()
  
  
  df = data.frame(
    Node = yl_tms_cluster,
    Total_count = total_count,
    Chemo_count = chemo_count,
    NoChemo_count = nochemo_count,
    Pvalue = os_pvalue_list,
    Chemo5ysr = chemo_5ysr_list,
    NoChemo5ysr = nochemo_5ysr_list
  )
  df$Diff = round(df$Chemo5ysr - df$NoChemo5ysr, 3)
  
  # All
  if(T){
    
    # *** KM curve for OS *** #
    survival_data = data.frame(
      OS_time = yl_tmc_data$OS_time,
      OS_status = yl_tmc_data$OS_status,
      Cluster = yl_tmc_data$Cluster
    )
    
    os_fit <- survfit(
      Surv(as.numeric(survival_data$OS_time), as.numeric(survival_data$OS_status)) ~ Cluster, data = survival_data
    )
    main_title = paste('Survival curve (OS) | ', 'N = ', nrow(survival_data), sep = '')
    os_fit_ggsurvplot <- ggsurvplot(
      os_fit, data = survival_data,
      title = main_title,
      surv.scale = 'percent',
      xlab = 'Years',
      legend = 'right', 
      legend.title = paste("NMF-", ' (', cluster_num, ')', sep = ''),
      legend.labs = c(
        # 'NoChemo',
        paste('', sort(unique(survival_data$Cluster)), sep = '')
      ),
      pval = TRUE,
      palette = c('grey', col.palette, '#000000')[-1],
      # risk.table.col = "Group",
      risk.table = TRUE, 
      risk.table.y.text.col = TRUE,
      surv.median.line = 'hv',
      break.time.by = 1,
      xlim = c(0,8)
    )
    pdf(
      normalizePath(
        file.path(TMC_DIR, paste('output', '_', nrow(clinical_data.d), '_all.pdf', sep = '')),
        mustWork = F
      )
      # height = 13, width = 13
    )
    print(os_fit_ggsurvplot)
    dev.off()
    os_pvalue = surv_pvalue(os_fit)$pval
    print(paste('os_p_value: ', os_pvalue, sep = ''))
    
    table_5ysr = get_n_year_surv_rate_using_clinical_data(
      na.omit(survival_data), n = 5, cluster_label = names(table(survival_data$Cluster))
    )
    print(table_5ysr)
    
    table_5ysr.match_index = match(df$Node, table_5ysr$cluster)
    df$All5ysr = table_5ysr$surv[table_5ysr.match_index]
    
    desc_order = order(df$All5ysr, decreasing = T)
    df = df[desc_order,]
    # df$Cluster = paste('C', seq(1, nrow(df)), sep = '')
    # write.xlsx(df, output_path, row.names = F, sheetName = paste('D564', sep = ''), append = T)
    write.csv(df, output_path, row.names = F)
    print(df)
  }
  
  
}

###############################
if(T){
  output_path = normalizePath(
    file.path(TMC_DIR, paste('output', '_', nrow(clinical_data.d), '_3ysr.csv', sep = '')),
    mustWork = F
  )
  
  yl_tmc_data = clinical_data.d
  yl_tms_cluster = sort(unique(yl_tmc_data$Cluster))
  yl_tms_cluster_count = length(yl_tms_cluster)
  cluster_num = yl_tms_cluster_count
  
  os_pvalue_list = NULL
  chemo_5ysr_list = NULL
  nochemo_5ysr_list = NULL
  total_count = NULL
  chemo_count = NULL
  nochemo_count = NULL
  pdf_path = normalizePath(
    file.path(TMC_DIR, paste('output', '_', nrow(clinical_data.d), '_3.pdf', sep = '')),
    mustWork = F
  )
  pdf(pdf_path)
  for (j in 1:yl_tms_cluster_count) {
    
    j.clinical_data = yl_tmc_data[yl_tmc_data$Cluster==yl_tms_cluster[j],]
    chemo_count = c(chemo_count, length(which(j.clinical_data$Chemotherapy_status==1)))
    nochemo_count = c(nochemo_count, length(which(j.clinical_data$Chemotherapy_status==0)))
    total_count = c(total_count, nrow(j.clinical_data))
    
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
      main_title = paste('Survival curve (OS) | ', yl_tms_cluster[j], sep = '')
      os_fit_ggsurvplot <- ggsurvplot(
        os_fit, data = survival_data,
        title = main_title,
        surv.scale = 'percent',
        xlab = 'Years',
        legend = 'right', 
        legend.title = "Chemo status",
        legend.labs = c('NoChemo', 'Chemo'),
        pval = TRUE,
        palette = c('grey', 'blue'),
        # risk.table.col = "Group",
        risk.table = TRUE, 
        risk.table.y.text.col = TRUE,
        surv.median.line = 'hv',
        break.time.by = 1,
        xlim = c(0,12)
      )
      print(os_fit_ggsurvplot)
      os_pvalue = surv_pvalue(os_fit)$pval
      os_pvalue = round(os_pvalue, 3)
      print(paste('os_p_value: ', os_pvalue, sep = ''))
      table_3ysr_each = get_n_year_surv_rate_using_clinical_data(
        na.omit(survival_data), n = 3, cluster_label = names(table(survival_data$Cluster))
      )
      print(table_3ysr_each)
      chemo_3ysr = table_3ysr_each$surv[table_3ysr_each$group==2]
      nochemo_3ysr = table_3ysr_each$surv[table_3ysr_each$group==1]
      
    }else{
      chemo_3ysr = NA
      nochemo_3ysr = NA
      os_pvalue = NA
    }
    
    
    
    chemo_3ysr_list = c(chemo_3ysr_list, chemo_3ysr)
    nochemo_3ysr_list = c(nochemo_3ysr_list, nochemo_3ysr)
    os_pvalue_list = c(os_pvalue_list, os_pvalue)
    
    
  }
  dev.off()
  
  
  df = data.frame(
    Node = yl_tms_cluster,
    Total_count = total_count,
    Chemo_count = chemo_count,
    NoChemo_count = nochemo_count,
    Pvalue = os_pvalue_list,
    Chemo3ysr = chemo_3ysr_list,
    NoChemo3ysr = nochemo_3ysr_list
  )
  df$Diff = round(df$Chemo3ysr - df$NoChemo3ysr, 3)
  
  # All
  if(T){
    
    # *** KM curve for OS *** #
    survival_data = data.frame(
      OS_time = yl_tmc_data$OS_time,
      OS_status = yl_tmc_data$OS_status,
      Cluster = yl_tmc_data$Cluster
    )
    
    os_fit <- survfit(
      Surv(as.numeric(survival_data$OS_time), as.numeric(survival_data$OS_status)) ~ Cluster, data = survival_data
    )
    main_title = paste('Survival curve (OS) | ', 'D', nrow(survival_data), sep = '')
    os_fit_ggsurvplot <- ggsurvplot(
      os_fit, data = survival_data,
      title = main_title,
      surv.scale = 'percent',
      xlab = 'Years',
      legend = 'right', 
      legend.title = paste("NMF-", ' (', cluster_num, ')', sep = ''),
      legend.labs = c(
        # 'NoChemo',
        paste('C', sort(unique(survival_data$Cluster)), sep = '')
      ),
      pval = TRUE,
      palette = c('grey', col.palette, '#000000')[-1],
      # risk.table.col = "Group",
      risk.table = TRUE, 
      risk.table.y.text.col = TRUE,
      surv.median.line = 'hv',
      break.time.by = 1,
      xlim = c(0,12)
    )
    pdf(
      normalizePath(
        file.path(TMC_DIR, paste('output', '_', nrow(clinical_data.d), '_all3.pdf', sep = '')),
        mustWork = F
      ),
      height = 13, width = 13
    )
    print(os_fit_ggsurvplot)
    dev.off()
    os_pvalue = surv_pvalue(os_fit)$pval
    print(paste('os_p_value: ', os_pvalue, sep = ''))
    
    table_3ysr = get_n_year_surv_rate_using_clinical_data(
      na.omit(survival_data), n = 5, cluster_label = names(table(survival_data$Cluster))
    )
    print(table_3ysr)
    
    table_3ysr.match_index = match(df$Node, table_3ysr$cluster)
    df$All3ysr = table_3ysr$surv[table_3ysr.match_index]
    
    desc_order = order(df$All3ysr, decreasing = T)
    df = df[desc_order,]
    # df$Cluster = paste('C', seq(1, nrow(df)), sep = '')
    # write.xlsx(df, output_path, row.names = F, sheetName = paste('D564', sep = ''), append = T)
    write.csv(df, output_path, row.names = F)
    print(df)
  }
  
  
}
