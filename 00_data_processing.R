# *** load library and path *** #
if(T){
  ### work dictory ###
  getwd()
  # BASE_DIR = "/home/ddzhan/SCLC/20201026"
  BASE_DIR = "/home/ddzhan/SCLC/20210726-SX"
  BASE_DIR = normalizePath(BASE_DIR)
  setwd(BASE_DIR)
  
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
  
  ### custom library ###
  custom_library_path = normalizePath(
    file.path(BASE_DIR, '..', 'library_directory', 'whole_gc_analysis_library.R')
  )
  source(custom_library_path, encoding = 'UTF-8', echo=TRUE)
  
  # metadata
  DATA_DIR = normalizePath(
    file.path(BASE_DIR, 'input')
  )
  clinical_exp_path = normalizePath(
    file.path(
      DATA_DIR, 'sclc_clinical_data_simple_20201013.xlsx'
    )
  )
  
  # quantification data
  dim_ifot_path = normalizePath(
    file.path(
      DATA_DIR, 'sclc_ifot_aft_qa.csv')
  )
}

# *** load data *** #
if(T){
  clinical_exp_data_0 = read_excel(clinical_exp_path)
  clinical_exp_data_0 = clinical_exp_data_0[clinical_exp_data_0$LDT!=1,]
  clinical_exp_data_0 = clinical_exp_data_0[clinical_exp_data_0$Proteome_fno_T!='Exp078306',]
  
  
  dim_ifot_data_0 = read.csv(dim_ifot_path)
  dim_ifot_data_0 = format_nairen_data(dim_ifot_data_0)
  
  dim_ifot_data_exps = colnames(dim_ifot_data_0)
  invalid_exps = setdiff(clinical_exp_data_0$Proteome_fno_T, dim_ifot_data_exps)
  if(length(invalid_exps)>0){
    invalid_exps_index = match(invalid_exps, clinical_exp_data_0$Proteome_fno_T)
    clinical_exp_data_0 = clinical_exp_data_0[-invalid_exps_index,]
  }
  
  dim_ifot_data_0 = get_fot_df_by_exps(dim_ifot_data_0, clinical_exp_data_0$Proteome_fno_T)
  dim_ifot_data_0 = fot_df_correction(dim_ifot_data_0)

  # 1-保留鉴定数大于2000
  # 2-去冗余蛋白
  # 3-取top3000
  # 4-重算ifot
}

# *** Batch info *** #
if(T){
  # *** ID > 2000 *** #
  ProtCount = apply(dim_ifot_data_0[,-1], 2, function(x){length(which(x>0))})
  quantile(ProtCount, seq(0, 1, 0.1))
  FOT_Median = apply(dim_ifot_data_0[,-1], 2, function(x){
    x = x[x>0]
    median(x)
  })
  clinical_exp_data_1 = data.frame(
    clinical_exp_data_0,
    ProtCount,
    FOT_Median
  )
  print(dim(clinical_exp_data_1))
  # clinical_exp_data_1 = clinical_exp_data_0[clinical_exp_data_0$Tumor_cells_percentage >=5,]
  print(dim(clinical_exp_data_1))
  # clinical_exp_data_1 = clinical_exp_data_1[clinical_exp_data_1$ProtCount > 2000,]
  print(dim(clinical_exp_data_1))
  
  # ***  *** #
  clinical_data_all = clinical_exp_data_1
}

if(T){
  # *** Global set *** #
  clinical_data_global = clinical_data_all
  print("Global set:")
  print(dim(clinical_data_global))
  print(quantile(clinical_data_global$OS_time))
  
  # *** Discovery set *** #
  clinical_data_discovery = clinical_data_global
  print("Discovery set:")
  print(dim(clinical_data_discovery))
  print(quantile(clinical_data_discovery$OS_time))
  
  # *** Prediction set *** #
  clinical_data_prediction = clinical_data_global
  print("Prediction set:")
  print(dim(clinical_data_prediction))
  print(quantile(clinical_data_prediction$OS_time))
}

if(T){
  if(T){
    # clear genes with duplication
    dset_ifot_data_0 = get_fot_df_by_exps(dim_ifot_data_0, clinical_data_discovery$Proteome_fno_T)
    pset_ifot_data_0 = get_fot_df_by_exps(dim_ifot_data_0, clinical_data_prediction$Proteome_fno_T)
    
    dset_ifot_data_1 = clear_duplicated_data(dset_ifot_data_0)
    duplicated_genes = setdiff(as.vector(dset_ifot_data_0$Symbol), as.vector(dset_ifot_data_1$Symbol))
    print(duplicated_genes)
    print(paste('duplicated_genes:', length(duplicated_genes)))
    
    # clear genes with high intensity 
    dset_ifot_data_1_genes = as.vector(unlist(dset_ifot_data_1[, 1]))
    dset_ifot_data_1_mat = dset_ifot_data_1[, -1]
    rownames(dset_ifot_data_1_mat) = dset_ifot_data_1_genes
    
    dset_ifot_data_1_mat_row_mean = apply(dset_ifot_data_1_mat, 1, mean)
    dset_ifot_data_1_mat_row_mean_sort_desc_cumsum = cumsum(
      sort(dset_ifot_data_1_mat_row_mean, decreasing = T)
    )
    dset_ifot_data_1_mat_row_mean_sort_desc_cumsum_probs = dset_ifot_data_1_mat_row_mean_sort_desc_cumsum/max(
      dset_ifot_data_1_mat_row_mean_sort_desc_cumsum
    )
    plot(
      seq(1, length(dset_ifot_data_1_mat_row_mean_sort_desc_cumsum_probs)), 
      dset_ifot_data_1_mat_row_mean_sort_desc_cumsum_probs,
      xlab = 'Gene rank',
      ylab = 'Cumulative probability'
    )
    cumsum_cutoff = 0.80
    abline(h = cumsum_cutoff, col = 2, lty = 3)
    high_intensity_genes = names(dset_ifot_data_1_mat_row_mean_sort_desc_cumsum_probs)[
      dset_ifot_data_1_mat_row_mean_sort_desc_cumsum_probs < cumsum_cutoff]
    print(sort(high_intensity_genes))
    print(paste('high_intensity_genes:', length(high_intensity_genes)))
    high_intensity_genes_index = match(high_intensity_genes, dset_ifot_data_1_genes)
    dset_ifot_data_2 = dset_ifot_data_1[-high_intensity_genes_index, ]
    
    redundant_genes = union(duplicated_genes, high_intensity_genes)
    print(paste('redundant_genes:', length(redundant_genes)))
  }
  
  if(T){
    DIM_DATA_DIR = normalizePath(
      file.path(BASE_DIR, 'input'), mustWork = F
    )
    if(!dir.exists(DIM_DATA_DIR)){
      dir.create(DIM_DATA_DIR)
    }
    
    # dset_ifot_data_0
    dset_ifot_data_0_path = normalizePath(
      file.path(DIM_DATA_DIR, 'dset_ifot_data_0.csv'), mustWork = F
    )
    write.csv(dset_ifot_data_0, dset_ifot_data_0_path, row.names = F)
    
    # pset_ifot_data_0
    pset_ifot_data_0_path = normalizePath(
      file.path(DIM_DATA_DIR, 'pset_ifot_data_0.csv'), mustWork = F
    )
    write.csv(pset_ifot_data_0, pset_ifot_data_0_path, row.names = F)
    
    clinical_data_global_path = normalizePath(
      file.path(DIM_DATA_DIR, 'clinical_data_global.csv'), mustWork = F
    )
    write.csv(clinical_data_global, clinical_data_global_path, row.names = F)
    
    clinical_data_discovery_path = normalizePath(
      file.path(DIM_DATA_DIR, 'clinical_data_discovery.csv'), mustWork = F
    )
    write.csv(clinical_data_discovery, clinical_data_discovery_path, row.names = F)
    
    clinical_data_prediction_path = normalizePath(
      file.path(DIM_DATA_DIR, 'clinical_data_prediction.csv'), mustWork = F
    )
    write.csv(clinical_data_prediction, clinical_data_prediction_path, row.names = F)
    
    write.csv(
      data.frame(duplicated_genes), 
      normalizePath(
        file.path(DIM_DATA_DIR, 'duplicated_genes.csv'), mustWork = F
      ), 
      row.names = F
    )
    write.csv(
      data.frame(high_intensity_genes), 
      normalizePath(
        file.path(DIM_DATA_DIR, 'high_intensity_genes.csv'), mustWork = F
      ),
      row.names = F
    )
    write.csv(
      data.frame(redundant_genes), 
      normalizePath(
        file.path(DIM_DATA_DIR, 'redundant_genes.csv'), mustWork = F
      ), 
      row.names = F
    )
  }
}



