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
}

# library(caret)
# set.seed(seed0)
# d.index = createDataPartition(clinical_data_all$OS_time, p = 0.5)
# d.index = unlist(d.index)



clinical_data_discovery = clinical_data_all[clinical_data_all$LDT==0,]
# *** dim_ifot_data_discovery *** #



dim_ifot_data_discovery = get_fot_df_by_exps(dim_ifot_data, clinical_data_discovery$Proteome_fno_T)

ProtCount = apply(dim_ifot_data_discovery[,-1], 2, function(x){length(which(x>0))})
print(quantile(ProtCount, seq(0, 1, 0.1)))

dim_ifot_data_discovery_1 = dim_ifot_data_discovery[-match(
  redundant_genes, dim_ifot_data_discovery$Symbol
),]

dim_ifot_data_discovery_topN = get_topx_df_1(dim_ifot_data_discovery_1, 1100)
dim_ifot_data_discovery_topN = fot_df_correction(dim_ifot_data_discovery_topN)

# *** Top1000 data *** #
all_data = data.frame(
  Symbol = dim_ifot_data_discovery_topN$Symbol,
  as.matrix(dim_ifot_data_discovery_topN[,-1])
)

write.csv(
  all_data, 'SCLC-75-Profile-3460-20220222.csv', row.names = F
)

write.csv(
  dim_ifot_data_discovery, 'SCLC-75-Profile-6415-20220222.csv', row.names = F
)

write.csv(
  dim_ifot_data_discovery_1, 'SCLC-75-Profile-6108-20220222.csv', row.names = F
)
