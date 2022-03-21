# 从当前目录下的Annotation文件中匹配EMT注释
get_EMT_anno <- function(GeneSymbol){
  emt_anno_path = normalizePath(
    file.path('.', 'Annotation', 'Epithelial-Mesenchymal Transition gene database(1).csv')
  )
  emt_data = read.csv(emt_anno_path)
  emt_data_genesymbol = as.vector(emt_data$Symbol)
  emt_anno = rep('', length(GeneSymbol))
  emt_match_index = match(GeneSymbol, emt_data_genesymbol)
  non_na_index = which(!is.na(emt_match_index))
  if(length(non_na_index) > 0){
    emt_anno[non_na_index] = GeneSymbol[non_na_index]
  }
  return(emt_anno)
}

# 从当前目录下的Annotation文件中匹配Tissue specifity注释
get_tissue_specific_info <- function(GeneSymbol){
  tsi_path = normalizePath(
    file.path('.', 'Annotation', 'gene_annotation.txt')
  )
  tsi_data = read.table(tsi_path, header = T, sep = '\t')
  tsi_data = tsi_data[,c('gene_symbol', 'elevated_tissues', 'tissue_specific_score')]
  
  tsi_genesymbol = as.vector(tsi_data$gene_symbol)
  tsi_tissue = as.vector(tsi_data$elevated_tissues)
  tsi_tissue_score = as.numeric(as.vector(tsi_data$tissue_specific_score))
  
  specific_tissue = rep('', length(GeneSymbol))
  specific_tissue_score = rep('', length(GeneSymbol))
  tsi_match_index = match(GeneSymbol, tsi_genesymbol)
  tsi_match_index_non_na = which(!is.na(tsi_match_index))
  
  if(length(tsi_match_index_non_na)>0){
    specific_tissue_index = tsi_match_index_non_na
    tsi_index = tsi_match_index[specific_tissue_index]
    specific_tissue[specific_tissue_index] = tsi_tissue[tsi_index]
    specific_tissue_score[specific_tissue_index] = tsi_tissue_score[tsi_index]
  }
  zero_index = which(specific_tissue == '0')
  if(length(zero_index)>0){
    specific_tissue[zero_index] = ''
  }
  df = data.frame(specific_tissue, specific_tissue_score)
  return(df)
  
}

# 从当前目录下的Annotation文件中匹配oncogene注释
get_oncogene <- function(GeneSymbol){
  onco_anno_path = normalizePath(
    file.path('.', 'Annotation', 'ongene_human.txt')
  )
  onco_data = read.table(onco_anno_path, sep = '\t', header = T, quote = "")
  onco_data_genesymbol = as.vector(onco_data$OncogeneName)
  onco_anno = rep('', length(GeneSymbol))
  onco_match_index = match(GeneSymbol, onco_data_genesymbol)
  non_na_index = which(!is.na(onco_match_index))
  if(length(non_na_index) > 0){
    onco_anno[non_na_index] = GeneSymbol[non_na_index]
  }
  return(onco_anno)
}


get_fda_human_target_info <- function(GeneSymbol){
  require(readxl)
  fda_human_target_path = normalizePath(
    file.path('.', 'Annotation', 'FDA_human_target_total.xlsx')
  )
  fda_human_target_data = read_excel(fda_human_target_path)
  fda_human_target_data_symbol = as.vector(fda_human_target_data$Symbol)
  parent_ref_name = NULL
  mechanism_of_action = NULL
  target_ChEMBL_ID = NULL
  GeneSymbol_count = length(GeneSymbol)
  for(i in 1:GeneSymbol_count){
    tmp_gene = GeneSymbol[i]
    index_of_match = which(fda_human_target_data_symbol == tmp_gene)
    index_of_match_count = length(index_of_match)
    if(index_of_match_count == 0){
      tmp_parent_ref_name = ''
      tmp_mechanism_of_action = ''
      tmp_target_ChEMBL_ID = ''
    }else{
      tmp_parent_ref_name = fda_human_target_data$PARENT_PREF_NAME[index_of_match]
      tmp_parent_ref_name = paste(tmp_parent_ref_name, collapse = ' | ')
      
      tmp_mechanism_of_action = fda_human_target_data$MECHANISM_OF_ACTION[index_of_match]
      tmp_mechanism_of_action_na_index = which(is.na(tmp_mechanism_of_action))
      if(length(tmp_mechanism_of_action_na_index)>0){
        tmp_mechanism_of_action[tmp_mechanism_of_action_na_index] = 'Unknown'
      }
      tmp_mechanism_of_action = paste(unique(tmp_mechanism_of_action), collapse = ' | ')
      
      tmp_target_ChEMBL_ID = fda_human_target_data$TARGET_CHEMBL_ID[index_of_match]
      tmp_target_ChEMBL_ID_na_index = which(is.na(tmp_target_ChEMBL_ID))
      if(length(tmp_target_ChEMBL_ID_na_index)>0){
        tmp_target_ChEMBL_ID[tmp_target_ChEMBL_ID_na_index] = 'Unknown'
      }
      tmp_target_ChEMBL_ID = paste(unique(tmp_target_ChEMBL_ID), collapse = ' | ')
      
    }
    parent_ref_name = c(parent_ref_name, tmp_parent_ref_name)
    mechanism_of_action = c(mechanism_of_action, tmp_mechanism_of_action)
    target_ChEMBL_ID = c(target_ChEMBL_ID, tmp_target_ChEMBL_ID)
  }
  fda_human_target_df = data.frame(parent_ref_name, mechanism_of_action, target_ChEMBL_ID)
  colnames(fda_human_target_df) = c('PARENT_PREF_NAME', 'Mechanism of action', 'Target ChEMBL ID')
  return(fda_human_target_df)
}


get_tf_tc_info <- function(GeneSymbol){
  require(readxl)
  tf_tc_path = normalizePath(
    file.path('.', 'Annotation', 'TF_TC.xlsx')
  )
  tf_data = read_excel(tf_tc_path, sheet = 'TF')
  tc_data = read_excel(tf_tc_path, sheet = 'TC')
  GeneSymbol_count = length(GeneSymbol)
  GeneSymbol_tf_flag = rep('', GeneSymbol_count)
  GeneSymbol_tc_flag = rep('', GeneSymbol_count)
  
  tf_match_index = match(GeneSymbol, as.vector(unlist(tf_data[,1])))
  tf_match_index_non_na = which(!is.na(tf_match_index))
  if(length(tf_match_index_non_na)>0){
    # tf_match_index_index = tf_match_index[tf_match_index_non_na]
    # GeneSymbol_tf_flag[tf_match_index_non_na] = as.vector(tf_data$Symbol)[tf_match_index_index]
    GeneSymbol_tf_flag[tf_match_index_non_na] = TRUE
  }
  
  tc_match_index = match(GeneSymbol, as.vector(as.vector(unlist(tc_data[,1]))))
  tc_match_index_non_na = which(!is.na(tc_match_index))
  if(length(tc_match_index_non_na)>0){
    # tc_match_index_index = tc_match_index[tc_match_index_non_na]
    # GeneSymbol_tc_flag[tc_match_index_non_na] = as.vector(tc_data$Symbol)[tc_match_index_index]
    GeneSymbol_tc_flag[tc_match_index_non_na] = TRUE
  }
  
  tf_tc_df = data.frame(GeneSymbol_tf_flag, GeneSymbol_tc_flag)
  colnames(tf_tc_df) = c('TF', 'TC')
  return(tf_tc_df)
  
}



# 输入GeneSymbol，t_mat和p_ma，t生成信息汇总表
get_data_statistic <- function(GeneSymbol, t_mat, p_mat){
  
  fda_human_target_df = get_fda_human_target_info(GeneSymbol)
  tf_tc_df = get_tf_tc_info(GeneSymbol)
  
  emt_anno = get_EMT_anno(GeneSymbol)
  tsi_df = get_tissue_specific_info(GeneSymbol)
  oncogene = get_oncogene(GeneSymbol)
  
  tp_ratio_mat = t_mat/p_mat
  pt_ratio_mat = p_mat/t_mat
  min_v = min(t_mat)
  row_count = nrow(t_mat)
  
  timesT <- apply(t_mat, 1, function(x, y){
    length(which(x>y))
  }, y = min_v)
  
  timesP <- apply(p_mat, 1, function(x,y){
    length(which(x>y))
  }, y = min_v)
  
  timesPatient <- NULL
  for(row_index in 1:row_count){
    t_x = as.vector(unlist(t_mat[row_index,]))
    p_x = as.vector(unlist(p_mat[row_index,]))
    tp_count = length(which(t_x>min_v | p_x>min_v))
    timesPatient = c(timesPatient, tp_count)
  }
  
  medianT = apply(t_mat, 1, function(x,y){
    index = which(x>y)
    if(length(index) > 0){
      return(median(x[index]))
    }else{
      return(y)
    }
  }, y = min_v)
  
  medianP = apply(p_mat, 1, function(x,y){
    index = which(x>y)
    if(length(index) > 0){
      return(median(x[index]))
    }else{
      return(y)
    }
  }, y = min_v)
  
  aveT = apply(t_mat, 1, function(x,y){
    index = which(x>y)
    if(length(index) > 0){
      return(mean(x[index]))
    }else{
      return(y)
    }
  }, y = min_v)
  
  aveP = apply(p_mat, 1, function(x,y){
    index = which(x>y)
    if(length(index) > 0){
      return(mean(x[index]))
    }else{
      return(y)
    }
  }, y = min_v)
  
  MaxT = apply(t_mat, 1, function(x, y){
    index = which(x>y)
    if(length(index) > 0){
      return(max(x[index]))
    }else{
      return(y)
    }
  }, y = min_v)
  
  MaxP = apply(p_mat, 1, function(x, y){
    index = which(x>y)
    if(length(index) > 0){
      return(max(x[index]))
    }else{
      return(y)
    }
  }, y = min_v)
  
  
  tp_ratio_3 = apply(tp_ratio_mat, 1, function(x){
    length(which(x>3))
  })
  
  tp_ratio_5 = apply(tp_ratio_mat, 1, function(x){
    length(which(x>5))
  })
  
  tp_ratio_10 = apply(tp_ratio_mat, 1, function(x){
    length(which(x>10))
  })
  
  
  
  pt_ratio_3 = apply(pt_ratio_mat, 1, function(x){
    length(which(x>3))
  })
  
  pt_ratio_5 = apply(pt_ratio_mat, 1, function(x){
    length(which(x>5))
  })
  
  pt_ratio_10 = apply(pt_ratio_mat, 1, function(x){
    length(which(x>10))
  })
  
  T_up_regu = tp_ratio_3/timesPatient
  T_up_regu_all = tp_ratio_3/ncol(t_mat)
  
  T_down_regu = pt_ratio_3/timesPatient
  T_down_regu_all = pt_ratio_3/ncol(t_mat)
  
  variation_ratio = (tp_ratio_3 - pt_ratio_3)/timesPatient
  
  df = data.frame(
    GeneSymbol,
    fda_human_target_df,
    tf_tc_df,
    emt_anno,
    tsi_df,
    oncogene,
    timesT,
    timesP,
    timesPatient,
    medianT,
    medianP,
    aveT,
    aveP,
    MaxT,
    MaxP,
    tp_ratio_3,
    tp_ratio_5,
    tp_ratio_10,
    pt_ratio_3,
    pt_ratio_5,
    pt_ratio_10,
    T_up_regu,
    T_up_regu_all,
    T_down_regu,
    T_down_regu_all,
    variation_ratio 
  )
  headers = c(
    'Symbol',
    'PARENT_PREF_NAME',
    'Mechanism of action', 
    'Target ChEMBL ID',
    'TF',
    'TC',
    'EMT',
    'Tissue specific',
    'Tissue specific score',	
    'oncogene'	,
    'timesT'	,
    'timesP'	,
    'timesPatient',	
    'medianT'	,
    'medianP'	,
    'aveT',
    'aveP',
    'MaxT',
    'MaxP',
    'T/P>3',
    'T/P>5',
    'T/P>10',
    'P/T>3'	,
    'P/T>5'	,
    'P/T>10'	,
    'T中上调比例'	,
    'T上调总人数占比'	,
    'T中下调比例'	,
    'T下调总人数占比',
    '变化比值数'
  )
  # headers = strsplit(headers, split = '[\t]')[[1]]
  colnames(df) = headers
  return(df)
}

# 内存回收
jgc <- function(){
  gc()
  .jcall("java/lang/System", method = "gc")
}  


# 按GeneSymbol合并列表
merge_list_by_genesymbol <- function(lst, by_name){
  lst_count = length(lst)
  merge_df = lst[[1]]
  for(i in 2:lst_count){
    merge_df = merge(merge_df, lst[[i]], by = by_name, all = TRUE)
  }
  Symbol = as.vector(unlist(merge_df[,1]))
  mat = as.matrix(merge_df[,-1])
  na_index = which(is.na(mat))
  if(length(na_index)>0){
    mat[na_index] = 0
  }
  df = data.frame(Symbol, mat)
  return(df)
}

# 删除全0行
remove_rows_with_all_zero <- function(df){
  df_value = df[,-1]
  df_rowsum = rowSums(df_value)
  zero_index = which(df_rowsum == 0)
  if(length(zero_index)>0){
    df = df[-zero_index,]
  }
  return(df)
}

# 获取gene symbol和OS，DFS的关系
get_continue_df_result <- function(data_path, ref_path, data_type, survival_type){
  # data_path = p_path
  data_df = read.csv(data_path)
  Symbol = as.vector(data_df$Symbol)
  Symbol_count = length(Symbol)
  Value = data_df[,-1]
  Value_colnames = colnames(Value)
  Value_exp_codes = apply(data.frame(Value_colnames), 1, function(x){
    strsplit(x, split = '_')[[1]][1]
  })
  
  # ref_path = ref_file_path
  ref_file_data_tp = read_excel(ref_path, sheet = 'proteome_fno')
  ref_file_data_t_clinical = read_excel(ref_path, sheet = 'proteome_T_fno')
  
  if(data_type == 'P'){
    p_index_of_match = match(Value_exp_codes, ref_file_data_tp$Proteome_fno)
    t_index_of_match = p_index_of_match - 1
    Value_exp_codes = ref_file_data_tp$Proteome_fno[t_index_of_match]
  }
  
  sample_index_of_match = match(Value_exp_codes, ref_file_data_t_clinical$Proteome_fno)
  ref_clinical_info = ref_file_data_t_clinical[sample_index_of_match,]
  
  if(data_type == 'T'){
    continue_df_colnames = c('PsigCountinCount-T', 'MinPvalue-T', 'HRminP-T', 'H_Expr_Count-T', 'Best_Cutpoint-T')
  }
  if(data_type == 'P'){
    continue_df_colnames = c('PsigCountinCount-P', 'MinPvalue-T', 'HRminP-T', 'H_Expr_Count-P', 'Best_Cutpoint-P')
  }
  if(data_type == 'TP'){
    continue_df_colnames = c('PsigCountinCount-RATIO', 'MinPvalue-RATIO', 'HRminP-RATIO', 'H_Expr_Count-RATIO', 'Best_Cutpoint-RATIO')
  }
  continue_df_colnames = paste(survival_type, continue_df_colnames, sep = '-')
  continue_df = NULL
  patient_number = ncol(Value)
  for(gene_index in 1:Symbol_count){
    gene_value = as.vector(unlist(Value[gene_index,]))
    gene_value_duplicated = duplicated(gene_value)
    if(length(which(!gene_value_duplicated)) == 1){
      continue_vector = rep('NA', 3)
    }else{
      attempt_number = 5
      gene_pvalue = NULL
      gene_HR = NULL
      na_break = F
      cutpoint_vector = NULL
      cutpoint_high_expr_count = NULL
      cutpoint_low_expr_count = NULL
      while(attempt_number <= (patient_number-5)){
        cutpoint = sort(gene_value)[attempt_number]
        label_flags = gene_value > cutpoint
        label_v = rep(1, patient_number)
        label_v[which(label_flags)] = 2 # high expression
        info_df = data.frame(ref_clinical_info, label_v)
        
        cutpoint_vector = c(cutpoint_vector, cutpoint)
        cutpoint_low_expr_count = c(cutpoint_low_expr_count, length(which(label_v==1)))
        cutpoint_high_expr_count = c(cutpoint_high_expr_count, length(which(label_v==2)))
        
        
        if(survival_type == 'OS'){
          coxph_reg =coxph(
            Surv(as.numeric(OS), as.numeric(OS_status))~ as.numeric(label_v)+as.numeric(Age)+as.factor(TNM_stage)+as.factor(Chemotherapy_status), 
            data = info_df
          )
        }else{
          coxph_reg =coxph(
            Surv(as.numeric(DFS), as.numeric(DFS_status))~ as.numeric(label_v)+as.numeric(Age)+as.factor(TNM_stage)+as.factor(Chemotherapy_status), 
            data = info_df
          )
        }
        coxph_reg_summary = summary(coxph_reg)
        coefficients = data.frame(coxph_reg_summary$coefficients)
        HR = coefficients[1,2]
        pvalue = coefficients[1,5]
        if(!is.na(pvalue)){
          gene_HR = c(gene_HR, HR)
          gene_pvalue = c(gene_pvalue, pvalue)
        }
        attempt_number = attempt_number + 1
      }
      if(length(gene_pvalue) == 0){
        continue_vector = rep('NA', 3)
      }else{
        continue_vector = get_pvalue_and_hr_in_max_continues(gene_pvalue, gene_HR)
      }
    }
    
    if('NA' %in% continue_vector){
      high_expr_count = 'NA'
      low_expr_count = 'NA'
      best_cutpoint = 'NA'
    }else{
      min_pvalue = continue_vector[2]
      min_HR = continue_vector[3]
      min_index = which(gene_pvalue == min_pvalue & gene_HR == min_HR)
      min_cutpoint = cutpoint_vector[min_index[1]]
      high_expr_count = cutpoint_high_expr_count[min_index[1]]
      low_expr_count = cutpoint_low_expr_count[min_index[1]]
      best_cutpoint = cutpoint_vector[min_index[1]]
    }
    continue_vector_extend =  c(continue_vector, high_expr_count, best_cutpoint)
    continue_df = rbind(continue_df, continue_vector_extend)
    if(gene_index %% 500 == 0 | gene_index == Symbol_count){
      cat('\n', data_type, '_', survival_type, '->', gene_index, '\n')
    }
  }
  
  colnames(continue_df) = continue_df_colnames
  continue_df_result = data.frame(Symbol, continue_df)
  return(continue_df_result)
}

# 获取最佳pvalue和HR
get_pvalue_and_hr_in_max_continues <- function(gene_pvalue, gene_HR){
  pvalue_count = length(gene_pvalue)
  pvalue_continue_index_set = NULL
  pvalue_continue_index = NULL
  continue_flag = F
  for(i in 1:pvalue_count){
    pvalue_i = gene_pvalue[i]
    if(pvalue_i<0.05 & i<pvalue_count){
      continue_flag = T
      pvalue_continue_index = c(pvalue_continue_index, i)
    }else{
      if(continue_flag){
        if(pvalue_i<0.05 & i==pvalue_count){
          pvalue_continue_index = c(pvalue_continue_index, i)
        }
        pvalue_continue_index_count = length(pvalue_continue_index)
        pvalue_continue_index_set_count = length(pvalue_continue_index_set)
        if(pvalue_continue_index_count > pvalue_continue_index_set_count){
          pvalue_continue_index_set = pvalue_continue_index
          pvalue_continue_index = NULL
        }
      }
      continue_flag = F
      
    }
  }
  if(pvalue_i<0.05 & i==pvalue_count & length(pvalue_continue_index_set) == 0){
    pvalue_continue_index_set = i
  }
  
  
  pvalue_continue_index_set_count = length(pvalue_continue_index_set)
  if(pvalue_continue_index_set_count>0){
    continue_sig_count = pvalue_continue_index_set_count
    gene_pvalue_subset = gene_pvalue[pvalue_continue_index_set]
    gene_HR_subset = gene_HR[pvalue_continue_index_set]
    min_pvalue = min(gene_pvalue_subset)
    min_HR = gene_HR_subset[which.min(gene_pvalue_subset)]
  }else{
    continue_sig_count = 'NA'
    min_pvalue = min(gene_pvalue)
    min_HR = gene_HR[which.min(gene_pvalue)]
  }
  continue_vector = c(continue_sig_count, min_pvalue, min_HR)
  return(continue_vector)
  
}


# 输出基本的汇总信息
write_summary_statistcs_file <- function(
  lauren_type, 
  csv_ifot_data, 
  csv_us_data, 
  csv_s_data,
  ref_file_data_tp,
  ref_exps_D_T
){
  index_of_match = match(ref_exps_D_T, ref_file_data_tp$Proteome_fno)
  t_index = index_of_match
  t_index_count = length(t_index)
  ref_exps = NULL
  for(i in 1:t_index_count){
    ref_exps = c(
      ref_exps, 
      ref_file_data_tp$Proteome_fno[t_index[i]],
      ref_file_data_tp$Proteome_fno[t_index[i]+1]
    )
  }
  csv_ifot_data_colnames = colnames(csv_ifot_data)
  csv_ifot_data_colnames_flag = apply(data.frame(csv_ifot_data_colnames), 1, function(x){
    x = strsplit(x, split = '_')[[1]]
    x[1]
  })
  index_of_match = match(ref_exps, csv_ifot_data_colnames_flag)
  index_of_match = c(1, index_of_match)
  ifot_df_of_match = remove_rows_with_all_zero(csv_ifot_data[,index_of_match])
  us_df_of_match = csv_us_data[match(ifot_df_of_match$Symbol, csv_us_data$Symbol),index_of_match]
  s_df_of_match = csv_s_data[match(ifot_df_of_match$Symbol, csv_s_data$Symbol),index_of_match]
  
  
  Symbol = as.vector(unlist(ifot_df_of_match[,1]))
  Values = as.matrix(ifot_df_of_match[,-1])
  index_of_zero = which(Values==0)
  replace_v = min(Values[-index_of_zero])/10
  Values[index_of_zero] = replace_v
  
  ### build T/P mat ###
  Values_col = ncol(Values)
  t_index = seq(1, Values_col, 2)
  p_index = seq(2, Values_col, 2)
  
  t_mat = Values[,t_index]
  p_mat = Values[,p_index]
  tp_mat = t_mat/p_mat
  
  t_output_name = paste(lauren_type, 'T_D1.csv', sep = '_')
  t_output_path = normalizePath(file.path(OUTPUT_DIR, t_output_name), mustWork = F)
  write.csv(data.frame(Symbol, t_mat), t_output_path, row.names = F)
  
  p_output_name = paste(lauren_type, 'P_D1.csv', sep = '_')
  p_output_path = normalizePath(file.path(OUTPUT_DIR, p_output_name), mustWork = F)
  write.csv(data.frame(Symbol, p_mat), p_output_path, row.names = F)
  
  tp_output_name = paste(lauren_type, 'TP_D1.csv', sep = '_')
  tp_output_path = normalizePath(file.path(OUTPUT_DIR, tp_output_name), mustWork = F)
  write.csv(data.frame(Symbol, tp_mat), tp_output_path, row.names = F)
  
  output_name = paste(lauren_type, length(t_index), 'pairs_summary_statistic.csv', sep = '_')
  output_path = normalizePath(file.path(OUTPUT_DIR, output_name), mustWork = F)
  ### process ###
  stat_df = get_data_statistic(Symbol, t_mat, p_mat)
  
  # uspept data
  us_pept_file_data = us_df_of_match
  us_pept_t_mat = us_pept_file_data[,-1][,t_index]
  us_pept_p_mat = us_pept_file_data[,-1][,p_index]
  us_pept_t_max = apply(us_pept_t_mat, 1, max)
  us_pept_p_max = apply(us_pept_p_mat, 1, max)
  
  # spept data
  s_pept_file_data = s_df_of_match
  s_pept_t_mat = s_pept_file_data[,-1][,t_index]
  s_pept_p_mat = s_pept_file_data[,-1][,p_index]
  s_pept_t_max = apply(s_pept_t_mat, 1, max)
  s_pept_p_max = apply(s_pept_p_mat, 1, max)
  
  pept_df = data.frame(us_pept_t_max, us_pept_p_max, s_pept_t_max, s_pept_p_max)
  colnames(pept_df) = c('Max_USPept_T', 'Max_USPept_P', 'Max_SPept_T', 'Max_SPept_P')
  
  stat_df = data.frame(stat_df, pept_df)
  write.csv(stat_df, output_path, row.names = F)
  message("Completed: ", output_path)
  # a = ref_file_data_tp[match(ref_exps, ref_file_data_tp$Proteome_fno),]
}

############################################################################
write_summary_statistcs_file_20190325 <- function(
  data_type, 
  csv_ifot_data, 
  csv_us_data, 
  csv_s_data,
  ref_exps_t,
  ref_exps_p
){
  ifot_df_of_match = csv_ifot_data
  us_df_of_match = csv_us_data
  s_df_of_match = csv_s_data
  Symbol = as.vector(unlist(ifot_df_of_match[,1]))
  Values = as.matrix(ifot_df_of_match[,-1])
  index_of_zero = which(Values==0)
  replace_v = min(Values[-index_of_zero])/10
  Values[index_of_zero] = replace_v
  ### build T/P mat ###
  csv_ifot_data_colnames = colnames(Values)
  csv_ifot_data_colnames_flag = apply(data.frame(csv_ifot_data_colnames), 1, function(x){
    x = strsplit(x, split = '_')[[1]]
    x[1]
  })
  t_index = match(ref_exps_t, csv_ifot_data_colnames_flag)
  p_index = match(ref_exps_p, csv_ifot_data_colnames_flag)
  
  t_mat = Values[,t_index]
  p_mat = Values[,p_index]
  tp_mat = t_mat/p_mat
  
  t_output_name = paste(data_type, 'T_D1.csv', sep = '_')
  t_output_path = normalizePath(file.path(OUTPUT_DIR, t_output_name), mustWork = F)
  write.csv(data.frame(Symbol, t_mat), t_output_path, row.names = F)
  
  p_output_name = paste(data_type, 'P_D1.csv', sep = '_')
  p_output_path = normalizePath(file.path(OUTPUT_DIR, p_output_name), mustWork = F)
  write.csv(data.frame(Symbol, p_mat), p_output_path, row.names = F)
  
  tp_output_name = paste(data_type, 'TP_D1.csv', sep = '_')
  tp_output_path = normalizePath(file.path(OUTPUT_DIR, tp_output_name), mustWork = F)
  write.csv(data.frame(Symbol, tp_mat), tp_output_path, row.names = F)
  
  output_name = paste(data_type, length(t_index), 'pairs_summary_statistic.csv', sep = '_')
  output_path = normalizePath(file.path(OUTPUT_DIR, output_name), mustWork = F)
  ### process ###
  stat_df = get_data_statistic(Symbol, t_mat, p_mat)
  stat_df_colnames = colnames(stat_df)
  
  # uspept data
  us_pept_file_data = us_df_of_match
  us_pept_t_mat = us_pept_file_data[,-1][,t_index]
  us_pept_p_mat = us_pept_file_data[,-1][,p_index]
  us_pept_t_max = apply(us_pept_t_mat, 1, max)
  us_pept_p_max = apply(us_pept_p_mat, 1, max)
  
  # spept data
  s_pept_file_data = s_df_of_match
  s_pept_t_mat = s_pept_file_data[,-1][,t_index]
  s_pept_p_mat = s_pept_file_data[,-1][,p_index]
  s_pept_t_max = apply(s_pept_t_mat, 1, max)
  s_pept_p_max = apply(s_pept_p_mat, 1, max)
  
  pept_df = data.frame(us_pept_t_max, us_pept_p_max, s_pept_t_max, s_pept_p_max)
  colnames(pept_df) = c('Max_USPept_T', 'Max_USPept_P', 'Max_SPept_T', 'Max_SPept_P')
  pept_df_colnames = colnames(pept_df)
  
  stat_df = data.frame(stat_df, pept_df)
  stat_df_colnames = c(stat_df_colnames, pept_df_colnames)
  colnames(stat_df) = stat_df_colnames
  # write.csv(stat_df, output_path, row.names = F, fileEncoding = 'GBK2312')
  require('readr')
  write_excel_csv(stat_df, output_path)
  message("Completed: ", output_path)
  # a = ref_file_data_tp[match(ref_exps, ref_file_data_tp$Proteome_fno),]
}





# 获取gene symbol和OS，DFS的关系
get_continue_df_result_20190327 <- function(data_path, ref_path, data_class, data_type, survival_type){
  cat('\n', 'get_continue_df_result_20190327', '\n')
  # data_path = p_path
  data_df = read.csv(data_path)
  Symbol = as.vector(data_df$Symbol)
  Symbol_count = length(Symbol)
  Value = data_df[,-1]
  Value_colnames = colnames(Value)
  Value_exp_codes = apply(data.frame(Value_colnames), 1, function(x){
    strsplit(x, split = '_')[[1]][1]
  })
  
  # ref_path = ref_file_path
  if(data_class == 'ALL'){
    clinical_data = read_excel(ref_path, sheet = 'filter_249')
  }else{
    clinical_data = read_excel(ref_path, sheet = 'filter_278')
  }
  
  
  if(data_type == 'P'){
    clinical_data_index = match(Value_exp_codes, clinical_data$Proteome_fno_P)
    clinical_data_subset = clinical_data[clinical_data_index,]
  }else{
    clinical_data_index = match(Value_exp_codes, clinical_data$Proteome_fno_T)
    clinical_data_subset = clinical_data[clinical_data_index,]
  }
  
  
  
  ref_clinical_info = clinical_data_subset
  
  if(data_type == 'T'){
    continue_df_colnames = c('PsigCountinCount-T', 'MinPvalue-T', 'HRminP-T', 'H_Expr_Count-T', 'Best_Cutpoint-T')
  }
  if(data_type == 'P'){
    continue_df_colnames = c('PsigCountinCount-P', 'MinPvalue-T', 'HRminP-T', 'H_Expr_Count-P', 'Best_Cutpoint-P')
  }
  if(data_type == 'TP'){
    continue_df_colnames = c('PsigCountinCount-RATIO', 'MinPvalue-RATIO', 'HRminP-RATIO', 'H_Expr_Count-RATIO', 'Best_Cutpoint-RATIO')
  }
  continue_df_colnames = paste(survival_type, continue_df_colnames, sep = '-')
  continue_df = NULL
  patient_number = ncol(Value)
  for(gene_index in 1:Symbol_count){
    gene_value = as.vector(unlist(Value[gene_index,]))
    gene_value_duplicated = duplicated(gene_value)
    if(length(which(!gene_value_duplicated)) == 1){
      continue_vector = rep('NA', 3)
    }else{
      attempt_number = 5
      gene_pvalue = NULL
      gene_HR = NULL
      na_break = F
      cutpoint_vector = NULL
      cutpoint_high_expr_count = NULL
      cutpoint_low_expr_count = NULL
      while(attempt_number <= (patient_number-5)){
        cutpoint = sort(gene_value)[attempt_number]
        label_flags = gene_value > cutpoint
        label_v = rep(1, patient_number)
        label_v[which(label_flags)] = 2 # high expression
        info_df = data.frame(ref_clinical_info, label_v)
        
        cutpoint_vector = c(cutpoint_vector, cutpoint)
        cutpoint_low_expr_count = c(cutpoint_low_expr_count, length(which(label_v==1)))
        cutpoint_high_expr_count = c(cutpoint_high_expr_count, length(which(label_v==2)))
        
        
        if(survival_type == 'OS'){
          coxph_reg =coxph(
            # Surv(as.numeric(OS_time), as.numeric(OS_status))~ as.numeric(label_v)+as.numeric(Age)+as.factor(TNM_stage)+as.factor(Chemotherapy_status),
            Surv(as.numeric(OS_time), as.numeric(OS_status))~ as.numeric(label_v)+as.numeric(Age)+as.factor(TNM_stage), 
            data = info_df
          )
        }else{
          coxph_reg =coxph(
            # Surv(as.numeric(DFS), as.numeric(DFS_status))~ as.numeric(label_v)+as.numeric(Age)+as.factor(TNM_stage2)+as.factor(Chemotherapy_status),
            Surv(as.numeric(DFS_time), as.numeric(DFS_status))~ as.numeric(label_v)+as.numeric(Age)+as.factor(TNM_stage), 
            data = info_df
          )
        }
        coxph_reg_summary = summary(coxph_reg)
        coefficients = data.frame(coxph_reg_summary$coefficients)
        HR = coefficients[1,2]
        pvalue = coefficients[1,5]
        if(!is.na(pvalue)){
          gene_HR = c(gene_HR, HR)
          gene_pvalue = c(gene_pvalue, pvalue)
        }
        attempt_number = attempt_number + 1
      }
      if(length(gene_pvalue) == 0){
        continue_vector = rep('NA', 3)
      }else{
        continue_vector = get_pvalue_and_hr_in_max_continues(gene_pvalue, gene_HR)
      }
    }
    
    if('NA' %in% continue_vector){
      high_expr_count = 'NA'
      low_expr_count = 'NA'
      best_cutpoint = 'NA'
    }else{
      min_pvalue = continue_vector[2]
      min_HR = continue_vector[3]
      min_index = which(gene_pvalue == min_pvalue & gene_HR == min_HR)
      min_cutpoint = cutpoint_vector[min_index[1]]
      high_expr_count = cutpoint_high_expr_count[min_index[1]]
      low_expr_count = cutpoint_low_expr_count[min_index[1]]
      best_cutpoint = cutpoint_vector[min_index[1]]
    }
    continue_vector_extend =  c(continue_vector, high_expr_count, best_cutpoint)
    continue_df = rbind(continue_df, continue_vector_extend)
    if(gene_index %% 500 == 0 | gene_index == Symbol_count){
      cat('\n', data_type, '_', survival_type, '->', gene_index, '\n')
    }
  }
  
  colnames(continue_df) = continue_df_colnames
  continue_df_result = data.frame(Symbol, continue_df)
  return(continue_df_result)
}


#################################### function zone ##############################################################
if(TRUE){
  format_nairen_data <- function(df){
    # df = dim_ifot_data
    Symbol = as.vector(unlist(df[,1]))
    Mat = as.matrix(df[,-1])
    na_index = which(is.na(Mat))
    if(length(na_index)>0){
      Mat[na_index] = 0
    }
    df_new = data.frame(Symbol, Mat)
    return(df_new)
  }
  
  clear_duplicated_data <- function(df){
    # df = dim_ifot_data
    Symbol = as.vector(unlist(df[,1]))
    Mat = as.matrix(df[,-1])
    Mat_col = ncol(Mat)
    for (i in 1:Mat_col) {
      x = Mat[,i]
      x_duplicated_index = which(
        (x > 0)
        & duplicated(x)
      )
      Mat[x_duplicated_index,i] = 0
    }
    Mat_rowsum = rowSums(Mat)
    df_new = df[Mat_rowsum>0,]
    return(df_new)
  }
  
  # clear genes with correlation = 1
  clear_duplicated_data_1 <- function(df){
    # df = dim_ifot_data
    Symbol = as.vector(unlist(df[,1]))
    Mat = as.matrix(df[,-1])
    rownames(Mat) = Symbol
    gene_cormat = cor(t(Mat))
    rownames(gene_cormat) = Symbol
    colnames(gene_cormat) = Symbol
    gene_cormat_row = nrow(gene_cormat)
    correlation_equal_index = NULL
    for (i in 1:gene_cormat_row) {
      x = gene_cormat[,i]
      x_cor1_index = which(x == 1)
      x_cor1_index_count = length(x_cor1_index)
      if(x_cor1_index_count>1){
        # x_cor1_index = setdiff(x_cor1_index, i)
        correlation_equal_index = union(correlation_equal_index, x_cor1_index)
      }
    }
    correlation_equal_mat = Mat[correlation_equal_index,]
    correlation_equal_mat_row = nrow(correlation_equal_mat)
    correlation_equal_mat_col = ncol(correlation_equal_mat)
    duplicated_index = NULL
    for (i in 1:correlation_equal_mat_row) {
      x_i = correlation_equal_mat[i,]
      for (j in 1:correlation_equal_mat_row) {
        if(i != j & !(j %in% duplicated_index)){
          x_j = correlation_equal_mat[j,]
          if(length(which(x_i == x_j)) == correlation_equal_mat_col){
            duplicated_index = c(duplicated_index, j)
            # print(duplicated_index)
          }
        }
      }
    }
    duplicated_index = unique(duplicated_index)
    duplicated_genes = rownames(correlation_equal_mat)[duplicated_index]
    
    duplicated_match_index = match(duplicated_genes, Symbol)
    df_new = df[-duplicated_match_index,]
    return(df_new)
  }
  
  
  get_bp_list = function(df){
    # df = dim_ifot_data
    Symbol = as.vector(unlist(df[,1]))
    Mat = as.matrix(df[,-1])
    Mat_col = ncol(Mat)
    bp_list = list()
    for (i in 1:Mat_col) {
      x = Mat[,i]
      x = x[x>0]
      bp_list[[i]] = log10(x)
    }
    return(bp_list)
  }
  
  # 1
  get_clinical_data_subset <- function(expnames, clinical_data){
    t_index = match(expnames, clinical_data$Proteome_fno_T)
    t_index_non_na_index = which(!is.na(t_index))
    t_index_subset = t_index[t_index_non_na_index]
    p_index = match(expnames, clinical_data$Proteome_fno_P)
    p_index_non_na_index = which(!is.na(p_index))
    p_index_subset = p_index[p_index_non_na_index]
    return(clinical_data[t_index_subset,])
  }
  
  # 2
  get_df_with_no_HB_ACT_KRT <- function(df){
    mat_subset_symbol = as.vector(unlist(df[,1]))
    delete_index_HB = which(grepl('^HB', mat_subset_symbol))
    delete_index_ACT = which(grepl('^ACT', mat_subset_symbol))
    delete_index_KRT = which(grepl('^KRT', mat_subset_symbol))
    delete_index = union(delete_index_HB, union(delete_index_ACT, delete_index_KRT))
    if(length(delete_index)>0){
      df = df[-delete_index,]
    }
    return(df)
  }
  
  # 2
  get_df_with_no_HB_ACT_KRT <- function(df){
    mat_subset_symbol = as.vector(unlist(df[,1]))
    delete_index_HB = which(grepl('^HB', mat_subset_symbol))
    delete_index_ACT = which(grepl('^ACT', mat_subset_symbol))
    delete_index_KRT = which(grepl('^KRT', mat_subset_symbol))
    delete_index = union(delete_index_HB, union(delete_index_ACT, delete_index_KRT))
    if(length(delete_index)>0){
      df = df[-delete_index,]
    }
    return(df)
  }
  
  get_df_with_no_HB_ACT_KRT_HIST <- function(df){
    mat_subset_symbol = as.vector(unlist(df[,1]))
    delete_index_HB = which(grepl('^HB', mat_subset_symbol))
    delete_index_ACT = which(grepl('^ACT', mat_subset_symbol))
    delete_index_KRT = which(grepl('^KRT', mat_subset_symbol))
    delete_index_HIST = which(grepl('^HIST', mat_subset_symbol))
    delete_index = union(delete_index_HB, union(delete_index_ACT, delete_index_KRT))
    delete_index = union(delete_index, delete_index_HIST)
    if(length(delete_index)>0){
      df = df[-delete_index,]
    }
    return(df)
  }
  
  
  # 3 
  remove_expdata_for_mat <- function(mat, expnames){
    mat_colnames = colnames(mat)
    del_index = match(expnames, mat_colnames)
    mat = mat[,-del_index]
    mat_rowsums = rowSums(mat)
    kept_index = which(mat_rowsums>0)
    return(mat[kept_index,])
  }
  
  # 4
  remove_expdata_for_df <- function(df, expnames){
    mat = df[,-1]
    mat_colnames = colnames(mat)
    del_index = match(expnames, mat_colnames)
    mat = mat[,-del_index]
    mat_rowsums = rowSums(mat)
    kept_index = which(mat_rowsums>0)
    return(df[kept_index,])
  }
  
  # 5
  get_fot_df_by_exps <- function(fot_df, exps){
    Symbol = as.vector(unlist(fot_df[,1]))
    mat = fot_df[,-1]
    mat_colnames = colnames(mat)
    match_index = match(exps, mat_colnames)
    mat_subset = mat[,match_index]
    mat_subset_rowsums = rowSums(mat_subset)
    kept_index = which(mat_subset_rowsums>0)
    result_df = data.frame(Symbol, mat_subset)[kept_index,]
    return(result_df)
  }
  
  # 6
  get_topx_df <- function(df, topx){
    mat_value = df[,-1]
    mat_value_col = ncol(mat_value)
    union_index = NULL
    topx_seq = seq(1, topx)
    for(i in 1:mat_value_col){
      x = as.vector(unlist(mat_value[,i]))
      x_order = order(x, decreasing = T)
      # x_order_index = match(topx_seq, x_order)
      x_order_index = x_order[topx_seq]
      # if(2 %in% x_order_index){
      #   print(i)
      # }
      union_index = union(union_index, x_order_index)
    }
    topx_df = df[union_index,]
    return(topx_df)
  }
  
  
  # 6
  get_topx_df_1 <- function(df, topx){
    mat_value = df[,-1]
    mat_value_col = ncol(mat_value)
    topx_seq = seq(1, topx)
    for(i in 1:mat_value_col){
      x = as.vector(unlist(mat_value[,i]))
      x_order = order(x, decreasing = T)
      x_order_index = x_order[topx_seq]
      y = x
      y[-x_order_index] = 0
      mat_value[,i] = y
    }
    kept_index = which(rowSums(mat_value)>0)
    topx_df = data.frame(
      Symbol = as.vector(unlist(df[,1])),
      mat_value
    )[kept_index,]
    return(topx_df)
  }
  
  # 8
  fot_df_correction <- function(df){
    Symbol = as.vector(unlist(df[,1]))
    Value = df[,-1]
    Value_row = nrow(Value)
    Value_col = ncol(Value)
    correction_value = matrix(0, Value_row, Value_col)
    colnames(correction_value) = colnames(Value)
    for(i in 1:Value_col){
      i_v = as.vector(unlist(Value[,i]))
      i_v_fot = i_v/sum(i_v)*1e5
      correction_value[,i] = i_v_fot
    }
    result_df = data.frame(Symbol, correction_value)
    return(result_df)
  }
  
  # 9
  get_xlim <- function(mat_2d){
    xlim <- c(floor(min(mat_2d[,1]))-1, ceiling(max(mat_2d[,1]))+1)
    return(xlim)
  }
  
  # 10
  get_ylim <- function(mat_2d){
    ylim <- c(floor(min(mat_2d[,2]))-1, ceiling(max(mat_2d[,2]))+1)
    return(ylim)
  }
  
  # 11  
  get_outliers_count <- function(x){
    x_scale = as.vector(scale(x))
    bp = boxplot(x_scale, plot = F)
    outliers_count = length(bp$out)
    return(outliers_count)
  }
  
  # 12
  get_max_cluster <- function(x){
    mydata = as.matrix(x)
    # print(x)
    wss_vector <- NULL
    wss_vector[1] <- (nrow(mydata)-1)*sum(apply(mydata,2,var))
    for(i in 2:12){
      wss <- tryCatch(
        {
          set.seed(1024)
          sum(kmeans(mydata, centers=i)$withinss)
        },
        error=function(cond) {
          # message("Here's the original error message:")
          # message(cond)
          return('error')
        },
        warning=function(cond) {
          # message("Here's the original warning message:")
          # message(cond)
          return('warning')
        }
      )
      if(typeof(wss) == "double"){
        wss_vector[i] = wss
      }else{
        break
      }
      
    }
    wss_cumsum = cumsum(wss_vector)
    for(i in 2:length(wss_cumsum)){
      x1 = wss_cumsum[i-1]
      x2 = wss_cumsum[i]
      if(1-x1/x2<0.05){
        break
      }
    }
    k = i
    return(k)
  }
  
  # 13
  get_rank_df <- function(df){
    Symbol = as.vector(df[,1])
    Value = df[,-1]
    Value_row = nrow(Value)
    Value_col = ncol(Value)
    rank_mat = matrix(0, Value_row, Value_col, dimnames = list(Symbol, colnames(Value)))
    for(i in 1:Value_col){
      x = as.vector(unlist(Value[,i]))
      index = which(x==0)
      x_order = order(x, decreasing = T)
      x_rank = rep(0, length(x))
      x_rank[x_order] = seq(1, length(x))
      x_rank[index] = NA
      rank_mat[,i] = x_rank
    }
    rank_df = data.frame(Symbol, rank_mat)
    return(rank_df)
  }
  
  # 14
  impute_df_with_random <- function(df){
    Symbol = as.vector(unlist(df[,1]))
    Mat = as.matrix(df[,-1])
    Mat_new = matrix(0, nrow(Mat), ncol(Mat), dimnames = list(Symbol, colnames(Mat)))
    Mat_col = ncol(Mat)
    for(i in 1:Mat_col){
      x = as.vector(unlist(Mat[,i]))
      zero_index = which(x==0)
      non_zero_index = which(x>0)
      x1 = x[non_zero_index]
      x1_log2 = log2(x1)
      x1_log2_mean = mean(x1_log2)
      x1_log2_sd = sd(x1_log2)
      x2_log2_mean = x1_log2_mean - 1.96*x1_log2_sd
      x2_log2_sd = 0.25*x1_log2_sd
      x2_log2 = rnorm(length(zero_index), x2_log2_mean, x2_log2_sd)
      Mat_new[non_zero_index,i] = x1_log2
      Mat_new[zero_index,i] = x2_log2
    }
    df_new = data.frame(Symbol, Mat_new)
    return(df_new)
  }
  
  scale_df_correction <- function(df){
    Symbol = as.vector(df$Symbol)
    Value = df[,-1]
    Value_scale = t(scale(t(Value), center = T, scale = T))
    return(data.frame(Symbol, Value_scale))
  }
  
  
  scale_df_correction_1 <- function(df){
    Symbol = as.vector(df$Symbol)
    Value = df[,-1]
    row_sd = apply(Value, 1, sd)
    Value_scale = Value/row_sd
    return(data.frame(Symbol, Value_scale))
  }
  
}




#################################### function zone ##############################################################




coxph_ggforest <- function(coxph_data, class_flag, data_flag, main){
  library(survival)
  library(survminer)
  # coxph_data = dim_summary_data
  coxph_data$OS_time = as.numeric(coxph_data$OS_time)
  coxph_data$OS_status = as.numeric(coxph_data$OS_status)
  coxph_data$Age_label = as.factor(coxph_data$Age_label)
  coxph_data$Gender = as.factor(coxph_data$Gender)
  coxph_data$TNM_stage_new = as.factor(coxph_data$TNM_stage_new)
  coxph_data$Chemotherapy_status = as.factor(coxph_data$Chemotherapy_status)
  
  if(data_flag==1){ # 分类变量
    coxph_data$Cluster = as.factor(coxph_data$Cluster)
    coxph_data$MCluster = as.factor(coxph_data$MCluster)
  }else{ # 连续变量
    coxph_data$Cluster = as.numeric(coxph_data$Cluster)
    coxph_data$MCluster = as.numeric(coxph_data$MCluster)
  }
  
  if(class_flag == 'Cluster'){
    cox_model = coxph(
      Surv(OS_time, OS_status) ~ Age_label + Gender + TNM_stage_new + Chemotherapy_status + Cluster, data = coxph_data
    )
  }else if(class_flag == 'MCluster'){
    cox_model = coxph(
      Surv(OS_time, OS_status) ~ Age_label + Gender + TNM_stage_new + Chemotherapy_status + MCluster, data = coxph_data
    )
  }else if(class_flag == 'ACT1'){
    cox_model = coxph(
      Surv(OS_time, OS_status) ~ as.numeric(Age_label) + Gender + TNM_stage_new + ACT1, data = coxph_data
    )
  }else if(class_flag == 'ACT2'){
    cox_model = coxph(
      Surv(OS_time, OS_status) ~ as.numeric(Age_label) + Gender + TNM_stage_new + ACT2, data = coxph_data
    )
  }else{
    cox_model = coxph(
      Surv(OS_time, OS_status) ~ Age_label + Gender + TNM_stage_new  + Chemotherapy_status, data = coxph_data
    )
  }
  plot(ggforest(cox_model, data = coxph_data, main = main))
}








coxph_ggforest_dfs <- function(coxph_data, class_flag, data_flag, main){
  library(survival)
  library(survminer)
  # coxph_data = dim_summary_data
  coxph_data$OS_time = as.numeric(coxph_data$DFS_time)
  coxph_data$OS_status = as.numeric(coxph_data$DFS_status)
  coxph_data$Age_label = as.factor(coxph_data$Age_label)
  coxph_data$Gender = as.factor(coxph_data$Gender)
  coxph_data$TNM_stage_new = as.factor(coxph_data$TNM_stage_new)
  coxph_data$Chemotherapy_status = as.factor(coxph_data$Chemotherapy_status)
  
  if(data_flag==1){ # 分类变量
    coxph_data$Cluster = as.factor(coxph_data$Cluster)
    coxph_data$MCluster = as.factor(coxph_data$MCluster)
  }else{ # 连续变量
    coxph_data$Cluster = as.numeric(coxph_data$Cluster)
    coxph_data$MCluster = as.numeric(coxph_data$MCluster)
  }
  
  if(class_flag == 'Cluster'){
    cox_model = coxph(
      Surv(OS_time, OS_status) ~ Age_label + Gender + TNM_stage_new + Chemotherapy_status + Cluster, data = coxph_data
    )
  }else if(class_flag == 'MCluster'){
    cox_model = coxph(
      Surv(OS_time, OS_status) ~ Age_label + Gender + TNM_stage_new + Chemotherapy_status + MCluster, data = coxph_data
    )
  }else if(class_flag == 'ACT1'){
    cox_model = coxph(
      Surv(OS_time, OS_status) ~ as.numeric(Age_label) + Gender + TNM_stage_new + ACT1, data = coxph_data
    )
  }else if(class_flag == 'ACT2'){
    cox_model = coxph(
      Surv(OS_time, OS_status) ~ as.numeric(Age_label) + Gender + TNM_stage_new + ACT2, data = coxph_data
    )
  }else{
    cox_model = coxph(
      Surv(OS_time, OS_status) ~ Age_label + Gender + TNM_stage_new + Chemotherapy_status, data = coxph_data
    )
  }
  plot(ggforest(cox_model, data = coxph_data, main = main))
}









get_predict_class <- function(ref_group, tmp_group){
  tmp_group_count = length(tmp_group)
  tmp_group_ref = tmp_group[-tmp_group_count]
  tmp_group_pred = tmp_group[tmp_group_count]
  map_table = table(ref_group, tmp_group_ref)
  map_table_v = map_table[,tmp_group_pred]
  predict_class = as.numeric(which.max(map_table_v))
  return(predict_class)
}


get_features_mat <- function(predict_data, reduncdant_genes, top=1000, features){
  # predict_data
  # features = features_all
  # reduncdant_genes
  
  # clear reduncdant genes
  genes = as.vector(unlist(predict_data[,1]))
  reduncdant_index = as.vector(na.omit(match(reduncdant_genes, genes)))
  if(length(reduncdant_index) > 0){
    predict_data = predict_data[-reduncdant_index,]
  }
  
  
  # get top data and fot correction
  # top = 1000
  predict_data_top = get_topx_df_1(predict_data, top)
  predict_data_top = fot_df_correction(predict_data_top)
  
  
  # get features matrix
  genes = as.vector(unlist(predict_data_top[,1]))
  mat = as.matrix(predict_data_top[,-1])
  rownames(mat) = genes
  mat_col = ncol(mat)
  
  features_index = match(features, genes)
  features_na_index = which(is.na(features_index))
  if(length(features_na_index)>0){
    print(features[features_na_index])
    stop('Having NA features')
  }
  features_matrix = mat[features_index,]
  
  return(features_matrix)
  # features_matrix = matrix(
  #   0, 
  #   nrow = length(features), ncol = ncol(mat), 
  #   dimnames = list(features, colnames(mat))
  # )
  # 
  # for (i in 1:mat_col) {
  #   x = as.vector(unlist(mat[,i]))
  #   x_features_index = match(features, genes)
  #   x_features_na_index = which(is.na(x_features_index))
  #   if(length(x_features_na_index)>0){
  #     features_matrix[-x_features_na_index,i] = x[x_features_index[-x_features_na_index]]
  #   }else{
  #     features_matrix[,i] = x[x_features_index]
  #   }
  # }
  
  
}


# kept genes occurred in different set
get_df_with_genes_occured_in_different_set <- function(df, clinical_data, occurrence = 0.1){
  # occurrence = 0.1
  # df = global_ifot_data_1
  # clinical_data = clinical_data_global
  df_mat = as.matrix(df[,-1])
  df_mat_row = nrow(df_mat)
  df_group = as.factor(as.vector(clinical_data$UID))
  df_group_num = length(table(df_group))
  df_group_min_id = as.vector(table(df_group)*occurrence)
  df_kept_index = NULL
  for (i in 1:df_mat_row) {
    x = as.vector(unlist(df_mat[i,]))
    x_id_v = tapply(x, df_group, function(y){length(which(y>0))})
    flag = x_id_v >= df_group_min_id
    if(length(which(flag))==df_group_num){
      df_kept_index = c(df_kept_index, i)
    }
  }
  df_subset = df[df_kept_index,]
  return(df_subset)
}

# kept genes occurred in different set
get_occurrence_freq_of_gene_in_different_set <- function(df, clinical_data, occurrence = 0.1){
  # occurrence = 0.1
  # df = global_ifot_data_1
  # clinical_data = clinical_data_global
  df_mat = as.matrix(df[,-1])
  df_mat_row = nrow(df_mat)
  df_group = as.factor(as.vector(clinical_data$UID))
  df_group_num = length(table(df_group))
  df_group_min_id = as.vector(table(df_group)*occurrence)
  occurrence_freq = NULL
  for (i in 1:df_mat_row) {
    x = as.vector(unlist(df_mat[i,]))
    x_id_v = tapply(x, df_group, function(y){length(which(y>0))})
    flag = x_id_v >= df_group_min_id
    occurrence_freq = c(occurrence_freq, length(which(flag)))
    # if(length(which(flag))==df_group_num){
    #   df_kept_index = c(df_kept_index, i)
    # }
  }
  # df_subset = df[df_kept_index,]
  return(occurrence_freq)
}


#################################### function zone ##############################################################


