#### function zone ####
if(T){
  # recalculate_fot_matrix
  recalculate_fot_matrix <-  function(df){
    symbol = as.vector(unlist(df[,1]))
    mat = as.matrix(df[,-1])
    mat.exp_sums = colSums(mat)
    mat_new = t(t(mat)/mat.exp_sums)*1e5
    mat_new.colnames = apply(data.frame(colnames(mat_new)), 1, function(x){
      if(grepl('MEI', x)){
        return(strsplit(x, '\\.')[[1]][1])
      }else{
        return(x)
      }
    })
    colnames(mat_new) = mat_new.colnames
    mat_new.rowsums = rowSums(mat_new)
    result_df = data.frame(
      Symbol = symbol, 
      mat_new
    )[mat_new.rowsums>0,]
    return(result_df)
  }
  
  # get_fot_df_by_exps
  get_fot_df_by_exps <- function(fot_df, exps){
    symbol = as.vector(unlist(fot_df[,1]))
    mat = fot_df[,-1]
    mat_colnames = colnames(mat)
    match_index = match(exps, mat_colnames)
    mat_subset = mat[,match_index]
    mat_subset_rowsums = rowSums(mat_subset)
    kept_index = which(mat_subset_rowsums>0)
    result_df = data.frame(
      Symbol = symbol, 
      mat_subset
    )[kept_index,]
    return(result_df)
  }
  
  # clear_duplicated_data
  # mean, sd
  clear_duplicated_data <- function(df){
    symbol = as.vector(unlist(df[,1]))
    mat = as.matrix(df[,-1])
    mat.row_sd = apply(mat, 1, sd)
    mat.row_mean = apply(mat, 1, mean)
    row_sd.dup_index = which(duplicated(mat.row_sd))
    row_mean.dup_index = which(duplicated(mat.row_mean))
    dup_index = intersect(row_sd.dup_index, row_mean.dup_index)
    if(length(dup_index)>0){
      print(paste0('Duplicated proteins: ', length(dup_index)))
      df_new = data.frame(
        Symbol = symbol,
        mat
      )[-dup_index,]
      dup_prots = symbol[dup_index]
    }else{
      df_new = df
      dup_prots = NA
    }
    
    results_lst = list(
      Data = df_new, Dup_prots = dup_prots
    )
    
    return(results_lst)
  }
  
  
  # clear_high_intensity_prots
  clear_high_intensity_prots <- function(df, cumsum_cutoff = 0.8){
    # clear genes with high intensity 
    df_genes = as.vector(unlist(df[, 1]))
    df_mat = df[, -1]
    rownames(df_mat) = df_genes
    
    df_mat_row_mean = apply(df_mat, 1, mean)
    df_mat_row_mean_sort_desc_cumsum = cumsum(
      sort(df_mat_row_mean, decreasing = T)
    )
    df_mat_row_mean_sort_desc_cumsum_probs = df_mat_row_mean_sort_desc_cumsum/max(
      df_mat_row_mean_sort_desc_cumsum
    )
    plot(
      seq(1, length(df_mat_row_mean_sort_desc_cumsum_probs)), 
      df_mat_row_mean_sort_desc_cumsum_probs,
      xlab = 'Gene rank',
      ylab = 'Cumulative probability'
    )
    abline(h = cumsum_cutoff, col = 2, lty = 3)
    high_intensity_genes = names(
      df_mat_row_mean_sort_desc_cumsum_probs
    )[df_mat_row_mean_sort_desc_cumsum_probs < cumsum_cutoff]
    
    high_intensity_genes_index = match(high_intensity_genes, df_genes)
    
    
    if(length(high_intensity_genes)>0){
      print(paste0('High_intensity_genes: ', length(high_intensity_genes)))
      df_new = data.frame(
        Symbol = df_genes,
        df_mat
      )[-high_intensity_genes_index,]
      high_intensity_prots = df_genes[high_intensity_genes_index]
    }else{
      df_new = df
      high_intensity_prots = NA
    }
    
    results_lst = list(
      Data = df_new, High_intensity_prots = high_intensity_prots
    )
    
    return(results_lst)
  }
  
  
  # get_topx_df
  get_topx_df <- function(df, topx){
    if(is.na(topx)){
      return(df)
    }
    mat_value = df[,-1]
    mat_value_col = ncol(mat_value)
    topx_seq = seq(1, topx)
    for(i in 1:mat_value_col){
      x = as.vector(unlist(mat_value[,i]))
      x_order = order(x, decreasing = T)
      x_order_index = x_order[topx_seq]
      y = rep(0, length(x))
      y[x_order_index] = x[x_order_index]
      mat_value[,i] = y
    }
    kept_index = which(rowSums(mat_value)>0)
    topx_df = data.frame(
      Symbol = as.vector(unlist(df[,1])),
      mat_value
    )[kept_index,]
    return(topx_df)
  }
  
  
  
}