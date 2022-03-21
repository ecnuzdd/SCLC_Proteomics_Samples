
all_data_mat = as.matrix(all_data[, -1])
row.names(all_data_mat) = all_data$Symbol
all_data_mat_cal = all_data_mat
all_data_mat_cal = log1p(all_data_mat)
range(all_data_mat_cal)
range(all_data_mat)

# cv ↓

all_data_id = apply(all_data_mat_cal, 1, function(x){length(which(x>0))})
all_data_mean = apply(all_data_mat_cal, 1, mean)
all_data_sd = apply(all_data_mat_cal, 1, sd)
all_data_cv = all_data_sd/all_data_mean
summary(all_data_cv)

# 1.8 497
# 1.9 449
if(T){
  # cv_features_ind = which(all_data_cv > 2 & all_data_id > ncol(all_data_mat)*0.15) # 152
  cv_features_ind = which(all_data_cv > 1.9 & all_data_id > ncol(all_data_mat)*0.1)
  cv_features = as.character(all_data$Symbol[cv_features_ind])
  features = cv_features
  
  features_id_cv_df = data.frame(
    Symbol = features,
    id = all_data_id[cv_features_ind],
    cv = all_data_cv[cv_features_ind]
  )
  print(dim(features_id_cv_df))
}
colnames(features_id_cv_df) = c('GeneSymbol', 'Count', 'CV')
write.csv(features_id_cv_df, 'SCLC-75-Count-CV-NMF-445-20220222.csv', row.names = F)
# write.xlsx(features_id_cv_df, "./object_RData/features_id_cv_df.xlsx", row.names = F)

