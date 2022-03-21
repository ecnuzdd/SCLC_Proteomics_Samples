if(T){
  #### library ####
  library('readxl')
  library('xlsx')
  library('stringr')
  library('stringi')
  library('seqinr')
  library("vioplot")
  library('ggpubr')
  library('ggsci')
  library('viridis')
  library('pheatmap')
  library('reshape2')
  library('survminer')
  library('ggbeeswarm')
  library('RColorBrewer')
  library('scales')
  
  #### directory ####
  # work 
  BASE_DIR = "/home/ddzhan/SCLC/20211022-plot"
  setwd(BASE_DIR)
  
  # input
  dat_path = normalizePath(file.path(BASE_DIR, 'Fig3_heatmap.xlsx'))
  
  
}





if(T){
  # annotation
  edata = read_excel(dat_path, sheet = 'edata')
  adata = read_excel(dat_path, sheet = 'adata')
  cdata = read_excel(dat_path, sheet = 'cdata')
  
  cdata$Subtype = factor(cdata$Subtype, levels = c('S-I', 'S-II', 'S-III'))
  
  # annotation
  annotation_col = data.frame(
    Subtype = cdata$Subtype
  )
  rownames(annotation_col) = cdata$ExpID
  
  ann_colors = list(
    Subtype = c("S-I"="#49A74B", "S-II"="#3279AE", "S-III"="#D42322")
  )
  
  ph_mat = as.matrix(edata[,-1])
  rownames(ph_mat) = edata$GeneSymbol
  
  ph_mat_scale = t(scale(t(ph_mat)))
  ph_mat_scale.minv = min(ph_mat_scale)
  ph_mat_scale.maxv = max(ph_mat_scale)
  
  breaks = seq(-9, 9, 0.05)
  breaks_n = length(breaks)
  breaks_cols = colorRampPalette(
    c('#11264f', '#145b7d', '#009ad6', 'white', '#FF6600', 'red', 'firebrick')
  )(breaks_n)
  
  
  ph_gap_col = cumsum(table(cdata$Subtype))
  # ph_gap_row = cumsum(table(adata$Subtype))
  p1 <- pheatmap(
    ph_mat,
    scale = 'row',
    cluster_cols = F,cluster_rows = F,
    annotation_col = annotation_col,
    # annotation_row = annotation_row,
    annotation_colors = ann_colors,
    # col = colorRampPalette(c( "blue2","blue1","blue","black","red","red1","red2"))(100),
    col = breaks_cols,
    breaks = breaks,
    # clustering_distance_rows = 'euclidean',
    # main = paste0('Features = ', nrow(nmfMatrix)),
    # cellwidth = 6,
    # cellheight = 3,
    # clustering_method = 'ward.D2',
    # border_color = "grey60",
    # fontsize_row = 2,
    show_rownames = F,
    show_colnames = F,
    # gaps_row = ph_gap_row,
    gaps_col = ph_gap_col
  )
  print(p1)
}

