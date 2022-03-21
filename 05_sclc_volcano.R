# setwd
BASE_DIR = normalizePath("/home/ddzhan/SCLC/20211022-plot")
setwd(BASE_DIR)

# loading packages
library(xlsx)
library(readxl)
library(stringi)
library(stringr)
library(ggplot2)
library(ggpubr)
library(ggsci)
library(DescTools)
library(RColorBrewer)
library(ggrepel)




# anno data
anno_fpath = normalizePath(
  file.path(BASE_DIR, 'sclc_wh_fot_20220214_zdd_v1.0.xlsx')
)
adata_volcano = read_excel(anno_fpath, sheet = 'adata_volcano')
adata_classfication = read_excel(anno_fpath, sheet = 'adata_classfication')
# adata_volcano = adata_volcano[which(!duplicated(adata_volcano$MeanFot)),]

df = adata_volcano
df$Label = 'No Signifcant'
df$Label[which(df$Pvalue<0.05 & df$FoldChange<0.5)] = 'Down-regulated'
df$Label[which(df$Pvalue<0.05 & df$FoldChange>2)] = 'Up-regulated'
df$Label[which(df$Pvalue<0.05 & df$FoldChange>10)] = 'Specifically overexpressed'
df$Label[match(adata_classfication$GeneSymbol, df$GeneSymbol)] = 'Druggable'
df$Label = factor(
  df$Label, levels = c(
    'No Signifcant', 'Down-regulated', 'Up-regulated', 'Specifically overexpressed', 'Druggable'
    )
)

classification_colors = pal_npg("nrc", alpha=0.7)(10)
show_col(classification_colors)
classification_colors_1 = pal_npg("nrc", alpha=1)(10)
show_col(classification_colors_1)

p <- ggplot(df, mapping=aes(x=log2(FoldChange), y=-log10(Pvalue), colour=Label))
p <- p + geom_point(size=2) 
p <- p + scale_color_manual(
  values = c(
    'No Signifcant' = 'grey60',
    'Down-regulated' = 'grey60', #classification_colors[2],
    'Up-regulated' = '#fdc9c9', #classification_colors[5],
    'Specifically overexpressed' = '#fa949e', #classification_colors[1],
    'Druggable' = 'firebrick' # classification_colors_1[8] 
  )
)
print(p)
chromatin_proteins = adata_classfication$GeneSymbol[adata_classfication$FunClassification=='Chromatin']
df1 <- df[df$Label=='Druggable',]
p <- p + geom_point(
  data=df1, 
  mapping = aes(x=log2(FoldChange), y=-log10(Pvalue)),
  color = 'firebrick',
  size = 2
) 

p <- p + geom_vline(
  xintercept = c(log2(0.5), log2(2), log2(10)), 
  linetype="dashed", 
  color = c("grey", "grey", classification_colors[1])
) + geom_hline(
  yintercept = c(-log10(0.05)), 
  linetype="dashed", 
  color = "grey"
)
p <- p + scale_x_continuous(
  breaks = c(-5, -3, -1, 0, 1, 3, 5, 10, 15, 20)
) + theme_classic2()

p <- p + labs(
  x = expression(paste(Log["2"], '(Fold Change)')),
  y = expression(paste('-', Log["10"], '(P-value)'))
)
print(p)


p + geom_text_repel(
  data = df1, 
  mapping = aes(x=log2(FoldChange), y=-log10(Pvalue), label=GeneSymbol),
  color = "black",
  force = 10,
  point.padding = 0.5,
  size = 3
)





