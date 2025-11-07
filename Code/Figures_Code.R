library(Seurat)
library(sc_utils)
library(ggplot2)
library(dplyr)
library(viridis)
library(stringr)
library(tidyr)
library(SingleR)
library(celldex)
library(BuenColors)
library(scProportionTest)

####### FIGURES HEART TRANSCRIPTOMICS ########

##### UMAP by cell type (2b) ##### 
color_palette <-  c("#1f77b4", "#d62728","#2ca02c",  "#ff7f0e", "darkblue","#bcbd22","lightskyblue2","#9467bd", "#e377c2" ,"#8c564b","#ad494a")

DimPlot(obj2, group.by = "annotation", reduction =  "umap", order=F,  label = F, repel = TRUE, raster=FALSE,  cols = color_palette) + 
  ggtitle(NULL) + NoAxes() + theme(legend.position="bottom",legend.text =   element_text(size=10), legend.key.size = unit(0.2, "cm"), plot.title = element_text(face = "bold", size=12)) +ggtitle("Cell Type Identity (heart)") 

##### UMAP by sample (Sup 3a) ##### 
sample.cols <- c("#E69F00", "darkblue", "#009E73")
DimPlot(obj2, order=T,  label = F, repel = TRUE, raster=FALSE,  cols = sample.cols, split.by = "Sample") + 
  ggtitle("Heart Tissue by Sample") + NoAxes() + theme(legend.text =   element_text(size=10), legend.key.size = unit(0.3, "cm"), plot.title = element_text(face = "bold", size=12))

##### Barplot with relative proportions of cell types by group (2c) ##### 
find_proportions_df <- function(seurat_obj, x, fill) {
  df <- seurat_obj@meta.data %>%
    dplyr::select(x, fill) %>%
    group_by(.data[[x]], .data[[fill]]) %>%
    summarise(n_cells = n()) %>%
    ungroup() %>%
    group_by(.data[[x]]) %>%
    mutate(n_cells_total = sum(n_cells)) %>%
    ungroup() %>%
    mutate(percentage_cells = round(n_cells / n_cells_total * 100, 3))
  df
}

plot_stacked_barplot <- function(df, x, fill, colors) {
  p <- df %>%
    ggplot(aes_string(x, "percentage_cells", fill = fill)) +
    geom_bar(stat="identity")+
    ggtitle("Cell type counts (Relative)") +
    labs(x = "Group", y = "Cells (%)", fill = "") + 
    scale_fill_manual(values = colors) +
    theme(
      axis.title.y = element_text(size = 10),
      axis.text=element_text(size=10),
      axis.text.x=element_text(angle=45, hjust=1,size=8),
      axis.title.x = element_text(vjust = 2, size=10),
      legend.position = "bottom",
      plot.title = element_text(face = "bold", size=11),
    ) 
  p
}

proportions_df_cells <- find_proportions_df(
  obj2,
  x = "Sample",
  fill = "annotation"
)


sample_levels <- (c("Sham", "MI", "MI+APAC"))

proportions_df_cells %>%
  ggplot(aes(x =  factor(Sample, levels=rev(sample_levels)), y = percentage_cells, group = annotation, fill = annotation)) +
  geom_col(position = position_stack(reverse = TRUE))+
  ggtitle("") + 
  coord_flip() +
  scale_fill_manual(values = color_palette) +
  theme( axis.text.x=element_text(angle=45, hjust=1,size=10),
         axis.text.y= element_text(colour="black", size=10),
         axis.title =  element_text(face = "bold"),
         legend.text = element_text(size=10),
         legend.position = "none",
         legend.title = element_text(size=12), 
         panel.background = element_rect(fill='transparent'), 
         panel.border = element_blank(), 
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(), 
         axis.line = element_line(colour = "black"))+
  labs(x="Group", y = "Cell counts (%)", fill = "Cell Type")


##### Cleveland plot (2d) #### 
obj2$orig.ident <- obj2@meta.data$Sample

prop_test <- sc_utils(obj2)

proptest_res <- list()
comp <- c("MI", "MI+APAC")
celltype <- names(table(obj2$annotation))

for (i in 1:length(celltype)){
  res1 <- permutation_test(
    prop_test, cluster_identity = "annotation", 
    sample_1 = "Sham", sample_2 =comp[i] ,
    sample_identity = "orig.ident")
  
  res <- c(res1@results$permutation)  
  proptest_res[[i]] <- res
}
names(proptest_res) <- comp

proptest_res_df <- data.frame(celltypes=rep(proptest_res[[1]]$clusters,2),
                              Pval= c(proptest_res[[1]]$FDR,proptest_res[[2]]$FDR),
                              FC=c(proptest_res[[1]]$obs_log2FD,proptest_res[[2]]$obs_log2FD),
                              Contrast= c(rep("MI vs Sham",8), rep("MI+APAC vs Sham",8)),
                              Color= c(rep("#6a9cd4",8), rep("#f0b72e",8) ))

res2 <- permutation_test(
  prop_test, cluster_identity = "annotation", 
  sample_1 = "MI", sample_2 = "MI+APAC",
  sample_identity = "orig.ident")

res2_ <- c(res2@results$permutation)  

proptest_res_df2 <- data.frame(celltypes=c(res2_$clusters),
                               Pval= c(res2_$FDR),
                               FC=c(res2_$obs_log2FD),
                               Contrast= c(rep("MI+APAC vs MI",8)),
                               Color=c(rep("#f183a1",8)))


proptest_res_df_2 <- rbind(proptest_res_df,proptest_res_df2)

proptest_res_df_2$celltypes <- as.factor(proptest_res_df_2$celltypes)
celltypes_levels <- rev(celltype)

ggplot(proptest_res_df_2, aes(FC, factor(celltypes, levels=celltypes_levels), color=Contrast)) +
  ggtitle("Differences in cell type proportion (Heart)") +
  geom_point(aes(size=Pval),alpha =0.8 ) +
  guides(size = guide_legend(reverse=F))+
  scale_size("P.value (FDR)",range=c(4,1), breaks=c(0.01, 0.05, 0.1)) +
  scale_colour_manual(values=setNames(proptest_res_df_2$Color,proptest_res_df_2$Contrast))+
  xlim(-1, 1)+ geom_vline(xintercept = c(-0.58,0.58), linetype="dotted", color = "black", size=0.5)+
  labs(x = "Log2FC", y = " ", fill = "") + 
  theme_bw()+
  theme(
    legend.position = "right",
    axis.title =element_text(size=10),
    axis.text.x=element_text(angle=45, hjust=1,size=8),
    axis.text.y=element_text(size=8),
    plot.title = element_text(face = "bold", size=11),
    legend.text =   element_text(size=8), legend.key.size = unit(0.4, "cm"), legend.title = element_text(size=9))


##### Volcano Plots (2 e) ####
Celltype_markers <- list() # Create a list to store markers for each Cell Type
for (i in 1:length(celltype)) {
  print(paste("Processing:", celltype[i]))
  cell_type <- subset(x = obj2, subset = annotation == celltype[i])
  Markers <- FindMarkers(cell_type, ident.1 = "MI", ident.2 = "MI+APAC", min.pct = 0.25)
  Markers$Gene <- rownames(Markers)
  Markers <- subset(Markers, !startsWith(Gene, "Rp"))
  Markers <- subset(Markers, !startsWith(Gene, "mt"))
  Markers <- subset(Markers, !startsWith(Gene, "Hb"))
  
  Markers$diffexpressed <- "Not sig."
  Markers$diffexpressed[(Markers$p_val_adj )< 0.05 ] <- "P-value"
  Markers$diffexpressed[abs(Markers$avg_log2FC )> 0.58 & Markers$p_val_adj>0.05] <- "Log(base 2) FC"
  Markers$diffexpressed[(Markers$avg_log2FC )> 0.58 & Markers$p_val_adj<0.05] <- "Adj.P-value & Log2 FC+"
  Markers$diffexpressed[Markers$avg_log2FC < -0.58 & Markers$p_val_adj<0.05] <- "Adj.P-value & Log2 FC-"
  Markers$plog10 <- -log10(Markers$p_val_adj)
  Markers$plog10[Markers$plog10=="Inf"] <- round(max(Markers$plog10[Markers$plog10 != Inf]))
  Celltype_markers[[i]] <- Markers  # Store markers for this patient
}

names(Celltype_markers) <- celltype

EnhancedVolcano(as.data.frame(Celltype_markers$atrial_cm), x="avg_log2FC", y="p_val_adj", lab=Celltype_markers$atrial_cm$Gene, title="Atrial CM (MI+APAC vs MI)",
                     subtitle="", titleLabSize = 10, legendPosition = "none",   # 👈 remove legend
                     subtitleLabSize = 0, captionLabSize = 0, pointSize = 1, labSize = 3,
                     legendLabSize = 8,axisLabSize = 10,cutoffLineWidth = 0.4,hlineWidth=0.1, vlineWidth=0.1, borderWidth=0.2)

EnhancedVolcano(as.data.frame(Celltype_markers$epicardial), x="avg_log2FC", y="p_val_adj", lab=Celltype_markers$epicardial$Gene, title="Epicardial CM (MI+APAC vs MI)",
                     subtitle="", titleLabSize = 10, subtitleLabSize = 0, captionLabSize = 0, pointSize = 1, labSize = 3, legendPosition = "none",
                     legendLabSize = 8,axisLabSize = 10,cutoffLineWidth = 0.4,hlineWidth=0.1, vlineWidth=0.1, borderWidth=0.2)

##### heatmap heart tissue annotation (sup 3b) #### 
Idents(obj2) <- "annotation"
genes <- c("Ttn", "Myh6" ,"Pecam1", "Npr3","Wt1","Pcdh15", "Muc16", "Ebf2", "Pdgfrb", "Hba-a1", "Hdac9", "F13a1","Fabp4", "Dach1")

obj2 <- ScaleData(obj2, features = genes, assay = "RNA")
DefaultAssay(obj2) <- "RNA"
DoHeatmap(obj2, features = genes, assay = "RNA",label = F, group.colors= color_palette) + 
  scale_fill_gradient2( low = (c("#3361A5", "#248AF3", "#14B3FF", "#88CEEF")), mid = "#C1D5DC", high =(c("#EAD397", "#FDB31A" ,"#E42A2A" ,"#A31D1D")), midpoint = 0, guide = "colourbar", aesthetics = "fill")


##### UMAP immune highlight (2f left) ####
obj2@meta.data$grey <- as.character(obj2$annotation)
obj2@meta.data$grey[obj2@meta.data$annotation=="Immune cells"] <- "Immune cells"
obj2@meta.data$grey[!obj2@meta.data$annotation=="Immune cells"] <- " "

table(obj2$grey)
Idents(obj2) <- "grey"

DimPlot(obj2, raster = F, cols=c("darkblue", "ivory4"), label = T) + ggtitle(NULL) + NoAxes() +
  theme( legend.position = "none")

##### UMAP immune subset (2f right) ####
heart_imm_cols <- c("#1f78b4","#b15928", "#33a02c","#e31a1c","#a6cee3", "#6a3d9a","#ff7f00","#b2df8a")

UMADimPlot(heart_imm, group.by="annotation",order=T,  label = F, repel = TRUE, raster=FALSE,  cols = heart_imm_cols) + 
  ggtitle("Immune cells (heart)") + NoAxes() + theme(legend.position="bottom",legend.text =   element_text(size=10), legend.key.size = unit(0.3, "cm"), plot.title = element_text(face = "bold", size=12))

##### Barplot with relative proportions of  immune subset (2g) ##### 
proportions_df_cells <- find_proportions_df(
  heart_imm_2,
  x = "Sample",
  fill = "new_annotation"
)

proportions_df_cells %>%
  ggplot(aes(x =  factor(Sample, levels=sample_levels), y = percentage_cells, group = new_annotation, fill = new_annotation)) +
  geom_col(position = position_stack(reverse = TRUE))+
  ggtitle("") + 
  coord_flip() +
  scale_fill_manual(values = heart_imm_cols) +
  theme( axis.text.x=element_text(angle=45, hjust=1,size=10),
         axis.text.y= element_text(colour="black", size=10),
         axis.title =  element_text(face = "bold"),
         legend.text = element_text(size=10),
         legend.position = "none",
         legend.title = element_text(size=12), 
         panel.background = element_rect(fill='transparent'), 
         panel.border = element_blank(), 
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(), 
         axis.line = element_line(colour = "black"))+
  labs(x="Group", y = "Cell counts (%)", fill = "Cell Type")

##### Heatmap (2h) ##### 
Idents(heart_imm) <- "annotation"

avg_exp <- AverageExpression(
  heart_imm,
  assays = "RNA",
  slot = "data"
)

mat <- avg_exp$RNA
mat <- as.matrix(mat)

markers <- list(
  A = c( "Cd19","Igkc","Iglc1","Iglc2","Iglc3","Ighm"),
  B = c("Cd3e","Cd247","Cd4","Cd8a","Cd8b1","Il7r","Gzma"),
  C = c( "Ccr2","Tlr2","Cd80","Nlrp3","Fcgr1"),
  D = c( "Csf1r","Spi1","Lyz2","Plac8","Cx3cr1"),
  E = c("H2-Ab1","H2-Aa","Itgax","Flt3","Zbtb46","Clec9a"),
  G = c( "Timd4","Lyve1","Folr2","Igf1","Mertk","Adgre1"),
  H = c("Lsamp","Tshz2","Sgip1","Nova2","Ebf2"),
  I = c("S100a8","S100a9","Lcn2","Cxcr2","Retnlg","Csf3r")
  
)

all_markers <- unlist(markers)


# Create a vector matching each marker to its lineage
row_groups <- rep(names(markers), times = sapply(markers, length))
names(row_groups) <- all_markers

# Reorder matrix according to marker list
mat <- mat[all_markers, ]
mat_markers_scaled <- t(scale(t(mat)))
mat_markers_scaled[is.na(mat_markers_scaled)] <- 0

# Define desired order for row splits
row_order <- c(
  "A",
  "B",
  "C",
  "D",
  "E",
  "G",
  "H",
  "I"
)

# Convert to factor with correct order
row_groups <- factor(row_groups, levels = row_order)
# Row annotation: which lineage each gene belongs to
row_ha <- rowAnnotation(
  Cell_types = row_groups,
  col = list(Cell_types = c(
    A = "cornflowerblue",
    B = "chartreuse4",
    C = "coral1",
    D = "tan2",
    E = "magenta4",
    G = "salmon4",
    H = "plum1",
    I = "darkolivegreen"
    
  ))
)

# Color scale
col_fun <- colorRamp2(c(min(mat_markers_scaled), 0, max(mat_markers_scaled)), c("navy", "white", "firebrick3"))

Heatmap(
  mat_markers_scaled[, c(3, 2,8,1,7,4,5,6)],
  name = "Expression",
  col = col_fun,
  cluster_rows = FALSE,
  cluster_columns = F,
  show_row_names = TRUE,
  show_column_names = TRUE,
  row_split = row_groups,
  left_annotation = row_ha,
  row_names_gp = gpar(fontsize = 8),
  column_names_gp = gpar(fontsize = 10),
  column_names_rot = 45,   # angled column names
  heatmap_legend_param = list(
    title = "Expression",
    legend_height = unit(3, "cm"),
    labels_gp = gpar(fontsize = 8),
    title_gp = gpar(fontsize = 9, fontface = "bold")
  )
)


##### myeloid_heart dotplot (2i) #### 
DotPlot(myeloid_heart, genes, scale=T)+
  RotatedAxis()+ scale_color_gradientn(colors = jdb_palette("solar_extra")) + 
  coord_flip()+
  xlab('Markers') +  ylab('Group')+ 
  theme(axis.text.x=element_text(size=8),axis.text.y=element_text(size=8)) +
  scale_size(breaks = c(0,25,50, 75, 100), limits=c(0,100), name = "Percentage") + labs(title="Myeloid lineage") +
  theme(legend.title=element_text(size = 10), legend.position = "bottom",
        axis.title = element_text(size = 10, face="bold"),
        title = element_text(size = 11, face="bold"),
        legend.text =   element_text(size=10), legend.key.size = unit(0.3, "cm"), plot.title = element_text(face = "bold", size=12))


##### cleveland immunes heart (Sup 3 c) ######
heart_imm@meta.data$Sample <-  as.character(heart_imm@meta.data$sample)

heart_imm@meta.data$Sample[heart_imm@meta.data$sample=="Sham heart"] <- "Sham"
heart_imm@meta.data$Sample[heart_imm@meta.data$sample=="MI heart"] <- "MI"
heart_imm@meta.data$Sample[heart_imm@meta.data$sample=="APAC MI heart"] <- "MI+APAC"

Idents(heart_imm) <- "Sample"
order_samples <- c( "Sham", "MI", "MI+APAC")
# Relevel the 'Sample' column to force the order
heart_imm@meta.data$Sample <- factor(heart_imm@meta.data$Sample, levels = order_samples)

# Set the active.ident with the correct order of levels
heart_imm@active.ident <- factor(heart_imm@active.ident, levels = order_samples)
levels(heart_imm) <- order_samples
heart_imm@active.ident <- factor(x = heart_imm@active.ident, levels = order_samples)



heart_imm$orig.ident <- heart_imm@meta.data$Sample
table(heart_imm$orig.ident)
(table(heart_imm$new_annotation))

prop_test <- sc_utils(heart_imm)


proptest_res <- list()
comp <- c("MI", "MI+APAC")
celltype <- names(table(heart_imm$new_annotation))

for (i in 1:length(celltype)){
  res1 <- permutation_test(
    prop_test, cluster_identity = "new_annotation", 
    sample_1 = "Sham", sample_2 =  comp[i],
    sample_identity = "orig.ident")
  
  res <- c(res1@results$permutation)  
  proptest_res[[i]] <- res
}
names(proptest_res) <- comp

proptest_res_df <- data.frame(celltypes=rep(proptest_res[[1]]$clusters,2),
                              Pval= c(proptest_res[[1]]$FDR,proptest_res[[2]]$FDR),
                              FC=c(proptest_res[[1]]$obs_log2FD,proptest_res[[2]]$obs_log2FD),
                              Contrast= c(rep("MI vs Sham",8), rep("MI+APAC vs Sham",8)),
                              Color= c(rep("#6a9cd4",8), rep("#f0b72e",8) ))

res2 <- permutation_test(
  prop_test, cluster_identity = "new_annotation", 
  sample_1 = "MI", sample_2 = "MI+APAC",
  sample_identity = "orig.ident")

res2_ <- c(res2@results$permutation)  

proptest_res_df2 <- data.frame(celltypes=c(res2_$clusters),
                               Pval= c(res2_$FDR),
                               FC=c(res2_$obs_log2FD),
                               Contrast= c(rep("MI+APAC vs MI",8)),
                               Color=c(rep("#f183a1",8)))


proptest_res_df_2 <- rbind(proptest_res_df,proptest_res_df2)

celltype_immune <- rev(res2_$clusters)

celltypes_levels <- (celltype_immune)

ggplot(proptest_res_df_2, aes(FC, factor(celltypes, levels=celltypes_levels), color=Contrast)) +
  ggtitle("Differences in cell type proportion (Heart Immune Cells)") +
  geom_point(aes(size=Pval),alpha =0.8 ) +
  guides(size = guide_legend(reverse=F))+
  scale_size("P.value (FDR)",range=c(4,1), breaks=c(0.01, 0.05, 0.1)) +
  scale_colour_manual(values=setNames(proptest_res_df_2$Color,proptest_res_df_2$Contrast))+
  xlim(-3, 3)+ geom_vline(xintercept = c(-0.58,0.58), linetype="dotted", color = "black", size=0.5)+
  labs(x = "Log2FC", y = " ", fill = "") + 
  theme(
    legend.position = "right",
    axis.title =element_text(size=10),
    axis.text.x=element_text(angle=45, hjust=1,size=8),
    axis.text.y=element_text(size=8),
    plot.title = element_text(face = "bold", size=11),
    legend.text =   element_text(size=8), legend.key.size = unit(0.4, "cm"), legend.title = element_text(size=9))+
  
  theme_classic()



####### FIGURES PBMCs TRANSCRIPTOMICS ########
##### UMAP samples (3b) ######
sample.cols <- c("#E69F00", "darkblue", "#009E73")

DimPlot(mouse.pbmc, reduction = "umap.rpca",group.by = "Sample", order=T,  label = F, repel = TRUE, raster=FALSE,  cols = sample.cols) + 
  ggtitle("scRNA-seq (PBMCs)") + NoAxes() + theme(legend.text =   element_text(size=10), legend.key.size = unit(0.5, "cm"), plot.title = element_text(face = "bold", size=12))

##### Dotplot Identity Markers (3c) ######
cell_type_markers <- c(top_expressed_markers$`B cells`, top_expressed_markers$Monocytes, top_expressed_markers$Neutrophils,
                       top_expressed_markers$`NK cells`, top_expressed_markers$`T cells`)
DotPlot(mouse.pbmc,cell_type_markers  , scale=T)+RotatedAxis()+ scale_color_gradientn(colors = jdb_palette("solar_extra")) + 
  xlab('Genes') +  ylab('Cell Type')+ 
  theme(axis.text.x=element_text(size=7),axis.text.y=element_text(size=7)) +
  scale_size(breaks = c(0,25,50, 75, 100), limits=c(0,100), name = "Percentage") + labs(title="RNA-Seq Markers") +
  theme(legend.text = element_text(size = 9),  
        legend.title=element_text(size = 11),
        axis.text.x = element_text(size = 9), # Adjust x-axis tick text size
        axis.text.y = element_text(size = 9))

##### UMAP PBMCs annotation (3d) ######
pbmc.cell.cols <- c("#7b7b7b", "#f183a1",  "#6a9cd4", "#08355a", "#866b1c")
DimPlot(mouse.pbmc, reduction = "umap.rpca", group.by = "annotation", order=T, label = F, repel = TRUE, raster=FALSE,  cols = pbmc.cell.cols) + 
  ggtitle("Cell Type Identity") + NoAxes() + theme(legend.text =   element_text(size=10), legend.key.size = unit(0.5, "cm"), plot.title = element_text(face = "bold", size=12))

##### Barplot with relative proportions of cell types by group (3e) ##### 
mouse.pbmc@meta.data$Sample <-  as.character(mouse.pbmc@meta.data$sample)

mouse.pbmc@meta.data$Sample[mouse.pbmc@meta.data$sample=="Sham PBMCs"] <- "Sham"
mouse.pbmc@meta.data$Sample[mouse.pbmc@meta.data$sample=="MI PBMCs"] <- "MI"
mouse.pbmc@meta.data$Sample[mouse.pbmc@meta.data$sample=="APAC MI PBMCs"] <- "MI+APAC"

Idents(mouse.pbmc) <- "Sample"
order_samples <- c( "Sham", "MI", "MI+APAC")
levels(mouse.pbmc) <- order_samples
mouse.pbmc@active.ident <- factor(x = mouse.pbmc@active.ident, levels = order_samples)

pmbc.integrated.data <- mouse.pbmc@meta.data %>% group_by(Sample, annotation) %>% 
  summarise(Nb = n()) %>%
  mutate(C = sum(Nb)) %>%
  mutate(percent = Nb/C*100)

pmbc.integrated.data$Sample <- factor(pmbc.integrated.data$Sample, levels=c( "Sham", "MI", "MI APAC"))

v_factor_levels <-  c("B cells", "Monocytes", "Neutrophils", "NK cells","T cells")

ggplot(pmbc.integrated.data, aes(x = factor(Sample, levels=sample_levels), y = percent, fill = factor(annotation, levels = v_factor_levels)))+
  geom_bar(stat = "identity") + 
  coord_flip()+
  geom_col(position = position_stack(reverse = TRUE))+ scale_fill_manual(values=pbmc.cell.cols)+
  theme( axis.text.x=element_text(angle=45, hjust=1,size=10),
         axis.text.y= element_text(colour="black", size=10),
         axis.title =  element_text(face = "bold"),
         legend.text = element_text(size=10),
         legend.position = "none",
         legend.title = element_text(size=12), 
         panel.background = element_rect(fill='transparent'), 
         panel.border = element_blank(), 
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(), 
         axis.line = element_line(colour = "black"))+
  labs(x="Group", y = "Cell counts (%)", fill = "Cell Type")

##### Heatmap (3f) ######
Markers_all <- FindAllMarkers(mouse.pbmc,  min.pct = 0.25, logfc.threshold=0.58)
Markers_all$genes <- rownames(Markers_all)
markers_use <- unique(Markers_all$genes)

mouse.pbmc <- ScaleData(mouse.pbmc, features = markers_use)

# Extract scaled expression matrix
mouse.pbmc <- ScaleData(mouse.pbmc, features = markers_use)
expr_mat <- as.data.frame(GetAssayData(mouse.pbmc, slot = "scale.data"))
expr_mat <- expr_mat[markers_use,]

expr_mat_clean <- expr_mat[complete.cases(expr_mat), ]

# Get row order from hierarchical clustering
row_order <- rownames(expr_mat_clean)[order.dendrogram(as.dendrogram(hclust(dist(expr_mat_clean))))]

DoHeatmap(mouse.pbmc, features = row_order,label = T , group.colors= sample.cols) + 
  scale_fill_gradient2( low = (c("#3361A5", "#248AF3", "#14B3FF", "#88CEEF")), mid = "#C1D5DC", high =(c("#EAD397", "#FDB31A" ,"#E42A2A" ,"#A31D1D")), midpoint  = 0, guide = "colourbar", aesthetics = "fill")


##### Dotplot (3g) ######
Idents(mouse.pbmc) <- "Sample"
table(Idents(mouse.pbmc))

highlight_genes <- c(
  
  "Cebpb", "S100a11", "S100a8", "S100a9","Tspo","Litaf", "Apoe", 
  "Cd3e", "Cd3g", "Bcl11a", "Mef2c","Card11"
  
)

DotPlot(mouse.pbmc,features=highlight_genes  , scale=T)+ RotatedAxis()+ scale_color_gradientn(colors = jdb_palette("solar_extra")) + 
  xlab('Genes') +  ylab('Group')+ coord_flip()+
  theme(axis.text.x=element_text(size=7),axis.text.y=element_text(size=7)) +
  scale_size(breaks = c(0,25,50, 75, 100), limits=c(0,100), name = "Percentage") + labs(title="") +
  theme(legend.text = element_text(size = 9),  
        legend.title=element_text(size = 10, face="bold"),
        axis.text.x = element_text(size = 9), # Adjust x-axis tick text size
        axis.text.y = element_text(size = 9),
        legend.position = "bottom",
        plot.title = element_text(face = "bold", size=11),
        axis.title = element_text(size=10, face="bold"),legend.key.size = unit(0.4, "cm"))


##### Barplot (3h) ######
Idents(mouse.pbmc) <- "sample"

order_cells <- names(table(mouse.pbmc@meta.data$annotation))
samples <- names(table(mouse.pbmc@meta.data$sample))[-1]

Celltype_markers <- list() # Create a list to store markers for each Cell Type
for (i in 1:length(order_cells)) {
  markers_for_contrast <- list()   # Create a list to store markers for each contrast
  for (j in 1:length(samples)) {
    print(paste("Processing:", order_cells[i], "for", samples[j]))
    cell_type <- subset(x = mouse.pbmc, subset = annotation == order_cells[i])
    Markers <- FindMarkers(cell_type, ident.1 = samples[j], ident.2 = "Sham PBMCs", min.pct = 0.25)
    Markers$Gene <- rownames(Markers)
    Markers <- subset(Markers, !startsWith(Gene, "Rp"))
    Markers <- subset(Markers, !startsWith(Gene, "Hb"))
    Markers <- subset(Markers, !startsWith(Gene, "mt-"))
    Markers$diffexpressed <- "Not sig."
    Markers$diffexpressed[(Markers$p_val_adj )< 0.05 ] <- "Adj.P-value"
    Markers$diffexpressed[abs(Markers$avg_log2FC )>= 0.58 & Markers$p_val_adj>0.05] <- "Log2FC"
    Markers$diffexpressed[(Markers$avg_log2FC )>= 0.58 & Markers$p_val_adj<0.05] <- "Adj.P-value & Log2FC +"
    Markers$diffexpressed[Markers$avg_log2FC <= -0.58 & Markers$p_val_adj<0.05] <- "Adj.P-value & Log2FC -"
    Markers$plog10 <- -log10(Markers$p_val_adj)
    Markers$plog10[Markers$plog10=="Inf"] <- round(max(Markers$plog10[Markers$plog10 != Inf]))
    markers_for_contrast[[j]] <- Markers  # Store markers for this patient
  }
  Celltype_markers[[i]] <- markers_for_contrast  # Store markers for this cell type
}
names(Celltype_markers) <- order_cells

samples <- names(table(mouse.pbmc@meta.data$sample))[-c(1,2)]
Celltype_markers2 <- list() # Create a list to store markers for each Cell Type
for (i in 1:length(order_cells)) {
  markers_for_contrast <- list()   # Create a list to store markers for each contrast
  for (j in 1:length(samples)) {
    print(paste("Processing:", order_cells[i], "for", samples[j]))
    cell_type <- subset(x = mouse.pbmc, subset = annotation == order_cells[i])
    Markers <- FindMarkers(cell_type, ident.1 = samples[j], ident.2 = "MI PBMCs", min.pct = 0.25)
    Markers$Gene <- rownames(Markers)
    Markers <- subset(Markers, !startsWith(Gene, "Rp"))
    Markers <- subset(Markers, !startsWith(Gene, "Hb"))
    Markers <- subset(Markers, !startsWith(Gene, "mt-"))
    Markers$diffexpressed <- "Not sig."
    Markers$diffexpressed[(Markers$p_val_adj )< 0.05 ] <- "Adj.P-value"
    Markers$diffexpressed[abs(Markers$avg_log2FC )>= 0.58 & Markers$p_val_adj>0.05] <- "Log2FC"
    Markers$diffexpressed[(Markers$avg_log2FC )>= 0.58 & Markers$p_val_adj<0.05] <- "Adj.P-value & Log2FC +"
    Markers$diffexpressed[Markers$avg_log2FC <= -0.58 & Markers$p_val_adj<0.05] <- "Adj.P-value & Log2FC -"
    Markers$plog10 <- -log10(Markers$p_val_adj)
    Markers$plog10[Markers$plog10=="Inf"] <- round(max(Markers$plog10[Markers$plog10 != Inf]))
    markers_for_contrast[[j]] <- Markers  # Store markers for this patient
  }
  Celltype_markers2[[i]] <- markers_for_contrast  # Store markers for this cell type
}
names(Celltype_markers2) <- order_cells


DEA_numbers_celltype <- list()
for(i in 1:length(order_cells)){
  Up <- list(a=length(unique(Celltype_markers[[i]][[1]]$Gene[Celltype_markers[[i]][[1]]$diffexpressed=="Adj.P-value & Log2FC +"])),
             b=length(unique(Celltype_markers[[i]][[2]]$Gene[Celltype_markers[[i]][[2]]$diffexpressed=="Adj.P-value & Log2FC +"])),
             c=length(unique(Celltype_markers2[[i]][[1]]$Gene[Celltype_markers2[[i]][[1]]$diffexpressed=="Adj.P-value & Log2FC +"])))
  Down <- list(a=length(unique(Celltype_markers[[i]][[1]]$Gene[Celltype_markers[[i]][[1]]$diffexpressed=="Adj.P-value & Log2FC -"])),
               b=length(unique(Celltype_markers[[i]][[2]]$Gene[Celltype_markers[[i]][[2]]$diffexpressed=="Adj.P-value & Log2FC -"])),
               c=length(unique(Celltype_markers2[[i]][[1]]$Gene[Celltype_markers2[[i]][[1]]$diffexpressed=="Adj.P-value & Log2FC -"])))
  DEA <- list(Up, Down)
  DEA_numbers_celltype[[i]] <- DEA
}


names(DEA_numbers_celltype) <- order_cells
df <- as.data.frame(do.call(rbind, (DEA_numbers_celltype)))
df <- as.data.frame(do.call(rbind, as.vector(DEA_numbers_celltype)))


markers_celltypes <- data.frame(UP= df[2,], 
                                DOWN= df[3,])


markers_celltypes$celltype <- rep(order_cells,each=3)
markers_celltypes$contrast <- rep(c("MI vs SHAM", "MI APAC vs SHAM", "APAC MI vs MI"), 5)

markers_celltypes$contrast <- factor(markers_celltypes$contrast , levels=c("MI vs SHAM", "MI APAC vs SHAM", "APAC MI vs MI"))


# Reshape to long so UP/DOWN become a single variable
df_long <- markers_celltypes %>%
  pivot_longer(cols = c(UP, DOWN), names_to = "direction", values_to = "count")


# Plot
 ggplot(df_long, aes(x = celltype, y = count, fill = direction)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~ contrast) +
  theme_minimal() +
  labs(x = "Cell type", y = "Number of DEGs", fill = "Regulation") +
  scale_fill_manual(values =  c("darkslateblue" , "red" ))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size=8),
        axis.title =element_text(size=10, face="bold"),
        axis.text.y=element_text(size=8),
        plot.title = element_text(face = "bold", size=11),
        legend.text =   element_text(size=8), legend.key.size = unit(0.4, "cm"), 
        legend.title = element_text(size=9, face="bold"), legend.position = "bottom")



##### Dotplot (4a) ######
genes <- c("Il1b",
           "Nlrp3",
           "Il18",
           "Tlr4",
           "C5ar1",
           "Itgam",
           "Selplg",
           "Il1rn",
           "Il18bp",
           "Socs3")

DotPlot(myeloids,features=genes , scale=T)+ RotatedAxis()+ scale_color_gradientn(colors = jdb_palette("solar_extra")) + 
  xlab('Genes') +  ylab('Cell Type')+ coord_flip()+
  theme(axis.text.x=element_text(size=7),axis.text.y=element_text(size=7)) +
  scale_size(breaks = c(0,25,50, 75, 100), limits=c(0,100), name = "Percentage") + labs(title="Inflammasome priming activity (PBMCs myeloids cells)") +
  theme(legend.text = element_text(size = 9),  
        legend.title=element_text(size = 10, face="bold"),
        axis.text.x = element_text(size = 9), # Adjust x-axis tick text size
        axis.text.y = element_text(size = 9),
        legend.position = "bottom",
        plot.title = element_text(face = "bold", size=11),
        axis.title = element_text(size=10, face="bold"),legend.key.size = unit(0.4, "cm"))


##### GSEA dotplot (4b) ######
ggplot(plot_df_myeloids,
       aes(x = cluster, y = Description, size = -log10(p.adjust), color = NES)) +
  geom_point(na.rm = FALSE) +
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, na.value = "grey80") +
  theme_minimal() +
  labs(
    title = "Top 10 GO Terms myeloids",
    x = "Cluster",
    y = "GO Term",
    size = "-log10(adj p-value)",
    color = "NES"
  ) +
  theme(axis.text.y = element_text(size = 8))

##### (Sup 4a) ######
plotScoreHeatmap(pred, labels.use = c("B cells"  ,"Neutrophils",
                                      "Monocytes" ,  "NK cells" , "T cells"), cluster_cols = F)
##### Cleveland Plot (Sup 4b) ######

mouse.pbmc$orig.ident <- mouse.pbmc@meta.data$Sample
table(mouse.pbmc$orig.ident)
(table(mouse.pbmc$annotation))

prop_test <- sc_utils(mouse.pbmc)


proptest_res <- list()
comp <- c("MI", "MI APAC")
celltype <- names(table(mouse.pbmc$annotation))

for (i in 1:length(celltype)){
  res1 <- permutation_test(
    prop_test, cluster_identity = "annotation", 
    sample_1 = "Sham", sample_2 = comp[i],
    sample_identity = "orig.ident")
  
  res <- c(res1@results$permutation)  
  proptest_res[[i]] <- res
}
names(proptest_res) <- comp

proptest_res_df <- data.frame(celltypes=rep(proptest_res[[1]]$clusters,2),
                              Pval= c(proptest_res[[1]]$FDR,proptest_res[[2]]$FDR),
                              FC=c(proptest_res[[1]]$obs_log2FD,proptest_res[[2]]$obs_log2FD),
                              Contrast= c(rep("MI vs Sham",5), rep("MI APAC vs Sham",5)),
                              Color= c(rep("#6a9cd4",5), rep("#f0b72e",5) ))

res2 <- permutation_test(
  prop_test, cluster_identity = "annotation", 
  sample_1 = "MI", sample_2 = "MI APAC",
  sample_identity = "orig.ident")

res2_ <- c(res2@results$permutation)  

proptest_res_df2 <- data.frame(celltypes=c(res2_$clusters),
                               Pval= c(res2_$FDR),
                               FC=c(res2_$obs_log2FD),
                               Contrast= c(rep("MI vs APAC",5)),
                               Color=c(rep("#f183a1",5)))


proptest_res_df_2 <- rbind(proptest_res_df,proptest_res_df2)

str(proptest_res_df_2)
proptest_res_df_2$celltypes <- as.factor(proptest_res_df_2$celltypes)


ggplot(proptest_res_df_2, aes(FC, celltypes, color=Contrast)) +
  ggtitle("Differences in cell type proportion") +
  geom_point(aes(size=Pval),alpha =0.8 ) +
  guides(size = guide_legend(reverse=F))+
  scale_size("P.value (FDR)",range=c(4,1), breaks=c(0.01, 0.05, 0.1)) +
  scale_colour_manual(values=setNames(proptest_res_df_2$Color,proptest_res_df_2$Contrast))+
  xlim(-8, 8)+ geom_vline(xintercept = c(-0.58,0.58), linetype="dotted", color = "black", size=0.5)+
  labs(x = "Log2FC", y = " ", fill = "") + 
  theme(
    axis.title =element_text(size=10),
    axis.text.x=element_text(angle=45, hjust=1,size=8),
    axis.text.y=element_text(size=8),
    plot.title = element_text(face = "bold", size=11),
    legend.text =   element_text(size=8), legend.key.size = unit(0.4, "cm"), legend.title = element_text(size=9))

##### (Sup 4c) ######
ggplot(plot_df_monocytes,
                               aes(x = cluster, y = Description, size = -log10(p.adjust), color = NES)) +
  geom_point(na.rm = FALSE) +
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, na.value = "grey80") +
  theme_minimal() +
  labs(
    title = "Top 10 GO Terms Monocytes",
    x = "Cluster",
    y = "GO Term",
    size = "-log10(adj p-value)",
    color = "NES"
  ) +
  theme(axis.text.y = element_text(size = 8))

##### (Sup 4d) ######
ggplot(plot_df_neutrophils,
                                 aes(x = cluster, y = Description, size = -log10(p.adjust), color = NES)) +
  geom_point(na.rm = FALSE) +
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, na.value = "grey80") +
  theme_minimal() +
  labs(
    title = "Top 10 GO Terms Neutrophils",
    x = "Cluster",
    y = "GO Term",
    size = "-log10(adj p-value)",
    color = "NES"
  ) +
  theme(axis.text.y = element_text(size = 8))

