
#Thymomstudy


#(1)_Base activation
{
  rainbow.colors<-palette(rainbow(25))
  library(Polychrome)
  library("viridis")
  color_level1 <-  as.character(c("#0F8299","#3D0F99","#FFA07A","#FF8800","#00E2FC","#967ACC","#54990F","#990F26","#CC7A88", "#666666"))
  
  color_TC <-  as.character(c("#0F8299"))
  color_levelTC <-  as.character(c("#39a0b3","#256875","#96bfe9","#74abe2","#5296da","#3182d3","#266eb6","#1f5994","#00942d","#147d39","#8aad34","#599717","#a7dca4","#009382"))
  color_TEC <-  as.character(c("#FF8800"))
  color_levelTEC <-  as.character(c("#ff006d","#9e0142","#f46d43","#fdae64","#fdc629","#fee08b","#599717","#abdda4","#66c2a5","#3298bd","#5e4fa2","white","#de0a26","#86cbd8","#c30010"))
  color_levelTETTEC <-  as.character(c("#01befe","#fdc629","#ff7d00","#ff006d","#adff02","#8f00ff","#2b2b3a"))
  color_BC <-  as.character(c("#3D0F99"))
  color_levelBC <-  as.character(c("#ce389c","#fc94af","#7434a4","#c19adf","#8b008b","#a9957b","#601a3e","#3D0F99","#8530d1"))
  color_MAC_Mono <-  as.character(c("#967ACC"))
  color_levelMAC_Mono <-  as.character(c("#a9957b","#BA68CB","#00EC9D","#a4de02","#4c9a2a","#00E2FC","#96bfe9","#8530d1","#3D0F99"))
  color_DC <-  as.character(c("#FF8800","#FF8800","#FF8800","#FF8800","#FF8800","#FF8800","#FF8800","#FF8800","#FF8800"))
  color_levelDC <-  as.character(c("#01befe","#fdc629","#ff7d00","#ff006d","#adff02","#8f00ff","#2b2b3a"))
  color_FB <-  as.character(c("#990F26"))
  color_levelFB <-  as.character(c("#9c6137","#e47200","red","yellow","#96bfe9","#990F26","#FFA07A","#adff02","#ff1493","#fc94af","#8b008b"))
  color_EC <-  as.character(c("#54990F"))
  color_levelEC <-  as.character(c("#9c6137","#e47200","#54990F","#a4de02","#22a384ff","#a9957b","#e69b00"))
  color_levelVSMC_Peri <- as.character(c("#cdcdcd","#9a9a9a","#808080","#778899","#2f4f4f","#4c5762","#404040"))
  color_levelPC <- as.character(c("#2a7b8eff","#22a384ff","#7ad151ff","#a4de02","#4c9a2a","#c19adf"))
  color_levelallcells <- as.character(c("#009dc4","#007e9d","#005e76","#003f4e","#000137","#a9957b","#FFA07A",
                                        "#BA68CB","#00EC9D","#00E2FC","#8530d1",
                                        "#990F26","#dc1c13","#ff1493","#CC7A88",
                                        "#54990F","#a4de02",
                                        "#FF8800",
                                        "#ce389c","#fc94af","#7434a4","#c19adf","#8b008b","#a9957b","#601a3e"))
  
  color_levelallcells_cch <- as.character(c("#009dc4","#007e9d","#005e76","#003f4e","#000137","#a9957b","#FFA07A",
                                            "#a9957b","#fc94af","#7434a4", "#c19adf","#8b008b",
                                            "#BA68CB","#00E2FC","#8530d1",
                                            "#FF8800",
                                            "#dc1c13","#ff1493","#CC7A88",
                                            "#54990F","#a4de02"
  ))
  
  color_levelallcells_cct <- as.character(c("#009dc4","#007e9d","#005e76","#003f4e","#000137","#a9957b","#FFA07A",
                                            "#a9957b","#fc94af","#7434a4", "#c19adf","#8b008b","#601a3e",
                                            "#BA68CB","#00E2FC","#8530d1",
                                            "#FF8800",
                                            "#990F26","#dc1c13","#ff1493","#CC7A88",
                                            "#54990F","#a4de02"))
  
  color_levelallcells_ccc <- as.character(c("#009dc4","#007e9d","#005e76","#003f4e","#000137","#a9957b","#FFA07A",
                                            "#7434a4", "#8b008b","#601a3e",
                                            "#BA68CB","#00E2FC","#00EC9D",
                                            "#FF8800",
                                            "#990F26","#ff1493","#CC7A88",
                                            "#54990F","#a4de02","#9c6137")) 
  
  color_maintissue <- as.character(c("#e83845","#ffce30","#288ba8"))
  
  color_maintissue.sample <- as.character(c("#e83845","#e83845","#e83845","#e83845","#ffce30","#ffce30","#ffce30","#ffce30","#288ba8","#288ba8"))
  
  color_UMAP <- as.character(c("#54990F", "#78B33E","#03AC13",
                               "#990F26", "#B33E52", "#CC7A88", "#E6B8BF", 
                               "#F9A602","#99600F", 
                               "#0F8299", "#3E9FB3", 
                               "#967ACC", "#B3823E",
                               "#666666", "#333333", "#999999", "#CCCCCC","#3D0F99","#3D0F99","#3D0F99","#3D0F99","#3D0F99","#3D0F99","#3D0F99"))
  
  color_UMAP2 <- as.character(c("#1F78C8", "#0000FF","#a6cee3",  "#ff0000", "#33a02c", 
                                "#6A33C2", "#ff7f00", "#FFD700", 
                                "#36648B", 
                                "#00E2E5", "#00FF00", "#778B00", "#BEBE00", 
                                "#8B3B00", "#A52A3C", "#FB6496", "#b2df8a", "#CAB2D6", 
                                "#FDBF6F", "#999999", "#EEE685", "#C8308C", 
                                "#FF83FA", "#C814FA"))
  color_th <- as.character(c("#009dc4","#008db0","#007e9d","#006e89","#005e76","#004f62","#003f4e","#015c5d","#016868","#F9A602","#c38452","#e3c099","#046307","#990F26","#E6B8BF","#888888"))
  color_DEG <- as.character(c("#54990F","#C51D34"))
  single_col <- as.character("#02025f")
  donor_col <- viridis(53)
  #yayon marker genes
  yayon_BC=c("MKI67","JCHAIN","PRDM1","XBP1","CLECL1","FCRL4","TNFRSF13B","IGHG1","TCL1A","MS4A1","IGHD","IGHM","SPIB","IDH2","MME","TNFRSF17","CD24","DNTT","CD27","CDC45","VPREB1","CD19","CD79A")
  yayon_DC_MAC=c("TOP2A","MKI67","CLEC4C","LILRA4","AIRE","TNFRSF11B","SYNPO2","EYA4","MT2A","SELENOM","LYPD3","CCL17","KCNN1","IL12B","CST7","IDO1","CD86","HLA-DRA","CD80","CCR7","LAMP3","CD1C","CLEC10A","XCR1","CLEC9A","GPNMB","MMP9","APOC2","TIMD4","SPIC","LILRB5","EGFL7","LYVE1","CD163","CD68")
  yayon_myeloid=c("CD34","SPINK2","PRSS57","PRTN3","AZU1","ELANE","DEFA4","LCN2","ORM1","CD52","S100A8","MS4A6A","CD14","CXCR4","CCR2","IL1B","GYPA","GYPB","HBA1","HBA2")
  yayon_FB=c("RSPO2","CES1","MYF5","PAX7","TNNT1","CXCL13","FDCSP","CCL21","LTBP1","HEY1","SBSPON","HLA-DRA","RGS5","CXCL10","CXCL9","TNFSF10","IL33","CCL19","DHRS3","WIF1","COL13A1","COL9A3","SEMA3C","MFAP5","PI16","FBN1","MKI67","GDF10","ALDH1A2","COLEC11","VIM","PDGFRB","PDGFRA")
  yayon_EC=c("CCL21","TFF3","PROX1","SELE","CCL2","ICAM1","ACKR1","MKI67","PLVAP","RGCC","ELN","SULF1","HEY1","SEMA3G","CXCL12","VWF","CDH5","PECAM1")
  yayon_VSMC=c("TOP2A","MKI67","SULF1","ELN","COL1A1","CCL21","CCL19","KCNJ8","ABCC9","PDGFRB","RGS5","CASQ2","RERGL","MYH11","ACTA2")
  yayon_TC=c("CD34","PECAM1","PCNA","RAG1","RAG2","TRDC","KLRB1","IL2RB","CD4","CD8A","CD8B","SATB1","GNG4","TOX2","KLF2","TIGIT","IL2RA","FOXP3","HSPH1","JUNB")
  yayon_TEC=c("PSMB11","LY75","CCL25","HLA-DRA","TBATA","TP53AIP1","DLL4","DLK2","IGFBP5","IGFBP6","CCN2","CCL2","KRT15","ITGA6","MKI67","EPCAM","ASCL1","CCL21","AIRE","FEZF2","CRIP1","SLPI","IVL","KRT10","CDKN2A","BEX1","NEUROD1","NEUROG1","NEUROD4","PCP4","FOXJ1","CHRNA1","MYOG","TTN","FOXI1","CFTR","POU2F3","PLCB2")
  
  library(Seurat)
  library(dplyr)
  library(magrittr)
  library(ggplot2)
  library(xlsx)
  library(patchwork)
  library(sctransform)
  library(SeuratWrappers)
  library(ggrepel)
  library(CellChat)
  library(patchwork)
  library(tidyr)
  library(clustree)
  library(ggsignif)
  library(limma)
  library(Matrix.utils)
  library(rstatix)
  library(readxl)
  library(Seurat)
  library(decoupleR)
  library(tibble)
  library(patchwork)
  library(ggplot2)
  library(pheatmap)
  library(celda)
  library(decontX)
  
  color_group <-  as.character(c("#bf2c34","#43a5be","#4fb06d","#be398d","#d49137"))
  color_condition <-  as.character(c("#fff100","#ff8c00","#e81123","#ec008c","#00188f","#009e49","#bad80a"))
  color_tissue <-  as.character(c("#fff100","#ff8c00","#cccccc","#e81123","#ec008c","#68217a","#00188f","#009e49","#bad80a"))
  color_classification <-  as.character(c("#fff100","#ff8c00","#e81123","#ec008c","#68217a","#00188f","#00bcf2","#00b294","#009e49","#bad80a"))
  color_maintissue <-  as.character(c("#76d7c4","#16a085","#1e8449"))
  color_sex <-  as.character(c("#34495e","#99a3a4","#bfc9ca"))
  color_development <-  as.character(c("#f4d03f","#f39c12","#ca6f1e","#a04000"))
  
  plot_integrated_clusters = function (srat) { 
    ## take an integrated Seurat object, plot distributions over orig.ident
    library(Seurat)
    library(patchwork)
    library(ggplot2)
    library(reshape2)
    library(RColorBrewer)
    
    
    count_table <- table(srat@meta.data$seurat_clusters, srat@meta.data$orig.ident)
    count_mtx   <- as.data.frame.matrix(count_table)
    count_mtx$cluster <- rownames(count_mtx)
    melt_mtx    <- melt(count_mtx)
    melt_mtx$cluster <- as.factor(melt_mtx$cluster)
    
    cluster_size   <- aggregate(value ~ cluster, data = melt_mtx, FUN = sum)
    
    sorted_labels <- paste(sort(as.integer(levels(cluster_size$cluster)),decreasing = T))
    cluster_size$cluster <- factor(cluster_size$cluster,levels = sorted_labels)
    melt_mtx$cluster <- factor(melt_mtx$cluster,levels = sorted_labels)
    
    colnames(melt_mtx)[2] <- "dataset"
    
    
    p1 <- ggplot(cluster_size, aes(y= cluster,x = value)) + geom_bar(position="dodge", stat="identity",fill = "grey60") + 
      theme_bw() + scale_x_log10() + xlab("Cells per cluster, log10 scale") + ylab("")
    p2 <- ggplot(melt_mtx,aes(x=cluster,y=value,fill=dataset)) + 
      geom_bar(position="fill", stat="identity") + theme_bw() + coord_flip() + 
      scale_fill_brewer(palette = "Set2") +
      ylab("Fraction of cells in each dataset") + xlab("Cluster number") + theme(legend.position="top")
    
    p2 + p1 + plot_layout(widths = c(3,1))
  }
  plot_integrated_clusters_group = function (srat) { 
    ## take an integrated Seurat object, plot distributions over orig.ident
    library(Seurat)
    library(patchwork)
    library(ggplot2)
    library(reshape2)
    library(RColorBrewer)
    
    
    count_table <- table(srat@meta.data$celltype.full, srat@meta.data$group)
    count_mtx   <- as.data.frame.matrix(count_table)
    count_mtx$celltype <- rownames(count_mtx)
    melt_mtx    <- melt(count_mtx)
    melt_mtx$celltype <- as.factor(melt_mtx$celltype)
    
    celltype_size   <- aggregate(value ~ celltype, data = melt_mtx, FUN = sum)
    
    sorted_labels <- levels(srat$celltype.full)
    celltype_size$celltype <- factor(celltype_size$celltype,levels = sorted_labels)
    melt_mtx$celltype <- factor(melt_mtx$celltype,levels = sorted_labels)
    
    colnames(melt_mtx)[2] <- "dataset"
    
    
    p1 <- ggplot(celltype_size, aes(y= celltype,x = value)) + geom_bar(position="dodge", stat="identity",fill = "grey60") + 
      theme_bw() + scale_x_log10() + xlab("Cells per celltype, log10 scale") + ylab("")
    p2 <- ggplot(melt_mtx,aes(x=celltype,y=value,fill=dataset)) + 
      geom_bar(position="fill", stat="identity") + theme_bw() + coord_flip() + 
      scale_fill_manual(values= color_group) +
      ylab("Fraction of cells in each dataset") + xlab("celltype number") + theme(legend.position="top")
    
    p2 + p1 + plot_layout(widths = c(3,1))
  }
  plot_integrated_clusters_Donor_ID = function (srat) { 
    ## take an integrated Seurat object, plot distributions over orig.ident
    library(Seurat)
    library(patchwork)
    library(ggplot2)
    library(reshape2)
    library(RColorBrewer)
    
    
    count_table <- table(srat@meta.data$celltype.full, srat@meta.data$Donor_ID)
    count_mtx   <- as.data.frame.matrix(count_table)
    count_mtx$celltype <- rownames(count_mtx)
    melt_mtx    <- melt(count_mtx)
    melt_mtx$celltype <- as.factor(melt_mtx$celltype)
    
    celltype_size   <- aggregate(value ~ celltype, data = melt_mtx, FUN = sum)
    
    sorted_labels <- levels(srat$celltype.full)
    celltype_size$celltype <- factor(celltype_size$celltype,levels = sorted_labels)
    melt_mtx$celltype <- factor(melt_mtx$celltype,levels = sorted_labels)
    
    colnames(melt_mtx)[2] <- "dataset"
    
    
    p1 <- ggplot(celltype_size, aes(y= celltype,x = value)) + geom_bar(position="dodge", stat="identity",fill = "grey60") + 
      theme_bw() + scale_x_log10() + xlab("Cells per celltype, log10 scale") + ylab("")
    p2 <- ggplot(melt_mtx,aes(x=celltype,y=value,fill=dataset)) + 
      geom_bar(position="fill", stat="identity") + theme_bw() + coord_flip() +
      ylab("Fraction of cells in each dataset") + xlab("celltype number") + theme(legend.position="top")
    
    p2 + p1 + plot_layout(widths = c(3,1))
  }
  plot_integrated_celltypes_condition = function (srat) { 
    ## take an integrated Seurat object, plot distributions over orig.ident
    library(Seurat)
    library(patchwork)
    library(ggplot2)
    library(reshape2)
    library(RColorBrewer)
    
    
    count_table <- table(srat@meta.data$celltype.full, srat@meta.data$condition)
    count_mtx   <- as.data.frame.matrix(count_table)
    count_mtx$celltype <- rownames(count_mtx)
    melt_mtx    <- melt(count_mtx)
    melt_mtx$celltype <- as.factor(melt_mtx$celltype)
    
    celltype_size   <- aggregate(value ~ celltype, data = melt_mtx, FUN = sum)
    
    sorted_labels <- levels(srat$celltype.full)
    celltype_size$celltype <- factor(celltype_size$celltype,levels = sorted_labels)
    melt_mtx$celltype <- factor(melt_mtx$celltype,levels = sorted_labels)
    
    colnames(melt_mtx)[2] <- "dataset"
    
    
    p1 <- ggplot(celltype_size, aes(y= celltype,x = value)) + geom_bar(position="dodge", stat="identity",fill = "grey60") + 
      theme_bw() + scale_x_log10() + xlab("Cells per celltype, log10 scale") + ylab("")
    p2 <- ggplot(melt_mtx,aes(x=celltype,y=value,fill=dataset)) + 
      geom_bar(position="fill", stat="identity") + theme_bw() + coord_flip() + 
      scale_fill_manual(values= color_condition) +
      ylab("Fraction of cells in each dataset") + xlab("celltype number") + theme(legend.position="top")
    
    p2 + p1 + plot_layout(widths = c(3,1))
  }
  plot_integrated_celltypes_classification = function (srat) { 
    ## take an integrated Seurat object, plot distributions over orig.ident
    library(Seurat)
    library(patchwork)
    library(ggplot2)
    library(reshape2)
    library(RColorBrewer)
    
    
    count_table <- table(srat@meta.data$celltype.full, srat@meta.data$classification)
    count_mtx   <- as.data.frame.matrix(count_table)
    count_mtx$celltype <- rownames(count_mtx)
    melt_mtx    <- melt(count_mtx)
    melt_mtx$celltype <- as.factor(melt_mtx$celltype)
    
    celltype_size   <- aggregate(value ~ celltype, data = melt_mtx, FUN = sum)
    
    sorted_labels <- levels(srat$celltype.full)
    celltype_size$celltype <- factor(celltype_size$celltype,levels = sorted_labels)
    melt_mtx$celltype <- factor(melt_mtx$celltype,levels = sorted_labels)
    
    colnames(melt_mtx)[2] <- "dataset"
    
    
    p1 <- ggplot(celltype_size, aes(y= celltype,x = value)) + geom_bar(position="dodge", stat="identity",fill = "grey60") + 
      theme_bw() + scale_x_log10() + xlab("Cells per celltype, log10 scale") + ylab("")
    p2 <- ggplot(melt_mtx,aes(x=celltype,y=value,fill=dataset)) + 
      geom_bar(position="fill", stat="identity") + theme_bw() + coord_flip() + 
      scale_fill_manual(values= color_classification) +
      ylab("Fraction of cells in each dataset") + xlab("celltype number") + theme(legend.position="top")
    
    p2 + p1 + plot_layout(widths = c(3,1))
  }
  plot_integrated_celltypes_tissue = function (srat) { 
    ## take an integrated Seurat object, plot distributions over orig.ident
    library(Seurat)
    library(patchwork)
    library(ggplot2)
    library(reshape2)
    library(RColorBrewer)
    
    
    count_table <- table(srat@meta.data$celltype.full, srat@meta.data$tissue)
    count_mtx   <- as.data.frame.matrix(count_table)
    count_mtx$celltype <- rownames(count_mtx)
    melt_mtx    <- melt(count_mtx)
    melt_mtx$celltype <- as.factor(melt_mtx$celltype)
    
    celltype_size   <- aggregate(value ~ celltype, data = melt_mtx, FUN = sum)
    
    sorted_labels <- levels(srat$celltype.full)
    celltype_size$celltype <- factor(celltype_size$celltype,levels = sorted_labels)
    melt_mtx$celltype <- factor(melt_mtx$celltype,levels = sorted_labels)
    
    colnames(melt_mtx)[2] <- "dataset"
    
    
    p1 <- ggplot(celltype_size, aes(y= celltype,x = value)) + geom_bar(position="dodge", stat="identity",fill = "grey60") + 
      theme_bw() + scale_x_log10() + xlab("Cells per celltype, log10 scale") + ylab("")
    p2 <- ggplot(melt_mtx,aes(x=celltype,y=value,fill=dataset)) + 
      geom_bar(position="fill", stat="identity") + theme_bw() + coord_flip() + 
      scale_fill_manual(values= color_tissue) +
      ylab("Fraction of cells in each dataset") + xlab("celltype number") + theme(legend.position="top")
    
    p2 + p1 + plot_layout(widths = c(3,1))
  }
  plot_integrated_celltypes_sex = function (srat) { 
    ## take an integrated Seurat object, plot distributions over orig.ident
    library(Seurat)
    library(patchwork)
    library(ggplot2)
    library(reshape2)
    library(RColorBrewer)
    
    
    count_table <- table(srat@meta.data$celltype.full, srat@meta.data$sex)
    count_mtx   <- as.data.frame.matrix(count_table)
    count_mtx$celltype <- rownames(count_mtx)
    melt_mtx    <- melt(count_mtx)
    melt_mtx$celltype <- as.factor(melt_mtx$celltype)
    
    celltype_size   <- aggregate(value ~ celltype, data = melt_mtx, FUN = sum)
    
    sorted_labels <- levels(srat$celltype.full)
    celltype_size$celltype <- factor(celltype_size$celltype,levels = sorted_labels)
    melt_mtx$celltype <- factor(melt_mtx$celltype,levels = sorted_labels)
    
    colnames(melt_mtx)[2] <- "dataset"
    
    
    p1 <- ggplot(celltype_size, aes(y= celltype,x = value)) + geom_bar(position="dodge", stat="identity",fill = "grey60") + 
      theme_bw() + scale_x_log10() + xlab("Cells per celltype, log10 scale") + ylab("")
    p2 <- ggplot(melt_mtx,aes(x=celltype,y=value,fill=dataset)) + 
      geom_bar(position="fill", stat="identity") + theme_bw() + coord_flip() + 
      scale_fill_manual(values= color_sex) +
      ylab("Fraction of cells in each dataset") + xlab("celltype number") + theme(legend.position="top")
    
    p2 + p1 + plot_layout(widths = c(3,1))
  }
  plot_integrated_celltypes_age = function (srat) { 
    ## take an integrated Seurat object, plot distributions over orig.ident
    library(Seurat)
    library(patchwork)
    library(ggplot2)
    library(reshape2)
    library(RColorBrewer)
    
    
    count_table <- table(srat@meta.data$celltype.full, srat@meta.data$age)
    count_mtx   <- as.data.frame.matrix(count_table)
    count_mtx$celltype <- rownames(count_mtx)
    melt_mtx    <- melt(count_mtx)
    melt_mtx$celltype <- as.factor(melt_mtx$celltype)
    
    celltype_size   <- aggregate(value ~ celltype, data = melt_mtx, FUN = sum)
    
    sorted_labels <- levels(srat$celltype.full)
    celltype_size$celltype <- factor(celltype_size$celltype,levels = sorted_labels)
    melt_mtx$celltype <- factor(melt_mtx$celltype,levels = sorted_labels)
    
    colnames(melt_mtx)[2] <- "dataset"
    
    
    p1 <- ggplot(celltype_size, aes(y= celltype,x = value)) + geom_bar(position="dodge", stat="identity",fill = "grey60") + 
      theme_bw() + scale_x_log10() + xlab("Cells per celltype, log10 scale") + ylab("")
    p2 <- ggplot(melt_mtx,aes(x=celltype,y=value,fill=dataset)) + 
      geom_bar(position="fill", stat="identity") + theme_bw() + coord_flip() +
      ylab("Fraction of cells in each dataset") + xlab("celltype number") + theme(legend.position="top")
    
    p2 + p1 + plot_layout(widths = c(3,1))
  }
  plot_integrated_celltypes_development = function (srat) { 
    ## take an integrated Seurat object, plot distributions over orig.ident
    library(Seurat)
    library(patchwork)
    library(ggplot2)
    library(reshape2)
    library(RColorBrewer)
    
    
    count_table <- table(srat@meta.data$celltype.full, srat@meta.data$development)
    count_mtx   <- as.data.frame.matrix(count_table)
    count_mtx$celltype <- rownames(count_mtx)
    melt_mtx    <- melt(count_mtx)
    melt_mtx$celltype <- as.factor(melt_mtx$celltype)
    
    celltype_size   <- aggregate(value ~ celltype, data = melt_mtx, FUN = sum)
    
    sorted_labels <- levels(srat$celltype.full)
    celltype_size$celltype <- factor(celltype_size$celltype,levels = sorted_labels)
    melt_mtx$celltype <- factor(melt_mtx$celltype,levels = sorted_labels)
    
    colnames(melt_mtx)[2] <- "dataset"
    
    
    p1 <- ggplot(celltype_size, aes(y= celltype,x = value)) + geom_bar(position="dodge", stat="identity",fill = "grey60") + 
      theme_bw() + scale_x_log10() + xlab("Cells per celltype, log10 scale") + ylab("")
    p2 <- ggplot(melt_mtx,aes(x=celltype,y=value,fill=dataset)) + 
      geom_bar(position="fill", stat="identity") + theme_bw() + coord_flip() + 
      scale_fill_manual(values= color_development) +
      ylab("Fraction of cells in each dataset") + xlab("celltype number") + theme(legend.position="top")
    
    p2 + p1 + plot_layout(widths = c(3,1))
  }
  plot_integrated_celltypes_donor = function (srat) { 
    ## take an integrated Seurat object, plot distributions over orig.ident
    library(Seurat)
    library(patchwork)
    library(ggplot2)
    library(reshape2)
    library(RColorBrewer)
    
    
    count_table <- table(srat@meta.data$celltype.full, srat@meta.data$donor)
    count_mtx   <- as.data.frame.matrix(count_table)
    count_mtx$celltype <- rownames(count_mtx)
    melt_mtx    <- melt(count_mtx)
    melt_mtx$celltype <- as.factor(melt_mtx$celltype)
    
    celltype_size   <- aggregate(value ~ celltype, data = melt_mtx, FUN = sum)
    
    sorted_labels <- levels(srat$celltype.full)
    celltype_size$celltype <- factor(celltype_size$celltype,levels = sorted_labels)
    melt_mtx$celltype <- factor(melt_mtx$celltype,levels = sorted_labels)
    
    colnames(melt_mtx)[2] <- "dataset"
    
    
    p1 <- ggplot(celltype_size, aes(y= celltype,x = value)) + geom_bar(position="dodge", stat="identity",fill = "grey60") + 
      theme_bw() + scale_x_log10() + xlab("Cells per celltype, log10 scale") + ylab("")
    p2 <- ggplot(melt_mtx,aes(x=celltype,y=value,fill=dataset)) + 
      geom_bar(position="fill", stat="identity") + theme_bw() + coord_flip() +
      ylab("Fraction of cells in each dataset") + xlab("celltype number") + theme(legend.position="top")
    
    p2 + p1 + plot_layout(widths = c(3,1))
  }
  #
  
  #Main
  plot_main_integrated_group_clusters = function (srat) { 
    ## take an integrated Seurat object, plot distributions over orig.ident
    library(Seurat)
    library(patchwork)
    library(ggplot2)
    library(reshape2)
    library(RColorBrewer)
    
    
    count_table <- table(srat@meta.data$group, srat@meta.data$celltype.full)
    count_mtx   <- as.data.frame.matrix(count_table)
    count_mtx$celltype <- rownames(count_mtx)
    melt_mtx    <- melt(count_mtx)
    melt_mtx$celltype <- as.factor(melt_mtx$celltype)
    
    celltype_size   <- aggregate(value ~ celltype, data = melt_mtx, FUN = sum)
    
    sorted_labels <- levels(srat$group)
    celltype_size$celltype <- factor(celltype_size$celltype,levels = sorted_labels)
    melt_mtx$celltype <- factor(melt_mtx$celltype,levels = sorted_labels)
    
    colnames(melt_mtx)[2] <- "dataset"
    
    
    p1 <- ggplot(celltype_size, aes(y= celltype,x = value)) + geom_bar(position="dodge", stat="identity",fill = "grey60") + 
      theme_bw() + scale_x_log10() + xlab("Cells per celltype, log10 scale") + ylab("")
    p2 <- ggplot(melt_mtx,aes(x=celltype,y=value,fill=dataset)) + 
      geom_bar(position="fill", stat="identity") + theme_bw() + coord_flip() + 
      scale_fill_manual(values= color_level1) +
      ylab("Fraction of cells in each dataset") + xlab("celltype number") + theme(legend.position="top")
    
    p2 + p1 + plot_layout(widths = c(3,1))
  }
  plot_main_integrated_Donor_ID_clusters = function (srat) { 
    ## take an integrated Seurat object, plot distributions over orig.ident
    library(Seurat)
    library(patchwork)
    library(ggplot2)
    library(reshape2)
    library(RColorBrewer)
    
    
    count_table <- table(srat@meta.data$Donor_ID, srat@meta.data$celltype.full)
    count_mtx   <- as.data.frame.matrix(count_table)
    count_mtx$celltype <- rownames(count_mtx)
    melt_mtx    <- melt(count_mtx)
    melt_mtx$celltype <- as.factor(melt_mtx$celltype)
    
    celltype_size   <- aggregate(value ~ celltype, data = melt_mtx, FUN = sum)
    
    sorted_labels <- levels(srat$Donor_ID)
    celltype_size$celltype <- factor(celltype_size$celltype,levels = sorted_labels)
    melt_mtx$celltype <- factor(melt_mtx$celltype,levels = sorted_labels)
    
    colnames(melt_mtx)[2] <- "dataset"
    
    
    p1 <- ggplot(celltype_size, aes(y= celltype,x = value)) + geom_bar(position="dodge", stat="identity",fill = "grey60") + 
      theme_bw() + scale_x_log10() + xlab("Cells per celltype, log10 scale") + ylab("")
    p2 <- ggplot(melt_mtx,aes(x=celltype,y=value,fill=dataset)) + 
      geom_bar(position="fill", stat="identity") + theme_bw() + coord_flip() + 
      scale_fill_manual(values= color_level1)+
      ylab("Fraction of cells in each dataset") + xlab("celltype number") + theme(legend.position="top")
    
    p2 + p1 + plot_layout(widths = c(3,1))
  }
  plot_main_integrated_condition_celltypes = function (srat) { 
    ## take an integrated Seurat object, plot distributions over orig.ident
    library(Seurat)
    library(patchwork)
    library(ggplot2)
    library(reshape2)
    library(RColorBrewer)
    
    
    count_table <- table(srat@meta.data$condition, srat@meta.data$celltype.full)
    count_mtx   <- as.data.frame.matrix(count_table)
    count_mtx$celltype <- rownames(count_mtx)
    melt_mtx    <- melt(count_mtx)
    melt_mtx$celltype <- as.factor(melt_mtx$celltype)
    
    celltype_size   <- aggregate(value ~ celltype, data = melt_mtx, FUN = sum)
    
    sorted_labels <- levels(srat$condition)
    celltype_size$celltype <- factor(celltype_size$celltype,levels = sorted_labels)
    melt_mtx$celltype <- factor(melt_mtx$celltype,levels = sorted_labels)
    
    colnames(melt_mtx)[2] <- "dataset"
    
    
    p1 <- ggplot(celltype_size, aes(y= celltype,x = value)) + geom_bar(position="dodge", stat="identity",fill = "grey60") + 
      theme_bw() + scale_x_log10() + xlab("Cells per celltype, log10 scale") + ylab("")
    p2 <- ggplot(melt_mtx,aes(x=celltype,y=value,fill=dataset)) + 
      geom_bar(position="fill", stat="identity") + theme_bw() + coord_flip() + 
      scale_fill_manual(values= color_level1) +
      ylab("Fraction of cells in each dataset") + xlab("celltype number") + theme(legend.position="top")
    
    p2 + p1 + plot_layout(widths = c(3,1))
  }
  plot_main_integrated_classification_celltypes = function (srat) { 
    ## take an integrated Seurat object, plot distributions over orig.ident
    library(Seurat)
    library(patchwork)
    library(ggplot2)
    library(reshape2)
    library(RColorBrewer)
    
    
    count_table <- table(srat@meta.data$classification, srat@meta.data$celltype.full)
    count_mtx   <- as.data.frame.matrix(count_table)
    count_mtx$celltype <- rownames(count_mtx)
    melt_mtx    <- melt(count_mtx)
    melt_mtx$celltype <- as.factor(melt_mtx$celltype)
    
    celltype_size   <- aggregate(value ~ celltype, data = melt_mtx, FUN = sum)
    
    sorted_labels <- levels(srat$classification)
    celltype_size$celltype <- factor(celltype_size$celltype,levels = sorted_labels)
    melt_mtx$celltype <- factor(melt_mtx$celltype,levels = sorted_labels)
    
    colnames(melt_mtx)[2] <- "dataset"
    
    
    p1 <- ggplot(celltype_size, aes(y= celltype,x = value)) + geom_bar(position="dodge", stat="identity",fill = "grey60") + 
      theme_bw() + scale_x_log10() + xlab("Cells per celltype, log10 scale") + ylab("")
    p2 <- ggplot(melt_mtx,aes(x=celltype,y=value,fill=dataset)) + 
      geom_bar(position="fill", stat="identity") + theme_bw() + coord_flip() + 
      scale_fill_manual(values= color_level1) +
      ylab("Fraction of cells in each dataset") + xlab("celltype number") + theme(legend.position="top")
    
    p2 + p1 + plot_layout(widths = c(3,1))
  }
  plot_main_integrated_tissue_celltypes = function (srat) { 
    ## take an integrated Seurat object, plot distributions over orig.ident
    library(Seurat)
    library(patchwork)
    library(ggplot2)
    library(reshape2)
    library(RColorBrewer)
    
    
    count_table <- table(srat@meta.data$tissue, srat@meta.data$celltype.full)
    count_mtx   <- as.data.frame.matrix(count_table)
    count_mtx$celltype <- rownames(count_mtx)
    melt_mtx    <- melt(count_mtx)
    melt_mtx$celltype <- as.factor(melt_mtx$celltype)
    
    celltype_size   <- aggregate(value ~ celltype, data = melt_mtx, FUN = sum)
    
    sorted_labels <- levels(srat$tissue)
    celltype_size$celltype <- factor(celltype_size$celltype,levels = sorted_labels)
    melt_mtx$celltype <- factor(melt_mtx$celltype,levels = sorted_labels)
    
    colnames(melt_mtx)[2] <- "dataset"
    
    
    p1 <- ggplot(celltype_size, aes(y= celltype,x = value)) + geom_bar(position="dodge", stat="identity",fill = "grey60") + 
      theme_bw() + scale_x_log10() + xlab("Cells per celltype, log10 scale") + ylab("")
    p2 <- ggplot(melt_mtx,aes(x=celltype,y=value,fill=dataset)) + 
      geom_bar(position="fill", stat="identity") + theme_bw() + coord_flip() + 
      scale_fill_manual(values= color_level1) +
      ylab("Fraction of cells in each dataset") + xlab("celltype number") + theme(legend.position="top")
    
    p2 + p1 + plot_layout(widths = c(3,1))
  }
  plot_main_integrated_development_celltypes = function (srat) { 
    ## take an integrated Seurat object, plot distributions over orig.ident
    library(Seurat)
    library(patchwork)
    library(ggplot2)
    library(reshape2)
    library(RColorBrewer)
    
    
    count_table <- table(srat@meta.data$development, srat@meta.data$celltype.full)
    count_mtx   <- as.data.frame.matrix(count_table)
    count_mtx$celltype <- rownames(count_mtx)
    melt_mtx    <- melt(count_mtx)
    melt_mtx$celltype <- as.factor(melt_mtx$celltype)
    
    celltype_size   <- aggregate(value ~ celltype, data = melt_mtx, FUN = sum)
    
    sorted_labels <- levels(srat$development)
    celltype_size$celltype <- factor(celltype_size$celltype,levels = sorted_labels)
    melt_mtx$celltype <- factor(melt_mtx$celltype,levels = sorted_labels)
    
    colnames(melt_mtx)[2] <- "dataset"
    
    
    p1 <- ggplot(celltype_size, aes(y= celltype,x = value)) + geom_bar(position="dodge", stat="identity",fill = "grey60") + 
      theme_bw() + scale_x_log10() + xlab("Cells per celltype, log10 scale") + ylab("")
    p2 <- ggplot(melt_mtx,aes(x=celltype,y=value,fill=dataset)) + 
      geom_bar(position="fill", stat="identity") + theme_bw() + coord_flip() + 
      scale_fill_manual(values= color_level1) +
      ylab("Fraction of cells in each dataset") + xlab("celltype number") + theme(legend.position="top")
    
    p2 + p1 + plot_layout(widths = c(3,1))
  }
  plot_main_integrated_donor_celltypes = function (srat) { 
    ## take an integrated Seurat object, plot distributions over orig.ident
    library(Seurat)
    library(patchwork)
    library(ggplot2)
    library(reshape2)
    library(RColorBrewer)
    
    
    count_table <- table(srat@meta.data$donor, srat@meta.data$celltype.full)
    count_mtx   <- as.data.frame.matrix(count_table)
    count_mtx$celltype <- rownames(count_mtx)
    melt_mtx    <- melt(count_mtx)
    melt_mtx$celltype <- as.factor(melt_mtx$celltype)
    
    celltype_size   <- aggregate(value ~ celltype, data = melt_mtx, FUN = sum)
    
    sorted_labels <- levels(srat$donor)
    celltype_size$celltype <- factor(celltype_size$celltype,levels = sorted_labels)
    melt_mtx$celltype <- factor(melt_mtx$celltype,levels = sorted_labels)
    
    colnames(melt_mtx)[2] <- "dataset"
    
    
    p1 <- ggplot(celltype_size, aes(y= celltype,x = value)) + geom_bar(position="dodge", stat="identity",fill = "grey60") + 
      theme_bw() + scale_x_log10() + xlab("Cells per celltype, log10 scale") + ylab("")
    p2 <- ggplot(melt_mtx,aes(x=celltype,y=value,fill=dataset)) + 
      geom_bar(position="fill", stat="identity") + theme_bw() + coord_flip() + 
      scale_fill_manual(values= color_level1) +
      ylab("Fraction of cells in each dataset") + xlab("celltype number") + theme(legend.position="top")
    
    p2 + p1 + plot_layout(widths = c(3,1))
  }
  plot_main_integrated_group_tissue_clusters = function (srat) { 
    ## take an integrated Seurat object, plot distributions over orig.ident
    library(Seurat)
    library(patchwork)
    library(ggplot2)
    library(reshape2)
    library(RColorBrewer)
    
    
    count_table <- table(srat@meta.data$group.tissue, srat@meta.data$celltype.full)
    count_mtx   <- as.data.frame.matrix(count_table)
    count_mtx$celltype <- rownames(count_mtx)
    melt_mtx    <- melt(count_mtx)
    melt_mtx$celltype <- as.factor(melt_mtx$celltype)
    
    celltype_size   <- aggregate(value ~ celltype, data = melt_mtx, FUN = sum)
    
    sorted_labels <- levels(srat$group.tissue)
    celltype_size$celltype <- factor(celltype_size$celltype,levels = sorted_labels)
    melt_mtx$celltype <- factor(melt_mtx$celltype,levels = sorted_labels)
    
    colnames(melt_mtx)[2] <- "dataset"
    
    
    p1 <- ggplot(celltype_size, aes(y= celltype,x = value)) + geom_bar(position="dodge", stat="identity",fill = "grey60") + 
      theme_bw() + scale_x_log10() + xlab("Cells per celltype, log10 scale") + ylab("")
    p2 <- ggplot(melt_mtx,aes(x=celltype,y=value,fill=dataset)) + 
      geom_bar(position="fill", stat="identity") + theme_bw() + coord_flip() + 
      scale_fill_manual(values= color_level1) +
      ylab("Fraction of datasets in each celltype") + xlab("celltype number") + theme(legend.position="top")
    
    p2 + p1 + plot_layout(widths = c(3,1))
  }
  #
  #TC
  plot_TC_integrated_group_clusters = function (srat) { 
    ## take an integrated Seurat object, plot distributions over orig.ident
    library(Seurat)
    library(patchwork)
    library(ggplot2)
    library(reshape2)
    library(RColorBrewer)
    
    
    count_table <- table(srat@meta.data$group, srat@meta.data$celltype.full)
    count_mtx   <- as.data.frame.matrix(count_table)
    count_mtx$celltype <- rownames(count_mtx)
    melt_mtx    <- melt(count_mtx)
    melt_mtx$celltype <- as.factor(melt_mtx$celltype)
    
    celltype_size   <- aggregate(value ~ celltype, data = melt_mtx, FUN = sum)
    
    sorted_labels <- levels(srat$group)
    celltype_size$celltype <- factor(celltype_size$celltype,levels = sorted_labels)
    melt_mtx$celltype <- factor(melt_mtx$celltype,levels = sorted_labels)
    
    colnames(melt_mtx)[2] <- "dataset"
    
    
    p1 <- ggplot(celltype_size, aes(y= celltype,x = value)) + geom_bar(position="dodge", stat="identity",fill = "grey60") + 
      theme_bw() + scale_x_log10() + xlab("Cells per celltype, log10 scale") + ylab("")
    p2 <- ggplot(melt_mtx,aes(x=celltype,y=value,fill=dataset)) + 
      geom_bar(position="fill", stat="identity") + theme_bw() + coord_flip() + 
      scale_fill_manual(values= color_levelTC) +
      ylab("Fraction of cells in each dataset") + xlab("celltype number") + theme(legend.position="top")
    
    p2 + p1 + plot_layout(widths = c(3,1))
  }
  plot_TC_integrated_Donor_ID_clusters = function (srat) { 
    ## take an integrated Seurat object, plot distributions over orig.ident
    library(Seurat)
    library(patchwork)
    library(ggplot2)
    library(reshape2)
    library(RColorBrewer)
    
    
    count_table <- table(srat@meta.data$Donor_ID, srat@meta.data$celltype.full)
    count_mtx   <- as.data.frame.matrix(count_table)
    count_mtx$celltype <- rownames(count_mtx)
    melt_mtx    <- melt(count_mtx)
    melt_mtx$celltype <- as.factor(melt_mtx$celltype)
    
    celltype_size   <- aggregate(value ~ celltype, data = melt_mtx, FUN = sum)
    
    sorted_labels <- levels(srat$Donor_ID)
    celltype_size$celltype <- factor(celltype_size$celltype,levels = sorted_labels)
    melt_mtx$celltype <- factor(melt_mtx$celltype,levels = sorted_labels)
    
    colnames(melt_mtx)[2] <- "dataset"
    
    
    p1 <- ggplot(celltype_size, aes(y= celltype,x = value)) + geom_bar(position="dodge", stat="identity",fill = "grey60") + 
      theme_bw() + scale_x_log10() + xlab("Cells per celltype, log10 scale") + ylab("")
    p2 <- ggplot(melt_mtx,aes(x=celltype,y=value,fill=dataset)) + 
      geom_bar(position="fill", stat="identity") + theme_bw() + coord_flip() + 
      scale_fill_manual(values= color_levelTC)+
      ylab("Fraction of cells in each dataset") + xlab("celltype number") + theme(legend.position="top")
    
    p2 + p1 + plot_layout(widths = c(3,1))
  }
  plot_TC_integrated_condition_celltypes = function (srat) { 
    ## take an integrated Seurat object, plot distributions over orig.ident
    library(Seurat)
    library(patchwork)
    library(ggplot2)
    library(reshape2)
    library(RColorBrewer)
    
    
    count_table <- table(srat@meta.data$condition, srat@meta.data$celltype.full)
    count_mtx   <- as.data.frame.matrix(count_table)
    count_mtx$celltype <- rownames(count_mtx)
    melt_mtx    <- melt(count_mtx)
    melt_mtx$celltype <- as.factor(melt_mtx$celltype)
    
    celltype_size   <- aggregate(value ~ celltype, data = melt_mtx, FUN = sum)
    
    sorted_labels <- levels(srat$condition)
    celltype_size$celltype <- factor(celltype_size$celltype,levels = sorted_labels)
    melt_mtx$celltype <- factor(melt_mtx$celltype,levels = sorted_labels)
    
    colnames(melt_mtx)[2] <- "dataset"
    
    
    p1 <- ggplot(celltype_size, aes(y= celltype,x = value)) + geom_bar(position="dodge", stat="identity",fill = "grey60") + 
      theme_bw() + scale_x_log10() + xlab("Cells per celltype, log10 scale") + ylab("")
    p2 <- ggplot(melt_mtx,aes(x=celltype,y=value,fill=dataset)) + 
      geom_bar(position="fill", stat="identity") + theme_bw() + coord_flip() + 
      scale_fill_manual(values= color_levelTC) +
      ylab("Fraction of cells in each dataset") + xlab("celltype number") + theme(legend.position="top")
    
    p2 + p1 + plot_layout(widths = c(3,1))
  }
  plot_TC_integrated_classification_celltypes = function (srat) { 
    ## take an integrated Seurat object, plot distributions over orig.ident
    library(Seurat)
    library(patchwork)
    library(ggplot2)
    library(reshape2)
    library(RColorBrewer)
    
    
    count_table <- table(srat@meta.data$classification, srat@meta.data$celltype.full)
    count_mtx   <- as.data.frame.matrix(count_table)
    count_mtx$celltype <- rownames(count_mtx)
    melt_mtx    <- melt(count_mtx)
    melt_mtx$celltype <- as.factor(melt_mtx$celltype)
    
    celltype_size   <- aggregate(value ~ celltype, data = melt_mtx, FUN = sum)
    
    sorted_labels <- levels(srat$classification)
    celltype_size$celltype <- factor(celltype_size$celltype,levels = sorted_labels)
    melt_mtx$celltype <- factor(melt_mtx$celltype,levels = sorted_labels)
    
    colnames(melt_mtx)[2] <- "dataset"
    
    
    p1 <- ggplot(celltype_size, aes(y= celltype,x = value)) + geom_bar(position="dodge", stat="identity",fill = "grey60") + 
      theme_bw() + scale_x_log10() + xlab("Cells per celltype, log10 scale") + ylab("")
    p2 <- ggplot(melt_mtx,aes(x=celltype,y=value,fill=dataset)) + 
      geom_bar(position="fill", stat="identity") + theme_bw() + coord_flip() + 
      scale_fill_manual(values= color_levelTC) +
      ylab("Fraction of cells in each dataset") + xlab("celltype number") + theme(legend.position="top")
    
    p2 + p1 + plot_layout(widths = c(3,1))
  }
  plot_TC_integrated_tissue_celltypes = function (srat) { 
    ## take an integrated Seurat object, plot distributions over orig.ident
    library(Seurat)
    library(patchwork)
    library(ggplot2)
    library(reshape2)
    library(RColorBrewer)
    
    
    count_table <- table(srat@meta.data$tissue, srat@meta.data$celltype.full)
    count_mtx   <- as.data.frame.matrix(count_table)
    count_mtx$celltype <- rownames(count_mtx)
    melt_mtx    <- melt(count_mtx)
    melt_mtx$celltype <- as.factor(melt_mtx$celltype)
    
    celltype_size   <- aggregate(value ~ celltype, data = melt_mtx, FUN = sum)
    
    sorted_labels <- levels(srat$tissue)
    celltype_size$celltype <- factor(celltype_size$celltype,levels = sorted_labels)
    melt_mtx$celltype <- factor(melt_mtx$celltype,levels = sorted_labels)
    
    colnames(melt_mtx)[2] <- "dataset"
    
    
    p1 <- ggplot(celltype_size, aes(y= celltype,x = value)) + geom_bar(position="dodge", stat="identity",fill = "grey60") + 
      theme_bw() + scale_x_log10() + xlab("Cells per celltype, log10 scale") + ylab("")
    p2 <- ggplot(melt_mtx,aes(x=celltype,y=value,fill=dataset)) + 
      geom_bar(position="fill", stat="identity") + theme_bw() + coord_flip() + 
      scale_fill_manual(values= color_levelTC) +
      ylab("Fraction of cells in each dataset") + xlab("celltype number") + theme(legend.position="top")
    
    p2 + p1 + plot_layout(widths = c(3,1))
  }
  plot_TC_integrated_development_celltypes = function (srat) { 
    ## take an integrated Seurat object, plot distributions over orig.ident
    library(Seurat)
    library(patchwork)
    library(ggplot2)
    library(reshape2)
    library(RColorBrewer)
    
    
    count_table <- table(srat@meta.data$development, srat@meta.data$celltype.full)
    count_mtx   <- as.data.frame.matrix(count_table)
    count_mtx$celltype <- rownames(count_mtx)
    melt_mtx    <- melt(count_mtx)
    melt_mtx$celltype <- as.factor(melt_mtx$celltype)
    
    celltype_size   <- aggregate(value ~ celltype, data = melt_mtx, FUN = sum)
    
    sorted_labels <- levels(srat$development)
    celltype_size$celltype <- factor(celltype_size$celltype,levels = sorted_labels)
    melt_mtx$celltype <- factor(melt_mtx$celltype,levels = sorted_labels)
    
    colnames(melt_mtx)[2] <- "dataset"
    
    
    p1 <- ggplot(celltype_size, aes(y= celltype,x = value)) + geom_bar(position="dodge", stat="identity",fill = "grey60") + 
      theme_bw() + scale_x_log10() + xlab("Cells per celltype, log10 scale") + ylab("")
    p2 <- ggplot(melt_mtx,aes(x=celltype,y=value,fill=dataset)) + 
      geom_bar(position="fill", stat="identity") + theme_bw() + coord_flip() + 
      scale_fill_manual(values= color_levelTC) +
      ylab("Fraction of cells in each dataset") + xlab("celltype number") + theme(legend.position="top")
    
    p2 + p1 + plot_layout(widths = c(3,1))
  }
  plot_TC_integrated_donor_celltypes = function (srat) { 
    ## take an integrated Seurat object, plot distributions over orig.ident
    library(Seurat)
    library(patchwork)
    library(ggplot2)
    library(reshape2)
    library(RColorBrewer)
    
    
    count_table <- table(srat@meta.data$donor, srat@meta.data$celltype.full)
    count_mtx   <- as.data.frame.matrix(count_table)
    count_mtx$celltype <- rownames(count_mtx)
    melt_mtx    <- melt(count_mtx)
    melt_mtx$celltype <- as.factor(melt_mtx$celltype)
    
    celltype_size   <- aggregate(value ~ celltype, data = melt_mtx, FUN = sum)
    
    sorted_labels <- levels(srat$donor)
    celltype_size$celltype <- factor(celltype_size$celltype,levels = sorted_labels)
    melt_mtx$celltype <- factor(melt_mtx$celltype,levels = sorted_labels)
    
    colnames(melt_mtx)[2] <- "dataset"
    
    
    p1 <- ggplot(celltype_size, aes(y= celltype,x = value)) + geom_bar(position="dodge", stat="identity",fill = "grey60") + 
      theme_bw() + scale_x_log10() + xlab("Cells per celltype, log10 scale") + ylab("")
    p2 <- ggplot(melt_mtx,aes(x=celltype,y=value,fill=dataset)) + 
      geom_bar(position="fill", stat="identity") + theme_bw() + coord_flip() + 
      scale_fill_manual(values= color_levelTC) +
      ylab("Fraction of cells in each dataset") + xlab("celltype number") + theme(legend.position="top")
    
    p2 + p1 + plot_layout(widths = c(3,1))
  }
  #
  #BC
  plot_BC_integrated_group_clusters = function (srat) { 
    ## take an integrated Seurat object, plot distributions over orig.ident
    library(Seurat)
    library(patchwork)
    library(ggplot2)
    library(reshape2)
    library(RColorBrewer)
    
    
    count_table <- table(srat@meta.data$group, srat@meta.data$celltype.full)
    count_mtx   <- as.data.frame.matrix(count_table)
    count_mtx$celltype <- rownames(count_mtx)
    melt_mtx    <- melt(count_mtx)
    melt_mtx$celltype <- as.factor(melt_mtx$celltype)
    
    celltype_size   <- aggregate(value ~ celltype, data = melt_mtx, FUN = sum)
    
    sorted_labels <- levels(srat$group)
    celltype_size$celltype <- factor(celltype_size$celltype,levels = sorted_labels)
    melt_mtx$celltype <- factor(melt_mtx$celltype,levels = sorted_labels)
    
    colnames(melt_mtx)[2] <- "dataset"
    
    
    p1 <- ggplot(celltype_size, aes(y= celltype,x = value)) + geom_bar(position="dodge", stat="identity",fill = "grey60") + 
      theme_bw() + scale_x_log10() + xlab("Cells per celltype, log10 scale") + ylab("")
    p2 <- ggplot(melt_mtx,aes(x=celltype,y=value,fill=dataset)) + 
      geom_bar(position="fill", stat="identity") + theme_bw() + coord_flip() + 
      scale_fill_manual(values= color_levelBC) +
      ylab("Fraction of cells in each dataset") + xlab("celltype number") + theme(legend.position="top")
    
    p2 + p1 + plot_layout(widths = c(3,1))
  }
  plot_BC_integrated_Donor_ID_clusters = function (srat) { 
    ## take an integrated Seurat object, plot distributions over orig.ident
    library(Seurat)
    library(patchwork)
    library(ggplot2)
    library(reshape2)
    library(RColorBrewer)
    
    
    count_table <- table(srat@meta.data$Donor_ID, srat@meta.data$celltype.full)
    count_mtx   <- as.data.frame.matrix(count_table)
    count_mtx$celltype <- rownames(count_mtx)
    melt_mtx    <- melt(count_mtx)
    melt_mtx$celltype <- as.factor(melt_mtx$celltype)
    
    celltype_size   <- aggregate(value ~ celltype, data = melt_mtx, FUN = sum)
    
    sorted_labels <- levels(srat$Donor_ID)
    celltype_size$celltype <- factor(celltype_size$celltype,levels = sorted_labels)
    melt_mtx$celltype <- factor(melt_mtx$celltype,levels = sorted_labels)
    
    colnames(melt_mtx)[2] <- "dataset"
    
    
    p1 <- ggplot(celltype_size, aes(y= celltype,x = value)) + geom_bar(position="dodge", stat="identity",fill = "grey60") + 
      theme_bw() + scale_x_log10() + xlab("Cells per celltype, log10 scale") + ylab("")
    p2 <- ggplot(melt_mtx,aes(x=celltype,y=value,fill=dataset)) + 
      geom_bar(position="fill", stat="identity") + theme_bw() + coord_flip() + 
      scale_fill_manual(values= color_levelBC)+
      ylab("Fraction of cells in each dataset") + xlab("celltype number") + theme(legend.position="top")
    
    p2 + p1 + plot_layout(widths = c(3,1))
  }
  plot_BC_integrated_condition_celltypes = function (srat) { 
    ## take an integrated Seurat object, plot distributions over orig.ident
    library(Seurat)
    library(patchwork)
    library(ggplot2)
    library(reshape2)
    library(RColorBrewer)
    
    
    count_table <- table(srat@meta.data$condition, srat@meta.data$celltype.full)
    count_mtx   <- as.data.frame.matrix(count_table)
    count_mtx$celltype <- rownames(count_mtx)
    melt_mtx    <- melt(count_mtx)
    melt_mtx$celltype <- as.factor(melt_mtx$celltype)
    
    celltype_size   <- aggregate(value ~ celltype, data = melt_mtx, FUN = sum)
    
    sorted_labels <- levels(srat$condition)
    celltype_size$celltype <- factor(celltype_size$celltype,levels = sorted_labels)
    melt_mtx$celltype <- factor(melt_mtx$celltype,levels = sorted_labels)
    
    colnames(melt_mtx)[2] <- "dataset"
    
    
    p1 <- ggplot(celltype_size, aes(y= celltype,x = value)) + geom_bar(position="dodge", stat="identity",fill = "grey60") + 
      theme_bw() + scale_x_log10() + xlab("Cells per celltype, log10 scale") + ylab("")
    p2 <- ggplot(melt_mtx,aes(x=celltype,y=value,fill=dataset)) + 
      geom_bar(position="fill", stat="identity") + theme_bw() + coord_flip() + 
      scale_fill_manual(values= color_levelBC) +
      ylab("Fraction of cells in each dataset") + xlab("celltype number") + theme(legend.position="top")
    
    p2 + p1 + plot_layout(widths = c(3,1))
  }
  plot_BC_integrated_classification_celltypes = function (srat) { 
    ## take an integrated Seurat object, plot distributions over orig.ident
    library(Seurat)
    library(patchwork)
    library(ggplot2)
    library(reshape2)
    library(RColorBrewer)
    
    
    count_table <- table(srat@meta.data$classification, srat@meta.data$celltype.full)
    count_mtx   <- as.data.frame.matrix(count_table)
    count_mtx$celltype <- rownames(count_mtx)
    melt_mtx    <- melt(count_mtx)
    melt_mtx$celltype <- as.factor(melt_mtx$celltype)
    
    celltype_size   <- aggregate(value ~ celltype, data = melt_mtx, FUN = sum)
    
    sorted_labels <- levels(srat$classification)
    celltype_size$celltype <- factor(celltype_size$celltype,levels = sorted_labels)
    melt_mtx$celltype <- factor(melt_mtx$celltype,levels = sorted_labels)
    
    colnames(melt_mtx)[2] <- "dataset"
    
    
    p1 <- ggplot(celltype_size, aes(y= celltype,x = value)) + geom_bar(position="dodge", stat="identity",fill = "grey60") + 
      theme_bw() + scale_x_log10() + xlab("Cells per celltype, log10 scale") + ylab("")
    p2 <- ggplot(melt_mtx,aes(x=celltype,y=value,fill=dataset)) + 
      geom_bar(position="fill", stat="identity") + theme_bw() + coord_flip() + 
      scale_fill_manual(values= color_levelBC) +
      ylab("Fraction of cells in each dataset") + xlab("celltype number") + theme(legend.position="top")
    
    p2 + p1 + plot_layout(widths = c(3,1))
  }
  plot_BC_integrated_tissue_celltypes = function (srat) { 
    ## take an integrated Seurat object, plot distributions over orig.ident
    library(Seurat)
    library(patchwork)
    library(ggplot2)
    library(reshape2)
    library(RColorBrewer)
    
    
    count_table <- table(srat@meta.data$tissue, srat@meta.data$celltype.full)
    count_mtx   <- as.data.frame.matrix(count_table)
    count_mtx$celltype <- rownames(count_mtx)
    melt_mtx    <- melt(count_mtx)
    melt_mtx$celltype <- as.factor(melt_mtx$celltype)
    
    celltype_size   <- aggregate(value ~ celltype, data = melt_mtx, FUN = sum)
    
    sorted_labels <- levels(srat$tissue)
    celltype_size$celltype <- factor(celltype_size$celltype,levels = sorted_labels)
    melt_mtx$celltype <- factor(melt_mtx$celltype,levels = sorted_labels)
    
    colnames(melt_mtx)[2] <- "dataset"
    
    
    p1 <- ggplot(celltype_size, aes(y= celltype,x = value)) + geom_bar(position="dodge", stat="identity",fill = "grey60") + 
      theme_bw() + scale_x_log10() + xlab("Cells per celltype, log10 scale") + ylab("")
    p2 <- ggplot(melt_mtx,aes(x=celltype,y=value,fill=dataset)) + 
      geom_bar(position="fill", stat="identity") + theme_bw() + coord_flip() + 
      scale_fill_manual(values= color_levelBC) +
      ylab("Fraction of cells in each dataset") + xlab("celltype number") + theme(legend.position="top")
    
    p2 + p1 + plot_layout(widths = c(3,1))
  }
  plot_BC_integrated_development_celltypes = function (srat) { 
    ## take an integrated Seurat object, plot distributions over orig.ident
    library(Seurat)
    library(patchwork)
    library(ggplot2)
    library(reshape2)
    library(RColorBrewer)
    
    
    count_table <- table(srat@meta.data$development, srat@meta.data$celltype.full)
    count_mtx   <- as.data.frame.matrix(count_table)
    count_mtx$celltype <- rownames(count_mtx)
    melt_mtx    <- melt(count_mtx)
    melt_mtx$celltype <- as.factor(melt_mtx$celltype)
    
    celltype_size   <- aggregate(value ~ celltype, data = melt_mtx, FUN = sum)
    
    sorted_labels <- levels(srat$development)
    celltype_size$celltype <- factor(celltype_size$celltype,levels = sorted_labels)
    melt_mtx$celltype <- factor(melt_mtx$celltype,levels = sorted_labels)
    
    colnames(melt_mtx)[2] <- "dataset"
    
    
    p1 <- ggplot(celltype_size, aes(y= celltype,x = value)) + geom_bar(position="dodge", stat="identity",fill = "grey60") + 
      theme_bw() + scale_x_log10() + xlab("Cells per celltype, log10 scale") + ylab("")
    p2 <- ggplot(melt_mtx,aes(x=celltype,y=value,fill=dataset)) + 
      geom_bar(position="fill", stat="identity") + theme_bw() + coord_flip() + 
      scale_fill_manual(values= color_levelBC) +
      ylab("Fraction of cells in each dataset") + xlab("celltype number") + theme(legend.position="top")
    
    p2 + p1 + plot_layout(widths = c(3,1))
  }
  plot_BC_integrated_donor_celltypes = function (srat) { 
    ## take an integrated Seurat object, plot distributions over orig.ident
    library(Seurat)
    library(patchwork)
    library(ggplot2)
    library(reshape2)
    library(RColorBrewer)
    
    
    count_table <- table(srat@meta.data$donor, srat@meta.data$celltype.full)
    count_mtx   <- as.data.frame.matrix(count_table)
    count_mtx$celltype <- rownames(count_mtx)
    melt_mtx    <- melt(count_mtx)
    melt_mtx$celltype <- as.factor(melt_mtx$celltype)
    
    celltype_size   <- aggregate(value ~ celltype, data = melt_mtx, FUN = sum)
    
    sorted_labels <- levels(srat$donor)
    celltype_size$celltype <- factor(celltype_size$celltype,levels = sorted_labels)
    melt_mtx$celltype <- factor(melt_mtx$celltype,levels = sorted_labels)
    
    colnames(melt_mtx)[2] <- "dataset"
    
    
    p1 <- ggplot(celltype_size, aes(y= celltype,x = value)) + geom_bar(position="dodge", stat="identity",fill = "grey60") + 
      theme_bw() + scale_x_log10() + xlab("Cells per celltype, log10 scale") + ylab("")
    p2 <- ggplot(melt_mtx,aes(x=celltype,y=value,fill=dataset)) + 
      geom_bar(position="fill", stat="identity") + theme_bw() + coord_flip() + 
      scale_fill_manual(values= color_levelBC) +
      ylab("Fraction of cells in each dataset") + xlab("celltype number") + theme(legend.position="top")
    
    p2 + p1 + plot_layout(widths = c(3,1))
  }
  ##
  #PC
  plot_PC_integrated_group_clusters = function (srat) { 
    ## take an integrated Seurat object, plot distributions over orig.ident
    library(Seurat)
    library(patchwork)
    library(ggplot2)
    library(reshape2)
    library(RColorBrewer)
    
    
    count_table <- table(srat@meta.data$group, srat@meta.data$celltype.full)
    count_mtx   <- as.data.frame.matrix(count_table)
    count_mtx$celltype <- rownames(count_mtx)
    melt_mtx    <- melt(count_mtx)
    melt_mtx$celltype <- as.factor(melt_mtx$celltype)
    
    celltype_size   <- aggregate(value ~ celltype, data = melt_mtx, FUN = sum)
    
    sorted_labels <- levels(srat$group)
    celltype_size$celltype <- factor(celltype_size$celltype,levels = sorted_labels)
    melt_mtx$celltype <- factor(melt_mtx$celltype,levels = sorted_labels)
    
    colnames(melt_mtx)[2] <- "dataset"
    
    
    p1 <- ggplot(celltype_size, aes(y= celltype,x = value)) + geom_bar(position="dodge", stat="identity",fill = "grey60") + 
      theme_bw() + scale_x_log10() + xlab("Cells per celltype, log10 scale") + ylab("")
    p2 <- ggplot(melt_mtx,aes(x=celltype,y=value,fill=dataset)) + 
      geom_bar(position="fill", stat="identity") + theme_bw() + coord_flip() + 
      scale_fill_manual(values= color_levelPC) +
      ylab("Fraction of cells in each dataset") + xlab("celltype number") + theme(legend.position="top")
    
    p2 + p1 + plot_layout(widths = c(3,1))
  }
  plot_PC_integrated_Donor_ID_clusters = function (srat) { 
    ## take an integrated Seurat object, plot distributions over orig.ident
    library(Seurat)
    library(patchwork)
    library(ggplot2)
    library(reshape2)
    library(RColorBrewer)
    
    
    count_table <- table(srat@meta.data$Donor_ID, srat@meta.data$celltype.full)
    count_mtx   <- as.data.frame.matrix(count_table)
    count_mtx$celltype <- rownames(count_mtx)
    melt_mtx    <- melt(count_mtx)
    melt_mtx$celltype <- as.factor(melt_mtx$celltype)
    
    celltype_size   <- aggregate(value ~ celltype, data = melt_mtx, FUN = sum)
    
    sorted_labels <- levels(srat$Donor_ID)
    celltype_size$celltype <- factor(celltype_size$celltype,levels = sorted_labels)
    melt_mtx$celltype <- factor(melt_mtx$celltype,levels = sorted_labels)
    
    colnames(melt_mtx)[2] <- "dataset"
    
    
    p1 <- ggplot(celltype_size, aes(y= celltype,x = value)) + geom_bar(position="dodge", stat="identity",fill = "grey60") + 
      theme_bw() + scale_x_log10() + xlab("Cells per celltype, log10 scale") + ylab("")
    p2 <- ggplot(melt_mtx,aes(x=celltype,y=value,fill=dataset)) + 
      geom_bar(position="fill", stat="identity") + theme_bw() + coord_flip() + 
      scale_fill_manual(values= color_levelPC)+
      ylab("Fraction of cells in each dataset") + xlab("celltype number") + theme(legend.position="top")
    
    p2 + p1 + plot_layout(widths = c(3,1))
  }
  plot_PC_integrated_condition_celltypes = function (srat) { 
    ## take an integrated Seurat object, plot distributions over orig.ident
    library(Seurat)
    library(patchwork)
    library(ggplot2)
    library(reshape2)
    library(RColorBrewer)
    
    
    count_table <- table(srat@meta.data$condition, srat@meta.data$celltype.full)
    count_mtx   <- as.data.frame.matrix(count_table)
    count_mtx$celltype <- rownames(count_mtx)
    melt_mtx    <- melt(count_mtx)
    melt_mtx$celltype <- as.factor(melt_mtx$celltype)
    
    celltype_size   <- aggregate(value ~ celltype, data = melt_mtx, FUN = sum)
    
    sorted_labels <- levels(srat$condition)
    celltype_size$celltype <- factor(celltype_size$celltype,levels = sorted_labels)
    melt_mtx$celltype <- factor(melt_mtx$celltype,levels = sorted_labels)
    
    colnames(melt_mtx)[2] <- "dataset"
    
    
    p1 <- ggplot(celltype_size, aes(y= celltype,x = value)) + geom_bar(position="dodge", stat="identity",fill = "grey60") + 
      theme_bw() + scale_x_log10() + xlab("Cells per celltype, log10 scale") + ylab("")
    p2 <- ggplot(melt_mtx,aes(x=celltype,y=value,fill=dataset)) + 
      geom_bar(position="fill", stat="identity") + theme_bw() + coord_flip() + 
      scale_fill_manual(values= color_levelPC) +
      ylab("Fraction of cells in each dataset") + xlab("celltype number") + theme(legend.position="top")
    
    p2 + p1 + plot_layout(widths = c(3,1))
  }
  plot_PC_integrated_classification_celltypes = function (srat) { 
    ## take an integrated Seurat object, plot distributions over orig.ident
    library(Seurat)
    library(patchwork)
    library(ggplot2)
    library(reshape2)
    library(RColorBrewer)
    
    
    count_table <- table(srat@meta.data$classification, srat@meta.data$celltype.full)
    count_mtx   <- as.data.frame.matrix(count_table)
    count_mtx$celltype <- rownames(count_mtx)
    melt_mtx    <- melt(count_mtx)
    melt_mtx$celltype <- as.factor(melt_mtx$celltype)
    
    celltype_size   <- aggregate(value ~ celltype, data = melt_mtx, FUN = sum)
    
    sorted_labels <- levels(srat$classification)
    celltype_size$celltype <- factor(celltype_size$celltype,levels = sorted_labels)
    melt_mtx$celltype <- factor(melt_mtx$celltype,levels = sorted_labels)
    
    colnames(melt_mtx)[2] <- "dataset"
    
    
    p1 <- ggplot(celltype_size, aes(y= celltype,x = value)) + geom_bar(position="dodge", stat="identity",fill = "grey60") + 
      theme_bw() + scale_x_log10() + xlab("Cells per celltype, log10 scale") + ylab("")
    p2 <- ggplot(melt_mtx,aes(x=celltype,y=value,fill=dataset)) + 
      geom_bar(position="fill", stat="identity") + theme_bw() + coord_flip() + 
      scale_fill_manual(values= color_levelPC) +
      ylab("Fraction of cells in each dataset") + xlab("celltype number") + theme(legend.position="top")
    
    p2 + p1 + plot_layout(widths = c(3,1))
  }
  plot_PC_integrated_tissue_celltypes = function (srat) { 
    ## take an integrated Seurat object, plot distributions over orig.ident
    library(Seurat)
    library(patchwork)
    library(ggplot2)
    library(reshape2)
    library(RColorBrewer)
    
    
    count_table <- table(srat@meta.data$tissue, srat@meta.data$celltype.full)
    count_mtx   <- as.data.frame.matrix(count_table)
    count_mtx$celltype <- rownames(count_mtx)
    melt_mtx    <- melt(count_mtx)
    melt_mtx$celltype <- as.factor(melt_mtx$celltype)
    
    celltype_size   <- aggregate(value ~ celltype, data = melt_mtx, FUN = sum)
    
    sorted_labels <- levels(srat$tissue)
    celltype_size$celltype <- factor(celltype_size$celltype,levels = sorted_labels)
    melt_mtx$celltype <- factor(melt_mtx$celltype,levels = sorted_labels)
    
    colnames(melt_mtx)[2] <- "dataset"
    
    
    p1 <- ggplot(celltype_size, aes(y= celltype,x = value)) + geom_bar(position="dodge", stat="identity",fill = "grey60") + 
      theme_bw() + scale_x_log10() + xlab("Cells per celltype, log10 scale") + ylab("")
    p2 <- ggplot(melt_mtx,aes(x=celltype,y=value,fill=dataset)) + 
      geom_bar(position="fill", stat="identity") + theme_bw() + coord_flip() + 
      scale_fill_manual(values= color_levelPC) +
      ylab("Fraction of cells in each dataset") + xlab("celltype number") + theme(legend.position="top")
    
    p2 + p1 + plot_layout(widths = c(3,1))
  }
  plot_PC_integrated_development_celltypes = function (srat) { 
    ## take an integrated Seurat object, plot distributions over orig.ident
    library(Seurat)
    library(patchwork)
    library(ggplot2)
    library(reshape2)
    library(RColorBrewer)
    
    
    count_table <- table(srat@meta.data$development, srat@meta.data$celltype.full)
    count_mtx   <- as.data.frame.matrix(count_table)
    count_mtx$celltype <- rownames(count_mtx)
    melt_mtx    <- melt(count_mtx)
    melt_mtx$celltype <- as.factor(melt_mtx$celltype)
    
    celltype_size   <- aggregate(value ~ celltype, data = melt_mtx, FUN = sum)
    
    sorted_labels <- levels(srat$development)
    celltype_size$celltype <- factor(celltype_size$celltype,levels = sorted_labels)
    melt_mtx$celltype <- factor(melt_mtx$celltype,levels = sorted_labels)
    
    colnames(melt_mtx)[2] <- "dataset"
    
    
    p1 <- ggplot(celltype_size, aes(y= celltype,x = value)) + geom_bar(position="dodge", stat="identity",fill = "grey60") + 
      theme_bw() + scale_x_log10() + xlab("Cells per celltype, log10 scale") + ylab("")
    p2 <- ggplot(melt_mtx,aes(x=celltype,y=value,fill=dataset)) + 
      geom_bar(position="fill", stat="identity") + theme_bw() + coord_flip() + 
      scale_fill_manual(values= color_levelPC) +
      ylab("Fraction of cells in each dataset") + xlab("celltype number") + theme(legend.position="top")
    
    p2 + p1 + plot_layout(widths = c(3,1))
  }
  plot_PC_integrated_donor_celltypes = function (srat) { 
    ## take an integrated Seurat object, plot distributions over orig.ident
    library(Seurat)
    library(patchwork)
    library(ggplot2)
    library(reshape2)
    library(RColorBrewer)
    
    
    count_table <- table(srat@meta.data$donor, srat@meta.data$celltype.full)
    count_mtx   <- as.data.frame.matrix(count_table)
    count_mtx$celltype <- rownames(count_mtx)
    melt_mtx    <- melt(count_mtx)
    melt_mtx$celltype <- as.factor(melt_mtx$celltype)
    
    celltype_size   <- aggregate(value ~ celltype, data = melt_mtx, FUN = sum)
    
    sorted_labels <- levels(srat$donor)
    celltype_size$celltype <- factor(celltype_size$celltype,levels = sorted_labels)
    melt_mtx$celltype <- factor(melt_mtx$celltype,levels = sorted_labels)
    
    colnames(melt_mtx)[2] <- "dataset"
    
    
    p1 <- ggplot(celltype_size, aes(y= celltype,x = value)) + geom_bar(position="dodge", stat="identity",fill = "grey60") + 
      theme_bw() + scale_x_log10() + xlab("Cells per celltype, log10 scale") + ylab("")
    p2 <- ggplot(melt_mtx,aes(x=celltype,y=value,fill=dataset)) + 
      geom_bar(position="fill", stat="identity") + theme_bw() + coord_flip() + 
      scale_fill_manual(values= color_levelPC) +
      ylab("Fraction of cells in each dataset") + xlab("celltype number") + theme(legend.position="top")
    
    p2 + p1 + plot_layout(widths = c(3,1))
  }
  ##
  #TEC
  plot_TEC_integrated_group_clusters = function (srat) { 
    ## take an integrated Seurat object, plot distributions over orig.ident
    library(Seurat)
    library(patchwork)
    library(ggplot2)
    library(reshape2)
    library(RColorBrewer)
    
    
    count_table <- table(srat@meta.data$group, srat@meta.data$celltype.full)
    count_mtx   <- as.data.frame.matrix(count_table)
    count_mtx$celltype <- rownames(count_mtx)
    melt_mtx    <- melt(count_mtx)
    melt_mtx$celltype <- as.factor(melt_mtx$celltype)
    
    celltype_size   <- aggregate(value ~ celltype, data = melt_mtx, FUN = sum)
    
    sorted_labels <- levels(srat$group)
    celltype_size$celltype <- factor(celltype_size$celltype,levels = sorted_labels)
    melt_mtx$celltype <- factor(melt_mtx$celltype,levels = sorted_labels)
    
    colnames(melt_mtx)[2] <- "dataset"
    
    
    p1 <- ggplot(celltype_size, aes(y= celltype,x = value)) + geom_bar(position="dodge", stat="identity",fill = "grey60") + 
      theme_bw() + scale_x_log10() + xlab("Cells per celltype, log10 scale") + ylab("")
    p2 <- ggplot(melt_mtx,aes(x=celltype,y=value,fill=dataset)) + 
      geom_bar(position="fill", stat="identity") + theme_bw() + coord_flip() + 
      scale_fill_manual(values= color_levelTEC) +
      ylab("Fraction of cells in each dataset") + xlab("celltype number") + theme(legend.position="top")
    
    p2 + p1 + plot_layout(widths = c(3,1))
  }
  plot_TEC_integrated_Donor_ID_clusters = function (srat) { 
    ## take an integrated Seurat object, plot distributions over orig.ident
    library(Seurat)
    library(patchwork)
    library(ggplot2)
    library(reshape2)
    library(RColorBrewer)
    
    
    count_table <- table(srat@meta.data$Donor_ID, srat@meta.data$celltype.full)
    count_mtx   <- as.data.frame.matrix(count_table)
    count_mtx$celltype <- rownames(count_mtx)
    melt_mtx    <- melt(count_mtx)
    melt_mtx$celltype <- as.factor(melt_mtx$celltype)
    
    celltype_size   <- aggregate(value ~ celltype, data = melt_mtx, FUN = sum)
    
    sorted_labels <- levels(srat$Donor_ID)
    celltype_size$celltype <- factor(celltype_size$celltype,levels = sorted_labels)
    melt_mtx$celltype <- factor(melt_mtx$celltype,levels = sorted_labels)
    
    colnames(melt_mtx)[2] <- "dataset"
    
    
    p1 <- ggplot(celltype_size, aes(y= celltype,x = value)) + geom_bar(position="dodge", stat="identity",fill = "grey60") + 
      theme_bw() + scale_x_log10() + xlab("Cells per celltype, log10 scale") + ylab("")
    p2 <- ggplot(melt_mtx,aes(x=celltype,y=value,fill=dataset)) + 
      geom_bar(position="fill", stat="identity") + theme_bw() + coord_flip() + 
      scale_fill_manual(values= color_levelTEC)+
      ylab("Fraction of cells in each dataset") + xlab("celltype number") + theme(legend.position="top")
    
    p2 + p1 + plot_layout(widths = c(3,1))
  }
  plot_TEC_integrated_condition_celltypes = function (srat) { 
    ## take an integrated Seurat object, plot distributions over orig.ident
    library(Seurat)
    library(patchwork)
    library(ggplot2)
    library(reshape2)
    library(RColorBrewer)
    
    
    count_table <- table(srat@meta.data$condition, srat@meta.data$celltype.full)
    count_mtx   <- as.data.frame.matrix(count_table)
    count_mtx$celltype <- rownames(count_mtx)
    melt_mtx    <- melt(count_mtx)
    melt_mtx$celltype <- as.factor(melt_mtx$celltype)
    
    celltype_size   <- aggregate(value ~ celltype, data = melt_mtx, FUN = sum)
    
    sorted_labels <- levels(srat$condition)
    celltype_size$celltype <- factor(celltype_size$celltype,levels = sorted_labels)
    melt_mtx$celltype <- factor(melt_mtx$celltype,levels = sorted_labels)
    
    colnames(melt_mtx)[2] <- "dataset"
    
    
    p1 <- ggplot(celltype_size, aes(y= celltype,x = value)) + geom_bar(position="dodge", stat="identity",fill = "grey60") + 
      theme_bw() + scale_x_log10() + xlab("Cells per celltype, log10 scale") + ylab("")
    p2 <- ggplot(melt_mtx,aes(x=celltype,y=value,fill=dataset)) + 
      geom_bar(position="fill", stat="identity") + theme_bw() + coord_flip() + 
      scale_fill_manual(values= color_levelTEC) +
      ylab("Fraction of cells in each dataset") + xlab("celltype number") + theme(legend.position="top")
    
    p2 + p1 + plot_layout(widths = c(3,1))
  }
  plot_TEC_integrated_classification_celltypes = function (srat) { 
    ## take an integrated Seurat object, plot distributions over orig.ident
    library(Seurat)
    library(patchwork)
    library(ggplot2)
    library(reshape2)
    library(RColorBrewer)
    
    
    count_table <- table(srat@meta.data$classification, srat@meta.data$celltype.full)
    count_mtx   <- as.data.frame.matrix(count_table)
    count_mtx$celltype <- rownames(count_mtx)
    melt_mtx    <- melt(count_mtx)
    melt_mtx$celltype <- as.factor(melt_mtx$celltype)
    
    celltype_size   <- aggregate(value ~ celltype, data = melt_mtx, FUN = sum)
    
    sorted_labels <- levels(srat$classification)
    celltype_size$celltype <- factor(celltype_size$celltype,levels = sorted_labels)
    melt_mtx$celltype <- factor(melt_mtx$celltype,levels = sorted_labels)
    
    colnames(melt_mtx)[2] <- "dataset"
    
    
    p1 <- ggplot(celltype_size, aes(y= celltype,x = value)) + geom_bar(position="dodge", stat="identity",fill = "grey60") + 
      theme_bw() + scale_x_log10() + xlab("Cells per celltype, log10 scale") + ylab("")
    p2 <- ggplot(melt_mtx,aes(x=celltype,y=value,fill=dataset)) + 
      geom_bar(position="fill", stat="identity") + theme_bw() + coord_flip() + 
      scale_fill_manual(values= color_levelTEC) +
      ylab("Fraction of cells in each dataset") + xlab("celltype number") + theme(legend.position="top")
    
    p2 + p1 + plot_layout(widths = c(3,1))
  }
  plot_TEC_integrated_tissue_celltypes = function (srat) { 
    ## take an integrated Seurat object, plot distributions over orig.ident
    library(Seurat)
    library(patchwork)
    library(ggplot2)
    library(reshape2)
    library(RColorBrewer)
    
    
    count_table <- table(srat@meta.data$tissue, srat@meta.data$celltype.full)
    count_mtx   <- as.data.frame.matrix(count_table)
    count_mtx$celltype <- rownames(count_mtx)
    melt_mtx    <- melt(count_mtx)
    melt_mtx$celltype <- as.factor(melt_mtx$celltype)
    
    celltype_size   <- aggregate(value ~ celltype, data = melt_mtx, FUN = sum)
    
    sorted_labels <- levels(srat$tissue)
    celltype_size$celltype <- factor(celltype_size$celltype,levels = sorted_labels)
    melt_mtx$celltype <- factor(melt_mtx$celltype,levels = sorted_labels)
    
    colnames(melt_mtx)[2] <- "dataset"
    
    
    p1 <- ggplot(celltype_size, aes(y= celltype,x = value)) + geom_bar(position="dodge", stat="identity",fill = "grey60") + 
      theme_bw() + scale_x_log10() + xlab("Cells per celltype, log10 scale") + ylab("")
    p2 <- ggplot(melt_mtx,aes(x=celltype,y=value,fill=dataset)) + 
      geom_bar(position="fill", stat="identity") + theme_bw() + coord_flip() + 
      scale_fill_manual(values= color_levelTEC) +
      ylab("Fraction of cells in each dataset") + xlab("celltype number") + theme(legend.position="top")
    
    p2 + p1 + plot_layout(widths = c(3,1))
  }
  plot_TEC_integrated_development_celltypes = function (srat) { 
    ## take an integrated Seurat object, plot distributions over orig.ident
    library(Seurat)
    library(patchwork)
    library(ggplot2)
    library(reshape2)
    library(RColorBrewer)
    
    
    count_table <- table(srat@meta.data$development, srat@meta.data$celltype.full)
    count_mtx   <- as.data.frame.matrix(count_table)
    count_mtx$celltype <- rownames(count_mtx)
    melt_mtx    <- melt(count_mtx)
    melt_mtx$celltype <- as.factor(melt_mtx$celltype)
    
    celltype_size   <- aggregate(value ~ celltype, data = melt_mtx, FUN = sum)
    
    sorted_labels <- levels(srat$development)
    celltype_size$celltype <- factor(celltype_size$celltype,levels = sorted_labels)
    melt_mtx$celltype <- factor(melt_mtx$celltype,levels = sorted_labels)
    
    colnames(melt_mtx)[2] <- "dataset"
    
    
    p1 <- ggplot(celltype_size, aes(y= celltype,x = value)) + geom_bar(position="dodge", stat="identity",fill = "grey60") + 
      theme_bw() + scale_x_log10() + xlab("Cells per celltype, log10 scale") + ylab("")
    p2 <- ggplot(melt_mtx,aes(x=celltype,y=value,fill=dataset)) + 
      geom_bar(position="fill", stat="identity") + theme_bw() + coord_flip() + 
      scale_fill_manual(values= color_levelTEC) +
      ylab("Fraction of cells in each dataset") + xlab("celltype number") + theme(legend.position="top")
    
    p2 + p1 + plot_layout(widths = c(3,1))
  }
  plot_TEC_integrated_donor_celltypes = function (srat) { 
    ## take an integrated Seurat object, plot distributions over orig.ident
    library(Seurat)
    library(patchwork)
    library(ggplot2)
    library(reshape2)
    library(RColorBrewer)
    
    
    count_table <- table(srat@meta.data$donor, srat@meta.data$celltype.full)
    count_mtx   <- as.data.frame.matrix(count_table)
    count_mtx$celltype <- rownames(count_mtx)
    melt_mtx    <- melt(count_mtx)
    melt_mtx$celltype <- as.factor(melt_mtx$celltype)
    
    celltype_size   <- aggregate(value ~ celltype, data = melt_mtx, FUN = sum)
    
    sorted_labels <- levels(srat$donor)
    celltype_size$celltype <- factor(celltype_size$celltype,levels = sorted_labels)
    melt_mtx$celltype <- factor(melt_mtx$celltype,levels = sorted_labels)
    
    colnames(melt_mtx)[2] <- "dataset"
    
    
    p1 <- ggplot(celltype_size, aes(y= celltype,x = value)) + geom_bar(position="dodge", stat="identity",fill = "grey60") + 
      theme_bw() + scale_x_log10() + xlab("Cells per celltype, log10 scale") + ylab("")
    p2 <- ggplot(melt_mtx,aes(x=celltype,y=value,fill=dataset)) + 
      geom_bar(position="fill", stat="identity") + theme_bw() + coord_flip() + 
      scale_fill_manual(values= color_levelTEC) +
      ylab("Fraction of cells in each dataset") + xlab("celltype number") + theme(legend.position="top")
    
    p2 + p1 + plot_layout(widths = c(3,1))
  }
  ##
  #DC
  plot_DC_integrated_group_clusters = function (srat) { 
    ## take an integrated Seurat object, plot distributions over orig.ident
    library(Seurat)
    library(patchwork)
    library(ggplot2)
    library(reshape2)
    library(RColorBrewer)
    
    
    count_table <- table(srat@meta.data$group, srat@meta.data$celltype.full)
    count_mtx   <- as.data.frame.matrix(count_table)
    count_mtx$celltype <- rownames(count_mtx)
    melt_mtx    <- melt(count_mtx)
    melt_mtx$celltype <- as.factor(melt_mtx$celltype)
    
    celltype_size   <- aggregate(value ~ celltype, data = melt_mtx, FUN = sum)
    
    sorted_labels <- levels(srat$group)
    celltype_size$celltype <- factor(celltype_size$celltype,levels = sorted_labels)
    melt_mtx$celltype <- factor(melt_mtx$celltype,levels = sorted_labels)
    
    colnames(melt_mtx)[2] <- "dataset"
    
    
    p1 <- ggplot(celltype_size, aes(y= celltype,x = value)) + geom_bar(position="dodge", stat="identity",fill = "grey60") + 
      theme_bw() + scale_x_log10() + xlab("Cells per celltype, log10 scale") + ylab("")
    p2 <- ggplot(melt_mtx,aes(x=celltype,y=value,fill=dataset)) + 
      geom_bar(position="fill", stat="identity") + theme_bw() + coord_flip() + 
      scale_fill_manual(values= color_levelDC) +
      ylab("Fraction of cells in each dataset") + xlab("celltype number") + theme(legend.position="top")
    
    p2 + p1 + plot_layout(widths = c(3,1))
  }
  plot_DC_integrated_Donor_ID_clusters = function (srat) { 
    ## take an integrated Seurat object, plot distributions over orig.ident
    library(Seurat)
    library(patchwork)
    library(ggplot2)
    library(reshape2)
    library(RColorBrewer)
    
    
    count_table <- table(srat@meta.data$Donor_ID, srat@meta.data$celltype.full)
    count_mtx   <- as.data.frame.matrix(count_table)
    count_mtx$celltype <- rownames(count_mtx)
    melt_mtx    <- melt(count_mtx)
    melt_mtx$celltype <- as.factor(melt_mtx$celltype)
    
    celltype_size   <- aggregate(value ~ celltype, data = melt_mtx, FUN = sum)
    
    sorted_labels <- levels(srat$Donor_ID)
    celltype_size$celltype <- factor(celltype_size$celltype,levels = sorted_labels)
    melt_mtx$celltype <- factor(melt_mtx$celltype,levels = sorted_labels)
    
    colnames(melt_mtx)[2] <- "dataset"
    
    
    p1 <- ggplot(celltype_size, aes(y= celltype,x = value)) + geom_bar(position="dodge", stat="identity",fill = "grey60") + 
      theme_bw() + scale_x_log10() + xlab("Cells per celltype, log10 scale") + ylab("")
    p2 <- ggplot(melt_mtx,aes(x=celltype,y=value,fill=dataset)) + 
      geom_bar(position="fill", stat="identity") + theme_bw() + coord_flip() + 
      scale_fill_manual(values= color_levelDC)+
      ylab("Fraction of cells in each dataset") + xlab("celltype number") + theme(legend.position="top")
    
    p2 + p1 + plot_layout(widths = c(3,1))
  }
  plot_DC_integrated_condition_celltypes = function (srat) { 
    ## take an integrated Seurat object, plot distributions over orig.ident
    library(Seurat)
    library(patchwork)
    library(ggplot2)
    library(reshape2)
    library(RColorBrewer)
    
    
    count_table <- table(srat@meta.data$condition, srat@meta.data$celltype.full)
    count_mtx   <- as.data.frame.matrix(count_table)
    count_mtx$celltype <- rownames(count_mtx)
    melt_mtx    <- melt(count_mtx)
    melt_mtx$celltype <- as.factor(melt_mtx$celltype)
    
    celltype_size   <- aggregate(value ~ celltype, data = melt_mtx, FUN = sum)
    
    sorted_labels <- levels(srat$condition)
    celltype_size$celltype <- factor(celltype_size$celltype,levels = sorted_labels)
    melt_mtx$celltype <- factor(melt_mtx$celltype,levels = sorted_labels)
    
    colnames(melt_mtx)[2] <- "dataset"
    
    
    p1 <- ggplot(celltype_size, aes(y= celltype,x = value)) + geom_bar(position="dodge", stat="identity",fill = "grey60") + 
      theme_bw() + scale_x_log10() + xlab("Cells per celltype, log10 scale") + ylab("")
    p2 <- ggplot(melt_mtx,aes(x=celltype,y=value,fill=dataset)) + 
      geom_bar(position="fill", stat="identity") + theme_bw() + coord_flip() + 
      scale_fill_manual(values= color_levelDC) +
      ylab("Fraction of cells in each dataset") + xlab("celltype number") + theme(legend.position="top")
    
    p2 + p1 + plot_layout(widths = c(3,1))
  }
  plot_DC_integrated_classification_celltypes = function (srat) { 
    ## take an integrated Seurat object, plot distributions over orig.ident
    library(Seurat)
    library(patchwork)
    library(ggplot2)
    library(reshape2)
    library(RColorBrewer)
    
    
    count_table <- table(srat@meta.data$classification, srat@meta.data$celltype.full)
    count_mtx   <- as.data.frame.matrix(count_table)
    count_mtx$celltype <- rownames(count_mtx)
    melt_mtx    <- melt(count_mtx)
    melt_mtx$celltype <- as.factor(melt_mtx$celltype)
    
    celltype_size   <- aggregate(value ~ celltype, data = melt_mtx, FUN = sum)
    
    sorted_labels <- levels(srat$classification)
    celltype_size$celltype <- factor(celltype_size$celltype,levels = sorted_labels)
    melt_mtx$celltype <- factor(melt_mtx$celltype,levels = sorted_labels)
    
    colnames(melt_mtx)[2] <- "dataset"
    
    
    p1 <- ggplot(celltype_size, aes(y= celltype,x = value)) + geom_bar(position="dodge", stat="identity",fill = "grey60") + 
      theme_bw() + scale_x_log10() + xlab("Cells per celltype, log10 scale") + ylab("")
    p2 <- ggplot(melt_mtx,aes(x=celltype,y=value,fill=dataset)) + 
      geom_bar(position="fill", stat="identity") + theme_bw() + coord_flip() + 
      scale_fill_manual(values= color_levelDC) +
      ylab("Fraction of cells in each dataset") + xlab("celltype number") + theme(legend.position="top")
    
    p2 + p1 + plot_layout(widths = c(3,1))
  }
  plot_DC_integrated_tissue_celltypes = function (srat) { 
    ## take an integrated Seurat object, plot distributions over orig.ident
    library(Seurat)
    library(patchwork)
    library(ggplot2)
    library(reshape2)
    library(RColorBrewer)
    
    
    count_table <- table(srat@meta.data$tissue, srat@meta.data$celltype.full)
    count_mtx   <- as.data.frame.matrix(count_table)
    count_mtx$celltype <- rownames(count_mtx)
    melt_mtx    <- melt(count_mtx)
    melt_mtx$celltype <- as.factor(melt_mtx$celltype)
    
    celltype_size   <- aggregate(value ~ celltype, data = melt_mtx, FUN = sum)
    
    sorted_labels <- levels(srat$tissue)
    celltype_size$celltype <- factor(celltype_size$celltype,levels = sorted_labels)
    melt_mtx$celltype <- factor(melt_mtx$celltype,levels = sorted_labels)
    
    colnames(melt_mtx)[2] <- "dataset"
    
    
    p1 <- ggplot(celltype_size, aes(y= celltype,x = value)) + geom_bar(position="dodge", stat="identity",fill = "grey60") + 
      theme_bw() + scale_x_log10() + xlab("Cells per celltype, log10 scale") + ylab("")
    p2 <- ggplot(melt_mtx,aes(x=celltype,y=value,fill=dataset)) + 
      geom_bar(position="fill", stat="identity") + theme_bw() + coord_flip() + 
      scale_fill_manual(values= color_levelDC) +
      ylab("Fraction of cells in each dataset") + xlab("celltype number") + theme(legend.position="top")
    
    p2 + p1 + plot_layout(widths = c(3,1))
  }
  plot_DC_integrated_development_celltypes = function (srat) { 
    ## take an integrated Seurat object, plot distributions over orig.ident
    library(Seurat)
    library(patchwork)
    library(ggplot2)
    library(reshape2)
    library(RColorBrewer)
    
    
    count_table <- table(srat@meta.data$development, srat@meta.data$celltype.full)
    count_mtx   <- as.data.frame.matrix(count_table)
    count_mtx$celltype <- rownames(count_mtx)
    melt_mtx    <- melt(count_mtx)
    melt_mtx$celltype <- as.factor(melt_mtx$celltype)
    
    celltype_size   <- aggregate(value ~ celltype, data = melt_mtx, FUN = sum)
    
    sorted_labels <- levels(srat$development)
    celltype_size$celltype <- factor(celltype_size$celltype,levels = sorted_labels)
    melt_mtx$celltype <- factor(melt_mtx$celltype,levels = sorted_labels)
    
    colnames(melt_mtx)[2] <- "dataset"
    
    
    p1 <- ggplot(celltype_size, aes(y= celltype,x = value)) + geom_bar(position="dodge", stat="identity",fill = "grey60") + 
      theme_bw() + scale_x_log10() + xlab("Cells per celltype, log10 scale") + ylab("")
    p2 <- ggplot(melt_mtx,aes(x=celltype,y=value,fill=dataset)) + 
      geom_bar(position="fill", stat="identity") + theme_bw() + coord_flip() + 
      scale_fill_manual(values= color_levelDC) +
      ylab("Fraction of cells in each dataset") + xlab("celltype number") + theme(legend.position="top")
    
    p2 + p1 + plot_layout(widths = c(3,1))
  }
  plot_DC_integrated_donor_celltypes = function (srat) { 
    ## take an integrated Seurat object, plot distributions over orig.ident
    library(Seurat)
    library(patchwork)
    library(ggplot2)
    library(reshape2)
    library(RColorBrewer)
    
    
    count_table <- table(srat@meta.data$donor, srat@meta.data$celltype.full)
    count_mtx   <- as.data.frame.matrix(count_table)
    count_mtx$celltype <- rownames(count_mtx)
    melt_mtx    <- melt(count_mtx)
    melt_mtx$celltype <- as.factor(melt_mtx$celltype)
    
    celltype_size   <- aggregate(value ~ celltype, data = melt_mtx, FUN = sum)
    
    sorted_labels <- levels(srat$donor)
    celltype_size$celltype <- factor(celltype_size$celltype,levels = sorted_labels)
    melt_mtx$celltype <- factor(melt_mtx$celltype,levels = sorted_labels)
    
    colnames(melt_mtx)[2] <- "dataset"
    
    
    p1 <- ggplot(celltype_size, aes(y= celltype,x = value)) + geom_bar(position="dodge", stat="identity",fill = "grey60") + 
      theme_bw() + scale_x_log10() + xlab("Cells per celltype, log10 scale") + ylab("")
    p2 <- ggplot(melt_mtx,aes(x=celltype,y=value,fill=dataset)) + 
      geom_bar(position="fill", stat="identity") + theme_bw() + coord_flip() + 
      scale_fill_manual(values= color_levelDC) +
      ylab("Fraction of cells in each dataset") + xlab("celltype number") + theme(legend.position="top")
    
    p2 + p1 + plot_layout(widths = c(3,1))
  }
  ##
  #MAC/Mono
  plot_MAC_Mono_integrated_group_clusters = function (srat) { 
    ## take an integrated Seurat object, plot distributions over orig.ident
    library(Seurat)
    library(patchwork)
    library(ggplot2)
    library(reshape2)
    library(RColorBrewer)
    
    
    count_table <- table(srat@meta.data$group, srat@meta.data$celltype.full)
    count_mtx   <- as.data.frame.matrix(count_table)
    count_mtx$celltype <- rownames(count_mtx)
    melt_mtx    <- melt(count_mtx)
    melt_mtx$celltype <- as.factor(melt_mtx$celltype)
    
    celltype_size   <- aggregate(value ~ celltype, data = melt_mtx, FUN = sum)
    
    sorted_labels <- levels(srat$group)
    celltype_size$celltype <- factor(celltype_size$celltype,levels = sorted_labels)
    melt_mtx$celltype <- factor(melt_mtx$celltype,levels = sorted_labels)
    
    colnames(melt_mtx)[2] <- "dataset"
    
    
    p1 <- ggplot(celltype_size, aes(y= celltype,x = value)) + geom_bar(position="dodge", stat="identity",fill = "grey60") + 
      theme_bw() + scale_x_log10() + xlab("Cells per celltype, log10 scale") + ylab("")
    p2 <- ggplot(melt_mtx,aes(x=celltype,y=value,fill=dataset)) + 
      geom_bar(position="fill", stat="identity") + theme_bw() + coord_flip() + 
      scale_fill_manual(values= color_levelMAC_Mono) +
      ylab("Fraction of cells in each dataset") + xlab("celltype number") + theme(legend.position="top")
    
    p2 + p1 + plot_layout(widths = c(3,1))
  }
  plot_MAC_Mono_integrated_Donor_ID_clusters = function (srat) { 
    ## take an integrated Seurat object, plot distributions over orig.ident
    library(Seurat)
    library(patchwork)
    library(ggplot2)
    library(reshape2)
    library(RColorBrewer)
    
    
    count_table <- table(srat@meta.data$Donor_ID, srat@meta.data$celltype.full)
    count_mtx   <- as.data.frame.matrix(count_table)
    count_mtx$celltype <- rownames(count_mtx)
    melt_mtx    <- melt(count_mtx)
    melt_mtx$celltype <- as.factor(melt_mtx$celltype)
    
    celltype_size   <- aggregate(value ~ celltype, data = melt_mtx, FUN = sum)
    
    sorted_labels <- levels(srat$Donor_ID)
    celltype_size$celltype <- factor(celltype_size$celltype,levels = sorted_labels)
    melt_mtx$celltype <- factor(melt_mtx$celltype,levels = sorted_labels)
    
    colnames(melt_mtx)[2] <- "dataset"
    
    
    p1 <- ggplot(celltype_size, aes(y= celltype,x = value)) + geom_bar(position="dodge", stat="identity",fill = "grey60") + 
      theme_bw() + scale_x_log10() + xlab("Cells per celltype, log10 scale") + ylab("")
    p2 <- ggplot(melt_mtx,aes(x=celltype,y=value,fill=dataset)) + 
      geom_bar(position="fill", stat="identity") + theme_bw() + coord_flip() + 
      scale_fill_manual(values= color_levelMAC_Mono)+
      ylab("Fraction of cells in each dataset") + xlab("celltype number") + theme(legend.position="top")
    
    p2 + p1 + plot_layout(widths = c(3,1))
  }
  plot_MAC_Mono_integrated_condition_celltypes = function (srat) { 
    ## take an integrated Seurat object, plot distributions over orig.ident
    library(Seurat)
    library(patchwork)
    library(ggplot2)
    library(reshape2)
    library(RColorBrewer)
    
    
    count_table <- table(srat@meta.data$condition, srat@meta.data$celltype.full)
    count_mtx   <- as.data.frame.matrix(count_table)
    count_mtx$celltype <- rownames(count_mtx)
    melt_mtx    <- melt(count_mtx)
    melt_mtx$celltype <- as.factor(melt_mtx$celltype)
    
    celltype_size   <- aggregate(value ~ celltype, data = melt_mtx, FUN = sum)
    
    sorted_labels <- levels(srat$condition)
    celltype_size$celltype <- factor(celltype_size$celltype,levels = sorted_labels)
    melt_mtx$celltype <- factor(melt_mtx$celltype,levels = sorted_labels)
    
    colnames(melt_mtx)[2] <- "dataset"
    
    
    p1 <- ggplot(celltype_size, aes(y= celltype,x = value)) + geom_bar(position="dodge", stat="identity",fill = "grey60") + 
      theme_bw() + scale_x_log10() + xlab("Cells per celltype, log10 scale") + ylab("")
    p2 <- ggplot(melt_mtx,aes(x=celltype,y=value,fill=dataset)) + 
      geom_bar(position="fill", stat="identity") + theme_bw() + coord_flip() + 
      scale_fill_manual(values= color_levelMAC_Mono) +
      ylab("Fraction of cells in each dataset") + xlab("celltype number") + theme(legend.position="top")
    
    p2 + p1 + plot_layout(widths = c(3,1))
  }
  plot_MAC_Mono_integrated_classification_celltypes = function (srat) { 
    ## take an integrated Seurat object, plot distributions over orig.ident
    library(Seurat)
    library(patchwork)
    library(ggplot2)
    library(reshape2)
    library(RColorBrewer)
    
    
    count_table <- table(srat@meta.data$classification, srat@meta.data$celltype.full)
    count_mtx   <- as.data.frame.matrix(count_table)
    count_mtx$celltype <- rownames(count_mtx)
    melt_mtx    <- melt(count_mtx)
    melt_mtx$celltype <- as.factor(melt_mtx$celltype)
    
    celltype_size   <- aggregate(value ~ celltype, data = melt_mtx, FUN = sum)
    
    sorted_labels <- levels(srat$classification)
    celltype_size$celltype <- factor(celltype_size$celltype,levels = sorted_labels)
    melt_mtx$celltype <- factor(melt_mtx$celltype,levels = sorted_labels)
    
    colnames(melt_mtx)[2] <- "dataset"
    
    
    p1 <- ggplot(celltype_size, aes(y= celltype,x = value)) + geom_bar(position="dodge", stat="identity",fill = "grey60") + 
      theme_bw() + scale_x_log10() + xlab("Cells per celltype, log10 scale") + ylab("")
    p2 <- ggplot(melt_mtx,aes(x=celltype,y=value,fill=dataset)) + 
      geom_bar(position="fill", stat="identity") + theme_bw() + coord_flip() + 
      scale_fill_manual(values= color_levelMAC_Mono) +
      ylab("Fraction of cells in each dataset") + xlab("celltype number") + theme(legend.position="top")
    
    p2 + p1 + plot_layout(widths = c(3,1))
  }
  plot_MAC_Mono_integrated_tissue_celltypes = function (srat) { 
    ## take an integrated Seurat object, plot distributions over orig.ident
    library(Seurat)
    library(patchwork)
    library(ggplot2)
    library(reshape2)
    library(RColorBrewer)
    
    
    count_table <- table(srat@meta.data$tissue, srat@meta.data$celltype.full)
    count_mtx   <- as.data.frame.matrix(count_table)
    count_mtx$celltype <- rownames(count_mtx)
    melt_mtx    <- melt(count_mtx)
    melt_mtx$celltype <- as.factor(melt_mtx$celltype)
    
    celltype_size   <- aggregate(value ~ celltype, data = melt_mtx, FUN = sum)
    
    sorted_labels <- levels(srat$tissue)
    celltype_size$celltype <- factor(celltype_size$celltype,levels = sorted_labels)
    melt_mtx$celltype <- factor(melt_mtx$celltype,levels = sorted_labels)
    
    colnames(melt_mtx)[2] <- "dataset"
    
    
    p1 <- ggplot(celltype_size, aes(y= celltype,x = value)) + geom_bar(position="dodge", stat="identity",fill = "grey60") + 
      theme_bw() + scale_x_log10() + xlab("Cells per celltype, log10 scale") + ylab("")
    p2 <- ggplot(melt_mtx,aes(x=celltype,y=value,fill=dataset)) + 
      geom_bar(position="fill", stat="identity") + theme_bw() + coord_flip() + 
      scale_fill_manual(values= color_levelMAC_Mono) +
      ylab("Fraction of cells in each dataset") + xlab("celltype number") + theme(legend.position="top")
    
    p2 + p1 + plot_layout(widths = c(3,1))
  }
  plot_MAC_Mono_integrated_development_celltypes = function (srat) { 
    ## take an integrated Seurat object, plot distributions over orig.ident
    library(Seurat)
    library(patchwork)
    library(ggplot2)
    library(reshape2)
    library(RColorBrewer)
    
    
    count_table <- table(srat@meta.data$development, srat@meta.data$celltype.full)
    count_mtx   <- as.data.frame.matrix(count_table)
    count_mtx$celltype <- rownames(count_mtx)
    melt_mtx    <- melt(count_mtx)
    melt_mtx$celltype <- as.factor(melt_mtx$celltype)
    
    celltype_size   <- aggregate(value ~ celltype, data = melt_mtx, FUN = sum)
    
    sorted_labels <- levels(srat$development)
    celltype_size$celltype <- factor(celltype_size$celltype,levels = sorted_labels)
    melt_mtx$celltype <- factor(melt_mtx$celltype,levels = sorted_labels)
    
    colnames(melt_mtx)[2] <- "dataset"
    
    
    p1 <- ggplot(celltype_size, aes(y= celltype,x = value)) + geom_bar(position="dodge", stat="identity",fill = "grey60") + 
      theme_bw() + scale_x_log10() + xlab("Cells per celltype, log10 scale") + ylab("")
    p2 <- ggplot(melt_mtx,aes(x=celltype,y=value,fill=dataset)) + 
      geom_bar(position="fill", stat="identity") + theme_bw() + coord_flip() + 
      scale_fill_manual(values= color_levelMAC_Mono) +
      ylab("Fraction of cells in each dataset") + xlab("celltype number") + theme(legend.position="top")
    
    p2 + p1 + plot_layout(widths = c(3,1))
  }
  plot_MAC_Mono_integrated_donor_celltypes = function (srat) { 
    ## take an integrated Seurat object, plot distributions over orig.ident
    library(Seurat)
    library(patchwork)
    library(ggplot2)
    library(reshape2)
    library(RColorBrewer)
    
    
    count_table <- table(srat@meta.data$donor, srat@meta.data$celltype.full)
    count_mtx   <- as.data.frame.matrix(count_table)
    count_mtx$celltype <- rownames(count_mtx)
    melt_mtx    <- melt(count_mtx)
    melt_mtx$celltype <- as.factor(melt_mtx$celltype)
    
    celltype_size   <- aggregate(value ~ celltype, data = melt_mtx, FUN = sum)
    
    sorted_labels <- levels(srat$donor)
    celltype_size$celltype <- factor(celltype_size$celltype,levels = sorted_labels)
    melt_mtx$celltype <- factor(melt_mtx$celltype,levels = sorted_labels)
    
    colnames(melt_mtx)[2] <- "dataset"
    
    
    p1 <- ggplot(celltype_size, aes(y= celltype,x = value)) + geom_bar(position="dodge", stat="identity",fill = "grey60") + 
      theme_bw() + scale_x_log10() + xlab("Cells per celltype, log10 scale") + ylab("")
    p2 <- ggplot(melt_mtx,aes(x=celltype,y=value,fill=dataset)) + 
      geom_bar(position="fill", stat="identity") + theme_bw() + coord_flip() + 
      scale_fill_manual(values= color_levelMAC_Mono) +
      ylab("Fraction of cells in each dataset") + xlab("celltype number") + theme(legend.position="top")
    
    p2 + p1 + plot_layout(widths = c(3,1))
  }
  ##
  #EC
  plot_EC_integrated_group_clusters = function (srat) { 
    ## take an integrated Seurat object, plot distributions over orig.ident
    library(Seurat)
    library(patchwork)
    library(ggplot2)
    library(reshape2)
    library(RColorBrewer)
    
    
    count_table <- table(srat@meta.data$group, srat@meta.data$celltype.full)
    count_mtx   <- as.data.frame.matrix(count_table)
    count_mtx$celltype <- rownames(count_mtx)
    melt_mtx    <- melt(count_mtx)
    melt_mtx$celltype <- as.factor(melt_mtx$celltype)
    
    celltype_size   <- aggregate(value ~ celltype, data = melt_mtx, FUN = sum)
    
    sorted_labels <- levels(srat$group)
    celltype_size$celltype <- factor(celltype_size$celltype,levels = sorted_labels)
    melt_mtx$celltype <- factor(melt_mtx$celltype,levels = sorted_labels)
    
    colnames(melt_mtx)[2] <- "dataset"
    
    
    p1 <- ggplot(celltype_size, aes(y= celltype,x = value)) + geom_bar(position="dodge", stat="identity",fill = "grey60") + 
      theme_bw() + scale_x_log10() + xlab("Cells per celltype, log10 scale") + ylab("")
    p2 <- ggplot(melt_mtx,aes(x=celltype,y=value,fill=dataset)) + 
      geom_bar(position="fill", stat="identity") + theme_bw() + coord_flip() + 
      scale_fill_manual(values= color_levelEC) +
      ylab("Fraction of cells in each dataset") + xlab("celltype number") + theme(legend.position="top")
    
    p2 + p1 + plot_layout(widths = c(3,1))
  }
  plot_EC_integrated_Donor_ID_clusters = function (srat) { 
    ## take an integrated Seurat object, plot distributions over orig.ident
    library(Seurat)
    library(patchwork)
    library(ggplot2)
    library(reshape2)
    library(RColorBrewer)
    
    
    count_table <- table(srat@meta.data$Donor_ID, srat@meta.data$celltype.full)
    count_mtx   <- as.data.frame.matrix(count_table)
    count_mtx$celltype <- rownames(count_mtx)
    melt_mtx    <- melt(count_mtx)
    melt_mtx$celltype <- as.factor(melt_mtx$celltype)
    
    celltype_size   <- aggregate(value ~ celltype, data = melt_mtx, FUN = sum)
    
    sorted_labels <- levels(srat$Donor_ID)
    celltype_size$celltype <- factor(celltype_size$celltype,levels = sorted_labels)
    melt_mtx$celltype <- factor(melt_mtx$celltype,levels = sorted_labels)
    
    colnames(melt_mtx)[2] <- "dataset"
    
    
    p1 <- ggplot(celltype_size, aes(y= celltype,x = value)) + geom_bar(position="dodge", stat="identity",fill = "grey60") + 
      theme_bw() + scale_x_log10() + xlab("Cells per celltype, log10 scale") + ylab("")
    p2 <- ggplot(melt_mtx,aes(x=celltype,y=value,fill=dataset)) + 
      geom_bar(position="fill", stat="identity") + theme_bw() + coord_flip() + 
      scale_fill_manual(values= color_levelEC)+
      ylab("Fraction of cells in each dataset") + xlab("celltype number") + theme(legend.position="top")
    
    p2 + p1 + plot_layout(widths = c(3,1))
  }
  plot_EC_integrated_condition_celltypes = function (srat) { 
    ## take an integrated Seurat object, plot distributions over orig.ident
    library(Seurat)
    library(patchwork)
    library(ggplot2)
    library(reshape2)
    library(RColorBrewer)
    
    
    count_table <- table(srat@meta.data$condition, srat@meta.data$celltype.full)
    count_mtx   <- as.data.frame.matrix(count_table)
    count_mtx$celltype <- rownames(count_mtx)
    melt_mtx    <- melt(count_mtx)
    melt_mtx$celltype <- as.factor(melt_mtx$celltype)
    
    celltype_size   <- aggregate(value ~ celltype, data = melt_mtx, FUN = sum)
    
    sorted_labels <- levels(srat$condition)
    celltype_size$celltype <- factor(celltype_size$celltype,levels = sorted_labels)
    melt_mtx$celltype <- factor(melt_mtx$celltype,levels = sorted_labels)
    
    colnames(melt_mtx)[2] <- "dataset"
    
    
    p1 <- ggplot(celltype_size, aes(y= celltype,x = value)) + geom_bar(position="dodge", stat="identity",fill = "grey60") + 
      theme_bw() + scale_x_log10() + xlab("Cells per celltype, log10 scale") + ylab("")
    p2 <- ggplot(melt_mtx,aes(x=celltype,y=value,fill=dataset)) + 
      geom_bar(position="fill", stat="identity") + theme_bw() + coord_flip() + 
      scale_fill_manual(values= color_levelEC) +
      ylab("Fraction of cells in each dataset") + xlab("celltype number") + theme(legend.position="top")
    
    p2 + p1 + plot_layout(widths = c(3,1))
  }
  plot_EC_integrated_classification_celltypes = function (srat) { 
    ## take an integrated Seurat object, plot distributions over orig.ident
    library(Seurat)
    library(patchwork)
    library(ggplot2)
    library(reshape2)
    library(RColorBrewer)
    
    
    count_table <- table(srat@meta.data$classification, srat@meta.data$celltype.full)
    count_mtx   <- as.data.frame.matrix(count_table)
    count_mtx$celltype <- rownames(count_mtx)
    melt_mtx    <- melt(count_mtx)
    melt_mtx$celltype <- as.factor(melt_mtx$celltype)
    
    celltype_size   <- aggregate(value ~ celltype, data = melt_mtx, FUN = sum)
    
    sorted_labels <- levels(srat$classification)
    celltype_size$celltype <- factor(celltype_size$celltype,levels = sorted_labels)
    melt_mtx$celltype <- factor(melt_mtx$celltype,levels = sorted_labels)
    
    colnames(melt_mtx)[2] <- "dataset"
    
    
    p1 <- ggplot(celltype_size, aes(y= celltype,x = value)) + geom_bar(position="dodge", stat="identity",fill = "grey60") + 
      theme_bw() + scale_x_log10() + xlab("Cells per celltype, log10 scale") + ylab("")
    p2 <- ggplot(melt_mtx,aes(x=celltype,y=value,fill=dataset)) + 
      geom_bar(position="fill", stat="identity") + theme_bw() + coord_flip() + 
      scale_fill_manual(values= color_levelEC) +
      ylab("Fraction of cells in each dataset") + xlab("celltype number") + theme(legend.position="top")
    
    p2 + p1 + plot_layout(widths = c(3,1))
  }
  plot_EC_integrated_tissue_celltypes = function (srat) { 
    ## take an integrated Seurat object, plot distributions over orig.ident
    library(Seurat)
    library(patchwork)
    library(ggplot2)
    library(reshape2)
    library(RColorBrewer)
    
    
    count_table <- table(srat@meta.data$tissue, srat@meta.data$celltype.full)
    count_mtx   <- as.data.frame.matrix(count_table)
    count_mtx$celltype <- rownames(count_mtx)
    melt_mtx    <- melt(count_mtx)
    melt_mtx$celltype <- as.factor(melt_mtx$celltype)
    
    celltype_size   <- aggregate(value ~ celltype, data = melt_mtx, FUN = sum)
    
    sorted_labels <- levels(srat$tissue)
    celltype_size$celltype <- factor(celltype_size$celltype,levels = sorted_labels)
    melt_mtx$celltype <- factor(melt_mtx$celltype,levels = sorted_labels)
    
    colnames(melt_mtx)[2] <- "dataset"
    
    
    p1 <- ggplot(celltype_size, aes(y= celltype,x = value)) + geom_bar(position="dodge", stat="identity",fill = "grey60") + 
      theme_bw() + scale_x_log10() + xlab("Cells per celltype, log10 scale") + ylab("")
    p2 <- ggplot(melt_mtx,aes(x=celltype,y=value,fill=dataset)) + 
      geom_bar(position="fill", stat="identity") + theme_bw() + coord_flip() + 
      scale_fill_manual(values= color_levelEC) +
      ylab("Fraction of cells in each dataset") + xlab("celltype number") + theme(legend.position="top")
    
    p2 + p1 + plot_layout(widths = c(3,1))
  }
  plot_EC_integrated_development_celltypes = function (srat) { 
    ## take an integrated Seurat object, plot distributions over orig.ident
    library(Seurat)
    library(patchwork)
    library(ggplot2)
    library(reshape2)
    library(RColorBrewer)
    
    
    count_table <- table(srat@meta.data$development, srat@meta.data$celltype.full)
    count_mtx   <- as.data.frame.matrix(count_table)
    count_mtx$celltype <- rownames(count_mtx)
    melt_mtx    <- melt(count_mtx)
    melt_mtx$celltype <- as.factor(melt_mtx$celltype)
    
    celltype_size   <- aggregate(value ~ celltype, data = melt_mtx, FUN = sum)
    
    sorted_labels <- levels(srat$development)
    celltype_size$celltype <- factor(celltype_size$celltype,levels = sorted_labels)
    melt_mtx$celltype <- factor(melt_mtx$celltype,levels = sorted_labels)
    
    colnames(melt_mtx)[2] <- "dataset"
    
    
    p1 <- ggplot(celltype_size, aes(y= celltype,x = value)) + geom_bar(position="dodge", stat="identity",fill = "grey60") + 
      theme_bw() + scale_x_log10() + xlab("Cells per celltype, log10 scale") + ylab("")
    p2 <- ggplot(melt_mtx,aes(x=celltype,y=value,fill=dataset)) + 
      geom_bar(position="fill", stat="identity") + theme_bw() + coord_flip() + 
      scale_fill_manual(values= color_levelEC) +
      ylab("Fraction of cells in each dataset") + xlab("celltype number") + theme(legend.position="top")
    
    p2 + p1 + plot_layout(widths = c(3,1))
  }
  plot_EC_integrated_donor_celltypes = function (srat) { 
    ## take an integrated Seurat object, plot distributions over orig.ident
    library(Seurat)
    library(patchwork)
    library(ggplot2)
    library(reshape2)
    library(RColorBrewer)
    
    
    count_table <- table(srat@meta.data$donor, srat@meta.data$celltype.full)
    count_mtx   <- as.data.frame.matrix(count_table)
    count_mtx$celltype <- rownames(count_mtx)
    melt_mtx    <- melt(count_mtx)
    melt_mtx$celltype <- as.factor(melt_mtx$celltype)
    
    celltype_size   <- aggregate(value ~ celltype, data = melt_mtx, FUN = sum)
    
    sorted_labels <- levels(srat$donor)
    celltype_size$celltype <- factor(celltype_size$celltype,levels = sorted_labels)
    melt_mtx$celltype <- factor(melt_mtx$celltype,levels = sorted_labels)
    
    colnames(melt_mtx)[2] <- "dataset"
    
    
    p1 <- ggplot(celltype_size, aes(y= celltype,x = value)) + geom_bar(position="dodge", stat="identity",fill = "grey60") + 
      theme_bw() + scale_x_log10() + xlab("Cells per celltype, log10 scale") + ylab("")
    p2 <- ggplot(melt_mtx,aes(x=celltype,y=value,fill=dataset)) + 
      geom_bar(position="fill", stat="identity") + theme_bw() + coord_flip() + 
      scale_fill_manual(values= color_levelEC) +
      ylab("Fraction of cells in each dataset") + xlab("celltype number") + theme(legend.position="top")
    
    p2 + p1 + plot_layout(widths = c(3,1))
  }
  ##
  #FB
  plot_FB_integrated_group_clusters = function (srat) { 
    ## take an integrated Seurat object, plot distributions over orig.ident
    library(Seurat)
    library(patchwork)
    library(ggplot2)
    library(reshape2)
    library(RColorBrewer)
    
    
    count_table <- table(srat@meta.data$group, srat@meta.data$celltype.full)
    count_mtx   <- as.data.frame.matrix(count_table)
    count_mtx$celltype <- rownames(count_mtx)
    melt_mtx    <- melt(count_mtx)
    melt_mtx$celltype <- as.factor(melt_mtx$celltype)
    
    celltype_size   <- aggregate(value ~ celltype, data = melt_mtx, FUN = sum)
    
    sorted_labels <- levels(srat$group)
    celltype_size$celltype <- factor(celltype_size$celltype,levels = sorted_labels)
    melt_mtx$celltype <- factor(melt_mtx$celltype,levels = sorted_labels)
    
    colnames(melt_mtx)[2] <- "dataset"
    
    
    p1 <- ggplot(celltype_size, aes(y= celltype,x = value)) + geom_bar(position="dodge", stat="identity",fill = "grey60") + 
      theme_bw() + scale_x_log10() + xlab("Cells per celltype, log10 scale") + ylab("")
    p2 <- ggplot(melt_mtx,aes(x=celltype,y=value,fill=dataset)) + 
      geom_bar(position="fill", stat="identity") + theme_bw() + coord_flip() + 
      scale_fill_manual(values= color_levelFB) +
      ylab("Fraction of cells in each dataset") + xlab("celltype number") + theme(legend.position="top")
    
    p2 + p1 + plot_layout(widths = c(3,1))
  }
  plot_FB_integrated_Donor_ID_clusters = function (srat) { 
    ## take an integrated Seurat object, plot distributions over orig.ident
    library(Seurat)
    library(patchwork)
    library(ggplot2)
    library(reshape2)
    library(RColorBrewer)
    
    
    count_table <- table(srat@meta.data$Donor_ID, srat@meta.data$celltype.full)
    count_mtx   <- as.data.frame.matrix(count_table)
    count_mtx$celltype <- rownames(count_mtx)
    melt_mtx    <- melt(count_mtx)
    melt_mtx$celltype <- as.factor(melt_mtx$celltype)
    
    celltype_size   <- aggregate(value ~ celltype, data = melt_mtx, FUN = sum)
    
    sorted_labels <- levels(srat$Donor_ID)
    celltype_size$celltype <- factor(celltype_size$celltype,levels = sorted_labels)
    melt_mtx$celltype <- factor(melt_mtx$celltype,levels = sorted_labels)
    
    colnames(melt_mtx)[2] <- "dataset"
    
    
    p1 <- ggplot(celltype_size, aes(y= celltype,x = value)) + geom_bar(position="dodge", stat="identity",fill = "grey60") + 
      theme_bw() + scale_x_log10() + xlab("Cells per celltype, log10 scale") + ylab("")
    p2 <- ggplot(melt_mtx,aes(x=celltype,y=value,fill=dataset)) + 
      geom_bar(position="fill", stat="identity") + theme_bw() + coord_flip() + 
      scale_fill_manual(values= color_levelFB)+
      ylab("Fraction of cells in each dataset") + xlab("celltype number") + theme(legend.position="top")
    
    p2 + p1 + plot_layout(widths = c(3,1))
  }
  plot_FB_integrated_condition_celltypes = function (srat) { 
    ## take an integrated Seurat object, plot distributions over orig.ident
    library(Seurat)
    library(patchwork)
    library(ggplot2)
    library(reshape2)
    library(RColorBrewer)
    
    
    count_table <- table(srat@meta.data$condition, srat@meta.data$celltype.full)
    count_mtx   <- as.data.frame.matrix(count_table)
    count_mtx$celltype <- rownames(count_mtx)
    melt_mtx    <- melt(count_mtx)
    melt_mtx$celltype <- as.factor(melt_mtx$celltype)
    
    celltype_size   <- aggregate(value ~ celltype, data = melt_mtx, FUN = sum)
    
    sorted_labels <- levels(srat$condition)
    celltype_size$celltype <- factor(celltype_size$celltype,levels = sorted_labels)
    melt_mtx$celltype <- factor(melt_mtx$celltype,levels = sorted_labels)
    
    colnames(melt_mtx)[2] <- "dataset"
    
    
    p1 <- ggplot(celltype_size, aes(y= celltype,x = value)) + geom_bar(position="dodge", stat="identity",fill = "grey60") + 
      theme_bw() + scale_x_log10() + xlab("Cells per celltype, log10 scale") + ylab("")
    p2 <- ggplot(melt_mtx,aes(x=celltype,y=value,fill=dataset)) + 
      geom_bar(position="fill", stat="identity") + theme_bw() + coord_flip() + 
      scale_fill_manual(values= color_levelFB) +
      ylab("Fraction of cells in each dataset") + xlab("celltype number") + theme(legend.position="top")
    
    p2 + p1 + plot_layout(widths = c(3,1))
  }
  plot_FB_integrated_classification_celltypes = function (srat) { 
    ## take an integrated Seurat object, plot distributions over orig.ident
    library(Seurat)
    library(patchwork)
    library(ggplot2)
    library(reshape2)
    library(RColorBrewer)
    
    
    count_table <- table(srat@meta.data$classification, srat@meta.data$celltype.full)
    count_mtx   <- as.data.frame.matrix(count_table)
    count_mtx$celltype <- rownames(count_mtx)
    melt_mtx    <- melt(count_mtx)
    melt_mtx$celltype <- as.factor(melt_mtx$celltype)
    
    celltype_size   <- aggregate(value ~ celltype, data = melt_mtx, FUN = sum)
    
    sorted_labels <- levels(srat$classification)
    celltype_size$celltype <- factor(celltype_size$celltype,levels = sorted_labels)
    melt_mtx$celltype <- factor(melt_mtx$celltype,levels = sorted_labels)
    
    colnames(melt_mtx)[2] <- "dataset"
    
    
    p1 <- ggplot(celltype_size, aes(y= celltype,x = value)) + geom_bar(position="dodge", stat="identity",fill = "grey60") + 
      theme_bw() + scale_x_log10() + xlab("Cells per celltype, log10 scale") + ylab("")
    p2 <- ggplot(melt_mtx,aes(x=celltype,y=value,fill=dataset)) + 
      geom_bar(position="fill", stat="identity") + theme_bw() + coord_flip() + 
      scale_fill_manual(values= color_levelFB) +
      ylab("Fraction of cells in each dataset") + xlab("celltype number") + theme(legend.position="top")
    
    p2 + p1 + plot_layout(widths = c(3,1))
  }
  plot_FB_integrated_tissue_celltypes = function (srat) { 
    ## take an integrated Seurat object, plot distributions over orig.ident
    library(Seurat)
    library(patchwork)
    library(ggplot2)
    library(reshape2)
    library(RColorBrewer)
    
    
    count_table <- table(srat@meta.data$tissue, srat@meta.data$celltype.full)
    count_mtx   <- as.data.frame.matrix(count_table)
    count_mtx$celltype <- rownames(count_mtx)
    melt_mtx    <- melt(count_mtx)
    melt_mtx$celltype <- as.factor(melt_mtx$celltype)
    
    celltype_size   <- aggregate(value ~ celltype, data = melt_mtx, FUN = sum)
    
    sorted_labels <- levels(srat$tissue)
    celltype_size$celltype <- factor(celltype_size$celltype,levels = sorted_labels)
    melt_mtx$celltype <- factor(melt_mtx$celltype,levels = sorted_labels)
    
    colnames(melt_mtx)[2] <- "dataset"
    
    
    p1 <- ggplot(celltype_size, aes(y= celltype,x = value)) + geom_bar(position="dodge", stat="identity",fill = "grey60") + 
      theme_bw() + scale_x_log10() + xlab("Cells per celltype, log10 scale") + ylab("")
    p2 <- ggplot(melt_mtx,aes(x=celltype,y=value,fill=dataset)) + 
      geom_bar(position="fill", stat="identity") + theme_bw() + coord_flip() + 
      scale_fill_manual(values= color_levelFB) +
      ylab("Fraction of cells in each dataset") + xlab("celltype number") + theme(legend.position="top")
    
    p2 + p1 + plot_layout(widths = c(3,1))
  }
  plot_FB_integrated_development_celltypes = function (srat) { 
    ## take an integrated Seurat object, plot distributions over orig.ident
    library(Seurat)
    library(patchwork)
    library(ggplot2)
    library(reshape2)
    library(RColorBrewer)
    
    
    count_table <- table(srat@meta.data$development, srat@meta.data$celltype.full)
    count_mtx   <- as.data.frame.matrix(count_table)
    count_mtx$celltype <- rownames(count_mtx)
    melt_mtx    <- melt(count_mtx)
    melt_mtx$celltype <- as.factor(melt_mtx$celltype)
    
    celltype_size   <- aggregate(value ~ celltype, data = melt_mtx, FUN = sum)
    
    sorted_labels <- levels(srat$development)
    celltype_size$celltype <- factor(celltype_size$celltype,levels = sorted_labels)
    melt_mtx$celltype <- factor(melt_mtx$celltype,levels = sorted_labels)
    
    colnames(melt_mtx)[2] <- "dataset"
    
    
    p1 <- ggplot(celltype_size, aes(y= celltype,x = value)) + geom_bar(position="dodge", stat="identity",fill = "grey60") + 
      theme_bw() + scale_x_log10() + xlab("Cells per celltype, log10 scale") + ylab("")
    p2 <- ggplot(melt_mtx,aes(x=celltype,y=value,fill=dataset)) + 
      geom_bar(position="fill", stat="identity") + theme_bw() + coord_flip() + 
      scale_fill_manual(values= color_levelFB) +
      ylab("Fraction of cells in each dataset") + xlab("celltype number") + theme(legend.position="top")
    
    p2 + p1 + plot_layout(widths = c(3,1))
  }
  plot_FB_integrated_donor_celltypes = function (srat) { 
    ## take an integrated Seurat object, plot distributions over orig.ident
    library(Seurat)
    library(patchwork)
    library(ggplot2)
    library(reshape2)
    library(RColorBrewer)
    
    
    count_table <- table(srat@meta.data$donor, srat@meta.data$celltype.full)
    count_mtx   <- as.data.frame.matrix(count_table)
    count_mtx$celltype <- rownames(count_mtx)
    melt_mtx    <- melt(count_mtx)
    melt_mtx$celltype <- as.factor(melt_mtx$celltype)
    
    celltype_size   <- aggregate(value ~ celltype, data = melt_mtx, FUN = sum)
    
    sorted_labels <- levels(srat$donor)
    celltype_size$celltype <- factor(celltype_size$celltype,levels = sorted_labels)
    melt_mtx$celltype <- factor(melt_mtx$celltype,levels = sorted_labels)
    
    colnames(melt_mtx)[2] <- "dataset"
    
    
    p1 <- ggplot(celltype_size, aes(y= celltype,x = value)) + geom_bar(position="dodge", stat="identity",fill = "grey60") + 
      theme_bw() + scale_x_log10() + xlab("Cells per celltype, log10 scale") + ylab("")
    p2 <- ggplot(melt_mtx,aes(x=celltype,y=value,fill=dataset)) + 
      geom_bar(position="fill", stat="identity") + theme_bw() + coord_flip() + 
      scale_fill_manual(values= color_levelFB) +
      ylab("Fraction of cells in each dataset") + xlab("celltype number") + theme(legend.position="top")
    
    p2 + p1 + plot_layout(widths = c(3,1))
  }
  ##
  #VSMC/Peri
  plot_VSMC_Peri_integrated_group_clusters = function (srat) { 
    ## take an integrated Seurat object, plot distributions over orig.ident
    library(Seurat)
    library(patchwork)
    library(ggplot2)
    library(reshape2)
    library(RColorBrewer)
    
    
    count_table <- table(srat@meta.data$group, srat@meta.data$celltype.full)
    count_mtx   <- as.data.frame.matrix(count_table)
    count_mtx$celltype <- rownames(count_mtx)
    melt_mtx    <- melt(count_mtx)
    melt_mtx$celltype <- as.factor(melt_mtx$celltype)
    
    celltype_size   <- aggregate(value ~ celltype, data = melt_mtx, FUN = sum)
    
    sorted_labels <- levels(srat$group)
    celltype_size$celltype <- factor(celltype_size$celltype,levels = sorted_labels)
    melt_mtx$celltype <- factor(melt_mtx$celltype,levels = sorted_labels)
    
    colnames(melt_mtx)[2] <- "dataset"
    
    
    p1 <- ggplot(celltype_size, aes(y= celltype,x = value)) + geom_bar(position="dodge", stat="identity",fill = "grey60") + 
      theme_bw() + scale_x_log10() + xlab("Cells per celltype, log10 scale") + ylab("")
    p2 <- ggplot(melt_mtx,aes(x=celltype,y=value,fill=dataset)) + 
      geom_bar(position="fill", stat="identity") + theme_bw() + coord_flip() + 
      scale_fill_manual(values= color_levelVSMC_Peri) +
      ylab("Fraction of cells in each dataset") + xlab("celltype number") + theme(legend.position="top")
    
    p2 + p1 + plot_layout(widths = c(3,1))
  }
  plot_VSMC_Peri_integrated_Donor_ID_clusters = function (srat) { 
    ## take an integrated Seurat object, plot distributions over orig.ident
    library(Seurat)
    library(patchwork)
    library(ggplot2)
    library(reshape2)
    library(RColorBrewer)
    
    
    count_table <- table(srat@meta.data$Donor_ID, srat@meta.data$celltype.full)
    count_mtx   <- as.data.frame.matrix(count_table)
    count_mtx$celltype <- rownames(count_mtx)
    melt_mtx    <- melt(count_mtx)
    melt_mtx$celltype <- as.factor(melt_mtx$celltype)
    
    celltype_size   <- aggregate(value ~ celltype, data = melt_mtx, FUN = sum)
    
    sorted_labels <- levels(srat$Donor_ID)
    celltype_size$celltype <- factor(celltype_size$celltype,levels = sorted_labels)
    melt_mtx$celltype <- factor(melt_mtx$celltype,levels = sorted_labels)
    
    colnames(melt_mtx)[2] <- "dataset"
    
    
    p1 <- ggplot(celltype_size, aes(y= celltype,x = value)) + geom_bar(position="dodge", stat="identity",fill = "grey60") + 
      theme_bw() + scale_x_log10() + xlab("Cells per celltype, log10 scale") + ylab("")
    p2 <- ggplot(melt_mtx,aes(x=celltype,y=value,fill=dataset)) + 
      geom_bar(position="fill", stat="identity") + theme_bw() + coord_flip() + 
      scale_fill_manual(values= color_levelVSMC_Peri)+
      ylab("Fraction of cells in each dataset") + xlab("celltype number") + theme(legend.position="top")
    
    p2 + p1 + plot_layout(widths = c(3,1))
  }
  plot_VSMC_Peri_integrated_condition_celltypes = function (srat) { 
    ## take an integrated Seurat object, plot distributions over orig.ident
    library(Seurat)
    library(patchwork)
    library(ggplot2)
    library(reshape2)
    library(RColorBrewer)
    
    
    count_table <- table(srat@meta.data$condition, srat@meta.data$celltype.full)
    count_mtx   <- as.data.frame.matrix(count_table)
    count_mtx$celltype <- rownames(count_mtx)
    melt_mtx    <- melt(count_mtx)
    melt_mtx$celltype <- as.factor(melt_mtx$celltype)
    
    celltype_size   <- aggregate(value ~ celltype, data = melt_mtx, FUN = sum)
    
    sorted_labels <- levels(srat$condition)
    celltype_size$celltype <- factor(celltype_size$celltype,levels = sorted_labels)
    melt_mtx$celltype <- factor(melt_mtx$celltype,levels = sorted_labels)
    
    colnames(melt_mtx)[2] <- "dataset"
    
    
    p1 <- ggplot(celltype_size, aes(y= celltype,x = value)) + geom_bar(position="dodge", stat="identity",fill = "grey60") + 
      theme_bw() + scale_x_log10() + xlab("Cells per celltype, log10 scale") + ylab("")
    p2 <- ggplot(melt_mtx,aes(x=celltype,y=value,fill=dataset)) + 
      geom_bar(position="fill", stat="identity") + theme_bw() + coord_flip() + 
      scale_fill_manual(values= color_levelVSMC_Peri) +
      ylab("Fraction of cells in each dataset") + xlab("celltype number") + theme(legend.position="top")
    
    p2 + p1 + plot_layout(widths = c(3,1))
  }
  plot_VSMC_Peri_integrated_classification_celltypes = function (srat) { 
    ## take an integrated Seurat object, plot distributions over orig.ident
    library(Seurat)
    library(patchwork)
    library(ggplot2)
    library(reshape2)
    library(RColorBrewer)
    
    
    count_table <- table(srat@meta.data$classification, srat@meta.data$celltype.full)
    count_mtx   <- as.data.frame.matrix(count_table)
    count_mtx$celltype <- rownames(count_mtx)
    melt_mtx    <- melt(count_mtx)
    melt_mtx$celltype <- as.factor(melt_mtx$celltype)
    
    celltype_size   <- aggregate(value ~ celltype, data = melt_mtx, FUN = sum)
    
    sorted_labels <- levels(srat$classification)
    celltype_size$celltype <- factor(celltype_size$celltype,levels = sorted_labels)
    melt_mtx$celltype <- factor(melt_mtx$celltype,levels = sorted_labels)
    
    colnames(melt_mtx)[2] <- "dataset"
    
    
    p1 <- ggplot(celltype_size, aes(y= celltype,x = value)) + geom_bar(position="dodge", stat="identity",fill = "grey60") + 
      theme_bw() + scale_x_log10() + xlab("Cells per celltype, log10 scale") + ylab("")
    p2 <- ggplot(melt_mtx,aes(x=celltype,y=value,fill=dataset)) + 
      geom_bar(position="fill", stat="identity") + theme_bw() + coord_flip() + 
      scale_fill_manual(values= color_levelVSMC_Peri) +
      ylab("Fraction of cells in each dataset") + xlab("celltype number") + theme(legend.position="top")
    
    p2 + p1 + plot_layout(widths = c(3,1))
  }
  plot_VSMC_Peri_integrated_tissue_celltypes = function (srat) { 
    ## take an integrated Seurat object, plot distributions over orig.ident
    library(Seurat)
    library(patchwork)
    library(ggplot2)
    library(reshape2)
    library(RColorBrewer)
    
    
    count_table <- table(srat@meta.data$tissue, srat@meta.data$celltype.full)
    count_mtx   <- as.data.frame.matrix(count_table)
    count_mtx$celltype <- rownames(count_mtx)
    melt_mtx    <- melt(count_mtx)
    melt_mtx$celltype <- as.factor(melt_mtx$celltype)
    
    celltype_size   <- aggregate(value ~ celltype, data = melt_mtx, FUN = sum)
    
    sorted_labels <- levels(srat$tissue)
    celltype_size$celltype <- factor(celltype_size$celltype,levels = sorted_labels)
    melt_mtx$celltype <- factor(melt_mtx$celltype,levels = sorted_labels)
    
    colnames(melt_mtx)[2] <- "dataset"
    
    
    p1 <- ggplot(celltype_size, aes(y= celltype,x = value)) + geom_bar(position="dodge", stat="identity",fill = "grey60") + 
      theme_bw() + scale_x_log10() + xlab("Cells per celltype, log10 scale") + ylab("")
    p2 <- ggplot(melt_mtx,aes(x=celltype,y=value,fill=dataset)) + 
      geom_bar(position="fill", stat="identity") + theme_bw() + coord_flip() + 
      scale_fill_manual(values= color_levelVSMC_Peri) +
      ylab("Fraction of cells in each dataset") + xlab("celltype number") + theme(legend.position="top")
    
    p2 + p1 + plot_layout(widths = c(3,1))
  }
  plot_VSMC_Peri_integrated_development_celltypes = function (srat) { 
    ## take an integrated Seurat object, plot distributions over orig.ident
    library(Seurat)
    library(patchwork)
    library(ggplot2)
    library(reshape2)
    library(RColorBrewer)
    
    
    count_table <- table(srat@meta.data$development, srat@meta.data$celltype.full)
    count_mtx   <- as.data.frame.matrix(count_table)
    count_mtx$celltype <- rownames(count_mtx)
    melt_mtx    <- melt(count_mtx)
    melt_mtx$celltype <- as.factor(melt_mtx$celltype)
    
    celltype_size   <- aggregate(value ~ celltype, data = melt_mtx, FUN = sum)
    
    sorted_labels <- levels(srat$development)
    celltype_size$celltype <- factor(celltype_size$celltype,levels = sorted_labels)
    melt_mtx$celltype <- factor(melt_mtx$celltype,levels = sorted_labels)
    
    colnames(melt_mtx)[2] <- "dataset"
    
    
    p1 <- ggplot(celltype_size, aes(y= celltype,x = value)) + geom_bar(position="dodge", stat="identity",fill = "grey60") + 
      theme_bw() + scale_x_log10() + xlab("Cells per celltype, log10 scale") + ylab("")
    p2 <- ggplot(melt_mtx,aes(x=celltype,y=value,fill=dataset)) + 
      geom_bar(position="fill", stat="identity") + theme_bw() + coord_flip() + 
      scale_fill_manual(values= color_levelVSMC_Peri) +
      ylab("Fraction of cells in each dataset") + xlab("celltype number") + theme(legend.position="top")
    
    p2 + p1 + plot_layout(widths = c(3,1))
  }
  plot_VSMC_Peri_integrated_donor_celltypes = function (srat) { 
    ## take an integrated Seurat object, plot distributions over orig.ident
    library(Seurat)
    library(patchwork)
    library(ggplot2)
    library(reshape2)
    library(RColorBrewer)
    
    
    count_table <- table(srat@meta.data$donor, srat@meta.data$celltype.full)
    count_mtx   <- as.data.frame.matrix(count_table)
    count_mtx$celltype <- rownames(count_mtx)
    melt_mtx    <- melt(count_mtx)
    melt_mtx$celltype <- as.factor(melt_mtx$celltype)
    
    celltype_size   <- aggregate(value ~ celltype, data = melt_mtx, FUN = sum)
    
    sorted_labels <- levels(srat$donor)
    celltype_size$celltype <- factor(celltype_size$celltype,levels = sorted_labels)
    melt_mtx$celltype <- factor(melt_mtx$celltype,levels = sorted_labels)
    
    colnames(melt_mtx)[2] <- "dataset"
    
    
    p1 <- ggplot(celltype_size, aes(y= celltype,x = value)) + geom_bar(position="dodge", stat="identity",fill = "grey60") + 
      theme_bw() + scale_x_log10() + xlab("Cells per celltype, log10 scale") + ylab("")
    p2 <- ggplot(melt_mtx,aes(x=celltype,y=value,fill=dataset)) + 
      geom_bar(position="fill", stat="identity") + theme_bw() + coord_flip() + 
      scale_fill_manual(values= color_levelVSMC_Peri) +
      ylab("Fraction of cells in each dataset") + xlab("celltype number") + theme(legend.position="top")
    
    p2 + p1 + plot_layout(widths = c(3,1))
  }
  #
  
  #####GO-Terms Enrichr
  library(enrichR)
  websiteLive <- getOption("enrichR.live")
  if (websiteLive) {
    listEnrichrSites()
    setEnrichrSite("Enrichr") # Human genes   
  }
  if (websiteLive) dbs <- listEnrichrDbs()
  if (websiteLive) head(dbs)
  dbs <- c("GO_Molecular_Function_2023", "GO_Cellular_Component_2023", "GO_Biological_Process_2023")
}
#(2)_data preproceeding####
##load datasets##

#Xin#
load(file="D:/Direder/Projekte/D003_10xThymom/Revision/Data_Matthias/Xin_ready_to_use_2.RData")
table(xin$name)
Idents(xin)<- xin$name
xin1 <- subset(xin,idents = "T01_1")
xin1 <- xin1@assays$RNA@counts
xin2 <- subset(xin,idents = "T02_1")
xin2 <- xin2@assays$RNA@counts
xin3 <- subset(xin,idents = "T03_1")
xin3 <- xin3@assays$RNA@counts
xin4 <- subset(xin,idents = "T04_1")
xin4 <- xin4@assays$RNA@counts
xin5 <- subset(xin,idents = "T05_1")
xin5 <- xin5@assays$RNA@counts
xin6 <- subset(xin,idents = "T06_1")
xin6 <- xin6@assays$RNA@counts

#yasumizu#
load(file="D:/Direder/Projekte/D003_10xThymom/Revision/Data_Matthias/yasmizu_ready_to_use.RData")
table(yasumizu_seurat$donor_id)
Idents(yasumizu_seurat)<- yasumizu_seurat$donor_id
yasumizu_MG03 <- subset(yasumizu_seurat,idents = "MG03")
yasumizu_MG03 <- yasumizu_MG03@assays$RNA@counts
yasumizu_MG21 <- subset(yasumizu_seurat,idents = "MG21")
yasumizu_MG21 <- yasumizu_MG21@assays$RNA@counts
yasumizu_MG22 <- subset(yasumizu_seurat,idents = "MG22")
yasumizu_MG22 <- yasumizu_MG22@assays$RNA@counts
yasumizu_MG23 <- subset(yasumizu_seurat,idents = "MG23")
yasumizu_MG23 <- yasumizu_MG23@assays$RNA@counts

#park#
load(file="D:/Direder/Projekte/D003_10xThymom/Revision/Data_Matthias/seurat_PARK_ready_to_use_FIXED.RData")
table(sce_park$donor)
Idents(sce_park)<- sce_park$donor
park_A16 <- subset(sce_park,idents = "A16")
park_A16 <- park_A16@assays$RNA@counts
park_A43 <- subset(sce_park,idents = "A43")
park_A43 <- park_A43@assays$RNA@counts
park_A45 <- subset(sce_park,idents = "A45")
park_A45 <- park_A45@assays$RNA@counts
park_C34 <- subset(sce_park,idents = "C34")
park_C34 <- park_C34@assays$RNA@counts
park_C40 <- subset(sce_park,idents = "C40")
park_C40 <- park_C40@assays$RNA@counts
park_C41 <- subset(sce_park,idents = "C41")
park_C41 <- park_C41@assays$RNA@counts
park_F21 <- subset(sce_park,idents = "F21")
park_F21 <- park_F21@assays$RNA@counts
park_F22 <- subset(sce_park,idents = "F22")
park_F22 <- park_F22@assays$RNA@counts
park_F23 <- subset(sce_park,idents = "F23")
park_F23 <- park_F23@assays$RNA@counts
park_F29 <- subset(sce_park,idents = "F29")
park_F29 <- park_F29@assays$RNA@counts
park_F30 <- subset(sce_park,idents = "F30")
park_F30 <- park_F30@assays$RNA@counts
park_F38 <- subset(sce_park,idents = "F38")
park_F38 <- park_F38@assays$RNA@counts
park_F41 <- subset(sce_park,idents = "F41")
park_F41 <- park_F41@assays$RNA@counts
park_F45 <- subset(sce_park,idents = "F45")
park_F45 <- park_F45@assays$RNA@counts
park_F64 <- subset(sce_park,idents = "F64")
park_F64 <- park_F64@assays$RNA@counts
park_F67 <- subset(sce_park,idents = "F67")
park_F67 <- park_F67@assays$RNA@counts
park_F74 <- subset(sce_park,idents = "F74")
park_F74 <- park_F74@assays$RNA@counts
park_F83 <- subset(sce_park,idents = "F83")
park_F83 <- park_F83@assays$RNA@counts
park_P1 <- subset(sce_park,idents = "P1")
park_P1 <- park_P1@assays$RNA@counts
park_P2 <- subset(sce_park,idents = "P2")
park_P2 <- park_P2@assays$RNA@counts
park_P3<- subset(sce_park,idents = "P3")
park_P3 <- park_P3@assays$RNA@counts
park_T03 <- subset(sce_park,idents = "T03")
park_T03 <- park_T03@assays$RNA@counts
park_T06 <- subset(sce_park,idents = "T06")
park_T06 <- park_T06@assays$RNA@counts
park_T07 <- subset(sce_park,idents = "T07")
park_T07 <- park_T07@assays$RNA@counts

#bautista#
bautista_a25<-Read10X("D:/Direder/Projekte/D003_10xThymom/Data/PMID_33597545_Bautista/Adult_25y")
bautista_f19<-Read10X("D:/Direder/Projekte/D003_10xThymom/Data/PMID_33597545_Bautista/Fetal_19w")
bautista_f231<-Read10X("D:/Direder/Projekte/D003_10xThymom/Data/PMID_33597545_Bautista/Fetal_23w_1")
bautista_f232<-Read10X("D:/Direder/Projekte/D003_10xThymom/Data/PMID_33597545_Bautista/Fetal_23w_2")
bautista_p6<-Read10X("D:/Direder/Projekte/D003_10xThymom/Data/PMID_33597545_Bautista/Postnatal_6d")
bautista_p101<-Read10X("D:/Direder/Projekte/D003_10xThymom/Data/PMID_33597545_Bautista/Postnatal_10m_1")
bautista_p102<-Read10X("D:/Direder/Projekte/D003_10xThymom/Data/PMID_33597545_Bautista/Postnatal_10m_2")

#direder#
direder_TH1<-Read10X("D:/Direder/Projekte/D003_10xThymom/Data/Thymus_hyperplasie/TH1")
direder_TH2<-Read10X("D:/Direder/Projekte/D003_10xThymom/Data/Thymus_hyperplasie/TH2")
direder_TH3<-Read10X("D:/Direder/Projekte/D003_10xThymom/Data/Thymus_hyperplasie/TH3")
direder_TH4<-Read10X("D:/Direder/Projekte/D003_10xThymom/Data/Thymus_hyperplasie/TH4")
direder_TB1<-Read10X("D:/Direder/Projekte/D003_10xThymom/Data/Thymom/TB1_1")
direder_TB2<-Read10X("D:/Direder/Projekte/D003_10xThymom/Data/Thymom/TB1_2")
direder_TB3<-Read10X("D:/Direder/Projekte/D003_10xThymom/Data/Thymom/TB2")
direder_TB4<-Read10X("D:/Direder/Projekte/D003_10xThymom/Data/Thymom/TB2B3")
direder_TC1<-Read10X("D:/Direder/Projekte/D003_10xThymom/Data/Thymus_carcinom/TC_1")
direder_TC2<-Read10X("D:/Direder/Projekte/D003_10xThymom/Data/Thymus_carcinom/TC_2")
direder_TH_MG1<-Read10X("D:/Direder/Projekte/D003_10xThymom/Data/Thymus_hyperplasie/TH_MG_1")
direder_TH_MG2<-Read10X("D:/Direder/Projekte/D003_10xThymom/Data/Thymus_hyperplasie/TH_MG_2")

#StartLoop_1
rawdat <- list(xin1, xin2, xin3, xin4, xin5, xin6,
               yasumizu_MG03, yasumizu_MG21, yasumizu_MG22, yasumizu_MG23,
               park_A16, park_A43, park_A45, park_C34, park_C40, park_C41, park_F21, park_F22, park_F23, park_F29, park_F30, park_F38,
               park_F41, park_F45, park_F64, park_F67, park_F74, park_F83, park_P1, park_P2, park_P3, park_T03, park_T06, park_T07,
               bautista_a25, bautista_f19, bautista_f231, bautista_f232, bautista_p6, bautista_p101, bautista_p102,
               direder_TH1, direder_TH2, direder_TH3, direder_TH4, direder_TB1, direder_TB2, direder_TB3, direder_TB4, direder_TC1, direder_TC2, direder_TH_MG1, direder_TH_MG2)

#convert to NCBI Symbol & remove duplicates
rawdat_1 <- lapply(rawdat, function(x){
  Signlist1=alias2SymbolUsingNCBI(as.character(rownames(x)),paste0(c("D:/Direder/SC-EnsembID/Homo_sapiens.gene_info")),required.columns = c("GeneID","Symbol"))
  Signlist2=as.data.frame(cbind(Signlist1,as.character(rownames(x))))
  colnames(Signlist2)[3]=c("rname")
  Signlist2$new = ifelse(is.na(Signlist2$Symbol), as.character(Signlist2$rname), Signlist2$Symbol)
  rownames(x)= as.character(Signlist2$new)
  x <- x[!duplicated(rownames(x)),]
})

#save(rawdat_1,file="rawdat_1.RData")
#remove surplus objects
rm(rawdat, 
   xin1, xin2, xin3, xin4, xin5, xin6,
   yasumizu_MG03, yasumizu_MG21, yasumizu_MG22, yasumizu_MG23,
   park_A16, park_A43, park_A45, park_C34, park_C40, park_C41, park_F21, park_F22, park_F23, park_F29, park_F30, park_F38,
   park_F41, park_F45, park_F64, park_F67, park_F74, park_F83, park_P1, park_P2, park_P3, park_T03, park_T06, park_T07,
   bautista_a25, bautista_f19, bautista_f231, bautista_f232, bautista_p6, bautista_p101, bautista_p102,
   direder_TH1, direder_TH2, direder_TH3, direder_TH4, direder_TB1, direder_TB2, direder_TB3, direder_TB4, direder_TC1, direder_TC2, direder_TH_MG1, direder_TH_MG2, sce_park, xin, yasumizu_seurat)


# keep only Genes detected in all dataset
ROWLIST1 <- rawdat_1[[1]]
ROWLIST2 <- rawdat_1[[7]]
ROWLIST3 <- rawdat_1[[11]]
ROWLIST4 <- rawdat_1[[35]]
ROWLIST5 <- rawdat_1[[42]]

dim(rawdat_1[[1]])
dim(rawdat_1[[7]])
dim(rawdat_1[[11]])
dim(rawdat_1[[35]])
dim(rawdat_1[[42]])


Temp1<-merge(ROWLIST1,ROWLIST2, by.x = rownames(ROWLIST1), by.y=rownames(ROWLIST2),all.x=F, all.y=F)
Temp1 <- Temp1[,1:2]
dim(ROWLIST1)
dim(ROWLIST2)
dim(Temp1)

Temp2<-merge(Temp1,ROWLIST3, by.x = rownames(Temp1), by.y=rownames(ROWLIST3),all.x=F, all.y=F)
Temp2 <- Temp2[,1:2]
dim(Temp1)
dim(ROWLIST3)
dim(Temp2)

Temp3<-merge(Temp2,ROWLIST4, by.x = rownames(Temp2), by.y=rownames(ROWLIST4),all.x=F, all.y=F)
Temp3 <- Temp3[,1:2]
dim(Temp2)
dim(ROWLIST4)
dim(Temp3)

Temp4<-merge(Temp3,ROWLIST5, by.x = rownames(Temp3), by.y=rownames(ROWLIST5),all.x=F, all.y=F)
Temp4 <- Temp3[,1:2]
dim(Temp3)
dim(ROWLIST5)
dim(Temp4)

rawdat_2 <- lapply(rawdat_1, function(x){
  Temp5<-join.Matrix(Temp4,x, by.x = rownames(Temp4), by.y=rownames(x), all.x = F, all.y = F)
  x <- Temp5[,3:ncol(Temp5)]
})

dim(rawdat_2[[1]])
dim(rawdat_2[[7]])
dim(rawdat_2[[11]])
dim(rawdat_2[[35]])
dim(rawdat_2[[42]])

#removal of background noise
# Create a SingleCellExperiment object and run decontX
rawdat_2a <- lapply(rawdat_2, function (x) { x<- SingleCellExperiment(list(counts = x))})
rawdat_2b <- lapply(rawdat_2a, function (x) { x <- decontX(x)})

#Create Sparcematrix from a SCE with decontX results
rawdat_3 <- lapply(rawdat_2b, function(x){
  x<-CreateSeuratObject(round(decontXcounts(x)))
})


#(3)_separate list for flags####
xin1_seurat <- rawdat_3[[1]] 
xin2_seurat <- rawdat_3[[2]] 
xin3_seurat <- rawdat_3[[3]] 
xin4_seurat <- rawdat_3[[4]] 
xin5_seurat <- rawdat_3[[5]] 
xin6_seurat <- rawdat_3[[6]] 
yasumizu_MG03_seurat <- rawdat_3[[7]] 
yasumizu_MG21_seurat <- rawdat_3[[8]] 
yasumizu_MG22_seurat <- rawdat_3[[9]] 
yasumizu_MG23_seurat <- rawdat_3[[10]] 
park_A16_seurat <- rawdat_3[[11]] 
park_A43_seurat <- rawdat_3[[12]] 
park_A45_seurat <- rawdat_3[[13]] 
park_C34_seurat <- rawdat_3[[14]] 
park_C40_seurat <- rawdat_3[[15]] 
park_C41_seurat <- rawdat_3[[16]] 
park_F21_seurat <- rawdat_3[[17]] 
park_F22_seurat <- rawdat_3[[18]] 
park_F23_seurat <- rawdat_3[[19]] 
park_F29_seurat <- rawdat_3[[20]] 
park_F30_seurat <- rawdat_3[[21]] 
park_F38_seurat <- rawdat_3[[22]] 
park_F41_seurat <- rawdat_3[[23]] 
park_F45_seurat <- rawdat_3[[24]] 
park_F64_seurat <- rawdat_3[[25]] 
park_F67_seurat <- rawdat_3[[26]] 
park_F74_seurat <- rawdat_3[[27]] 
park_F83_seurat <- rawdat_3[[28]] 
park_P1_seurat <- rawdat_3[[29]] 
park_P2_seurat <- rawdat_3[[30]] 
park_P3_seurat <- rawdat_3[[31]] 
park_T03_seurat <- rawdat_3[[32]] 
park_T06_seurat <- rawdat_3[[33]] 
park_T07_seurat <- rawdat_3[[34]] 
bautista_a25_seurat <- rawdat_3[[35]] 
bautista_f19_seurat <- rawdat_3[[36]] 
bautista_f231_seurat <- rawdat_3[[37]] 
bautista_f232_seurat <- rawdat_3[[38]] 
bautista_p6_seurat <- rawdat_3[[39]] 
bautista_p101_seurat <- rawdat_3[[40]] 
bautista_p102_seurat <- rawdat_3[[41]] 
direder_TH1_seurat <- rawdat_3[[42]] 
direder_TH2_seurat <- rawdat_3[[43]] 
direder_TH3_seurat <- rawdat_3[[44]] 
direder_TH4_seurat <- rawdat_3[[45]] 
direder_TB1_seurat <- rawdat_3[[46]] 
direder_TB2_seurat <- rawdat_3[[47]] 
direder_TB3_seurat <- rawdat_3[[48]] 
direder_TB4_seurat <- rawdat_3[[49]] 
direder_TC1_seurat <- rawdat_3[[50]] 
direder_TC2_seurat <- rawdat_3[[51]] 
direder_TH_MG1_seurat <- rawdat_3[[52]] 
direder_TH_MG2_seurat <- rawdat_3[[53]] 

#remove surplus objects
rm(ROWLIST1,ROWLIST2, ROWLIST3, ROWLIST4,ROWLIST5,Temp1, Temp2, Temp3, Temp4, rawdat_1,rawdat_2, rawdat_2a,rawdat_2b,rawdat_3)

#(4)_flag-1_group####
xin1_seurat <- AddMetaData(xin1_seurat, "Xin", col.name = "group")
xin2_seurat <- AddMetaData(xin2_seurat, "Xin", col.name = "group")
xin3_seurat <- AddMetaData(xin3_seurat, "Xin", col.name = "group")
xin4_seurat <- AddMetaData(xin4_seurat, "Xin", col.name = "group")
xin5_seurat <- AddMetaData(xin5_seurat, "Xin", col.name = "group")
xin6_seurat <- AddMetaData(xin6_seurat, "Xin", col.name = "group")
yasumizu_MG03_seurat <- AddMetaData(yasumizu_MG03_seurat, "Yasumizu", col.name = "group")
yasumizu_MG21_seurat <- AddMetaData(yasumizu_MG21_seurat, "Yasumizu", col.name = "group")
yasumizu_MG22_seurat <- AddMetaData(yasumizu_MG22_seurat, "Yasumizu", col.name = "group")
yasumizu_MG23_seurat <- AddMetaData(yasumizu_MG23_seurat, "Yasumizu", col.name = "group")
park_A16_seurat <- AddMetaData(park_A16_seurat, "Park", col.name = "group")
park_A43_seurat <- AddMetaData(park_A43_seurat, "Park", col.name = "group") 
park_A45_seurat <- AddMetaData(park_A45_seurat, "Park", col.name = "group")
park_C34_seurat <- AddMetaData(park_C34_seurat, "Park", col.name = "group")
park_C40_seurat <- AddMetaData(park_C40_seurat, "Park", col.name = "group") 
park_C41_seurat <- AddMetaData(park_C41_seurat, "Park", col.name = "group")
park_F21_seurat <- AddMetaData(park_F21_seurat, "Park", col.name = "group")
park_F22_seurat <- AddMetaData(park_F22_seurat, "Park", col.name = "group")
park_F23_seurat <- AddMetaData(park_F23_seurat, "Park", col.name = "group")
park_F29_seurat <- AddMetaData(park_F29_seurat, "Park", col.name = "group")
park_F30_seurat <- AddMetaData(park_F30_seurat, "Park", col.name = "group")
park_F38_seurat <- AddMetaData(park_F38_seurat, "Park", col.name = "group") 
park_F41_seurat <- AddMetaData(park_F41_seurat, "Park", col.name = "group")
park_F45_seurat <- AddMetaData(park_F45_seurat, "Park", col.name = "group")
park_F64_seurat <- AddMetaData(park_F64_seurat, "Park", col.name = "group")
park_F67_seurat <- AddMetaData(park_F67_seurat, "Park", col.name = "group") 
park_F74_seurat <- AddMetaData(park_F74_seurat, "Park", col.name = "group")
park_F83_seurat <- AddMetaData(park_F83_seurat, "Park", col.name = "group")
park_P1_seurat <- AddMetaData(park_P1_seurat, "Park", col.name = "group")
park_P2_seurat <- AddMetaData(park_P2_seurat, "Park", col.name = "group")
park_P3_seurat <- AddMetaData(park_P3_seurat, "Park", col.name = "group")
park_T03_seurat <- AddMetaData(park_T03_seurat, "Park", col.name = "group")
park_T06_seurat <- AddMetaData(park_T06_seurat, "Park", col.name = "group")
park_T07_seurat <- AddMetaData(park_T07_seurat, "Park", col.name = "group")
bautista_a25_seurat <- AddMetaData(bautista_a25_seurat, "Bautista", col.name = "group")
bautista_f19_seurat <- AddMetaData(bautista_f19_seurat, "Bautista", col.name = "group")
bautista_f231_seurat <- AddMetaData(bautista_f231_seurat, "Bautista", col.name = "group")
bautista_f232_seurat <- AddMetaData(bautista_f232_seurat, "Bautista", col.name = "group")
bautista_p6_seurat <- AddMetaData(bautista_p6_seurat, "Bautista", col.name = "group")
bautista_p101_seurat <- AddMetaData(bautista_p101_seurat, "Bautista", col.name = "group")
bautista_p102_seurat <- AddMetaData(bautista_p102_seurat, "Bautista", col.name = "group")
direder_TH1_seurat <- AddMetaData(direder_TH1_seurat, "Direder", col.name = "group")
direder_TH2_seurat <- AddMetaData(direder_TH2_seurat, "Direder", col.name = "group")
direder_TH3_seurat <- AddMetaData(direder_TH3_seurat, "Direder", col.name = "group")
direder_TH4_seurat <- AddMetaData(direder_TH4_seurat, "Direder", col.name = "group")
direder_TB1_seurat <- AddMetaData(direder_TB1_seurat, "Direder", col.name = "group")
direder_TB2_seurat <- AddMetaData(direder_TB2_seurat, "Direder", col.name = "group")
direder_TB3_seurat <- AddMetaData(direder_TB3_seurat, "Direder", col.name = "group")
direder_TB4_seurat <- AddMetaData(direder_TB4_seurat, "Direder", col.name = "group")
direder_TC1_seurat <- AddMetaData(direder_TC1_seurat, "Direder", col.name = "group")
direder_TC2_seurat <- AddMetaData(direder_TC2_seurat, "Direder", col.name = "group")
direder_TH_MG1_seurat <- AddMetaData(direder_TH_MG1_seurat, "Direder", col.name = "group")
direder_TH_MG2_seurat <- AddMetaData(direder_TH_MG2_seurat, "Direder", col.name = "group")

#(5)_flag-2_Donor_ID####
xin1_seurat <- AddMetaData(xin1_seurat, "T01_1", col.name = "Donor_ID")
xin2_seurat <- AddMetaData(xin2_seurat, "T02_1", col.name = "Donor_ID")
xin3_seurat <- AddMetaData(xin3_seurat, "T03_1", col.name = "Donor_ID")
xin4_seurat <- AddMetaData(xin4_seurat, "T04_1", col.name = "Donor_ID")
xin5_seurat <- AddMetaData(xin5_seurat, "T05_1", col.name = "Donor_ID")
xin6_seurat <- AddMetaData(xin6_seurat, "T06_1", col.name = "Donor_ID")
yasumizu_MG03_seurat <- AddMetaData(yasumizu_MG03_seurat, "MG03", col.name = "Donor_ID")
yasumizu_MG21_seurat <- AddMetaData(yasumizu_MG21_seurat, "MG21", col.name = "Donor_ID")
yasumizu_MG22_seurat <- AddMetaData(yasumizu_MG22_seurat, "MG22", col.name = "Donor_ID")
yasumizu_MG23_seurat <- AddMetaData(yasumizu_MG23_seurat, "MG23", col.name = "Donor_ID")
park_A16_seurat <- AddMetaData(park_A16_seurat, "A16", col.name = "Donor_ID")
park_A43_seurat <- AddMetaData(park_A43_seurat, "A43", col.name = "Donor_ID") 
park_A45_seurat <- AddMetaData(park_A45_seurat, "A45", col.name = "Donor_ID")
park_C34_seurat <- AddMetaData(park_C34_seurat, "C34", col.name = "Donor_ID")
park_C40_seurat <- AddMetaData(park_C40_seurat, "C40", col.name = "Donor_ID") 
park_C41_seurat <- AddMetaData(park_C41_seurat, "C41", col.name = "Donor_ID")
park_F21_seurat <- AddMetaData(park_F21_seurat, "F21", col.name = "Donor_ID")
park_F22_seurat <- AddMetaData(park_F22_seurat, "F22", col.name = "Donor_ID")
park_F23_seurat <- AddMetaData(park_F23_seurat, "F23", col.name = "Donor_ID")
park_F29_seurat <- AddMetaData(park_F29_seurat, "F29", col.name = "Donor_ID")
park_F30_seurat <- AddMetaData(park_F30_seurat, "F30", col.name = "Donor_ID")
park_F38_seurat <- AddMetaData(park_F38_seurat, "F38", col.name = "Donor_ID") 
park_F41_seurat <- AddMetaData(park_F41_seurat, "F41", col.name = "Donor_ID")
park_F45_seurat <- AddMetaData(park_F45_seurat, "F45", col.name = "Donor_ID")
park_F64_seurat <- AddMetaData(park_F64_seurat, "F64", col.name = "Donor_ID")
park_F67_seurat <- AddMetaData(park_F67_seurat, "F67", col.name = "Donor_ID") 
park_F74_seurat <- AddMetaData(park_F74_seurat, "F74", col.name = "Donor_ID")
park_F83_seurat <- AddMetaData(park_F83_seurat, "F83", col.name = "Donor_ID")
park_P1_seurat <- AddMetaData(park_P1_seurat, "P1", col.name = "Donor_ID")
park_P2_seurat <- AddMetaData(park_P2_seurat, "P2", col.name = "Donor_ID")
park_P3_seurat <- AddMetaData(park_P3_seurat, "P3", col.name = "Donor_ID")
park_T03_seurat <- AddMetaData(park_T03_seurat, "T03", col.name = "Donor_ID")
park_T06_seurat <- AddMetaData(park_T06_seurat, "T06", col.name = "Donor_ID")
park_T07_seurat <- AddMetaData(park_T07_seurat, "T07", col.name = "Donor_ID")
bautista_a25_seurat <- AddMetaData(bautista_a25_seurat, "A25", col.name = "Donor_ID")
bautista_f19_seurat <- AddMetaData(bautista_f19_seurat, "F19", col.name = "Donor_ID")
bautista_f231_seurat <- AddMetaData(bautista_f231_seurat, "F234", col.name = "Donor_ID")
bautista_f232_seurat <- AddMetaData(bautista_f232_seurat, "F232", col.name = "Donor_ID")
bautista_p6_seurat <- AddMetaData(bautista_p6_seurat, "P6", col.name = "Donor_ID")
bautista_p101_seurat <- AddMetaData(bautista_p101_seurat, "P101", col.name = "Donor_ID")
bautista_p102_seurat <- AddMetaData(bautista_p102_seurat, "P102", col.name = "Donor_ID")
direder_TH1_seurat <- AddMetaData(direder_TH1_seurat, "TH1", col.name = "Donor_ID")
direder_TH2_seurat <- AddMetaData(direder_TH2_seurat, "TH2", col.name = "Donor_ID")
direder_TH3_seurat <- AddMetaData(direder_TH3_seurat, "TH3", col.name = "Donor_ID")
direder_TH4_seurat <- AddMetaData(direder_TH4_seurat, "TH4", col.name = "Donor_ID")
direder_TB1_seurat <- AddMetaData(direder_TB1_seurat, "TB1", col.name = "Donor_ID")
direder_TB2_seurat <- AddMetaData(direder_TB2_seurat, "TB2", col.name = "Donor_ID")
direder_TB3_seurat <- AddMetaData(direder_TB3_seurat, "TB3", col.name = "Donor_ID")
direder_TB4_seurat <- AddMetaData(direder_TB4_seurat, "TB4", col.name = "Donor_ID")
direder_TC1_seurat <- AddMetaData(direder_TC1_seurat, "TC1", col.name = "Donor_ID")
direder_TC2_seurat <- AddMetaData(direder_TC2_seurat, "TC2", col.name = "Donor_ID")
direder_TH_MG1_seurat <- AddMetaData(direder_TH_MG1_seurat, "TH_MG1", col.name = "Donor_ID")
direder_TH_MG2_seurat <- AddMetaData(direder_TH_MG2_seurat, "TH_MG2", col.name = "Donor_ID")

#(6)_flag-3_condition####
xin1_seurat <- AddMetaData(xin1_seurat, "TET_AB", col.name = "condition")
xin2_seurat <- AddMetaData(xin2_seurat, "TET_C", col.name = "condition")
xin3_seurat <- AddMetaData(xin3_seurat, "TET_A", col.name = "condition")
xin4_seurat <- AddMetaData(xin4_seurat, "TET_AB", col.name = "condition")
xin5_seurat <- AddMetaData(xin5_seurat, "TET_MNT", col.name = "condition")
xin6_seurat <- AddMetaData(xin6_seurat, "TET_B", col.name = "condition")
yasumizu_MG03_seurat <- AddMetaData(yasumizu_MG03_seurat, "TET_AB", col.name = "condition")
yasumizu_MG21_seurat <- AddMetaData(yasumizu_MG21_seurat, "TET_B", col.name = "condition")
yasumizu_MG22_seurat <- AddMetaData(yasumizu_MG22_seurat, "TET_B", col.name = "condition")
yasumizu_MG23_seurat <- AddMetaData(yasumizu_MG23_seurat, "TET_AB", col.name = "condition")
park_A16_seurat <- AddMetaData(park_A16_seurat, "Thymus", col.name = "condition")
park_A43_seurat <- AddMetaData(park_A43_seurat, "Thymus", col.name = "condition") 
park_A45_seurat <- AddMetaData(park_A45_seurat, "Thymus", col.name = "condition")
park_C34_seurat <- AddMetaData(park_C34_seurat, "Thymus", col.name = "condition")
park_C40_seurat <- AddMetaData(park_C40_seurat, "Thymus", col.name = "condition") 
park_C41_seurat <- AddMetaData(park_C41_seurat, "Thymus", col.name = "condition")
park_F21_seurat <- AddMetaData(park_F21_seurat, "Thymus", col.name = "condition")
park_F22_seurat <- AddMetaData(park_F22_seurat, "Thymus", col.name = "condition")
park_F23_seurat <- AddMetaData(park_F23_seurat, "Thymus", col.name = "condition")
park_F29_seurat <- AddMetaData(park_F29_seurat, "Thymus", col.name = "condition")
park_F30_seurat <- AddMetaData(park_F30_seurat, "Thymus", col.name = "condition")
park_F38_seurat <- AddMetaData(park_F38_seurat, "Thymus", col.name = "condition") 
park_F41_seurat <- AddMetaData(park_F41_seurat, "Thymus", col.name = "condition")
park_F45_seurat <- AddMetaData(park_F45_seurat, "Thymus", col.name = "condition")
park_F64_seurat <- AddMetaData(park_F64_seurat, "Thymus", col.name = "condition")
park_F67_seurat <- AddMetaData(park_F67_seurat, "Thymus", col.name = "condition") 
park_F74_seurat <- AddMetaData(park_F74_seurat, "Thymus", col.name = "condition")
park_F83_seurat <- AddMetaData(park_F83_seurat, "Thymus", col.name = "condition")
park_P1_seurat <- AddMetaData(park_P1_seurat, "Thymus", col.name = "condition")
park_P2_seurat <- AddMetaData(park_P2_seurat, "Thymus", col.name = "condition")
park_P3_seurat <- AddMetaData(park_P3_seurat, "Thymus", col.name = "condition")
park_T03_seurat <- AddMetaData(park_T03_seurat, "Thymus", col.name = "condition")
park_T06_seurat <- AddMetaData(park_T06_seurat, "Thymus", col.name = "condition")
park_T07_seurat <- AddMetaData(park_T07_seurat, "Thymus", col.name = "condition")
bautista_a25_seurat <- AddMetaData(bautista_a25_seurat, "Thymus", col.name = "condition")
bautista_f19_seurat <- AddMetaData(bautista_f19_seurat, "Thymus", col.name = "condition")
bautista_f231_seurat <- AddMetaData(bautista_f231_seurat, "Thymus", col.name = "condition")
bautista_f232_seurat <- AddMetaData(bautista_f232_seurat, "Thymus", col.name = "condition")
bautista_p6_seurat <- AddMetaData(bautista_p6_seurat, "Thymus", col.name = "condition")
bautista_p101_seurat <- AddMetaData(bautista_p101_seurat, "Thymus", col.name = "condition")
bautista_p102_seurat <- AddMetaData(bautista_p102_seurat, "Thymus", col.name = "condition")
direder_TH1_seurat <- AddMetaData(direder_TH1_seurat, "TTH", col.name = "condition")
direder_TH2_seurat <- AddMetaData(direder_TH2_seurat, "TTH", col.name = "condition")
direder_TH3_seurat <- AddMetaData(direder_TH3_seurat, "TTH", col.name = "condition")
direder_TH4_seurat <- AddMetaData(direder_TH4_seurat, "TTH", col.name = "condition")
direder_TB1_seurat <- AddMetaData(direder_TB1_seurat, "TET_B", col.name = "condition")
direder_TB2_seurat <- AddMetaData(direder_TB2_seurat, "TET_B", col.name = "condition")
direder_TB3_seurat <- AddMetaData(direder_TB3_seurat, "TET_B", col.name = "condition")
direder_TB4_seurat <- AddMetaData(direder_TB4_seurat, "TET_B", col.name = "condition")
direder_TC1_seurat <- AddMetaData(direder_TC1_seurat, "TET_C", col.name = "condition")
direder_TC2_seurat <- AddMetaData(direder_TC2_seurat, "TET_C", col.name = "condition")
direder_TH_MG1_seurat <- AddMetaData(direder_TH_MG1_seurat, "TTH", col.name = "condition")
direder_TH_MG2_seurat <- AddMetaData(direder_TH_MG2_seurat, "TTH", col.name = "condition")

#(7)_flag-4_classification####
xin1_seurat <- AddMetaData(xin1_seurat, "TET_AB", col.name = "classification")
xin2_seurat <- AddMetaData(xin2_seurat, "TET_C", col.name = "classification")
xin3_seurat <- AddMetaData(xin3_seurat, "TET_A", col.name = "classification")
xin4_seurat <- AddMetaData(xin4_seurat, "TET_AB", col.name = "classification")
xin5_seurat <- AddMetaData(xin5_seurat, "TET_MNT", col.name = "classification")
xin6_seurat <- AddMetaData(xin6_seurat, "TET_B3", col.name = "classification")
yasumizu_MG03_seurat <- AddMetaData(yasumizu_MG03_seurat, "TET_AB", col.name = "classification")
yasumizu_MG21_seurat <- AddMetaData(yasumizu_MG21_seurat, "TET_B1", col.name = "classification")
yasumizu_MG22_seurat <- AddMetaData(yasumizu_MG22_seurat, "TET_B2", col.name = "classification")
yasumizu_MG23_seurat <- AddMetaData(yasumizu_MG23_seurat, "TET_AB", col.name = "classification")
park_A16_seurat <- AddMetaData(park_A16_seurat, "Thymus", col.name = "classification")
park_A43_seurat <- AddMetaData(park_A43_seurat, "Thymus", col.name = "classification") 
park_A45_seurat <- AddMetaData(park_A45_seurat, "Thymus", col.name = "classification")
park_C34_seurat <- AddMetaData(park_C34_seurat, "Thymus", col.name = "classification")
park_C40_seurat <- AddMetaData(park_C40_seurat, "Thymus", col.name = "classification") 
park_C41_seurat <- AddMetaData(park_C41_seurat, "Thymus", col.name = "classification")
park_F21_seurat <- AddMetaData(park_F21_seurat, "Thymus", col.name = "classification")
park_F22_seurat <- AddMetaData(park_F22_seurat, "Thymus", col.name = "classification")
park_F23_seurat <- AddMetaData(park_F23_seurat, "Thymus", col.name = "classification")
park_F29_seurat <- AddMetaData(park_F29_seurat, "Thymus", col.name = "classification")
park_F30_seurat <- AddMetaData(park_F30_seurat, "Thymus", col.name = "classification")
park_F38_seurat <- AddMetaData(park_F38_seurat, "Thymus", col.name = "classification") 
park_F41_seurat <- AddMetaData(park_F41_seurat, "Thymus", col.name = "classification")
park_F45_seurat <- AddMetaData(park_F45_seurat, "Thymus", col.name = "classification")
park_F64_seurat <- AddMetaData(park_F64_seurat, "Thymus", col.name = "classification")
park_F67_seurat <- AddMetaData(park_F67_seurat, "Thymus", col.name = "classification") 
park_F74_seurat <- AddMetaData(park_F74_seurat, "Thymus", col.name = "classification")
park_F83_seurat <- AddMetaData(park_F83_seurat, "Thymus", col.name = "classification")
park_P1_seurat <- AddMetaData(park_P1_seurat, "Thymus", col.name = "classification")
park_P2_seurat <- AddMetaData(park_P2_seurat, "Thymus", col.name = "classification")
park_P3_seurat <- AddMetaData(park_P3_seurat, "Thymus", col.name = "classification")
park_T03_seurat <- AddMetaData(park_T03_seurat, "Thymus", col.name = "classification")
park_T06_seurat <- AddMetaData(park_T06_seurat, "Thymus", col.name = "classification")
park_T07_seurat <- AddMetaData(park_T07_seurat, "Thymus", col.name = "classification")
bautista_a25_seurat <- AddMetaData(bautista_a25_seurat, "Thymus", col.name = "classification")
bautista_f19_seurat <- AddMetaData(bautista_f19_seurat, "Thymus", col.name = "classification")
bautista_f231_seurat <- AddMetaData(bautista_f231_seurat, "Thymus", col.name = "classification")
bautista_f232_seurat <- AddMetaData(bautista_f232_seurat, "Thymus", col.name = "classification")
bautista_p6_seurat <- AddMetaData(bautista_p6_seurat, "Thymus", col.name = "classification")
bautista_p101_seurat <- AddMetaData(bautista_p101_seurat, "Thymus", col.name = "classification")
bautista_p102_seurat <- AddMetaData(bautista_p102_seurat, "Thymus", col.name = "classification")
direder_TH1_seurat <- AddMetaData(direder_TH1_seurat, "TTH", col.name = "classification")
direder_TH2_seurat <- AddMetaData(direder_TH2_seurat, "TTH", col.name = "classification")
direder_TH3_seurat <- AddMetaData(direder_TH3_seurat, "TTH", col.name = "classification")
direder_TH4_seurat <- AddMetaData(direder_TH4_seurat, "TTH", col.name = "classification")
direder_TB1_seurat <- AddMetaData(direder_TB1_seurat, "TET_B1", col.name = "classification")
direder_TB2_seurat <- AddMetaData(direder_TB2_seurat, "TET_B1", col.name = "classification")
direder_TB3_seurat <- AddMetaData(direder_TB3_seurat, "TET_B2", col.name = "classification")
direder_TB4_seurat <- AddMetaData(direder_TB4_seurat, "TET_B2_B3", col.name = "classification")
direder_TC1_seurat <- AddMetaData(direder_TC1_seurat, "TET_C", col.name = "classification")
direder_TC2_seurat <- AddMetaData(direder_TC2_seurat, "TET_C", col.name = "classification")
direder_TH_MG1_seurat <- AddMetaData(direder_TH_MG1_seurat, "TTH", col.name = "classification")
direder_TH_MG2_seurat <- AddMetaData(direder_TH_MG2_seurat, "TTH", col.name = "classification")

#(8)_flag-5_tissue####
xin1_seurat <- AddMetaData(xin1_seurat, "TET_AB", col.name = "tissue")
xin2_seurat <- AddMetaData(xin2_seurat, "TET_C", col.name = "tissue")
xin3_seurat <- AddMetaData(xin3_seurat, "TET_A", col.name = "tissue")
xin4_seurat <- AddMetaData(xin4_seurat, "TET_AB", col.name = "tissue")
xin5_seurat <- AddMetaData(xin5_seurat, "TET_MNT", col.name = "tissue")
xin6_seurat <- AddMetaData(xin6_seurat, "TET_B", col.name = "tissue")
yasumizu_MG03_seurat <- AddMetaData(yasumizu_MG03_seurat, "TET_AB", col.name = "tissue")
yasumizu_MG21_seurat <- AddMetaData(yasumizu_MG21_seurat, "TET_B", col.name = "tissue")
yasumizu_MG22_seurat <- AddMetaData(yasumizu_MG22_seurat, "TET_B", col.name = "tissue")
yasumizu_MG23_seurat <- AddMetaData(yasumizu_MG23_seurat, "TET_AB", col.name = "tissue")
park_A16_seurat <- AddMetaData(park_A16_seurat, "adult_Thymus", col.name = "tissue")
park_A43_seurat <- AddMetaData(park_A43_seurat, "adult_Thymus", col.name = "tissue") 
park_A45_seurat <- AddMetaData(park_A45_seurat, "adult_Thymus", col.name = "tissue")
park_C34_seurat <- AddMetaData(park_C34_seurat, "prenatal_Thymus", col.name = "tissue")
park_C40_seurat <- AddMetaData(park_C40_seurat, "prenatal_Thymus", col.name = "tissue") 
park_C41_seurat <- AddMetaData(park_C41_seurat, "prenatal_Thymus", col.name = "tissue")
park_F21_seurat <- AddMetaData(park_F21_seurat, "prenatal_Thymus", col.name = "tissue")
park_F22_seurat <- AddMetaData(park_F22_seurat, "prenatal_Thymus", col.name = "tissue")
park_F23_seurat <- AddMetaData(park_F23_seurat, "prenatal_Thymus", col.name = "tissue")
park_F29_seurat <- AddMetaData(park_F29_seurat, "prenatal_Thymus", col.name = "tissue")
park_F30_seurat <- AddMetaData(park_F30_seurat, "prenatal_Thymus", col.name = "tissue")
park_F38_seurat <- AddMetaData(park_F38_seurat, "prenatal_Thymus", col.name = "tissue") 
park_F41_seurat <- AddMetaData(park_F41_seurat, "prenatal_Thymus", col.name = "tissue")
park_F45_seurat <- AddMetaData(park_F45_seurat, "prenatal_Thymus", col.name = "tissue")
park_F64_seurat <- AddMetaData(park_F64_seurat, "prenatal_Thymus", col.name = "tissue")
park_F67_seurat <- AddMetaData(park_F67_seurat, "prenatal_Thymus", col.name = "tissue") 
park_F74_seurat <- AddMetaData(park_F74_seurat, "prenatal_Thymus", col.name = "tissue")
park_F83_seurat <- AddMetaData(park_F83_seurat, "prenatal_Thymus", col.name = "tissue")
park_P1_seurat <- AddMetaData(park_P1_seurat, "pediatric_Thymus", col.name = "tissue")
park_P2_seurat <- AddMetaData(park_P2_seurat, "adult_Thymus", col.name = "tissue")
park_P3_seurat <- AddMetaData(park_P3_seurat, "pediatric_Thymus", col.name = "tissue")
park_T03_seurat <- AddMetaData(park_T03_seurat, "pediatric_Thymus", col.name = "tissue")
park_T06_seurat <- AddMetaData(park_T06_seurat, "pediatric_Thymus", col.name = "tissue")
park_T07_seurat <- AddMetaData(park_T07_seurat, "pediatric_Thymus", col.name = "tissue")
bautista_a25_seurat <- AddMetaData(bautista_a25_seurat, "adult_Thymus", col.name = "tissue")
bautista_f19_seurat <- AddMetaData(bautista_f19_seurat, "prenatal_Thymus", col.name = "tissue")
bautista_f231_seurat <- AddMetaData(bautista_f231_seurat, "prenatal_Thymus", col.name = "tissue")
bautista_f232_seurat <- AddMetaData(bautista_f232_seurat, "prenatal_Thymus", col.name = "tissue")
bautista_p6_seurat <- AddMetaData(bautista_p6_seurat, "pediatric_Thymus", col.name = "tissue")
bautista_p101_seurat <- AddMetaData(bautista_p101_seurat, "pediatric_Thymus", col.name = "tissue")
bautista_p102_seurat <- AddMetaData(bautista_p102_seurat, "pediatric_Thymus", col.name = "tissue")
direder_TH1_seurat <- AddMetaData(direder_TH1_seurat, "TTH", col.name = "tissue")
direder_TH2_seurat <- AddMetaData(direder_TH2_seurat, "TTH", col.name = "tissue")
direder_TH3_seurat <- AddMetaData(direder_TH3_seurat, "TTH", col.name = "tissue")
direder_TH4_seurat <- AddMetaData(direder_TH4_seurat, "TTH", col.name = "tissue")
direder_TB1_seurat <- AddMetaData(direder_TB1_seurat, "TET_B", col.name = "tissue")
direder_TB2_seurat <- AddMetaData(direder_TB2_seurat, "TET_B", col.name = "tissue")
direder_TB3_seurat <- AddMetaData(direder_TB3_seurat, "TET_B", col.name = "tissue")
direder_TB4_seurat <- AddMetaData(direder_TB4_seurat, "TET_B", col.name = "tissue")
direder_TC1_seurat <- AddMetaData(direder_TC1_seurat, "TET_C", col.name = "tissue")
direder_TC2_seurat <- AddMetaData(direder_TC2_seurat, "TET_C", col.name = "tissue")
direder_TH_MG1_seurat <- AddMetaData(direder_TH_MG1_seurat, "TTH", col.name = "tissue")
direder_TH_MG2_seurat <- AddMetaData(direder_TH_MG2_seurat, "TTH", col.name = "tissue")

#(9)_flag-6_sex####
xin1_seurat <- AddMetaData(xin1_seurat, "male", col.name = "sex")
xin2_seurat <- AddMetaData(xin2_seurat, "male", col.name = "sex")
xin3_seurat <- AddMetaData(xin3_seurat, "male", col.name = "sex")
xin4_seurat <- AddMetaData(xin4_seurat, "female", col.name = "sex")
xin5_seurat <- AddMetaData(xin5_seurat, "male", col.name = "sex")
xin6_seurat <- AddMetaData(xin6_seurat, "female", col.name = "sex")
yasumizu_MG03_seurat <- AddMetaData(yasumizu_MG03_seurat, "male", col.name = "sex")
yasumizu_MG21_seurat <- AddMetaData(yasumizu_MG21_seurat, "female", col.name = "sex")
yasumizu_MG22_seurat <- AddMetaData(yasumizu_MG22_seurat, "female", col.name = "sex")
yasumizu_MG23_seurat <- AddMetaData(yasumizu_MG23_seurat, "female", col.name = "sex")
park_A16_seurat <- AddMetaData(park_A16_seurat, "female", col.name = "sex")
park_A43_seurat <- AddMetaData(park_A43_seurat, "male", col.name = "sex") 
park_A45_seurat <- AddMetaData(park_A45_seurat, "female", col.name = "sex")
park_C34_seurat <- AddMetaData(park_C34_seurat, "female", col.name = "sex")
park_C40_seurat <- AddMetaData(park_C40_seurat, "female", col.name = "sex") 
park_C41_seurat <- AddMetaData(park_C41_seurat, "male", col.name = "sex")
park_F21_seurat <- AddMetaData(park_F21_seurat, "male", col.name = "sex")
park_F22_seurat <- AddMetaData(park_F22_seurat, "female", col.name = "sex")
park_F23_seurat <- AddMetaData(park_F23_seurat, "male", col.name = "sex")
park_F29_seurat <- AddMetaData(park_F29_seurat, "female", col.name = "sex")
park_F30_seurat <- AddMetaData(park_F30_seurat, "male", col.name = "sex")
park_F38_seurat <- AddMetaData(park_F38_seurat, "male", col.name = "sex") 
park_F41_seurat <- AddMetaData(park_F41_seurat, "female", col.name = "sex")
park_F45_seurat <- AddMetaData(park_F45_seurat, "female", col.name = "sex")
park_F64_seurat <- AddMetaData(park_F64_seurat, "male", col.name = "sex")
park_F67_seurat <- AddMetaData(park_F67_seurat, "female", col.name = "sex") 
park_F74_seurat <- AddMetaData(park_F74_seurat, "male", col.name = "sex")
park_F83_seurat <- AddMetaData(park_F83_seurat, "female", col.name = "sex")
park_P1_seurat <- AddMetaData(park_P1_seurat, "male", col.name = "sex")
park_P2_seurat <- AddMetaData(park_P2_seurat, "male", col.name = "sex")
park_P3_seurat <- AddMetaData(park_P3_seurat, "male", col.name = "sex")
park_T03_seurat <- AddMetaData(park_T03_seurat, "female", col.name = "sex")
park_T06_seurat <- AddMetaData(park_T06_seurat, "male", col.name = "sex")
park_T07_seurat <- AddMetaData(park_T07_seurat, "male", col.name = "sex")
bautista_a25_seurat <- AddMetaData(bautista_a25_seurat, "n.a.", col.name = "sex")
bautista_f19_seurat <- AddMetaData(bautista_f19_seurat, "n.a.", col.name = "sex")
bautista_f231_seurat <- AddMetaData(bautista_f231_seurat, "n.a.", col.name = "sex")
bautista_f232_seurat <- AddMetaData(bautista_f232_seurat, "n.a.", col.name = "sex")
bautista_p6_seurat <- AddMetaData(bautista_p6_seurat, "n.a.", col.name = "sex")
bautista_p101_seurat <- AddMetaData(bautista_p101_seurat, "n.a.", col.name = "sex")
bautista_p102_seurat <- AddMetaData(bautista_p102_seurat, "n.a.", col.name = "sex")
direder_TH1_seurat <- AddMetaData(direder_TH1_seurat, "male", col.name = "sex")
direder_TH2_seurat <- AddMetaData(direder_TH2_seurat, "male", col.name = "sex")
direder_TH3_seurat <- AddMetaData(direder_TH3_seurat, "male", col.name = "sex")
direder_TH4_seurat <- AddMetaData(direder_TH4_seurat, "female", col.name = "sex")
direder_TB1_seurat <- AddMetaData(direder_TB1_seurat, "male", col.name = "sex")
direder_TB2_seurat <- AddMetaData(direder_TB2_seurat, "male", col.name = "sex")
direder_TB3_seurat <- AddMetaData(direder_TB3_seurat, "male", col.name = "sex")
direder_TB4_seurat <- AddMetaData(direder_TB4_seurat, "male", col.name = "sex")
direder_TC1_seurat <- AddMetaData(direder_TC1_seurat, "female", col.name = "sex")
direder_TC2_seurat <- AddMetaData(direder_TC2_seurat, "male", col.name = "sex")
direder_TH_MG1_seurat <- AddMetaData(direder_TH_MG1_seurat, "male", col.name = "sex")
direder_TH_MG2_seurat <- AddMetaData(direder_TH_MG2_seurat, "female", col.name = "sex")

#(10)_flag-7_age####
xin1_seurat <- AddMetaData(xin1_seurat, "56y", col.name = "age")
xin2_seurat <- AddMetaData(xin2_seurat, "69y", col.name = "age")
xin3_seurat <- AddMetaData(xin3_seurat, "72y", col.name = "age")
xin4_seurat <- AddMetaData(xin4_seurat, "55y", col.name = "age")
xin5_seurat <- AddMetaData(xin5_seurat, "67y", col.name = "age")
xin6_seurat <- AddMetaData(xin6_seurat, "44y", col.name = "age")
yasumizu_MG03_seurat <- AddMetaData(yasumizu_MG03_seurat, "46y", col.name = "age")
yasumizu_MG21_seurat <- AddMetaData(yasumizu_MG21_seurat, "55y", col.name = "age")
yasumizu_MG22_seurat <- AddMetaData(yasumizu_MG22_seurat, "49y", col.name = "age")
yasumizu_MG23_seurat <- AddMetaData(yasumizu_MG23_seurat, "35y", col.name = "age")
park_A16_seurat <- AddMetaData(park_A16_seurat, "20-25y", col.name = "age")
park_A43_seurat <- AddMetaData(park_A43_seurat, "35-40y", col.name = "age") 
park_A45_seurat <- AddMetaData(park_A45_seurat, "15-20y", col.name = "age")
park_C34_seurat <- AddMetaData(park_C34_seurat, "9w", col.name = "age")
park_C40_seurat <- AddMetaData(park_C40_seurat, "7w", col.name = "age") 
park_C41_seurat <- AddMetaData(park_C41_seurat, "8w", col.name = "age")
park_F21_seurat <- AddMetaData(park_F21_seurat, "16w", col.name = "age")
park_F22_seurat <- AddMetaData(park_F22_seurat, "9w", col.name = "age")
park_F23_seurat <- AddMetaData(park_F23_seurat, "11w", col.name = "age")
park_F29_seurat <- AddMetaData(park_F29_seurat, "17w", col.name = "age")
park_F30_seurat <- AddMetaData(park_F30_seurat, "14w", col.name = "age")
park_F38_seurat <- AddMetaData(park_F38_seurat, "13w", col.name = "age") 
park_F41_seurat <- AddMetaData(park_F41_seurat, "16w", col.name = "age")
park_F45_seurat <- AddMetaData(park_F45_seurat, "12w", col.name = "age")
park_F64_seurat <- AddMetaData(park_F64_seurat, "11w", col.name = "age")
park_F67_seurat <- AddMetaData(park_F67_seurat, "12w", col.name = "age") 
park_F74_seurat <- AddMetaData(park_F74_seurat, "10w", col.name = "age")
park_F83_seurat <- AddMetaData(park_F83_seurat, "17w", col.name = "age")
park_P1_seurat <- AddMetaData(park_P1_seurat, "15m", col.name = "age")
park_P2_seurat <- AddMetaData(park_P2_seurat, "10-15y", col.name = "age")
park_P3_seurat <- AddMetaData(park_P3_seurat, "6m", col.name = "age")
park_T03_seurat <- AddMetaData(park_T03_seurat, "10m", col.name = "age")
park_T06_seurat <- AddMetaData(park_T06_seurat, "30m", col.name = "age")
park_T07_seurat <- AddMetaData(park_T07_seurat, "3m", col.name = "age")
bautista_a25_seurat <- AddMetaData(bautista_a25_seurat, "25y", col.name = "age")
bautista_f19_seurat <- AddMetaData(bautista_f19_seurat, "19m", col.name = "age")
bautista_f231_seurat <- AddMetaData(bautista_f231_seurat, "23m", col.name = "age")
bautista_f232_seurat <- AddMetaData(bautista_f232_seurat, "23m", col.name = "age")
bautista_p6_seurat <- AddMetaData(bautista_p6_seurat, "6m", col.name = "age")
bautista_p101_seurat <- AddMetaData(bautista_p101_seurat, "10m", col.name = "age")
bautista_p102_seurat <- AddMetaData(bautista_p102_seurat, "10m", col.name = "age")
direder_TH1_seurat <- AddMetaData(direder_TH1_seurat, "70y", col.name = "age")
direder_TH2_seurat <- AddMetaData(direder_TH2_seurat, "40y", col.name = "age")
direder_TH3_seurat <- AddMetaData(direder_TH3_seurat, "25y", col.name = "age")
direder_TH4_seurat <- AddMetaData(direder_TH4_seurat, "54y", col.name = "age")
direder_TB1_seurat <- AddMetaData(direder_TB1_seurat, "70y", col.name = "age")
direder_TB2_seurat <- AddMetaData(direder_TB2_seurat, "40y", col.name = "age")
direder_TB3_seurat <- AddMetaData(direder_TB3_seurat, "68y", col.name = "age")
direder_TB4_seurat <- AddMetaData(direder_TB4_seurat, "77y", col.name = "age")
direder_TC1_seurat <- AddMetaData(direder_TC1_seurat, "54y", col.name = "age")
direder_TC2_seurat <- AddMetaData(direder_TC2_seurat, "76y", col.name = "age")
direder_TH_MG1_seurat <- AddMetaData(direder_TH_MG1_seurat, "23y", col.name = "age")
direder_TH_MG2_seurat <- AddMetaData(direder_TH_MG2_seurat, "55y", col.name = "age")

#(11)_flag-8_development####
xin1_seurat <- AddMetaData(xin1_seurat, "adult", col.name = "development")
xin2_seurat <- AddMetaData(xin2_seurat, "adult", col.name = "development")
xin3_seurat <- AddMetaData(xin3_seurat, "adult", col.name = "development")
xin4_seurat <- AddMetaData(xin4_seurat, "adult", col.name = "development")
xin5_seurat <- AddMetaData(xin5_seurat, "adult", col.name = "development")
xin6_seurat <- AddMetaData(xin6_seurat, "adult", col.name = "development")
yasumizu_MG03_seurat <- AddMetaData(yasumizu_MG03_seurat, "adult", col.name = "development")
yasumizu_MG21_seurat <- AddMetaData(yasumizu_MG21_seurat, "adult", col.name = "development")
yasumizu_MG22_seurat <- AddMetaData(yasumizu_MG22_seurat, "adult", col.name = "development")
yasumizu_MG23_seurat <- AddMetaData(yasumizu_MG23_seurat, "adult", col.name = "development")
park_A16_seurat <- AddMetaData(park_A16_seurat, "adult", col.name = "development")
park_A43_seurat <- AddMetaData(park_A43_seurat, "adult", col.name = "development") 
park_A45_seurat <- AddMetaData(park_A45_seurat, "adult", col.name = "development")
park_C34_seurat <- AddMetaData(park_C34_seurat, "fetal", col.name = "development")
park_C40_seurat <- AddMetaData(park_C40_seurat, "embryo", col.name = "development") 
park_C41_seurat <- AddMetaData(park_C41_seurat, "embryo", col.name = "development")
park_F21_seurat <- AddMetaData(park_F21_seurat, "fetal", col.name = "development")
park_F22_seurat <- AddMetaData(park_F22_seurat, "fetal", col.name = "development")
park_F23_seurat <- AddMetaData(park_F23_seurat, "fetal", col.name = "development")
park_F29_seurat <- AddMetaData(park_F29_seurat, "fetal", col.name = "development")
park_F30_seurat <- AddMetaData(park_F30_seurat, "fetal", col.name = "development")
park_F38_seurat <- AddMetaData(park_F38_seurat, "fetal", col.name = "development") 
park_F41_seurat <- AddMetaData(park_F41_seurat, "fetal", col.name = "development")
park_F45_seurat <- AddMetaData(park_F45_seurat, "fetal", col.name = "development")
park_F64_seurat <- AddMetaData(park_F64_seurat, "fetal", col.name = "development")
park_F67_seurat <- AddMetaData(park_F67_seurat, "fetal", col.name = "development") 
park_F74_seurat <- AddMetaData(park_F74_seurat, "fetal", col.name = "development")
park_F83_seurat <- AddMetaData(park_F83_seurat, "fetal", col.name = "development")
park_P1_seurat <- AddMetaData(park_P1_seurat, "pediatric", col.name = "development")
park_P2_seurat <- AddMetaData(park_P2_seurat, "adult", col.name = "development")
park_P3_seurat <- AddMetaData(park_P3_seurat, "pediatric", col.name = "development")
park_T03_seurat <- AddMetaData(park_T03_seurat, "pediatric", col.name = "development")
park_T06_seurat <- AddMetaData(park_T06_seurat, "pediatric", col.name = "development")
park_T07_seurat <- AddMetaData(park_T07_seurat, "pediatric", col.name = "development")
bautista_a25_seurat <- AddMetaData(bautista_a25_seurat, "adult", col.name = "development")
bautista_f19_seurat <- AddMetaData(bautista_f19_seurat, "fetal", col.name = "development")
bautista_f231_seurat <- AddMetaData(bautista_f231_seurat, "fetal", col.name = "development")
bautista_f232_seurat <- AddMetaData(bautista_f232_seurat, "fetal", col.name = "development")
bautista_p6_seurat <- AddMetaData(bautista_p6_seurat, "pediatric", col.name = "development")
bautista_p101_seurat <- AddMetaData(bautista_p101_seurat, "pediatric", col.name = "development")
bautista_p102_seurat <- AddMetaData(bautista_p102_seurat, "pediatric", col.name = "development")
direder_TH1_seurat <- AddMetaData(direder_TH1_seurat, "adult", col.name = "development")
direder_TH2_seurat <- AddMetaData(direder_TH2_seurat, "adult", col.name = "development")
direder_TH3_seurat <- AddMetaData(direder_TH3_seurat, "adult", col.name = "development")
direder_TH4_seurat <- AddMetaData(direder_TH4_seurat, "adult", col.name = "development")
direder_TB1_seurat <- AddMetaData(direder_TB1_seurat, "adult", col.name = "development")
direder_TB2_seurat <- AddMetaData(direder_TB2_seurat, "adult", col.name = "development")
direder_TB3_seurat <- AddMetaData(direder_TB3_seurat, "adult", col.name = "development")
direder_TB4_seurat <- AddMetaData(direder_TB4_seurat, "adult", col.name = "development")
direder_TC1_seurat <- AddMetaData(direder_TC1_seurat, "adult", col.name = "development")
direder_TC2_seurat <- AddMetaData(direder_TC2_seurat, "adult", col.name = "development")
direder_TH_MG1_seurat <- AddMetaData(direder_TH_MG1_seurat, "adult", col.name = "development")
direder_TH_MG2_seurat <- AddMetaData(direder_TH_MG2_seurat, "adult", col.name = "development")

#(12)_flag-9_donor####
xin1_seurat <- AddMetaData(xin1_seurat, "Xin1", col.name = "donor")
xin2_seurat <- AddMetaData(xin2_seurat, "Xin2", col.name = "donor")
xin3_seurat <- AddMetaData(xin3_seurat, "Xin3", col.name = "donor")
xin4_seurat <- AddMetaData(xin4_seurat, "Xin4", col.name = "donor")
xin5_seurat <- AddMetaData(xin5_seurat, "Xin5", col.name = "donor")
xin6_seurat <- AddMetaData(xin6_seurat, "Xin6", col.name = "donor")
yasumizu_MG03_seurat <- AddMetaData(yasumizu_MG03_seurat, "Yasumizu1", col.name = "donor")
yasumizu_MG21_seurat <- AddMetaData(yasumizu_MG21_seurat, "Yasumizu2", col.name = "donor")
yasumizu_MG22_seurat <- AddMetaData(yasumizu_MG22_seurat, "Yasumizu3", col.name = "donor")
yasumizu_MG23_seurat <- AddMetaData(yasumizu_MG23_seurat, "Yasumizu4", col.name = "donor")
park_A16_seurat <- AddMetaData(park_A16_seurat, "Park1", col.name = "donor")
park_A43_seurat <- AddMetaData(park_A43_seurat, "Park2", col.name = "donor") 
park_A45_seurat <- AddMetaData(park_A45_seurat, "Park3", col.name = "donor")
park_C34_seurat <- AddMetaData(park_C34_seurat, "Park4", col.name = "donor")
park_C40_seurat <- AddMetaData(park_C40_seurat, "Park5", col.name = "donor") 
park_C41_seurat <- AddMetaData(park_C41_seurat, "Park6", col.name = "donor")
park_F21_seurat <- AddMetaData(park_F21_seurat, "Park7", col.name = "donor")
park_F22_seurat <- AddMetaData(park_F22_seurat, "Park8", col.name = "donor")
park_F23_seurat <- AddMetaData(park_F23_seurat, "Park9", col.name = "donor")
park_F29_seurat <- AddMetaData(park_F29_seurat, "Park10", col.name = "donor")
park_F30_seurat <- AddMetaData(park_F30_seurat, "Park11", col.name = "donor")
park_F38_seurat <- AddMetaData(park_F38_seurat, "Park12", col.name = "donor") 
park_F41_seurat <- AddMetaData(park_F41_seurat, "Park13", col.name = "donor")
park_F45_seurat <- AddMetaData(park_F45_seurat, "Park14", col.name = "donor")
park_F64_seurat <- AddMetaData(park_F64_seurat, "Park15", col.name = "donor")
park_F67_seurat <- AddMetaData(park_F67_seurat, "Park16", col.name = "donor") 
park_F74_seurat <- AddMetaData(park_F74_seurat, "Park17", col.name = "donor")
park_F83_seurat <- AddMetaData(park_F83_seurat, "Park18", col.name = "donor")
park_P1_seurat <- AddMetaData(park_P1_seurat, "Park19", col.name = "donor")
park_P2_seurat <- AddMetaData(park_P2_seurat, "Park20", col.name = "donor")
park_P3_seurat <- AddMetaData(park_P3_seurat, "Park21", col.name = "donor")
park_T03_seurat <- AddMetaData(park_T03_seurat, "Park22", col.name = "donor")
park_T06_seurat <- AddMetaData(park_T06_seurat, "Park23", col.name = "donor")
park_T07_seurat <- AddMetaData(park_T07_seurat, "Park24", col.name = "donor")
bautista_a25_seurat <- AddMetaData(bautista_a25_seurat, "Bautista1", col.name = "donor")
bautista_f19_seurat <- AddMetaData(bautista_f19_seurat, "Bautista2", col.name = "donor")
bautista_f231_seurat <- AddMetaData(bautista_f231_seurat, "Bautista3", col.name = "donor")
bautista_f232_seurat <- AddMetaData(bautista_f232_seurat, "Bautista4", col.name = "donor")
bautista_p6_seurat <- AddMetaData(bautista_p6_seurat, "Bautista5", col.name = "donor")
bautista_p101_seurat <- AddMetaData(bautista_p101_seurat, "Bautista6", col.name = "donor")
bautista_p102_seurat <- AddMetaData(bautista_p102_seurat, "Bautista7", col.name = "donor")
direder_TH1_seurat <- AddMetaData(direder_TH1_seurat, "Direder1", col.name = "donor")
direder_TH2_seurat <- AddMetaData(direder_TH2_seurat, "Direder2", col.name = "donor")
direder_TH3_seurat <- AddMetaData(direder_TH3_seurat, "Direder3", col.name = "donor")
direder_TH4_seurat <- AddMetaData(direder_TH4_seurat, "Direder4", col.name = "donor")
direder_TB1_seurat <- AddMetaData(direder_TB1_seurat, "Direder5", col.name = "donor")
direder_TB2_seurat <- AddMetaData(direder_TB2_seurat, "Direder6", col.name = "donor")
direder_TB3_seurat <- AddMetaData(direder_TB3_seurat, "Direder7", col.name = "donor")
direder_TB4_seurat <- AddMetaData(direder_TB4_seurat, "Direder8", col.name = "donor")
direder_TC1_seurat <- AddMetaData(direder_TC1_seurat, "Direder9", col.name = "donor")
direder_TC2_seurat <- AddMetaData(direder_TC2_seurat, "Direder10", col.name = "donor")
direder_TH_MG1_seurat <- AddMetaData(direder_TH_MG1_seurat, "Direder11", col.name = "donor")
direder_TH_MG2_seurat <- AddMetaData(direder_TH_MG2_seurat, "Direder12", col.name = "donor")


#(13)_flag-10_sample####
xin1_seurat <- AddMetaData(xin1_seurat, "TET_AB01", col.name = "sample")
xin2_seurat <- AddMetaData(xin2_seurat, "TET_C01", col.name = "sample")
xin3_seurat <- AddMetaData(xin3_seurat, "TET_A01", col.name = "sample")
xin4_seurat <- AddMetaData(xin4_seurat, "TET_AB02", col.name = "sample")
xin5_seurat <- AddMetaData(xin5_seurat, "TET_MNT01", col.name = "sample")
xin6_seurat <- AddMetaData(xin6_seurat, "TET_B01", col.name = "sample")
yasumizu_MG03_seurat <- AddMetaData(yasumizu_MG03_seurat, "TET_AB03", col.name = "sample")
yasumizu_MG21_seurat <- AddMetaData(yasumizu_MG21_seurat, "TET_B02", col.name = "sample")
yasumizu_MG22_seurat <- AddMetaData(yasumizu_MG22_seurat, "TET_B03", col.name = "sample")
yasumizu_MG23_seurat <- AddMetaData(yasumizu_MG23_seurat, "TET_AB04", col.name = "sample")
park_A16_seurat <- AddMetaData(park_A16_seurat, "Thymus01", col.name = "sample")
park_A43_seurat <- AddMetaData(park_A43_seurat, "Thymus02", col.name = "sample") 
park_A45_seurat <- AddMetaData(park_A45_seurat, "Thymus03", col.name = "sample")
park_C34_seurat <- AddMetaData(park_C34_seurat, "Thymus04", col.name = "sample")
park_C40_seurat <- AddMetaData(park_C40_seurat, "Thymus05", col.name = "sample") 
park_C41_seurat <- AddMetaData(park_C41_seurat, "Thymus06", col.name = "sample")
park_F21_seurat <- AddMetaData(park_F21_seurat, "Thymus07", col.name = "sample")
park_F22_seurat <- AddMetaData(park_F22_seurat, "Thymus08", col.name = "sample")
park_F23_seurat <- AddMetaData(park_F23_seurat, "Thymus09", col.name = "sample")
park_F29_seurat <- AddMetaData(park_F29_seurat, "Thymus10", col.name = "sample")
park_F30_seurat <- AddMetaData(park_F30_seurat, "Thymus11", col.name = "sample")
park_F38_seurat <- AddMetaData(park_F38_seurat, "Thymus12", col.name = "sample") 
park_F41_seurat <- AddMetaData(park_F41_seurat, "Thymus13", col.name = "sample")
park_F45_seurat <- AddMetaData(park_F45_seurat, "Thymus14", col.name = "sample")
park_F64_seurat <- AddMetaData(park_F64_seurat, "Thymus15", col.name = "sample")
park_F67_seurat <- AddMetaData(park_F67_seurat, "Thymus16", col.name = "sample") 
park_F74_seurat <- AddMetaData(park_F74_seurat, "Thymus17", col.name = "sample")
park_F83_seurat <- AddMetaData(park_F83_seurat, "Thymus18", col.name = "sample")
park_P1_seurat <- AddMetaData(park_P1_seurat, "Thymus19", col.name = "sample")
park_P2_seurat <- AddMetaData(park_P2_seurat, "Thymus20", col.name = "sample")
park_P3_seurat <- AddMetaData(park_P3_seurat, "Thymus21", col.name = "sample")
park_T03_seurat <- AddMetaData(park_T03_seurat, "Thymus22", col.name = "sample")
park_T06_seurat <- AddMetaData(park_T06_seurat, "Thymus23", col.name = "sample")
park_T07_seurat <- AddMetaData(park_T07_seurat, "Thymus24", col.name = "sample")
bautista_a25_seurat <- AddMetaData(bautista_a25_seurat, "Thymus25", col.name = "sample")
bautista_f19_seurat <- AddMetaData(bautista_f19_seurat, "Thymus26", col.name = "sample")
bautista_f231_seurat <- AddMetaData(bautista_f231_seurat, "Thymus27", col.name = "sample")
bautista_f232_seurat <- AddMetaData(bautista_f232_seurat, "Thymus28", col.name = "sample")
bautista_p6_seurat <- AddMetaData(bautista_p6_seurat, "Thymus29", col.name = "sample")
bautista_p101_seurat <- AddMetaData(bautista_p101_seurat, "Thymus30", col.name = "sample")
bautista_p102_seurat <- AddMetaData(bautista_p102_seurat, "Thymus31", col.name = "sample")
direder_TH1_seurat <- AddMetaData(direder_TH1_seurat, "TTH01", col.name = "sample")
direder_TH2_seurat <- AddMetaData(direder_TH2_seurat, "TTH02", col.name = "sample")
direder_TH3_seurat <- AddMetaData(direder_TH3_seurat, "TTH03", col.name = "sample")
direder_TH4_seurat <- AddMetaData(direder_TH4_seurat, "TTH04", col.name = "sample")
direder_TB1_seurat <- AddMetaData(direder_TB1_seurat, "TET_B04", col.name = "sample")
direder_TB2_seurat <- AddMetaData(direder_TB2_seurat, "TET_B05", col.name = "sample")
direder_TB3_seurat <- AddMetaData(direder_TB3_seurat, "TET_B06", col.name = "sample")
direder_TB4_seurat <- AddMetaData(direder_TB4_seurat, "TET_B07", col.name = "sample")
direder_TC1_seurat <- AddMetaData(direder_TC1_seurat, "TET_C02", col.name = "sample")
direder_TC2_seurat <- AddMetaData(direder_TC2_seurat, "TET_C03", col.name = "sample")
direder_TH_MG1_seurat <- AddMetaData(direder_TH_MG1_seurat, "TTH05", col.name = "sample")
direder_TH_MG2_seurat <- AddMetaData(direder_TH_MG2_seurat, "TTH06", col.name = "sample")

#(14)_flag-11_maintissue####
xin1_seurat <- AddMetaData(xin1_seurat, "TET", col.name = "maintissue")
xin2_seurat <- AddMetaData(xin2_seurat, "TET", col.name = "maintissue")
xin3_seurat <- AddMetaData(xin3_seurat, "TET", col.name = "maintissue")
xin4_seurat <- AddMetaData(xin4_seurat, "TET", col.name = "maintissue")
xin5_seurat <- AddMetaData(xin5_seurat, "TET", col.name = "maintissue")
xin6_seurat <- AddMetaData(xin6_seurat, "TET", col.name = "maintissue")
yasumizu_MG03_seurat <- AddMetaData(yasumizu_MG03_seurat, "TET", col.name = "maintissue")
yasumizu_MG21_seurat <- AddMetaData(yasumizu_MG21_seurat, "TET", col.name = "maintissue")
yasumizu_MG22_seurat <- AddMetaData(yasumizu_MG22_seurat, "TET", col.name = "maintissue")
yasumizu_MG23_seurat <- AddMetaData(yasumizu_MG23_seurat, "TET", col.name = "maintissue")
park_A16_seurat <- AddMetaData(park_A16_seurat, "Thymus", col.name = "maintissue")
park_A43_seurat <- AddMetaData(park_A43_seurat, "Thymus", col.name = "maintissue") 
park_A45_seurat <- AddMetaData(park_A45_seurat, "Thymus", col.name = "maintissue")
park_C34_seurat <- AddMetaData(park_C34_seurat, "Thymus", col.name = "maintissue")
park_C40_seurat <- AddMetaData(park_C40_seurat, "Thymus", col.name = "maintissue") 
park_C41_seurat <- AddMetaData(park_C41_seurat, "Thymus", col.name = "maintissue")
park_F21_seurat <- AddMetaData(park_F21_seurat, "Thymus", col.name = "maintissue")
park_F22_seurat <- AddMetaData(park_F22_seurat, "Thymus", col.name = "maintissue")
park_F23_seurat <- AddMetaData(park_F23_seurat, "Thymus", col.name = "maintissue")
park_F29_seurat <- AddMetaData(park_F29_seurat, "Thymus", col.name = "maintissue")
park_F30_seurat <- AddMetaData(park_F30_seurat, "Thymus", col.name = "maintissue")
park_F38_seurat <- AddMetaData(park_F38_seurat, "Thymus", col.name = "maintissue") 
park_F41_seurat <- AddMetaData(park_F41_seurat, "Thymus", col.name = "maintissue")
park_F45_seurat <- AddMetaData(park_F45_seurat, "Thymus", col.name = "maintissue")
park_F64_seurat <- AddMetaData(park_F64_seurat, "Thymus", col.name = "maintissue")
park_F67_seurat <- AddMetaData(park_F67_seurat, "Thymus", col.name = "maintissue") 
park_F74_seurat <- AddMetaData(park_F74_seurat, "Thymus", col.name = "maintissue")
park_F83_seurat <- AddMetaData(park_F83_seurat, "Thymus", col.name = "maintissue")
park_P1_seurat <- AddMetaData(park_P1_seurat, "Thymus", col.name = "maintissue")
park_P2_seurat <- AddMetaData(park_P2_seurat, "Thymus", col.name = "maintissue")
park_P3_seurat <- AddMetaData(park_P3_seurat, "Thymus", col.name = "maintissue")
park_T03_seurat <- AddMetaData(park_T03_seurat, "Thymus", col.name = "maintissue")
park_T06_seurat <- AddMetaData(park_T06_seurat, "Thymus", col.name = "maintissue")
park_T07_seurat <- AddMetaData(park_T07_seurat, "Thymus", col.name = "maintissue")
bautista_a25_seurat <- AddMetaData(bautista_a25_seurat, "Thymus", col.name = "maintissue")
bautista_f19_seurat <- AddMetaData(bautista_f19_seurat, "Thymus", col.name = "maintissue")
bautista_f231_seurat <- AddMetaData(bautista_f231_seurat, "Thymus", col.name = "maintissue")
bautista_f232_seurat <- AddMetaData(bautista_f232_seurat, "Thymus", col.name = "maintissue")
bautista_p6_seurat <- AddMetaData(bautista_p6_seurat, "Thymus", col.name = "maintissue")
bautista_p101_seurat <- AddMetaData(bautista_p101_seurat, "Thymus", col.name = "maintissue")
bautista_p102_seurat <- AddMetaData(bautista_p102_seurat, "Thymus", col.name = "maintissue")
direder_TH1_seurat <- AddMetaData(direder_TH1_seurat, "TTH", col.name = "maintissue")
direder_TH2_seurat <- AddMetaData(direder_TH2_seurat, "TTH", col.name = "maintissue")
direder_TH3_seurat <- AddMetaData(direder_TH3_seurat, "TTH", col.name = "maintissue")
direder_TH4_seurat <- AddMetaData(direder_TH4_seurat, "TTH", col.name = "maintissue")
direder_TB1_seurat <- AddMetaData(direder_TB1_seurat, "TET", col.name = "maintissue")
direder_TB2_seurat <- AddMetaData(direder_TB2_seurat, "TET", col.name = "maintissue")
direder_TB3_seurat <- AddMetaData(direder_TB3_seurat, "TET", col.name = "maintissue")
direder_TB4_seurat <- AddMetaData(direder_TB4_seurat, "TET", col.name = "maintissue")
direder_TC1_seurat <- AddMetaData(direder_TC1_seurat, "TET", col.name = "maintissue")
direder_TC2_seurat <- AddMetaData(direder_TC2_seurat, "TET", col.name = "maintissue")
direder_TH_MG1_seurat <- AddMetaData(direder_TH_MG1_seurat, "TTH", col.name = "maintissue")
direder_TH_MG2_seurat <- AddMetaData(direder_TH_MG2_seurat, "TTH", col.name = "maintissue")
#(15)_Quality ctrl####
#StartLoop_2
rawdat_4 <- list(xin1_seurat, xin2_seurat, xin3_seurat, xin4_seurat, xin5_seurat, xin6_seurat, 
                 yasumizu_MG03_seurat, yasumizu_MG21_seurat, yasumizu_MG22_seurat, yasumizu_MG23_seurat, 
                 park_A16_seurat, park_A43_seurat, park_A45_seurat, park_C34_seurat, park_C40_seurat, park_C41_seurat, park_F21_seurat, park_F22_seurat, park_F23_seurat, park_F29_seurat, park_F30_seurat, park_F38_seurat, 
                 park_F41_seurat, park_F45_seurat, park_F64_seurat, park_F67_seurat, park_F74_seurat, park_F83_seurat, park_P1_seurat, park_P2_seurat, park_P3_seurat, park_T03_seurat, park_T06_seurat, park_T07_seurat, 
                 bautista_a25_seurat, bautista_f19_seurat, bautista_f231_seurat, bautista_f232_seurat, bautista_p6_seurat, bautista_p101_seurat, bautista_p102_seurat, 
                 direder_TH1_seurat, direder_TH2_seurat, direder_TH3_seurat, direder_TH4_seurat, direder_TB1_seurat, direder_TB2_seurat, direder_TB3_seurat, direder_TB4_seurat, direder_TC1_seurat, direder_TC2_seurat, direder_TH_MG1_seurat, direder_TH_MG2_seurat)

rawdat_5 <- lapply(rawdat_4, function(x){
  x[["percent.mt"]] <- PercentageFeatureSet(x, pattern = "^MT-")
  x[["percent.ERY"]] <- PercentageFeatureSet(x, pattern = "HBB")
  x <- subset(x, subset = percent.ERY < 5 )
})

#QC_Graph#
pdf("Graphs/Quality control/QC.pdf")
VlnPlot(rawdat_5[[1]], features = c("nFeature_RNA", "nCount_RNA","percent.ERY","percent.mt"), ncol = 4) & ggtitle("percent.mt_1")
FeatureScatter(rawdat_5[[1]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + ggtitle("1")
FeatureScatter(rawdat_5[[1]], feature1 = "nCount_RNA", feature2 = "percent.mt") + ggtitle("1")
VlnPlot(rawdat_5[[2]], features = c("nFeature_RNA", "nCount_RNA","percent.ERY","percent.mt"), ncol = 3) + ggtitle("percent.mt_2")
FeatureScatter(rawdat_5[[2]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + ggtitle("2")
FeatureScatter(rawdat_5[[2]], feature1 = "nCount_RNA", feature2 = "percent.mt") + ggtitle("2")
VlnPlot(rawdat_5[[3]], features = c("nFeature_RNA", "nCount_RNA","percent.ERY","percent.mt"), ncol = 3) + ggtitle("percent.mt_3")
FeatureScatter(rawdat_5[[3]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + ggtitle("3")
FeatureScatter(rawdat_5[[3]], feature1 = "nCount_RNA", feature2 = "percent.mt") + ggtitle("3")
VlnPlot(rawdat_5[[4]], features = c("nFeature_RNA", "nCount_RNA","percent.ERY","percent.mt"), ncol = 3) + ggtitle("percent.mt_4")
FeatureScatter(rawdat_5[[4]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + ggtitle("4")
FeatureScatter(rawdat_5[[4]], feature1 = "nCount_RNA", feature2 = "percent.mt") + ggtitle("4")
VlnPlot(rawdat_5[[5]], features = c("nFeature_RNA", "nCount_RNA","percent.ERY","percent.mt"), ncol = 3) + ggtitle("percent.mt_5")
FeatureScatter(rawdat_5[[5]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + ggtitle("5")
FeatureScatter(rawdat_5[[5]], feature1 = "nCount_RNA", feature2 = "percent.mt") + ggtitle("5")
VlnPlot(rawdat_5[[6]], features = c("nFeature_RNA", "nCount_RNA","percent.ERY","percent.mt"), ncol = 3) + ggtitle("percent.mt_6")
FeatureScatter(rawdat_5[[6]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + ggtitle("6")
FeatureScatter(rawdat_5[[6]], feature1 = "nCount_RNA", feature2 = "percent.mt") + ggtitle("6")
VlnPlot(rawdat_5[[7]], features = c("nFeature_RNA", "nCount_RNA","percent.ERY","percent.mt"), ncol = 3) + ggtitle("percent.mt_7")
FeatureScatter(rawdat_5[[7]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + ggtitle("7")
FeatureScatter(rawdat_5[[7]], feature1 = "nCount_RNA", feature2 = "percent.mt") + ggtitle("7")
VlnPlot(rawdat_5[[8]], features = c("nFeature_RNA", "nCount_RNA","percent.ERY","percent.mt"), ncol = 3) + ggtitle("percent.mt_8")
FeatureScatter(rawdat_5[[8]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + ggtitle("8")
FeatureScatter(rawdat_5[[8]], feature1 = "nCount_RNA", feature2 = "percent.mt") + ggtitle("8")
VlnPlot(rawdat_5[[9]], features = c("nFeature_RNA", "nCount_RNA","percent.ERY","percent.mt"), ncol = 3) + ggtitle("percent.mt_9")
FeatureScatter(rawdat_5[[9]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + ggtitle("9")
FeatureScatter(rawdat_5[[9]], feature1 = "nCount_RNA", feature2 = "percent.mt") + ggtitle("9")
VlnPlot(rawdat_5[[10]], features = c("nFeature_RNA", "nCount_RNA","percent.ERY","percent.mt"), ncol = 3) + ggtitle("percent.mt_10")
FeatureScatter(rawdat_5[[10]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + ggtitle("10")
FeatureScatter(rawdat_5[[10]], feature1 = "nCount_RNA", feature2 = "percent.mt") + ggtitle("10")
VlnPlot(rawdat_5[[11]], features = c("nFeature_RNA", "nCount_RNA","percent.ERY","percent.mt"), ncol = 3) + ggtitle("percent.mt_11")
FeatureScatter(rawdat_5[[11]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + ggtitle("11")
FeatureScatter(rawdat_5[[11]], feature1 = "nCount_RNA", feature2 = "percent.mt") + ggtitle("11")
VlnPlot(rawdat_5[[12]], features = c("nFeature_RNA", "nCount_RNA","percent.ERY","percent.mt"), ncol = 3) + ggtitle("percent.mt_12")
FeatureScatter(rawdat_5[[12]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + ggtitle("12")
FeatureScatter(rawdat_5[[12]], feature1 = "nCount_RNA", feature2 = "percent.mt") + ggtitle("12")
VlnPlot(rawdat_5[[13]], features = c("nFeature_RNA", "nCount_RNA","percent.ERY","percent.mt"), ncol = 3) + ggtitle("percent.mt_13")
FeatureScatter(rawdat_5[[13]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + ggtitle("13")
FeatureScatter(rawdat_5[[13]], feature1 = "nCount_RNA", feature2 = "percent.mt") + ggtitle("13")
VlnPlot(rawdat_5[[14]], features = c("nFeature_RNA", "nCount_RNA","percent.ERY","percent.mt"), ncol = 3) + ggtitle("percent.mt_14")
FeatureScatter(rawdat_5[[14]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + ggtitle("14")
FeatureScatter(rawdat_5[[14]], feature1 = "nCount_RNA", feature2 = "percent.mt") + ggtitle("14")
VlnPlot(rawdat_5[[15]], features = c("nFeature_RNA", "nCount_RNA","percent.ERY","percent.mt"), ncol = 3) + ggtitle("percent.mt_15")
FeatureScatter(rawdat_5[[15]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + ggtitle("15")
FeatureScatter(rawdat_5[[15]], feature1 = "nCount_RNA", feature2 = "percent.mt") + ggtitle("15")
VlnPlot(rawdat_5[[16]], features = c("nFeature_RNA", "nCount_RNA","percent.ERY","percent.mt"), ncol = 3) + ggtitle("percent.mt_16")
FeatureScatter(rawdat_5[[16]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + ggtitle("16")
FeatureScatter(rawdat_5[[16]], feature1 = "nCount_RNA", feature2 = "percent.mt") + ggtitle("16")
VlnPlot(rawdat_5[[17]], features = c("nFeature_RNA", "nCount_RNA","percent.ERY","percent.mt"), ncol = 3) + ggtitle("percent.mt_17")
FeatureScatter(rawdat_5[[17]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + ggtitle("17")
FeatureScatter(rawdat_5[[17]], feature1 = "nCount_RNA", feature2 = "percent.mt") + ggtitle("17")
VlnPlot(rawdat_5[[18]], features = c("nFeature_RNA", "nCount_RNA","percent.ERY","percent.mt"), ncol = 3) + ggtitle("percent.mt_18")
FeatureScatter(rawdat_5[[18]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + ggtitle("18")
FeatureScatter(rawdat_5[[18]], feature1 = "nCount_RNA", feature2 = "percent.mt") + ggtitle("18")
VlnPlot(rawdat_5[[19]], features = c("nFeature_RNA", "nCount_RNA","percent.ERY","percent.mt"), ncol = 3) + ggtitle("percent.mt_19")
FeatureScatter(rawdat_5[[19]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + ggtitle("19")
FeatureScatter(rawdat_5[[19]], feature1 = "nCount_RNA", feature2 = "percent.mt") + ggtitle("19")
VlnPlot(rawdat_5[[20]], features = c("nFeature_RNA", "nCount_RNA","percent.ERY","percent.mt"), ncol = 3) + ggtitle("percent.mt_20")
FeatureScatter(rawdat_5[[20]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + ggtitle("20")
FeatureScatter(rawdat_5[[20]], feature1 = "nCount_RNA", feature2 = "percent.mt") + ggtitle("20")
VlnPlot(rawdat_5[[21]], features = c("nFeature_RNA", "nCount_RNA","percent.ERY","percent.mt"), ncol = 3) + ggtitle("percent.mt_21")
FeatureScatter(rawdat_5[[21]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + ggtitle("21")
FeatureScatter(rawdat_5[[21]], feature1 = "nCount_RNA", feature2 = "percent.mt") + ggtitle("21")
VlnPlot(rawdat_5[[22]], features = c("nFeature_RNA", "nCount_RNA","percent.ERY","percent.mt"), ncol = 3) + ggtitle("percent.mt_22")
FeatureScatter(rawdat_5[[22]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + ggtitle("22")
FeatureScatter(rawdat_5[[22]], feature1 = "nCount_RNA", feature2 = "percent.mt") + ggtitle("22")
VlnPlot(rawdat_5[[23]], features = c("nFeature_RNA", "nCount_RNA","percent.ERY","percent.mt"), ncol = 3) + ggtitle("percent.mt_23")
FeatureScatter(rawdat_5[[23]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + ggtitle("23")
FeatureScatter(rawdat_5[[23]], feature1 = "nCount_RNA", feature2 = "percent.mt") + ggtitle("23")
VlnPlot(rawdat_5[[24]], features = c("nFeature_RNA", "nCount_RNA","percent.ERY","percent.mt"), ncol = 3) + ggtitle("percent.mt_24")
FeatureScatter(rawdat_5[[24]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + ggtitle("24")
FeatureScatter(rawdat_5[[24]], feature1 = "nCount_RNA", feature2 = "percent.mt") + ggtitle("24")
VlnPlot(rawdat_5[[25]], features = c("nFeature_RNA", "nCount_RNA","percent.ERY","percent.mt"), ncol = 3) + ggtitle("percent.mt_25")
FeatureScatter(rawdat_5[[25]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + ggtitle("25")
FeatureScatter(rawdat_5[[25]], feature1 = "nCount_RNA", feature2 = "percent.mt") + ggtitle("25")
VlnPlot(rawdat_5[[26]], features = c("nFeature_RNA", "nCount_RNA","percent.ERY","percent.mt"), ncol = 3) + ggtitle("percent.mt_26")
FeatureScatter(rawdat_5[[26]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + ggtitle("26")
FeatureScatter(rawdat_5[[26]], feature1 = "nCount_RNA", feature2 = "percent.mt") + ggtitle("26")
VlnPlot(rawdat_5[[27]], features = c("nFeature_RNA", "nCount_RNA","percent.ERY","percent.mt"), ncol = 3) + ggtitle("percent.mt_27")
FeatureScatter(rawdat_5[[27]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + ggtitle("27")
FeatureScatter(rawdat_5[[27]], feature1 = "nCount_RNA", feature2 = "percent.mt") + ggtitle("27")
VlnPlot(rawdat_5[[28]], features = c("nFeature_RNA", "nCount_RNA","percent.ERY","percent.mt"), ncol = 3) + ggtitle("percent.mt_28")
FeatureScatter(rawdat_5[[28]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + ggtitle("28")
FeatureScatter(rawdat_5[[28]], feature1 = "nCount_RNA", feature2 = "percent.mt") + ggtitle("28")
VlnPlot(rawdat_5[[29]], features = c("nFeature_RNA", "nCount_RNA","percent.ERY","percent.mt"), ncol = 3) + ggtitle("percent.mt_29")
FeatureScatter(rawdat_5[[29]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + ggtitle("29")
FeatureScatter(rawdat_5[[29]], feature1 = "nCount_RNA", feature2 = "percent.mt") + ggtitle("29")
VlnPlot(rawdat_5[[30]], features = c("nFeature_RNA", "nCount_RNA","percent.ERY","percent.mt"), ncol = 3) + ggtitle("percent.mt_30")
FeatureScatter(rawdat_5[[30]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + ggtitle("30")
FeatureScatter(rawdat_5[[30]], feature1 = "nCount_RNA", feature2 = "percent.mt") + ggtitle("30")
VlnPlot(rawdat_5[[31]], features = c("nFeature_RNA", "nCount_RNA","percent.ERY","percent.mt"), ncol = 3) + ggtitle("percent.mt_31")
FeatureScatter(rawdat_5[[31]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + ggtitle("31")
FeatureScatter(rawdat_5[[31]], feature1 = "nCount_RNA", feature2 = "percent.mt") + ggtitle("31")
VlnPlot(rawdat_5[[32]], features = c("nFeature_RNA", "nCount_RNA","percent.ERY","percent.mt"), ncol = 3) + ggtitle("percent.mt_32")
FeatureScatter(rawdat_5[[32]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + ggtitle("32")
FeatureScatter(rawdat_5[[32]], feature1 = "nCount_RNA", feature2 = "percent.mt") + ggtitle("32")
VlnPlot(rawdat_5[[33]], features = c("nFeature_RNA", "nCount_RNA","percent.ERY","percent.mt"), ncol = 3) + ggtitle("percent.mt_33")
FeatureScatter(rawdat_5[[33]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + ggtitle("33")
FeatureScatter(rawdat_5[[33]], feature1 = "nCount_RNA", feature2 = "percent.mt") + ggtitle("33")
VlnPlot(rawdat_5[[34]], features = c("nFeature_RNA", "nCount_RNA","percent.ERY","percent.mt"), ncol = 3) + ggtitle("percent.mt_34")
FeatureScatter(rawdat_5[[34]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + ggtitle("34")
FeatureScatter(rawdat_5[[34]], feature1 = "nCount_RNA", feature2 = "percent.mt") + ggtitle("34")
VlnPlot(rawdat_5[[35]], features = c("nFeature_RNA", "nCount_RNA","percent.ERY","percent.mt"), ncol = 3) + ggtitle("percent.mt_35")
FeatureScatter(rawdat_5[[35]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + ggtitle("35")
FeatureScatter(rawdat_5[[35]], feature1 = "nCount_RNA", feature2 = "percent.mt") + ggtitle("35")
VlnPlot(rawdat_5[[36]], features = c("nFeature_RNA", "nCount_RNA","percent.ERY","percent.mt"), ncol = 3) + ggtitle("percent.mt_36")
FeatureScatter(rawdat_5[[36]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + ggtitle("36")
FeatureScatter(rawdat_5[[36]], feature1 = "nCount_RNA", feature2 = "percent.mt") + ggtitle("36")
VlnPlot(rawdat_5[[37]], features = c("nFeature_RNA", "nCount_RNA","percent.ERY","percent.mt"), ncol = 3) + ggtitle("percent.mt_37")
FeatureScatter(rawdat_5[[37]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + ggtitle("37")
FeatureScatter(rawdat_5[[37]], feature1 = "nCount_RNA", feature2 = "percent.mt") + ggtitle("37")
VlnPlot(rawdat_5[[38]], features = c("nFeature_RNA", "nCount_RNA","percent.ERY","percent.mt"), ncol = 3) + ggtitle("percent.mt_38")
FeatureScatter(rawdat_5[[38]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + ggtitle("38")
FeatureScatter(rawdat_5[[38]], feature1 = "nCount_RNA", feature2 = "percent.mt") + ggtitle("38")
VlnPlot(rawdat_5[[39]], features = c("nFeature_RNA", "nCount_RNA","percent.ERY","percent.mt"), ncol = 3) + ggtitle("percent.mt_39")
FeatureScatter(rawdat_5[[39]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + ggtitle("39")
FeatureScatter(rawdat_5[[39]], feature1 = "nCount_RNA", feature2 = "percent.mt") + ggtitle("39")
VlnPlot(rawdat_5[[40]], features = c("nFeature_RNA", "nCount_RNA","percent.ERY","percent.mt"), ncol = 3) + ggtitle("percent.mt_40")
FeatureScatter(rawdat_5[[40]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + ggtitle("40")
FeatureScatter(rawdat_5[[40]], feature1 = "nCount_RNA", feature2 = "percent.mt") + ggtitle("40")
VlnPlot(rawdat_5[[41]], features = c("nFeature_RNA", "nCount_RNA","percent.ERY","percent.mt"), ncol = 3) + ggtitle("percent.mt_41")
FeatureScatter(rawdat_5[[41]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + ggtitle("41")
FeatureScatter(rawdat_5[[41]], feature1 = "nCount_RNA", feature2 = "percent.mt") + ggtitle("41")
VlnPlot(rawdat_5[[42]], features = c("nFeature_RNA", "nCount_RNA","percent.ERY","percent.mt"), ncol = 3) + ggtitle("percent.mt_42")
FeatureScatter(rawdat_5[[42]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + ggtitle("42")
FeatureScatter(rawdat_5[[42]], feature1 = "nCount_RNA", feature2 = "percent.mt") + ggtitle("42")
VlnPlot(rawdat_5[[43]], features = c("nFeature_RNA", "nCount_RNA","percent.ERY","percent.mt"), ncol = 3) + ggtitle("percent.mt_43")
FeatureScatter(rawdat_5[[43]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + ggtitle("43")
FeatureScatter(rawdat_5[[43]], feature1 = "nCount_RNA", feature2 = "percent.mt") + ggtitle("43")
VlnPlot(rawdat_5[[44]], features = c("nFeature_RNA", "nCount_RNA","percent.ERY","percent.mt"), ncol = 3) + ggtitle("percent.mt_44")
FeatureScatter(rawdat_5[[44]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + ggtitle("44")
FeatureScatter(rawdat_5[[44]], feature1 = "nCount_RNA", feature2 = "percent.mt") + ggtitle("44")
VlnPlot(rawdat_5[[45]], features = c("nFeature_RNA", "nCount_RNA","percent.ERY","percent.mt"), ncol = 3) + ggtitle("percent.mt_45")
FeatureScatter(rawdat_5[[45]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + ggtitle("45")
FeatureScatter(rawdat_5[[45]], feature1 = "nCount_RNA", feature2 = "percent.mt") + ggtitle("45")
VlnPlot(rawdat_5[[46]], features = c("nFeature_RNA", "nCount_RNA","percent.ERY","percent.mt"), ncol = 3) + ggtitle("percent.mt_46")
FeatureScatter(rawdat_5[[46]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + ggtitle("46")
FeatureScatter(rawdat_5[[46]], feature1 = "nCount_RNA", feature2 = "percent.mt") + ggtitle("46")
VlnPlot(rawdat_5[[47]], features = c("nFeature_RNA", "nCount_RNA","percent.ERY","percent.mt"), ncol = 3) + ggtitle("percent.mt_47")
FeatureScatter(rawdat_5[[47]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + ggtitle("47")
FeatureScatter(rawdat_5[[47]], feature1 = "nCount_RNA", feature2 = "percent.mt") + ggtitle("47")
VlnPlot(rawdat_5[[48]], features = c("nFeature_RNA", "nCount_RNA","percent.ERY","percent.mt"), ncol = 3) + ggtitle("percent.mt_48")
FeatureScatter(rawdat_5[[48]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + ggtitle("48")
FeatureScatter(rawdat_5[[48]], feature1 = "nCount_RNA", feature2 = "percent.mt") + ggtitle("48")
VlnPlot(rawdat_5[[49]], features = c("nFeature_RNA", "nCount_RNA","percent.ERY","percent.mt"), ncol = 3) + ggtitle("percent.mt_49")
FeatureScatter(rawdat_5[[49]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + ggtitle("49")
FeatureScatter(rawdat_5[[49]], feature1 = "nCount_RNA", feature2 = "percent.mt") + ggtitle("49")
VlnPlot(rawdat_5[[50]], features = c("nFeature_RNA", "nCount_RNA","percent.ERY","percent.mt"), ncol = 3) + ggtitle("percent.mt_50")
FeatureScatter(rawdat_5[[50]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + ggtitle("50")
FeatureScatter(rawdat_5[[50]], feature1 = "nCount_RNA", feature2 = "percent.mt") + ggtitle("50")
VlnPlot(rawdat_5[[51]], features = c("nFeature_RNA", "nCount_RNA","percent.ERY","percent.mt"), ncol = 3) + ggtitle("percent.mt_51")
FeatureScatter(rawdat_5[[51]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + ggtitle("51")
FeatureScatter(rawdat_5[[51]], feature1 = "nCount_RNA", feature2 = "percent.mt") + ggtitle("51")
VlnPlot(rawdat_5[[52]], features = c("nFeature_RNA", "nCount_RNA","percent.ERY","percent.mt"), ncol = 3) + ggtitle("percent.mt_52")
FeatureScatter(rawdat_5[[52]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + ggtitle("52")
FeatureScatter(rawdat_5[[52]], feature1 = "nCount_RNA", feature2 = "percent.mt") + ggtitle("52")
VlnPlot(rawdat_5[[53]], features = c("nFeature_RNA", "nCount_RNA","percent.ERY","percent.mt"), ncol = 3) + ggtitle("percent.mt_53")
FeatureScatter(rawdat_5[[53]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + ggtitle("53")
FeatureScatter(rawdat_5[[53]], feature1 = "nCount_RNA", feature2 = "percent.mt") + ggtitle("53")
dev.off()


#remove surplus objects
rm(xin1_seurat, xin2_seurat, xin3_seurat, xin4_seurat, xin5_seurat, xin6_seurat, 
   yasumizu_MG03_seurat, yasumizu_MG21_seurat, yasumizu_MG22_seurat, yasumizu_MG23_seurat, 
   park_A16_seurat, park_A43_seurat, park_A45_seurat, park_C34_seurat, park_C40_seurat, park_C41_seurat, park_F21_seurat, park_F22_seurat, park_F23_seurat, park_F29_seurat, park_F30_seurat, park_F38_seurat, 
   park_F41_seurat, park_F45_seurat, park_F64_seurat, park_F67_seurat, park_F74_seurat, park_F83_seurat, park_P1_seurat, park_P2_seurat, park_P3_seurat, park_T03_seurat, park_T06_seurat, park_T07_seurat, 
   bautista_a25_seurat, bautista_f19_seurat, bautista_f231_seurat, bautista_f232_seurat, bautista_p6_seurat, bautista_p101_seurat, bautista_p102_seurat, 
   direder_TH1_seurat, direder_TH2_seurat, direder_TH3_seurat, direder_TH4_seurat, direder_TB1_seurat, direder_TB2_seurat, direder_TB3_seurat, direder_TB4_seurat, direder_TC1_seurat, direder_TC2_seurat, direder_TH_MG1_seurat, direder_TH_MG2_seurat, rawdat_4)

#save rawdat_5
#save(rawdat_5,file="rawdat_5.RData")

rawdat_6 <- lapply(rawdat_5, function(x){
  x <- subset(x, subset = nFeature_RNA > 200 & nFeature_RNA < 3000 & percent.mt < 7)
})

#(16)_separate list to merge####
xin1_seurat <- rawdat_6[[1]] 
xin2_seurat <- rawdat_6[[2]] 
xin3_seurat <- rawdat_6[[3]] 
xin4_seurat <- rawdat_6[[4]] 
xin5_seurat <- rawdat_6[[5]] 
xin6_seurat <- rawdat_6[[6]] 
yasumizu_MG03_seurat <- rawdat_6[[7]] 
yasumizu_MG21_seurat <- rawdat_6[[8]] 
yasumizu_MG22_seurat <- rawdat_6[[9]] 
yasumizu_MG23_seurat <- rawdat_6[[10]] 
park_A16_seurat <- rawdat_6[[11]] 
park_A43_seurat <- rawdat_6[[12]] 
park_A45_seurat <- rawdat_6[[13]] 
park_C34_seurat <- rawdat_6[[14]] 
park_C40_seurat <- rawdat_6[[15]] 
park_C41_seurat <- rawdat_6[[16]] 
park_F21_seurat <- rawdat_6[[17]] 
park_F22_seurat <- rawdat_6[[18]] 
park_F23_seurat <- rawdat_6[[19]] 
park_F29_seurat <- rawdat_6[[20]] 
park_F30_seurat <- rawdat_6[[21]] 
park_F38_seurat <- rawdat_6[[22]] 
park_F41_seurat <- rawdat_6[[23]] 
park_F45_seurat <- rawdat_6[[24]] 
park_F64_seurat <- rawdat_6[[25]] 
park_F67_seurat <- rawdat_6[[26]] 
park_F74_seurat <- rawdat_6[[27]] 
park_F83_seurat <- rawdat_6[[28]] 
park_P1_seurat <- rawdat_6[[29]] 
park_P2_seurat <- rawdat_6[[30]] 
park_P3_seurat <- rawdat_6[[31]] 
park_T03_seurat <- rawdat_6[[32]] 
park_T06_seurat <- rawdat_6[[33]] 
park_T07_seurat <- rawdat_6[[34]] 
bautista_a25_seurat <- rawdat_6[[35]] 
bautista_f19_seurat <- rawdat_6[[36]] 
bautista_f231_seurat <- rawdat_6[[37]] 
bautista_f232_seurat <- rawdat_6[[38]] 
bautista_p6_seurat <- rawdat_6[[39]] 
bautista_p101_seurat <- rawdat_6[[40]] 
bautista_p102_seurat <- rawdat_6[[41]] 
direder_TH1_seurat <- rawdat_6[[42]] 
direder_TH2_seurat <- rawdat_6[[43]] 
direder_TH3_seurat <- rawdat_6[[44]] 
direder_TH4_seurat <- rawdat_6[[45]] 
direder_TB1_seurat <- rawdat_6[[46]] 
direder_TB2_seurat <- rawdat_6[[47]] 
direder_TB3_seurat <- rawdat_6[[48]] 
direder_TB4_seurat <- rawdat_6[[49]] 
direder_TC1_seurat <- rawdat_6[[50]] 
direder_TC2_seurat <- rawdat_6[[51]] 
direder_TH_MG1_seurat <- rawdat_6[[52]] 
direder_TH_MG2_seurat <- rawdat_6[[53]] 

rawdat_7 <- merge( xin1_seurat, y= c(xin2_seurat, xin3_seurat, xin4_seurat, xin5_seurat, xin6_seurat, 
                                     yasumizu_MG03_seurat, yasumizu_MG21_seurat, yasumizu_MG22_seurat, yasumizu_MG23_seurat, 
                                     park_A16_seurat, park_A43_seurat, park_A45_seurat, park_C34_seurat, park_C40_seurat, park_C41_seurat, park_F21_seurat, park_F22_seurat, park_F23_seurat, park_F29_seurat, park_F30_seurat, park_F38_seurat, 
                                     park_F41_seurat, park_F45_seurat, park_F64_seurat, park_F67_seurat, park_F74_seurat, park_F83_seurat, park_P1_seurat, park_P2_seurat, park_P3_seurat, park_T03_seurat, park_T06_seurat, park_T07_seurat, 
                                     bautista_a25_seurat, bautista_f19_seurat, bautista_f231_seurat, bautista_f232_seurat, bautista_p6_seurat, bautista_p101_seurat, bautista_p102_seurat, 
                                     direder_TH1_seurat, direder_TH2_seurat, direder_TH3_seurat, direder_TH4_seurat, direder_TB1_seurat, direder_TB2_seurat, direder_TB3_seurat, direder_TB4_seurat, direder_TC1_seurat, direder_TC2_seurat, direder_TH_MG1_seurat, direder_TH_MG2_seurat), add.cell.ids = c("Xin1","Xin2","Xin3","Xin4","Xin5","Xin6","Yasumizu1","Yasumizu2","Yasumizu3","Yasumizu4","Park1","Park2","Park3","Park4","Park5","Park6","Park7","Park8","Park9","Park10","Park11","Park12","Park13","Park14","Park15","Park16","Park17","Park18","Park19","Park20","Park21","Park22","Park23","Park24","Bautista1","Bautista2","Bautista3","Bautista4","Bautista5","Bautista6","Bautista7","Direder1","Direder2","Direder3","Direder4","Direder5","Direder6","Direder7","Direder8","Direder9","Direder10","Direder11","Direder12"), project = "TET")
rawdat_7[["RNA"]] <- JoinLayers(rawdat_7[["RNA"]])

#remove surplus objects
rm(rawdat_5,rawdat_6, xin1_seurat, xin2_seurat, xin3_seurat, xin4_seurat, xin5_seurat, xin6_seurat,yasumizu_MG03_seurat, yasumizu_MG21_seurat, yasumizu_MG22_seurat, yasumizu_MG23_seurat, 
   park_A16_seurat, park_A43_seurat, park_A45_seurat, park_C34_seurat, park_C40_seurat, park_C41_seurat, park_F21_seurat, park_F22_seurat, park_F23_seurat, park_F29_seurat, park_F30_seurat, park_F38_seurat, 
   park_F41_seurat, park_F45_seurat, park_F64_seurat, park_F67_seurat, park_F74_seurat, park_F83_seurat, park_P1_seurat, park_P2_seurat, park_P3_seurat, park_T03_seurat, park_T06_seurat, park_T07_seurat, 
   bautista_a25_seurat, bautista_f19_seurat, bautista_f231_seurat, bautista_f232_seurat, bautista_p6_seurat, bautista_p101_seurat, bautista_p102_seurat, 
   direder_TH1_seurat, direder_TH2_seurat, direder_TH3_seurat, direder_TH4_seurat, direder_TB1_seurat, direder_TB2_seurat, direder_TB3_seurat, direder_TB4_seurat, direder_TC1_seurat, direder_TC2_seurat, direder_TH_MG1_seurat, direder_TH_MG2_seurat, rawdat_4)


#(17)_integrate all objects ####
rawdat_7[["RNA"]] <- split(rawdat_7[["RNA"]], f = rawdat_7$donor)
rawdat_8 <- NormalizeData(rawdat_7)
rawdat_8 <- FindVariableFeatures(rawdat_8, verbose =F)
rawdat_8 <- SketchData(rawdat_8, ncells = 5000, method = "LeverageScore", sketched.assay = "sketch")
DefaultAssay(rawdat_8) <- "sketch"
rawdat_8 <- FindVariableFeatures(rawdat_8, verbose =T)
rawdat_8 <- ScaleData(rawdat_8, verbose =T)
rawdat_8 <- RunPCA (rawdat_8, verbose =T)
#save rawdat_8
#save(rawdat_8,file="rawdat_8.RData")

#RPCa
rawdat_9 <- IntegrateLayers(rawdat_8, method = RPCAIntegration,  orig.reduction = "pca",
                            new.reduction = "integrated.rpca",
                            dims = 1:30, k.anchor = 20, reference = which(Layers(rawdat_8, search = "data") %in% c("data.Park1")), verbose = T)
# cluster the integrated data
pdf("Graphs/Main/Elbowplot.pdf")
ElbowPlot(rawdat_9, ndims=30)
dev.off()
rawdat_9 <- FindNeighbors(rawdat_9, reduction = "integrated.rpca", dims = 1:18)
#DefaultAssay(rawdat_9)<- "integrated"

rawdat_9 <- FindClusters(rawdat_9, resolution = 0.6)
plot_integrated_clusters(rawdat_9) 
rawdat_9 <- RunUMAP(rawdat_9, reduction = "integrated.rpca", dims = 1:18, return.model = T, verbose = F)

pdf("Graphs/Quality control/UMAP_test_rpca.pdf")
UMAPPlot(rawdat_9) + ggtitle("Basic-UMAP")
dev.off()

#visualize cellcluster
plot.s1 <- DimPlot(rawdat_9, group.by = "donor", reduction = "umap")
plot.s2 <- DimPlot(rawdat_9, group.by = "seurat_clusters", reduction = "umap", label = T)
plot.s3 <- DimPlot(rawdat_9, group.by = "group", reduction = "umap", label = T)
plot.s4 <- DimPlot(rawdat_9, group.by = "tissue", reduction = "umap", label = T)
pdf("Graphs/Quality control/Integrationcheck_test_rpca.pdf", height=15, width=20)
plot.s1 + plot.s2 + plot.s3 + plot.s4 + plot_layout(ncol = 2)
dev.off()

test_rpca <- rawdat_9

rawdat_9[["sketch"]] <- JoinLayers(rawdat_9[["sketch"]])

#define Ident
rawdat_9$celltype<-Idents(rawdat_9)

Idents(rawdat_9)<-rawdat_9$celltype

# resplit the sketched cell assay into layers this is required to project the integration onto all cells
rawdat_10 <- rawdat_9
rawdat_10[["sketch"]] <- split(rawdat_10[["sketch"]], f = rawdat_10$donor)

rawdat_10 <- ProjectIntegration(rawdat_10, sketched.assay = "sketch", assay = "RNA", reduction = "integrated.rpca")


rawdat_10 <- ProjectData(rawdat_10, sketched.assay = "sketch", assay = "RNA", sketched.reduction = "integrated.rpca.full",
                         full.reduction = "integrated.rpca.full", dims = 1:30, refdata = list(celltype.full = "celltype"))

rawdat_10 <- RunUMAP(rawdat_10, reduction = "integrated.rpca.full", dims = 1:30, reduction.name = "umap.full",
                     reduction.key = "UMAP_full_")



p1 <- DimPlot(rawdat_10, reduction = "umap.full", group.by = "donor", alpha = 0.1)
p2 <- DimPlot(rawdat_10, reduction = "umap.full", group.by = "celltype.full", alpha = 0.1)
p1 + p2 + plot_layout(ncol = 1)

plot_integrated_clusters(rawdat_10) 


#save rawdat_10
save(rawdat_10,file="rawdat_10.RData")
rm(rawdat_7,rawdat_8,rawdat_9)


DefaultAssay(rawdat_10)<- "RNA"
th.total.uncut<-rawdat_10
th.total.uncut <- JoinLayers(th.total.uncut)
th.total.uncut$orig.cluster<-th.total.uncut$celltype.full

clusters_ordered<-c("0","1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28")
th.total.uncut$celltype.full<- factor(th.total.uncut$celltype.full, levels = clusters_ordered)
Idents(th.total.uncut)<-th.total.uncut$celltype.full
table(th.total.uncut$celltype.full)
table(th.total.uncut$orig.cluster)
table(th.total.uncut$celltype)

#(18)_Module scores based on Xin, Yasumizu, Bautista, Park####
DefaultAssay(th.total.uncut)<-"RNA"
Idents(th.total.uncut)<-"celltype.full"
#Xin: logfc>1.5
TC_Xin_signature <- list(c("CD1E","CD3E","CD3D","CD8B","TRBC2","TCF7","TRBC1","CD8A","DNTT","CD3G","LEF1","CD1A","ELOVL4","CD1B","SATB1","RAG1","STMN1","CD7","LCK","ARPP21"))
th.total.uncut <- AddModuleScore(th.total.uncut, features = TC_Xin_signature, name = "TC_Xin module score")
BC_Xin_signature <- list(c("IGKV3-20","IGKV3-15","IGKV4-1","IGKV3-11","IGKC","IGHV4-39","IGLV2-14","IGHV3-23","IGHG1","IGLV2-8","IGKV1-5","IGHV3-30","IGHV3-21","IGHV1-3","IGHV4-34",
                           "IGKV2-30","IGLV1-47","MS4A1","IGLV3-1","IGHV3-7","CD79A","IGHG2","IGHA1","IGLC2","IGHM","IGHG4","CD83","IGHG3","LY9","JCHAIN","CD74","HLA-DRA"))
th.total.uncut <- AddModuleScore(th.total.uncut, features = BC_Xin_signature, name = "BC_Xin module score")
TEC_Xin_signature <- list(c("KRT19","ELF3","CLDN4","WFDC2","SLPI","LGALS3","KRT17","KRT5","CD9","CLU","CHI3L1","KRT18","KRT8","MDK","CALML3","KRT15","PERP","CD24",
                            "GSTP1","CKB","CSTB","CDC42EP4","ATF3","CALML5","TSPO","S100A2","SFN","S100A14","TACSTD2","SPINT2","PLTP","ANXA4","S100A8","HSPB1",
                            "TFF3"))
th.total.uncut <- AddModuleScore(th.total.uncut, features = TEC_Xin_signature, name = "TEC_Xin module score")
DC_Xin_signature <- list(c("GZMB","LYZ","IRF8","HLA-DPB1","IRF7","PLD4","LILRA4","CPVL","GPR183","CD74","HLA-DQB1","HLA-DRB5","HLA-DPA1","TXN",
                           "HLA-DRA","CST3","HLA-DRB1","HLA-DQA1"))
th.total.uncut <- AddModuleScore(th.total.uncut, features = DC_Xin_signature, name = "DC_Xin module score")
MAC_Xin_signature <- list(c("C1QC","CCL3","CCL4L2","C1QB","C1QA","CCL4","RNASE1","CD14","IL1B","TYROBP","HLA-DPA1","HLA-DRB1","HLA-DRA","MS4A7",
                            "APOE","FCGR3A","HLA-DQA1","MS4A6A","HLA-DRB5","HLA-DPB1","SELENOP","SLC40A1","CD74","IER3","FCER1G","CXCL8","CST3","C3","FTL",
                            "HLA-DQB1","FCGR2A","RGS1","CTSB","HLA-DMA","CSF1R","APOC1","AIF1","FCGRT","SGK1","TNF","NPC2","MS4A4A","FCGBP",
                            "PLAU","SAT1","GPR34","LGMN","FOLR2","HMOX1","GRN","CX3CR1","TREM2","BCL2A1","MAFB"))
th.total.uncut <- AddModuleScore(th.total.uncut, features = MAC_Xin_signature, name = "MAC_Xin module score")
Mono_Xin_signature <- list(c("G0S2","SPP1","PLAUR","LYZ","TYROBP","C15orf48","IL1B","APOE","CXCL8","CTSL","FTL","FCER1G","CXCL2","APOC1","S100A9","CTSB",
                             "BCL2A1","CXCL3","SERPINA1","TIMP1","CCL3","IL1RN","IER3","EREG","RNASE1","CTSD","HLA-DPB1","CST3","HLA-DRA","HLA-DRB5","FTH1","HLA-DRB1",
                             "PLEK","S100A8","TYMP","HLA-DQB1","CCL2","GRN","HLA-DPA1","SOD2","PPIF","SAT1","SDS","HLA-DQA1","LGALS1","IFITM3","PLIN2",
                             "NFKBIA"))
th.total.uncut <- AddModuleScore(th.total.uncut, features = Mono_Xin_signature, name = "Mono_Xin module score")
EC_Xin_signature <- list(c("IGFBP7","SPARCL1","PLVAP","TM4SF1","IFI27","RAMP2","SPARC","IGFBP3","A2M","HSPG2","GNG11","EGFL7","VWF","ACKR1","MGP","IFITM3","CRIP2","PLPP1","CAV1",
                           "ENG","RAMP3","RBP7","AQP1","CLEC14A","IGFBP4","FABP4","NPDC1","SLC9A3R2","ESAM","CLDN5","GSN","COL4A1","PECAM1","CD34","TIMP3","COL4A2","RNASE1","CLEC3B",
                           "ENPP2","BCAM","EMCN","EPAS1","CALCRL","ADGRL4","ID3","CAVIN1","FLT1","EMP1"))
th.total.uncut <- AddModuleScore(th.total.uncut, features = EC_Xin_signature, name = "EC_Xin module score")
FB_Xin_signature <- list(c("SFRP2","COL1A1","TIMP1","MGP","COL1A2","SFRP4","DCN","LUM","COL3A1","FN1","SPARC","RARRES2","APOD","BGN","CCN2","TAGLN","IGFBP4","CXCL14","SERPINF1","CCN1","PLA2G2A",
                           "IGFBP7","CFD","FBLN1","MT2A","CLU","CTHRC1","C1S","IGFBP6","NNMT","CALD1","COL6A2","SERPING1","IFITM3","LGALS1","CXCL12","POSTN",
                           "COL6A3","TPM1","MFAP4","TIMP3","CCDC80","MMP2","MYL9","CCL19","CCL2","IGFBP5","MXRA8","COL6A1","VCAN","PCOLCE","IFI27","TPM2","RARRES1","CTSK",
                           "THBS1","CXCL2","ISLR","EFEMP1","C1R","S100A6","THY1","IGF1","C3","CCN5","CD63","PRSS23","SERPINE1","CST3","CYP1B1","C7"))
th.total.uncut <- AddModuleScore(th.total.uncut, features = FB_Xin_signature, name = "FB_Xin module score")
VSMC_Xin_signature <- list(c("TAGLN","ACTA2","RGS5","IGFBP7","MYL9","CALD1","SPARCL1","TPM2","SPARC","ADIRF","TPM1","MGP","NDUFA4L2","IFITM3","TIMP1","IGFBP4","PPP1R14A","CAV1","MYLK","SOD3","BGN",
                             "CCL2","LHFPL6","SERPING1","PLAC9","COL18A1","COL4A2","CSRP2","THY1","C11orf96","COL4A1","CSRP1","COL1A2","IGFBP5","CCN1","DSTN","MT2A",
                             "MYH11","COL1A1","FN1","COL6A2","COX4I2","CRIP2","MFGE8","GJA4","PDGFRB","ID3","A2M","ACTG2","NOTCH3","RGS16","IFI27","ADAMTS1","FRZB",
                             "CD151","NNMT","COL3A1","LGALS1","CAVIN3","SELENOM","TGFB1I1","CAVIN1","TINAGL1","PTP4A3","MAP1B","CNN3","PHLDA1","HIGD1B","APOE","PCOLCE",
                             "CXCL2","GSN"))
th.total.uncut <- AddModuleScore(th.total.uncut, features = VSMC_Xin_signature, name = "VSMC_Xin module score")

module_DEG_Xin_signature <- c("TC_Xin module score1","BC_Xin module score1","TEC_Xin module score1","DC_Xin module score1","MAC_Xin module score1","Mono_Xin module score1","EC_Xin module score1","FB_Xin module score1","VSMC_Xin module score1")

#Yasumizu upper 20 supp data 10
EC_DEG_yasumizu_signature <- list(c("GSN","IGFBP7","CAV1","SPARCL1","A2M","IFI27","TIMP3","SPARC","VWF","ENG","EGFL7","AQP1","EPAS1",
                                    "ADGRF5","IGFBP4","HSPG2","PECAM1","IFITM3","GNG11","PTPRB"))
th.total.uncut <- AddModuleScore(th.total.uncut, features = EC_DEG_yasumizu_signature, name = "EC_DEG_yasumizu module score")
tumo_assoc_fb_DEG_yasumizu_signature <- list(c("DCN","COL1A2","GSN","CCDC80","CST3","MGP","FBLN1","CFD","C1R","C1S","COL6A2","CCN1",
                                               "MFAP5","FSTL1","SERPING1","C3","CPE","MFAP4","IGFBP5","IGFBP6"))
th.total.uncut <- AddModuleScore(th.total.uncut, features = tumo_assoc_fb_DEG_yasumizu_signature, name = "tumo_assoc_fb_DEG_yasumizu module score")
norm_fb_DEG_yasumizu_signature <- list(c("FN1","CST3","CPE","TSC22D1","CCN2","CD63","SERPINF1","CLU","LGALS1","GAS6","COL6A2","MDK","PRRX1","IGFBP5","FBLN1","IGFBP7","SLPI",
                                         "IFITM3","IGFBP6","MXRA8"))
th.total.uncut <- AddModuleScore(th.total.uncut, features = norm_fb_DEG_yasumizu_signature, name = "norm_fb_DEG_yasumizu module score")
cTEC_DEG_yasumizu_signature <- list(c("PLTP","PRSS16","CCL25","GAS6","CST3","IGFBP2","APP","KRT5","IGFBP7","MDK","SAA1","KRT19","DSP",
                                      "PAX1","CYP1B1","COL6A1","WFDC2","COL6A2","SERPINF1"))
th.total.uncut <- AddModuleScore(th.total.uncut, features = cTEC_DEG_yasumizu_signature, name = "cTEC_DEG_yasumizu module score")
mTEC_I_DEG_yasumizu_signature <- list(c("KRT19","S100A14","IFI27","KRT5","CALML3","CTSV","MT1X","CST3","MDK","BST2",
                                        "KRT15","ANXA2","IFITM3","TSC22D1","HSPB1","LGALS1","IFI6","SAA1","CCL25"))
th.total.uncut <- AddModuleScore(th.total.uncut, features = mTEC_I_DEG_yasumizu_signature, name = "mTEC_I_DEG_yasumizu module score")
mTEC_II_DEG_yasumizu_signature <- list(c("CLU","APP","NFIA","SPINT2","BCAM","KRT19","KRT8","CLDN4","SCPEP1","S100A14","WFDC2","SLPI",
                                         "CD63","CST3","MDK","SDC4","TSC22D1","KRT18","CLDN7","AQP5"))
th.total.uncut <- AddModuleScore(th.total.uncut, features = mTEC_II_DEG_yasumizu_signature, name = "mTEC_II_DEG_yasumizu module score")
nmTEC_DEG_yasumizu_signature <- list(c("KRT19","S100A14","NFIB","CNN3","LGALS3","CLU","MDK","RBMS3","IFI27","KRT5","APLP2","NFIA","CEBPD",
                                       "BST2","CALML3","HSPB1","IFITM3","PERP","CD63","KRT8"))
th.total.uncut <- AddModuleScore(th.total.uncut, features = nmTEC_DEG_yasumizu_signature, name = "nmTEC_DEG_yasumizu module score")
NK_DEG_yasumizu_signature <- list(c("NKG7","GNLY","CST7","KLRD1","GZMB","PRF1","GZMA","CTSW","FGFBP2","TYROBP","CCL5","KLRF1","HOPX",
                                    "FCGR3A","KLRB1","CCL4","GZMH","HCST","GZMM","EFHD2"))
th.total.uncut <- AddModuleScore(th.total.uncut, features = NK_DEG_yasumizu_signature, name = "NK_DEG_yasumizu module score")
NKT_thym_DEG_yasumizu_signature <- list(c("CTSW","FOSB","XCL1","FOS","HCST","TRGC2","SRGN","HSPA1A","ZNF683","HSPA1B","STK17A","KLRK1",
                                          "MATK","JUNB","XCL2","JUN","DUSP1","CD69","CD7","CD63"))
th.total.uncut <- AddModuleScore(th.total.uncut, features = NKT_thym_DEG_yasumizu_signature, name = "NKT_thym_DEG_yasumizu module score")
NKT_peri_DEG_yasumizu_signature <- list(c("LTB","ZBTB16","TPT1","CRIP1","IL32","GIMAP4","IL7R","ANXA1","KLF2","PIM1",
                                          "B2M","RPL10","RPS29","S100A4","EMP3","RPS27A","RPS27","RPL34","RPS14"))
th.total.uncut <- AddModuleScore(th.total.uncut, features = NKT_peri_DEG_yasumizu_signature, name = "NKT_peri_DEG_yasumizu module score")
gdTC_DEG_yasumizu_signature <- list(c("RTKN2","TRDC","TRGC2","TRG-AS1","LIME1","PPP2R5C","SLFN5","LEF1","TOX2","IKZF2","CD7","GYPC","MYOM2","ACTN1",
                                      "HCST","PRKCH","CD27","TRGC1","LAT"))
th.total.uncut <- AddModuleScore(th.total.uncut, features = gdTC_DEG_yasumizu_signature, name = "gdTC_DEG_yasumizu module score")
ILC_DEG_yasumizu_signature <- list(c("KLRB1","CTSW","XCL1","SPINK2","NFKBIA","CD63","XCL2","CD69","BST2","TNFRSF18","FXYD5","LTC4S",
                                     "ID2","LINC00299","FXYD7","FOS","TNFRSF4","MDFIC","KRT81","TRDC"))
th.total.uncut <- AddModuleScore(th.total.uncut, features = ILC_DEG_yasumizu_signature, name = "ILC_DEG_yasumizu module score")
DN_TC_DEG_yasumizu_signature <- list(c("CDK6","GAPDH","GSTP1","IGLL1","FXYD2","HSPD1","YBX3",
                                       "FABP5","PPP1R14B","UHRF1","HSP90AB1","CD99","NPM1","SLC25A3","SMIM24","GSTO1"))
th.total.uncut <- AddModuleScore(th.total.uncut, features = DN_TC_DEG_yasumizu_signature, name = "DN_TC_DEG_yasumizu module score")
DP_TC_DEG_yasumizu_signature <- list(c("CD1E","MZB1","ARPP21","CD1B","RAG1","CD8B","TFDP2","CD8A","CD99","DNTT","ELOVL4","MDM4","PTP4A2",
                                       "SOX4","CD1A","RAG2","CD3D","CD38","CD3G"))
th.total.uncut <- AddModuleScore(th.total.uncut, features = DP_TC_DEG_yasumizu_signature, name = "DP_TC_DEG_yasumizu module score")
cyc_DNDP_TC_DEG_yasumizu_signature <- list(c("HMGB2","TUBA1B","PCLAF","TYMS","TUBB","MKI67","HMGB1","HMGN2","NUSAP1",
                                             "STMN1","TOP2A","SMC2","CENPF","RRM2","SMC4","BIRC5","DEK","NUCKS1"))
th.total.uncut <- AddModuleScore(th.total.uncut, features = cyc_DNDP_TC_DEG_yasumizu_signature, name = "cyc_DNDP_TC_DEG_yasumizu module score")
aaCD8_I_TC_DEG_yasumizu_signature <- list(c("LEF1","TRAC","DUSP2","TOX2","CD2","ITM2A","CD3D","NUCB2","ACTB","PTPN7","TOX","SLAMF1","SRGN",
                                            "SATB1","ACTG1","IKZF2","ITGA4","MALAT1","ID3","MIR181A1HG"))
th.total.uncut <- AddModuleScore(th.total.uncut, features = aaCD8_I_TC_DEG_yasumizu_signature, name = "aaCD8_I_TC_DEG_yasumizu module score")
aaCD8_II_TC_DEG_yasumizu_signature <- list(c("CTSW","HCST","PPP2R5C","B2M","TRGC2","TXNIP","KLRK1","IL32","GZMK","CD7","HLA-B",
                                             "CD27","LYAR","IFITM2","DUSP2","HLA-C","TMSB4X","GIMAP4","NKG7","MBP"))
th.total.uncut <- AddModuleScore(th.total.uncut, features = aaCD8_II_TC_DEG_yasumizu_signature, name = "aaCD8_II_TC_DEG_yasumizu module score")
thym_CD8_TC_DEG_yasumizu_signature <- list(c("ARMH1","GZMM","SIRPG","CLEC2D","TOX2","CHI3L2","ABLIM1","CD3D","SATB1","TRAC","LCP2","LAT",
                                             "ACTN1","STK17B","CD3E","PTPRC","CD2","DENND2D","OXNAD1","CORO1A"))
th.total.uncut <- AddModuleScore(th.total.uncut, features = thym_CD8_TC_DEG_yasumizu_signature, name = "thym_CD8_TC_DEG_yasumizu module score")
naiveTC_TET1_TC_DEG_yasumizu_signature <- list(c("RPS10","RPS29","RPL34","RPS3A","RPS27","RPS12","RPL39","NOSIP","RPS18","RPL13",
                                                 "RPS4X","RPS27A","RPL21","RPS15A","RPS28","RPS21","RPS14","RPL9","RPL38"))
th.total.uncut <- AddModuleScore(th.total.uncut, features = naiveTC_TET1_TC_DEG_yasumizu_signature, name = "naiveTC_TET1_TC_DEG_yasumizu module score")
CD8_Tem_DEG_yasumizu_signature <- list(c("GZMK","CCL5","GZMA","CST7","IL32","B2M","HLA-B","HLA-A","HLA-C","HCST","S100A4","GZMM","NKG7",
                                         "GIMAP4","S100A10","LYAR","S100A6","HLA-E","GIMAP7","CALM1"))
th.total.uncut <- AddModuleScore(th.total.uncut, features = CD8_Tem_DEG_yasumizu_signature, name = "CD8_Tem_DEG_yasumizu module score")
CD8_Trm_DEG_yasumizu_signature <- list(c("CCL5","GZMK","GZMA","CCL4","CST7","SRGN","B2M","IL32","NKG7","GZMH","CCL4L2","HLA-C","HLA-A",
                                         "HLA-B","APOBEC3G","RGS1","HCST","CD69","CD2","DUSP2"))
th.total.uncut <- AddModuleScore(th.total.uncut, features = CD8_Trm_DEG_yasumizu_signature, name = "CD8_Trm_DEG_yasumizu module score")
CD8_Temra_DEG_yasumizu_signature <- list(c("NKG7","CCL5","CST7","GNLY","GZMH","GZMA","GZMB","KLRD1","CTSW","PRF1","S100A4","HCST",
                                           "CCL4","GZMM","HOPX","B2M","IL32","FGFBP2","IFITM2","PPP2R5C"))
th.total.uncut <- AddModuleScore(th.total.uncut, features = CD8_Temra_DEG_yasumizu_signature, name = "CD8_Temra_DEG_yasumizu module score")
thym_CD4_I_TC_DEG_yasumizu_signature <- list(c("ITM2A","SATB1","TOX2","CCR9","BCL11B","TRAC","CD3D","CHI3L2","ARMH1","SOX4","PTPRC",
                                               "ACTB","FKBP5","LZTFL1","MALAT1","LAT","CD3E","CD3G","NDFIP1","CD2"))
th.total.uncut <- AddModuleScore(th.total.uncut, features = thym_CD4_I_TC_DEG_yasumizu_signature, name = "thym_CD4_I_TC_DEG_yasumizu module score")
thym_CD4_II_TC_DEG_yasumizu_signature <- list(c("ITM2A","ARMH1","IL6ST","TMSB10","NDFIP1","SARAF","CD40LG","TRAC","TPT1","HLA-E",
                                                "OXNAD1","ID2","CYTIP","RPL23","CLDN1","SCML4","RCAN3","CD27","SERPINE2","RPS24"))
th.total.uncut <- AddModuleScore(th.total.uncut, features = thym_CD4_II_TC_DEG_yasumizu_signature, name = "thym_CD4_II_TC_DEG_yasumizu module score")
CD4_Tnaive_DEG_yasumizu_signature <- list(c("GIMAP7","CCR7","IFITM1","NOSIP","KLF2","GIMAP5","LINC00861","SELL","HLA-B","RPL30",
                                            "SARAF","GIMAP1","RPL32","HLA-C","RPS14","RPS12","LDLRAP1","RPS25","EEF1A1"))
th.total.uncut <- AddModuleScore(th.total.uncut, features = CD4_Tnaive_DEG_yasumizu_signature, name = "CD4_Tnaive_DEG_yasumizu module score")
CD4_Tcm_Th0_DEG_yasumizu_signature <- list(c("GIMAP7","HLA-B","IL7R","IL32","GIMAP5","GIMAP4","B2M","HLA-C","FXYD5","KLF2","HLA-A","RPS12","RPL32","RPL30","RPS14","HLA-E","GPR183",
                                             "RPS25","ANXA1","RPS27A"))
th.total.uncut <- AddModuleScore(th.total.uncut, features = CD4_Tcm_Th0_DEG_yasumizu_signature, name = "CD4_Tcm_Th0_DEG_yasumizu module score")
CD4_Tcm_Th2_DEG_yasumizu_signature <- list(c("CRIP1","S100A4","ANXA1","IL32","FXYD5","S100A10","GSTK1","VIM","EMP3","GIMAP7","AHNAK",
                                             "HLA-C","KLF2","B2M","ITGB1","HLA-A","HLA-B","S100A6","IFITM1","COTL1"))
th.total.uncut <- AddModuleScore(th.total.uncut, features = CD4_Tcm_Th2_DEG_yasumizu_signature, name = "CD4_Tcm_Th2_DEG_yasumizu module score")
T_ago_DEG_yasumizu_signature <- list(c("SARAF","NDFIP1","IFITM1","CD5","HLA-E","ARMH1","CD6","LIMD2","SIRPG","CORO1A","CCR7","NR4A1",
                                       "B2M","CD27","ACTB","CD3E","DENND2D","TSHZ2","PTPRC","HLA-A"))
th.total.uncut <- AddModuleScore(th.total.uncut, features = T_ago_DEG_yasumizu_signature, name = "T_ago_DEG_yasumizu module score")
CD4_Tcm_Tfh_DEG_yasumizu_signature <- list(c("GIMAP7","IFITM1","HLA-A","GIMAP5","FXYD5","HLA-B","HLA-C","IL32","B2M","GIMAP4",
                                             "GSTK1","LIMS1","CRIP1","HLA-E","CD27","S100A4","FYB1","SPOCK2","KLF2","GPR183"))
th.total.uncut <- AddModuleScore(th.total.uncut, features = CD4_Tcm_Tfh_DEG_yasumizu_signature, name = "CD4_Tcm_Tfh_DEG_yasumizu module score")
CD4_Tcm_Th17_DEG_yasumizu_signature <- list(c("S100A4","CRIP1","IL32","S100A10","AHNAK","IFITM1","ANXA1","FXYD5","EMP3","VIM",
                                              "S100A6","GIMAP7","HLA-C","HLA-A","HLA-B","GSTK1","B2M","KLF2","FLT3LG","CALM1"))
th.total.uncut <- AddModuleScore(th.total.uncut, features = CD4_Tcm_Th17_DEG_yasumizu_signature, name = "CD4_Tcm_Th17_DEG_yasumizu module score")
CD4_Tem_Th1_17_DEG_yasumizu_signature <- list(c("S100A4","IL32","AHNAK","ANXA1","IFITM1","IL7R","HLA-C","HLA-A","HLA-B","KLRB1",
                                                "CCL5","GZMK","B2M","S100A10","KLRG1","GIMAP7","TGFB1","GZMA","CRIP1","FXYD5"))
th.total.uncut <- AddModuleScore(th.total.uncut, features = CD4_Tem_Th1_17_DEG_yasumizu_signature, name = "CD4_Tem_Th1_17_DEG_yasumizu module score")
CD4_Tem_Th1_DEG_yasumizu_signature <- list(c("CCL5","GZMA","CST7","NKG7","PRF1","IL32","S100A4","GZMK","GZMH","B2M","GZMM",
                                             "HLA-C","KLRG1","HLA-A","AHNAK","ANXA1","HLA-B","EMP3","S100A10","HCST"))
th.total.uncut <- AddModuleScore(th.total.uncut, features = CD4_Tem_Th1_DEG_yasumizu_signature, name = "CD4_Tem_Th1_DEG_yasumizu module score")
CD4_Temra_Th1_DEG_yasumizu_signature <- list(c("GZMH","GNLY","NKG7","CCL5","CST7","GZMA","PRF1","FGFBP2","S100A4","CX3CR1",
                                               "EFHD2","KLRG1","IL32","GZMB","IFITM1","AHNAK","HOPX","FLNA","S100A10","ANXA1"))
th.total.uncut <- AddModuleScore(th.total.uncut, features = CD4_Temra_Th1_DEG_yasumizu_signature, name = "CD4_Temra_Th1_DEG_yasumizu module score")
naive_Treg_DEG_yasumizu_signature <- list(c("IL32","HLA-A","HLA-B","TIGIT","HLA-C","FOXP3","S100A4","B2M","IL10RA","FXYD5","CD27","IFITM1",
                                            "HLA-E","GBP5","SAMHD1","GIMAP7","SELL","TNFRSF1B","S100A10","RPS26"))
th.total.uncut <- AddModuleScore(th.total.uncut, features = naive_Treg_DEG_yasumizu_signature, name = "naive_Treg_DEG_yasumizu module score")
act_Treg_DEG_yasumizu_signature <- list(c("IL32","S100A4","S100A10","HLA-A","SAMHD1","B2M","CRIP1","EMP3","GSTK1","HLA-C","ANXA2",
                                          "AHNAK","ITGB1","HLA-B","IL10RA","S100A6","GBP5","CYTOR","KLF2","FXYD5"))
th.total.uncut <- AddModuleScore(th.total.uncut, features = act_Treg_DEG_yasumizu_signature, name = "act_Treg_DEG_yasumizu module score")
cDC1_DEG_yasumizu_signature <- list(c("CST3","LYZ","CPVL","HLA-DPA1","HLA-DPB1","HLA-DQB1","LGALS1","HLA-DRA","IRF8","BASP1","S100A10",
                                      "HLA-DRB1","FGL2","HLA-DQA1","C1orf54","GSTP1","CD74","COTL1","BATF3","CTSZ"))
th.total.uncut <- AddModuleScore(th.total.uncut, features = cDC1_DEG_yasumizu_signature, name = "cDC1_DEG_yasumizu module score")
cDC2_DEG_yasumizu_signature <- list(c("CST3","TYROBP","FCER1G","HLA-DRB1","HLA-DPB1","HLA-DQB1","HLA-DPA1","LYZ","FGL2","HLA-DRA",
                                      "LST1","HLA-DQA1","CPVL","HLA-DMA","GRN","MNDA","FCGRT","IGSF6","LGALS1","JAML"))
th.total.uncut <- AddModuleScore(th.total.uncut, features = cDC2_DEG_yasumizu_signature, name = "cDC2_DEG_yasumizu module score")
pDC_DEG_yasumizu_signature <- list(c("CCDC50","TCF4","IRF7","IRF8","LILRA4","UGCG","TPM2","JCHAIN","PLD4","IGKC","SERPINF1","GZMB",
                                     "IL3RA","TYROBP","CST3","HERPUD1","SCT","PLEK","MPEG1","CLEC4C"))
th.total.uncut <- AddModuleScore(th.total.uncut, features = pDC_DEG_yasumizu_signature, name = "pDC_DEG_yasumizu module score")
mono_CD14_DEG_yasumizu_signature <- list(c("TYROBP","LYZ","CST3","S100A9","MNDA","FCN1","S100A8","CTSS","FCER1G","TYMP","VCAN",
                                           "CD14","FTL","LGALS1","LST1","CEBPD","IFI30","SERPINA1","AIF1","S100A12"))
th.total.uncut <- AddModuleScore(th.total.uncut, features = mono_CD14_DEG_yasumizu_signature, name = "mono_CD14_DEG_yasumizu module score")
mono_CD16_DEG_yasumizu_signature <- list(c("LST1","FCER1G","TYROBP","CTSS","FCN1","IFI30","CST3","SERPINA1","AIF1","TYMP","COTL1",
                                           "FCGR3A","FTL","SAT1","TNFRSF1B","FTH1","S100A4","BRI3","MS4A7"))
th.total.uncut <- AddModuleScore(th.total.uncut, features = mono_CD16_DEG_yasumizu_signature, name = "mono_CD16_DEG_yasumizu module score")
Mac_DEG_yasumizu_signature <- list(c("C1QA","C1QB","C1QC","MS4A7","MS4A6A","CST3","TYROBP","A2M","CD14","SAT1","NPC2","SLC40A1","FCGRT",
                                     "CTSB","HLA-DRA","HLA-DPA1","SELENOP","MAF","APOE","HLA-DPB1"))
th.total.uncut <- AddModuleScore(th.total.uncut, features = Mac_DEG_yasumizu_signature, name = "Mac_DEG_yasumizu module score")
naive_BC_DEG_yasumizu_signature <- list(c("CD74","MS4A1","TCL1A","HLA-DRB5","CD79A","HLA-DRA","FCER2","IGHM","CD79B","HLA-DPA1",
                                          "HLA-DQA1","HLA-DQA2","HVCN1","HLA-DPB1","LINC00926","HLA-DRB1","BANK1","NCF1","CD37","HLA-DMA"))
th.total.uncut <- AddModuleScore(th.total.uncut, features = naive_BC_DEG_yasumizu_signature, name = "naive_BC_DEG_yasumizu module score")
pre_GC_BC_DEG_yasumizu_signature <- list(c("MS4A1","HLA-DRA","HLA-DRB5","CD74","CD79A","HLA-DQA2","HLA-DPA1","HLA-DPB1","HLA-DQA1",
                                           "CD79B","TCL1A","HLA-DRB1","IGHM","NCF1","HLA-DQB1","HLA-DMA","FCER2","HVCN1","BANK1","HLA-DMB"))
th.total.uncut <- AddModuleScore(th.total.uncut, features = pre_GC_BC_DEG_yasumizu_signature, name = "pre_GC_BC_DEG_yasumizu module score")
GC_BC_DEG_yasumizu_signature <- list(c("MS4A1","CD79A","HLA-DRA","CD74","HLA-DRB5","CD79B","HLA-DQA1","TCL1A","HLA-DQA2","HLA-DPB1",
                                       "HLA-DPA1","HLA-DMB","HLA-DMA","NCF1","HLA-DRB1","MEF2C","HLA-DQB1","CD37","CD22","IGHM"))
th.total.uncut <- AddModuleScore(th.total.uncut, features = GC_BC_DEG_yasumizu_signature, name = "GC_BC_DEG_yasumizu module score")
PB_DEG_yasumizu_signature <- list(c("CD74","CD79A","MS4A1","HLA-DRB5","HLA-DRA","CD79B","HLA-DPA1","HLA-DQA2","IGHG3","BANK1","HLA-DQA1",
                                    "LINC00926","HLA-DPB1","HLA-DRB1","IGHM","RALGPS2","NCF1","NIBAN3","HLA-DMA","HVCN1"))
th.total.uncut <- AddModuleScore(th.total.uncut, features = PB_DEG_yasumizu_signature, name = "PB_DEG_yasumizu module score")
unswi_mem_BC_DEG_yasumizu_signature <- list(c("MS4A1","BANK1","CD79A","HLA-DRB5","CD74","CD79B","HLA-DRA","HLA-DPA1","HLA-DPB1","HLA-DQA2",
                                              "HLA-DQA1","HLA-DRB1","RALGPS2","CD37","HLA-DQB1","NIBAN3","FCRLA","HLA-DMA","LINC00926","NCF1"))
th.total.uncut <- AddModuleScore(th.total.uncut, features = unswi_mem_BC_DEG_yasumizu_signature, name = "unswi_mem_BC_DEG_yasumizu module score")
thym_mem_BC_DEG_yasumizu_signature <- list(c("MS4A1","HLA-DRA","HLA-DPB1","CD74","HLA-DPA1","HLA-DRB5","CD79A","HLA-DQA2","BANK1",
                                             "HLA-DQB1","CD37","HLA-DRB1","HLA-DQA1","HLA-DMA","CD79B","VPREB3","BLK","TNFRSF13C","CD24"))
th.total.uncut <- AddModuleScore(th.total.uncut, features = thym_mem_BC_DEG_yasumizu_signature, name = "thym_mem_BC_DEG_yasumizu module score")
mem_BC_I_DEG_yasumizu_signature <- list(c("MS4A1","BANK1","CD74","CD79A","HLA-DRA","HLA-DRB5","CD79B","RALGPS2","HLA-DPA1","HLA-DPB1",
                                          "HLA-DQA1","CD37","HLA-DQA2","HLA-DRB1","LINC00926","NIBAN3","ALOX5","SPIB","CD24","NCF1"))
th.total.uncut <- AddModuleScore(th.total.uncut, features = mem_BC_I_DEG_yasumizu_signature, name = "mem_BC_I_DEG_yasumizu module score")
mem_BC_II_DEG_yasumizu_signature <- list(c("MS4A1","CD74","HLA-DRA","BANK1","HLA-DRB5","HLA-DPB1","CD79A","HLA-DPA1","HLA-DQA2","HLA-DQA1","HLA-DRB1",
                                           "BLK","HLA-DQB1","SYNGR2","CD79B","CRIP1","CD37","POU2AF1","SPIB","HLA-DMA"))
th.total.uncut <- AddModuleScore(th.total.uncut, features = mem_BC_II_DEG_yasumizu_signature, name = "mem_BC_II_DEG_yasumizu module score")

module_DEG_yasumizu_signature <- c("EC_DEG_yasumizu module score1","tumo_assoc_fb_DEG_yasumizu module score1","norm_fb_DEG_yasumizu module score1","cTEC_DEG_yasumizu module score1","mTEC_I_DEG_yasumizu module score1","mTEC_II_DEG_yasumizu module score1","nmTEC_DEG_yasumizu module score1","NK_DEG_yasumizu module score1","NKT_thym_DEG_yasumizu module score1","NKT_peri_DEG_yasumizu module score1","gdTC_DEG_yasumizu module score1","ILC_DEG_yasumizu module score1","DN_TC_DEG_yasumizu module score1","DP_TC_DEG_yasumizu module score1","cyc_DNDP_TC_DEG_yasumizu module score1","aaCD8_I_TC_DEG_yasumizu module score1","aaCD8_II_TC_DEG_yasumizu module score1","thym_CD8_TC_DEG_yasumizu module score1","naiveTC_TET1_TC_DEG_yasumizu module score1","CD8_Tem_DEG_yasumizu module score1","CD8_Trm_DEG_yasumizu module score1","CD8_Temra_DEG_yasumizu module score1","thym_CD4_I_TC_DEG_yasumizu module score1","thym_CD4_II_TC_DEG_yasumizu module score1","CD4_Tnaive_DEG_yasumizu module score1","CD4_Tcm_Th0_DEG_yasumizu module score1","CD4_Tcm_Th2_DEG_yasumizu module score1","T_ago_DEG_yasumizu module score1","CD4_Tcm_Tfh_DEG_yasumizu module score1","CD4_Tcm_Th17_DEG_yasumizu module score1", "CD4_Tem_Th1_17_DEG_yasumizu module score1","CD4_Tem_Th1_DEG_yasumizu module score1","CD4_Temra_Th1_DEG_yasumizu module score1","naive_Treg_DEG_yasumizu module score1","act_Treg_DEG_yasumizu module score1","cDC1_DEG_yasumizu module score1","cDC2_DEG_yasumizu module score1","pDC_DEG_yasumizu module score1","mono_CD14_DEG_yasumizu module score1","mono_CD16_DEG_yasumizu module score1","Mac_DEG_yasumizu module score1","naive_BC_DEG_yasumizu module score1","pre_GC_BC_DEG_yasumizu module score1","GC_BC_DEG_yasumizu module score1","PB_DEG_yasumizu module score1", "unswi_mem_BC_DEG_yasumizu module score1","thym_mem_BC_DEG_yasumizu module score1","mem_BC_I_DEG_yasumizu module score1","mem_BC_II_DEG_yasumizu module score1")


#Yasumizu marker genes
EC_yasumizu_signature <- list(c("PECAM1","VWF","A2M","ADGRF5","AQP1"))
th.total.uncut <- AddModuleScore(th.total.uncut, features = EC_yasumizu_signature, name = "EC_yasumizu module score")
tumo_assoc_fb_yasumizu_signature <- list(c("PDGFRA","ADH1B","FBN1","COL1A1","S100A4"))
th.total.uncut <- AddModuleScore(th.total.uncut, features = tumo_assoc_fb_yasumizu_signature, name = "tumo_assoc_fb_yasumizu module score")
norm_fb_yasumizu_signature <- list(c("FN1","EGFL6","MATN4","PTN","VIM"))
th.total.uncut <- AddModuleScore(th.total.uncut, features = norm_fb_yasumizu_signature, name = "norm_fb_yasumizu module score")
cTEC_yasumizu_signature <- list(c("CCL25","PSMB11","LY75","PRSS16","PMP22"))
th.total.uncut <- AddModuleScore(th.total.uncut, features = cTEC_yasumizu_signature, name = "cTEC_yasumizu module score")
mTEC_I_yasumizu_signature <- list(c("KRT15","IFI27","MT1X","CBR3","GPX2"))
th.total.uncut <- AddModuleScore(th.total.uncut, features = mTEC_I_yasumizu_signature, name = "mTEC_I_yasumizu module score")
mTEC_II_yasumizu_signature <- list(c("CLDN4","SDC4","CLU","STEAP4","ALPI"))
th.total.uncut <- AddModuleScore(th.total.uncut, features = mTEC_II_yasumizu_signature, name = "mTEC_II_yasumizu module score")
nmTEC_yasumizu_signature <- list(c("GABRA5","NEFL","KRT15","KRT6C","IFI6"))
th.total.uncut <- AddModuleScore(th.total.uncut, features = nmTEC_yasumizu_signature, name = "nmTEC_yasumizu module score")
NK_yasumizu_signature <- list(c("GNLY","NKG7","TYROBP","KLRF1"))
th.total.uncut <- AddModuleScore(th.total.uncut, features = NK_yasumizu_signature, name = "NK_yasumizu module score")
NKT_thym_yasumizu_signature <- list(c("XCL1","FOSB","JUN","XCL2","HSPA1A"))
th.total.uncut <- AddModuleScore(th.total.uncut, features = NKT_thym_yasumizu_signature, name = "NKT_thym_yasumizu module score")
NKT_peri_yasumizu_signature <- list(c("ZBTB16","ZNF462","SOCS3","TNFSF10","TRAT1"))
th.total.uncut <- AddModuleScore(th.total.uncut, features = NKT_peri_yasumizu_signature, name = "NKT_peri_yasumizu module score")
gdTC_yasumizu_signature <- list(c("TRDC","RTKN2","SMC4","MME","TOX2"))
th.total.uncut <- AddModuleScore(th.total.uncut, features = gdTC_yasumizu_signature, name = "gdTC_yasumizu module score")
ILC_yasumizu_signature <- list(c("SPINK2","KLRB1","CD63","FCER1G","BST2"))
th.total.uncut <- AddModuleScore(th.total.uncut, features = ILC_yasumizu_signature, name = "ILC_yasumizu module score")
DN_TC_yasumizu_signature <- list(c("CDK6","IGLL1","DNTT","STMN1","NPM1"))
th.total.uncut <- AddModuleScore(th.total.uncut, features = DN_TC_yasumizu_signature, name = "DN_TC_yasumizu module score")
DP_TC_yasumizu_signature <- list(c("RAG1","CD1B","RAG2","ARPP21","MDM4"))
th.total.uncut <- AddModuleScore(th.total.uncut, features = DP_TC_yasumizu_signature, name = "DP_TC_yasumizu module score")
cyc_DNDP_TC_yasumizu_signature <- list(c("HMGB2","TUBA1B","STMN1","HMGN2","TUBB"))
th.total.uncut <- AddModuleScore(th.total.uncut, features = cyc_DNDP_TC_yasumizu_signature, name = "cyc_DNDP_TC_yasumizu module score")
aaCD8_I_TC_yasumizu_signature <- list(c("CD8A","GNG4","SLAMF1","MYB","SH2D1A"))
th.total.uncut <- AddModuleScore(th.total.uncut, features = aaCD8_I_TC_yasumizu_signature, name = "aaCD8_I_TC_yasumizu module score")
aaCD8_II_TC_yasumizu_signature <- list(c("ZNF683","KLRK1","LYAR","GZMK","NCR3"))
th.total.uncut <- AddModuleScore(th.total.uncut, features = aaCD8_II_TC_yasumizu_signature, name = "aaCD8_II_TC_yasumizu module score")
thym_CD8_TC_yasumizu_signature <- list(c("ARMH1","CHI3L1","SATB1","SIRPG","CD8A"))
th.total.uncut <- AddModuleScore(th.total.uncut, features = thym_CD8_TC_yasumizu_signature, name = "thym_CD8_TC_yasumizu module score")
naiveTC_TET1_TC_yasumizu_signature <- list(c("CD8B","KLRK1","S100B","LRRN3","PASK"))
th.total.uncut <- AddModuleScore(th.total.uncut, features = naiveTC_TET1_TC_yasumizu_signature, name = "naiveTC_TET1_TC_yasumizu module score")
CD8_Tem_yasumizu_signature <- list(c("EOMES","GZMK","CXCR4","CXCR3","CD84"))
th.total.uncut <- AddModuleScore(th.total.uncut, features = CD8_Tem_yasumizu_signature, name = "CD8_Tem_yasumizu module score")
CD8_Trm_yasumizu_signature <- list(c("RGS1","CCL4L2","CXCR6","ITGA1"))
th.total.uncut <- AddModuleScore(th.total.uncut, features = CD8_Trm_yasumizu_signature, name = "CD8_Trm_yasumizu module score")
CD8_Temra_yasumizu_signature <- list(c("FCGR3A","NKG7","GNLY","KLRD1","CCL4"))
th.total.uncut <- AddModuleScore(th.total.uncut, features = CD8_Temra_yasumizu_signature, name = "CD8_Temra_yasumizu module score")
thym_CD4_I_TC_yasumizu_signature <- list(c("SATB1","CCR9","TOX2","LZTFL1","FOS"))
th.total.uncut <- AddModuleScore(th.total.uncut, features = thym_CD4_I_TC_yasumizu_signature, name = "thym_CD4_I_TC_yasumizu module score")
thym_CD4_II_TC_yasumizu_signature <- list(c("PRKD3","ARMH1","CLDN1","ID2","CD69"))
th.total.uncut <- AddModuleScore(th.total.uncut, features = thym_CD4_II_TC_yasumizu_signature, name = "thym_CD4_II_TC_yasumizu module score")
CD4_Tnaive_yasumizu_signature <- list(c("CCR7","NOG","SELL","AK5","TCF7"))
th.total.uncut <- AddModuleScore(th.total.uncut, features = CD4_Tnaive_yasumizu_signature, name = "CD4_Tnaive_yasumizu module score")
CD4_Tcm_Th0_yasumizu_signature <- list(c("FOSB","FOS","CD69","CD84","GCNT4"))
th.total.uncut <- AddModuleScore(th.total.uncut, features = CD4_Tcm_Th0_yasumizu_signature, name = "CD4_Tcm_Th0_yasumizu module score")
CD4_Tcm_Th2_yasumizu_signature <- list(c("GATA3","CCR4","CD84","ITGB1","ANXA1"))
th.total.uncut <- AddModuleScore(th.total.uncut, features = CD4_Tcm_Th2_yasumizu_signature, name = "CD4_Tcm_Th2_yasumizu module score")
T_ago_yasumizu_signature <- list(c("NR4A1","NFKBID","BCL2L11","REL","BACH2"))
th.total.uncut <- AddModuleScore(th.total.uncut, features = T_ago_yasumizu_signature, name = "T_ago_yasumizu module score")
CD4_Tcm_Tfh_yasumizu_signature <- list(c("CXCR5","MAF","TIGIT","IKZF3","PDCD1"))
th.total.uncut <- AddModuleScore(th.total.uncut, features = CD4_Tcm_Tfh_yasumizu_signature, name = "CD4_Tcm_Tfh_yasumizu module score")
CD4_Tcm_Th17_yasumizu_signature <- list(c("RORC","TNFRSF4","CCR6","LGALS1","ANXA2"))
th.total.uncut <- AddModuleScore(th.total.uncut, features = CD4_Tcm_Th17_yasumizu_signature, name = "CD4_Tcm_Th17_yasumizu module score")
CD4_Tem_Th1_17_yasumizu_signature <- list(c("GZMK","KLRB1","CCL5","CEBPD","CXCR3"))
th.total.uncut <- AddModuleScore(th.total.uncut, features = CD4_Tem_Th1_17_yasumizu_signature, name = "CD4_Tem_Th1_17_yasumizu module score")
CD4_Tem_Th1_yasumizu_signature <- list(c("GZMH","CST7","EOMES","S1PR5","CXCR3"))
th.total.uncut <- AddModuleScore(th.total.uncut, features = CD4_Tem_Th1_yasumizu_signature, name = "CD4_Tem_Th1_yasumizu module score")
CD4_Temra_Th1_yasumizu_signature <- list(c("TBX21","CX3CR1","GNLY","FGFBP2","CCL5"))
th.total.uncut <- AddModuleScore(th.total.uncut, features = CD4_Temra_Th1_yasumizu_signature, name = "CD4_Temra_Th1_yasumizu module score")
naive_Treg_yasumizu_signature <- list(c("FOXP3","IKZF2","FCRL3","IL2RA","CTLA4"))
th.total.uncut <- AddModuleScore(th.total.uncut, features = naive_Treg_yasumizu_signature, name = "naive_Treg_yasumizu module score")
act_Treg_yasumizu_signature <- list(c("CYTOR","CCR4","IKZF2","FANK1","FOXP3"))
th.total.uncut <- AddModuleScore(th.total.uncut, features = act_Treg_yasumizu_signature, name = "act_Treg_yasumizu module score")
cDC1_yasumizu_signature <- list(c("XCR1","CLEC9A","BATF3","RGCC","IRF8"))
th.total.uncut <- AddModuleScore(th.total.uncut, features = cDC1_yasumizu_signature, name = "cDC1_yasumizu module score")
cDC2_yasumizu_signature <- list(c("SIRPA","FCER1A","CST3","CD74"))
th.total.uncut <- AddModuleScore(th.total.uncut, features = cDC2_yasumizu_signature, name = "cDC2_yasumizu module score")
pDC_yasumizu_signature <- list(c("CLEC4C","IL3RA","IRF7","TCF4","JCHAIN"))
th.total.uncut <- AddModuleScore(th.total.uncut, features = pDC_yasumizu_signature, name = "pDC_yasumizu module score")
mono_CD14_yasumizu_signature <- list(c("CD14","S100A8","S100A9","VCAN","FCN1"))
th.total.uncut <- AddModuleScore(th.total.uncut, features = mono_CD14_yasumizu_signature, name = "mono_CD14_yasumizu module score")
mono_CD16_yasumizu_signature <- list(c("FCGR3A","CD52","IFITM2","LST1","TCF7L2"))
th.total.uncut <- AddModuleScore(th.total.uncut, features = mono_CD16_yasumizu_signature, name = "mono_CD16_yasumizu module score")
Mac_yasumizu_signature <- list(c("MAFB","C1QC","C1QA","C1QB","A2M"))
th.total.uncut <- AddModuleScore(th.total.uncut, features = Mac_yasumizu_signature, name = "Mac_yasumizu module score")
naive_BC_yasumizu_signature <- list(c("FCER2","IL4R","IGHD","KCNG1","NT5E"))
th.total.uncut <- AddModuleScore(th.total.uncut, features = naive_BC_yasumizu_signature, name = "naive_BC_yasumizu module score")
pre_GC_BC_yasumizu_signature <- list(c("STMN1","TYMS","TCL1A","PCLAF","ZWINT"))
th.total.uncut <- AddModuleScore(th.total.uncut, features = pre_GC_BC_yasumizu_signature, name = "pre_GC_BC_yasumizu module score")
GC_BC_yasumizu_signature <- list(c("BCL6","MEF2B","BACH2","RGS13"))
th.total.uncut <- AddModuleScore(th.total.uncut, features = GC_BC_yasumizu_signature, name = "GC_BC_yasumizu module score")
PB_yasumizu_signature <- list(c("IGHG1","XBP1","MZB1","TXNDC5","FKBP11"))
th.total.uncut <- AddModuleScore(th.total.uncut, features = PB_yasumizu_signature, name = "PB_yasumizu module score")
unswi_mem_BC_yasumizu_signature <- list(c("FCRL2","SYK","LY6E","FGR"))
th.total.uncut <- AddModuleScore(th.total.uncut, features = unswi_mem_BC_yasumizu_signature, name = "unswi_mem_BC_yasumizu module score")
thym_mem_BC_yasumizu_signature <- list(c("XIST","CD27","AIM2","TEX9","IGHA1"))
th.total.uncut <- AddModuleScore(th.total.uncut, features = thym_mem_BC_yasumizu_signature, name = "thym_mem_BC_yasumizu module score")
mem_BC_I_yasumizu_signature <- list(c("JCHAIN","PLD4","ITM2C","ALOX5","CCR7"))
th.total.uncut <- AddModuleScore(th.total.uncut, features = mem_BC_I_yasumizu_signature, name = "mem_BC_I_yasumizu module score")
mem_BC_II_yasumizu_signature <- list(c("S100A10","ANXA2","ITGB1","CRIP2"))
th.total.uncut <- AddModuleScore(th.total.uncut, features = mem_BC_II_yasumizu_signature, name = "mem_BC_II_yasumizu module score")

module_yasumizu_signature <- c("EC_yasumizu module score1","tumo_assoc_fb_yasumizu module score1","norm_fb_yasumizu module score1","cTEC_yasumizu module score1","mTEC_I_yasumizu module score1","mTEC_II_yasumizu module score1","nmTEC_yasumizu module score1","NK_yasumizu module score1","NKT_thym_yasumizu module score1","NKT_peri_yasumizu module score1","gdTC_yasumizu module score1","ILC_yasumizu module score1","DN_TC_yasumizu module score1","DP_TC_yasumizu module score1","cyc_DNDP_TC_yasumizu module score1","aaCD8_I_TC_yasumizu module score1","aaCD8_II_TC_yasumizu module score1","thym_CD8_TC_yasumizu module score1","naiveTC_TET1_TC_yasumizu module score1","CD8_Tem_yasumizu module score1","CD8_Trm_yasumizu module score1","CD8_Temra_yasumizu module score1","thym_CD4_I_TC_yasumizu module score1","thym_CD4_II_TC_yasumizu module score1","CD4_Tnaive_yasumizu module score1","CD4_Tcm_Th0_yasumizu module score1","CD4_Tcm_Th2_yasumizu module score1","T_ago_yasumizu module score1","CD4_Tcm_Tfh_yasumizu module score1","CD4_Tcm_Th17_yasumizu module score1", "CD4_Tem_Th1_17_yasumizu module score1","CD4_Tem_Th1_yasumizu module score1","CD4_Temra_Th1_yasumizu module score1","naive_Treg_yasumizu module score1","act_Treg_yasumizu module score1","cDC1_yasumizu module score1","cDC2_yasumizu module score1","pDC_yasumizu module score1","mono_CD14_yasumizu module score1","mono_CD16_yasumizu module score1","Mac_yasumizu module score1","naive_BC_yasumizu module score1","pre_GC_BC_yasumizu module score1","GC_BC_yasumizu module score1","PB_yasumizu module score1", "unswi_mem_BC_yasumizu module score1","thym_mem_BC_yasumizu module score1","mem_BC_I_yasumizu module score1","mem_BC_II_yasumizu module score1")

#Bautista: upper 20
TEC1_bautista_signature <- list(c("S100A14","KRT19","KRT5","ZBED2","CTSV","PSMB11","CCL25","SLC46A2","PAX1","KREMEN2","COL17A1",
                                  "TBATA","RAB42","CALML3","KRT15","TDGF1","PRSS16","KRT17","PLEK2"))
th.total.uncut <- AddModuleScore(th.total.uncut, features = TEC1_bautista_signature, name = "TEC1_bautista module score")
TEC2_bautista_signature <- list(c("TAGLN3","MRLN","NEUROD1","NHLH1","NEB","MYOG","KLHL41","SMPX","CHGB","GNG8","KIF1A","CLDN9",
                                  "DUSP26","RTN1","HES6","CKM","KIF19","CHRNA1","PCSK1N"))
th.total.uncut <- AddModuleScore(th.total.uncut, features = TEC2_bautista_signature, name = "TEC2_bautista module score")
TEC3_bautista_signature <- list(c("PLA2G2F","A2ML1","SERPINB13","CALML5","ASPG","NCCRP1","KRT23","SERPINB3","PYDC1","SERPINB4",
                                  "IVL","DEFB1","SBSN","SPRR1B","IL1RN","SPINK5","KRT6A"))
th.total.uncut <- AddModuleScore(th.total.uncut, features = TEC3_bautista_signature, name = "TEC3_bautista module score")
MES_bautista_signature <- list(c("LUM","DPT","PID1","DCN","ABCA9","LAMA2","SFRP2","PDGFRA","GPC3","SCN7A","ABCA6","ABCA8","COL1A1","COL3A1","IGF1","COL1A2",
                                 "ASPN","OGN","OLFML3","SRPX"))
th.total.uncut <- AddModuleScore(th.total.uncut, features = MES_bautista_signature, name = "MES_bautista module score")
Peri_bautista_signature <- list(c("PLN","RERGL","TAGLN","CNN1","KCNAB1","LMOD1","CASQ2","RGS5","RRAD","MYH11","NDUFA4L2","AVPR1A","EFHD1","COX4I2","ACTA2",
                                  "HIGD1B","TBX2-AS1","LGI4","PHLDA2","CARMN"))
th.total.uncut <- AddModuleScore(th.total.uncut, features = Peri_bautista_signature, name = "Peri_bautista module score")
EC1_bautista_signature <- list(c("ACKR1","LCN6","RAMP3","NOSTRIN","PLVAP","CCL14","TPD52L1","TSPAN7","HLA-DRB5","HLA-DRB1","HYAL2","GIMAP7","PECAM1","VWF",
                                 "AQP1","FAM110D","CDC42EP3","DUSP23"))
th.total.uncut <- AddModuleScore(th.total.uncut, features = EC1_bautista_signature, name = "EC1_bautista module score")
EC2_bautista_signature <- list(c("SEMA3G","GJA5","SRGN","GJA4","PLLP","HEY1","SOX17","MECOM","SLC9A3R2","ACE","PODXL","ANXA3","ATP13A3","CXCL12","UNC5B",
                                 "ADAMTS6","ARL15","PREX1","TM4SF1","EDN1"))
th.total.uncut <- AddModuleScore(th.total.uncut, features = EC2_bautista_signature, name = "EC2_bautista module score")
EC3_bautista_signature <- list(c("APLNR","COL4A1","FABP4","COL15A1","RBP7","TM4SF18","PRND","COL4A2","SPP1","PCDH12","IL32","RGCC","KDR","PRCP",
                                 "CD34","DYSF","ROBO4","SOX7","GNG11"))
th.total.uncut <- AddModuleScore(th.total.uncut, features = EC3_bautista_signature, name = "EC3_bautista module score")
EC4_bautista_signature <- list(c("NTS","TFF3","MRC1","LYVE1","STAB2","CCL21","PROX1","TBX1","MMRN1","DPP4","MPP7","GPM6A","SEMA3D","LAPTM5","PTX3","SBSPON",
                                 "SCN3B","KLHL4","PARD6G","TSPAN5"))
th.total.uncut <- AddModuleScore(th.total.uncut, features = EC4_bautista_signature, name = "EC4_bautista module score")
RBCs_bautista_signature <- list(c("HBG2","HBA1","GYPA","GYPB","HBA2","HEMGN","HBM","AHSP","RHAG","ALAS2","TSPO2","HBG1","HBB","TRIM10","KLF1",
                                  "SLC4A1","SPTA1","GATA1","HBQ1"))
th.total.uncut <- AddModuleScore(th.total.uncut, features = RBCs_bautista_signature, name = "RBCs_bautista module score")
Immune_bautista_signature <- list(c("SMIM24","IGLL1","CD53","CD7","CD27","CD3E","RHOH","WAS","CD48","TRAT1","AC002454.1","CD79A","ITGB2","CD37","CD8A","SELL",
                                    "SASH3","TRAF3IP3","ALOX5AP","PRKCB"))
th.total.uncut <- AddModuleScore(th.total.uncut, features = Immune_bautista_signature, name = "Immune_bautista module score")
Meso_bautista_signature <- list(c("ITLN1","PRG4","MSLN","UPK3B","RSPO1","NPY","CPB1","LGALS2","TFPI2","CPA4","HHIP","TGM1","PLA2G2A","CST6","PDZK1IP1","GPM6A",
                                  "NSG1","SBSPON","KLK11","CCDC80"))
th.total.uncut <- AddModuleScore(th.total.uncut, features = Meso_bautista_signature, name = "Meso_bautista module score")

module_bautista_signature <- c("TEC1_bautista module score1","TEC2_bautista module score1","TEC3_bautista module score1","MES_bautista module score1","Peri_bautista module score1","EC1_bautista module score1","EC2_bautista module score1","EC3_bautista module score1","EC4_bautista module score1","RBCs_bautista module score1","Immune_bautista module score1","Meso_bautista module score1")

module_th.total_signature <- c(module_DEG_Xin_signature, module_DEG_yasumizu_signature,module_bautista_signature,"IGHA1", "MZB1", "XBP1")

#PC: "IGHA1", "MZB1", "XBP1" [35589692]

DotPlot(object = th.total.uncut, features = module_DEG_Xin_signature,assay="RNA") + coord_flip()+ scale_color_viridis(option = "H")+theme(axis.text.x = element_text(angle = 45)) + ggtitle("Main_modulescore") + scale_y_discrete(guide = guide_axis(angle = 90))
DotPlot(object = th.total.uncut, features = module_DEG_yasumizu_signature,assay="RNA") + coord_flip()+ scale_color_viridis(option = "H")+theme(axis.text.x = element_text(angle = 45)) + ggtitle("Main_modulescore") + scale_y_discrete(guide = guide_axis(angle = 90))
DotPlot(object = th.total.uncut, features = module_yasumizu_signature,assay="RNA") + coord_flip()+ scale_color_viridis(option = "H")+theme(axis.text.x = element_text(angle = 45)) + ggtitle("Main_modulescore") + scale_y_discrete(guide = guide_axis(angle = 90))
DotPlot(object = th.total.uncut, features = module_bautista_signature,assay="RNA") + coord_flip()+ scale_color_viridis(option = "H")+theme(axis.text.x = element_text(angle = 45)) + ggtitle("Main_modulescore") + scale_y_discrete(guide = guide_axis(angle = 90))
DotPlot(object = th.total.uncut, features =c("IGHA1", "MZB1", "XBP1"),assay="RNA") + coord_flip()+ scale_color_viridis(option = "H")+theme(axis.text.x = element_text(angle = 45)) + ggtitle("Main_modulescore") + scale_y_discrete(guide = guide_axis(angle = 90))

#Markerpanel-Yayon
DotPlot(object = th.total.uncut, features =yayon_BC,assay="RNA") + coord_flip()+ scale_color_viridis(option = "H")+theme(axis.text.x = element_text(angle = 45)) + ggtitle("Main_modulescore") + scale_y_discrete(guide = guide_axis(angle = 90))
DotPlot(object = th.total.uncut, features =yayon_DC_MAC,assay="RNA") + coord_flip()+ scale_color_viridis(option = "H")+theme(axis.text.x = element_text(angle = 45)) + ggtitle("Main_modulescore") + scale_y_discrete(guide = guide_axis(angle = 90))
DotPlot(object = th.total.uncut, features =yayon_myeloid,assay="RNA") + coord_flip()+ scale_color_viridis(option = "H")+theme(axis.text.x = element_text(angle = 45)) + ggtitle("Main_modulescore") + scale_y_discrete(guide = guide_axis(angle = 90))
DotPlot(object = th.total.uncut, features =yayon_FB,assay="RNA") + coord_flip()+ scale_color_viridis(option = "H")+theme(axis.text.x = element_text(angle = 45)) + ggtitle("Main_modulescore") + scale_y_discrete(guide = guide_axis(angle = 90))
DotPlot(object = th.total.uncut, features =yayon_EC,assay="RNA") + coord_flip()+ scale_color_viridis(option = "H")+theme(axis.text.x = element_text(angle = 45)) + ggtitle("Main_modulescore") + scale_y_discrete(guide = guide_axis(angle = 90))
DotPlot(object = th.total.uncut, features =yayon_VSMC,assay="RNA") + coord_flip()+ scale_color_viridis(option = "H")+theme(axis.text.x = element_text(angle = 45)) + ggtitle("Main_modulescore") + scale_y_discrete(guide = guide_axis(angle = 90))
DotPlot(object = th.total.uncut, features =yayon_TC,assay="RNA") + coord_flip()+ scale_color_viridis(option = "H")+theme(axis.text.x = element_text(angle = 45)) + ggtitle("Main_modulescore") + scale_y_discrete(guide = guide_axis(angle = 90))
DotPlot(object = th.total.uncut, features =yayon_TEC,assay="RNA") + coord_flip()+ scale_color_viridis(option = "H")+theme(axis.text.x = element_text(angle = 45)) + ggtitle("Main_modulescore") + scale_y_discrete(guide = guide_axis(angle = 90))

#Modulescore imaging
pdf("Graphs/Main/Modulescore_Clusteridentification.pdf", height=20, width=15)
DotPlot(object = th.total.uncut, features = module_th.total.uncut_signature,assay="RNA") + coord_flip()+ scale_color_viridis(option = "H")+theme(axis.text.x = element_text(angle = 45)) + ggtitle("Main_modulescore") + scale_y_discrete(guide = guide_axis(angle = 90))
dev.off()


#(19)_Name Cluster and contamination check####
table(th.total.uncut$celltype.full)
#"TC","BC","PC","TEC","DC", "MAC/Mono","EC","FB","VSMC/Peri", "RBCs"
th.total.uncut<-RenameIdents(th.total.uncut,
                             `0`="TC", `1`="TC", `2`="TC", `3`="TC", `4`="TC", `5`="TC", 
                             `6`="FB", `7`="EC", `8`="TC", `9`="TC", `10`="TC",
                             `11`="TEC", `12`="BC", `13`="VSMC/Peri", `14`="MAC/Mono", `15`="TEC", 
                             `16`="FB", `17`="TC", `18`="MAC/Mono", `19`="TC", `20`="PC", 
                             `21`="MAC/Mono", `22`="DC", `23`="FB", `24`="RBCs", `25`="MAC/Mono", `26`="TC", `27`="BC", `28`="EC")
th.total.uncut$celltype.full<-Idents(th.total.uncut)


#Define cluster levels
table(th.total.uncut$celltype.full)
clusters_ordered<-c("TC","BC","PC","TEC","DC", "MAC/Mono","EC","FB","VSMC/Peri", "RBCs")
th.total.uncut$celltype.full<- factor(th.total.uncut$celltype.full, levels = clusters_ordered)
Idents(th.total.uncut)<-th.total.uncut$celltype.full

p1 <- DimPlot(th.total.uncut, reduction = "umap.full", group.by = "group", alpha = 0.1)
p2 <- DimPlot(th.total.uncut, reduction = "umap.full", group.by = "celltype.full", alpha = 0.1)
#pdf("Graphs/Main/UMAP_main.pdf")
p1 + p2 + plot_layout(ncol = 1)
#dev.off()

#Set group in tissue identy 
Idents(th.total.uncut) <- th.total.uncut$group
th.total.uncut$group.tissue <- paste(Idents(th.total.uncut), th.total.uncut$tissue, sep = "_")
Idents(th.total.uncut)<-th.total.uncut$group.tissue
table(th.total.uncut$group.tissue)


#order conditions
#Name Cluster
th.total.uncut$group <- factor(x = th.total.uncut$group, levels = c("Park", "Bautista","Xin", "Yasumizu", "Direder"))
th.total.uncut$condition <- factor(x = th.total.uncut$condition, levels = c("Thymus", "TTH", "TET_A","TET_AB","TET_B","TET_C","TET_MNT"))
th.total.uncut$classification <- factor(x = th.total.uncut$classification, levels = c("Thymus", "TTH", "TET_A","TET_AB","TET_B1","TET_B2","TET_B2_B3","TET_B3","TET_C","TET_MNT"))
th.total.uncut$maintissue <- factor(x = th.total.uncut$maintissue, levels = c("Thymus","TTH","TET"))
th.total.uncut$tissue <- factor(x = th.total.uncut$tissue, levels = c("prenatal_Thymus", "pediatric_Thymus", "adult_Thymus", "TTH", "TET_A","TET_AB","TET_B","TET_C","TET_MNT"))
th.total.uncut$sex <- factor(x = th.total.uncut$sex, levels = c("female","male","n.a."))
th.total.uncut$age <- factor(x = th.total.uncut$age, levels = c("7w","8w","9w","10w","11w","12w","13w","14w","16w","17w","3m","6m","10m","15m","19m","23m","30m","10-15y","15-20y","20-25y","23y","25y","35y","35-40y","40y","44y","46y","49y","54y","55y","56y","67y","68y","69y","70y","72y","76y","77y"))
th.total.uncut$development <- factor(x = th.total.uncut$development, levels = c("embryo","fetal","pediatric","adult"))
th.total.uncut$group.tissue <- factor(x = th.total.uncut$group.tissue, levels = c("Park_prenatal_Thymus", "Park_pediatric_Thymus", "Park_adult_Thymus", "Bautista_prenatal_Thymus", "Bautista_pediatric_Thymus","Bautista_adult_Thymus","Xin_TET_A","Xin_TET_AB","Xin_TET_B","Xin_TET_C","Xin_TET_MNT","Yasumizu_TET_AB","Yasumizu_TET_B","Direder_TTH","Direder_TET_B","Direder_TET_C"))
th.total.uncut$sample <- factor(x = th.total.uncut$sample, levels = c("Thymus01","Thymus02","Thymus03","Thymus04","Thymus05","Thymus06","Thymus07","Thymus08","Thymus09","Thymus10","Thymus11","Thymus12","Thymus13","Thymus14","Thymus15","Thymus16","Thymus17","Thymus18","Thymus19","Thymus20","Thymus21","Thymus22","Thymus23","Thymus24","Thymus25","Thymus26","Thymus27","Thymus28","Thymus29","Thymus30","Thymus31",
                                                                      "TTH01","TTH02","TTH03","TTH04","TTH05","TTH06",
                                                                      "TET_A01",
                                                                      "TET_AB01","TET_AB02","TET_AB03","TET_AB04",
                                                                      "TET_B01","TET_B02","TET_B03","TET_B04","TET_B05","TET_B06","TET_B07",
                                                                      "TET_C01", "TET_C02","TET_C03",
                                                                      "TET_MNT01"
))
levels(th.total.uncut$sample)


DimPlot(th.total.uncut, reduction = "umap.full", group.by = "celltype.full", alpha = 0.4)+ ggtitle("Basic-UMAP")+ scale_color_manual(values=color_level1)

DimPlot(th.total.uncut, reduction = "umap.full",  group.by="group", split.by="tissue", ncol=3)+ ggtitle("Basic-UMAP")
DimPlot(th.total.uncut, reduction = "umap.full",  group.by="tissue", alpha = 0.1)+ ggtitle("Basic-UMAP")


DefaultAssay(th.total.uncut)<-"RNA"
Idents(th.total.uncut)<-"celltype.full"
save(th.total.uncut,file="th.total.uncut.RData")
th.total <- th.total.uncut

#Doublet identification
FeaturePlot(th.total,c("KRT19","CD3D"), blend=T, order=T, label=T, cols= c("lightgrey","green","blue"), max.cutoff = 2)
FeaturePlot(th.total,c("KRT19","CD3E"), blend=T, order=T, label=T, cols= c("lightgrey","green","blue"), max.cutoff = 2)
FeaturePlot(th.total,c("KRT19","CD3G"), blend=T, order=T, label=T, cols= c("lightgrey","green","blue"), max.cutoff = 2)

table(th.total$celltype)
th.total$celltype<-th.total$celltype.full
table(th.total$celltype)
table(th.total$celltype.full)

#(20)_Subset TEC####
Idents(th.total)<- th.total$celltype.full
TEC_sub <- subset(th.total,idents = "TEC")


DefaultAssay(TEC_sub)<-"RNA"
TEC_sub[["RNA"]]$data <- as(TEC_sub[["RNA"]]$data, Class = "dgCMatrix")

TEC_sub <- FindVariableFeatures(TEC_sub)
TEC_sub <- ScaleData(TEC_sub)
TEC_sub <- RunPCA(TEC_sub)
ElbowPlot(TEC_sub, ndims=20)
TEC_sub <- RunUMAP(TEC_sub, dims = 1:10)
TEC_sub <- FindNeighbors(TEC_sub, graph.name = "test",dims = 1:10)
TEC_sub <- FindClusters(TEC_sub, graph.name = "test",resolution = 0.1)
DimPlot(TEC_sub, label = T, label.size = 3, alpha = 1)
DimPlot(TEC_sub, label=T, split.by="condition", alpha = 1)
DimPlot(TEC_sub, label=T, split.by="development", alpha = 1)
DimPlot(TEC_sub, label=T, split.by="group", alpha = 1)
DimPlot(TEC_sub, label=T, group.by="seurat_clusters", alpha = 1)
DefaultAssay(TEC_sub)<- "RNA"

#Modulscore:
#Xin >1.5
CCL25_cTEC_Xin_signature <- list(c("CCL19","KRT15","CCL25","MT1X","CTSV","S100A14","IFI27","GPNMB","KRT5","CALML3",
                                   "KRT6A","IFI6","SFN","SORBS2","S100A9","CHCHD10","EGR1","HOPX","COX7A1","ASS1",
                                   "CXCL12","RPL37A","FBXO2","NEAT1","CCL2","PIK3R1","RPL37","IFIT3","OLFM3",
                                   "MGST3","RPL13A","NDUFA1"))
TEC_sub <- AddModuleScore(TEC_sub, features = CCL25_cTEC_Xin_signature, name = "CCL25_cTEC_Xin module score")
CHGA_TEC_Xin_signature <- list(c("TUBA1A","CHGA","NTHL1","ASCL1","TESC","KIF19","BTBD17","SINHCAF","GNG8","CHGB","DEPP1","BTG2",
                                 "CYBA","SMIM10L1","SOX4","FAIM2","CD69","IFI27"))
TEC_sub <- AddModuleScore(TEC_sub, features = CHGA_TEC_Xin_signature, name = "CHGA_TEC_Xin module score")
CHI3L1_mTEC_Xin_signature <- list(c("IGKV3-15","CHI3L1","CALML5","TSPO","LGALS3","CSF3","ELF3","CDC42EP4","SYNGR2","MSLN",
                                    "PIM3","CD24","PLAC8","RPS4Y1","PI3","MYBPC1","SELENOM","LYZ","CD9","S100A2","TFF3","GSTP1","FXYD3",
                                    "DHRS9","TNFAIP2","RARRES1","TNFRSF4","SLC9A3R1","ANXA1","AGR2","DLX5","GPX2","ATF5","PYCARD","S100A8",
                                    "ANXA4","C15orf48"))
TEC_sub <- AddModuleScore(TEC_sub, features = CHI3L1_mTEC_Xin_signature, name = "CHI3L1_mTEC_Xin module score")
GNB3_mTEC_Xin_signature <- list(c("AZGP1","PRH2","GNG13","GNB3","HEPACAM2","CXCL14","ADH1B","S100A1","ADH1C","GNAT3",
                                  "RARRES2","PLCB2","CRYM","AVIL","OVOL3","MCUB","FXYD6","TPM1","PCSK1N","TRPM5",
                                  "SMIM10L1","DCN","SH2D6","RTN1","SNHG18","CHGB","NEUROD1","CPB1","TAGLN3","RGS21","S100A13",
                                  "TMEM176B","PRH1","HLA-B","QPCT","RBP1","HAP1","PDCD5","CXCL12","HLA-C","MAP1B","PCBP4","TMEM176A",
                                  "TAS1R3","CFAP298","FMO1","HES6","SPON2","C11orf53","LCN1","GNG8","DHRS2","FBLN1","CAT",
                                  "AKR1B1","EFHD1","MUC20","STMN1","TCF4","RAMP1","NAAA","CPE","FXYD2","ANXA1","RAB26","AKAP12","GPX3",
                                  "ASCL3","ESPL1","ALKAL1","CALM2"))
TEC_sub <- AddModuleScore(TEC_sub, features = GNB3_mTEC_Xin_signature, name = "GNB3_mTEC_Xin module score")
KRT14_mTEC_Xin_signature <- list(c("CLU","SCGB3A1","CST3","CCN2","KRT14","IGFBP7","SAA1","FN1","ITM2C","SAA2","BCAM","APP","ID3",
                                   "VIM","GAS6","IGFL1","PLTP","GMDS","TGM2","CXCL1","CRABP2","MYL9","SERPINF1","SCPEP1","CCN1",
                                   "SLPI","MEG3","IGFBP5","MT1X","IGFBP6","TUSC3","XIST","KRT7","GADD45B","RARRES2","HLA-B","FBLN1",
                                   "MAOA","MYADM","FOSB","ARHGAP29","HCAR2","SPINT2"))
TEC_sub <- AddModuleScore(TEC_sub, features = KRT14_mTEC_Xin_signature, name = "KRT14_mTEC_Xin module score")
mcTEC_like_Xin_signature <- list(c("PRSS16","CCL25","PLTP","CD1E","IGFBP7","TRBC2","CD8B","CYP1B1","CD3D","TBATA","NBPF20","COL6A1",
                                   "PASK","GADD45B","C20orf204","TRBC1","CD3E","KRT14","COL6A2","POTEI","PMP22","GAS6",
                                   "MS4A1","MZB1","ELOVL4"))
TEC_sub <- AddModuleScore(TEC_sub, features = mcTEC_like_Xin_signature, name = "mcTEC_like_Xin module score")
MYOG_TEC_Xin_signature <- list(c("ACTA1","TNNT3","MYLPF","CKM","TNNC2","TPM2","MYL1","MYL5","TTN","ACTC1","TNNI1","MYOG",
                                 "SLN","TNNT1","DES","PDLIM3","TNNI2","NEB","MYH3","KLHL41","TPM1","TNNT2",
                                 "RGS16","CELA3A","CD81","RGS2","SIRT2","HSPB3","SELENOW","CYB5R1","PDE4DIP","SPARC","MYBPH","RASSF4",
                                 "CNN3","NPNT","TRIM55","MYMK","TMEM38A","ID1","BIN1","MRLN","FITM1","ANKRD36C","COX6A2","ATP2A1",
                                 "DUSP26","BLCAP","NEURL1","PTP4A3","CAV3"))
TEC_sub <- AddModuleScore(TEC_sub, features = MYOG_TEC_Xin_signature, name = "MYOG_TEC_Xin module score")

module_TEC_Xin_signature <- c("CCL25_cTEC_Xin module score1", "CHGA_TEC_Xin module score1", "CHI3L1_mTEC_Xin module score1", "GNB3_mTEC_Xin module score1", "KRT14_mTEC_Xin module score1", "mcTEC_like_Xin module score1", "MYOG_TEC_Xin module score1")

#Bautista upper 20
Immature_TEC_Bautista_signature <- list(c("SERTM1","LINC01564","CLCN1","OLAH","NNMT","IGFBP5","KRT14","ATP8B4","MAOA","LRP1B",
                                          "DPYS","KCNQ5","ACADL","ALDH1L1-AS2","MLXIPL"))
TEC_sub <- AddModuleScore(TEC_sub, features = Immature_TEC_Bautista_signature, name = "Immature_TEC_Bautista module score")
cTEC_lo_Bautista_signature <- list(c("TCF24","COL26A1","GJB6","CCL25","SCUBE1","FBN3","GALNT9","COL9A1","ENO3","UBE2C","CENPA","GLB1L3","KIF4A","CDCA3","DEPDC1",
                                     "LINC00618"))
TEC_sub <- AddModuleScore(TEC_sub, features = cTEC_lo_Bautista_signature, name = "cTEC_lo_Bautista module score")
cTEC_hi_Bautista_signature <- list(c("CFC1","TNFRSF17","CFHR1","PRSS16","RTBDN","FOXR1","PGAM2","TBATA","CLEC2L","CCL25",
                                     "PNCK","LHCGR","C4orf50","NPHS1"))
TEC_sub <- AddModuleScore(TEC_sub, features = cTEC_hi_Bautista_signature, name = "cTEC_hi_Bautista module score")
mTEC_lo_Bautista_signature <- list(c("LINC00922","CXCL9","LYPD1","LRFN2","CCL19","MROH2A","KRT15","FOXQ1","CXCL10","MMP9","GPR88","CXCL11","LINC00839","FXYD2",
                                     "IL22RA2","HTR1F","CHI3L1"))
TEC_sub <- AddModuleScore(TEC_sub, features = mTEC_lo_Bautista_signature, name = "mTEC_lo_Bautista module score")
mTEC_hi_Bautista_signature <- list(c("EDDM3B","BANF2","UMOD","AMTN","SPATA16","CSN1S1","CLDN18","DNMT3L","S100A5","CLEC4G","GIP","CD209","NT5DC4","SERPINA6",
                                     "HPD","COLEC10"))
TEC_sub <- AddModuleScore(TEC_sub, features = mTEC_hi_Bautista_signature, name = "mTEC_hi_Bautista module score")
Keratinocyte_like_Bautista_signature <- list(c("LIPK","ARHGAP40","VTCN1","TMPRSS11D","CEACAM5","SCEL","DSC1","SERPINB3","SDR9C7","SERPINB12","SERPINB13","FETUB",
                                               "PIGR","SPINK5","IVL","A2ML1","CLCA2"))
TEC_sub <- AddModuleScore(TEC_sub, features = Keratinocyte_like_Bautista_signature, name = "Keratinocyte_like_Bautista module score")
Neuroendo_Bautista_signature <- list(c("ADAD1","ATOH1","BARHL1","NEUROD6","CABP2","CLRN1","NHLH1","NEUROD1","STMN2","NEUROG1","DYDC2","SMIM18","KCNH6",
                                       "GRXCR1","TEX26","CLEC18C","CCER2"))
TEC_sub <- AddModuleScore(TEC_sub, features = Neuroendo_Bautista_signature, name = "Neuroendo_Bautista module score")
Myoid_Bautista_signature <- list(c("ASB5","DHRS7C","FGF5","ATP1B4","LMOD2","MIR1-1HG","HSPB3","GJD2","MIR133A1HG","TTN","MYH3","MYOZ2","LANCL1-AS1","CKM","COX6A2","SOHLH2","TNNI1",
                                   "ACTA1"))
TEC_sub <- AddModuleScore(TEC_sub, features = Myoid_Bautista_signature, name = "Myoid_Bautista module score")
Myelin_Bautista_signature <- list(c("TMEM215","LINC01608","CDH19","SOX10","S100B","SFRP5","PMP2","MPZ","GJB1","GJC3","C16orf82","MAG","CMTM5","GLDN","GRIK3","KCNJ10","TMPRSS5",
                                    "CNTF"))
TEC_sub <- AddModuleScore(TEC_sub, features = Myelin_Bautista_signature, name = "Myelin_Bautista module score")

module_TEC_Bautista_signature <- c("Immature_TEC_Bautista module score1","cTEC_lo_Bautista module score1","cTEC_hi_Bautista module score1","mTEC_lo_Bautista module score1","mTEC_hi_Bautista module score1","Keratinocyte_like_Bautista module score1","Neuroendo_Bautista module score1","Myoid_Bautista module score1","Myelin_Bautista module score1")


#Bautista_subset upper 20
Immature_TECI_Bautista_subset_signature <- list(c("NFKBID","RPS4Y1","EIF1AY","IRF1","EGR2","TK1","SOCS3","TRIB1","CCL25","KRT17","PMAIP1","CCL2","HLA-DQA2","CXCL14","MAFF","TMC8","PSMB11","MYADM","ATF3"))
TEC_sub <- AddModuleScore(TEC_sub, features = Immature_TECI_Bautista_subset_signature, name = "Immature_TECI_Bautista_subset module score")
Immature_TECII_Bautista_subset_signature <- list(c("IGFBP5","MYBPC1","MT1X","NPTX1","IGKC","NNMT","SAA1","GAL","IGLC3","IGLC2","FAM107A","LRP1B","IGHG3","MAOA","CORO6","HSD17B2","TXNIP","SLC16A11","DPYS","ADHFE1","XIST"))
TEC_sub <- AddModuleScore(TEC_sub, features = Immature_TECII_Bautista_subset_signature, name = "Immature_TECII_Bautista_subset module score")
mTEC_lo_Bautista_subset_signature <- list(c("KRT15","CCL19","SAA1","KRT14","CH25H","CCDC3","LAMB3","CTSV","CXCL14","HCAR2","CLU","C1S","KRT5",
                                            "DCN","ANGPT2","LAYN","KRT19","GADD45B","TGFBI","KRT17"))
TEC_sub <- AddModuleScore(TEC_sub, features = mTEC_lo_Bautista_subset_signature, name = "mTEC_lo_Bautista_subset module score")
mTEC_hi_Bautista_subset_signature <- list(c("FCER2","AIRE","FABP6","MARCO","CFP","TNFRSF9","DEFB4A","C17orf99","IL4I1","CD70","MT3","NPB","NCF1","PLEK","JCHAIN",
                                            "FEZF2","CCL27","SNX8","CLEC7A","SPIB","SLC15A3"))
TEC_sub <- AddModuleScore(TEC_sub, features = mTEC_hi_Bautista_subset_signature, name = "mTEC_hi_Bautista_subset module score")
Keratinocyte_like_Bautista_subset_signature <- list(c("SCEL","IVL","SPINK5","SERPINB3","SERPINB13","KRT23","LYPD2","PSORS1C2","A2ML1","SPRR1B",
                                                      "LCN2","PI3","SERPINB4","LY6G6C","KRT1","CYSRT1","CST6","PLA2G2F","SBSN","ASPG","CXCL17"))
TEC_sub <- AddModuleScore(TEC_sub, features = Keratinocyte_like_Bautista_subset_signature, name = "Keratinocyte_like_Bautista_subset module score")
Tuft_ionocyte_Bautista_subset_signature <- list(c("GFRA4","TRPM2","FOXI1","ZACN","MUC6","NOL4","CAMK2B","RNF186","GHRH","TMEM61","KCNH2","SPON2","USHBP1",
                                                  "NWD1","NEURL1","PROX1","PCSK1N","TOX3","KLK12","NKAIN4","ESPL1","SCG5"))
TEC_sub <- AddModuleScore(TEC_sub, features = Tuft_ionocyte_Bautista_subset_signature, name = "Tuft_ionocyte_Bautista_subset module score")
Neuroendo_Bautista_subset_signature <- list(c("NEUROD1","STMN2","GAP43","ARHGAP36","RTN1","HIGD1B","GNG8","NHLH1","GNG3","FABP3","POU4F1","PPP1R17","KLHL35",
                                              "KRT81","NEUROG1","SHOX2","PTX3","JPH4","GKAP1","PCSK1N","SOX11","STX1A"))
TEC_sub <- AddModuleScore(TEC_sub, features = Neuroendo_Bautista_subset_signature, name = "Neuroendo_Bautista_subset module score")
ciliated_Bautista_subset_signature <- list(c("USH2A","GRXCR1","BARHL1","CABP2","ATOH1","LHX3","S100A1","STRC","PCP4","CLEC18A",
                                             "CLEC18B","CCER2","KCNH6","CFAP77","CRIP3","CFAP126","C1orf194","FAM183A","PAX2","PKIB"))
TEC_sub <- AddModuleScore(TEC_sub, features = ciliated_Bautista_subset_signature, name = "ciliated_Bautista_subset module score")
Myoid_Bautista_subset_signature <- list(c("HSPB3","TTN","MYH3","COX6A2","SCT","CKM","ACTA1","MYL4","TNNI1","DES","UNC45B","SLN","ACTC1",
                                          "MYOD1","TNNT2","STAC3","FITM1","MYLPF","PITX2","KLHL41","MYOG"))
TEC_sub <- AddModuleScore(TEC_sub, features = Myoid_Bautista_subset_signature, name = "Myoid_Bautista_subset module score")
Myelin_Bautista_subset_signature <- list(c("CDH19","SOX10","S100B","SFRP5","PMP2","MAG","GJC3","MPZ","PLP1","GLDN","CMTM5","PMP22","GFRA3","RELN","LGI4","PRX","TTYH1","APOA1",
                                           "AZGP1","MLIP","AATK","NKAIN3"))
TEC_sub <- AddModuleScore(TEC_sub, features = Myelin_Bautista_subset_signature, name = "Myelin_Bautista_subset module score")
mcTEC_Yayon_signature <- list(c("DLK2","IGFBP5","IGFBP6","CCN2","CCL2","KRT15","ITGA6","MKI67"))
TEC_sub <- AddModuleScore(TEC_sub, features = mcTEC_Yayon_signature, name = "mcTEC_Yayon module score")


module_TEC_Bautista_subset_signature <- c("Immature_TECI_Bautista_subset module score1","Immature_TECII_Bautista_subset module score1","mTEC_lo_Bautista_subset module score1","mTEC_hi_Bautista_subset module score1",
                                          "Keratinocyte_like_Bautista_subset module score1","Tuft_ionocyte_Bautista_subset module score1","Neuroendo_Bautista_subset module score1","ciliated_Bautista_subset module score1",
                                          "Myoid_Bautista_subset module score1","Myelin_Bautista_subset module score1")


module_TEC_DEG_yasumizu_signature <- c("cTEC_DEG_yasumizu module score1","mTEC_I_DEG_yasumizu module score1","mTEC_II_DEG_yasumizu module score1","nmTEC_DEG_yasumizu module score1")

TEC_module <- c(module_TEC_DEG_yasumizu_signature, module_TEC_Bautista_signature, module_TEC_Bautista_subset_signature, module_TEC_Xin_signature)

DotPlot(object = TEC_sub, features = module_TEC_Xin_signature,assay="RNA") + coord_flip()+ scale_color_viridis(option = "H")+theme(axis.text.x = element_text(angle = 45)) + ggtitle("TEC_modulescore") + scale_y_discrete(guide = guide_axis(angle = 90))
DotPlot(object = TEC_sub, features = module_TEC_Bautista_signature,assay="RNA") + coord_flip()+ scale_color_viridis(option = "H")+theme(axis.text.x = element_text(angle = 45)) + ggtitle("TEC_modulescore") + scale_y_discrete(guide = guide_axis(angle = 90))
DotPlot(object = TEC_sub, features = module_TEC_DEG_yasumizu_signature,assay="RNA") + coord_flip()+ scale_color_viridis(option = "H")+theme(axis.text.x = element_text(angle = 45)) + ggtitle("TEC_modulescore") + scale_y_discrete(guide = guide_axis(angle = 90))
DotPlot(object = TEC_sub, features = module_TEC_Bautista_subset_signature,assay="RNA") + coord_flip()+ scale_color_viridis(option = "H")+theme(axis.text.x = element_text(angle = 45)) + ggtitle("TEC_modulescore") + scale_y_discrete(guide = guide_axis(angle = 90))
DotPlot(object = TEC_sub, features = module_DEG_yasumizu_signature,assay="RNA") + coord_flip()+ scale_color_viridis(option = "H")+theme(axis.text.x = element_text(angle = 45)) + ggtitle("TEC_modulescore") + scale_y_discrete(guide = guide_axis(angle = 90))
DotPlot(object = TEC_sub, features = yayon_TEC,assay="RNA") + coord_flip()+ scale_color_viridis(option = "H")+theme(axis.text.x = element_text(angle = 45)) + ggtitle("TEC_modulescore") + scale_y_discrete(guide = guide_axis(angle = 90))

FeaturePlot(TEC_sub,"TEC_Xin module score1", order=T)

#Park: cTEC:"PSMB11", mcTEC: "DLK2","KRT14", mTEC(I):"KRT14","CCL19",mTEC(II): "CCL19","AIRE", mTEC(III):"KRT1","TEC(myo):"MYOD1", TEC(neuro): NEUROG1
DotPlot(TEC_sub,features=c("PSMB11", "DLK2","KRT14", "CCL19","AIRE", "KRT1","MYOD1", "NEUROG1")) + scale_color_viridis(option = "H")+coord_flip()#+ theme(panel.background = element_rect(fill = 'lightgrey'))
#Yasumiz: cTEC: "CCL25","PSMB11", mTEC(I):"KRT15","IFI27",mTEC(II): "CLDN4","SDC4", nmTEC: "GABRA5","NEFL"
DotPlot(TEC_sub,features=c("CCL25","PSMB11", "KRT15","IFI27","CLDN4","SDC4", "GABRA5","NEFL")) + scale_color_viridis(option = "H")+coord_flip()#+ theme(panel.background = element_rect(fill = 'lightgrey'))

# Bautista: cTEC: "FOXN1","PAX9","SIX1", "PSMB11","HLA-DQB1","PRSS16","CCL25", Immature TEC: "CCL25",Aire+mTEC hi: "HLA-DQB1","CLDN4","SPIB","AIRE","FEZF2",Corneo-like mTEC: "KRT1","Neuroendocrine: "BEX1","NEUROD1","Myoid: "MYOD1","DES",Myelin+: SOX10","MPZ"
DotPlot(TEC_sub, features= c("FOXN1","PAX9","SIX1", "PSMB11","HLA-DQB1","PRSS16","CCL25", "CLDN4","SPIB","AIRE","FEZF2","KRT1","BEX1","NEUROD1","MYOD1","DES","SOX10","MPZ"))+coord_flip()
FeaturePlot(TEC_sub,"SPIB", order=T, label=T)

#Name Cluster
Idents(TEC_sub)<-TEC_sub$seurat_clusters
TEC_sub<-RenameIdents(TEC_sub,
                      `0`="mTEC_lo", `1`="cTEC", `2`="Carc_TEC_1", `3`="TEC_myo",`4`="mixed_TEC",`5`="Carc_TEC_2",
                      `6`="TET_TEC",`7`="DP",`8`="mTEC_hi")
TEC_sub$celltype<-Idents(TEC_sub)

levels(TEC_sub)

#Define cluster levels
clusters_ordered<-c("cTEC", "mTEC_lo","mTEC_hi", "mixed_TEC","TEC_myo", "TET_TEC", "Carc_TEC_1", "Carc_TEC_2","DP")
TEC_sub$celltype<- factor(TEC_sub$celltype, levels = clusters_ordered)
Idents(TEC_sub)<-TEC_sub$celltype
TEC_sub_safe <- TEC_sub

#Doublet and contamination check
DimPlot(TEC_sub, label = TRUE, label.size = 6)
FeaturePlot(TEC_sub,c("KRT19","CD3D"), blend=T, order=T, label=T, cols= c("lightgrey","green","blue"), max.cutoff = 2)

cells.of.interest_1 <- WhichCells(
  TEC_sub,
  idents = "DP"
)

TEC_sub <- subset(TEC_sub,idents = c("cTEC", "mTEC_lo","mTEC_hi", "mixed_TEC","TEC_myo", "TET_TEC", "Carc_TEC_1", "Carc_TEC_2"))

table(TEC_sub$celltype)
table(TEC_sub$celltype.full)


DefaultAssay(TEC_sub)<-"RNA"
Idents(TEC_sub)<-"celltype"

#subcluster_mixed_TEC
TEC_sub <- FindSubCluster(TEC_sub, "mixed_TEC", "test", subcluster.name = "subcelltype",  resolution = 0.2, algorithm = 1)
DimPlot(TEC_sub,  group.by = "subcelltype", label = TRUE, label.size = 6)

DotPlot(object = TEC_sub, features = TEC_module,assay="RNA") + coord_flip()+ scale_color_viridis(option = "H")+theme(axis.text.x = element_text(angle = 45)) + ggtitle("TEC_modulescore") + scale_y_discrete(guide = guide_axis(angle = 90))
DotPlot(object = TEC_sub, features = yayon_TEC,assay="RNA") + coord_flip()+ scale_color_viridis(option = "H")+theme(axis.text.x = element_text(angle = 45)) + ggtitle("TEC_modulescore") + scale_y_discrete(guide = guide_axis(angle = 90))

#Name Cluster
Idents(TEC_sub)<-TEC_sub$subcelltype
table(TEC_sub$subcelltype)
TEC_sub<-RenameIdents(TEC_sub,
                      `mTEC_lo`="mTEC_lo", `Carc_TEC_1`="Carc_TEC_1", `cTEC`="cTEC", `TEC_myo`="TEC_myo_1",
                      `mixed_TEC_0`="TEC_tuft", `mixed_TEC_1`="TEC_neuro_ciliated", `mixed_TEC_2`="TEC_tuft", `mixed_TEC_3`="TEC_myo_2",
                      `Carc_TEC_2`="Carc_TEC_2",`TET_TEC`="TET_TEC",`TC`="TC",`mTEC_hi`="mTEC_hi")
TEC_sub$subcelltype<-Idents(TEC_sub)

TEC_sub$celltype<-TEC_sub$subcelltype

#identification mcTEC
TEC_sub <- FindSubCluster(TEC_sub, "mTEC_lo", "test", subcluster.name = "subcelltype",  resolution = 0.1, algorithm = 1)
DimPlot(TEC_sub,  group.by = "subcelltype", label = TRUE, label.size = 6)

DotPlot(object = TEC_sub,  group.by = "subcelltype", features = c(TEC_module,"mcTEC_Yayon module score1"),assay="RNA") + coord_flip()+ scale_color_viridis(option = "H")+theme(axis.text.x = element_text(angle = 45)) + ggtitle("TEC_modulescore") + scale_y_discrete(guide = guide_axis(angle = 90))
DotPlot(object = TEC_sub,  group.by = "subcelltype", features = yayon_TEC,assay="RNA") + coord_flip()+ scale_color_viridis(option = "H")+theme(axis.text.x = element_text(angle = 45)) + ggtitle("TEC_modulescore") + scale_y_discrete(guide = guide_axis(angle = 90))
DotPlot(object = TEC_sub,  group.by = "subcelltype", features = yayon_TEC,assay="RNA") + coord_flip()+ scale_color_viridis(option = "H")+theme(axis.text.x = element_text(angle = 45)) + ggtitle("TEC_modulescore") + scale_y_discrete(guide = guide_axis(angle = 90))

Idents(TEC_sub)<-TEC_sub$subcelltype
table(TEC_sub$subcelltype)
TEC_sub<-RenameIdents(TEC_sub,
                      `mTEC_lo_0`="mTEC_lo", `mTEC_lo_1`="mTEC_lo", `mTEC_lo_2`="mTEC_lo", `mTEC_lo_3`="mcTEC", `Carc_TEC_1`="Carc_TEC_1", `cTEC`="cTEC", `TEC_myo_1`="TEC_myo_1",
                      `TEC_tuft`="TEC_tuft", `TEC_neuro_ciliated`="TEC_neuro_ciliated", `TEC_tuft`="TEC_tuft", `TEC_myo_2`="TEC_myo_2",
                      `Carc_TEC_2`="Carc_TEC_2",`TET_TEC`="TET_TEC",`TC`="TC",`mTEC_hi`="mTEC_hi")
TEC_sub$subcelltype<-Idents(TEC_sub)


#Define cluster levels
table(TEC_sub$subcelltype)
Idents(TEC_sub)<-TEC_sub$subcelltype
clusters_ordered<-c("cTEC", "mcTEC","mTEC_lo","mTEC_hi","TEC_tuft", "TEC_neuro_ciliated","TEC_myo_1","TEC_myo_2", "TET_TEC", "Carc_TEC_1", "Carc_TEC_2","TC")
TEC_sub$subcelltype<- factor(TEC_sub$subcelltype, levels = clusters_ordered)
Idents(TEC_sub)<-TEC_sub$subcelltype
table(TEC_sub$subcelltype)
TEC_sub$celltype<-TEC_sub$subcelltype
TEC_sub$celltype.full<-TEC_sub$subcelltype
levels(TEC_sub)
TEC_sub_uncut<-TEC_sub

#Modulescore imaging
pdf("Graphs/Main/Modulescore_TEC_Clusteridentification.pdf", width=15, height=7)
DotPlot(object = TEC_sub, features = c(TEC_module, yayon_TEC),assay="RNA")+ scale_color_viridis(option = "H")+theme(axis.text.x = element_text(angle = 90)) + ggtitle("TEC_modulescore") 
dev.off()

pdf("Graphs/Main/Modulescore_TEC_Clusteridentification_yayon.pdf", width=9, height=9)
DotPlot(object = TEC_sub, features = yayon_TEC,assay="RNA") + coord_flip()+ scale_color_viridis(option = "H")+theme(axis.text.x = element_text(angle = 45)) + ggtitle("TEC_modulescore") + scale_y_discrete(guide = guide_axis(angle = 90))
dev.off()
pdf("Graphs/Main/TEC_UMAPPlots.pdf", width=5, height=5)
DimPlot(TEC_sub, label=T, label.color = "black",label.size = 5, repel=T, alpha=1) + ggtitle("TEC-subset")+ scale_color_manual(values=color_levelTEC) + NoLegend()
dev.off()

pdf("Graphs/Main/TEC_split_celltype_UMAPPlots.pdf", width=30, height=5)
DimPlot(TEC_sub, label=T, label.color = "black",label.size = 5, repel=T, alpha=1, split.by="celltype.full") + ggtitle("TEC-subset")+ scale_color_manual(values=color_levelTEC) + NoLegend()
dev.off()


pdf("Graphs/Main/UMAPPlots_TEC_sub_QC.pdf", width=5, height=5)
DimPlot(TEC_sub, group.by="group")+ ggtitle("Basic-UMAP")
DimPlot(TEC_sub, group.by="tissue")+ ggtitle("Basic-UMAP")
dev.off()

pdf("Graphs/Main/TEC_sub_APTEC_Clusteridentification.pdf", width=10, height=5)
FeaturePlot(TEC_sub,c("KRT19","HLA-DRA"), blend=T, order=T, cols= c("lightgrey","red","blue"))
dev.off()

#Cell_cluster_count
Idents(TEC_sub) <- TEC_sub$subcelltype
TEC_sub$celltype.tissue <- paste(Idents(TEC_sub), TEC_sub$tissue, sep = "_")
Idents(TEC_sub)<-TEC_sub$celltype.tissue
levels(TEC_sub)
Idents(TEC_sub) <- TEC_sub$celltype.tissue
TEC_sub$celltype.tissue.donor <- paste(Idents(TEC_sub), TEC_sub$donor, sep = "_")
Idents(TEC_sub)<-TEC_sub$celltype.tissue.donor
levels(TEC_sub)
table(TEC_sub$celltype.tissue.donor)
write.xlsx(table(TEC_sub$celltype.tissue.donor), "Lists/TEC_sub_ctd.xlsx")

pdf("Graphs/Main/TEC_Barplots.pdf")
plot_integrated_clusters_group (TEC_sub)
plot_integrated_clusters_Donor_ID (TEC_sub)
plot_integrated_celltypes_condition (TEC_sub)
plot_integrated_celltypes_classification (TEC_sub)
plot_integrated_celltypes_tissue (TEC_sub)
plot_integrated_celltypes_sex (TEC_sub)
plot_integrated_celltypes_age (TEC_sub) 
plot_integrated_celltypes_development (TEC_sub)
plot_integrated_celltypes_donor (TEC_sub)
plot_TEC_integrated_condition_celltypes(TEC_sub)
plot_TEC_integrated_development_celltypes(TEC_sub)
plot_TEC_integrated_group_clusters(TEC_sub)
plot_TEC_integrated_tissue_celltypes(TEC_sub)
dev.off()



#add identity to analyze groups
#Set group in identy 
Idents(TEC_sub) <- TEC_sub$celltype.full
TEC_sub$celltype.group <- paste(Idents(TEC_sub), TEC_sub$group, sep = "_")
Idents(TEC_sub)<-TEC_sub$celltype.group
table(TEC_sub$celltype.group)

Idents(TEC_sub)<- TEC_sub$celltype.group
levels(TEC_sub)
clusters_ordered = c("cTEC_Park", "cTEC_Xin", "cTEC_Yasumizu", "cTEC_Bautista", "cTEC_Direder",
                     "mcTEC_Park", "mcTEC_Bautista", "mcTEC_Direder",
                     "mTEC_lo_Park", "mTEC_lo_Xin","mTEC_lo_Yasumizu","mTEC_lo_Bautista","mTEC_lo_Direder",
                     "mTEC_hi_Park", "mTEC_hi_Xin", "mTEC_hi_Yasumizu", "mTEC_hi_Bautista","mTEC_hi_Direder",
                     "TEC_tuft_Park", "TEC_tuft_Xin", "TEC_tuft_Yasumizu", "TEC_tuft_Bautista", "TEC_tuft_Direder",       
                     "TEC_neuro_ciliated_Park", "TEC_neuro_ciliated_Bautista",         
                     "TEC_myo_1_Park", "TEC_myo_1_Bautista", "TEC_myo_1_Xin", "TEC_myo_1_Direder",                                   
                     "TEC_myo_2_Park", "TEC_myo_2_Bautista", 
                     "TET_TEC_Park", "TET_TEC_Xin", "TET_TEC_Yasumizu", "TET_TEC_Bautista", "TET_TEC_Direder",
                     "Carc_TEC_1_Direder", 
                     "Carc_TEC_2_Xin","Carc_TEC_2_Direder")
TEC_sub$celltype.group<- factor(TEC_sub$celltype.group, levels = clusters_ordered)

#Modulescore imaging
pdf("Graphs/Main/Modulescore_TEC_Clusteridentification_groupspecific.pdf", width=12, height=15)
DotPlot(object = TEC_sub, features = TEC_module,assay="RNA") + scale_color_viridis(option = "H") + ggtitle("TEC_modulescore") + scale_x_discrete(guide = guide_axis(angle = 90))
dev.off()

#add identity to analyze tissue/health diffgenes
#Set tissue in identy 
Idents(TEC_sub) <- TEC_sub$celltype.full
TEC_sub$celltype.condition <- paste(Idents(TEC_sub), TEC_sub$condition, sep = "_")
Idents(TEC_sub)<-TEC_sub$celltype.condition
table(TEC_sub$celltype.condition)
#Set tissue in identy 
Idents(TEC_sub) <- TEC_sub$celltype.condition
TEC_sub$celltype.condition.development <- paste(Idents(TEC_sub), TEC_sub$development, sep = "_")
Idents(TEC_sub)<-TEC_sub$celltype.condition.development
table(TEC_sub$celltype.condition.development)

#Clustermarker
Idents(TEC_sub)<-"celltype.full"
TEC_sub.clustermarker<-FindAllMarkers(TEC_sub, assay = "RNA", min.pct = 0.25, logfc.threshold = 0.25)
TEC_sub.clustermarker$Foldchange_UP <- 2^(TEC_sub.clustermarker$avg_log2FC)
TEC_sub.clustermarker$Foldchange_DOWN <- 2^(-TEC_sub.clustermarker$avg_log2FC)
TEC_sub.clustermarker$Ratio_pct1_pct2 <- (TEC_sub.clustermarker$pct.1)/(TEC_sub.clustermarker$pct.2)
write.xlsx(TEC_sub.clustermarker, "Lists/Clustermarker_TEC_sub_clustermarker.xlsx")

#DotPlot Top up and down regulated genes 
#up
TEC_sub.clustermarker_UP <- subset(TEC_sub.clustermarker, Foldchange_UP >=2 )
top10_UP_clustermarker <- TEC_sub.clustermarker_UP %>% group_by(cluster) %>% top_n(n=10,wt=Foldchange_UP)
top10_UP_clustermarker<-arrange(top10_UP_clustermarker,cluster,Foldchange_UP)
top10_UP_clustermarker <- as.character(top10_UP_clustermarker$gene)
top10_UP_clustermarker<-top10_UP_clustermarker[!duplicated(top10_UP_clustermarker)]

pdf("Graphs/Main/DotPlot_Top_10UP_Clustermarker_TEC_FCmin2.pdf", width = 7, height = 20)
DotPlot(TEC_sub, features=top10_UP_clustermarker, group.by = "celltype.full",assay="RNA") + coord_flip()+ scale_color_viridis(option = "H")+theme(axis.text.x = element_text(angle = 45)) + ggtitle("TEC Cluster_DEG") + scale_y_discrete(guide = guide_axis(angle = 90))
dev.off()

pdf("Graphs/Main/DotPlot_Clustermarker_TET_TEC_FCmin2.pdf", width = 8, height = 8)
DotPlot(TEC_sub, features=c("ERBB4","CCDC144NL-AS1","CASC9","PDGFRL","KHDRBS2","GNGT1","CFH","AUTS2","DLEU1","MIR100HG","ADGRV1","PLCB1","ROBO2","RARRES2","CNTN4",
                            "SCPEP1","FTX","GPC4","FARP1","GAS6","IMMP2L","ARHGAP29","IGFBP7","PDZD2","RBMS3","TGM2","APBB2","ZBTB20","FNDC3B","MEIS1","RORA","NEO1",
                            "PRSS16","GMDS","TBC1D5","KIAA1217","EYA1","MAOA"), group.by = "celltype.full",assay="RNA") + coord_flip()+ scale_color_viridis(option = "H")+theme(axis.text.x = element_text(angle = 45)) + ggtitle("TEC Cluster_DEG") + scale_y_discrete(guide = guide_axis(angle = 90))
dev.off()


#Heatmap#
## Define Gene lists that you want to plot as Heatmap
TEC_sub.clustermarker_UP<-arrange(TEC_sub.clustermarker_UP,cluster,Foldchange_UP)
genes <- as.character(TEC_sub.clustermarker_UP$gene)
avgexp <- TEC_sub

## Plot Heatmaps
## Calculate Average gene expressions (pseudobulked) for each of the celltypes in each condition
levels(avgexp)
DefaultAssay(avgexp) <- "RNA"
avgexp <- NormalizeData(avgexp)
avgexp <- AggregateExpression(avgexp, return.seurat =T, assays = "RNA")


# Modify the heatmap plot
library(RColorBrewer)
heatmap_plot <- DoHeatmap(
  avgexp,
  features = genes,
  size = 6,
  group.bar.height = 0.02,
  draw.lines = FALSE,angle = 90,
  group.colors = color_levelTEC
) +
  scale_fill_gradientn(colours = rev(brewer.pal(8, name = "RdYlBu")))

# Modify the theme to make rownames cursive
heatmap_plot 
ggsave("Graphs/Main/Heatmap_TEC_clustermarker.pdf", width = 40, height = 20, units = c("cm"))


#Carc-TEC Heatmap TOP 10 genes both cluster#
## Define Gene lists that you want to plot as Heatmap

Idents(TEC_sub)<- TEC_sub$celltype.full
genes <- as.character(c("FCN2","SOX18","CCL4","PITX1","CCL20","FXYD2","PARD6G-AS1","ATP2A3","TAMALIN","DUSP2","CD177","CA9",
                        "KRTAP19-1","MUC4","MUC1","MUC16","CEL","LGALS7","CXCL14","PLAU"))
avgexp <- TEC_sub

## Plot Heatmaps
## Calculate Average gene expressions (pseudobulked) for each of the celltypes in each condition
levels(avgexp)
DefaultAssay(avgexp) <- "RNA"
avgexp <- NormalizeData(avgexp)
avgexp <- AggregateExpression(avgexp, return.seurat =T, assays = "RNA")


# Modify the heatmap plot
library(RColorBrewer)
heatmap_plot <- DoHeatmap(
  avgexp,
  features = genes,
  size = 6,
  group.bar.height = 0.02,
  draw.lines = FALSE,angle = 90,
  group.colors = color_levelTEC
) +
  scale_fill_gradientn(colours = rev(brewer.pal(8, name = "RdYlBu")))

# Modify the theme to make rownames cursive
heatmap_plot 
ggsave("Graphs/Main/Heatmap_TECsub_Carc.pdf", width = 10, height = 10, units = c("cm"))


pdf("Graphs/Main/Featureblend_TET_A_APP_CD74.pdf", width=10, height=5)
FeaturePlot(TET_A_pathway,c("APP","CD74"), blend=T, order=T, label=T, cols= c("lightgrey","green","blue"), max.cutoff = 2)
dev.off()


#Feature Marker TETEC
pdf("Graphs/Main/Feature_KRT14_TET_TEC.pdf", width=5, height=5)
FeaturePlot(TEC_sub, c("KRT14"), order=T, label=T,pt.size=0.01, max.cutoff = 3)
dev.off()

save(TEC_sub,file="TEC_sub.RData")

#(21)_Subset TC#####


Idents(th.total)<- th.total$celltype.full
TC_sub <- subset(th.total,idents = "TC")

DefaultAssay(TC_sub)<-"RNA"
TC_sub[["RNA"]]$data <- as(TC_sub[["RNA"]]$data, Class = "dgCMatrix")

TC_sub <- FindVariableFeatures(TC_sub)
TC_sub <- ScaleData(TC_sub)
TC_sub <- RunPCA(TC_sub)
ElbowPlot(TC_sub, ndims=20)
TC_sub <- RunUMAP(TC_sub, dims = 1:10)
TC_sub <- FindNeighbors(TC_sub, dims = 1:10)
TC_sub <- FindClusters(TC_sub, resolution = 0.1)
DimPlot(TC_sub, label = T, label.size = 3, alpha = 0.1)
DimPlot(TC_sub, label=T, split.by="condition", alpha = 1)
DimPlot(TC_sub, label=T, split.by="development", alpha = 1)
DimPlot(TC_sub, label=T, group.by="seurat_clusters", alpha = 1)
DefaultAssay(TC_sub)<- "RNA"

#Modulscore:
#XIN logfc >1
DN_Xin_signature <- list(c("IGKV3-15","IGKC","IGHG1","IGKV4-1","IGKV1-5","HLA-DRA","CD74","IGKV2-30","MS4A1","IGHG2","HLA-DRB5",
                           "RRM2","HLA-DRB1","HLA-DPA1","IGHV4-39","TYMS","TUBA1B","HLA-DQA1","HMGB2","UBE2C","HLA-DPB1","TUBB","TK1","IGHM",
                           "PKM","MKI67","PTTG1","IGHG4","CD79B"))
TC_sub <- AddModuleScore(TC_sub, features = DN_Xin_signature, name = "DN_Xin module score")
DP_P_Xin_signature <- list(c("TUBA1B","TYMS","HMGB2","TUBB","HMGN2","RRM2","UBE2C","STMN1","TOP2A","MKI67","CENPF","UHRF1",
                             "HMGB1","TPX2","NUSAP1","MCM7","DUT","CCNB2","GAPDH","CDC20","LMNB1","ASPM",
                             "CCNA2","SMC4","ID1","CDK1","PTTG1","PCNA","DEK","UBE2T","TMPO","TRDC","BIRC5","RANBP1","KIF22","CDCA7",
                             "NUCKS1","HMGA1","HMGB3","CCNB1","CENPU","RAN","MCM3","ANP32B","SMC2","PPIA","NEIL3","LSM4","MAD2L1","CKS2"))
TC_sub <- AddModuleScore(TC_sub, features = DP_P_Xin_signature, name = "DP_P_Xin module score")
DP_Q_Xin_signature <- list(c("RAG1","CD1E","MZB1","ELOVL4","ARPP21","CD8B","DNTT","CYP2U1","CD1B","SOX4","CD1A","AQP3","CD38","RASSF6",
                             "PTCRA","TFDP2","SH3TC1","SATB1","TCF7","PTP4A2","RAG2","GLUL","MDM4","BCL2L1","CD99","AEBP1","XBP1",
                             "TRBC2","HDAC7","TCF12","SSBP2","MIR181A1HG","LDLRAD4","CD8A","ADA","APBB1IP","EPHB6",
                             "MAL","SMPD3","NFATC3","TRBC1","PITPNM2","PSAP","CAMK4","LEF1","RORC","LAIR1","VOPP1","SCAI","SH2D1A","DAPK1","CHRNA3","GALNT7","CD3G","CD3D","MAP1A","DGKA","NEIL1","MTA3","DCK","CCR9","SLAMF1","NUCB2"))
TC_sub <- AddModuleScore(TC_sub, features = DP_Q_Xin_signature, name = "DP_Q_Xin module score")
Memory_CD4_Xin_signature <- list(c("CXCL13","COTL1","GEM","FAAH2","THADA","BATF","TIGIT","KLRB1","MAGEH1","TSHZ2","SLC2A3","RNF19A"))
TC_sub <- AddModuleScore(TC_sub, features = Memory_CD4_Xin_signature, name = "Memory_CD4__Xin module score")
Memory_CD8_Xin_signature <- list(c("GZMK","CCL5","TNFSF9","CST7","HSPA1B","GZMA","TUBA4A","HSPA1A","TNF","CCL4","EGR1","FOSB"))
TC_sub <- AddModuleScore(TC_sub, features = Memory_CD8_Xin_signature, name = "Memory_CD8_Xin module score")
CD8_Trm_Xin_signature <- list(c("CCL4","NKG7","CCL5","GZMB","CCL4L2","GZMH","GNLY","IFNG","GZMA","CCL3","CST7","GZMK","TNFSF9","RGS1",
                                "CXCL13","HLA-DRB1","APOBEC3G","DUSP4","LGALS1","HLA-DPA1","PRF1","HLA-DRB5","VCAM1","CD74","LAG3","XCL1",
                                "PHLDA1","HSPA1B","COTL1","DUSP1","HLA-DPB1","CLEC2B","ZNF683","HSPA1A","IFITM2","ATF3",
                                "IFITM1","OASL","PSMB9","ISG15","S100A6","CTSC"))
TC_sub <- AddModuleScore(TC_sub, features = CD8_Trm_Xin_signature, name = "CD8_Trm_Xin module score")
naiveTC_TET2_Xin_signature <- list(c("LMNA","ANXA1","CCR7","GPR183","FTH1"))
TC_sub <- AddModuleScore(TC_sub, features = naiveTC_TET2_Xin_signature, name = "naiveTC_TET2_Xin module score")
naiveTC_TET1_Xin_signature <- list(c("CCR7","LYPD3"))
TC_sub <- AddModuleScore(TC_sub, features = naiveTC_TET1_Xin_signature, name = "naiveTC_TET1_Xin module score")
NKT_Xin_signature <- list(c("GNLY","XCL2","XCL1","TYROBP","NKG7","CTSW","KLRB1","GZMB","CCL4","PRF1","FCER1G","KLRD1",
                            "HOPX","IFITM2","CCL5","CST7","FCGR3A","BHLHE40","KRT86","CCL3","IL4I1","FGFBP2","CD83",
                            "KLRC1","IFITM3","LMNA","MAP3K8","GZMH","GZMA","ADGRE5","NFKBIA","IL2RB","CLIC3","SPON2"))
TC_sub <- AddModuleScore(TC_sub, features = NKT_Xin_signature, name = "NKT_Xin module score")
TH_like_Xin_signature <- list(c("TPSAB1","NEAT1","MT-CYB","MT-ND4L","MT-ND5","MT-ND3","MT-CO3","LINC-PINT","MT-CO1","MT-ND4","XIST",
                                "MT-ATP6","EML4","PDE4B","SLC2A3","MT-ND1","MT-ND2","CD44","ADGRE5","IGKV2-30",
                                "KRT17","CLU","NFKBIZ","CST3","RALGAPA1","FOXP1","AHNAK","MT-ATP8","STAT4"))
TC_sub <- AddModuleScore(TC_sub, features = TH_like_Xin_signature, name = "TH_like_Xin module score")
Treg_Xin_signature <- list(c("TNFRSF4","FOXP3","TNFRSF18","DUSP4","BATF","PHLDA1","S100A4","TIGIT","IL2RA","IL32","CTSC","TNFRSF9",
                             "CARD16","SAT1","TYMP","TNFRSF1B","CTLA4","LGALS1","ICOS","GBP5","SYNGR2","ISG15","RGS1",
                             "IFI27","DNPH1","LAIR2","CXCR6","SH2D2A"))
TC_sub <- AddModuleScore(TC_sub, features = Treg_Xin_signature, name = "Treg_Xin module score")
aßT_entry_Xin_signature <- list(c("TRDV1","ITM2A","TOX2"))
TC_sub <- AddModuleScore(TC_sub, features = aßT_entry_Xin_signature, name = "aßT_entry_Xin module score")

module_TC_XIN_signature <- c("DN_Xin module score1","DP_P_Xin module score1","DP_Q_Xin module score1","Memory_CD4_Xin module score1","Memory_CD8_Xin module score1","CD8_Trm_Xin module score1","naiveTC_TET2_Xin module score1","naiveTC_TET1_Xin module score1","NKT_Xin module score1","TH_like_Xin module score1","Treg_Xin module score1","aßT_entry_Xin module score1")
module_TC_yasumizu_signature <- c("gdTC_yasumizu module score1","ILC_yasumizu module score1","DN_TC_yasumizu module score1","DP_TC_yasumizu module score1","cyc_DNDP_TC_yasumizu module score1","aaCD8_I_TC_yasumizu module score1","aaCD8_II_TC_yasumizu module score1","thym_CD8_TC_yasumizu module score1","naiveTC_TET1_TC_yasumizu module score1","CD8_Tem_yasumizu module score1","CD8_Trm_yasumizu module score1","CD8_Temra_yasumizu module score1","thym_CD4_I_TC_yasumizu module score1","thym_CD4_II_TC_yasumizu module score1","CD4_Tnaive_yasumizu module score1","CD4_Tcm_Th0_yasumizu module score1","CD4_Tcm_Th2_yasumizu module score1","T_ago_yasumizu module score1","CD4_Tcm_Tfh_yasumizu module score1","CD4_Tcm_Th17_yasumizu module score1", "CD4_Tem_Th1_17_yasumizu module score1","CD4_Tem_Th1_yasumizu module score1","CD4_Temra_Th1_yasumizu module score1","naive_Treg_yasumizu module score1","act_Treg_yasumizu module score1")
module_TC_DEG_yasumizu_signature <- c("gdTC_DEG_yasumizu module score1","ILC_DEG_yasumizu module score1","DN_TC_DEG_yasumizu module score1","DP_TC_DEG_yasumizu module score1","cyc_DNDP_TC_DEG_yasumizu module score1","aaCD8_I_TC_DEG_yasumizu module score1","aaCD8_II_TC_DEG_yasumizu module score1","thym_CD8_TC_DEG_yasumizu module score1","naiveTC_TET1_TC_DEG_yasumizu module score1","CD8_Tem_DEG_yasumizu module score1","CD8_Trm_DEG_yasumizu module score1","CD8_Temra_DEG_yasumizu module score1","thym_CD4_I_TC_DEG_yasumizu module score1","thym_CD4_II_TC_DEG_yasumizu module score1","CD4_Tnaive_DEG_yasumizu module score1","CD4_Tcm_Th0_DEG_yasumizu module score1","CD4_Tcm_Th2_DEG_yasumizu module score1","T_ago_DEG_yasumizu module score1","CD4_Tcm_Tfh_DEG_yasumizu module score1","CD4_Tcm_Th17_DEG_yasumizu module score1", "CD4_Tem_Th1_17_DEG_yasumizu module score1","CD4_Tem_Th1_DEG_yasumizu module score1","CD4_Temra_Th1_DEG_yasumizu module score1","naive_Treg_DEG_yasumizu module score1","act_Treg_DEG_yasumizu module score1")
TC_module_list <- c(module_TC_DEG_yasumizu_signature,module_TC_XIN_signature)

DotPlot(object = TC_sub, features = module_TC_XIN_signature,assay="RNA") + coord_flip()+ scale_color_viridis(option = "H")+theme(axis.text.x = element_text(angle = 45)) + ggtitle("TC_modulescore") + scale_y_discrete(guide = guide_axis(angle = 90))
DotPlot(object = TC_sub, features = c("ST18","CD34","IGLL1", "TRGC2","TRDC",  "HIVEP3","RGPD3", "PTCRA", "TP53INP1", "SMPD3","PCNA","CDK1","MKI67","TRBC2","CD4","CD8A","CD8B",  "RAG1","RAG2","AQP3","RORC","SATB1","CCND3", "TOX2","CCR9","TRAC", "CCR7","CD5","CD27","CCND2"),assay="RNA") + coord_flip()+ scale_color_viridis(option = "H")+theme(axis.text.x = element_text(angle = 45)) + ggtitle("TC_modulescore") + scale_y_discrete(guide = guide_axis(angle = 90))
DotPlot(object = TC_sub, features = module_TC_yasumizu_signature,assay="RNA") + coord_flip()+ scale_color_viridis(option = "H")+theme(axis.text.x = element_text(angle = 45)) + ggtitle("TC_modulescore") + scale_y_discrete(guide = guide_axis(angle = 90))
DotPlot(object = TC_sub, features = module_TC_DEG_yasumizu_signature,assay="RNA") + coord_flip()+ scale_color_viridis(option = "H")+theme(axis.text.x = element_text(angle = 45)) + ggtitle("TC_modulescore") + scale_y_discrete(guide = guide_axis(angle = 90))
DotPlot(object = TC_sub, features = yayon_TC,assay="RNA") + coord_flip()+ scale_color_viridis(option = "H")+theme(axis.text.x = element_text(angle = 45)) + ggtitle("TC_modulescore") + scale_y_discrete(guide = guide_axis(angle = 90))

#Modulescore imaging
pdf("Graphs/Main/Modulescore_TC_Clusteridentification.pdf", width=10, height=15)
DotPlot(object = TC_sub, features = TC_module_list,assay="RNA") + coord_flip()+ scale_color_viridis(option = "H")+theme(axis.text.x = element_text(angle = 45)) + ggtitle("TC_modulescore") + scale_y_discrete(guide = guide_axis(angle = 90))
DotPlot(object = TC_sub, features = yayon_TC,assay="RNA") + coord_flip()+ scale_color_viridis(option = "H")+theme(axis.text.x = element_text(angle = 45)) + ggtitle("TC_modulescore") + scale_y_discrete(guide = guide_axis(angle = 90))
dev.off()



#park:  ETP:HOXA9  DN: CD3G, HES1, TRDC, DP:PTPRC, CD3G, RORC, CD8A, CD8B, CD4  aßT(entry):PTPRC, CD3G, CCR9,   CD8aa:PTPRC, CD3G, PDCD1, CD8A  CD8+T:CD8A, CD8B, PTPRC, CD3G,  CD8+Tmem: ANXA1, CRTAM, CD8A, CD8B, CD3G, PTPRC  CD4+T:PTPRC, CD3G, CD40LG, DHRS3  CD4+Tmem:ANXA1, CD40LG, PTPRC, CD3G,  Treg: FOXP3, PTPRC, CD3G  T(agonist):DHRS3, CD40LG, PTPRC, CD3G  ydT:TRDC, ANXA1, CD8B, CCR9, HES1, CD3G, PTPRC  NKT:KLRD1, EOMES, TRDC, CD8A, CD3G, PTPRC
#park track: DN(early):"ST18""CD34","IGLL1", DN(P):"TRGC2","TRDC",  DN(Q):"HIVEP3","RGPD3","TRDC", "PTCRA", "TP53INP1", DP(P):"SMPD3","PCNA","CDK1","MKI67","TRBC2","CD4","CD8A","CD8B",  DP(Q): "RAG1","RAG2","AQP3","RORC","SATB1","TRBC2","CD4","CD8A","CD8B","CCND3", aßT(entry):"TOX2","SATB1","CCR9","TRAC",  CD8+T/CD4+T:"CCR7","CD5","CD27","CCND2",



#Name Cluster
Idents(TC_sub)<-TC_sub$seurat_clusters
TC_sub<-RenameIdents(TC_sub,
                     `0`="CD4_CD8", `1`="DP_Q_1", `2`="DP_P_1", `3`="aßT_entry",`4`="DN_P_1",`5`="DP_Q_2",`6`="DP_P_2",`7`="DN_P_2",`8`="naiveTC_TET1",`9`="naiveTC_TET2",`10`="DP",
                     `11`="DP",`12`="DNearly_T_ETP",`13`="TH_like",`14`="DP")
TC_sub$celltype.full<-Idents(TC_sub)

TC_sub_backup<-TC_sub

table(TC_sub$celltype.full)

#Doublet and contamination check
DimPlot(TC_sub, label = TRUE, label.size = 6)
FeaturePlot(TC_sub,c("KRT19","CD3D"), blend=T, order=T, label=T, cols= c("lightgrey","green","blue"), max.cutoff = 2)

cells.of.interest_2 <- WhichCells(
  TC_sub,
  idents = "DP"
)

TC_sub <- subset(TC_sub,idents = c("DNearly_T_ETP", "DN_P_1","DN_P_2", "DP_Q_1","DP_Q_2", "DP_P_1", "DP_P_2","aßT_entry","CD4_CD8", "TH_like","naiveTC_TET1","naiveTC_TET2"))


#Define cluster levels
clusters_ordered<-c("DNearly_T_ETP", "DN_P_1","DN_P_2", "DP_Q_1","DP_Q_2", "DP_P_1", "DP_P_2","aßT_entry","CD4_CD8", "TH_like","naiveTC_TET1","naiveTC_TET2")
TC_sub$celltype.full<- factor(TC_sub$celltype.full, levels = clusters_ordered)
Idents(TC_sub)<-TC_sub$celltype.full


DefaultAssay(TC_sub)<-"RNA"
Idents(TC_sub)<-"celltype.full"

pdf("Graphs/Main/TC_UMAPPlots.pdf", width=5, height=5)
DimPlot(TC_sub, label=T, label.color = "black",label.size = 5, repel=T, alpha=0.4) + ggtitle("TC-subset")+ scale_color_manual(values=color_levelTC) + NoLegend()
dev.off()

pdf("Graphs/Main/UMAPPlots_TC_sub_QC.pdf", width=5, height=5)
DimPlot(TC_sub, group.by="group")+ ggtitle("Basic-UMAP")
DimPlot(TC_sub, group.by="tissue")+ ggtitle("Basic-UMAP")
dev.off()

pdf("Graphs/Main/Modulescore_TC_Clusteridentification.pdf", width=10, height=15)
DotPlot(object = TC_sub, features = TC_module_list,assay="RNA") + coord_flip()+ scale_color_viridis(option = "H")+theme(axis.text.x = element_text(angle = 45)) + ggtitle("TC_modulescore") + scale_y_discrete(guide = guide_axis(angle = 90))
dev.off()

pdf("Graphs/Main/Modulescore_TC_Clusteridentification_Yayon.pdf", width=7, height=7)
DotPlot(object = TC_sub, features = yayon_TC,assay="RNA") + coord_flip()+ scale_color_viridis(option = "H")+theme(axis.text.x = element_text(angle = 45)) + ggtitle("TC_modulescore") + scale_y_discrete(guide = guide_axis(angle = 90))
dev.off()

pdf("Graphs/Main/TC_Barplots.pdf")
plot_integrated_clusters_group (TC_sub)
plot_integrated_clusters_Donor_ID (TC_sub)
plot_integrated_celltypes_condition (TC_sub)
plot_integrated_celltypes_classification (TC_sub)
plot_integrated_celltypes_tissue (TC_sub)
plot_integrated_celltypes_sex (TC_sub)
plot_integrated_celltypes_age (TC_sub)
plot_integrated_celltypes_development (TC_sub)
plot_integrated_celltypes_donor (TC_sub)
plot_TC_integrated_condition_celltypes(TC_sub)
plot_TC_integrated_development_celltypes(TC_sub)
plot_TC_integrated_group_clusters(TC_sub)
plot_TC_integrated_tissue_celltypes(TC_sub)
dev.off()

#add identity to analyze tissue/health diffgenes
#Set tissue in identy
Idents(TC_sub) <- TC_sub$celltype.full
TC_sub$celltype.condition <- paste(Idents(TC_sub), TC_sub$condition, sep = "_")
Idents(TC_sub)<-TC_sub$celltype.condition
table(TC_sub$celltype.condition)


#Clustermarker
Idents(TC_sub)<-"celltype.full"
TC_sub.clustermarker<-FindAllMarkers(TC_sub, assay = "RNA", min.pct = 0.25, logfc.threshold = 0.25)
TC_sub.clustermarker$Foldchange_UP <- 2^(TC_sub.clustermarker$avg_log2FC)
TC_sub.clustermarker$Foldchange_DOWN <- 2^(-TC_sub.clustermarker$avg_log2FC)
TC_sub.clustermarker$Ratio_pct1_pct2 <- (TC_sub.clustermarker$pct.1)/(TC_sub.clustermarker$pct.2)
write.xlsx(TC_sub.clustermarker, "Lists/Clustermarker_TC_sub_clustermarker.xlsx")

#DotPlot Top up and down regulated genes
#up
TC_sub.clustermarker_UP <- subset(TC_sub.clustermarker, Foldchange_UP >=2 )
top10_UP_clustermarker <- TC_sub.clustermarker_UP %>% group_by(cluster) %>% top_n(n=10,wt=Foldchange_UP)
top10_UP_clustermarker<-arrange(top10_UP_clustermarker,cluster,Foldchange_UP)
top10_UP_clustermarker <- as.character(top10_UP_clustermarker$gene)
top10_UP_clustermarker<-top10_UP_clustermarker[!duplicated(top10_UP_clustermarker)]

pdf("Graphs/Main/DotPlot_Top_10UP_Clustermarker_TC_FCmin2.pdf", width = 7, height = 20)
DotPlot(TC_sub, features=top10_UP_clustermarker,group.by = "celltype.full", assay="RNA") + coord_flip()+ scale_color_viridis(option = "H")+theme(axis.text.x = element_text(angle = 45)) + ggtitle("TC Cluster_DEG") + scale_y_discrete(guide = guide_axis(angle = 90))
dev.off()

#Heatmap#
## Define Gene lists that you want to plot as Heatmap
TC_sub.clustermarker_UP<-arrange(TC_sub.clustermarker_UP,cluster,Foldchange_UP)
genes <- as.character(TC_sub.clustermarker_UP$gene)
avgexp <- TC_sub

## Plot Heatmaps
## Calculate Average gene expressions (pseudobulked) for each of the celltypes in each condition
levels(avgexp)
DefaultAssay(avgexp) <- "RNA"
avgexp <- NormalizeData(avgexp)
avgexp <- AggregateExpression(avgexp, return.seurat =T, assays = "RNA")


# Modify the heatmap plot
library(RColorBrewer)
heatmap_plot <- DoHeatmap(
  avgexp,
  features = genes,
  size = 6,
  group.bar.height = 0.02,
  draw.lines = FALSE,angle = 90,
  group.colors = color_levelTC
) +
  scale_fill_gradientn(colours = rev(brewer.pal(8, name = "RdYlBu")))

# Modify the theme to make rownames cursive
heatmap_plot 
ggsave("Graphs/Main/Heatmap_TC_clustermarker.pdf", width = 40, height = 20, units = c("cm"))


#GO-Terms_DNearly_T_ETP
table(TC_sub.clustermarker_UP$cluster)
DNearly_T_ETP_sub.clustermarker_UP <- subset(TC_sub.clustermarker_UP, cluster=="DNearly_T_ETP")
DNearly_T_ETP_sub.clustermarker_UP <- as.character(DNearly_T_ETP_sub.clustermarker_UP$gene)

if (websiteLive) {
  enriched <- enrichr(DNearly_T_ETP_sub.clustermarker_UP, dbs)
}

if (websiteLive) enriched[["GO_Biological_Process_2023"]]
pdf("Graphs/Main/UP_Enrichment_DNearly_T_ETP.pdf")
if (websiteLive) {
  plotEnrich(enriched[[3]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")
}
dev.off()

#GO-Terms_DN_P_1
table(TC_sub.clustermarker_UP$cluster)
DN_P_1_sub.clustermarker_UP <- subset(TC_sub.clustermarker_UP, cluster=="DN_P_1")
DN_P_1_sub.clustermarker_UP <- as.character(DN_P_1_sub.clustermarker_UP$gene)

if (websiteLive) {
  enriched <- enrichr(DN_P_1_sub.clustermarker_UP, dbs)
}

if (websiteLive) enriched[["GO_Biological_Process_2023"]]
pdf("Graphs/Main/UP_Enrichment_DN_P_1.pdf")
if (websiteLive) {
  plotEnrich(enriched[[3]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")
}
dev.off()

#GO-Terms_DN_P_2
table(TC_sub.clustermarker_UP$cluster)
DN_P_2_sub.clustermarker_UP <- subset(TC_sub.clustermarker_UP, cluster=="DN_P_2")
DN_P_2_sub.clustermarker_UP <- as.character(DN_P_2_sub.clustermarker_UP$gene)

if (websiteLive) {
  enriched <- enrichr(DN_P_2_sub.clustermarker_UP, dbs)
}

if (websiteLive) enriched[["GO_Biological_Process_2023"]]
pdf("Graphs/Main/UP_Enrichment_DN_P_2.pdf")
if (websiteLive) {
  plotEnrich(enriched[[3]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")
}
dev.off()

#GO-Terms_DP_Q_1 
table(TC_sub.clustermarker_UP$cluster)
DP_Q_1_sub.clustermarker_UP <- subset(TC_sub.clustermarker_UP, cluster=="DP_Q_1")
DP_Q_1_sub.clustermarker_UP <- as.character(DP_Q_1_sub.clustermarker_UP$gene)

if (websiteLive) {
  enriched <- enrichr(DP_Q_1_sub.clustermarker_UP, dbs)
}

if (websiteLive) enriched[["GO_Biological_Process_2023"]]
pdf("Graphs/Main/UP_Enrichment_DP_Q_1.pdf")
if (websiteLive) {
  plotEnrich(enriched[[3]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")
}
dev.off()

#GO-Terms_DP_Q_2
table(TC_sub.clustermarker_UP$cluster)
DP_Q_2_sub.clustermarker_UP <- subset(TC_sub.clustermarker_UP, cluster=="DP_Q_2")
DP_Q_2_sub.clustermarker_UP <- as.character(DP_Q_2_sub.clustermarker_UP$gene)

if (websiteLive) {
  enriched <- enrichr(DP_Q_2_sub.clustermarker_UP, dbs)
}

if (websiteLive) enriched[["GO_Biological_Process_2023"]]
pdf("Graphs/Main/UP_Enrichment_DP_Q_2.pdf")
if (websiteLive) {
  plotEnrich(enriched[[3]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")
}
dev.off()

#GO-Terms_DP_P_1
table(TC_sub.clustermarker_UP$cluster)
DP_P_1_sub.clustermarker_UP <- subset(TC_sub.clustermarker_UP, cluster=="DP_P_1")
DP_P_1_sub.clustermarker_UP <- as.character(DP_P_1_sub.clustermarker_UP$gene)

if (websiteLive) {
  enriched <- enrichr(DP_P_1_sub.clustermarker_UP, dbs)
}

if (websiteLive) enriched[["GO_Biological_Process_2023"]]
pdf("Graphs/Main/UP_Enrichment_DP_P_1.pdf")
if (websiteLive) {
  plotEnrich(enriched[[3]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")
}
dev.off()

#GO-Terms_DP_P_2
table(TC_sub.clustermarker_UP$cluster)
DP_P_2_sub.clustermarker_UP <- subset(TC_sub.clustermarker_UP, cluster=="DP_P_2")
DP_P_2_sub.clustermarker_UP <- as.character(DP_P_2_sub.clustermarker_UP$gene)

if (websiteLive) {
  enriched <- enrichr(DP_P_2_sub.clustermarker_UP, dbs)
}

if (websiteLive) enriched[["GO_Biological_Process_2023"]]
pdf("Graphs/Main/UP_Enrichment_DP_P_2.pdf")
if (websiteLive) {
  plotEnrich(enriched[[3]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")
}
dev.off()

#GO-Terms_aßT_entry
table(TC_sub.clustermarker_UP$cluster)
aßT_entry_sub.clustermarker_UP <- subset(TC_sub.clustermarker_UP, cluster=="aßT_entry")
aßT_entry_sub.clustermarker_UP <- as.character(aßT_entry_sub.clustermarker_UP$gene)

if (websiteLive) {
  enriched <- enrichr(aßT_entry_sub.clustermarker_UP, dbs)
}

if (websiteLive) enriched[["GO_Biological_Process_2023"]]
pdf("Graphs/Main/UP_Enrichment_aßT_entry.pdf")
if (websiteLive) {
  plotEnrich(enriched[[3]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")
}
dev.off()

#GO-Terms_CD4_CD8
table(TC_sub.clustermarker_UP$cluster)
CD4_CD8_sub.clustermarker_UP <- subset(TC_sub.clustermarker_UP, cluster=="CD4_CD8")
CD4_CD8_sub.clustermarker_UP <- as.character(CD4_CD8_sub.clustermarker_UP$gene)

if (websiteLive) {
  enriched <- enrichr(CD4_CD8_sub.clustermarker_UP, dbs)
}

if (websiteLive) enriched[["GO_Biological_Process_2023"]]
pdf("Graphs/Main/UP_Enrichment_CD4_CD8.pdf")
if (websiteLive) {
  plotEnrich(enriched[[3]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")
}
dev.off()

#GO-Terms_TH_like
table(TC_sub.clustermarker_UP$cluster)
TH_like_sub.clustermarker_UP <- subset(TC_sub.clustermarker_UP, cluster=="TH_like")
TH_like_sub.clustermarker_UP <- as.character(TH_like_sub.clustermarker_UP$gene)

if (websiteLive) {
  enriched <- enrichr(TH_like_sub.clustermarker_UP, dbs)
}

if (websiteLive) enriched[["GO_Biological_Process_2023"]]
pdf("Graphs/Main/UP_Enrichment_TH_like.pdf")
if (websiteLive) {
  plotEnrich(enriched[[3]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")
}
dev.off()

#GO-Terms_naiveTC_TET1
table(TC_sub.clustermarker_UP$cluster)
naiveTC_TET1_sub.clustermarker_UP <- subset(TC_sub.clustermarker_UP, cluster=="naiveTC_TET1")
naiveTC_TET1_sub.clustermarker_UP <- as.character(naiveTC_TET1_sub.clustermarker_UP$gene)

if (websiteLive) {
  enriched <- enrichr(naiveTC_TET1_sub.clustermarker_UP, dbs)
}

if (websiteLive) enriched[["GO_Biological_Process_2023"]]
pdf("Graphs/Main/UP_Enrichment_naiveTC_TET1.pdf")
if (websiteLive) {
  plotEnrich(enriched[[3]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")
}
dev.off()

#GO-Terms_naiveTC_TET2
table(TC_sub.clustermarker_UP$cluster)
naiveTC_TET2_sub.clustermarker_UP <- subset(TC_sub.clustermarker_UP, cluster=="naiveTC_TET2")
naiveTC_TET2_sub.clustermarker_UP <- as.character(naiveTC_TET2_sub.clustermarker_UP$gene)

if (websiteLive) {
  enriched <- enrichr(naiveTC_TET2_sub.clustermarker_UP, dbs)
}

if (websiteLive) enriched[["GO_Biological_Process_2023"]]
pdf("Graphs/Main/UP_Enrichment_naiveTC_TET2.pdf")
if (websiteLive) {
  plotEnrich(enriched[[3]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")
}
dev.off()

save(TC_sub,file="TC_sub.RData")

#(22)_Subset FB####
Idents(th.total)<- th.total$celltype.full
FB_sub <- subset(th.total,idents = "FB")

DefaultAssay(FB_sub)<-"RNA"
FB_sub[["RNA"]]$data <- as(FB_sub[["RNA"]]$data, Class = "dgCMatrix")

FB_sub <- FindVariableFeatures(FB_sub)
FB_sub <- ScaleData(FB_sub)
FB_sub <- RunPCA(FB_sub)
ElbowPlot(FB_sub, ndims=20)
FB_sub <- RunUMAP(FB_sub, dims = 1:9)
FB_sub <- FindNeighbors(FB_sub, dims = 1:9)
FB_sub <- FindClusters(FB_sub, resolution = 0.1)
DimPlot(FB_sub, label = T, label.size = 3, alpha = 1)
DimPlot(FB_sub, label=T, group.by="condition", alpha = 1)
DimPlot(FB_sub, label=T, split.by="development", alpha = 1)
DimPlot(FB_sub, label=T, split.by="group", alpha = 1)
DimPlot(FB_sub, label=T, group.by="seurat_clusters", alpha = 1)
DefaultAssay(FB_sub)<- "RNA"

#Modulscore:
#Ascension PMID: 33385399
A1_Ascension_signature <- list(c("ANGPTL5","C1QTNF3","CD151","CD55","CD99","CPE","CTSB","CYBRD1","DCN","FBLN1","FGL2","GPX3","GSN","IGFBP6","LOX","MFAP5",
                                 "MGST1","MMP2","OLFML3","PDGFRL","PI16","PIGT","PODN","REXO2","SCARA5","SERPINF1","SLPI","TSPAN8","XG"))
FB_sub <- AddModuleScore(FB_sub, features = A1_Ascension_signature, name = "A1_Ascension module score")
A2_Ascension_signature <- list(c("APCDD1","AXIN2","C1orf198","CLEC2A","COL13A1","COL18A1","COL23A1","COMP","CTSC","EMX2","F13A1","GNG11","GREM2",
                                 "HSPB3","ID1","LAMC3","NKD2","NPTX2","PTK7","RGS2","RGS3","RSPO1","SPRY1","STC2","TGFBI"))
FB_sub <- AddModuleScore(FB_sub, features = A2_Ascension_signature, name = "A2_Ascension module score")
A3_Ascension_signature <- list(c("CD9","COL6A1","ELN","LEPR","RGCC","SGCA","WIF1"))
FB_sub <- AddModuleScore(FB_sub, features = A3_Ascension_signature, name = "A3_Ascension module score")
A4_Ascension_signature <- list(c("C1QTNF3","FBN1","FSTL1","HSD3B7","IGFBP6","ISLR","MFAP5","PCOLCE2","PRG4","PRSS23","SFRP4","TNXB"))
FB_sub <- AddModuleScore(FB_sub, features = A4_Ascension_signature, name = "A4_Ascension module score")
B1_Ascension_signature <- list(c("ARL6IP1","CCL2","CXCL1","CXCL2","CXCL3","DNAJA1","ERRFI1","GPC3","ITM2A","MCL1","MYC","PLA2G2A","PLIN2","SOD2","UAP1"))
FB_sub <- AddModuleScore(FB_sub, features = B1_Ascension_signature, name = "B1_Ascension module score")
B2_Ascension_signature <- list(c("BIRC3","BTG1","C3","CCDC146","CCL19","CD74","CTSH","IGFBP3","OLFM2","PSME2","TNFSF13B"))
FB_sub <- AddModuleScore(FB_sub, features = B2_Ascension_signature, name = "B2_Ascension module score")
C1_Ascension_signature <- list(c("CCND1","CDH11","COL11A1","COL5A2","DPEP1","EDNRA","GPC3","LAMC3","MEF2C","MME","POSTN","SPARC","STMN1"))
FB_sub <- AddModuleScore(FB_sub, features = C1_Ascension_signature, name = "C1_Ascension module score")
C2_Ascension_signature <- list(c("ADAMTS6","ARHGAP15","CADM2","CCK","CHADL","CLEC14A","COCH","CRABP1","DKK2","EMID1","FIBIN","FZD1","GAP43","HSPA2",
                                 "LIPA","MEIS2","MEOX2","MKX","NDNF","NECAB1","OGN","PCSK1N","PLXDC1","PRRG4","RHPN1","RSPO4","SLC22A16","SLITRK6","SYTL2"))
FB_sub <- AddModuleScore(FB_sub, features = C2_Ascension_signature, name = "C2_Ascension module score")
C3_Ascension_signature <- list(c("ASPN","BGN","DIO2","DKK2","F2R","FIBIN","GPM6B","LRRC15","LTBP2","MARCKS","PLEKHH2","PMEPA1","POSTN"))
FB_sub <- AddModuleScore(FB_sub, features = C3_Ascension_signature, name = "C3_Ascension module score")
C4_Ascension_signature <- list(c("ANGPTL7","APOD","CLDN1","CYP1B1","EBF2","EIF4A3","FGFBP2","IFI27","KLK1","PODNL1","SCN7A","SFRP4","TAGLN","TENM2","TM4SF1"))
FB_sub <- AddModuleScore(FB_sub, features = C4_Ascension_signature, name = "C4_Ascension module score")

tumo_assoc_fb_yasumizu_signature <- list(c("PDGFRA","ADH1B","FBN1","COL1A1","S100A4"))
FB_sub <- AddModuleScore(FB_sub, features = tumo_assoc_fb_yasumizu_signature, name = "tumo_assoc_fb_yasumizu module score")
norm_fb_yasumizu_signature <- list(c("FN1","EGFL6","MATN4","PTN","VIM"))
FB_sub <- AddModuleScore(FB_sub, features = norm_fb_yasumizu_signature, name = "norm_fb_yasumizu module score")

module_FB_Ascension_signature <- c("A1_Ascension module score1","A2_Ascension module score1","A3_Ascension module score1","A4_Ascension module score1","B1_Ascension module score1","B2_Ascension module score1","C1_Ascension module score1","C2_Ascension module score1","C3_Ascension module score1","C4_Ascension module score1")

module_FB_DEG_yasumizu_signature <- c("tumo_assoc_fb_DEG_yasumizu module score1","norm_fb_DEG_yasumizu module score1")


#Modulescore imaging
pdf("Graphs/Main/Modulescore_FB_Clusteridentification.pdf", width=10, height=10)
DotPlot(object = FB_sub, features = module_FB_Ascension_signature,assay="RNA") + coord_flip()+ scale_color_viridis(option = "H")+theme(axis.text.x = element_text(angle = 45)) + ggtitle("FB_modulescore") + scale_y_discrete(guide = guide_axis(angle = 90))
DotPlot(object = FB_sub, features = module_FB_DEG_yasumizu_signature,assay="RNA") + scale_color_viridis(option = "H")+theme(axis.text.x = element_text(angle = 45)) + ggtitle("FB_modulescore") + scale_y_discrete(guide = guide_axis(angle = 0))
DotPlot(object = FB_sub, features = yayon_FB,assay="RNA") + coord_flip()+ scale_color_viridis(option = "H")+theme(axis.text.x = element_text(angle = 45)) + ggtitle("FB_modulescore") + scale_y_discrete(guide = guide_axis(angle = 90))
dev.off()

DotPlot(FB_sub,features=c("ELN","QPCT","MMP2","SFRP2",
                          "DKK3","TNMD","SFRP1","TNN",
                          "APOE","C7","CYGB","IGFBP7"))
DotPlot(FB_sub,features=c("FBN1","PRG4","PCOLCE2","SFRP4","IGFBP6","SLPI","PI16","ELN","RGCC","SGCA","WIF1","APCDD1","COMP","COL18A1","NKD2"))
DotPlot(FB_sub,features=c("ASPN","F2R","GPM6B","POSTN","COCH","FIBIN","CRABP1","RSPO4","COL11A1","DPEP1","TNMD","WFDC1","ANGPTL7","APOD","TM4SF1"))
DotPlot(FB_sub,features=c("CCL2","ITM2A","SPSB1","TNFAIP6","CCDC146","CCL19","CD74"))
DotPlot(FB_sub,c("LRRC15","ARL4C","COL11A1","MMP11","EDIL3","TOMM6","SDC1"))
FeaturePlot(FB_sub,"ACTA2", order=T, split.by="tissue")

#Name Cluster
Idents(FB_sub)<-FB_sub$seurat_clusters
FB_sub<-RenameIdents(FB_sub,
                     `0`="fetal_FB_SFRP1", `1`="medFB_CCL2", `2`="fetal_FB_POSTN", `3`="TET_FB",`4`="PeriloFB",`5`="medFB_MHCIIh",
                     `6`="InterloFB_COL9A3",`7`="Perilo_Interlo_medFB",`8`="InterloFB",`9`="PeriloFB_Prolif",`10`="fetal_FB_SFRP1",`11`="medFB_CCL2",`12`="medFB_MHCIIh")
FB_sub$celltype.full<-Idents(FB_sub)
table(FB_sub$celltype.full)

#Doublet and contamination check
DimPlot(FB_sub, label = TRUE, label.size = 6)
FeaturePlot(FB_sub,c("KRT19","CD3D"), blend=T, order=T, label=T, cols= c("lightgrey","green","blue"), max.cutoff = 2)

#Define cluster levels
clusters_ordered<-c("PeriloFB_Prolif","PeriloFB", "InterloFB","InterloFB_COL9A3", "Perilo_Interlo_medFB", "medFB_MHCIIh", "medFB_CCL2", "TET_FB","fetal_FB_POSTN","fetal_FB_SFRP1")
FB_sub$celltype.full<- factor(FB_sub$celltype.full, levels = clusters_ordered)
Idents(FB_sub)<-FB_sub$celltype.full

DefaultAssay(FB_sub)<-"RNA"
Idents(FB_sub)<-"celltype.full"

pdf("Graphs/Main/FB_UMAPPlots.pdf", width=5, height=5)
DimPlot(FB_sub, label=T, label.color = "black",label.size = 5, repel=T, alpha=1) + ggtitle("FB-subset")+ scale_color_manual(values=color_levelFB) + NoLegend()
dev.off()

pdf("Graphs/Main/FB_split_celltype_UMAPPlots.pdf", width=30, height=5)
DimPlot(FB_sub, label=T, label.color = "black",label.size = 5, repel=T, alpha=1, split.by="celltype.full") + ggtitle("FB-subset")+ scale_color_manual(values=color_levelFB) + NoLegend()
dev.off()


pdf("Graphs/Main/UMAPPlots_FB_sub_QC.pdf", width=5, height=5)
DimPlot(FB_sub, group.by="group")+ ggtitle("Basic-UMAP")
DimPlot(FB_sub, group.by="tissue")+ ggtitle("Basic-UMAP")
dev.off()



pdf("Graphs/Main/FB_Barplots.pdf")
plot_integrated_clusters_group (FB_sub)
plot_integrated_clusters_Donor_ID (FB_sub)
plot_integrated_celltypes_condition (FB_sub)
plot_integrated_celltypes_classification (FB_sub)
plot_integrated_celltypes_tissue (FB_sub)
plot_integrated_celltypes_sex (FB_sub)
plot_integrated_celltypes_age (FB_sub) 
plot_integrated_celltypes_development (FB_sub)
plot_integrated_celltypes_donor (FB_sub)
plot_FB_integrated_condition_celltypes(FB_sub)
plot_FB_integrated_development_celltypes(FB_sub)
plot_FB_integrated_group_clusters(FB_sub)
plot_FB_integrated_tissue_celltypes(FB_sub)
dev.off()

#Cell_cluster_count
Idents(FB_sub) <- FB_sub$celltype.full
FB_sub$celltype.tissue <- paste(Idents(FB_sub), FB_sub$tissue, sep = "_")
Idents(FB_sub)<-FB_sub$celltype.tissue
levels(FB_sub)
Idents(FB_sub) <- FB_sub$celltype.tissue
FB_sub$celltype.tissue.donor <- paste(Idents(FB_sub), FB_sub$donor, sep = "_")
Idents(FB_sub)<-FB_sub$celltype.tissue.donor
levels(FB_sub)
table(FB_sub$celltype.tissue.donor)
write.xlsx(table(FB_sub$celltype.tissue.donor), "Lists/FB_sub_ctd.xlsx")

#add identity to analyze tissue/health diffgenes
#Set tissue in identy 
Idents(FB_sub) <- FB_sub$celltype.full
FB_sub$celltype.condition <- paste(Idents(FB_sub), FB_sub$condition, sep = "_")
Idents(FB_sub)<-FB_sub$celltype.condition
table(FB_sub$celltype.condition)
#Set tissue in identy 
Idents(FB_sub) <- FB_sub$celltype.condition
FB_sub$celltype.condition.development <- paste(Idents(FB_sub), FB_sub$development, sep = "_")
Idents(FB_sub)<-FB_sub$celltype.condition.development
table(FB_sub$celltype.condition.development)

#Clustermarker
Idents(FB_sub)<-"celltype.full"
FB_sub.clustermarker<-FindAllMarkers(FB_sub, assay = "RNA", min.pct = 0.25, logfc.threshold = 0.25)
FB_sub.clustermarker$Foldchange_UP <- 2^(FB_sub.clustermarker$avg_log2FC)
FB_sub.clustermarker$Foldchange_DOWN <- 2^(-FB_sub.clustermarker$avg_log2FC)
FB_sub.clustermarker$Ratio_pct1_pct2 <- (FB_sub.clustermarker$pct.1)/(FB_sub.clustermarker$pct.2)
write.xlsx(FB_sub.clustermarker, "Lists/Clustermarker_FB_sub_clustermarker.xlsx")

#DotPlot Top up and down regulated genes 
#up
FB_sub.clustermarker_UP <- subset(FB_sub.clustermarker, Foldchange_UP >=2 )
top10_UP_clustermarker <- FB_sub.clustermarker_UP %>% group_by(cluster) %>% top_n(n=10,wt=Foldchange_UP)
top10_UP_clustermarker<-arrange(top10_UP_clustermarker,cluster,Foldchange_UP)
top10_UP_clustermarker <- as.character(top10_UP_clustermarker$gene)
top10_UP_clustermarker<-top10_UP_clustermarker[!duplicated(top10_UP_clustermarker)]

pdf("Graphs/Main/DotPlot_Top_10UP_Clustermarker_FB_FCmin2.pdf", width = 7, height = 20)
DotPlot(FB_sub, features=top10_UP_clustermarker, group.by = "celltype.full",assay="RNA") + coord_flip()+ scale_color_viridis(option = "H")+theme(axis.text.x = element_text(angle = 45)) + ggtitle("FB Cluster_DEG") + scale_y_discrete(guide = guide_axis(angle = 90))
dev.off()


pdf("Graphs/Main/DotPlot_Top_20UP_Clustermarker_TET_FB_FCmin2.pdf", width = 6, height = 8)
DotPlot(FB_sub, features=c("LRRC15","MMP11","COL11A1","EDIL3","TXNDC5","HMCN1","HTRA1","CTHRC1","PMEPA1","ITGA11","FN1","ITGB5","MEGF6","TWIST1","ARL4C","SPON1",
                           "APCDD1","C17orf49","SFRP4","NBL1", "SDC1"), group.by = "celltype.full",assay="RNA") + coord_flip()+ scale_color_viridis(option = "H")+theme(axis.text.x = element_text(angle = 45)) + ggtitle("TEC Cluster_DEG") + scale_y_discrete(guide = guide_axis(angle = 90))
dev.off()

#Feature Marker TETFB
pdf("Graphs/Main/Feature_KRT14_TET_FB.pdf", width=5, height=5)
FeaturePlot(FB_sub, c("THY1"), order=T, label=T,pt.size=0.01, max.cutoff = 3)
dev.off()


#Heatmap#
## Define Gene lists that you want to plot as Heatmap
FB_sub.clustermarker_UP<-arrange(FB_sub.clustermarker_UP,cluster,Foldchange_UP)
genes <- as.character(FB_sub.clustermarker_UP$gene)
avgexp <- FB_sub

## Plot Heatmaps
## Calculate Average gene expressions (pseudobulked) for each of the celltypes in each condition
levels(avgexp)
DefaultAssay(avgexp) <- "RNA"
avgexp <- NormalizeData(avgexp)
avgexp <- AggregateExpression(avgexp, return.seurat =T, assays = "RNA")


# Modify the heatmap plot
library(RColorBrewer)
heatmap_plot <- DoHeatmap(
  avgexp,
  features = genes,
  size = 6,
  group.bar.height = 0.02,
  draw.lines = FALSE,angle = 90,
  group.colors = color_levelFB
) +
  scale_fill_gradientn(colours = rev(brewer.pal(8, name = "RdYlBu")))

# Modify the theme to make rownames cursive
heatmap_plot 
ggsave("Graphs/Main/Heatmap_FB_clustermarker.pdf", width = 40, height = 20, units = c("cm"))

#Enrichment analysis FB_sub
#GO-Terms_PeriloFB_Prolif
table(FB_sub.clustermarker_UP$cluster)
PeriloFB_Prolif_sub.clustermarker_UP <- subset(FB_sub.clustermarker_UP, cluster=="PeriloFB_Prolif")
PeriloFB_Prolif_sub.clustermarker_UP <- as.character(PeriloFB_Prolif_sub.clustermarker_UP$gene)

if (websiteLive) {
  enriched <- enrichr(PeriloFB_Prolif_sub.clustermarker_UP, dbs)
}

if (websiteLive) enriched[["GO_Biological_Process_2023"]]
pdf("Graphs/Main/UP_Enrichment_PeriloFB_Prolif.pdf")
if (websiteLive) {
  plotEnrich(enriched[[3]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")
}
dev.off()

#GO-Terms_PeriloFB
table(FB_sub.clustermarker_UP$cluster)
PeriloFB_sub.clustermarker_UP <- subset(FB_sub.clustermarker_UP, cluster=="PeriloFB")
PeriloFB_sub.clustermarker_UP <- as.character(PeriloFB_sub.clustermarker_UP$gene)

if (websiteLive) {
  enriched <- enrichr(PeriloFB_sub.clustermarker_UP, dbs)
}

if (websiteLive) enriched[["GO_Biological_Process_2023"]]
pdf("Graphs/Main/UP_Enrichment_PeriloFB.pdf")
if (websiteLive) {
  plotEnrich(enriched[[3]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")
}
dev.off()

#GO-Terms_InterloFB
table(FB_sub.clustermarker_UP$cluster)
InterloFB_sub.clustermarker_UP <- subset(FB_sub.clustermarker_UP, cluster=="InterloFB")
InterloFB_sub.clustermarker_UP <- as.character(InterloFB_sub.clustermarker_UP$gene)

if (websiteLive) {
  enriched <- enrichr(InterloFB_sub.clustermarker_UP, dbs)
}

if (websiteLive) enriched[["GO_Biological_Process_2023"]]
pdf("Graphs/Main/UP_Enrichment_InterloFB.pdf")
if (websiteLive) {
  plotEnrich(enriched[[3]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")
}
dev.off()

#GO-Terms_InterloFB_COL9A3
table(FB_sub.clustermarker_UP$cluster)
InterloFB_COL9A3_sub.clustermarker_UP <- subset(FB_sub.clustermarker_UP, cluster=="InterloFB_COL9A3")
InterloFB_COL9A3_sub.clustermarker_UP <- as.character(InterloFB_COL9A3_sub.clustermarker_UP$gene)

if (websiteLive) {
  enriched <- enrichr(InterloFB_COL9A3_sub.clustermarker_UP, dbs)
}

if (websiteLive) enriched[["GO_Biological_Process_2023"]]
pdf("Graphs/Main/UP_Enrichment_InterloFB_COL9A3.pdf")
if (websiteLive) {
  plotEnrich(enriched[[3]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")
}
dev.off()

#GO-Terms_Perilo_Interlo_medFB
table(FB_sub.clustermarker_UP$cluster)
Perilo_Interlo_medFB_sub.clustermarker_UP <- subset(FB_sub.clustermarker_UP, cluster=="Perilo_Interlo_medFB")
Perilo_Interlo_medFB_sub.clustermarker_UP <- as.character(Perilo_Interlo_medFB_sub.clustermarker_UP$gene)

if (websiteLive) {
  enriched <- enrichr(Perilo_Interlo_medFB_sub.clustermarker_UP, dbs)
}

if (websiteLive) enriched[["GO_Biological_Process_2023"]]
pdf("Graphs/Main/UP_Enrichment_Perilo_Interlo_medFB.pdf")
if (websiteLive) {
  plotEnrich(enriched[[3]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")
}
dev.off()

#GO-Terms_medFB_MHCIIh
table(FB_sub.clustermarker_UP$cluster)
medFB_MHCIIh_sub.clustermarker_UP <- subset(FB_sub.clustermarker_UP, cluster=="medFB_MHCIIh")
medFB_MHCIIh_sub.clustermarker_UP <- as.character(medFB_MHCIIh_sub.clustermarker_UP$gene)

if (websiteLive) {
  enriched <- enrichr(medFB_MHCIIh_sub.clustermarker_UP, dbs)
}

if (websiteLive) enriched[["GO_Biological_Process_2023"]]
pdf("Graphs/Main/UP_Enrichment_medFB_MHCIIh.pdf")
if (websiteLive) {
  plotEnrich(enriched[[3]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")
}
dev.off()

#GO-Terms_medFB_CCL2
table(FB_sub.clustermarker_UP$cluster)
medFB_CCL2_sub.clustermarker_UP <- subset(FB_sub.clustermarker_UP, cluster=="medFB_CCL2")
medFB_CCL2_sub.clustermarker_UP <- as.character(medFB_CCL2_sub.clustermarker_UP$gene)

if (websiteLive) {
  enriched <- enrichr(medFB_CCL2_sub.clustermarker_UP, dbs)
}

if (websiteLive) enriched[["GO_Biological_Process_2023"]]
pdf("Graphs/Main/UP_Enrichment_medFB_CCL2.pdf")
if (websiteLive) {
  plotEnrich(enriched[[3]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")
}
dev.off()

#GO-Terms_TET_FB
table(FB_sub.clustermarker_UP$cluster)
TET_FB_sub.clustermarker_UP <- subset(FB_sub.clustermarker_UP, cluster=="TET_FB")
TET_FB_sub.clustermarker_UP <- as.character(TET_FB_sub.clustermarker_UP$gene)

if (websiteLive) {
  enriched <- enrichr(TET_FB_sub.clustermarker_UP, dbs)
}

if (websiteLive) enriched[["GO_Biological_Process_2023"]]
pdf("Graphs/Main/UP_Enrichment_TET_FB.pdf")
if (websiteLive) {
  plotEnrich(enriched[[3]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")
}
dev.off()

#GO-Terms_fetal_FB_POSTN
table(FB_sub.clustermarker_UP$cluster)
fetal_FB_POSTN_sub.clustermarker_UP <- subset(FB_sub.clustermarker_UP, cluster=="fetal_FB_POSTN")
fetal_FB_POSTN_sub.clustermarker_UP <- as.character(fetal_FB_POSTN_sub.clustermarker_UP$gene)

if (websiteLive) {
  enriched <- enrichr(fetal_FB_POSTN_sub.clustermarker_UP, dbs)
}

if (websiteLive) enriched[["GO_Biological_Process_2023"]]
pdf("Graphs/Main/UP_Enrichment_fetal_FB_POSTN.pdf")
if (websiteLive) {
  plotEnrich(enriched[[3]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")
}
dev.off()

#GO-Terms_fetal_FB_SFRP1 
table(FB_sub.clustermarker_UP$cluster)
fetal_FB_SFRP1_sub.clustermarker_UP <- subset(FB_sub.clustermarker_UP, cluster=="fetal_FB_SFRP1")
fetal_FB_SFRP1_sub.clustermarker_UP <- as.character(fetal_FB_SFRP1_sub.clustermarker_UP$gene)

if (websiteLive) {
  enriched <- enrichr(fetal_FB_SFRP1_sub.clustermarker_UP, dbs)
}

if (websiteLive) enriched[["GO_Biological_Process_2023"]]
pdf("Graphs/Main/UP_Enrichment_fetal_FB_SFRP1.pdf")
if (websiteLive) {
  plotEnrich(enriched[[3]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")
}
dev.off()

save(FB_sub,file="FB_sub.RData")


#(23)_Subset BC#####

Idents(th.total)<- th.total$celltype.full
BC_sub <- subset(th.total,idents = "BC")

DefaultAssay(BC_sub)<-"RNA"
BC_sub[["RNA"]]$data <- as(BC_sub[["RNA"]]$data, Class = "dgCMatrix")

BC_sub <- FindVariableFeatures(BC_sub)
BC_sub <- ScaleData(BC_sub)
BC_sub <- RunPCA(BC_sub)
ElbowPlot(BC_sub, ndims=20)
BC_sub <- RunUMAP(BC_sub, dims = 1:10)
BC_sub <- FindNeighbors(BC_sub, dims = 1:10)
BC_sub <- FindClusters(BC_sub, resolution = 0.1)
DimPlot(BC_sub, label = T, label.size = 3, alpha = 1)
DimPlot(BC_sub, label=T, split.by="condition", alpha = 1)
DimPlot(BC_sub, label=T, split.by="development", alpha = 1)
DimPlot(BC_sub, label=T, group.by="seurat_clusters", alpha = 1)
DefaultAssay(BC_sub)<- "RNA"


module_BC_XIN_signature <- c("naive_BC_DEG_yasumizu module score1","pre_GC_BC_DEG_yasumizu module score1","GC_BC_DEG_yasumizu module score1","PB_DEG_yasumizu module score1", "unswi_mem_BC_DEG_yasumizu module score1","thym_mem_BC_DEG_yasumizu module score1","mem_BC_I_DEG_yasumizu module score1","mem_BC_II_DEG_yasumizu module score1")

#Modulescore imaging
pdf("Graphs/Main/Modulescore_BC_Clusteridentification.pdf", width=10, height=10)
DotPlot(object = BC_sub, features = module_BC_XIN_signature,assay="RNA") + coord_flip()+ scale_color_viridis(option = "H")+theme(axis.text.x = element_text(angle = 45)) + ggtitle("BC_modulescore") + scale_y_discrete(guide = guide_axis(angle = 90))
DotPlot(object = BC_sub, features = c("FCER2","IL4R","STMN1","TYMS","BCL6","MEF2B","CDC1","FCRL2","XIST","CD27","JCHAIN","PLD4","S100A10","ANXA2"),assay="RNA") + coord_flip()+ scale_color_viridis(option = "H")+theme(axis.text.x = element_text(angle = 45)) + ggtitle("BC_modulescore") + scale_y_discrete(guide = guide_axis(angle = 90))
DotPlot(object = BC_sub, features = c("CD19","MS4A1","CD27","IGHD","IGHM","MME"),assay="RNA") + coord_flip()+ scale_color_viridis(option = "H")+theme(axis.text.x = element_text(angle = 45)) + ggtitle("BC_modulescore") + scale_y_discrete(guide = guide_axis(angle = 90))
VlnPlot(BC_sub,features=c("CD19","MS4A1","CD27","IGHD","IGHM","MME"),pt.size=0, cols=color_levelBC)
DotPlot(object = BC_sub, features = yayon_BC,assay="RNA") + coord_flip()+ scale_color_viridis(option = "H")+theme(axis.text.x = element_text(angle = 45)) + ggtitle("BC_modulescore") + scale_y_discrete(guide = guide_axis(angle = 90))
DotPlot(object = BC_sub, features = c("TNFSF9","CD70","LMNA","NR4A2","CDKN1A","CD2","ELOVL4","CD1E","CD8A","AQP3","RGS13","MEF2B","SUGCT","PDGFD",
                                      "AC023590.1","VPREB1","IGLL1","DEPP1","AKAP12","SMAD1","IGHG2","PTPRJ","CSGALNACT1","ST14","NAPSA","KCNG1","CNTNAP2","GCNT1","IL4R",
                                      "ADARB1","CLEC2B","CHMP1B","BORCS5","SIPA1L1","SNX9","DONSON","ZNF90","IGLC3","IGLC2","DNTT"), assay="RNA") + coord_flip()+ scale_color_viridis(option = "H")+theme(axis.text.x = element_text(angle = 45)) + ggtitle("BC_modulescore") + scale_y_discrete(guide = guide_axis(angle = 90))
dev.off()

#Modulescore_PMID:33815362
#Transitional (CD19+IgD+CD27-CD10+), Naive (CD19+IgD+CD27-CD10-), IgM Memory (CD19+IgD+CD27+), Classical Memory (CD19+IgD-CD27+) and Double Negative (CD19+IgD-CD27-) 
DotPlot(BC_sub,features=c("FCER2","IL4R","STMN1","TYMS","BCL6","MEF2B","CDC1","FCRL2","XIST","CD27","JCHAIN","PLD4","S100A10","ANXA2"))

#Name Cluster
Idents(BC_sub)<-BC_sub$seurat_clusters
BC_sub<-RenameIdents(BC_sub,
                     `0`="Naive_BC_1", `1`="classical_memory_BC_1", `2`="IgM_memory_BC", `3`="Naive_BC_2",`4`="Naive_BC_3",`5`="classical_memory_BC_3",`6`="DN_BC_1",`7`="classical_memory_BC_2",`8`="DN_BC_2")
BC_sub$celltype.full<-Idents(BC_sub)
table(BC_sub$celltype.full)

#Doublet and contamination check
DimPlot(BC_sub, label = TRUE, label.size = 6)
FeaturePlot(BC_sub,c("KRT19","CD3D"), blend=T, order=T, label=T, cols= c("lightgrey","green","blue"), max.cutoff = 2)

#Define cluster levels
clusters_ordered<-c("Naive_BC_1","Naive_BC_2","Naive_BC_3","classical_memory_BC_1","classical_memory_BC_2","classical_memory_BC_3","IgM_memory_BC","DN_BC_1","DN_BC_2")
BC_sub$celltype.full<- factor(BC_sub$celltype.full, levels = clusters_ordered)
Idents(BC_sub)<-BC_sub$celltype.full

pdf("Graphs/Main/BC_Modulescores_Vlnplot.pdf")
VlnPlot(BC_sub,features=c("CD19","MS4A1","CD27","IGHD","IGHM","MME"),pt.size=0, cols=color_levelBC)
dev.off()

DefaultAssay(BC_sub)<-"RNA"
Idents(BC_sub)<-"celltype.full"

pdf("Graphs/Main/BC_UMAPPlots.pdf", width=5, height=5)
DimPlot(BC_sub, label=T, label.color = "black",label.size = 5, repel=T, alpha=1) + ggtitle("BC-subset")+ scale_color_manual(values=color_levelBC) + NoLegend()
dev.off()


pdf("Graphs/Main/UMAPPlots_BC_sub_QC.pdf", width=5, height=5)
DimPlot(BC_sub, group.by="group")+ ggtitle("Basic-UMAP")
DimPlot(BC_sub, group.by="tissue")+ ggtitle("Basic-UMAP")
dev.off()

pdf("Graphs/Main/BC_Barplots.pdf")
plot_integrated_clusters_group (BC_sub)
plot_integrated_clusters_Donor_ID (BC_sub)
plot_integrated_celltypes_condition (BC_sub)
plot_integrated_celltypes_classification (BC_sub)
plot_integrated_celltypes_tissue (BC_sub)
plot_integrated_celltypes_sex (BC_sub)
plot_integrated_celltypes_age (BC_sub) 
plot_integrated_celltypes_development (BC_sub)
plot_integrated_celltypes_donor (BC_sub)
plot_BC_integrated_condition_celltypes(BC_sub)
plot_BC_integrated_development_celltypes(BC_sub)
plot_BC_integrated_group_clusters(BC_sub)
plot_BC_integrated_tissue_celltypes(BC_sub)
dev.off()

#add identity to analyze tissue/health diffgenes
#Set tissue in identy 
Idents(BC_sub) <- BC_sub$celltype.full
BC_sub$celltype.condition <- paste(Idents(BC_sub), BC_sub$condition, sep = "_")
Idents(BC_sub)<-BC_sub$celltype.condition
table(BC_sub$celltype.condition)

#Clustermarker
Idents(BC_sub)<-"celltype.full"
BC_sub.clustermarker<-FindAllMarkers(BC_sub, assay = "RNA", min.pct = 0.25, logfc.threshold = 0.25)
BC_sub.clustermarker$Foldchange_UP <- 2^(BC_sub.clustermarker$avg_log2FC)
BC_sub.clustermarker$Foldchange_DOWN <- 2^(-BC_sub.clustermarker$avg_log2FC)
BC_sub.clustermarker$Ratio_pct1_pct2 <- (BC_sub.clustermarker$pct.1)/(BC_sub.clustermarker$pct.2)
write.xlsx(BC_sub.clustermarker, "Lists/Clustermarker_BC_sub_clustermarker.xlsx")

#DotPlot Top up and down regulated genes 
#up
BC_sub.clustermarker_UP <- subset(BC_sub.clustermarker, Foldchange_UP >=2 )
top10_UP_clustermarker <- BC_sub.clustermarker_UP %>% group_by(cluster) %>% top_n(n=10,wt=Foldchange_UP)
top10_UP_clustermarker<-arrange(top10_UP_clustermarker,cluster,Foldchange_UP)
top10_UP_clustermarker <- as.character(top10_UP_clustermarker$gene)
top10_UP_clustermarker<-top10_UP_clustermarker[!duplicated(top10_UP_clustermarker)]

pdf("Graphs/Main/DotPlot_Top_10UP_Clustermarker_BC_FCmin2.pdf", width = 7, height = 20)
DotPlot(BC_sub, features=top10_UP_clustermarker, group.by = "celltype.full",assay="RNA") + coord_flip()+ scale_color_viridis(option = "H")+theme(axis.text.x = element_text(angle = 45)) + ggtitle("BC Cluster_DEG") + scale_y_discrete(guide = guide_axis(angle = 90))
dev.off()

#Heatmap#
## Define Gene lists that you want to plot as Heatmap
BC_sub.clustermarker_UP<-arrange(BC_sub.clustermarker_UP,cluster,Foldchange_UP)
genes <- as.character(BC_sub.clustermarker_UP$gene)
avgexp <- BC_sub

## Plot Heatmaps
## Calculate Average gene expressions (pseudobulked) for each of the celltypes in each condition
levels(avgexp)
DefaultAssay(avgexp) <- "RNA"
avgexp <- NormalizeData(avgexp)
avgexp <- AggregateExpression(avgexp, return.seurat =T, assays = "RNA")


# Modify the heatmap plot
library(RColorBrewer)
heatmap_plot <- DoHeatmap(
  avgexp,
  features = genes,
  size = 6,
  group.bar.height = 0.02,
  draw.lines = FALSE,angle = 90,
  group.colors = color_levelBC
) +
  scale_fill_gradientn(colours = rev(brewer.pal(8, name = "RdYlBu")))

# Modify the theme to make rownames cursive
heatmap_plot 
ggsave("Graphs/Main/Heatmap_BC_clustermarker.pdf", width = 40, height = 20, units = c("cm"))



#Enrichment analysis BC_sub
#GO-Terms_Naive_BC_1
table(BC_sub.clustermarker_UP$cluster)
Naive_BC_1_sub.clustermarker_UP <- subset(BC_sub.clustermarker_UP, cluster=="Naive_BC_1")
Naive_BC_1_sub.clustermarker_UP <- as.character(Naive_BC_1_sub.clustermarker_UP$gene)

if (websiteLive) {
  enriched <- enrichr(Naive_BC_1_sub.clustermarker_UP, dbs)
}

if (websiteLive) enriched[["GO_Biological_Process_2023"]]
pdf("Graphs/Main/UP_Enrichment_Naive_BC_1.pdf")
if (websiteLive) {
  plotEnrich(enriched[[3]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")
}
dev.off()

#GO-Terms_Naive_BC_2
table(BC_sub.clustermarker_UP$cluster)
Naive_BC_2_sub.clustermarker_UP <- subset(BC_sub.clustermarker_UP, cluster=="Naive_BC_2")
Naive_BC_2_sub.clustermarker_UP <- as.character(Naive_BC_2_sub.clustermarker_UP$gene)

if (websiteLive) {
  enriched <- enrichr(Naive_BC_2_sub.clustermarker_UP, dbs)
}

if (websiteLive) enriched[["GO_Biological_Process_2023"]]
pdf("Graphs/Main/UP_Enrichment_Naive_BC_2.pdf")
if (websiteLive) {
  plotEnrich(enriched[[3]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")
}
dev.off()

#GO-Terms_Naive_BC_3
table(BC_sub.clustermarker_UP$cluster)
Naive_BC_3_sub.clustermarker_UP <- subset(BC_sub.clustermarker_UP, cluster=="Naive_BC_3")
Naive_BC_3_sub.clustermarker_UP <- as.character(Naive_BC_3_sub.clustermarker_UP$gene)

if (websiteLive) {
  enriched <- enrichr(Naive_BC_3_sub.clustermarker_UP, dbs)
}

if (websiteLive) enriched[["GO_Biological_Process_2023"]]
pdf("Graphs/Main/UP_Enrichment_Naive_BC_3.pdf")
if (websiteLive) {
  plotEnrich(enriched[[3]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")
}
dev.off()

#GO-Terms_DN_BC_2
table(BC_sub.clustermarker_UP$cluster)
DN_BC_2_sub.clustermarker_UP <- subset(BC_sub.clustermarker_UP, cluster=="DN_BC_2")
DN_BC_2_sub.clustermarker_UP <- as.character(DN_BC_2_sub.clustermarker_UP$gene)

if (websiteLive) {
  enriched <- enrichr(DN_BC_2_sub.clustermarker_UP, dbs)
}

if (websiteLive) enriched[["GO_Biological_Process_2023"]]
pdf("Graphs/Main/UP_Enrichment_DN_BC_2.pdf")
if (websiteLive) {
  plotEnrich(enriched[[3]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")
}
dev.off()

#GO-Terms_classical_memory_BC_1
table(BC_sub.clustermarker_UP$cluster)
classical_memory_BC_1_sub.clustermarker_UP <- subset(BC_sub.clustermarker_UP, cluster=="classical_memory_BC_1")
classical_memory_BC_1_sub.clustermarker_UP <- as.character(classical_memory_BC_1_sub.clustermarker_UP$gene)

if (websiteLive) {
  enriched <- enrichr(classical_memory_BC_1_sub.clustermarker_UP, dbs)
}

if (websiteLive) enriched[["GO_Biological_Process_2023"]]
pdf("Graphs/Main/UP_Enrichment_classical_memory_BC_1.pdf")
if (websiteLive) {
  plotEnrich(enriched[[3]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")
}
dev.off()

#GO-Terms_classical_memory_BC_2
table(BC_sub.clustermarker_UP$cluster)
classical_memory_BC_2_sub.clustermarker_UP <- subset(BC_sub.clustermarker_UP, cluster=="classical_memory_BC_2")
classical_memory_BC_2_sub.clustermarker_UP <- as.character(classical_memory_BC_2_sub.clustermarker_UP$gene)

if (websiteLive) {
  enriched <- enrichr(classical_memory_BC_2_sub.clustermarker_UP, dbs)
}

if (websiteLive) enriched[["GO_Biological_Process_2023"]]
pdf("Graphs/Main/UP_Enrichment_classical_memory_BC_2.pdf")
if (websiteLive) {
  plotEnrich(enriched[[3]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")
}
dev.off()

#GO-Terms_classical_memory_BC_3
table(BC_sub.clustermarker_UP$cluster)
classical_memory_BC_3_sub.clustermarker_UP <- subset(BC_sub.clustermarker_UP, cluster=="classical_memory_BC_3")
classical_memory_BC_3_sub.clustermarker_UP <- as.character(classical_memory_BC_3_sub.clustermarker_UP$gene)

if (websiteLive) {
  enriched <- enrichr(classical_memory_BC_3_sub.clustermarker_UP, dbs)
}

if (websiteLive) enriched[["GO_Biological_Process_2023"]]
pdf("Graphs/Main/UP_Enrichment_classical_memory_BC_3.pdf")
if (websiteLive) {
  plotEnrich(enriched[[3]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")
}
dev.off()

#GO-Terms_IgM_memory_BC
table(BC_sub.clustermarker_UP$cluster)
IgM_memory_BC_sub.clustermarker_UP <- subset(BC_sub.clustermarker_UP, cluster=="IgM_memory_BC")
IgM_memory_BC_sub.clustermarker_UP <- as.character(IgM_memory_BC_sub.clustermarker_UP$gene)

if (websiteLive) {
  enriched <- enrichr(IgM_memory_BC_sub.clustermarker_UP, dbs)
}

if (websiteLive) enriched[["GO_Biological_Process_2023"]]
pdf("Graphs/Main/UP_Enrichment_IgM_memory_BC.pdf")
if (websiteLive) {
  plotEnrich(enriched[[3]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")
}
dev.off()


#GO-Terms_DN_BC_1
table(BC_sub.clustermarker_UP$cluster)
DN_BC_1_sub.clustermarker_UP <- subset(BC_sub.clustermarker_UP, cluster=="DN_BC_1")
DN_BC_1_sub.clustermarker_UP <- as.character(DN_BC_1_sub.clustermarker_UP$gene)

if (websiteLive) {
  enriched <- enrichr(DN_BC_1_sub.clustermarker_UP, dbs)
}

if (websiteLive) enriched[["GO_Biological_Process_2023"]]
pdf("Graphs/Main/UP_Enrichment_DN_BC_1.pdf")
if (websiteLive) {
  plotEnrich(enriched[[3]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")
}
dev.off()

save(BC_sub,file="BC_sub.RData")

#(24)_Subset PC####
Idents(th.total)<- th.total$celltype.full
PC_sub <- subset(th.total,idents = "PC")

DefaultAssay(PC_sub)<-"RNA"
PC_sub[["RNA"]]$data <- as(PC_sub[["RNA"]]$data, Class = "dgCMatrix")

PC_sub <- FindVariableFeatures(PC_sub)
PC_sub <- ScaleData(PC_sub)
PC_sub <- RunPCA(PC_sub)
ElbowPlot(PC_sub, ndims=20)
PC_sub <- RunUMAP(PC_sub, dims = 1:8)
PC_sub <- FindNeighbors(PC_sub, dims = 1:8)
PC_sub <- FindClusters(PC_sub, resolution = 0.1)
DimPlot(PC_sub, label = T, label.size = 3, alpha = 1)
DimPlot(PC_sub, label=T, split.by="condition", alpha = 1)
DimPlot(PC_sub, label=T, split.by="development", alpha = 1)
DimPlot(PC_sub, label=T, split.by="group", alpha = 1)
DimPlot(PC_sub, label=T, group.by="seurat_clusters", alpha = 1)
DefaultAssay(PC_sub)<- "RNA"

#Modulescore imaging  PMID:36098652
pdf("Graphs/Main/Modulescore_PC_Clusteridentification.pdf")
DotPlot(object = PC_sub, features = c("IGHM","IGHD","IGHG1","IGHG2","IGHG3","IGHG4","IGHA1","IGHA2"),assay="RNA") + coord_flip()+ scale_color_viridis(option = "H")+theme(axis.text.x = element_text(angle = 45)) + ggtitle("PC_modulescore") + scale_y_discrete(guide = guide_axis(angle = 90))
DotPlot(object = PC_sub, features = c("IGHA1", "MZB1", "XBP1"),assay="RNA") + coord_flip()+ scale_color_viridis(option = "H")+theme(axis.text.x = element_text(angle = 45)) + ggtitle("PC_modulescore") + scale_y_discrete(guide = guide_axis(angle = 90))
DotPlot(object = PC_sub, features = c("IGLC1","STMN1","CD52","HLA-DRA","HLA-DPB1","SOX5","TEX14","LINC00910","TRIO","BFSP2","MTRNR2L8",
                                      "MTRNR2L12","PTP4A1","SMOC1","ABCB9","DLX5","CCL19","FABP4","SPARCL1","TPM2"),assay="RNA") + coord_flip()+ scale_color_viridis(option = "H")+theme(axis.text.x = element_text(angle = 45)) + ggtitle("PC_modulescore") + scale_y_discrete(guide = guide_axis(angle = 90))
FeaturePlot(object = PC_sub, features = c("IGHA1", "MZB1", "XBP1"), order=T)
dev.off()

#Name Cluster
Idents(PC_sub)<-PC_sub$seurat_clusters
PC_sub<-RenameIdents(PC_sub,
                     `0`="PC2", `1`="PC1", `2`="PC4", `3`="PC5", `4`="PC6", `5`="PC3")
PC_sub$celltype.full<-Idents(PC_sub)

#Doublet and contamination check
DimPlot(PC_sub, label = TRUE, label.size = 6)
FeaturePlot(PC_sub,c("KRT19","CD3D"), blend=T, order=T, label=T, cols= c("lightgrey","green","blue"), max.cutoff = 2)

#Define cluster levels
clusters_ordered<-c("PC1","PC2","PC3","PC4","PC5","PC6")
PC_sub$celltype.full<- factor(PC_sub$celltype.full, levels = clusters_ordered)
Idents(PC_sub)<-PC_sub$celltype.full

DefaultAssay(PC_sub)<-"RNA"
Idents(PC_sub)<-"celltype.full"

pdf("Graphs/Main/PC_UMAPPlots.pdf", width=5, height=5)
DimPlot(PC_sub, label=T, label.color = "black",label.size = 5, repel=T, alpha=1) + ggtitle("PC-subset")+ scale_color_manual(values=color_levelPC) + NoLegend()
dev.off()

pdf("Graphs/Main/UMAPPlots_PC_sub_QC.pdf", width=5, height=5)
DimPlot(PC_sub, group.by="group")+ ggtitle("Basic-UMAP")
DimPlot(PC_sub, group.by="tissue")+ ggtitle("Basic-UMAP")
dev.off()


pdf("Graphs/Main/PC_Barplots.pdf")
plot_integrated_clusters_group (PC_sub)
plot_integrated_clusters_Donor_ID (PC_sub)
plot_integrated_celltypes_condition (PC_sub)
plot_integrated_celltypes_classification (PC_sub)
plot_integrated_celltypes_tissue (PC_sub)
plot_integrated_celltypes_sex (PC_sub)
plot_integrated_celltypes_age (PC_sub) 
plot_integrated_celltypes_development (PC_sub)
plot_integrated_celltypes_donor (PC_sub)
plot_PC_integrated_condition_celltypes(PC_sub)
plot_PC_integrated_development_celltypes(PC_sub)
plot_PC_integrated_group_clusters(PC_sub)
plot_PC_integrated_tissue_celltypes(PC_sub)
dev.off()


pdf("Graphs/Main/PC_sub_PC4_Clusteridentification.pdf", width=10, height=5)
FeaturePlot(PC_sub,c("MZB1","SPARCL1"), blend=T, order=T, cols= c("lightgrey","red","blue"))
dev.off()

#add identity to analyze tissue/health diffgenes
#Set tissue in identy 
Idents(PC_sub) <- PC_sub$celltype.full
PC_sub$celltype.condition <- paste(Idents(PC_sub), PC_sub$condition, sep = "_")
Idents(PC_sub)<-PC_sub$celltype.condition
table(PC_sub$celltype.condition)
#Set tissue in identy 
Idents(PC_sub) <- PC_sub$celltype.condition
PC_sub$celltype.condition.development <- paste(Idents(PC_sub), PC_sub$development, sep = "_")
Idents(PC_sub)<-PC_sub$celltype.condition.development
table(PC_sub$celltype.condition.development)

#Clustermarker
Idents(PC_sub)<-"celltype.full"
PC_sub.clustermarker<-FindAllMarkers(PC_sub, assay = "RNA", min.pct = 0.25, logfc.threshold = 0.25)
PC_sub.clustermarker$Foldchange_UP <- 2^(PC_sub.clustermarker$avg_log2FC)
PC_sub.clustermarker$Foldchange_DOWN <- 2^(-PC_sub.clustermarker$avg_log2FC)
PC_sub.clustermarker$Ratio_pct1_pct2 <- (PC_sub.clustermarker$pct.1)/(PC_sub.clustermarker$pct.2)
write.xlsx(PC_sub.clustermarker, "Lists/Clustermarker_PC_sub_clustermarker.xlsx")

#DotPlot Top up and down regulated genes 
#up
PC_sub.clustermarker_UP <- subset(PC_sub.clustermarker, Foldchange_UP >=2 )
top10_UP_clustermarker <- PC_sub.clustermarker_UP %>% group_by(cluster) %>% top_n(n=10,wt=Foldchange_UP)
top10_UP_clustermarker<-arrange(top10_UP_clustermarker,cluster,Foldchange_UP)
top10_UP_clustermarker <- as.character(top10_UP_clustermarker$gene)
top10_UP_clustermarker<-top10_UP_clustermarker[!duplicated(top10_UP_clustermarker)]

pdf("Graphs/Main/DotPlot_Top_10UP_Clustermarker_PC_FCmin2.pdf", width = 7, height = 12)
DotPlot(PC_sub, features=top10_UP_clustermarker, group.by = "celltype.full",assay="RNA") + coord_flip()+ scale_color_viridis(option = "H")+theme(axis.text.x = element_text(angle = 45)) + ggtitle("PC Cluster_DEG") + scale_y_discrete(guide = guide_axis(angle = 90))
dev.off()

#Heatmap#
## Define Gene lists that you want to plot as Heatmap
PC_sub.clustermarker_UP<-arrange(PC_sub.clustermarker_UP,cluster,Foldchange_UP)
genes <- as.character(PC_sub.clustermarker_UP$gene)
avgexp <- PC_sub

## Plot Heatmaps
## Calculate Average gene expressions (pseudobulked) for each of the celltypes in each condition
levels(avgexp)
DefaultAssay(avgexp) <- "RNA"
avgexp <- NormalizeData(avgexp)
avgexp <- AggregateExpression(avgexp, return.seurat =T, assays = "RNA")


# Modify the heatmap plot
library(RColorBrewer)
heatmap_plot <- DoHeatmap(
  avgexp,
  features = genes,
  size = 6,
  group.bar.height = 0.02,
  draw.lines = FALSE,angle = 90,
  group.colors = color_levelPC
) +
  scale_fill_gradientn(colours = rev(brewer.pal(8, name = "RdYlBu")))

# Modify the theme to make rownames cursive
heatmap_plot 
ggsave("Graphs/Main/Heatmap_PC_clustermarker.pdf", width = 40, height = 20, units = c("cm"))

#Enrichment analysis PC_sub
#GO-Terms_PC1
table(PC_sub.clustermarker_UP$cluster)
PC1_sub.clustermarker_UP <- subset(PC_sub.clustermarker_UP, cluster=="PC1")
PC1_sub.clustermarker_UP <- as.character(PC1_sub.clustermarker_UP$gene)

if (websiteLive) {
  enriched <- enrichr(PC1_sub.clustermarker_UP, dbs)
}

if (websiteLive) enriched[["GO_Biological_Process_2023"]]
pdf("Graphs/Main/UP_Enrichment_PC1.pdf")
if (websiteLive) {
  plotEnrich(enriched[[3]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")
}
dev.off()

#GO-Terms_PC2
table(PC_sub.clustermarker_UP$cluster)
PC2_sub.clustermarker_UP <- subset(PC_sub.clustermarker_UP, cluster=="PC2")
PC2_sub.clustermarker_UP <- as.character(PC2_sub.clustermarker_UP$gene)

if (websiteLive) {
  enriched <- enrichr(PC2_sub.clustermarker_UP, dbs)
}

if (websiteLive) enriched[["GO_Biological_Process_2023"]]
pdf("Graphs/Main/UP_Enrichment_PC2.pdf")
if (websiteLive) {
  plotEnrich(enriched[[3]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")
}
dev.off()

#GO-Terms_PC3
table(PC_sub.clustermarker_UP$cluster)
PC3_sub.clustermarker_UP <- subset(PC_sub.clustermarker_UP, cluster=="PC3")
PC3_sub.clustermarker_UP <- as.character(PC3_sub.clustermarker_UP$gene)

if (websiteLive) {
  enriched <- enrichr(PC3_sub.clustermarker_UP, dbs)
}

if (websiteLive) enriched[["GO_Biological_Process_2023"]]
pdf("Graphs/Main/UP_Enrichment_PC3.pdf")
if (websiteLive) {
  plotEnrich(enriched[[3]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")
}
dev.off()

#GO-Terms_PC4
table(PC_sub.clustermarker_UP$cluster)
PC4_sub.clustermarker_UP <- subset(PC_sub.clustermarker_UP, cluster=="PC4")
PC4_sub.clustermarker_UP <- as.character(PC4_sub.clustermarker_UP$gene)

if (websiteLive) {
  enriched <- enrichr(PC4_sub.clustermarker_UP, dbs)
}

if (websiteLive) enriched[["GO_Biological_Process_2023"]]
pdf("Graphs/Main/UP_Enrichment_PC4.pdf")
if (websiteLive) {
  plotEnrich(enriched[[3]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")
}
dev.off()

#GO-Terms_PC5
table(PC_sub.clustermarker_UP$cluster)
PC5_sub.clustermarker_UP <- subset(PC_sub.clustermarker_UP, cluster=="PC5")
PC5_sub.clustermarker_UP <- as.character(PC5_sub.clustermarker_UP$gene)

if (websiteLive) {
  enriched <- enrichr(PC5_sub.clustermarker_UP, dbs)
}

if (websiteLive) enriched[["GO_Biological_Process_2023"]]
pdf("Graphs/Main/UP_Enrichment_PC5.pdf")
if (websiteLive) {
  plotEnrich(enriched[[3]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")
}
dev.off()

#GO-Terms_PC6
table(PC_sub.clustermarker_UP$cluster)
PC6_sub.clustermarker_UP <- subset(PC_sub.clustermarker_UP, cluster=="PC6")
PC6_sub.clustermarker_UP <- as.character(PC6_sub.clustermarker_UP$gene)

if (websiteLive) {
  enriched <- enrichr(PC6_sub.clustermarker_UP, dbs)
}

if (websiteLive) enriched[["GO_Biological_Process_2023"]]
pdf("Graphs/Main/UP_Enrichment_PC6.pdf")
if (websiteLive) {
  plotEnrich(enriched[[3]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")
}
dev.off()


save(PC_sub,file="PC_sub.RData")

#(25)_Subset DC####
Idents(th.total)<- th.total$celltype.full
DC_sub <- subset(th.total,idents = "DC")

DefaultAssay(DC_sub)<-"RNA"
DC_sub[["RNA"]]$data <- as(DC_sub[["RNA"]]$data, Class = "dgCMatrix")

DC_sub <- FindVariableFeatures(DC_sub)
DC_sub <- ScaleData(DC_sub)
DC_sub <- RunPCA(DC_sub)
ElbowPlot(DC_sub, ndims=20)
DC_sub <- RunUMAP(DC_sub, dims = 1:7)
DC_sub <- FindNeighbors(DC_sub, dims = 1:7)
DC_sub <- FindClusters(DC_sub, resolution = 0.05)
DimPlot(DC_sub, label = T, label.size = 3, alpha = 1)
DimPlot(DC_sub, label=T, split.by="condition", alpha = 1)
DimPlot(DC_sub, label=T, split.by="development", alpha = 1)
DimPlot(DC_sub, label=T, split.by="group", alpha = 1)
DimPlot(DC_sub, label=T, group.by="seurat_clusters", alpha = 1)
DefaultAssay(DC_sub)<- "RNA"

#Modulscore:
module_DC_DEG_yasumizu_signature <- c("cDC1_DEG_yasumizu module score1","cDC2_DEG_yasumizu module score1","pDC_DEG_yasumizu module score1")

#Modulescore imaging + Genelist Park
pdf("Graphs/Main/Modulescore_DC_Clusteridentification.pdf")
DotPlot(object = DC_sub, features = module_DC_DEG_yasumizu_signature, assay="RNA") + coord_flip()+ scale_color_viridis(option = "H")+theme(axis.text.x = element_text(angle = 45)) + ggtitle("DC_modulescore") + scale_y_discrete(guide = guide_axis(angle = 90))
DotPlot(object = DC_sub, features = c("XCR1","CLEC9A","SIRPA","CLEC10A","IL3RA","CLEC4C","LAMP3","CCR7"), assay="RNA") + coord_flip()+ scale_color_viridis(option = "H")+theme(axis.text.x = element_text(angle = 45)) + ggtitle("DC_modulescore") + scale_y_discrete(guide = guide_axis(angle = 90))
DotPlot(object = DC_sub, features = c("XCR1","CLEC9A","ANPEP","FBXO27","SIRPA","CLEC10A","DENND3","IL3RA","JCHAIN","CLEC4C","LAMP3","CCR7","AIRE","FOXD4","TNFRSF11A","ST7","CST7","TNFRSF11B","NEK6","MT2A","CXCL9","CXCL10","CCL22","CCL17","CCL19","CD40","CD80","CD274","CD74","HLA-DRA"), assay="RNA") + coord_flip()+ scale_color_viridis(option = "H")+theme(axis.text.x = element_text(angle = 45)) + ggtitle("DC_modulescore") + scale_y_discrete(guide = guide_axis(angle = 90))
DotPlot(object = DC_sub, features =yayon_DC_MAC, assay="RNA") + coord_flip()+ scale_color_viridis(option = "H")+theme(axis.text.x = element_text(angle = 45)) + ggtitle("DC_modulescore") + scale_y_discrete(guide = guide_axis(angle = 90))
DotPlot(object = DC_sub, features =c("ENTHD1","CCL17","EBI3","XCR1","CLEC9A","SLAMF8","ASIP","GZMB","LINC00996"), assay="RNA") + coord_flip()+ scale_color_viridis(option = "H")+theme(axis.text.x = element_text(angle = 45)) + ggtitle("DC_modulescore") + scale_y_discrete(guide = guide_axis(angle = 90))
dev.off()

#Name Cluster
Idents(DC_sub)<-DC_sub$seurat_clusters
DC_sub<-RenameIdents(DC_sub,
                     `0`="pDC1", `1`="pDC2", `2`="pDC3")
DC_sub$celltype.full<-Idents(DC_sub)

#Doublet and contamination check
DimPlot(DC_sub, label = TRUE, label.size = 6)
FeaturePlot(DC_sub,c("KRT19","CD3D"), blend=T, order=T, label=T, cols= c("lightgrey","green","blue"), max.cutoff = 2)

#Define cluster levels
clusters_ordered<-c("pDC1","pDC2","pDC3")
DC_sub$celltype.full<- factor(DC_sub$celltype.full, levels = clusters_ordered)
Idents(DC_sub)<-DC_sub$celltype.full

DefaultAssay(DC_sub)<-"RNA"
Idents(DC_sub)<-"celltype.full"

pdf("Graphs/Main/DC_UMAPPlots.pdf", width=5, height=5)
DimPlot(DC_sub, label=T, label.color = "black",label.size = 5, repel=T, alpha=1) + ggtitle("DC-subset")+ scale_color_manual(values=color_levelDC) + NoLegend()
dev.off()

pdf("Graphs/Main/UMAPPlots_DC_sub_QC.pdf", width=5, height=5)
DimPlot(DC_sub, group.by="group")+ ggtitle("Basic-UMAP")
DimPlot(DC_sub, group.by="tissue")+ ggtitle("Basic-UMAP")
dev.off()

pdf("Graphs/Main/DC_Barplots.pdf")
plot_integrated_clusters_group (DC_sub)
plot_integrated_clusters_Donor_ID (DC_sub)
plot_integrated_celltypes_condition (DC_sub)
plot_integrated_celltypes_classification (DC_sub)
plot_integrated_celltypes_tissue (DC_sub)
plot_integrated_celltypes_sex (DC_sub)
plot_integrated_celltypes_age (DC_sub) 
plot_integrated_celltypes_development (DC_sub)
plot_integrated_celltypes_donor (DC_sub)
plot_DC_integrated_condition_celltypes(DC_sub)
plot_DC_integrated_development_celltypes(DC_sub)
plot_DC_integrated_group_clusters(DC_sub)
plot_DC_integrated_tissue_celltypes(DC_sub)
dev.off()

#add identity to analyze tissue/health diffgenes
#Set tissue in identy 
Idents(DC_sub) <- DC_sub$celltype.full
DC_sub$celltype.condition <- paste(Idents(DC_sub), DC_sub$condition, sep = "_")
Idents(DC_sub)<-DC_sub$celltype.condition
table(DC_sub$celltype.condition)
#Set tissue in identy 
Idents(DC_sub) <- DC_sub$celltype.condition
DC_sub$celltype.condition.development <- paste(Idents(DC_sub), DC_sub$development, sep = "_")
Idents(DC_sub)<-DC_sub$celltype.condition.development
table(DC_sub$celltype.condition.development)

#Clustermarker
Idents(DC_sub)<-"celltype.full"
DC_sub.clustermarker<-FindAllMarkers(DC_sub, assay = "RNA", min.pct = 0.25, logfc.threshold = 0.25)
DC_sub.clustermarker$Foldchange_UP <- 2^(DC_sub.clustermarker$avg_log2FC)
DC_sub.clustermarker$Foldchange_DOWN <- 2^(-DC_sub.clustermarker$avg_log2FC)
DC_sub.clustermarker$Ratio_pct1_pct2 <- (DC_sub.clustermarker$pct.1)/(DC_sub.clustermarker$pct.2)
write.xlsx(DC_sub.clustermarker, "Lists/Clustermarker_DC_sub_clustermarker.xlsx")

#DotPlot Top up and down regulated genes 
#up
DC_sub.clustermarker_UP <- subset(DC_sub.clustermarker, Foldchange_UP >=2 )
top10_UP_clustermarker <- DC_sub.clustermarker_UP %>% group_by(cluster) %>% top_n(n=10,wt=Foldchange_UP)
top10_UP_clustermarker<-arrange(top10_UP_clustermarker,cluster,Foldchange_UP)
top10_UP_clustermarker <- as.character(top10_UP_clustermarker$gene)
top10_UP_clustermarker<-top10_UP_clustermarker[!duplicated(top10_UP_clustermarker)]

pdf("Graphs/Main/DotPlot_Top_10UP_Clustermarker_DC_FCmin2.pdf", width = 7, height = 12)
DotPlot(DC_sub, features=top10_UP_clustermarker, group.by = "celltype.full",assay="RNA") + coord_flip()+ scale_color_viridis(option = "H")+theme(axis.text.x = element_text(angle = 45)) + ggtitle("DC Cluster_DEG") + scale_y_discrete(guide = guide_axis(angle = 90))
dev.off()


#Heatmap#
## Define Gene lists that you want to plot as Heatmap
DC_sub.clustermarker_UP<-arrange(DC_sub.clustermarker_UP,cluster,Foldchange_UP)
genes <- as.character(DC_sub.clustermarker_UP$gene)
avgexp <- DC_sub

## Plot Heatmaps
## Calculate Average gene expressions (pseudobulked) for each of the celltypes in each condition
levels(avgexp)
DefaultAssay(avgexp) <- "RNA"
avgexp <- NormalizeData(avgexp)
avgexp <- AggregateExpression(avgexp, return.seurat =T, assays = "RNA")


# Modify the heatmap plot
library(RColorBrewer)
heatmap_plot <- DoHeatmap(
  avgexp,
  features = genes,
  size = 6,
  group.bar.height = 0.02,
  draw.lines = FALSE,angle = 90,
  group.colors = color_levelDC
) +
  scale_fill_gradientn(colours = rev(brewer.pal(8, name = "RdYlBu")))

# Modify the theme to make rownames cursive
heatmap_plot 
ggsave("Graphs/Main/Heatmap_DC_clustermarker.pdf", width = 40, height = 20, units = c("cm"))


#Clustermarker
Idents(DC_sub)<-"celltype.condition"
DC_cc_sub.clustermarker<-FindAllMarkers(DC_sub, assay = "RNA", min.pct = 0.25, logfc.threshold = 0.25)
DC_cc_sub.clustermarker$Foldchange_UP <- 2^(DC_cc_sub.clustermarker$avg_log2FC)
DC_cc_sub.clustermarker$Foldchange_DOWN <- 2^(-DC_cc_sub.clustermarker$avg_log2FC)
DC_cc_sub.clustermarker$Ratio_pct1_pct2 <- (DC_cc_sub.clustermarker$pct.1)/(DC_cc_sub.clustermarker$pct.2)
write.xlsx(DC_cc_sub.clustermarker, "Lists/Clustermarker_DC_sub_clustermarker_celltype_condition.xlsx")

#DotPlot Top up and down regulated genes 
#up
DC_cc_sub.clustermarker_UP <- subset(DC_cc_sub.clustermarker, Foldchange_UP >=2 )
top10_UP_clustermarker <- DC_cc_sub.clustermarker_UP %>% group_by(cluster) %>% top_n(n=10,wt=Foldchange_UP)
top10_UP_clustermarker<-arrange(top10_UP_clustermarker,cluster,Foldchange_UP)
top10_UP_clustermarker <- as.character(top10_UP_clustermarker$gene)
top10_UP_clustermarker<-top10_UP_clustermarker[!duplicated(top10_UP_clustermarker)]

pdf("Graphs/Main/DotPlot_Top_10UP_Clustermarker_DC_cc_FCmin2.pdf", width = 7, height = 12)
DotPlot(DC_sub, features=top10_UP_clustermarker, group.by = "celltype.condition",assay="RNA") + coord_flip()+ scale_color_viridis(option = "H")+theme(axis.text.x = element_text(angle = 45)) + ggtitle("DC Cluster_DEG") + scale_y_discrete(guide = guide_axis(angle = 90))
dev.off()


#Enrichment analysis DC_sub
#GO-Terms_pDC1
table(DC_sub.clustermarker_UP$cluster)
pDC1_sub.clustermarker_UP <- subset(DC_sub.clustermarker_UP, cluster=="pDC1")
pDC1_sub.clustermarker_UP <- as.character(pDC1_sub.clustermarker_UP$gene)

if (websiteLive) {
  enriched <- enrichr(pDC1_sub.clustermarker_UP, dbs)
}

if (websiteLive) enriched[["GO_Biological_Process_2023"]]
pdf("Graphs/Main/UP_Enrichment_pDC1.pdf")
if (websiteLive) {
  plotEnrich(enriched[[3]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")
}
dev.off()

#Enrichment analysis DC_sub
#GO-Terms_pDC2
table(DC_sub.clustermarker_UP$cluster)
pDC2_sub.clustermarker_UP <- subset(DC_sub.clustermarker_UP, cluster=="pDC2")
pDC2_sub.clustermarker_UP <- as.character(pDC2_sub.clustermarker_UP$gene)

if (websiteLive) {
  enriched <- enrichr(pDC2_sub.clustermarker_UP, dbs)
}

if (websiteLive) enriched[["GO_Biological_Process_2023"]]
pdf("Graphs/Main/UP_Enrichment_pDC2.pdf")
if (websiteLive) {
  plotEnrich(enriched[[3]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")
}
dev.off()

#Enrichment analysis DC_sub
#GO-Terms_pDC3
table(DC_sub.clustermarker_UP$cluster)
pDC3_sub.clustermarker_UP <- subset(DC_sub.clustermarker_UP, cluster=="pDC3")
pDC3_sub.clustermarker_UP <- as.character(pDC3_sub.clustermarker_UP$gene)

if (websiteLive) {
  enriched <- enrichr(pDC3_sub.clustermarker_UP, dbs)
}

if (websiteLive) enriched[["GO_Biological_Process_2023"]]
pdf("Graphs/Main/UP_Enrichment_pDC3.pdf")
if (websiteLive) {
  plotEnrich(enriched[[3]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")
}
dev.off()


save(DC_sub,file="DC_sub.RData")


#(26)_Subset MAC/Mono####
Idents(th.total)<- th.total$celltype.full
MAC_Mono_sub <- subset(th.total,idents = "MAC/Mono")

DefaultAssay(MAC_Mono_sub)<-"RNA"
MAC_Mono_sub[["RNA"]]$data <- as(MAC_Mono_sub[["RNA"]]$data, Class = "dgCMatrix")

MAC_Mono_sub <- FindVariableFeatures(MAC_Mono_sub)
MAC_Mono_sub <- ScaleData(MAC_Mono_sub)
MAC_Mono_sub <- RunPCA(MAC_Mono_sub)
ElbowPlot(MAC_Mono_sub, ndims=20)
MAC_Mono_sub <- RunUMAP(MAC_Mono_sub, dims = 1:10)
MAC_Mono_sub <- FindNeighbors(MAC_Mono_sub, dims = 1:10)
MAC_Mono_sub <- FindClusters(MAC_Mono_sub, resolution = 0.1)
DimPlot(MAC_Mono_sub, label = T, label.size = 3, alpha = 1)
DimPlot(MAC_Mono_sub, label=T, split.by="condition", alpha = 1)
DimPlot(MAC_Mono_sub, label=T, split.by="development", alpha = 1)
DimPlot(MAC_Mono_sub, label=T, split.by="group", alpha = 1)
DimPlot(MAC_Mono_sub, label=T, group.by="seurat_clusters", alpha = 1)
DefaultAssay(MAC_Mono_sub)<- "RNA"

#Modulscore:
module_MAC_Mono_DEG_yasumizu_subset_signature <- c("mono_CD14_DEG_yasumizu module score1","mono_CD16_DEG_yasumizu module score1","Mac_DEG_yasumizu module score1")

#Modulescore imaging
pdf("Graphs/Main/Modulescore_MAC_Mono_Clusteridentification.pdf")
DotPlot(object = MAC_Mono_sub, features = c(module_MAC_Mono_DEG_yasumizu_subset_signature,"CD8A","CD8B","CD3D","LUM","DCN","COL1A1"),assay="RNA") + coord_flip()+ scale_color_viridis(option = "H")+theme(axis.text.x = element_text(angle = 45)) + ggtitle("MAC_Mono_modulescore") + scale_y_discrete(guide = guide_axis(angle = 90))
DotPlot(MAC_Mono_sub,features=c("AIF1","IL1B","CXCL2","CD163","MRC1","MARCO","MSR1","CD8A","CD8B","CD3D","CD4","CEBPE","HDC","MS4A2", "LUM","DCN","COL1A1"))+coord_flip() + scale_color_viridis(option = "H")#+ theme(panel.background = element_rect(fill = 'lightgrey'))
FeaturePlot(MAC_Mono_sub,c("CD68","AIF1","IL1B","CXCL2","CD163","MRC1"), order=T, label=T)
DotPlot(MAC_Mono_sub,features=yayon_DC_MAC)+coord_flip() + scale_color_viridis(option = "H")#+ theme(panel.background = element_rect(fill = 'lightgrey'))
DotPlot(MAC_Mono_sub,features=c("ELOVL4","RAG2","CD8A","EGFL6","C7","HAS2","IFI27","APOC1","PLTP","EREG","G0S2","AREG","TEX14",
                                "CCL4L2","AC009313.1","PF4","PPBP","NRGN")
)+coord_flip() + scale_color_viridis(option = "H")#+ theme(panel.background = element_rect(fill = 'lightgrey'))
dev.off()

FeaturePlot(MAC_Mono_sub,c("CD68","AIF1","IL1B","CXCL2","CD163","MRC1"), order=T)
DotPlot(MAC_Mono_sub,features=c("AIF1","IL1B","CXCL2","CD163","MRC1","MARCO","MSR1","CD8A","CD8B","CD3D","CD4","CEBPE","HDC","MS4A2"))+coord_flip() + scale_color_viridis(option = "H")#+ theme(panel.background = element_rect(fill = 'lightgrey'))
DotPlot(MAC_Mono_sub,features=c("CD14","FCGR3A","ICAM1","CD11B","CCR2","CALR","LGALS2","LGALS1","LGALS3","LGALS4","LGALS6","LGALS7"))+coord_flip() + scale_color_viridis(option = "H")#+ theme(panel.background = element_rect(fill = 'lightgrey'))
DotPlot(MAC_Mono_sub,features=c("DCN","LUM"))+coord_flip() + scale_color_viridis(option = "H")#+ theme(panel.background = element_rect(fill = 'lightgrey'))

pdf("Graphs/Main/MAC_FB_Clusteridentification.pdf", width=10, height=5)
FeaturePlot(MAC_Mono_sub,c("CD68","DCN"), blend=T, order=T, max.cutoff = 2, cols= c("lightgrey","red","blue"))
dev.off()


#Name Cluster
Idents(MAC_Mono_sub)<-MAC_Mono_sub$seurat_clusters
MAC_Mono_sub<-RenameIdents(MAC_Mono_sub,
                           `0`="M1_M2_Mac", `1`="M2_Mac", `2`="M1_Mac", `3`="Mono",`4`="cDC",`5`="aDC",`6`="eMac",`7`="Mast",`8`="FB_Mac",`9`="M2_Mac")
MAC_Mono_sub$celltype.full<-Idents(MAC_Mono_sub)
table(MAC_Mono_sub$celltype.full)

#Doublet and contamination check
DimPlot(MAC_Mono_sub, label = TRUE, label.size = 6)
FeaturePlot(MAC_Mono_sub,c("KRT19","CD3D"), blend=T, order=T, label=T, cols= c("lightgrey","green","blue"), max.cutoff = 2)

#Define cluster levels
clusters_ordered<-c("Mono","M1_Mac","M1_M2_Mac","M2_Mac","eMac","FB_Mac","Mast" ,"cDC", "aDC")
MAC_Mono_sub$celltype.full<- factor(MAC_Mono_sub$celltype.full, levels = clusters_ordered)
Idents(MAC_Mono_sub)<-MAC_Mono_sub$celltype.full

DefaultAssay(MAC_Mono_sub)<-"RNA"
Idents(MAC_Mono_sub)<-"celltype.full"

pdf("Graphs/Main/MAC_Mono_UMAPPlots.pdf", width=5, height=5)
DimPlot(MAC_Mono_sub, label=T, label.color = "black",label.size = 5, repel=T, alpha=1) + ggtitle("MAC_Mono-subset")+ scale_color_manual(values=color_levelMAC_Mono) + NoLegend()
dev.off()

pdf("Graphs/Main/UMAPPlots_MAC_Mono_sub_QC.pdf", width=5, height=5)
DimPlot(MAC_Mono_sub, group.by="group")+ ggtitle("Basic-UMAP")
DimPlot(MAC_Mono_sub, group.by="tissue")+ ggtitle("Basic-UMAP")
dev.off()


pdf("Graphs/Main/MAC_Mono_Barplots.pdf")
plot_integrated_clusters_group (MAC_Mono_sub)
plot_integrated_clusters_Donor_ID (MAC_Mono_sub)
plot_integrated_celltypes_condition (MAC_Mono_sub)
plot_integrated_celltypes_classification (MAC_Mono_sub)
plot_integrated_celltypes_tissue (MAC_Mono_sub)
plot_integrated_celltypes_sex (MAC_Mono_sub)
plot_integrated_celltypes_age (MAC_Mono_sub) 
plot_integrated_celltypes_development (MAC_Mono_sub)
plot_integrated_celltypes_donor (MAC_Mono_sub)
plot_MAC_Mono_integrated_condition_celltypes(MAC_Mono_sub)
plot_MAC_Mono_integrated_development_celltypes(MAC_Mono_sub)
plot_MAC_Mono_integrated_group_clusters(MAC_Mono_sub)
plot_MAC_Mono_integrated_tissue_celltypes(MAC_Mono_sub)
dev.off()

#add identity to analyze tissue/health diffgenes
#Set tissue in identy 
Idents(MAC_Mono_sub) <- MAC_Mono_sub$celltype.full
MAC_Mono_sub$celltype.condition <- paste(Idents(MAC_Mono_sub), MAC_Mono_sub$condition, sep = "_")
Idents(MAC_Mono_sub)<-MAC_Mono_sub$celltype.condition
table(MAC_Mono_sub$celltype.condition)
#Set tissue in identy 
Idents(MAC_Mono_sub) <- MAC_Mono_sub$celltype.condition
MAC_Mono_sub$celltype.condition.development <- paste(Idents(MAC_Mono_sub), MAC_Mono_sub$development, sep = "_")
Idents(MAC_Mono_sub)<-MAC_Mono_sub$celltype.condition.development
table(MAC_Mono_sub$celltype.condition.development)

#Clustermarker
Idents(MAC_Mono_sub)<-"celltype.full"
MAC_Mono_sub.clustermarker<-FindAllMarkers(MAC_Mono_sub, assay = "RNA", min.pct = 0.25, logfc.threshold = 0.25)
MAC_Mono_sub.clustermarker$Foldchange_UP <- 2^(MAC_Mono_sub.clustermarker$avg_log2FC)
MAC_Mono_sub.clustermarker$Foldchange_DOWN <- 2^(-MAC_Mono_sub.clustermarker$avg_log2FC)
MAC_Mono_sub.clustermarker$Ratio_pct1_pct2 <- (MAC_Mono_sub.clustermarker$pct.1)/(MAC_Mono_sub.clustermarker$pct.2)
write.xlsx(MAC_Mono_sub.clustermarker, "Lists/Clustermarker_MAC_Mono_sub_clustermarker.xlsx")

#DotPlot Top up and down regulated genes 
#up
MAC_Mono_sub.clustermarker_UP <- subset(MAC_Mono_sub.clustermarker, Foldchange_UP >=2 )
top10_UP_clustermarker <- MAC_Mono_sub.clustermarker_UP %>% group_by(cluster) %>% top_n(n=10,wt=Foldchange_UP)
top10_UP_clustermarker<-arrange(top10_UP_clustermarker,cluster,Foldchange_UP)
top10_UP_clustermarker <- as.character(top10_UP_clustermarker$gene)
top10_UP_clustermarker<-top10_UP_clustermarker[!duplicated(top10_UP_clustermarker)]

pdf("Graphs/Main/DotPlot_Top_10UP_Clustermarker_MAC_Mono_FCmin2.pdf", width = 7, height = 15)
DotPlot(MAC_Mono_sub, features=top10_UP_clustermarker, group.by = "celltype.full",assay="RNA") + coord_flip()+ scale_color_viridis(option = "H")+theme(axis.text.x = element_text(angle = 45)) + ggtitle("MAC_Mono Cluster_DEG") + scale_y_discrete(guide = guide_axis(angle = 90))
dev.off()

#Heatmap#
## Define Gene lists that you want to plot as Heatmap
MAC_Mono_sub.clustermarker_UP<-arrange(MAC_Mono_sub.clustermarker_UP,cluster,Foldchange_UP)
genes <- as.character(MAC_Mono_sub.clustermarker_UP$gene)
avgexp <- MAC_Mono_sub

## Plot Heatmaps
## Calculate Average gene expressions (pseudobulked) for each of the celltypes in each condition
levels(avgexp)
DefaultAssay(avgexp) <- "RNA"
avgexp <- NormalizeData(avgexp)
avgexp <- AggregateExpression(avgexp, return.seurat =T, assays = "RNA")


# Modify the heatmap plot
library(RColorBrewer)
heatmap_plot <- DoHeatmap(
  avgexp,
  features = genes,
  size = 6,
  group.bar.height = 0.02,
  draw.lines = FALSE,angle = 90,
  group.colors = color_levelMAC_Mono
) +
  scale_fill_gradientn(colours = rev(brewer.pal(8, name = "RdYlBu")))

# Modify the theme to make rownames cursive
heatmap_plot 
ggsave("Graphs/Main/Heatmap_MAC_Mono_clustermarker.pdf", width = 40, height = 20, units = c("cm"))


#Enrichment analysis MAC_Mono_sub
#GO-Terms_Mono
table(MAC_Mono_sub.clustermarker_UP$cluster)
Mono_sub.clustermarker_UP <- subset(MAC_Mono_sub.clustermarker_UP, cluster=="Mono")
Mono_sub.clustermarker_UP <- as.character(Mono_sub.clustermarker_UP$gene)

if (websiteLive) {
  enriched <- enrichr(Mono_sub.clustermarker_UP, dbs)
}

if (websiteLive) enriched[["GO_Biological_Process_2023"]]
pdf("Graphs/Main/UP_Enrichment_Mono.pdf")
if (websiteLive) {
  plotEnrich(enriched[[3]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")
}
dev.off()

#GO-Terms_M1_Mac
table(MAC_Mono_sub.clustermarker_UP$cluster)
M1_Mac_sub.clustermarker_UP <- subset(MAC_Mono_sub.clustermarker_UP, cluster=="M1_Mac")
M1_Mac_sub.clustermarker_UP <- as.character(M1_Mac_sub.clustermarker_UP$gene)

if (websiteLive) {
  enriched <- enrichr(M1_Mac_sub.clustermarker_UP, dbs)
}

if (websiteLive) enriched[["GO_Biological_Process_2023"]]
pdf("Graphs/Main/UP_Enrichment_M1_Mac.pdf")
if (websiteLive) {
  plotEnrich(enriched[[3]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")
}
dev.off()

#GO-Terms_M1_M2_Mac
table(MAC_Mono_sub.clustermarker_UP$cluster)
M1_M2_Mac_sub.clustermarker_UP <- subset(MAC_Mono_sub.clustermarker_UP, cluster=="M1_M2_Mac")
M1_M2_Mac_sub.clustermarker_UP <- as.character(M1_M2_Mac_sub.clustermarker_UP$gene)

if (websiteLive) {
  enriched <- enrichr(M1_M2_Mac_sub.clustermarker_UP, dbs)
}

if (websiteLive) enriched[["GO_Biological_Process_2023"]]
pdf("Graphs/Main/UP_Enrichment_M1_M2_Mac.pdf")
if (websiteLive) {
  plotEnrich(enriched[[3]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")
}
dev.off()
#GO-Terms_M2_Mac
table(MAC_Mono_sub.clustermarker_UP$cluster)
M2_Mac_sub.clustermarker_UP <- subset(MAC_Mono_sub.clustermarker_UP, cluster=="M2_Mac")
M2_Mac_sub.clustermarker_UP <- as.character(M2_Mac_sub.clustermarker_UP$gene)

if (websiteLive) {
  enriched <- enrichr(M2_Mac_sub.clustermarker_UP, dbs)
}

if (websiteLive) enriched[["GO_Biological_Process_2023"]]
pdf("Graphs/Main/UP_Enrichment_M2_Mac.pdf")
if (websiteLive) {
  plotEnrich(enriched[[3]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")
}
dev.off()

#GO-Terms_eMac
table(MAC_Mono_sub.clustermarker_UP$cluster)
eMac_sub.clustermarker_UP <- subset(MAC_Mono_sub.clustermarker_UP, cluster=="eMac")
eMac_sub.clustermarker_UP <- as.character(eMac_sub.clustermarker_UP$gene)

if (websiteLive) {
  enriched <- enrichr(eMac_sub.clustermarker_UP, dbs)
}

if (websiteLive) enriched[["GO_Biological_Process_2023"]]
pdf("Graphs/Main/UP_Enrichment_eMac.pdf")
if (websiteLive) {
  plotEnrich(enriched[[3]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")
}
dev.off()

#GO-Terms_FB_Mac
table(MAC_Mono_sub.clustermarker_UP$cluster)
FB_Mac_sub.clustermarker_UP <- subset(MAC_Mono_sub.clustermarker_UP, cluster=="FB_Mac")
FB_Mac_sub.clustermarker_UP <- as.character(FB_Mac_sub.clustermarker_UP$gene)

if (websiteLive) {
  enriched <- enrichr(FB_Mac_sub.clustermarker_UP, dbs)
}

if (websiteLive) enriched[["GO_Biological_Process_2023"]]
pdf("Graphs/Main/UP_Enrichment_FB_Mac.pdf")
if (websiteLive) {
  plotEnrich(enriched[[3]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")
}
dev.off()

#GO-Terms_Mast
table(MAC_Mono_sub.clustermarker_UP$cluster)
Mast_sub.clustermarker_UP <- subset(MAC_Mono_sub.clustermarker_UP, cluster=="Mast")
Mast_sub.clustermarker_UP <- as.character(Mast_sub.clustermarker_UP$gene)

if (websiteLive) {
  enriched <- enrichr(Mast_sub.clustermarker_UP, dbs)
}

if (websiteLive) enriched[["GO_Biological_Process_2023"]]
pdf("Graphs/Main/UP_Enrichment_Mast.pdf")
if (websiteLive) {
  plotEnrich(enriched[[3]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")
}
dev.off()

#GO-Terms_cDC
table(MAC_Mono_sub.clustermarker_UP$cluster)
cDC_sub.clustermarker_UP <- subset(MAC_Mono_sub.clustermarker_UP, cluster=="cDC")
cDC_sub.clustermarker_UP <- as.character(cDC_sub.clustermarker_UP$gene)

if (websiteLive) {
  enriched <- enrichr(cDC_sub.clustermarker_UP, dbs)
}

if (websiteLive) enriched[["GO_Biological_Process_2023"]]
pdf("Graphs/Main/UP_Enrichment_cDC.pdf")
if (websiteLive) {
  plotEnrich(enriched[[3]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")
}
dev.off()

#GO-Terms_aDC
table(MAC_Mono_sub.clustermarker_UP$cluster)
aDC_sub.clustermarker_UP <- subset(MAC_Mono_sub.clustermarker_UP, cluster=="aDC")
aDC_sub.clustermarker_UP <- as.character(aDC_sub.clustermarker_UP$gene)

if (websiteLive) {
  enriched <- enrichr(aDC_sub.clustermarker_UP, dbs)
}

if (websiteLive) enriched[["GO_Biological_Process_2023"]]
pdf("Graphs/Main/UP_Enrichment_aDC.pdf")
if (websiteLive) {
  plotEnrich(enriched[[3]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")
}
dev.off()

save(MAC_Mono_sub,file="MAC_Mono_sub.RData")

#(27)_Subset EC####
Idents(th.total)<- th.total$celltype.full
EC_sub <- subset(th.total,idents = "EC")

DefaultAssay(EC_sub)<-"RNA"
EC_sub[["RNA"]]$data <- as(EC_sub[["RNA"]]$data, Class = "dgCMatrix")

EC_sub <- FindVariableFeatures(EC_sub)
EC_sub <- ScaleData(EC_sub)
EC_sub <- RunPCA(EC_sub)
ElbowPlot(EC_sub, ndims=30)
EC_sub <- RunUMAP(EC_sub, dims = 1:6)
EC_sub <- FindNeighbors(EC_sub, dims = 1:6)
EC_sub <- FindClusters(EC_sub, resolution = 0.1)
DimPlot(EC_sub, label = T, label.size = 3, alpha = 1)
DimPlot(EC_sub, label=T, split.by="condition", alpha = 1)
DimPlot(EC_sub, label=T, split.by="development", alpha = 1)
DimPlot(EC_sub, label=T, split.by="group", alpha = 1)
DimPlot(EC_sub, label=T, group.by="seurat_clusters", alpha = 1)
DefaultAssay(EC_sub)<- "RNA"

#Modulscore:
#PMID: 34030460
EC_aero_signature <- list(c("SPARCL1","CLDN5","HLA-E","HPGD","EDNRB","EMCN","IFI27","EPAS1","TIMP3","ADGRL2","EGFL7","RAMP2","PECAM1","AQP1","CA4",
                            "CAVIN2","IFITM3","APP","GNG11","GALNT18","SPARC","ITM2B","ICAM2"))
EC_sub <- AddModuleScore(EC_sub, features = EC_aero_signature, name = "EC_aero module score")

EC_art_signature <- list(c("EPAS1","TIMP3","LDB2","SPARCL1","CLDN5","IFITM3","IFI27","CALCRL","PTPRB","MT2A","MGP","VWF","ID1","PECAM1","RAMP2","A2M",
                           "TM4SF1","CAV1","TFPI"))
EC_sub <- AddModuleScore(EC_sub, features = EC_art_signature, name = "EC_art module score")

EC_cap_signature <- list(c("FCN3","EPAS1","TIMP3","CLDN5","LDB2","SPARCL1","IFITM3","IFI27","VWF","PECAM1","TMEM100","MT2A","RAMP2","HLA-E","EGFL7",
                           "PTPRB","BTNL9","SLCO2A1","CALCRL"))
EC_sub <- AddModuleScore(EC_sub, features = EC_cap_signature, name = "EC_cap module score")

EC_lym_signature <- list(c("CCL21","IGFBP7","TFF3","GNG11","MMRN1","TFPI","PPFIBP1","AKAP12","TIMP3","ADIRF","CLDN5"))
EC_sub <- AddModuleScore(EC_sub, features = EC_lym_signature, name = "EC_lym module score")

EC_ven_signature <- list(c("VWF","TIMP3","CLU","SPARCL1","IGFBP7","ACKR1","CALCRL","EPAS1","PTPRB","IFITM3","SLCO2A1","PRSS23","PECAM1","MGP",
                           "GPX3","C7","RAMP3","LIFR"))
EC_sub <- AddModuleScore(EC_sub, features = EC_ven_signature, name = "EC_ven module score")

EC_systven_signature <- list(c("SPARCL1","ACKR1","IGFBP7","LDB2","ZNF385D","EMCN","IFI27"))
EC_sub <- AddModuleScore(EC_sub, features = EC_systven_signature, name = "EC_systven module score")

cTEC_lo_Bautista_signature <- list(c("TCF24","COL26A1","GJB6","CCL25","SCUBE1","FBN3","GALNT9","COL9A1","ENO3","UBE2C","CENPA","GLB1L3","KIF4A","CDCA3","DEPDC1",
                                     "LINC00618"))
EC_sub <- AddModuleScore(EC_sub, features = cTEC_lo_Bautista_signature, name = "cTEC_lo_Bautista module score")
cTEC_hi_Bautista_signature <- list(c("CFC1","TNFRSF17","CFHR1","PRSS16","RTBDN","FOXR1","PGAM2","TBATA","CLEC2L","CCL25",
                                     "PNCK","LHCGR","C4orf50","NPHS1"))
EC_sub <- AddModuleScore(EC_sub, features = cTEC_hi_Bautista_signature, name = "cTEC_hi_Bautista module score")
mTEC_lo_Bautista_signature <- list(c("LINC00922","CXCL9","LYPD1","LRFN2","CCL19","MROH2A","KRT15","FOXQ1","CXCL10","MMP9","GPR88","CXCL11","LINC00839","FXYD2",
                                     "IL22RA2","HTR1F","CHI3L1"))
EC_sub <- AddModuleScore(EC_sub, features = mTEC_lo_Bautista_signature, name = "mTEC_lo_Bautista module score")
mTEC_hi_Bautista_signature <- list(c("EDDM3B","BANF2","UMOD","AMTN","SPATA16","CSN1S1","CLDN18","DNMT3L","S100A5","CLEC4G","GIP","CD209","NT5DC4","SERPINA6",
                                     "HPD","COLEC10"))
EC_sub <- AddModuleScore(EC_sub, features = mTEC_hi_Bautista_signature, name = "mTEC_hi_Bautista module score")


module_EC_subset_signature <- c("EC_art module score1","EC_cap module score1","EC_lym module score1","EC_ven module score1","EC_systven module score1")
module_EC_subset_bautista_signature <- c("EC1_bautista module score1","EC2_bautista module score1","EC3_bautista module score1","EC4_bautista module score1")


#Modulescore imaging
pdf("Graphs/Main/Modulescore_EC_Clusteridentification.pdf")
DotPlot(object = EC_sub, features = module_EC_subset_bautista_signature,assay="RNA") + coord_flip()+ scale_color_viridis(option = "H")+theme(axis.text.x = element_text(angle = 45)) + ggtitle("EC_modulescore") + scale_y_discrete(guide = guide_axis(angle = 90))
DotPlot(object = EC_sub, features = module_EC_subset_signature,assay="RNA") + coord_flip()+ scale_color_viridis(option = "H")+theme(axis.text.x = element_text(angle = 45)) + ggtitle("EC_modulescore PMID:34030460") + scale_y_discrete(guide = guide_axis(angle = 90))
FeaturePlot(EC_sub,"LYVE1",order=T)
FeaturePlot(object = EC_sub, features = "EC_Xin module score1")
VlnPlot(object = EC_sub, features = "EC_Xin module score1", pt.size=0)
DotPlot(object = EC_sub, features = yayon_EC,assay="RNA") + coord_flip()+ scale_color_viridis(option = "H")+theme(axis.text.x = element_text(angle = 45)) + ggtitle("EC_modulescore PMID:34030460") + scale_y_discrete(guide = guide_axis(angle = 90))
DotPlot(object = EC_sub, features = c("FCN2","HCK","DHRS9","MMP11","HMCN1","CTHRC1","ZBED2","CALN1","PSMB11","PRAP1","COL2A1","OLFM3","CCL21","MT1A","CFD","KRTDAP","SBSN","PSORS1C2",
                                      "DPYS","CH25H","KRT14","GPC3","NR2F1","OGN","PTPRCAP","XCL1","CD27","MYL4","MYOD1","MYH3","TRIM63","ASB5","CSRP3","ATOH1","CCER2","STMN2","PCSK1N",                                      
                                      "TMEM61","CLDN9","ERBB4","CCDC144NL-AS1","CASC9","ACKR1","SOX17","RAMP3","CCL14"),assay="RNA") + coord_flip()+ scale_color_viridis(option = "H")+theme(axis.text.x = element_text(angle = 45)) + ggtitle("EC_modulescore PMID:34030460") + scale_y_discrete(guide = guide_axis(angle = 90))
dev.off()

#Name Cluster
Idents(EC_sub)<-EC_sub$seurat_clusters
EC_sub<-RenameIdents(EC_sub,
                     `0`="EC_Cap_1", `1`="EC_Cap_2", `2`="EC_Ven_1", `3`="EC_ART",`4`="DP",`5`="EC_Ven_2",
                     `6`="DP")
EC_sub$celltype.full<-Idents(EC_sub)
table(EC_sub$celltype.full)

#Doublet and contamination check
DimPlot(EC_sub, label = TRUE, label.size = 6)
FeaturePlot(EC_sub,c("KRT19","CD3D"), blend=T, order=T, label=T, cols= c("lightgrey","green","blue"), max.cutoff = 2)

cells.of.interest_3 <- WhichCells(
  EC_sub,
  idents = "DP"
)

EC_sub <- subset(EC_sub,idents = c("EC_Cap_1","EC_Cap_2","EC_Ven_1","EC_Ven_2","EC_ART"))


#Define cluster levels
clusters_ordered<-c("EC_Cap_1","EC_Cap_2","EC_Ven_1","EC_Ven_2","EC_ART")
EC_sub$celltype.full<- factor(EC_sub$celltype.full, levels = clusters_ordered)
Idents(EC_sub)<-EC_sub$celltype.full

DefaultAssay(EC_sub)<-"RNA"
Idents(EC_sub)<-"celltype.full"

pdf("Graphs/Main/EC_UMAPPlots.pdf", width=5, height=5)
DimPlot(EC_sub, label=T, label.color = "black",label.size = 5, repel=T, alpha=1) + ggtitle("EC-subset")+ scale_color_manual(values=color_levelEC) + NoLegend()
dev.off()

pdf("Graphs/Main/UMAPPlots_EC_sub_QC.pdf", width=5, height=5)
DimPlot(EC_sub, group.by="group")+ ggtitle("Basic-UMAP")
DimPlot(EC_sub, group.by="tissue")+ ggtitle("Basic-UMAP")
dev.off()

pdf("Graphs/Main/EC_Barplots.pdf")
plot_integrated_clusters_group (EC_sub)
plot_integrated_clusters_Donor_ID (EC_sub)
plot_integrated_celltypes_condition (EC_sub)
plot_integrated_celltypes_classification (EC_sub)
plot_integrated_celltypes_tissue (EC_sub)
plot_integrated_celltypes_sex (EC_sub)
plot_integrated_celltypes_age (EC_sub) 
plot_integrated_celltypes_development (EC_sub)
plot_integrated_celltypes_donor (EC_sub)
plot_EC_integrated_condition_celltypes(EC_sub)
plot_EC_integrated_development_celltypes(EC_sub)
plot_EC_integrated_group_clusters(EC_sub)
plot_EC_integrated_tissue_celltypes(EC_sub)
dev.off()

#add identity to analyze tissue/health diffgenes
#Set tissue in identy 
Idents(EC_sub) <- EC_sub$celltype.full
EC_sub$celltype.condition <- paste(Idents(EC_sub), EC_sub$condition, sep = "_")
Idents(EC_sub)<-EC_sub$celltype.condition
table(EC_sub$celltype.condition)
#Set tissue in identy 
Idents(EC_sub) <- EC_sub$celltype.condition
EC_sub$celltype.condition.development <- paste(Idents(EC_sub), EC_sub$development, sep = "_")
Idents(EC_sub)<-EC_sub$celltype.condition.development
table(EC_sub$celltype.condition.development)

#Clustermarker
Idents(EC_sub)<-"celltype.full"
EC_sub.clustermarker<-FindAllMarkers(EC_sub, assay = "RNA", min.pct = 0.25, logfc.threshold = 0.25)
EC_sub.clustermarker$Foldchange_UP <- 2^(EC_sub.clustermarker$avg_log2FC)
EC_sub.clustermarker$Foldchange_DOWN <- 2^(-EC_sub.clustermarker$avg_log2FC)
EC_sub.clustermarker$Ratio_pct1_pct2 <- (EC_sub.clustermarker$pct.1)/(EC_sub.clustermarker$pct.2)
write.xlsx(EC_sub.clustermarker, "Lists/Clustermarker_EC_sub_clustermarker.xlsx")

#DotPlot Top up and down regulated genes 
#up
EC_sub.clustermarker_UP <- subset(EC_sub.clustermarker, Foldchange_UP >=2 )
top10_UP_clustermarker <- EC_sub.clustermarker_UP %>% group_by(cluster) %>% top_n(n=10,wt=Foldchange_UP)
top10_UP_clustermarker<-arrange(top10_UP_clustermarker,cluster,Foldchange_UP)
top10_UP_clustermarker <- as.character(top10_UP_clustermarker$gene)
top10_UP_clustermarker<-top10_UP_clustermarker[!duplicated(top10_UP_clustermarker)]

pdf("Graphs/Main/DotPlot_Top_10UP_Clustermarker_EC_FCmin2.pdf", width = 7, height = 12)
DotPlot(EC_sub, features=top10_UP_clustermarker, group.by = "celltype.full",assay="RNA") + coord_flip()+ scale_color_viridis(option = "H")+theme(axis.text.x = element_text(angle = 45)) + ggtitle("EC Cluster_DEG") + scale_y_discrete(guide = guide_axis(angle = 90))
dev.off()

#Heatmap#
## Define Gene lists that you want to plot as Heatmap
EC_sub.clustermarker_UP<-arrange(EC_sub.clustermarker_UP,cluster,Foldchange_UP)
genes <- as.character(EC_sub.clustermarker_UP$gene)
avgexp <- EC_sub

## Plot Heatmaps
## Calculate Average gene expressions (pseudobulked) for each of the celltypes in each condition
levels(avgexp)
DefaultAssay(avgexp) <- "RNA"
avgexp <- NormalizeData(avgexp)
avgexp <- AggregateExpression(avgexp, return.seurat =T, assays = "RNA")


# Modify the heatmap plot
library(RColorBrewer)
heatmap_plot <- DoHeatmap(
  avgexp,
  features = genes,
  size = 6,
  group.bar.height = 0.02,
  draw.lines = FALSE,angle = 90,
  group.colors = color_levelEC
) +
  scale_fill_gradientn(colours = rev(brewer.pal(8, name = "RdYlBu")))

# Modify the theme to make rownames cursive
heatmap_plot 
ggsave("Graphs/Main/Heatmap_EC_clustermarker.pdf", width = 40, height = 20, units = c("cm"))

#Enrichment analysis EC_sub
#GO-Terms_EC_Cap_1
table(EC_sub.clustermarker_UP$cluster)
EC_Cap_1_sub.clustermarker_UP <- subset(EC_sub.clustermarker_UP, cluster=="EC_Cap_1")
EC_Cap_1_sub.clustermarker_UP <- as.character(EC_Cap_1_sub.clustermarker_UP$gene)

if (websiteLive) {
  enriched <- enrichr(EC_Cap_1_sub.clustermarker_UP, dbs)
}

if (websiteLive) enriched[["GO_Biological_Process_2023"]]
pdf("Graphs/Main/UP_Enrichment_EC_Cap_1.pdf")
if (websiteLive) {
  plotEnrich(enriched[[3]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")
}
dev.off()

#GO-Terms_EC_Cap_2
table(EC_sub.clustermarker_UP$cluster)
EC_Cap_2_sub.clustermarker_UP <- subset(EC_sub.clustermarker_UP, cluster=="EC_Cap_2")
EC_Cap_2_sub.clustermarker_UP <- as.character(EC_Cap_2_sub.clustermarker_UP$gene)

if (websiteLive) {
  enriched <- enrichr(EC_Cap_2_sub.clustermarker_UP, dbs)
}

if (websiteLive) enriched[["GO_Biological_Process_2023"]]
pdf("Graphs/Main/UP_Enrichment_EC_Cap_2.pdf")
if (websiteLive) {
  plotEnrich(enriched[[3]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")
}
dev.off()

#GO-Terms_EC_Ven_1
table(EC_sub.clustermarker_UP$cluster)
EC_Ven_1_sub.clustermarker_UP <- subset(EC_sub.clustermarker_UP, cluster=="EC_Ven_1")
EC_Ven_1_sub.clustermarker_UP <- as.character(EC_Ven_1_sub.clustermarker_UP$gene)

if (websiteLive) {
  enriched <- enrichr(EC_Ven_1_sub.clustermarker_UP, dbs)
}

if (websiteLive) enriched[["GO_Biological_Process_2023"]]
pdf("Graphs/Main/UP_Enrichment_EC_Ven_1.pdf")
if (websiteLive) {
  plotEnrich(enriched[[3]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")
}
dev.off()

#GO-Terms_EC_Ven_2
table(EC_sub.clustermarker_UP$cluster)
EC_Ven_2_sub.clustermarker_UP <- subset(EC_sub.clustermarker_UP, cluster=="EC_Ven_2")
EC_Ven_2_sub.clustermarker_UP <- as.character(EC_Ven_2_sub.clustermarker_UP$gene)

if (websiteLive) {
  enriched <- enrichr(EC_Ven_2_sub.clustermarker_UP, dbs)
}

if (websiteLive) enriched[["GO_Biological_Process_2023"]]
pdf("Graphs/Main/UP_Enrichment_EC_Ven_2.pdf")
if (websiteLive) {
  plotEnrich(enriched[[3]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")
}
dev.off()

#GO-Terms_EC_ART
table(EC_sub.clustermarker_UP$cluster)
EC_ART_sub.clustermarker_UP <- subset(EC_sub.clustermarker_UP, cluster=="EC_ART")
EC_ART_sub.clustermarker_UP <- as.character(EC_ART_sub.clustermarker_UP$gene)

if (websiteLive) {
  enriched <- enrichr(EC_ART_sub.clustermarker_UP, dbs)
}

if (websiteLive) enriched[["GO_Biological_Process_2023"]]
pdf("Graphs/Main/UP_Enrichment_EC_ART.pdf")
if (websiteLive) {
  plotEnrich(enriched[[3]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")
}
dev.off()

save(EC_sub,file="EC_sub.RData")

#(28)_Subset VSMC/Peri####
Idents(th.total)<- th.total$celltype.full
VSMC_Peri_sub <- subset(th.total,idents = "VSMC/Peri")

DefaultAssay(VSMC_Peri_sub)<-"RNA"
VSMC_Peri_sub[["RNA"]]$data <- as(VSMC_Peri_sub[["RNA"]]$data, Class = "dgCMatrix")

VSMC_Peri_sub <- FindVariableFeatures(VSMC_Peri_sub)
VSMC_Peri_sub <- ScaleData(VSMC_Peri_sub)
VSMC_Peri_sub <- RunPCA(VSMC_Peri_sub)
ElbowPlot(VSMC_Peri_sub, ndims=20)
VSMC_Peri_sub <- RunUMAP(VSMC_Peri_sub, dims = 1:5)
VSMC_Peri_sub <- FindNeighbors(VSMC_Peri_sub, dims = 1:5)
VSMC_Peri_sub <- FindClusters(VSMC_Peri_sub, resolution = 0.1)
DimPlot(VSMC_Peri_sub, label = T, label.size = 3, alpha = 1)
DimPlot(VSMC_Peri_sub, label=T, split.by="condition", alpha = 1)
DimPlot(VSMC_Peri_sub, label=T, split.by="development", alpha = 1)
DimPlot(VSMC_Peri_sub, label=T, split.by="group", alpha = 1)
DimPlot(VSMC_Peri_sub, label=T, group.by="seurat_clusters", alpha = 1)
DefaultAssay(VSMC_Peri_sub)<- "RNA"

#Modulscore:
#PMID: 37189183

#Modulescore imaging
pdf("Graphs/Main/Modulescore_VSMC_Peri_Clusteridentification.pdf")
DotPlot(object = VSMC_Peri_sub, features = c("ACTA2","MYH11","COL1A1","COL1A2","PDGFRA","FABP4","EBF2","CD14","CD68","CD3D","CD3G","ENG"),assay="RNA") + coord_flip()+ scale_color_viridis(option = "H")+theme(axis.text.x = element_text(angle = 45)) + ggtitle("VSMC_Peri_modulescore") + scale_y_discrete(guide = guide_axis(angle = 90))
DotPlot(object = VSMC_Peri_sub, features = yayon_VSMC,assay="RNA") + coord_flip()+ scale_color_viridis(option = "H")+theme(axis.text.x = element_text(angle = 45)) + ggtitle("VSMC_Peri_modulescore") + scale_y_discrete(guide = guide_axis(angle = 90))
FeaturePlot(object = VSMC_Peri_sub, features = "Peri_bautista module score1", order=T)
FeaturePlot(object = VSMC_Peri_sub, features = "VSMC_Xin module score1", order=T)
FeaturePlot(VSMC_Peri_sub,"ACTA2")
dev.off()

#Name Cluster
Idents(VSMC_Peri_sub)<-VSMC_Peri_sub$seurat_clusters
VSMC_Peri_sub<-RenameIdents(VSMC_Peri_sub,
                            `0`="Pericyte_COL1A1", `1`="SMC1", `2`="SMC2", `3`="Pericyte_CCL19",`4`="FetalPericyte",`5`="SMC3",`6`="ProlifPericyte")
VSMC_Peri_sub$celltype.full<-Idents(VSMC_Peri_sub)
table(VSMC_Peri_sub$celltype.full)

#Doublet and contamination check
DimPlot(VSMC_Peri_sub, label = TRUE, label.size = 6)
FeaturePlot(VSMC_Peri_sub,c("KRT19","CD3D"), blend=T, order=T, label=T, cols= c("lightgrey","green","blue"), max.cutoff = 2)

#Define cluster levels
clusters_ordered<-c("Pericyte_CCL19","Pericyte_COL1A1","ProlifPericyte","FetalPericyte","SMC1","SMC2","SMC3")
VSMC_Peri_sub$celltype.full<- factor(VSMC_Peri_sub$celltype.full, levels = clusters_ordered)
Idents(VSMC_Peri_sub)<-VSMC_Peri_sub$celltype.full

DefaultAssay(VSMC_Peri_sub)<-"RNA"
Idents(VSMC_Peri_sub)<-"celltype.full"

pdf("Graphs/Main/VSMC_Peri_UMAPPlots.pdf", width=5, height=5)
DimPlot(VSMC_Peri_sub, label=T, label.color = "black",label.size = 5, repel=T, alpha=1) + ggtitle("VSMC_Peri-subset")+ scale_color_manual(values=color_levelVSMC_Peri) + NoLegend()
dev.off()

pdf("Graphs/Main/UMAPPlots_VSMC_Peri_sub_QC.pdf", width=5, height=5)
DimPlot(VSMC_Peri_sub, group.by="group")+ ggtitle("Basic-UMAP")
DimPlot(VSMC_Peri_sub, group.by="tissue")+ ggtitle("Basic-UMAP")
dev.off()


pdf("Graphs/Main/VSMC_Peri_Barplots.pdf")
plot_integrated_clusters_group (VSMC_Peri_sub)
plot_integrated_clusters_Donor_ID (VSMC_Peri_sub)
plot_integrated_celltypes_condition (VSMC_Peri_sub)
plot_integrated_celltypes_classification (VSMC_Peri_sub)
plot_integrated_celltypes_tissue (VSMC_Peri_sub)
plot_integrated_celltypes_sex (VSMC_Peri_sub)
plot_integrated_celltypes_age (VSMC_Peri_sub) 
plot_integrated_celltypes_development (VSMC_Peri_sub)
plot_integrated_celltypes_donor (VSMC_Peri_sub)
plot_VSMC_Peri_integrated_condition_celltypes(VSMC_Peri_sub)
plot_VSMC_Peri_integrated_development_celltypes(VSMC_Peri_sub)
plot_VSMC_Peri_integrated_group_clusters(VSMC_Peri_sub)
plot_VSMC_Peri_integrated_tissue_celltypes(VSMC_Peri_sub)
dev.off()

#add identity to analyze tissue/health diffgenes
#Set tissue in identy 
Idents(VSMC_Peri_sub) <- VSMC_Peri_sub$celltype.full
VSMC_Peri_sub$celltype.condition <- paste(Idents(VSMC_Peri_sub), VSMC_Peri_sub$condition, sep = "_")
Idents(VSMC_Peri_sub)<-VSMC_Peri_sub$celltype.condition
table(VSMC_Peri_sub$celltype.condition)
#Set tissue in identy 
Idents(VSMC_Peri_sub) <- VSMC_Peri_sub$celltype.condition
VSMC_Peri_sub$celltype.condition.development <- paste(Idents(VSMC_Peri_sub), VSMC_Peri_sub$development, sep = "_")
Idents(VSMC_Peri_sub)<-VSMC_Peri_sub$celltype.condition.development
table(VSMC_Peri_sub$celltype.condition.development)

#Clustermarker
Idents(VSMC_Peri_sub)<-"celltype.full"
VSMC_Peri_sub.clustermarker<-FindAllMarkers(VSMC_Peri_sub, assay = "RNA", min.pct = 0.25, logfc.threshold = 0.25)
VSMC_Peri_sub.clustermarker$Foldchange_UP <- 2^(VSMC_Peri_sub.clustermarker$avg_log2FC)
VSMC_Peri_sub.clustermarker$Foldchange_DOWN <- 2^(-VSMC_Peri_sub.clustermarker$avg_log2FC)
VSMC_Peri_sub.clustermarker$Ratio_pct1_pct2 <- (VSMC_Peri_sub.clustermarker$pct.1)/(VSMC_Peri_sub.clustermarker$pct.2)
write.xlsx(VSMC_Peri_sub.clustermarker, "Lists/Clustermarker_VSMC_Peri_sub_clustermarker.xlsx")

#DotPlot Top up and down regulated genes 
#up
VSMC_Peri_sub.clustermarker_UP <- subset(VSMC_Peri_sub.clustermarker, Foldchange_UP >=2 )
top10_UP_clustermarker <- VSMC_Peri_sub.clustermarker_UP %>% group_by(cluster) %>% top_n(n=10,wt=Foldchange_UP)
top10_UP_clustermarker<-arrange(top10_UP_clustermarker,cluster,Foldchange_UP)
top10_UP_clustermarker <- as.character(top10_UP_clustermarker$gene)
top10_UP_clustermarker<-top10_UP_clustermarker[!duplicated(top10_UP_clustermarker)]

pdf("Graphs/Main/DotPlot_Top_10UP_Clustermarker_VSMC_Peri_FCmin2.pdf", width = 7, height = 12)
DotPlot(VSMC_Peri_sub, features=top10_UP_clustermarker, group.by = "celltype.full",assay="RNA") + coord_flip()+ scale_color_viridis(option = "H")+theme(axis.text.x = element_text(angle = 45)) + ggtitle("VSMC_Peri Cluster_DEG") + scale_y_discrete(guide = guide_axis(angle = 90))
dev.off()

#Heatmap#
## Define Gene lists that you want to plot as Heatmap
VSMC_Peri_sub.clustermarker_UP<-arrange(VSMC_Peri_sub.clustermarker_UP,cluster,Foldchange_UP)
genes <- as.character(VSMC_Peri_sub.clustermarker_UP$gene)
avgexp <- VSMC_Peri_sub

## Plot Heatmaps
## Calculate Average gene expressions (pseudobulked) for each of the celltypes in each condition
levels(avgexp)
DefaultAssay(avgexp) <- "RNA"
avgexp <- NormalizeData(avgexp)
avgexp <- AggregateExpression(avgexp, return.seurat =T, assays = "RNA")


# Modify the heatmap plot
library(RColorBrewer)
heatmap_plot <- DoHeatmap(
  avgexp,
  features = genes,
  size = 6,
  group.bar.height = 0.02,
  draw.lines = FALSE,angle = 90,
  group.colors = color_levelVSMC_Peri
) +
  scale_fill_gradientn(colours = rev(brewer.pal(8, name = "RdYlBu")))

# Modify the theme to make rownames cursive
heatmap_plot 
ggsave("Graphs/Main/Heatmap_VSMC_Peri_clustermarker.pdf", width = 40, height = 20, units = c("cm"))

#Enrichment analysis VSMC_Peri_sub
#GO-Terms_Pericyte_CCL19
table(VSMC_Peri_sub.clustermarker_UP$cluster)
Pericyte_CCL19_sub.clustermarker_UP <- subset(VSMC_Peri_sub.clustermarker_UP, cluster=="Pericyte_CCL19")
Pericyte_CCL19_sub.clustermarker_UP <- as.character(Pericyte_CCL19_sub.clustermarker_UP$gene)

if (websiteLive) {
  enriched <- enrichr(Pericyte_CCL19_sub.clustermarker_UP, dbs)
}

if (websiteLive) enriched[["GO_Biological_Process_2023"]]
pdf("Graphs/Main/UP_Enrichment_Pericyte_CCL19.pdf")
if (websiteLive) {
  plotEnrich(enriched[[3]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")
}
dev.off()

#GO-Terms_Pericyte_COL1A1
table(VSMC_Peri_sub.clustermarker_UP$cluster)
Pericyte_COL1A1_sub.clustermarker_UP <- subset(VSMC_Peri_sub.clustermarker_UP, cluster=="Pericyte_COL1A1")
Pericyte_COL1A1_sub.clustermarker_UP <- as.character(Pericyte_COL1A1_sub.clustermarker_UP$gene)

if (websiteLive) {
  enriched <- enrichr(Pericyte_COL1A1_sub.clustermarker_UP, dbs)
}

if (websiteLive) enriched[["GO_Biological_Process_2023"]]
pdf("Graphs/Main/UP_Enrichment_Pericyte_COL1A1.pdf")
if (websiteLive) {
  plotEnrich(enriched[[3]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")
}
dev.off()

#GO-Terms_ProlifPericyte
table(VSMC_Peri_sub.clustermarker_UP$cluster)
ProlifPericyte_sub.clustermarker_UP <- subset(VSMC_Peri_sub.clustermarker_UP, cluster=="ProlifPericyte")
ProlifPericyte_sub.clustermarker_UP <- as.character(ProlifPericyte_sub.clustermarker_UP$gene)

if (websiteLive) {
  enriched <- enrichr(ProlifPericyte_sub.clustermarker_UP, dbs)
}

if (websiteLive) enriched[["GO_Biological_Process_2023"]]
pdf("Graphs/Main/UP_Enrichment_ProlifPericyte.pdf")
if (websiteLive) {
  plotEnrich(enriched[[3]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")
}
dev.off()

#GO-Terms_FetalPericyte
table(VSMC_Peri_sub.clustermarker_UP$cluster)
FetalPericyte_sub.clustermarker_UP <- subset(VSMC_Peri_sub.clustermarker_UP, cluster=="FetalPericyte")
FetalPericyte_sub.clustermarker_UP <- as.character(FetalPericyte_sub.clustermarker_UP$gene)

if (websiteLive) {
  enriched <- enrichr(FetalPericyte_sub.clustermarker_UP, dbs)
}

if (websiteLive) enriched[["GO_Biological_Process_2023"]]
pdf("Graphs/Main/UP_Enrichment_FetalPericyte.pdf")
if (websiteLive) {
  plotEnrich(enriched[[3]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")
}
dev.off()

#GO-Terms_SMC1
table(VSMC_Peri_sub.clustermarker_UP$cluster)
SMC1_sub.clustermarker_UP <- subset(VSMC_Peri_sub.clustermarker_UP, cluster=="SMC1")
SMC1_sub.clustermarker_UP <- as.character(SMC1_sub.clustermarker_UP$gene)

if (websiteLive) {
  enriched <- enrichr(SMC1_sub.clustermarker_UP, dbs)
}

if (websiteLive) enriched[["GO_Biological_Process_2023"]]
pdf("Graphs/Main/UP_Enrichment_SMC1.pdf")
if (websiteLive) {
  plotEnrich(enriched[[3]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")
}
dev.off()

#GO-Terms_SMC2
table(VSMC_Peri_sub.clustermarker_UP$cluster)
SMC2_sub.clustermarker_UP <- subset(VSMC_Peri_sub.clustermarker_UP, cluster=="SMC2")
SMC2_sub.clustermarker_UP <- as.character(SMC2_sub.clustermarker_UP$gene)

if (websiteLive) {
  enriched <- enrichr(SMC2_sub.clustermarker_UP, dbs)
}

if (websiteLive) enriched[["GO_Biological_Process_2023"]]
pdf("Graphs/Main/UP_Enrichment_SMC2.pdf")
if (websiteLive) {
  plotEnrich(enriched[[3]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")
}
dev.off()

#GO-Terms_SMC3
table(VSMC_Peri_sub.clustermarker_UP$cluster)
SMC3_sub.clustermarker_UP <- subset(VSMC_Peri_sub.clustermarker_UP, cluster=="SMC3")
SMC3_sub.clustermarker_UP <- as.character(SMC3_sub.clustermarker_UP$gene)

if (websiteLive) {
  enriched <- enrichr(SMC3_sub.clustermarker_UP, dbs)
}

if (websiteLive) enriched[["GO_Biological_Process_2023"]]
pdf("Graphs/Main/UP_Enrichment_SMC3.pdf")
if (websiteLive) {
  plotEnrich(enriched[[3]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")
}
dev.off()


save(VSMC_Peri_sub,file="VSMC_Peri_sub.RData")




#(29)_Doublet & Contamination Removal####
FeaturePlot(th.total,c("KRT19","CD3D"), blend=T, order=T, label=T, cols= c("lightgrey","green","blue"), max.cutoff = 2)
FeaturePlot(th.total,c("KRT19","CD3E"), blend=T, order=T, label=T, cols= c("lightgrey","green","blue"), max.cutoff = 2)
FeaturePlot(th.total,c("KRT19","CD3G"), blend=T, order=T, label=T, cols= c("lightgrey","green","blue"), max.cutoff = 2)

# DP1
th.total$cont1 <- ifelse(
  Cells(th.total) %in% cells.of.interest_1,
  "DP",
  "other"
)


DimPlot(
  th.total,
  reduction = "umap.full",
  group.by = "cont1",
  cols = c("grey80", "red")
)

Idents(th.total)<- th.total$cont1
th.total <- subset(th.total,idents = "other")
DimPlot(th.total, reduction = "umap.full", group.by = "celltype.full", alpha = 0.4)+ ggtitle("Basic-UMAP")+ scale_color_manual(values=color_level1)

# DP2
th.total$cont2 <- ifelse(
  Cells(th.total) %in% cells.of.interest_2,
  "DP",
  "other"
)

DimPlot(
  th.total,
  reduction = "umap.full",
  group.by = "cont2",
  cols = c( "red","grey80")
)

Idents(th.total)<- th.total$cont2
th.total <- subset(th.total,idents = "other")
DimPlot(th.total, reduction = "umap.full", group.by = "celltype.full", alpha = 0.4)+ ggtitle("Basic-UMAP")+ scale_color_manual(values=color_level1)

# DP3
th.total$cont3 <- ifelse(
  Cells(th.total) %in% cells.of.interest_3,
  "DP",
  "other"
)

DimPlot(
  th.total,
  reduction = "umap.full",
  group.by = "cont3",
  cols = c("grey80", "red")
)

Idents(th.total)<- th.total$cont3
th.total <- subset(th.total,idents = "other")
DimPlot(th.total, reduction = "umap.full", group.by = "celltype.full", alpha = 0.4)+ ggtitle("Basic-UMAP")+ scale_color_manual(values=color_level1)

FeaturePlot(th.total,c("KRT19","CD3D"), blend=T, order=T, label=T, cols= c("lightgrey","green","blue"), max.cutoff = 2)
FeaturePlot(th.total,c("KRT19","CD3E"), blend=T, order=T, label=T, cols= c("lightgrey","green","blue"), max.cutoff = 2)
FeaturePlot(th.total,c("KRT19","CD3G"), blend=T, order=T, label=T, cols= c("lightgrey","green","blue"), max.cutoff = 2)

#Remove remaining RBCs
Idents(th.total)<-"celltype.full"
th.total <- subset(th.total,idents = c("TC","BC","PC","TEC","DC", "MAC/Mono","EC","FB","VSMC/Peri"))




table(th.total$celltype)
th.total$celltype<-th.total$celltype.full
table(th.total$celltype)
table(th.total$celltype.full)

#

#(30)_UMAPPlots, named cellgroups####

clusters_ordered<-c("TC","BC","PC","TEC","DC", "MAC/Mono","EC","FB","VSMC/Peri")
th.total$celltype.full<- factor(th.total$celltype.full, levels = clusters_ordered)
clusters_ordered<-c("TC","BC","PC","TEC","DC", "MAC/Mono","EC","FB","VSMC/Peri")
th.total$celltype<- factor(th.total$celltype, levels = clusters_ordered)


pdf("Graphs/Main/UMAPPlots_th.total.pdf", width=5, height=5)
DimPlot(th.total, reduction = "umap.full", group.by = "celltype.full", alpha = 0.4)+ ggtitle("Basic-UMAP")+ scale_color_manual(values=color_level1)
dev.off()

pdf("Graphs/Main/UMAPPlots_th.total_QC.pdf", width=5, height=5)
DimPlot(th.total, reduction = "umap.full",  group.by="group", split.by="tissue", ncol=3)+ ggtitle("Basic-UMAP")
DimPlot(th.total, reduction = "umap.full",  group.by="tissue", alpha = 0.1)+ ggtitle("Basic-UMAP")
dev.off()


pdf("Graphs/Main/UMAPPlots_th.total.identified.pdf")
DimPlot(th.total, reduction = "umap.full",  split.by="group", ncol = 3, alpha = 0.1)+ ggtitle("Basic-UMAP")+ scale_color_manual(values=color_level1)
DimPlot(th.total, reduction = "umap.full",  split.by="condition", ncol = 3, alpha = 0.1)+ ggtitle("Basic-UMAP")+ scale_color_manual(values=color_level1)
DimPlot(th.total, reduction = "umap.full",  split.by="classification", ncol = 3, alpha = 0.1)+ ggtitle("Basic-UMAP")+ scale_color_manual(values=color_level1)
DimPlot(th.total, reduction = "umap.full",  split.by="sex", ncol = 3, alpha = 0.1)+ ggtitle("Basic-UMAP")+ scale_color_manual(values=color_level1)
DimPlot(th.total, reduction = "umap.full",  split.by="age", ncol = 3, alpha = 0.1)+ ggtitle("Basic-UMAP")+ scale_color_manual(values=color_level1)
DimPlot(th.total, reduction = "umap.full",  split.by="development", ncol = 3, alpha = 0.1)+ ggtitle("Basic-UMAP")+ scale_color_manual(values=color_level1)
DimPlot(th.total, reduction = "umap.full",  split.by="maintissue", ncol = 3, alpha = 0.4)+ ggtitle("Basic-UMAP")+ scale_color_manual(values=color_level1)
dev.off()

pdf("Graphs/Main/UMAPPlots_th.total.identified_maintissue.pdf", height=6, width=15)
DimPlot(th.total, reduction = "umap.full",  split.by="maintissue", group.by = "celltype.full", ncol = 3, alpha = 0.1, raster=F)+ ggtitle("Basic-UMAP")+ scale_color_manual(values=color_level1)
dev.off()

#Modulescore imaging
pdf("Graphs/Main/Modulescore_Clusteridentification.pdf", height=20, width=15)
DotPlot(object = th.total, features = module_th.total_signature,assay="RNA") + coord_flip()+ scale_color_viridis(option = "H")+theme(axis.text.x = element_text(angle = 45)) + ggtitle("Main_modulescore") + scale_y_discrete(guide = guide_axis(angle = 90))
dev.off()

#Clustermarker
th.total.clustermarker<-FindAllMarkers(th.total,verbose = T, min.pct = 0.25, logfc.threshold = 0.25)
head(th.total.clustermarker)
th.total.clustermarker$Foldchange_UP <- 2^(th.total.clustermarker$avg_log2FC)
th.total.clustermarker$Foldchange_DOWN <- 2^(-th.total.clustermarker$avg_log2FC)
th.total.clustermarker$Ratio_pct1_pct2 <- (th.total.clustermarker$pct.1)/(th.total.clustermarker$pct.2)
write.xlsx(th.total.clustermarker, "Lists/Clustermarker_th.total_clustermarker.xlsx")
#th.total.clustermarker <- read_excel("Lists/Clustermarker_th.total_clustermarker.xlsx")

#DotPlot Top up and down regulated genes 
#up
th.total.clustermarker_UP <- subset(th.total.clustermarker, Foldchange_UP >=2 )
top10_UP_clustermarker <- th.total.clustermarker_UP %>% group_by(cluster) %>% top_n(n=10,wt=Foldchange_UP)
top10_UP_clustermarker<-arrange(top10_UP_clustermarker,cluster,Foldchange_UP)
top10_UP_clustermarker <- as.character(top10_UP_clustermarker$gene)
top10_UP_clustermarker<-top10_UP_clustermarker[!duplicated(top10_UP_clustermarker)]

pdf("Graphs/Main/DotPlot_Top_10UP_Clustermarker_th.total_FCmin2.pdf", width = 7, height = 20)
DotPlot(th.total, features=top10_UP_clustermarker,assay="RNA") + coord_flip()+ scale_color_viridis(option = "H")+theme(axis.text.x = element_text(angle = 45)) + ggtitle("DEG_th.total_clustermarker") + scale_y_discrete(guide = guide_axis(angle = 90))
dev.off()

#Clustermarker
Idents(th.total)<-"tissue"
th.total.clustermarker.t<-FindAllMarkers(th.total, assay = "RNA", min.pct = 0.25, logfc.threshold = 0.25)
th.total.clustermarker.t$Foldchange_UP <- 2^(th.total.clustermarker.t$avg_log2FC)
th.total.clustermarker.t$Foldchange_DOWN <- 2^(-th.total.clustermarker.t$avg_log2FC)
th.total.clustermarker.t$Ratio_pct1_pct2 <- (th.total.clustermarker.t$pct.1)/(th.total.clustermarker.t$pct.2)
write.xlsx(th.total.clustermarker.t, "Lists/Clustermarker_th.total_tissue.xlsx")
#th.total.clustermarker.t <- read_excel("Lists/Clustermarker_th.total_tissue.xlsx")

#DotPlot Top up and down regulated genes 
#up
th.total.clustermarker.t_UP <- subset(th.total.clustermarker.t, Foldchange_UP >=2 )
top10_UP_clustermarker <- th.total.clustermarker.t_UP %>% group_by(cluster) %>% top_n(n=10,wt=Foldchange_UP)
top10_UP_clustermarker<-arrange(top10_UP_clustermarker,cluster,Foldchange_UP)
top10_UP_clustermarker <- as.character(top10_UP_clustermarker$gene)
top10_UP_clustermarker<-top10_UP_clustermarker[!duplicated(top10_UP_clustermarker)]

pdf("Graphs/Main/DotPlot_Top_10UP_Clustermarker_th.total_tissue.pdf", width = 7, height = 12)
DotPlot(th.total, features=top10_UP_clustermarker, group.by = "tissue",assay="RNA") + coord_flip()+ scale_color_viridis(option = "H")+theme(axis.text.x = element_text(angle = 45)) + ggtitle("DEG_th.total_tissue") + scale_y_discrete(guide = guide_axis(angle = 90))
dev.off()


#Clustermarker
Idents(th.total)<-"condition"
th.total.clustermarker.cond<-FindAllMarkers(th.total, assay = "RNA", min.pct = 0.25, logfc.threshold = 0.25)
th.total.clustermarker.cond$Foldchange_UP <- 2^(th.total.clustermarker.cond$avg_log2FC)
th.total.clustermarker.cond$Foldchange_DOWN <- 2^(-th.total.clustermarker.cond$avg_log2FC)
th.total.clustermarker.cond$Ratio_pct1_pct2 <- (th.total.clustermarker.cond$pct.1)/(th.total.clustermarker.cond$pct.2)
write.xlsx(th.total.clustermarker.cond, "Lists/Clustermarker_th.total_condition.xlsx")

#DotPlot Top up and down regulated genes 
#up
th.total.clustermarker.t_UP <- subset(th.total.clustermarker.t, Foldchange_UP >=2 )
top10_UP_clustermarker <- th.total.clustermarker.t_UP %>% group_by(cluster) %>% top_n(n=10,wt=Foldchange_UP)
top10_UP_clustermarker<-arrange(top10_UP_clustermarker,cluster,Foldchange_UP)
top10_UP_clustermarker <- as.character(top10_UP_clustermarker$gene)
top10_UP_clustermarker<-top10_UP_clustermarker[!duplicated(top10_UP_clustermarker)]

pdf("Graphs/Main/DotPlot_Top_10UP_Clustermarker_th.total_condition.pdf", width = 7, height = 12)
DotPlot(th.total, features=top10_UP_clustermarker, group.by = "condition",assay="RNA") + coord_flip()+ scale_color_viridis(option = "H")+theme(axis.text.x = element_text(angle = 45)) + ggtitle("DEG_th.total_condition") + scale_y_discrete(guide = guide_axis(angle = 90))
dev.off()


#Clustermarker celltype_tissue

#add identity to analyze tissue/health diffgenes
#Set tissue in identy 
Idents(th.total) <- th.total$celltype.full
th.total$celltype.tissue <- paste(Idents(th.total), th.total$tissue, sep = "_")
Idents(th.total)<-th.total$celltype.tissue

Idents(th.total)<-"celltype.tissue"
th.total.clustermarker.c.t<-FindAllMarkers(th.total, assay = "RNA", min.pct = 0.7, logfc.threshold = 1.5)
th.total.clustermarker.c.t$Foldchange_UP <- 2^(th.total.clustermarker.c.t$avg_log2FC)
th.total.clustermarker.c.t$Foldchange_DOWN <- 2^(-th.total.clustermarker.c.t$avg_log2FC)
th.total.clustermarker.c.t$Ratio_pct1_pct2 <- (th.total.clustermarker.c.t$pct.1)/(th.total.clustermarker.c.t$pct.2)
write.xlsx(th.total.clustermarker.c.t, "Lists/Clustermarker_th.total_celltype_condition.xlsx")

#DotPlot Top up and down regulated genes 
#up
th.total.clustermarker.c.t_UP <- subset(th.total.clustermarker.c.t, Foldchange_UP >=2 )
top10_UP_clustermarker <- th.total.clustermarker.c.t_UP %>% group_by(cluster) %>% top_n(n=10,wt=Foldchange_UP)
top10_UP_clustermarker<-arrange(top10_UP_clustermarker,cluster,Foldchange_UP)
top10_UP_clustermarker <- as.character(top10_UP_clustermarker$gene)
top10_UP_clustermarker<-top10_UP_clustermarker[!duplicated(top10_UP_clustermarker)]

pdf("Graphs/Main/DotPlot_Top_10UP_Clustermarker.c.t_th.total_FCmin2_SC.pdf", width = 10, height = 35)
DotPlot(th.total, features=top10_UP_clustermarker,  cols = c("blue","orange"), dot.scale = 8, assay="RNA") + coord_flip()+theme(axis.text.x = element_text(angle = 45)) + ggtitle("celltype.condition_up DEG") + scale_y_discrete(guide = guide_axis(angle = 90))
dev.off()

#down
th.total.clustermarker.c.t_DOWN <- subset(th.total.clustermarker.c.t, Foldchange_DOWN >=2 )
top10_DOWN_clustermarker <- th.total.clustermarker.c.t_DOWN %>% group_by(cluster) %>% top_n(n=10,wt=Foldchange_DOWN)
top10_DOWN_clustermarker<-arrange(top10_DOWN_clustermarker,cluster,Foldchange_DOWN)
top10_DOWN_clustermarker <- as.character(top10_DOWN_clustermarker$gene)
top10_DOWN_clustermarker<-top10_DOWN_clustermarker[!duplicated(top10_DOWN_clustermarker)]

pdf("Graphs/Main/DotPlot_Top_10DOWN_Clustermarker.c.t_th.total_FCmin2_SC.pdf", width = 10, height = 35)
DotPlot(th.total, features=top10_DOWN_clustermarker,  cols = c("blue","orange"), dot.scale = 8, assay="RNA") + coord_flip()+theme(axis.text.x = element_text(angle = 45)) + ggtitle("celltype.condition_down DEG") + scale_y_discrete(guide = guide_axis(angle = 90))
dev.off()


#(31)_Bar_plots_identifications####
pdf("Graphs/Main/Barplots.pdf", width=13, height=10)
plot_integrated_clusters_group (th.total)
plot_integrated_clusters_Donor_ID (th.total)
plot_integrated_celltypes_condition (th.total)
plot_integrated_celltypes_classification (th.total)
plot_integrated_celltypes_tissue (th.total)
plot_integrated_celltypes_sex (th.total)
plot_integrated_celltypes_age (th.total) 
plot_integrated_celltypes_development (th.total)
plot_integrated_celltypes_donor (th.total)
plot_main_integrated_condition_celltypes(th.total)
plot_main_integrated_development_celltypes(th.total)
plot_main_integrated_group_clusters(th.total)
plot_main_integrated_tissue_celltypes(th.total)
plot_main_integrated_group_tissue_clusters(th.total)
dev.off()

#save th.total
save(th.total,file="th.total.RData")


#(32)_cellchat####
#https://htmlpreview.github.io/?https://github.com/jinworks/CellChat/blob/master/tutorial/Comparison_analysis_of_multiple_datasets.html
#https://htmlpreview.github.io/?https://github.com/jinworks/CellChat/blob/master/tutorial/Comparison_analysis_of_multiple_datasets_with_different_cellular_compositions.html

Idents(th.total)<-th.total$tissue
levels(th.total)
table(th.total$tissue)
table(th.total$celltype)


#(33)_cc adult Thymus####
aduthymcellchat <- subset(th.total, idents = "adult_Thymus")
Idents(aduthymcellchat)<- aduthymcellchat$celltype
table(aduthymcellchat$celltype)

# Create cell Chat object, use RNA assay as negative values can be found in the integrated data
aduthymcellchat <- createCellChat(aduthymcellchat, group.by = "ident",assay = "RNA")

# load CellChat Database is a manually curated database of literature-supported ligand-receptor interactions
#in both human and mouse. CellChatDB in mouse contains 2,021 validated molecular interactions, including 60%
#of secrete autocrine/paracrine signaling interactions, 21% of extracellular matrix (ECM)-receptor interactions
#and 19% of cell-cell contact interactions. CellChatDB in human contains 1,939 validated molecular interactions,
#including 61.8% of paracrine/autocrine signaling interactions, 21.7% of extracellular matrix (ECM)-receptor interactions
#and 16.5% of cell-cell contact interactions.

#If needed one can manually add ligand receptor pairs to the database, for tutorial see to
#https://htmlpreview.github.io/?https://github.com/sqjin/CellChat/blob/master/tutorial/CellChat-vignette.html#load-the-required-libraries

CellChatDB <- CellChatDB.human

CellChatDB.use <- CellChatDB
aduthymcellchat@DB<- CellChatDB.use

aduthymcellchat<- subsetData(aduthymcellchat)


# identify overexpressed genes
aduthymcellchat <- identifyOverExpressedGenes(aduthymcellchat)
aduthymcellchat <- identifyOverExpressedInteractions(aduthymcellchat)

# compute communication probability
aduthymcellchat <- computeCommunProb(aduthymcellchat)
aduthymcellchat <- filterCommunication(aduthymcellchat, min.cells = 20)

aduthymcellchat <- computeCommunProbPathway(aduthymcellchat)
aduthymcellchat <- aggregateNet(aduthymcellchat)

df.net <- subsetCommunication(aduthymcellchat)
df.netP <- subsetCommunication(aduthymcellchat,slot.name = "netP")
write.xlsx(df.net, "Lists/df_net_aduthymus.xlsx")
write.xlsx(df.netP, "Lists/df_netP_aduthymus.xlsx")

groupSize <- as.numeric(table(aduthymcellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(aduthymcellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(aduthymcellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

netVisual_aggregate(aduthymcellchat, signaling =c("COLLAGEN"), layout = "chord")

save(aduthymcellchat,file="aduthymcellchat.RData")

####
pedthymcellchat <- subset(th.total, idents = "pediatric_Thymus")
Idents(pedthymcellchat)<- pedthymcellchat$celltype
table(pedthymcellchat$celltype)

# Create cell Chat object, use RNA assay as negative values can be found in the integrated data
pedthymcellchat <- createCellChat(pedthymcellchat, group.by = "ident",assay = "RNA")

# load CellChat Database is a manually curated database of literature-supported ligand-receptor interactions
#in both human and mouse. CellChatDB in mouse contains 2,021 validated molecular interactions, including 60%
#of secrete autocrine/paracrine signaling interactions, 21% of extracellular matrix (ECM)-receptor interactions
#and 19% of cell-cell contact interactions. CellChatDB in human contains 1,939 validated molecular interactions,
#including 61.8% of paracrine/autocrine signaling interactions, 21.7% of extracellular matrix (ECM)-receptor interactions
#and 16.5% of cell-cell contact interactions.

#If needed one can manually add ligand receptor pairs to the database, for tutorial see to
#https://htmlpreview.github.io/?https://github.com/sqjin/CellChat/blob/master/tutorial/CellChat-vignette.html#load-the-required-libraries

CellChatDB <- CellChatDB.human

CellChatDB.use <- CellChatDB
pedthymcellchat@DB<- CellChatDB.use

pedthymcellchat<- subsetData(pedthymcellchat)


# identify overexpressed genes
pedthymcellchat <- identifyOverExpressedGenes(pedthymcellchat)
pedthymcellchat <- identifyOverExpressedInteractions(pedthymcellchat)

# compute communication probability
pedthymcellchat <- computeCommunProb(pedthymcellchat)
pedthymcellchat <- filterCommunication(pedthymcellchat, min.cells = 20)

pedthymcellchat <- computeCommunProbPathway(pedthymcellchat)
pedthymcellchat <- aggregateNet(pedthymcellchat)

df.net <- subsetCommunication(pedthymcellchat)
df.netP <- subsetCommunication(pedthymcellchat,slot.name = "netP")
write.xlsx(df.net, "Lists/df_net_pedthymus.xlsx")
write.xlsx(df.netP, "Lists/df_netP_pedthymus.xlsx")

groupSize <- as.numeric(table(pedthymcellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(pedthymcellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(pedthymcellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

netVisual_aggregate(pedthymcellchat, signaling =c("COLLAGEN"), layout = "chord")

save(pedthymcellchat,file="pedthymcellchat.RData")


#(34)_cc prenatal Thymus vs adult Thymus####
Idents(th.total)<- th.total$tissue
table(th.total$tissue)
prethym_cellchat <- subset(th.total, idents = "prenatal_Thymus")
Idents(prethym_cellchat)<- prethym_cellchat$celltype
table(prethym_cellchat$celltype)

# Create cell Chat object, use RNA assay as negative values can be found in the integrated data
prethym_cellchat <- createCellChat(prethym_cellchat, group.by = "ident",assay = "RNA")

# load CellChat Database is a manually curated database of literature-supported ligand-receptor interactions
#in both human and mouse. CellChatDB in mouse contains 2,021 validated molecular interactions, including 60%
#of secrete autocrine/paracrine signaling interactions, 21% of extracellular matrix (ECM)-receptor interactions
#and 19% of cell-cell contact interactions. CellChatDB in human contains 1,939 validated molecular interactions,
#including 61.8% of paracrine/autocrine signaling interactions, 21.7% of extracellular matrix (ECM)-receptor interactions
#and 16.5% of cell-cell contact interactions.

#If needed one can manually add ligand receptor pairs to the database, for tutorial see to
#https://htmlpreview.github.io/?https://github.com/sqjin/CellChat/blob/master/tutorial/CellChat-vignette.html#load-the-required-libraries

CellChatDB <- CellChatDB.human

CellChatDB.use <- CellChatDB
prethym_cellchat@DB<- CellChatDB.use

prethym_cellchat<- subsetData(prethym_cellchat)


# identify overexpressed genes
prethym_cellchat <- identifyOverExpressedGenes(prethym_cellchat)
prethym_cellchat <- identifyOverExpressedInteractions(prethym_cellchat)

# compute communication probability
prethym_cellchat <- computeCommunProb(prethym_cellchat)
prethym_cellchat <- filterCommunication(prethym_cellchat, min.cells = 0)

prethym_cellchat <- computeCommunProbPathway(prethym_cellchat)
prethym_cellchat <- aggregateNet(prethym_cellchat)

df.net <- subsetCommunication(prethym_cellchat)
df.netP <- subsetCommunication(prethym_cellchat,slot.name = "netP")
write.xlsx(df.net, "Lists/df_net_prethym.xlsx")
write.xlsx(df.netP, "Lists/df_netP_prethym.xlsx")

save(prethym_cellchat,file="prethym_cellchat.RData")


###Comparative analysis#
#https://htmlpreview.github.io/?https://github.com/jinworks/CellChat/blob/master/tutorial/Comparison_analysis_of_multiple_datasets_with_different_cellular_compositions.html
#180324

#Load CellChat object of each dataset

#Lift up CellChat objects and merge them together
# Define the cell labels to lift up by combining both cell labels from the conditions
# Please note that the order of cell groups in `group.new` will affect the appearance order when visualizing the cell-cell communication. If there are unique cell groups in both the conditions, you should define `group.new = union(levels(cellchat.E14@idents), levels(cellchat.E13@idents))`
levels(prethym_cellchat@idents)
levels(aduthymcellchat@idents)

prethym_cellchat <- netAnalysis_computeCentrality(prethym_cellchat)
aduthymcellchat <- netAnalysis_computeCentrality(aduthymcellchat)

prethym.object.list <- list(thym = aduthymcellchat, prethym = prethym_cellchat)
cellchat_prethym <- mergeCellChat(prethym.object.list, add.names = names(prethym.object.list))
cellchat_prethym

pdf("Graphs/Main/cellchat/NoInteractions_prethymvsThym.pdf")
gg1 <- compareInteractions(cellchat_prethym, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat_prethym, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2
dev.off

pdf("Graphs/Main/cellchat/Heatmap_prethymvsThym.pdf")
gg1 <- netVisual_heatmap(cellchat_prethym)
#> Do heatmap based on a merged object
gg2 <- netVisual_heatmap(cellchat_prethym, measure = "weight")
#> Do heatmap based on a merged object
gg1 + gg2
dev.off()

pdf("Graphs/Main/cellchat/Circle_Diffno_prethymvsThym.pdf")
par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat_prethym, weight.scale = T)
netVisual_diffInteraction(cellchat_prethym, weight.scale = T, measure = "weight")
dev.off()

weight.max <- getMaxWeight(prethym.object.list, attribute = c("idents","count"))
par(mfrow = c(1,2), xpd=TRUE)
pdf("Graphs/Main/cellchat/Circle_Interactions_prethymvsThym.pdf")
for (i in 1:length(prethym.object.list)) {
  netVisual_circle(prethym.object.list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(prethym.object.list)[i]))
}
dev.off()

num.link <- sapply(prethym.object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(prethym.object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(prethym.object.list[[i]], title = names(prethym.object.list)[i], weight.MinMax = weight.MinMax)
}

pdf("Graphs/Main/cellchat/Scatterplot_prethymvsThym.pdf")
patchwork::wrap_plots(plots = gg)
dev.off()

gg1 <- netAnalysis_signalingChanges_scatter(cellchat_prethym, idents.use = "TEC", signaling.exclude = "MIF")
gg2 <- netAnalysis_signalingChanges_scatter(cellchat_prethym, idents.use = "TEC", signaling.exclude = c("MIF"))

pdf("Graphs/Main/cellchat/Specificchanges_TEC_prethymvsThym.pdf")
patchwork::wrap_plots(plots = list(gg1,gg2))
dev.off()


library(ComplexHeatmap)
i = 1
# combining all the identified signaling pathways from different datasets 
pathway.union <- union(prethym.object.list[[i]]@netP$pathways, prethym.object.list[[i+1]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(prethym.object.list[[i]], pattern = "outgoing", signaling = pathway.union, title = names(prethym.object.list)[i], width = 5, height = 15)
ht2 = netAnalysis_signalingRole_heatmap(prethym.object.list[[i+1]], pattern = "outgoing", signaling = pathway.union, title = names(prethym.object.list)[i+1], width = 5, height = 15)
pdf("Graphs/Main/cellchat/Outgoing_signaling_pattern_prethymvsThym.pdf")
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
dev.off()

ht1 = netAnalysis_signalingRole_heatmap(prethym.object.list[[i]], pattern = "incoming", signaling = pathway.union, title = names(prethym.object.list)[i], width = 5, height = 15, color.heatmap = "GnBu")
ht2 = netAnalysis_signalingRole_heatmap(prethym.object.list[[i+1]], pattern = "incoming", signaling = pathway.union, title = names(prethym.object.list)[i+1], width = 5, height = 15, color.heatmap = "GnBu")
pdf("Graphs/Main/cellchat/Incoming_signaling_pattern_prethymvsThym.pdf")
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
dev.off()

ht1 = netAnalysis_signalingRole_heatmap(prethym.object.list[[i]], pattern = "all", signaling = pathway.union, title = names(prethym.object.list)[i], width = 5, height = 15, color.heatmap = "OrRd")
ht2 = netAnalysis_signalingRole_heatmap(prethym.object.list[[i+1]], pattern = "all", signaling = pathway.union, title = names(prethym.object.list)[i+1], width = 5, height = 15, color.heatmap = "OrRd")
pdf("Graphs/Main/cellchat/Overall_signaling_pattern_prethymvsThym.pdf")
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
dev.off()



pdf("Graphs/Main/cellchat/Dotplot_cellchat_prethymvsThym.pdf", height = 30, width = 30)
netVisual_bubble(cellchat_prethym, comparison = c(1, 2), angle.x = 45)
dev.off()
#> Comparing communications on a merged object

gg1 <- netVisual_bubble(cellchat_prethym,  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
gg2 <- netVisual_bubble(cellchat_prethym,  comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
pdf("Graphs/Main/cellchat/Signaling_dotplot_cellchat_prethymvsThym.pdf", height = 30, width = 50)
gg1 + gg2
dev.off()


#(35)_cc pediatric Thymus vs adult Thymus####
Idents(th.total)<- th.total$tissue
pedthym_cellchat <- subset(th.total, idents = "pediatric_Thymus")
Idents(pedthym_cellchat)<- pedthym_cellchat$celltype
table(pedthym_cellchat$celltype)

# Create cell Chat object, use RNA assay as negative values can be found in the integrated data
pedthym_cellchat <- createCellChat(pedthym_cellchat, group.by = "ident",assay = "RNA")

# load CellChat Database is a manually curated database of literature-supported ligand-receptor interactions
#in both human and mouse. CellChatDB in mouse contains 2,021 validated molecular interactions, including 60%
#of secrete autocrine/paracrine signaling interactions, 21% of extracellular matrix (ECM)-receptor interactions
#and 19% of cell-cell contact interactions. CellChatDB in human contains 1,939 validated molecular interactions,
#including 61.8% of paracrine/autocrine signaling interactions, 21.7% of extracellular matrix (ECM)-receptor interactions
#and 16.5% of cell-cell contact interactions.

#If needed one can manually add ligand receptor pairs to the database, for tutorial see to
#https://htmlpreview.github.io/?https://github.com/sqjin/CellChat/blob/master/tutorial/CellChat-vignette.html#load-the-required-libraries

CellChatDB <- CellChatDB.human

CellChatDB.use <- CellChatDB
pedthym_cellchat@DB<- CellChatDB.use

pedthym_cellchat<- subsetData(pedthym_cellchat)


# identify overexpressed genes
pedthym_cellchat <- identifyOverExpressedGenes(pedthym_cellchat)
pedthym_cellchat <- identifyOverExpressedInteractions(pedthym_cellchat)

# compute communication probability
pedthym_cellchat <- computeCommunProb(pedthym_cellchat)
pedthym_cellchat <- filterCommunication(pedthym_cellchat, min.cells = 0)

pedthym_cellchat <- computeCommunProbPathway(pedthym_cellchat)
pedthym_cellchat <- aggregateNet(pedthym_cellchat)

df.net <- subsetCommunication(pedthym_cellchat)
df.netP <- subsetCommunication(pedthym_cellchat,slot.name = "netP")
write.xlsx(df.net, "Lists/df_net_pedthym.xlsx")
write.xlsx(df.netP, "Lists/df_netP_pedthym.xlsx")

save(pedthym_cellchat,file="pedthym_cellchat.RData")


###Comparative analysis#
#https://htmlpreview.github.io/?https://github.com/jinworks/CellChat/blob/master/tutorial/Comparison_analysis_of_multiple_datasets_with_different_cellular_compositions.html
#180324

#Load CellChat object of each dataset

#Lift up CellChat objects and merge them together
# Define the cell labels to lift up by combining both cell labels from the conditions
# Please note that the order of cell groups in `group.new` will affect the appearance order when visualizing the cell-cell communication. If there are unique cell groups in both the conditions, you should define `group.new = union(levels(cellchat.E14@idents), levels(cellchat.E13@idents))`
levels(pedthym_cellchat@idents)
levels(aduthymcellchat@idents)

pedthym_cellchat <- netAnalysis_computeCentrality(pedthym_cellchat)
aduthymcellchat <- netAnalysis_computeCentrality(aduthymcellchat)

pedthym.object.list <- list(thym = aduthymcellchat, pedthym = pedthym_cellchat)
cellchat_pedthym <- mergeCellChat(pedthym.object.list, add.names = names(pedthym.object.list))
cellchat_pedthym

pdf("Graphs/Main/cellchat/NoInteractions_pedthymvsThym.pdf")
gg1 <- compareInteractions(cellchat_pedthym, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat_pedthym, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2
dev.off

pdf("Graphs/Main/cellchat/Heatmap_pedthymvsThym.pdf")
gg1 <- netVisual_heatmap(cellchat_pedthym)
#> Do heatmap based on a merged object
gg2 <- netVisual_heatmap(cellchat_pedthym, measure = "weight")
#> Do heatmap based on a merged object
gg1 + gg2
dev.off()

pdf("Graphs/Main/cellchat/Circle_Diffno_pedthymvsThym.pdf")
par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat_pedthym, weight.scale = T)
netVisual_diffInteraction(cellchat_pedthym, weight.scale = T, measure = "weight")
dev.off()

weight.max <- getMaxWeight(pedthym.object.list, attribute = c("idents","count"))
par(mfrow = c(1,2), xpd=TRUE)
pdf("Graphs/Main/cellchat/Circle_Interactions_pedthymvsThym.pdf")
for (i in 1:length(pedthym.object.list)) {
  netVisual_circle(pedthym.object.list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(pedthym.object.list)[i]))
}
dev.off()

num.link <- sapply(pedthym.object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(pedthym.object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(pedthym.object.list[[i]], title = names(pedthym.object.list)[i], weight.MinMax = weight.MinMax)
}

pdf("Graphs/Main/cellchat/Scatterplot_pedthymvsThym.pdf")
patchwork::wrap_plots(plots = gg)
dev.off()

gg1 <- netAnalysis_signalingChanges_scatter(cellchat_pedthym, idents.use = "TEC", signaling.exclude = "MIF")
gg2 <- netAnalysis_signalingChanges_scatter(cellchat_pedthym, idents.use = "TEC", signaling.exclude = c("MIF"))

pdf("Graphs/Main/cellchat/Specificchanges_TEC_pedthymvsThym.pdf")
patchwork::wrap_plots(plots = list(gg1,gg2))
dev.off()


library(ComplexHeatmap)
i = 1
# combining all the identified signaling pathways from different datasets 
pathway.union <- union(pedthym.object.list[[i]]@netP$pathways, pedthym.object.list[[i+1]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(pedthym.object.list[[i]], pattern = "outgoing", signaling = pathway.union, title = names(pedthym.object.list)[i], width = 5, height = 15)
ht2 = netAnalysis_signalingRole_heatmap(pedthym.object.list[[i+1]], pattern = "outgoing", signaling = pathway.union, title = names(pedthym.object.list)[i+1], width = 5, height = 15)
pdf("Graphs/Main/cellchat/Outgoing_signaling_pattern_pedthymvsThym.pdf")
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
dev.off()

ht1 = netAnalysis_signalingRole_heatmap(pedthym.object.list[[i]], pattern = "incoming", signaling = pathway.union, title = names(pedthym.object.list)[i], width = 5, height = 15, color.heatmap = "GnBu")
ht2 = netAnalysis_signalingRole_heatmap(pedthym.object.list[[i+1]], pattern = "incoming", signaling = pathway.union, title = names(pedthym.object.list)[i+1], width = 5, height = 15, color.heatmap = "GnBu")
pdf("Graphs/Main/cellchat/Incoming_signaling_pattern_pedthymvsThym.pdf")
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
dev.off()

ht1 = netAnalysis_signalingRole_heatmap(pedthym.object.list[[i]], pattern = "all", signaling = pathway.union, title = names(pedthym.object.list)[i], width = 5, height = 15, color.heatmap = "OrRd")
ht2 = netAnalysis_signalingRole_heatmap(pedthym.object.list[[i+1]], pattern = "all", signaling = pathway.union, title = names(pedthym.object.list)[i+1], width = 5, height = 15, color.heatmap = "OrRd")
pdf("Graphs/Main/cellchat/Overall_signaling_pattern_pedthymvsThym.pdf")
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
dev.off()



pdf("Graphs/Main/cellchat/Dotplot_cellchat_pedthymvsThym.pdf", height = 30, width = 30)
netVisual_bubble(cellchat_pedthym, comparison = c(1, 2), angle.x = 45)
dev.off()
#> Comparing communications on a merged object

gg1 <- netVisual_bubble(cellchat_pedthym,  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
gg2 <- netVisual_bubble(cellchat_pedthym,  comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
pdf("Graphs/Main/cellchat/Signaling_dotplot_cellchat_pedthymvsThym.pdf", height = 30, width = 50)
gg1 + gg2
dev.off()


#(36)_cc TTH vs adult Thymus####
Idents(th.total)<- th.total$tissue
TTH_cellchat <- subset(th.total, idents = "TTH")
Idents(TTH_cellchat)<- TTH_cellchat$celltype
table(TTH_cellchat$celltype)

# Create cell Chat object, use RNA assay as negative values can be found in the integrated data
TTH_cellchat <- createCellChat(TTH_cellchat, group.by = "ident",assay = "RNA")

# load CellChat Database is a manually curated database of literature-supported ligand-receptor interactions
#in both human and mouse. CellChatDB in mouse contains 2,021 validated molecular interactions, including 60%
#of secrete autocrine/paracrine signaling interactions, 21% of extracellular matrix (ECM)-receptor interactions
#and 19% of cell-cell contact interactions. CellChatDB in human contains 1,939 validated molecular interactions,
#including 61.8% of paracrine/autocrine signaling interactions, 21.7% of extracellular matrix (ECM)-receptor interactions
#and 16.5% of cell-cell contact interactions.

#If needed one can manually add ligand receptor pairs to the database, for tutorial see to
#https://htmlpreview.github.io/?https://github.com/sqjin/CellChat/blob/master/tutorial/CellChat-vignette.html#load-the-required-libraries

CellChatDB <- CellChatDB.human

CellChatDB.use <- CellChatDB
TTH_cellchat@DB<- CellChatDB.use

TTH_cellchat<- subsetData(TTH_cellchat)


# identify overexpressed genes
TTH_cellchat <- identifyOverExpressedGenes(TTH_cellchat)
TTH_cellchat <- identifyOverExpressedInteractions(TTH_cellchat)

# compute communication probability
TTH_cellchat <- computeCommunProb(TTH_cellchat)
TTH_cellchat <- filterCommunication(TTH_cellchat, min.cells = 0)

TTH_cellchat <- computeCommunProbPathway(TTH_cellchat)
TTH_cellchat <- aggregateNet(TTH_cellchat)

df.net <- subsetCommunication(TTH_cellchat)
df.netP <- subsetCommunication(TTH_cellchat,slot.name = "netP")
write.xlsx(df.net, "Lists/df_net_TTH.xlsx")
write.xlsx(df.netP, "Lists/df_netP_TTH.xlsx")

save(TTH_cellchat,file="TTH_cellchat.RData")


###Comparative analysis#
#https://htmlpreview.github.io/?https://github.com/jinworks/CellChat/blob/master/tutorial/Comparison_analysis_of_multiple_datasets_with_different_cellular_compositions.html
#180324

#Load CellChat object of each dataset

#Lift up CellChat objects and merge them together
# Define the cell labels to lift up by combining both cell labels from the conditions
# Please note that the order of cell groups in `group.new` will affect the appearance order when visualizing the cell-cell communication. If there are unique cell groups in both the conditions, you should define `group.new = union(levels(cellchat.E14@idents), levels(cellchat.E13@idents))`
levels(TTH_cellchat@idents)
levels(aduthymcellchat@idents)

TTH_cellchat <- netAnalysis_computeCentrality(TTH_cellchat)
aduthymcellchat <- netAnalysis_computeCentrality(aduthymcellchat)

TTH.object.list <- list(thym = aduthymcellchat, TTH = TTH_cellchat)
cellchat_TTH <- mergeCellChat(TTH.object.list, add.names = names(TTH.object.list))
cellchat_TTH

pdf("Graphs/Main/cellchat/NoInteractions_TTHvsThym.pdf")
gg1 <- compareInteractions(cellchat_TTH, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat_TTH, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2
dev.off

pdf("Graphs/Main/cellchat/Heatmap_TTHvsThym.pdf")
gg1 <- netVisual_heatmap(cellchat_TTH)
#> Do heatmap based on a merged object
gg2 <- netVisual_heatmap(cellchat_TTH, measure = "weight")
#> Do heatmap based on a merged object
gg1 + gg2
dev.off()

pdf("Graphs/Main/cellchat/Circle_Diffno_TTHvsThym.pdf")
par(mfrow = c(1,2), xpd=TRUE)
#netVisual_diffInteraction(cellchat_TTH, weight.scale = T)
netVisual_diffInteraction(cellchat_TTH, weight.scale = T, measure = "weight")

dev.off()

weight.max <- getMaxWeight(TTH.object.list, attribute = c("idents","count"))
par(mfrow = c(1,2), xpd=TRUE)
pdf("Graphs/Main/cellchat/Circle_Interactions_TTHvsThym.pdf")
for (i in 1:length(TTH.object.list)) {
  netVisual_circle(TTH.object.list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(TTH.object.list)[i]))
}
dev.off()

num.link <- sapply(TTH.object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(TTH.object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(TTH.object.list[[i]], title = names(TTH.object.list)[i], weight.MinMax = weight.MinMax)
}

pdf("Graphs/Main/cellchat/Scatterplot_TTHvsThym.pdf")
patchwork::wrap_plots(plots = gg)
dev.off()

gg1 <- netAnalysis_signalingChanges_scatter(cellchat_TTH, idents.use = "TEC", signaling.exclude = "MIF")
gg2 <- netAnalysis_signalingChanges_scatter(cellchat_TTH, idents.use = "TEC", signaling.exclude = c("MIF"))

pdf("Graphs/Main/cellchat/Specificchanges_TEC_TTHvsThym.pdf")
patchwork::wrap_plots(plots = list(gg1,gg2))
dev.off()


library(ComplexHeatmap)
i = 1
# combining all the identified signaling pathways from different datasets 
pathway.union <- union(TTH.object.list[[i]]@netP$pathways, TTH.object.list[[i+1]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(TTH.object.list[[i]], pattern = "outgoing", signaling = pathway.union, title = names(TTH.object.list)[i], width = 5, height = 15)
ht2 = netAnalysis_signalingRole_heatmap(TTH.object.list[[i+1]], pattern = "outgoing", signaling = pathway.union, title = names(TTH.object.list)[i+1], width = 5, height = 15)
pdf("Graphs/Main/cellchat/Outgoing_signaling_pattern_TTHvsThym.pdf")
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
dev.off()

ht1 = netAnalysis_signalingRole_heatmap(TTH.object.list[[i]], pattern = "incoming", signaling = pathway.union, title = names(TTH.object.list)[i], width = 5, height = 15, color.heatmap = "GnBu")
ht2 = netAnalysis_signalingRole_heatmap(TTH.object.list[[i+1]], pattern = "incoming", signaling = pathway.union, title = names(TTH.object.list)[i+1], width = 5, height = 15, color.heatmap = "GnBu")
pdf("Graphs/Main/cellchat/Incoming_signaling_pattern_TTHvsThym.pdf")
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
dev.off()

ht1 = netAnalysis_signalingRole_heatmap(TTH.object.list[[i]], pattern = "all", signaling = pathway.union, title = names(TTH.object.list)[i], width = 5, height = 15, color.heatmap = "OrRd")
ht2 = netAnalysis_signalingRole_heatmap(TTH.object.list[[i+1]], pattern = "all", signaling = pathway.union, title = names(TTH.object.list)[i+1], width = 5, height = 15, color.heatmap = "OrRd")
pdf("Graphs/Main/cellchat/Overall_signaling_pattern_TTHvsThym.pdf")
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
dev.off()



pdf("Graphs/Main/cellchat/Dotplot_cellchat_TTHvsThym.pdf", height = 30, width = 30)
netVisual_bubble(cellchat_TTH, comparison = c(1, 2), angle.x = 45)
dev.off()
#> Comparing communications on a merged object

gg1 <- netVisual_bubble(cellchat_TTH,  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
gg2 <- netVisual_bubble(cellchat_TTH,  comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
pdf("Graphs/Main/cellchat/Signaling_dotplot_cellchat_TTHvsThym.pdf", height = 30, width = 50)
gg1 + gg2
dev.off()


#(37)_cc TET_A vs adult Thymus####
Idents(th.total)<- th.total$tissue
TET_A_cellchat <- subset(th.total, idents = "TET_A")
Idents(TET_A_cellchat)<- TET_A_cellchat$celltype
table(TET_A_cellchat$celltype)

# Create cell Chat object, use RNA assay as negative values can be found in the integrated data
TET_A_cellchat <- createCellChat(TET_A_cellchat, group.by = "ident",assay = "RNA")

# load CellChat Database is a manually curated database of literature-supported ligand-receptor interactions
#in both human and mouse. CellChatDB in mouse contains 2,021 validated molecular interactions, including 60%
#of secrete autocrine/paracrine signaling interactions, 21% of extracellular matrix (ECM)-receptor interactions
#and 19% of cell-cell contact interactions. CellChatDB in human contains 1,939 validated molecular interactions,
#including 61.8% of paracrine/autocrine signaling interactions, 21.7% of extracellular matrix (ECM)-receptor interactions
#and 16.5% of cell-cell contact interactions.

#If needed one can manually add ligand receptor pairs to the database, for tutorial see to
#https://htmlpreview.github.io/?https://github.com/sqjin/CellChat/blob/master/tutorial/CellChat-vignette.html#load-the-required-libraries

CellChatDB <- CellChatDB.human

CellChatDB.use <- CellChatDB
TET_A_cellchat@DB<- CellChatDB.use

TET_A_cellchat<- subsetData(TET_A_cellchat)


# identify overexpressed genes
TET_A_cellchat <- identifyOverExpressedGenes(TET_A_cellchat)
TET_A_cellchat <- identifyOverExpressedInteractions(TET_A_cellchat)

# compute communication probability
TET_A_cellchat <- computeCommunProb(TET_A_cellchat)
TET_A_cellchat <- filterCommunication(TET_A_cellchat, min.cells = 0)

TET_A_cellchat <- computeCommunProbPathway(TET_A_cellchat)
TET_A_cellchat <- aggregateNet(TET_A_cellchat)

df.net <- subsetCommunication(TET_A_cellchat)
df.netP <- subsetCommunication(TET_A_cellchat,slot.name = "netP")
write.xlsx(df.net, "Lists/df_net_TET_A.xlsx")
write.xlsx(df.netP, "Lists/df_netP_TET_A.xlsx")

save(TET_A_cellchat,file="TET_A_cellchatt.RData")

###Comparative analysis#
#https://htmlpreview.github.io/?https://github.com/jinworks/CellChat/blob/master/tutorial/Comparison_analysis_of_multiple_datasets_with_different_cellular_compositions.html
#180324

#Load CellChat object of each dataset

#Lift up CellChat objects and merge them together
# Define the cell labels to lift up by combining both cell labels from the conditions
# Please note that the order of cell groups in `group.new` will affect the appearance order when visualizing the cell-cell communication. If there are unique cell groups in both the conditions, you should define `group.new = union(levels(cellchat.E14@idents), levels(cellchat.E13@idents))`
levels(TET_A_cellchat@idents)
levels(aduthymcellchat@idents)

TET_A_cellchat <- netAnalysis_computeCentrality(TET_A_cellchat)
aduthymcellchat <- netAnalysis_computeCentrality(aduthymcellchat)

TET_A.object.list <- list(thym = aduthymcellchat, TET_A = TET_A_cellchat)
cellchat_TET_A <- mergeCellChat(TET_A.object.list, add.names = names(TET_A.object.list))
cellchat_TET_A

pdf("Graphs/Main/cellchat/NoInteractions_TET_AvsThym.pdf")
gg1 <- compareInteractions(cellchat_TET_A, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat_TET_A, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2
dev.off()

pdf("Graphs/Main/cellchat/Heatmap_TET_AvsThym.pdf")
gg1 <- netVisual_heatmap(cellchat_TET_A)
#> Do heatmap based on a merged object
gg2 <- netVisual_heatmap(cellchat_TET_A, measure = "weight")
#> Do heatmap based on a merged object
gg1 + gg2
dev.off()

pdf("Graphs/Main/cellchat/Circle_Diffno_TET_AvsThym.pdf")
par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat_TET_A, weight.scale = T)
netVisual_diffInteraction(cellchat_TET_A, weight.scale = T, measure = "weight")
dev.off()

weight.max <- getMaxWeight(TET_A.object.list, attribute = c("idents","count"))
par(mfrow = c(1,2), xpd=TRUE)
pdf("Graphs/Main/cellchat/Circle_Interactions_TET_AvsThym.pdf")
for (i in 1:length(TET_A.object.list)) {
  netVisual_circle(TET_A.object.list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(TET_A.object.list)[i]))
}
dev.off()

num.link <- sapply(TET_A.object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(TET_A.object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(TET_A.object.list[[i]], title = names(TET_A.object.list)[i], weight.MinMax = weight.MinMax)
}

pdf("Graphs/Main/cellchat/Scatterplot_TET_AvsThym.pdf")
patchwork::wrap_plots(plots = gg)
dev.off()

gg1 <- netAnalysis_signalingChanges_scatter(cellchat_TET_A, idents.use = "TEC", signaling.exclude = "MIF")
gg2 <- netAnalysis_signalingChanges_scatter(cellchat_TET_A, idents.use = "TEC", signaling.exclude = c("MIF"))

pdf("Graphs/Main/cellchat/Specificchanges_TEC_TET_AvsThym.pdf")
patchwork::wrap_plots(plots = list(gg1,gg2))
dev.off()


library(ComplexHeatmap)
i = 1
# combining all the identified signaling pathways from different datasets 
pathway.union <- union(TET_A.object.list[[i]]@netP$pathways, TET_A.object.list[[i+1]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(TET_A.object.list[[i]], pattern = "outgoing", signaling = pathway.union, title = names(TET_A.object.list)[i], width = 5, height = 15)
ht2 = netAnalysis_signalingRole_heatmap(TET_A.object.list[[i+1]], pattern = "outgoing", signaling = pathway.union, title = names(TET_A.object.list)[i+1], width = 5, height = 15)
pdf("Graphs/Main/cellchat/Outgoing_signaling_pattern_TET_AvsThym.pdf")
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
dev.off()

ht1 = netAnalysis_signalingRole_heatmap(TET_A.object.list[[i]], pattern = "incoming", signaling = pathway.union, title = names(TET_A.object.list)[i], width = 5, height = 15, color.heatmap = "GnBu")
ht2 = netAnalysis_signalingRole_heatmap(TET_A.object.list[[i+1]], pattern = "incoming", signaling = pathway.union, title = names(TET_A.object.list)[i+1], width = 5, height = 15, color.heatmap = "GnBu")
pdf("Graphs/Main/cellchat/Incoming_signaling_pattern_TET_AvsThym.pdf")
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
dev.off()

ht1 = netAnalysis_signalingRole_heatmap(TET_A.object.list[[i]], pattern = "all", signaling = pathway.union, title = names(TET_A.object.list)[i], width = 5, height = 15, color.heatmap = "OrRd")
ht2 = netAnalysis_signalingRole_heatmap(TET_A.object.list[[i+1]], pattern = "all", signaling = pathway.union, title = names(TET_A.object.list)[i+1], width = 5, height = 15, color.heatmap = "OrRd")
pdf("Graphs/Main/cellchat/Overall_signaling_pattern_TET_AvsThym.pdf")
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
dev.off()



pdf("Graphs/Main/cellchat/Dotplot_cellchat_TET_AvsThym.pdf", height = 30, width = 30)
netVisual_bubble(cellchat_TET_A, comparison = c(1, 2), angle.x = 45)
dev.off()
#> Comparing communications on a merged object

gg1 <- netVisual_bubble(cellchat_TET_A,  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
gg2 <- netVisual_bubble(cellchat_TET_A,  comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
pdf("Graphs/Main/cellchat/Signaling_dotplot_cellchat_TET_AvsThym.pdf", height = 30, width = 50)
gg1 + gg2
dev.off()







#(38)_cc TET_AB vs adult Thymus####
Idents(th.total)<- th.total$tissue
TET_AB_cellchat <- subset(th.total, idents = "TET_AB")
Idents(TET_AB_cellchat)<- TET_AB_cellchat$celltype
table(TET_AB_cellchat$celltype)

# Create cell Chat object, use RNA assay as negative values can be found in the integrated data
TET_AB_cellchat <- createCellChat(TET_AB_cellchat, group.by = "ident",assay = "RNA")

# load CellChat Database is a manually curated database of literature-supported ligand-receptor interactions
#in both human and mouse. CellChatDB in mouse contains 2,021 validated molecular interactions, including 60%
#of secrete autocrine/paracrine signaling interactions, 21% of extracellular matrix (ECM)-receptor interactions
#and 19% of cell-cell contact interactions. CellChatDB in human contains 1,939 validated molecular interactions,
#including 61.8% of paracrine/autocrine signaling interactions, 21.7% of extracellular matrix (ECM)-receptor interactions
#and 16.5% of cell-cell contact interactions.

#If needed one can manually add ligand receptor pairs to the database, for tutorial see to
#https://htmlpreview.github.io/?https://github.com/sqjin/CellChat/blob/master/tutorial/CellChat-vignette.html#load-the-required-libraries

CellChatDB <- CellChatDB.human

CellChatDB.use <- CellChatDB
TET_AB_cellchat@DB<- CellChatDB.use

TET_AB_cellchat<- subsetData(TET_AB_cellchat)


# identify overexpressed genes
TET_AB_cellchat <- identifyOverExpressedGenes(TET_AB_cellchat)
TET_AB_cellchat <- identifyOverExpressedInteractions(TET_AB_cellchat)

# compute communication probability
TET_AB_cellchat <- computeCommunProb(TET_AB_cellchat)
TET_AB_cellchat <- filterCommunication(TET_AB_cellchat, min.cells = 0)

TET_AB_cellchat <- computeCommunProbPathway(TET_AB_cellchat)
TET_AB_cellchat <- aggregateNet(TET_AB_cellchat)

df.net <- subsetCommunication(TET_AB_cellchat)
df.netP <- subsetCommunication(TET_AB_cellchat,slot.name = "netP")
write.xlsx(df.net, "Lists/df_net_TET_AB.xlsx")
write.xlsx(df.netP, "Lists/df_netP_TET_AB.xlsx")

save(TET_AB_cellchat,file="TET_AB_cellchat.RData")

###Comparative analysis#
#https://htmlpreview.github.io/?https://github.com/jinworks/CellChat/blob/master/tutorial/Comparison_analysis_of_multiple_datasets_with_different_cellular_compositions.html
#180324

#Load CellChat object of each dataset

#Lift up CellChat objects and merge them together
# Define the cell labels to lift up by combining both cell labels from the conditions
# Please note that the order of cell groups in `group.new` will affect the appearance order when visualizing the cell-cell communication. If there are unique cell groups in both the conditions, you should define `group.new = union(levels(cellchat.E14@idents), levels(cellchat.E13@idents))`
levels(TET_AB_cellchat@idents)
levels(aduthymcellchat@idents)

TET_AB_cellchat <- netAnalysis_computeCentrality(TET_AB_cellchat)
aduthymcellchat <- netAnalysis_computeCentrality(aduthymcellchat)

TET_AB.object.list <- list(thym = aduthymcellchat, TET_AB = TET_AB_cellchat)
cellchat_TET_AB <- mergeCellChat(TET_AB.object.list, add.names = names(TET_AB.object.list))
cellchat_TET_AB

pdf("Graphs/Main/cellchat/NoInteractions_TET_ABvsThym.pdf")
gg1 <- compareInteractions(cellchat_TET_AB, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat_TET_AB, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2
dev.off()

pdf("Graphs/Main/cellchat/Heatmap_TET_ABvsThym.pdf")
gg1 <- netVisual_heatmap(cellchat_TET_AB)
#> Do heatmap based on a merged object
gg2 <- netVisual_heatmap(cellchat_TET_AB, measure = "weight")
#> Do heatmap based on a merged object
gg1 + gg2
dev.off()

pdf("Graphs/Main/cellchat/Circle_Diffno_TET_ABvsThym.pdf")
par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat_TET_AB, weight.scale = T)
netVisual_diffInteraction(cellchat_TET_AB, weight.scale = T, measure = "weight")
dev.off()

weight.max <- getMaxWeight(TET_AB.object.list, attribute = c("idents","count"))
par(mfrow = c(1,2), xpd=TRUE)
pdf("Graphs/Main/cellchat/Circle_Interactions_TET_ABvsThym.pdf")
for (i in 1:length(TET_AB.object.list)) {
  netVisual_circle(TET_AB.object.list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(TET_AB.object.list)[i]))
}
dev.off()

num.link <- sapply(TET_AB.object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(TET_AB.object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(TET_AB.object.list[[i]], title = names(TET_AB.object.list)[i], weight.MinMax = weight.MinMax)
}

pdf("Graphs/Main/cellchat/Scatterplot_TET_ABvsThym.pdf")
patchwork::wrap_plots(plots = gg)
dev.off()

gg1 <- netAnalysis_signalingChanges_scatter(cellchat_TET_AB, idents.use = "TEC", signaling.exclude = "MIF")
gg2 <- netAnalysis_signalingChanges_scatter(cellchat_TET_AB, idents.use = "TEC", signaling.exclude = c("MIF"))

pdf("Graphs/Main/cellchat/Specificchanges_TEC_TET_ABvsThym.pdf")
patchwork::wrap_plots(plots = list(gg1,gg2))
dev.off()


library(ComplexHeatmap)
i = 1
# combining all the identified signaling pathways from different datasets 
pathway.union <- union(TET_AB.object.list[[i]]@netP$pathways, TET_AB.object.list[[i+1]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(TET_AB.object.list[[i]], pattern = "outgoing", signaling = pathway.union, title = names(TET_AB.object.list)[i], width = 5, height = 15)
ht2 = netAnalysis_signalingRole_heatmap(TET_AB.object.list[[i+1]], pattern = "outgoing", signaling = pathway.union, title = names(TET_AB.object.list)[i+1], width = 5, height = 15)
pdf("Graphs/Main/cellchat/Outgoing_signaling_pattern_TET_ABvsThym.pdf")
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
dev.off()

ht1 = netAnalysis_signalingRole_heatmap(TET_AB.object.list[[i]], pattern = "incoming", signaling = pathway.union, title = names(TET_AB.object.list)[i], width = 5, height = 15, color.heatmap = "GnBu")
ht2 = netAnalysis_signalingRole_heatmap(TET_AB.object.list[[i+1]], pattern = "incoming", signaling = pathway.union, title = names(TET_AB.object.list)[i+1], width = 5, height = 15, color.heatmap = "GnBu")
pdf("Graphs/Main/cellchat/Incoming_signaling_pattern_TET_ABvsThym.pdf")
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
dev.off()

ht1 = netAnalysis_signalingRole_heatmap(TET_AB.object.list[[i]], pattern = "all", signaling = pathway.union, title = names(TET_AB.object.list)[i], width = 5, height = 15, color.heatmap = "OrRd")
ht2 = netAnalysis_signalingRole_heatmap(TET_AB.object.list[[i+1]], pattern = "all", signaling = pathway.union, title = names(TET_AB.object.list)[i+1], width = 5, height = 15, color.heatmap = "OrRd")
pdf("Graphs/Main/cellchat/Overall_signaling_pattern_TET_ABvsThym.pdf")
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
dev.off()



pdf("Graphs/Main/cellchat/Dotplot_cellchat_TET_ABvsThym.pdf", height = 30, width = 30)
netVisual_bubble(cellchat_TET_AB, comparison = c(1, 2), angle.x = 45)
dev.off()
#> Comparing communications on a merged object

gg1 <- netVisual_bubble(cellchat_TET_AB,  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
gg2 <- netVisual_bubble(cellchat_TET_AB,  comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
pdf("Graphs/Main/cellchat/Signaling_dotplot_cellchat_TET_ABvsThym.pdf", height = 30, width = 50)
gg1 + gg2
dev.off()





#(39)_cc TET_B vs adult Thymus####
Idents(th.total)<- th.total$tissue
TET_B_cellchat <- subset(th.total, idents = "TET_B")
Idents(TET_B_cellchat)<- TET_B_cellchat$celltype
table(TET_B_cellchat$celltype)

# Create cell Chat object, use RNA assay as negative values can be found in the integrated data
TET_B_cellchat <- createCellChat(TET_B_cellchat, group.by = "ident",assay = "RNA")

# load CellChat Database is a manually curated database of literature-supported ligand-receptor interactions
#in both human and mouse. CellChatDB in mouse contains 2,021 validated molecular interactions, including 60%
#of secrete autocrine/paracrine signaling interactions, 21% of extracellular matrix (ECM)-receptor interactions
#and 19% of cell-cell contact interactions. CellChatDB in human contains 1,939 validated molecular interactions,
#including 61.8% of paracrine/autocrine signaling interactions, 21.7% of extracellular matrix (ECM)-receptor interactions
#and 16.5% of cell-cell contact interactions.

#If needed one can manually add ligand receptor pairs to the database, for tutorial see to
#https://htmlpreview.github.io/?https://github.com/sqjin/CellChat/blob/master/tutorial/CellChat-vignette.html#load-the-required-libraries

CellChatDB <- CellChatDB.human

CellChatDB.use <- CellChatDB
TET_B_cellchat@DB<- CellChatDB.use

TET_B_cellchat<- subsetData(TET_B_cellchat)


# identify overexpressed genes
TET_B_cellchat <- identifyOverExpressedGenes(TET_B_cellchat)
TET_B_cellchat <- identifyOverExpressedInteractions(TET_B_cellchat)

# compute communication probability
TET_B_cellchat <- computeCommunProb(TET_B_cellchat)
TET_B_cellchat <- filterCommunication(TET_B_cellchat, min.cells = 0)

TET_B_cellchat <- computeCommunProbPathway(TET_B_cellchat)
TET_B_cellchat <- aggregateNet(TET_B_cellchat)

df.net <- subsetCommunication(TET_B_cellchat)
df.netP <- subsetCommunication(TET_B_cellchat,slot.name = "netP")
write.xlsx(df.net, "Lists/df_net_TET_B.xlsx")
write.xlsx(df.netP, "Lists/df_netP_TET_B.xlsx")

save(TET_B_cellchat,file="TET_B_cellchat.RData")

###Comparative analysis#
#https://htmlpreview.github.io/?https://github.com/jinworks/CellChat/blob/master/tutorial/Comparison_analysis_of_multiple_datasets_with_different_cellular_compositions.html
#180324

#Load CellChat object of each dataset

#Lift up CellChat objects and merge them together
# Define the cell labels to lift up by combining both cell labels from the conditions
# Please note that the order of cell groups in `group.new` will affect the appearance order when visualizing the cell-cell communication. If there are unique cell groups in both the conditions, you should define `group.new = union(levels(cellchat.E14@idents), levels(cellchat.E13@idents))`
levels(TET_B_cellchat@idents)
levels(aduthymcellchat@idents)

TET_B_cellchat <- netAnalysis_computeCentrality(TET_B_cellchat)
aduthymcellchat <- netAnalysis_computeCentrality(aduthymcellchat)

TET_B.object.list <- list(thym = aduthymcellchat, TET_B = TET_B_cellchat)
cellchat_TET_B <- mergeCellChat(TET_B.object.list, add.names = names(TET_B.object.list))
cellchat_TET_B

pdf("Graphs/Main/cellchat/NoInteractions_TET_BvsThym.pdf")
gg1 <- compareInteractions(cellchat_TET_B, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat_TET_B, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2
dev.off()

pdf("Graphs/Main/cellchat/Heatmap_TET_BvsThym.pdf")
gg1 <- netVisual_heatmap(cellchat_TET_B)
#> Do heatmap based on a merged object
gg2 <- netVisual_heatmap(cellchat_TET_B, measure = "weight")
#> Do heatmap based on a merged object
gg1 + gg2
dev.off()

pdf("Graphs/Main/cellchat/Circle_Diffno_TET_BvsThym.pdf")
par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat_TET_B, weight.scale = T)
netVisual_diffInteraction(cellchat_TET_B, weight.scale = T, measure = "weight")
dev.off()

weight.max <- getMaxWeight(TET_B.object.list, attribute = c("idents","count"))
par(mfrow = c(1,2), xpd=TRUE)
pdf("Graphs/Main/cellchat/Circle_Interactions_TET_BvsThym.pdf")
for (i in 1:length(TET_B.object.list)) {
  netVisual_circle(TET_B.object.list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(TET_B.object.list)[i]))
}
dev.off()

num.link <- sapply(TET_B.object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(TET_B.object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(TET_B.object.list[[i]], title = names(TET_B.object.list)[i], weight.MinMax = weight.MinMax)
}

pdf("Graphs/Main/cellchat/Scatterplot_TET_BvsThym.pdf")
patchwork::wrap_plots(plots = gg)
dev.off()

gg1 <- netAnalysis_signalingChanges_scatter(cellchat_TET_B, idents.use = "TEC", signaling.exclude = "MIF")
gg2 <- netAnalysis_signalingChanges_scatter(cellchat_TET_B, idents.use = "TEC", signaling.exclude = c("MIF"))

pdf("Graphs/Main/cellchat/Specificchanges_TEC_TET_BvsThym.pdf")
patchwork::wrap_plots(plots = list(gg1,gg2))
dev.off()


library(ComplexHeatmap)
i = 1
# combining all the identified signaling pathways from different datasets 
pathway.union <- union(TET_B.object.list[[i]]@netP$pathways, TET_B.object.list[[i+1]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(TET_B.object.list[[i]], pattern = "outgoing", signaling = pathway.union, title = names(TET_B.object.list)[i], width = 5, height = 15)
ht2 = netAnalysis_signalingRole_heatmap(TET_B.object.list[[i+1]], pattern = "outgoing", signaling = pathway.union, title = names(TET_B.object.list)[i+1], width = 5, height = 15)
pdf("Graphs/Main/cellchat/Outgoing_signaling_pattern_TET_BvsThym.pdf")
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
dev.off()

ht1 = netAnalysis_signalingRole_heatmap(TET_B.object.list[[i]], pattern = "incoming", signaling = pathway.union, title = names(TET_B.object.list)[i], width = 5, height = 15, color.heatmap = "GnBu")
ht2 = netAnalysis_signalingRole_heatmap(TET_B.object.list[[i+1]], pattern = "incoming", signaling = pathway.union, title = names(TET_B.object.list)[i+1], width = 5, height = 15, color.heatmap = "GnBu")
pdf("Graphs/Main/cellchat/Incoming_signaling_pattern_TET_BvsThym.pdf")
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
dev.off()

ht1 = netAnalysis_signalingRole_heatmap(TET_B.object.list[[i]], pattern = "all", signaling = pathway.union, title = names(TET_B.object.list)[i], width = 5, height = 15, color.heatmap = "OrRd")
ht2 = netAnalysis_signalingRole_heatmap(TET_B.object.list[[i+1]], pattern = "all", signaling = pathway.union, title = names(TET_B.object.list)[i+1], width = 5, height = 15, color.heatmap = "OrRd")
pdf("Graphs/Main/cellchat/Overall_signaling_pattern_TET_BvsThym.pdf")
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
dev.off()



pdf("Graphs/Main/cellchat/Dotplot_cellchat_TET_BvsThym.pdf", height = 30, width = 30)
netVisual_bubble(cellchat_TET_B, comparison = c(1, 2), angle.x = 45)
dev.off()
#> Comparing communications on a merged object

gg1 <- netVisual_bubble(cellchat_TET_B,  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
gg2 <- netVisual_bubble(cellchat_TET_B,  comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
pdf("Graphs/Main/cellchat/Signaling_dotplot_cellchat_TET_BvsThym.pdf", height = 30, width = 50)
gg1 + gg2
dev.off()

#(40)_cc TET_C vs adult Thymus####
Idents(th.total)<- th.total$tissue
TET_C_cellchat <- subset(th.total, idents = "TET_C")
Idents(TET_C_cellchat)<- TET_C_cellchat$celltype
table(TET_C_cellchat$celltype)

# Create cell Chat object, use RNA assay as negative values can be found in the integrated data
TET_C_cellchat <- createCellChat(TET_C_cellchat, group.by = "ident",assay = "RNA")

# load CellChat Database is a manually curated database of literature-supported ligand-receptor interactions
#in both human and mouse. CellChatDB in mouse contains 2,021 validated molecular interactions, including 60%
#of secrete autocrine/paracrine signaling interactions, 21% of extracellular matrix (ECM)-receptor interactions
#and 19% of cell-cell contact interactions. CellChatDB in human contains 1,939 validated molecular interactions,
#including 61.8% of paracrine/autocrine signaling interactions, 21.7% of extracellular matrix (ECM)-receptor interactions
#and 16.5% of cell-cell contact interactions.

#If needed one can manually add ligand receptor pairs to the database, for tutorial see to
#https://htmlpreview.github.io/?https://github.com/sqjin/CellChat/blob/master/tutorial/CellChat-vignette.html#load-the-required-libraries

CellChatDB <- CellChatDB.human

CellChatDB.use <- CellChatDB
TET_C_cellchat@DB<- CellChatDB.use

TET_C_cellchat<- subsetData(TET_C_cellchat)


# identify overexpressed genes
TET_C_cellchat <- identifyOverExpressedGenes(TET_C_cellchat)
TET_C_cellchat <- identifyOverExpressedInteractions(TET_C_cellchat)

# compute communication probability
TET_C_cellchat <- computeCommunProb(TET_C_cellchat)
TET_C_cellchat <- filterCommunication(TET_C_cellchat, min.cells = 0)

TET_C_cellchat <- computeCommunProbPathway(TET_C_cellchat)
TET_C_cellchat <- aggregateNet(TET_C_cellchat)

df.net <- subsetCommunication(TET_C_cellchat)
df.netP <- subsetCommunication(TET_C_cellchat,slot.name = "netP")
write.xlsx(df.net, "Lists/df_net_TET_C.xlsx")
write.xlsx(df.netP, "Lists/df_netP_TET_C.xlsx")

save(TET_C_cellchat,file="TET_C_cellchat.RData")

###Comparative analysis#
#https://htmlpreview.github.io/?https://github.com/jinworks/CellChat/blob/master/tutorial/Comparison_analysis_of_multiple_datasets_with_different_cellular_compositions.html
#180324

#Load CellChat object of each dataset

#Lift up CellChat objects and merge them together
# Define the cell labels to lift up by combining both cell labels from the conditions
# Please note that the order of cell groups in `group.new` will affect the appearance order when visualizing the cell-cell communication. If there are unique cell groups in both the conditions, you should define `group.new = union(levels(cellchat.E14@idents), levels(cellchat.E13@idents))`
levels(TET_C_cellchat@idents)
levels(aduthymcellchat@idents)

TET_C_cellchat <- netAnalysis_computeCentrality(TET_C_cellchat)
aduthymcellchat <- netAnalysis_computeCentrality(aduthymcellchat)

TET_C.object.list <- list(thym = aduthymcellchat, TET_C = TET_C_cellchat)
cellchat_TET_C <- mergeCellChat(TET_C.object.list, add.names = names(TET_C.object.list))
cellchat_TET_C

pdf("Graphs/Main/cellchat/NoInteractions_TET_CvsThym.pdf")
gg1 <- compareInteractions(cellchat_TET_C, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat_TET_C, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2
dev.off()

pdf("Graphs/Main/cellchat/Heatmap_TET_CvsThym.pdf")
gg1 <- netVisual_heatmap(cellchat_TET_C)
#> Do heatmap based on a merged object
gg2 <- netVisual_heatmap(cellchat_TET_C, measure = "weight")
#> Do heatmap based on a merged object
gg1 + gg2
dev.off()

pdf("Graphs/Main/cellchat/Circle_Diffno_TET_CvsThym.pdf")
par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat_TET_C, weight.scale = T)
netVisual_diffInteraction(cellchat_TET_C, weight.scale = T, measure = "weight")
dev.off()

weight.max <- getMaxWeight(TET_C.object.list, attribute = c("idents","count"))
par(mfrow = c(1,2), xpd=TRUE)
pdf("Graphs/Main/cellchat/Circle_Interactions_TET_CvsThym.pdf")
for (i in 1:length(TET_C.object.list)) {
  netVisual_circle(TET_C.object.list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(TET_C.object.list)[i]))
}
dev.off()

num.link <- sapply(TET_C.object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(TET_C.object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(TET_C.object.list[[i]], title = names(TET_C.object.list)[i], weight.MinMax = weight.MinMax)
}

pdf("Graphs/Main/cellchat/Scatterplot_TET_CvsThym.pdf")
patchwork::wrap_plots(plots = gg)
dev.off()

gg1 <- netAnalysis_signalingChanges_scatter(cellchat_TET_C, idents.use = "TEC", signaling.exclude = "MIF")
gg2 <- netAnalysis_signalingChanges_scatter(cellchat_TET_C, idents.use = "TEC", signaling.exclude = c("MIF"))

pdf("Graphs/Main/cellchat/Specificchanges_TEC_TET_CvsThym.pdf")
patchwork::wrap_plots(plots = list(gg1,gg2))
dev.off()


library(ComplexHeatmap)
i = 1
# combining all the identified signaling pathways from different datasets 
pathway.union <- union(TET_C.object.list[[i]]@netP$pathways, TET_C.object.list[[i+1]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(TET_C.object.list[[i]], pattern = "outgoing", signaling = pathway.union, title = names(TET_C.object.list)[i], width = 5, height = 15)
ht2 = netAnalysis_signalingRole_heatmap(TET_C.object.list[[i+1]], pattern = "outgoing", signaling = pathway.union, title = names(TET_C.object.list)[i+1], width = 5, height = 15)
pdf("Graphs/Main/cellchat/Outgoing_signaling_pattern_TET_CvsThym.pdf")
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
dev.off()

ht1 = netAnalysis_signalingRole_heatmap(TET_C.object.list[[i]], pattern = "incoming", signaling = pathway.union, title = names(TET_C.object.list)[i], width = 5, height = 15, color.heatmap = "GnBu")
ht2 = netAnalysis_signalingRole_heatmap(TET_C.object.list[[i+1]], pattern = "incoming", signaling = pathway.union, title = names(TET_C.object.list)[i+1], width = 5, height = 15, color.heatmap = "GnBu")
pdf("Graphs/Main/cellchat/Incoming_signaling_pattern_TET_CvsThym.pdf")
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
dev.off()

ht1 = netAnalysis_signalingRole_heatmap(TET_C.object.list[[i]], pattern = "all", signaling = pathway.union, title = names(TET_C.object.list)[i], width = 5, height = 15, color.heatmap = "OrRd")
ht2 = netAnalysis_signalingRole_heatmap(TET_C.object.list[[i+1]], pattern = "all", signaling = pathway.union, title = names(TET_C.object.list)[i+1], width = 5, height = 15, color.heatmap = "OrRd")
pdf("Graphs/Main/cellchat/Overall_signaling_pattern_TET_CvsThym.pdf")
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
dev.off()



pdf("Graphs/Main/cellchat/Dotplot_cellchat_TET_CvsThym.pdf", height = 30, width = 30)
netVisual_bubble(cellchat_TET_C, comparison = c(1, 2), angle.x = 45)
dev.off()
#> Comparing communications on a merged object

gg1 <- netVisual_bubble(cellchat_TET_C,  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
gg2 <- netVisual_bubble(cellchat_TET_C,  comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
pdf("Graphs/Main/cellchat/Signaling_dotplot_cellchat_TET_CvsThym.pdf", height = 30, width = 50)
gg1 + gg2
dev.off()



#(41)_cc TET_MNT vs adult Thymus####
Idents(th.total)<- th.total$tissue
TET_MNT_cellchat <- subset(th.total, idents = "TET_MNT")
Idents(TET_MNT_cellchat)<- TET_MNT_cellchat$celltype
table(TET_MNT_cellchat$celltype)

# Create cell Chat object, use RNA assay as negative values can be found in the integrated data
TET_MNT_cellchat <- createCellChat(TET_MNT_cellchat, group.by = "ident",assay = "RNA")

# load CellChat Database is a manually curated database of literature-supported ligand-receptor interactions
#in both human and mouse. CellChatDB in mouse contains 2,021 validated molecular interactions, including 60%
#of secrete autocrine/paracrine signaling interactions, 21% of extracellular matrix (ECM)-receptor interactions
#and 19% of cell-cell contact interactions. CellChatDB in human contains 1,939 validated molecular interactions,
#including 61.8% of paracrine/autocrine signaling interactions, 21.7% of extracellular matrix (ECM)-receptor interactions
#and 16.5% of cell-cell contact interactions.

#If needed one can manually add ligand receptor pairs to the database, for tutorial see to
#https://htmlpreview.github.io/?https://github.com/sqjin/CellChat/blob/master/tutorial/CellChat-vignette.html#load-the-required-libraries

CellChatDB <- CellChatDB.human

CellChatDB.use <- CellChatDB
TET_MNT_cellchat@DB<- CellChatDB.use

TET_MNT_cellchat<- subsetData(TET_MNT_cellchat)


# identify overexpressed genes
TET_MNT_cellchat <- identifyOverExpressedGenes(TET_MNT_cellchat)
TET_MNT_cellchat <- identifyOverExpressedInteractions(TET_MNT_cellchat)

# compute communication probability
TET_MNT_cellchat <- computeCommunProb(TET_MNT_cellchat)
TET_MNT_cellchat <- filterCommunication(TET_MNT_cellchat, min.cells = 0)

TET_MNT_cellchat <- computeCommunProbPathway(TET_MNT_cellchat)
TET_MNT_cellchat <- aggregateNet(TET_MNT_cellchat)

df.net <- subsetCommunication(TET_MNT_cellchat)
df.netP <- subsetCommunication(TET_MNT_cellchat,slot.name = "netP")
write.xlsx(df.net, "Lists/df_net_TET_MNT.xlsx")
write.xlsx(df.netP, "Lists/df_netP_TET_MNT.xlsx")

save(TET_MNT_cellchat,file="TET_MNT_cellchat.RData")

###Comparative analysis#
#https://htmlpreview.github.io/?https://github.com/jinworks/CellChat/blob/master/tutorial/Comparison_analysis_of_multiple_datasets_with_different_cellular_compositions.html
#180324

#Load CellChat object of each dataset

#Lift up CellChat objects and merge them together
# Define the cell labels to lift up by combining both cell labels from the conditions
# Please note that the order of cell groups in `group.new` will affect the appearance order when visualizing the cell-cell communication. If there are unique cell groups in both the conditions, you should define `group.new = union(levels(cellchat.E14@idents), levels(cellchat.E13@idents))`
levels(TET_MNT_cellchat@idents)
levels(aduthymcellchat@idents)

TET_MNT_cellchat <- netAnalysis_computeCentrality(TET_MNT_cellchat)
aduthymcellchat <- netAnalysis_computeCentrality(aduthymcellchat)

TET_MNT.object.list <- list(thym = aduthymcellchat, TET_MNT = TET_MNT_cellchat)
cellchat_TET_MNT <- mergeCellChat(TET_MNT.object.list, add.names = names(TET_MNT.object.list))
cellchat_TET_MNT

pdf("Graphs/Main/cellchat/NoInteractions_TET_MNTvsThym.pdf")
gg1 <- compareInteractions(cellchat_TET_MNT, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat_TET_MNT, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2
dev.off()

pdf("Graphs/Main/cellchat/Heatmap_TET_MNTvsThym.pdf")
gg1 <- netVisual_heatmap(cellchat_TET_MNT)
#> Do heatmap based on a merged object
gg2 <- netVisual_heatmap(cellchat_TET_MNT, measure = "weight")
#> Do heatmap based on a merged object
gg1 + gg2
dev.off()

pdf("Graphs/Main/cellchat/Circle_Diffno_TET_MNTvsThym.pdf")
par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat_TET_MNT, weight.scale = T)
netVisual_diffInteraction(cellchat_TET_MNT, weight.scale = T, measure = "weight")
dev.off()

weight.max <- getMaxWeight(TET_MNT.object.list, attribute = c("idents","count"))
par(mfrow = c(1,2), xpd=TRUE)
pdf("Graphs/Main/cellchat/Circle_Interactions_TET_MNTvsThym.pdf")
for (i in 1:length(TET_MNT.object.list)) {
  netVisual_circle(TET_MNT.object.list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(TET_MNT.object.list)[i]))
}
dev.off()

num.link <- sapply(TET_MNT.object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(TET_MNT.object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(TET_MNT.object.list[[i]], title = names(TET_MNT.object.list)[i], weight.MinMax = weight.MinMax)
}

pdf("Graphs/Main/cellchat/Scatterplot_TET_MNTvsThym.pdf")
patchwork::wrap_plots(plots = gg)
dev.off()

gg1 <- netAnalysis_signalingChanges_scatter(cellchat_TET_MNT, idents.use = "TEC", signaling.exclude = "MIF")
gg2 <- netAnalysis_signalingChanges_scatter(cellchat_TET_MNT, idents.use = "TEC", signaling.exclude = c("MIF"))

pdf("Graphs/Main/cellchat/Specificchanges_TEC_TET_MNTvsThym.pdf")
patchwork::wrap_plots(plots = list(gg1,gg2))
dev.off()


library(ComplexHeatmap)
i = 1
# combining all the identified signaling pathways from different datasets 
pathway.union <- union(TET_MNT.object.list[[i]]@netP$pathways, TET_MNT.object.list[[i+1]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(TET_MNT.object.list[[i]], pattern = "outgoing", signaling = pathway.union, title = names(TET_MNT.object.list)[i], width = 5, height = 15)
ht2 = netAnalysis_signalingRole_heatmap(TET_MNT.object.list[[i+1]], pattern = "outgoing", signaling = pathway.union, title = names(TET_MNT.object.list)[i+1], width = 5, height = 15)
pdf("Graphs/Main/cellchat/Outgoing_signaling_pattern_TET_MNTvsThym.pdf")
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
dev.off()

ht1 = netAnalysis_signalingRole_heatmap(TET_MNT.object.list[[i]], pattern = "incoming", signaling = pathway.union, title = names(TET_MNT.object.list)[i], width = 5, height = 15, color.heatmap = "GnBu")
ht2 = netAnalysis_signalingRole_heatmap(TET_MNT.object.list[[i+1]], pattern = "incoming", signaling = pathway.union, title = names(TET_MNT.object.list)[i+1], width = 5, height = 15, color.heatmap = "GnBu")
pdf("Graphs/Main/cellchat/Incoming_signaling_pattern_TET_MNTvsThym.pdf")
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
dev.off()

ht1 = netAnalysis_signalingRole_heatmap(TET_MNT.object.list[[i]], pattern = "all", signaling = pathway.union, title = names(TET_MNT.object.list)[i], width = 5, height = 15, color.heatmap = "OrRd")
ht2 = netAnalysis_signalingRole_heatmap(TET_MNT.object.list[[i+1]], pattern = "all", signaling = pathway.union, title = names(TET_MNT.object.list)[i+1], width = 5, height = 15, color.heatmap = "OrRd")
pdf("Graphs/Main/cellchat/Overall_signaling_pattern_TET_MNTvsThym.pdf")
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
dev.off()



pdf("Graphs/Main/cellchat/Dotplot_cellchat_TET_MNTvsThym.pdf", height = 30, width = 30)
netVisual_bubble(cellchat_TET_MNT, comparison = c(1, 2), angle.x = 45)
dev.off()
#> Comparing communications on a merged object

gg1 <- netVisual_bubble(cellchat_TET_MNT,  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
gg2 <- netVisual_bubble(cellchat_TET_MNT,  comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
pdf("Graphs/Main/cellchat/Signaling_dotplot_cellchat_TET_MNTvsThym.pdf", height = 30, width = 50)
gg1 + gg2
dev.off()


#(42)_cellchat_comparison_all####
#Load CellChat object of each dataset

#Lift up CellChat objects and merge them together
# Define the cell labels to lift up by combining both cell labels from the conditions
# Please note that the order of cell groups in `group.new` will affect the appearance order when visualizing the cell-cell communication. If there are unique cell groups in both the conditions, you should define `group.new = union(levels(cellchat.E14@idents), levels(cellchat.E13@idents))`
levels(TET_A_cellchat@idents)
levels(TET_AB_cellchat@idents)
levels(TET_B_cellchat@idents)
levels(TET_C_cellchat@idents)
levels(TET_MNT_cellchat@idents)
levels(TTH_cellchat@idents)
levels(prethym_cellchat@idents)
levels(pedthym_cellchat@idents)
levels(aduthymcellchat@idents)

TET_A_cellchat <- netAnalysis_computeCentrality(TET_A_cellchat)
TET_AB_cellchat <- netAnalysis_computeCentrality(TET_AB_cellchat)
TET_B_cellchat <- netAnalysis_computeCentrality(TET_B_cellchat)
TET_C_cellchat <- netAnalysis_computeCentrality(TET_C_cellchat)
TET_MNT_cellchat <- netAnalysis_computeCentrality(TET_MNT_cellchat)
TTH_cellchat <- netAnalysis_computeCentrality(TTH_cellchat)
prethym_cellchat <- netAnalysis_computeCentrality(prethym_cellchat)
pedthym_cellchat <- netAnalysis_computeCentrality(pedthym_cellchat)
aduthymcellchat <- netAnalysis_computeCentrality(aduthymcellchat)

TET_all.object.list <- list(aduthym = aduthymcellchat, pedthym = pedthym_cellchat, prethym = prethym_cellchat, TTH = TTH_cellchat, TET_A = TET_A_cellchat, TET_AB = TET_AB_cellchat, TET_B = TET_B_cellchat, TET_C = TET_C_cellchat, TET_MNT = TET_MNT_cellchat)
cellchat_TET_all <- mergeCellChat(TET_all.object.list, add.names = names(TET_all.object.list))
cellchat_TET_all

pdf("Graphs/Main/cellchat/NoInteractions_TET_allvsThym.pdf", height=3, width=3)
gg1 <- compareInteractions(cellchat_TET_all, show.legend = F, group = c(1:9), color.use = color_tissue )
gg2 <- compareInteractions(cellchat_TET_all, show.legend = F, group = c(1:9), measure = "weight", color.use = color_tissue)
gg1 + gg2
dev.off()

pdf("Graphs/Main/cellchat/Heatmap_cellchat_TET_allvsThym_nointeract.pdf", width=15, height=3)
gg1 <- netVisual_heatmap(cellchat_TET_all,
                         comparison = c(1, 2))
gg2 <- netVisual_heatmap(cellchat_TET_all,
                         comparison = c(1, 3))
gg3 <- netVisual_heatmap(cellchat_TET_all,
                         comparison = c(1, 4))
gg4 <- netVisual_heatmap(cellchat_TET_all,
                         comparison = c(1, 5))
gg5 <- netVisual_heatmap(cellchat_TET_all,
                         comparison = c(1, 6))
gg6 <- netVisual_heatmap(cellchat_TET_all,
                         comparison = c(1, 7))
gg7 <- netVisual_heatmap(cellchat_TET_all,
                         comparison = c(1, 8))
gg8 <- netVisual_heatmap(cellchat_TET_all,
                         comparison = c(1, 9))

#> Do heatmap based on a merged object
gg1 + gg2+gg3+gg4+gg5+gg6+gg7+gg8
dev.off()

pdf("Graphs/Main/cellchat/Heatmap_cellchat_TET_allvsThym_strength_interact1.pdf", width=15, height=3)
#> Do heatmap based on a merged object
gg7 <- netVisual_heatmap(cellchat_TET_all, measure = "weight",
                         comparison = c(1, 2))
gg8 <- netVisual_heatmap(cellchat_TET_all, measure = "weight",
                         comparison = c(1, 3))
gg9 <- netVisual_heatmap(cellchat_TET_all, measure = "weight",
                         comparison = c(1, 4))
gg10 <- netVisual_heatmap(cellchat_TET_all, measure = "weight",
                          comparison = c(1, 5))
gg11 <- netVisual_heatmap(cellchat_TET_all, measure = "weight",
                          comparison = c(1, 6))
gg12 <- netVisual_heatmap(cellchat_TET_all, measure = "weight",
                          comparison = c(1, 7))
gg13 <- netVisual_heatmap(cellchat_TET_all, measure = "weight",
                          comparison = c(1, 8))
gg14 <- netVisual_heatmap(cellchat_TET_all, measure = "weight",
                          comparison = c(1, 9))
gg7+gg8+gg9+gg10+gg11+gg12+gg13+gg14
dev.off()


#pdf("Graphs/Main/cellchat/Circle_Diffno_cellchat_TET_allvsThym.pdf")
#par(mfrow = c(1,2), xpd=TRUE)
#netVisual_diffInteraction(cellchat_TET_all, weight.scale = T)
#netVisual_diffInteraction(cellchat_TET_all, weight.scale = T, measure = "weight")
#dev.off()

weight.max <- getMaxWeight(TET_all.object.list, attribute = c("idents","count"))
par(mfrow = c(1,9), xpd=TRUE)
pdf("Graphs/Main/cellchat/Circle_Interactions_cellchat_TET_allvsThym_labelF.pdf", width=4, height = 6)
for (i in 1:length(TET_all.object.list)) {
  netVisual_circle(TET_all.object.list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(TET_all.object.list)[i]), color.use = color_level1)
}
dev.off()

num.link <- sapply(TET_all.object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(TET_all.object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(TET_all.object.list[[i]], color.use = color_level1, title = names(TET_all.object.list)[i], weight.MinMax = weight.MinMax)+ xlim(0, 10)+ylim(0, 10)
}

pdf("Graphs/Main/cellchat/Scatterplot_cellchat_TET_allvsThym.pdf", width=10, height=10)
patchwork::wrap_plots(plots = gg, nrows=4, ncols=2)
dev.off()


library(ComplexHeatmap)
i = 1
# combining all the identified signaling pathways from different datasets 
pathway.union <- union(TET_all.object.list[[i]]@netP$pathways, TET_all.object.list[[i+1]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(TET_all.object.list[[i]], pattern = "outgoing", signaling = pathway.union, title = names(TET_all.object.list)[i], width = 5, height = 15)
ht2 = netAnalysis_signalingRole_heatmap(TET_all.object.list[[i+1]], pattern = "outgoing", signaling = pathway.union, title = names(TET_all.object.list)[i+1], width = 5, height = 15)
ht3 = netAnalysis_signalingRole_heatmap(TET_all.object.list[[i+2]], pattern = "outgoing", signaling = pathway.union, title = names(TET_all.object.list)[i+2], width = 5, height = 15)
ht4 = netAnalysis_signalingRole_heatmap(TET_all.object.list[[i+3]], pattern = "outgoing", signaling = pathway.union, title = names(TET_all.object.list)[i+3], width = 5, height = 15)
ht5 = netAnalysis_signalingRole_heatmap(TET_all.object.list[[i+4]], pattern = "outgoing", signaling = pathway.union, title = names(TET_all.object.list)[i+4], width = 5, height = 15)
ht6 = netAnalysis_signalingRole_heatmap(TET_all.object.list[[i+5]], pattern = "outgoing", signaling = pathway.union, title = names(TET_all.object.list)[i+5], width = 5, height = 15)
ht7 = netAnalysis_signalingRole_heatmap(TET_all.object.list[[i+6]], pattern = "outgoing", signaling = pathway.union, title = names(TET_all.object.list)[i+6], width = 5, height = 15)
ht8 = netAnalysis_signalingRole_heatmap(TET_all.object.list[[i+7]], pattern = "outgoing", signaling = pathway.union, title = names(TET_all.object.list)[i+7], width = 5, height = 15)
ht9 = netAnalysis_signalingRole_heatmap(TET_all.object.list[[i+8]], pattern = "outgoing", signaling = pathway.union, title = names(TET_all.object.list)[i+8], width = 5, height = 15)
pdf("Graphs/Main/cellchat/Outgoing_signaling_pattern_cellchat_TET_allvsThym.pdf", width=25, height=10)
draw(ht1 + ht2 + ht3 + ht4 + ht5 + ht6 + ht7 + ht8 + ht9, ht_gap = unit(0.5, "cm"))
dev.off()

ht1 = netAnalysis_signalingRole_heatmap(TET_all.object.list[[i]], pattern = "incoming", signaling = pathway.union, title = names(TET_all.object.list)[i], width = 5, height = 15, color.heatmap = "GnBu")
ht2 = netAnalysis_signalingRole_heatmap(TET_all.object.list[[i+1]], pattern = "incoming", signaling = pathway.union, title = names(TET_all.object.list)[i+1], width = 5, height = 15, color.heatmap = "GnBu")
ht3 = netAnalysis_signalingRole_heatmap(TET_all.object.list[[i+2]], pattern = "incoming", signaling = pathway.union, title = names(TET_all.object.list)[i+2], width = 5, height = 15, color.heatmap = "GnBu")
ht4 = netAnalysis_signalingRole_heatmap(TET_all.object.list[[i+3]], pattern = "incoming", signaling = pathway.union, title = names(TET_all.object.list)[i+3], width = 5, height = 15, color.heatmap = "GnBu")
ht5 = netAnalysis_signalingRole_heatmap(TET_all.object.list[[i+4]], pattern = "incoming", signaling = pathway.union, title = names(TET_all.object.list)[i+4], width = 5, height = 15, color.heatmap = "GnBu")
ht6 = netAnalysis_signalingRole_heatmap(TET_all.object.list[[i+5]], pattern = "incoming", signaling = pathway.union, title = names(TET_all.object.list)[i+5], width = 5, height = 15, color.heatmap = "GnBu")
ht7 = netAnalysis_signalingRole_heatmap(TET_all.object.list[[i+6]], pattern = "incoming", signaling = pathway.union, title = names(TET_all.object.list)[i+6], width = 5, height = 15, color.heatmap = "GnBu")
ht8 = netAnalysis_signalingRole_heatmap(TET_all.object.list[[i+7]], pattern = "incoming", signaling = pathway.union, title = names(TET_all.object.list)[i+7], width = 5, height = 15, color.heatmap = "GnBu")
ht9 = netAnalysis_signalingRole_heatmap(TET_all.object.list[[i+8]], pattern = "incoming", signaling = pathway.union, title = names(TET_all.object.list)[i+8], width = 5, height = 15, color.heatmap = "GnBu")
pdf("Graphs/Main/cellchat/Incoming_signaling_pattern_cellchat_TET_allvsThym.pdf", width=25, height=10)
draw(ht1 + ht2 + ht3 + ht4 + ht5 + ht6 + ht7 + ht8 + ht9, ht_gap = unit(0.5, "cm"))
dev.off()

ht1 = netAnalysis_signalingRole_heatmap(TET_all.object.list[[i]], pattern = "all", signaling = pathway.union, title = names(TET_all.object.list)[i], width = 5, height = 15, color.heatmap = "OrRd")
ht2 = netAnalysis_signalingRole_heatmap(TET_all.object.list[[i+1]], pattern = "all", signaling = pathway.union, title = names(TET_all.object.list)[i+1], width = 5, height = 15, color.heatmap = "OrRd")
ht3 = netAnalysis_signalingRole_heatmap(TET_all.object.list[[i+2]], pattern = "all", signaling = pathway.union, title = names(TET_all.object.list)[i+2], width = 5, height = 15, color.heatmap = "OrRd")
ht4 = netAnalysis_signalingRole_heatmap(TET_all.object.list[[i+3]], pattern = "all", signaling = pathway.union, title = names(TET_all.object.list)[i+3], width = 5, height = 15, color.heatmap = "OrRd")
ht5 = netAnalysis_signalingRole_heatmap(TET_all.object.list[[i+4]], pattern = "all", signaling = pathway.union, title = names(TET_all.object.list)[i+4], width = 5, height = 15, color.heatmap = "OrRd")
ht6 = netAnalysis_signalingRole_heatmap(TET_all.object.list[[i+5]], pattern = "all", signaling = pathway.union, title = names(TET_all.object.list)[i+5], width = 5, height = 15, color.heatmap = "OrRd")
ht7 = netAnalysis_signalingRole_heatmap(TET_all.object.list[[i+6]], pattern = "all", signaling = pathway.union, title = names(TET_all.object.list)[i+6], width = 5, height = 15, color.heatmap = "OrRd")
ht8 = netAnalysis_signalingRole_heatmap(TET_all.object.list[[i+7]], pattern = "all", signaling = pathway.union, title = names(TET_all.object.list)[i+7], width = 5, height = 15, color.heatmap = "OrRd")
ht9 = netAnalysis_signalingRole_heatmap(TET_all.object.list[[i+8]], pattern = "all", signaling = pathway.union, title = names(TET_all.object.list)[i+8], width = 5, height = 15, color.heatmap = "OrRd")
pdf("Graphs/Main/cellchat/Overall_signaling_pattern_cellchat_TET_allvsThym.pdf", width=25, height=10)
draw(ht1 + ht2 + ht3 + ht4 + ht5 + ht6 + ht7 + ht8 + ht9, ht_gap = unit(0.5, "cm"))
dev.off()



pdf("Graphs/Main/cellchat/Dotplot_cellchat_TET_allvsThym.pdf", height = 60, width = 100)
netVisual_bubble(cellchat_TET_all, comparison = c(1:9),color.text = color_level1, angle.x = 90)
dev.off()


pdf("Graphs/Main/cellchat/Signaling_dotplot_cellchat_TETallvsThym_cd99_APP.pdf", height = 3, width = 20)
netVisual_bubble(cellchat_TET_all,   comparison = c(1:9), signaling=c("CD99"), color.text = color_level1, title.name = "Increased signaling", angle.x = 90, remove.isolate = T)
netVisual_bubble(cellchat_TET_all,   comparison = c(1:9), signaling=c("APP"),max.dataset = 5, color.text = color_level1, title.name = "Increased signaling", angle.x = 90, remove.isolate = T)
dev.off()


#(43)_Stripplot_DEG_subset_tissue_specific####
Idents(th.total) <- th.total$celltype.full
th.total$celltype.tissue <- paste(Idents(th.total), th.total$tissue, sep = "_")
Idents(th.total)<-th.total$celltype.tissue

######prenatal_Thymus## #
#Clustermarker juxtaposition celltype.tissue against adult counterpart###
Idents(th.total)<-"celltype.tissue"
table(th.total$celltype.tissue)

TC_prenatal_Thymusvsadult_DEG       <- FindMarkers(th.total, ident.1= "TC_prenatal_Thymus", ident.2="TC_adult_Thymus",assay = "RNA", min.pct = 0.1)
BC_prenatal_Thymusvsadult_DEG       <- FindMarkers(th.total, ident.1= "BC_prenatal_Thymus", ident.2="BC_adult_Thymus",assay = "RNA", min.pct = 0.1)
PC_prenatal_Thymusvsadult_DEG       <- FindMarkers(th.total, ident.1= "PC_prenatal_Thymus", ident.2="PC_adult_Thymus",assay = "RNA", min.pct = 0.1)
TEC_prenatal_Thymusvsadult_DEG       <- FindMarkers(th.total, ident.1= "TEC_prenatal_Thymus", ident.2="TEC_adult_Thymus",assay = "RNA", min.pct = 0.1)
DC_prenatal_Thymusvsadult_DEG       <- FindMarkers(th.total, ident.1= "DC_prenatal_Thymus", ident.2="DC_adult_Thymus",assay = "RNA", min.pct = 0.1)
MAC_Mono_prenatal_Thymusvsadult_DEG       <- FindMarkers(th.total, ident.1= "MAC/Mono_prenatal_Thymus", ident.2="MAC/Mono_adult_Thymus",assay = "RNA", min.pct = 0.1)
EC_prenatal_Thymusvsadult_DEG       <- FindMarkers(th.total, ident.1= "EC_prenatal_Thymus", ident.2="EC_adult_Thymus",assay = "RNA", min.pct = 0.1)
FB_prenatal_Thymusvsadult_DEG       <- FindMarkers(th.total, ident.1= "FB_prenatal_Thymus", ident.2="FB_adult_Thymus",assay = "RNA", min.pct = 0.1)
VSMC_Peri_prenatal_Thymusvsadult_DEG       <- FindMarkers(th.total, ident.1= "VSMC/Peri_prenatal_Thymus", ident.2="VSMC/Peri_adult_Thymus",assay = "RNA", min.pct = 0.1)


TC_prenatal_Thymusvsadult_DEG$celltype <- "TC"
BC_prenatal_Thymusvsadult_DEG$celltype <- "BC"
PC_prenatal_Thymusvsadult_DEG$celltype <- "PC"
TEC_prenatal_Thymusvsadult_DEG$celltype <- "TEC"
DC_prenatal_Thymusvsadult_DEG$celltype <- "DC"
MAC_Mono_prenatal_Thymusvsadult_DEG$celltype <- "MAC/Mono"
EC_prenatal_Thymusvsadult_DEG$celltype <- "EC"
FB_prenatal_Thymusvsadult_DEG$celltype <- "FB"
VSMC_Peri_prenatal_Thymusvsadult_DEG$celltype <- "VSMC/Peri"


write.xlsx(TC_prenatal_Thymusvsadult_DEG, "Lists/Stripplot_lists/TC_prenatal_Thymusvsadult_DEG.xlsx")
write.xlsx(BC_prenatal_Thymusvsadult_DEG, "Lists/Stripplot_lists/BC_prenatal_Thymusvsadult_DEG.xlsx")
write.xlsx(PC_prenatal_Thymusvsadult_DEG, "Lists/Stripplot_lists/PC_prenatal_Thymusvsadult_DEG.xlsx")
write.xlsx(TEC_prenatal_Thymusvsadult_DEG, "Lists/Stripplot_lists/TEC_prenatal_Thymusvsadult_DEG.xlsx")
write.xlsx(DC_prenatal_Thymusvsadult_DEG, "Lists/Stripplot_lists/DC_prenatal_Thymusvsadult_DEG.xlsx")
write.xlsx(MAC_Mono_prenatal_Thymusvsadult_DEG, "Lists/Stripplot_lists/MAC_Mono_prenatal_Thymusvsadult_DEG.xlsx")
write.xlsx(EC_prenatal_Thymusvsadult_DEG, "Lists/Stripplot_lists/EC_prenatal_Thymusvsadult_DEG.xlsx")
write.xlsx(FB_prenatal_Thymusvsadult_DEG, "Lists/Stripplot_lists/FB_prenatal_Thymusvsadult_DEG.xlsx")
write.xlsx(VSMC_Peri_prenatal_Thymusvsadult_DEG, "Lists/Stripplot_lists/VSMC_Peri_prenatal_Thymusvsadult_DEG.xlsx")


###### Strip plot for DEG###

#1 Give each dataframe an additional column with the gene names the rownames might change later this will help you keep track

TC_prenatal_Thymusvsadult_DEG$gene<- rownames(TC_prenatal_Thymusvsadult_DEG)
BC_prenatal_Thymusvsadult_DEG$gene<- rownames(BC_prenatal_Thymusvsadult_DEG)
PC_prenatal_Thymusvsadult_DEG$gene<- rownames(PC_prenatal_Thymusvsadult_DEG)
TEC_prenatal_Thymusvsadult_DEG$gene<- rownames(TEC_prenatal_Thymusvsadult_DEG)
DC_prenatal_Thymusvsadult_DEG$gene<- rownames(DC_prenatal_Thymusvsadult_DEG)
MAC_Mono_prenatal_Thymusvsadult_DEG$gene<- rownames(MAC_Mono_prenatal_Thymusvsadult_DEG)
EC_prenatal_Thymusvsadult_DEG$gene<- rownames(EC_prenatal_Thymusvsadult_DEG)
FB_prenatal_Thymusvsadult_DEG$gene<- rownames(FB_prenatal_Thymusvsadult_DEG)
VSMC_Peri_prenatal_Thymusvsadult_DEG$gene<- rownames(VSMC_Peri_prenatal_Thymusvsadult_DEG)


#2 Create a list of dataframes with differentially expressed genes that you want to combine
DEG_strip_list <- list(TC_prenatal_Thymusvsadult_DEG, BC_prenatal_Thymusvsadult_DEG, PC_prenatal_Thymusvsadult_DEG, TEC_prenatal_Thymusvsadult_DEG, DC_prenatal_Thymusvsadult_DEG, MAC_Mono_prenatal_Thymusvsadult_DEG, EC_prenatal_Thymusvsadult_DEG, FB_prenatal_Thymusvsadult_DEG, VSMC_Peri_prenatal_Thymusvsadult_DEG)



# combine the dataframes into a single dataframe with an added column indicating the cell population
DEG_strip_df <- bind_rows(DEG_strip_list, .id = "celltype.tissue")


DEG_strip_df <- DEG_strip_df %>% mutate(Threshold_met =
                                          ifelse(abs(avg_log2FC) > 0.58 & p_val_adj < 0.05, 'p-value_adj & log2FC',
                                                 ifelse(abs(avg_log2FC) < 0.58 & p_val_adj < 0.05, 'p-value_adj',
                                                        ifelse(abs(avg_log2FC) > 0.58 & p_val_adj > 0.05, 'log2FC', 'Not sig.'))))




celltype_order <- c("TC","BC","PC","TEC","DC", "MAC/Mono","EC","FB","VSMC/Peri")



pdf("Graphs/Main/DEG_prenatal_Thymus_vs_adult_labeled.pdf", height= 5, width= 25)
ggplot(DEG_strip_df, aes(x = celltype, y = avg_log2FC, color = Threshold_met)) +
  geom_jitter(width = 0.25) +
  scale_color_manual(values = c('p-value_adj & log2FC'="#9c2162",
                                'p-value_adj' ="#0b9ed2",
                                'log2FC' ="#b2bd5b",
                                'Not sig.'="#A3A3A3"))  +
  theme_classic() +
  geom_text_repel(aes(label = ifelse(p_val_adj < 0.05 & abs(avg_log2FC) >= 8, gene, "")),
                  size = 4, max.overlaps = 8, fontface = "bold", point.padding = 0, box.padding = 0.2, min.segment.length = Inf,
                  segment.color = "grey", color = "black")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size=14),
        axis.text.y = element_text(face = "bold", size = 16),
        axis.line = element_line(size = 1.0),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold", margin = margin(r = 45)),
        plot.title = element_text(size = 16, face = "bold"),
        legend.text = element_text(size = 12, face = "bold"),
        legend.title = element_text(size = 14, face = "bold")) +
  ggtitle("DEG: prenatal_Thymus vs. Thym") +
  theme(plot.title = element_text(size = 16, face = "bold")) +
  scale_x_discrete(limits = celltype_order)
dev.off()

######pediatric_Thymus## #
#Clustermarker juxtaposition celltype.tissue against adult counterpart###
Idents(th.total)<-"celltype.tissue"
table(th.total$celltype.tissue)

TC_pediatric_Thymusvsadult_DEG       <- FindMarkers(th.total, ident.1= "TC_pediatric_Thymus", ident.2="TC_adult_Thymus",assay = "RNA", min.pct = 0.1)
BC_pediatric_Thymusvsadult_DEG       <- FindMarkers(th.total, ident.1= "BC_pediatric_Thymus", ident.2="BC_adult_Thymus",assay = "RNA", min.pct = 0.1)
PC_pediatric_Thymusvsadult_DEG       <- FindMarkers(th.total, ident.1= "PC_pediatric_Thymus", ident.2="PC_adult_Thymus",assay = "RNA", min.pct = 0.1)
TEC_pediatric_Thymusvsadult_DEG       <- FindMarkers(th.total, ident.1= "TEC_pediatric_Thymus", ident.2="TEC_adult_Thymus",assay = "RNA", min.pct = 0.1)
DC_pediatric_Thymusvsadult_DEG       <- FindMarkers(th.total, ident.1= "DC_pediatric_Thymus", ident.2="DC_adult_Thymus",assay = "RNA", min.pct = 0.1)
MAC_Mono_pediatric_Thymusvsadult_DEG       <- FindMarkers(th.total, ident.1= "MAC/Mono_pediatric_Thymus", ident.2="MAC/Mono_adult_Thymus",assay = "RNA", min.pct = 0.1)
EC_pediatric_Thymusvsadult_DEG       <- FindMarkers(th.total, ident.1= "EC_pediatric_Thymus", ident.2="EC_adult_Thymus",assay = "RNA", min.pct = 0.1)
FB_pediatric_Thymusvsadult_DEG       <- FindMarkers(th.total, ident.1= "FB_pediatric_Thymus", ident.2="FB_adult_Thymus",assay = "RNA", min.pct = 0.1)
VSMC_Peri_pediatric_Thymusvsadult_DEG       <- FindMarkers(th.total, ident.1= "VSMC/Peri_pediatric_Thymus", ident.2="VSMC/Peri_adult_Thymus",assay = "RNA", min.pct = 0.1)


TC_pediatric_Thymusvsadult_DEG$celltype <- "TC"
BC_pediatric_Thymusvsadult_DEG$celltype <- "BC"
PC_pediatric_Thymusvsadult_DEG$celltype <- "PC"
TEC_pediatric_Thymusvsadult_DEG$celltype <- "TEC"
DC_pediatric_Thymusvsadult_DEG$celltype <- "DC"
MAC_Mono_pediatric_Thymusvsadult_DEG$celltype <- "MAC/Mono"
EC_pediatric_Thymusvsadult_DEG$celltype <- "EC"
FB_pediatric_Thymusvsadult_DEG$celltype <- "FB"
VSMC_Peri_pediatric_Thymusvsadult_DEG$celltype <- "VSMC/Peri"


write.xlsx(TC_pediatric_Thymusvsadult_DEG, "Lists/Stripplot_lists/TC_pediatric_Thymusvsadult_DEG.xlsx")
write.xlsx(BC_pediatric_Thymusvsadult_DEG, "Lists/Stripplot_lists/BC_pediatric_Thymusvsadult_DEG.xlsx")
write.xlsx(PC_pediatric_Thymusvsadult_DEG, "Lists/Stripplot_lists/PC_pediatric_Thymusvsadult_DEG.xlsx")
write.xlsx(TEC_pediatric_Thymusvsadult_DEG, "Lists/Stripplot_lists/TEC_pediatric_Thymusvsadult_DEG.xlsx")
write.xlsx(DC_pediatric_Thymusvsadult_DEG, "Lists/Stripplot_lists/DC_pediatric_Thymusvsadult_DEG.xlsx")
write.xlsx(MAC_Mono_pediatric_Thymusvsadult_DEG, "Lists/Stripplot_lists/MAC_Mono_pediatric_Thymusvsadult_DEG.xlsx")
write.xlsx(EC_pediatric_Thymusvsadult_DEG, "Lists/Stripplot_lists/EC_pediatric_Thymusvsadult_DEG.xlsx")
write.xlsx(FB_pediatric_Thymusvsadult_DEG, "Lists/Stripplot_lists/FB_pediatric_Thymusvsadult_DEG.xlsx")
write.xlsx(VSMC_Peri_pediatric_Thymusvsadult_DEG, "Lists/Stripplot_lists/VSMC_Peri_pediatric_Thymusvsadult_DEG.xlsx")


###### Strip plot for DEG###

#1 Give each dataframe an additional column with the gene names the rownames might change later this will help you keep track

TC_pediatric_Thymusvsadult_DEG$gene<- rownames(TC_pediatric_Thymusvsadult_DEG)
BC_pediatric_Thymusvsadult_DEG$gene<- rownames(BC_pediatric_Thymusvsadult_DEG)
PC_pediatric_Thymusvsadult_DEG$gene<- rownames(PC_pediatric_Thymusvsadult_DEG)
TEC_pediatric_Thymusvsadult_DEG$gene<- rownames(TEC_pediatric_Thymusvsadult_DEG)
DC_pediatric_Thymusvsadult_DEG$gene<- rownames(DC_pediatric_Thymusvsadult_DEG)
MAC_Mono_pediatric_Thymusvsadult_DEG$gene<- rownames(MAC_Mono_pediatric_Thymusvsadult_DEG)
EC_pediatric_Thymusvsadult_DEG$gene<- rownames(EC_pediatric_Thymusvsadult_DEG)
FB_pediatric_Thymusvsadult_DEG$gene<- rownames(FB_pediatric_Thymusvsadult_DEG)
VSMC_Peri_pediatric_Thymusvsadult_DEG$gene<- rownames(VSMC_Peri_pediatric_Thymusvsadult_DEG)


#2 Create a list of dataframes with differentially expressed genes that you want to combine
DEG_strip_list <- list(TC_pediatric_Thymusvsadult_DEG, BC_pediatric_Thymusvsadult_DEG, PC_pediatric_Thymusvsadult_DEG, TEC_pediatric_Thymusvsadult_DEG, DC_pediatric_Thymusvsadult_DEG, MAC_Mono_pediatric_Thymusvsadult_DEG, EC_pediatric_Thymusvsadult_DEG, FB_pediatric_Thymusvsadult_DEG, VSMC_Peri_pediatric_Thymusvsadult_DEG)



# combine the dataframes into a single dataframe with an added column indicating the cell population
DEG_strip_df <- bind_rows(DEG_strip_list, .id = "celltype.tissue")


DEG_strip_df <- DEG_strip_df %>% mutate(Threshold_met =
                                          ifelse(abs(avg_log2FC) > 0.58 & p_val_adj < 0.05, 'p-value_adj & log2FC',
                                                 ifelse(abs(avg_log2FC) < 0.58 & p_val_adj < 0.05, 'p-value_adj',
                                                        ifelse(abs(avg_log2FC) > 0.58 & p_val_adj > 0.05, 'log2FC', 'Not sig.'))))




celltype_order <- c("TC","BC","PC","TEC","DC", "MAC/Mono","EC","FB","VSMC/Peri")



pdf("Graphs/Main/DEG_pediatric_Thymus_vs_adult_labeled.pdf", height= 5, width= 25)
ggplot(DEG_strip_df, aes(x = celltype, y = avg_log2FC, color = Threshold_met)) +
  geom_jitter(width = 0.25) +
  scale_color_manual(values = c('p-value_adj & log2FC'="#9c2162",
                                'p-value_adj' ="#0b9ed2",
                                'log2FC' ="#b2bd5b",
                                'Not sig.'="#A3A3A3"))  +
  theme_classic() +
  geom_text_repel(aes(label = ifelse(p_val_adj < 0.05 & abs(avg_log2FC) >= 8, gene, "")),
                  size = 4, max.overlaps = 8, fontface = "bold", point.padding = 0, box.padding = 0.2, min.segment.length = Inf,
                  segment.color = "grey", color = "black")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size=14),
        axis.text.y = element_text(face = "bold", size = 16),
        axis.line = element_line(size = 1.0),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold", margin = margin(r = 45)),
        plot.title = element_text(size = 16, face = "bold"),
        legend.text = element_text(size = 12, face = "bold"),
        legend.title = element_text(size = 14, face = "bold")) +
  ggtitle("DEG: pediatric_Thymus vs. Thym") +
  theme(plot.title = element_text(size = 16, face = "bold")) +
  scale_x_discrete(limits = celltype_order)
dev.off()

######TTH## #
#Clustermarker juxtaposition celltype.tissue against adult counterpart###
Idents(th.total)<-"celltype.tissue"
table(th.total$celltype.tissue)

TC_TTHvsadult_DEG       <- FindMarkers(th.total, ident.1= "TC_TTH", ident.2="TC_adult_Thymus",assay = "RNA", min.pct = 0.1)
BC_TTHvsadult_DEG       <- FindMarkers(th.total, ident.1= "BC_TTH", ident.2="BC_adult_Thymus",assay = "RNA", min.pct = 0.1)
PC_TTHvsadult_DEG       <- FindMarkers(th.total, ident.1= "PC_TTH", ident.2="PC_adult_Thymus",assay = "RNA", min.pct = 0.1)
TEC_TTHvsadult_DEG       <- FindMarkers(th.total, ident.1= "TEC_TTH", ident.2="TEC_adult_Thymus",assay = "RNA", min.pct = 0.1)
DC_TTHvsadult_DEG       <- FindMarkers(th.total, ident.1= "DC_TTH", ident.2="DC_adult_Thymus",assay = "RNA", min.pct = 0.1)
MAC_Mono_TTHvsadult_DEG       <- FindMarkers(th.total, ident.1= "MAC/Mono_TTH", ident.2="MAC/Mono_adult_Thymus",assay = "RNA", min.pct = 0.1)
EC_TTHvsadult_DEG       <- FindMarkers(th.total, ident.1= "EC_TTH", ident.2="EC_adult_Thymus",assay = "RNA", min.pct = 0.1)
FB_TTHvsadult_DEG       <- FindMarkers(th.total, ident.1= "FB_TTH", ident.2="FB_adult_Thymus",assay = "RNA", min.pct = 0.1)
VSMC_Peri_TTHvsadult_DEG       <- FindMarkers(th.total, ident.1= "VSMC/Peri_TTH", ident.2="VSMC/Peri_adult_Thymus",assay = "RNA", min.pct = 0.1)


TC_TTHvsadult_DEG$celltype <- "TC"
BC_TTHvsadult_DEG$celltype <- "BC"
PC_TTHvsadult_DEG$celltype <- "PC"
TEC_TTHvsadult_DEG$celltype <- "TEC"
DC_TTHvsadult_DEG$celltype <- "DC"
MAC_Mono_TTHvsadult_DEG$celltype <- "MAC/Mono"
EC_TTHvsadult_DEG$celltype <- "EC"
FB_TTHvsadult_DEG$celltype <- "FB"
VSMC_Peri_TTHvsadult_DEG$celltype <- "VSMC/Peri"


write.xlsx(TC_TTHvsadult_DEG, "Lists/Stripplot_lists/TC_TTHvsadult_DEG.xlsx")
write.xlsx(BC_TTHvsadult_DEG, "Lists/Stripplot_lists/BC_TTHvsadult_DEG.xlsx")
write.xlsx(PC_TTHvsadult_DEG, "Lists/Stripplot_lists/PC_TTHvsadult_DEG.xlsx")
write.xlsx(TEC_TTHvsadult_DEG, "Lists/Stripplot_lists/TEC_TTHvsadult_DEG.xlsx")
write.xlsx(DC_TTHvsadult_DEG, "Lists/Stripplot_lists/DC_TTHvsadult_DEG.xlsx")
write.xlsx(MAC_Mono_TTHvsadult_DEG, "Lists/Stripplot_lists/MAC_Mono_TTHvsadult_DEG.xlsx")
write.xlsx(EC_TTHvsadult_DEG, "Lists/Stripplot_lists/EC_TTHvsadult_DEG.xlsx")
write.xlsx(FB_TTHvsadult_DEG, "Lists/Stripplot_lists/FB_TTHvsadult_DEG.xlsx")
write.xlsx(VSMC_Peri_TTHvsadult_DEG, "Lists/Stripplot_lists/VSMC_Peri_TTHvsadult_DEG.xlsx")


###### Strip plot for DEG###

#1 Give each dataframe an additional column with the gene names the rownames might change later this will help you keep track

TC_TTHvsadult_DEG$gene<- rownames(TC_TTHvsadult_DEG)
BC_TTHvsadult_DEG$gene<- rownames(BC_TTHvsadult_DEG)
PC_TTHvsadult_DEG$gene<- rownames(PC_TTHvsadult_DEG)
TEC_TTHvsadult_DEG$gene<- rownames(TEC_TTHvsadult_DEG)
DC_TTHvsadult_DEG$gene<- rownames(DC_TTHvsadult_DEG)
MAC_Mono_TTHvsadult_DEG$gene<- rownames(MAC_Mono_TTHvsadult_DEG)
EC_TTHvsadult_DEG$gene<- rownames(EC_TTHvsadult_DEG)
FB_TTHvsadult_DEG$gene<- rownames(FB_TTHvsadult_DEG)
VSMC_Peri_TTHvsadult_DEG$gene<- rownames(VSMC_Peri_TTHvsadult_DEG)


#2 Create a list of dataframes with differentially expressed genes that you want to combine
DEG_strip_list <- list(TC_TTHvsadult_DEG, BC_TTHvsadult_DEG, PC_TTHvsadult_DEG, TEC_TTHvsadult_DEG, DC_TTHvsadult_DEG, MAC_Mono_TTHvsadult_DEG, EC_TTHvsadult_DEG, FB_TTHvsadult_DEG, VSMC_Peri_TTHvsadult_DEG)



# combine the dataframes into a single dataframe with an added column indicating the cell population
DEG_strip_df <- bind_rows(DEG_strip_list, .id = "celltype.tissue")


DEG_strip_df <- DEG_strip_df %>% mutate(Threshold_met =
                                          ifelse(abs(avg_log2FC) > 0.58 & p_val_adj < 0.05, 'p-value_adj & log2FC',
                                                 ifelse(abs(avg_log2FC) < 0.58 & p_val_adj < 0.05, 'p-value_adj',
                                                        ifelse(abs(avg_log2FC) > 0.58 & p_val_adj > 0.05, 'log2FC', 'Not sig.'))))




celltype_order <- c("TC","BC","PC","TEC","DC", "MAC/Mono","EC","FB","VSMC/Peri")



pdf("Graphs/Main/DEG_TTH_vs_adult_labeled.pdf", height= 5, width= 25)
ggplot(DEG_strip_df, aes(x = celltype, y = avg_log2FC, color = Threshold_met)) +
  geom_jitter(width = 0.25) +
  scale_color_manual(values = c('p-value_adj & log2FC'="#9c2162",
                                'p-value_adj' ="#0b9ed2",
                                'log2FC' ="#b2bd5b",
                                'Not sig.'="#A3A3A3"))  +
  theme_classic() +
  geom_text_repel(aes(label = ifelse(p_val_adj < 0.05 & abs(avg_log2FC) >= 8, gene, "")),
                  size = 4, max.overlaps = 8, fontface = "bold", point.padding = 0, box.padding = 0.2, min.segment.length = Inf,
                  segment.color = "grey", color = "black")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size=14),
        axis.text.y = element_text(face = "bold", size = 16),
        axis.line = element_line(size = 1.0),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold", margin = margin(r = 45)),
        plot.title = element_text(size = 16, face = "bold"),
        legend.text = element_text(size = 12, face = "bold"),
        legend.title = element_text(size = 14, face = "bold")) +
  ggtitle("DEG: TTH vs. Thym") +
  theme(plot.title = element_text(size = 16, face = "bold")) +
  scale_x_discrete(limits = celltype_order)
dev.off()


######TET_A## #
#Clustermarker juxtaposition celltype.tissue against adult counterpart###
Idents(th.total)<-"celltype.tissue"
TC_TET_Avsadult_DEG       <- FindMarkers(th.total, ident.1= "TC_TET_A", ident.2="TC_adult_Thymus",assay = "RNA", min.pct = 0.1)
BC_TET_Avsadult_DEG       <- FindMarkers(th.total, ident.1= "BC_TET_A", ident.2="BC_adult_Thymus",assay = "RNA", min.pct = 0.1)
PC_TET_Avsadult_DEG       <- FindMarkers(th.total, ident.1= "PC_TET_A", ident.2="PC_adult_Thymus",assay = "RNA", min.pct = 0.1)
TEC_TET_Avsadult_DEG       <- FindMarkers(th.total, ident.1= "TEC_TET_A", ident.2="TEC_adult_Thymus",assay = "RNA", min.pct = 0.1)
DC_TET_Avsadult_DEG       <- FindMarkers(th.total, ident.1= "DC_TET_A", ident.2="DC_adult_Thymus",assay = "RNA", min.pct = 0.1)
MAC_Mono_TET_Avsadult_DEG       <- FindMarkers(th.total, ident.1= "MAC/Mono_TET_A", ident.2="MAC/Mono_adult_Thymus",assay = "RNA", min.pct = 0.1)
EC_TET_Avsadult_DEG       <- FindMarkers(th.total, ident.1= "EC_TET_A", ident.2="EC_adult_Thymus",assay = "RNA", min.pct = 0.1)
FB_TET_Avsadult_DEG       <- FindMarkers(th.total, ident.1= "FB_TET_A", ident.2="FB_adult_Thymus",assay = "RNA", min.pct = 0.1)
VSMC_Peri_TET_Avsadult_DEG       <- FindMarkers(th.total, ident.1= "VSMC/Peri_TET_A", ident.2="VSMC/Peri_adult_Thymus",assay = "RNA", min.pct = 0.1)


TC_TET_Avsadult_DEG$celltype <- "TC"
BC_TET_Avsadult_DEG$celltype <- "BC"
PC_TET_Avsadult_DEG$celltype <- "PC"
TEC_TET_Avsadult_DEG$celltype <- "TEC"
DC_TET_Avsadult_DEG$celltype <- "DC"
MAC_Mono_TET_Avsadult_DEG$celltype <- "MAC/Mono"
EC_TET_Avsadult_DEG$celltype <- "EC"
FB_TET_Avsadult_DEG$celltype <- "FB"
VSMC_Peri_TET_Avsadult_DEG$celltype <- "VSMC/Peri"

write.xlsx(TC_TET_Avsadult_DEG, "Lists/Stripplot_lists/TC_TET_Avsadult_DEG.xlsx")
write.xlsx(BC_TET_Avsadult_DEG, "Lists/Stripplot_lists/BC_TET_Avsadult_DEG.xlsx")
write.xlsx(PC_TET_Avsadult_DEG, "Lists/Stripplot_lists/PC_TET_Avsadult_DEG.xlsx")
write.xlsx(TEC_TET_Avsadult_DEG, "Lists/Stripplot_lists/TEC_TET_Avsadult_DEG.xlsx")
write.xlsx(DC_TET_Avsadult_DEG, "Lists/Stripplot_lists/DC_TET_Avsadult_DEG.xlsx")
write.xlsx(MAC_Mono_TET_Avsadult_DEG, "Lists/Stripplot_lists/MAC_Mono_TET_Avsadult_DEG.xlsx")
write.xlsx(EC_TET_Avsadult_DEG, "Lists/Stripplot_lists/EC_TET_Avsadult_DEG.xlsx")
write.xlsx(FB_TET_Avsadult_DEG, "Lists/Stripplot_lists/FB_TET_Avsadult_DEG.xlsx")
write.xlsx(VSMC_Peri_TET_Avsadult_DEG, "Lists/Stripplot_lists/VSMC_Peri_TET_Avsadult_DEG.xlsx")


###### Strip plot for DEG###

#1 Give each dataframe an additional column with the gene names the rownames might change later this will help you keep track

TC_TET_Avsadult_DEG$gene<- rownames(TC_TET_Avsadult_DEG)
BC_TET_Avsadult_DEG$gene<- rownames(BC_TET_Avsadult_DEG)
PC_TET_Avsadult_DEG$gene<- rownames(PC_TET_Avsadult_DEG)
TEC_TET_Avsadult_DEG$gene<- rownames(TEC_TET_Avsadult_DEG)
DC_TET_Avsadult_DEG$gene<- rownames(DC_TET_Avsadult_DEG)
MAC_Mono_TET_Avsadult_DEG$gene<- rownames(MAC_Mono_TET_Avsadult_DEG)
EC_TET_Avsadult_DEG$gene<- rownames(EC_TET_Avsadult_DEG)
FB_TET_Avsadult_DEG$gene<- rownames(FB_TET_Avsadult_DEG)
VSMC_Peri_TET_Avsadult_DEG$gene<- rownames(VSMC_Peri_TET_Avsadult_DEG)


#2 Create a list of dataframes with differentially expressed genes that you want to combine
DEG_strip_list <- list(TC_TET_Avsadult_DEG, BC_TET_Avsadult_DEG, PC_TET_Avsadult_DEG, TEC_TET_Avsadult_DEG, DC_TET_Avsadult_DEG, MAC_Mono_TET_Avsadult_DEG, EC_TET_Avsadult_DEG, FB_TET_Avsadult_DEG, VSMC_Peri_TET_Avsadult_DEG)



# combine the dataframes into a single dataframe with an added column indicating the cell population
DEG_strip_df <- bind_rows(DEG_strip_list, .id = "celltype.tissue")


DEG_strip_df <- DEG_strip_df %>% mutate(Threshold_met =
                                          ifelse(abs(avg_log2FC) > 0.58 & p_val_adj < 0.05, 'p-value_adj & log2FC',
                                                 ifelse(abs(avg_log2FC) < 0.58 & p_val_adj < 0.05, 'p-value_adj',
                                                        ifelse(abs(avg_log2FC) > 0.58 & p_val_adj > 0.05, 'log2FC', 'Not sig.'))))




celltype_order <- c("TC","BC","PC","TEC","DC", "MAC/Mono","EC","FB","VSMC/Peri")

pdf("Graphs/Main/DEG_TET_A_vs_adult_labeled.pdf", height= 5, width= 25)
ggplot(DEG_strip_df, aes(x = celltype, y = avg_log2FC, color = Threshold_met)) +
  geom_jitter(width = 0.25) +
  scale_color_manual(values = c('p-value_adj & log2FC'="#9c2162",
                                'p-value_adj' ="#0b9ed2",
                                'log2FC' ="#b2bd5b",
                                'Not sig.'="#A3A3A3"))  +
  labs(x = "Cluster", y = "Average log2 fold change", color = "Threshold_met") +
  theme_classic() +
  geom_text_repel(aes(label = ifelse(p_val_adj < 0.05 & abs(avg_log2FC) >= 8, gene, "")),
                  size = 4, max.overlaps = 8, fontface = "bold", point.padding = 0, box.padding = 0.2, min.segment.length = Inf,
                  segment.color = "grey", color = "black")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size=14),
        axis.text.y = element_text(face = "bold", size = 16),
        axis.line = element_line(size = 1.0),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold", margin = margin(r = 45)),
        plot.title = element_text(size = 16, face = "bold"),
        legend.text = element_text(size = 12, face = "bold"),
        legend.title = element_text(size = 14, face = "bold")) +
  ggtitle("DEG: TET-A vs. Thym") +
  theme(plot.title = element_text(size = 16, face = "bold")) +
  scale_x_discrete(limits = celltype_order)
dev.off()

######TET_AB## #
#Clustermarker juxtaposition celltype.tissue against adult counterpart###
Idents(th.total)<-"celltype.tissue"
TC_TET_ABvsadult_DEG       <- FindMarkers(th.total, ident.1= "TC_TET_AB", ident.2="TC_adult_Thymus",assay = "RNA", min.pct = 0.1)
BC_TET_ABvsadult_DEG       <- FindMarkers(th.total, ident.1= "BC_TET_AB", ident.2="BC_adult_Thymus",assay = "RNA", min.pct = 0.1)
PC_TET_ABvsadult_DEG       <- FindMarkers(th.total, ident.1= "PC_TET_AB", ident.2="PC_adult_Thymus",assay = "RNA", min.pct = 0.1)
TEC_TET_ABvsadult_DEG       <- FindMarkers(th.total, ident.1= "TEC_TET_AB", ident.2="TEC_adult_Thymus",assay = "RNA", min.pct = 0.1)
DC_TET_ABvsadult_DEG       <- FindMarkers(th.total, ident.1= "DC_TET_AB", ident.2="DC_adult_Thymus",assay = "RNA", min.pct = 0.1)
MAC_Mono_TET_ABvsadult_DEG       <- FindMarkers(th.total, ident.1= "MAC/Mono_TET_AB", ident.2="MAC/Mono_adult_Thymus",assay = "RNA", min.pct = 0.1)
EC_TET_ABvsadult_DEG       <- FindMarkers(th.total, ident.1= "EC_TET_AB", ident.2="EC_adult_Thymus",assay = "RNA", min.pct = 0.1)
FB_TET_ABvsadult_DEG       <- FindMarkers(th.total, ident.1= "FB_TET_AB", ident.2="FB_adult_Thymus",assay = "RNA", min.pct = 0.1)
VSMC_Peri_TET_ABvsadult_DEG       <- FindMarkers(th.total, ident.1= "VSMC/Peri_TET_AB", ident.2="VSMC/Peri_adult_Thymus",assay = "RNA", min.pct = 0.1)

TC_TET_ABvsadult_DEG$celltype <- "TC"
BC_TET_ABvsadult_DEG$celltype <- "BC"
PC_TET_ABvsadult_DEG$celltype <- "PC"
TEC_TET_ABvsadult_DEG$celltype <- "TEC"
DC_TET_ABvsadult_DEG$celltype <- "DC"
MAC_Mono_TET_ABvsadult_DEG$celltype <- "MAC/Mono"
EC_TET_ABvsadult_DEG$celltype <- "EC"
FB_TET_ABvsadult_DEG$celltype <- "FB"
VSMC_Peri_TET_ABvsadult_DEG$celltype <- "VSMC/Peri"

write.xlsx(TC_TET_ABvsadult_DEG, "Lists/Stripplot_lists/TC_TET_ABvsadult_DEG.xlsx")
write.xlsx(BC_TET_ABvsadult_DEG, "Lists/Stripplot_lists/BC_TET_ABvsadult_DEG.xlsx")
write.xlsx(PC_TET_ABvsadult_DEG, "Lists/Stripplot_lists/PC_TET_ABvsadult_DEG.xlsx")
write.xlsx(TEC_TET_ABvsadult_DEG, "Lists/Stripplot_lists/TEC_TET_ABvsadult_DEG.xlsx")
write.xlsx(DC_TET_ABvsadult_DEG, "Lists/Stripplot_lists/DC_TET_ABvsadult_DEG.xlsx")
write.xlsx(MAC_Mono_TET_ABvsadult_DEG, "Lists/Stripplot_lists/MAC_Mono_TET_ABvsadult_DEG.xlsx")
write.xlsx(EC_TET_ABvsadult_DEG, "Lists/Stripplot_lists/EC_TET_ABvsadult_DEG.xlsx")
write.xlsx(FB_TET_ABvsadult_DEG, "Lists/Stripplot_lists/FB_TET_ABvsadult_DEG.xlsx")
write.xlsx(VSMC_Peri_TET_ABvsadult_DEG, "Lists/Stripplot_lists/VSMC_Peri_TET_ABvsadult_DEG.xlsx")


###### Strip plot for DEG###

#1 Give each dataframe an additional column with the gene names the rownames might change later this will help you keep track

TC_TET_ABvsadult_DEG$gene<- rownames(TC_TET_ABvsadult_DEG)
BC_TET_ABvsadult_DEG$gene<- rownames(BC_TET_ABvsadult_DEG)
PC_TET_ABvsadult_DEG$gene<- rownames(PC_TET_ABvsadult_DEG)
TEC_TET_ABvsadult_DEG$gene<- rownames(TEC_TET_ABvsadult_DEG)
DC_TET_ABvsadult_DEG$gene<- rownames(DC_TET_ABvsadult_DEG)
MAC_Mono_TET_ABvsadult_DEG$gene<- rownames(MAC_Mono_TET_ABvsadult_DEG)
EC_TET_ABvsadult_DEG$gene<- rownames(EC_TET_ABvsadult_DEG)
FB_TET_ABvsadult_DEG$gene<- rownames(FB_TET_ABvsadult_DEG)
VSMC_Peri_TET_ABvsadult_DEG$gene<- rownames(VSMC_Peri_TET_ABvsadult_DEG)


#2 Create a list of dataframes with differentially expressed genes that you want to combine
DEG_strip_list <- list(TC_TET_ABvsadult_DEG, BC_TET_ABvsadult_DEG, PC_TET_ABvsadult_DEG, TEC_TET_ABvsadult_DEG, DC_TET_ABvsadult_DEG, MAC_Mono_TET_ABvsadult_DEG, EC_TET_ABvsadult_DEG, FB_TET_ABvsadult_DEG, VSMC_Peri_TET_ABvsadult_DEG)



# combine the dataframes into a single dataframe with an added column indicating the cell population
DEG_strip_df <- bind_rows(DEG_strip_list, .id = "celltype.tissue")


DEG_strip_df <- DEG_strip_df %>% mutate(Threshold_met =
                                          ifelse(abs(avg_log2FC) > 0.58 & p_val_adj < 0.05, 'p-value_adj & log2FC',
                                                 ifelse(abs(avg_log2FC) < 0.58 & p_val_adj < 0.05, 'p-value_adj',
                                                        ifelse(abs(avg_log2FC) > 0.58 & p_val_adj > 0.05, 'log2FC', 'Not sig.'))))




celltype_order <- c("TC","BC","PC","TEC","DC", "MAC/Mono","EC","FB","VSMC/Peri")



pdf("Graphs/Main/DEG_TET_AB_vs_adult_labeled.pdf", height= 5, width= 25)
ggplot(DEG_strip_df, aes(x = celltype, y = avg_log2FC, color = Threshold_met)) +
  geom_jitter(width = 0.25) +
  scale_color_manual(values = c('p-value_adj & log2FC'="#9c2162",
                                'p-value_adj' ="#0b9ed2",
                                'log2FC' ="#b2bd5b",
                                'Not sig.'="#A3A3A3"))  +
  theme_classic() +
  geom_text_repel(aes(label = ifelse(p_val_adj < 0.05 & abs(avg_log2FC) >= 8, gene, "")),
                  size = 4, max.overlaps = 8, fontface = "bold", point.padding = 0, box.padding = 0.2, min.segment.length = Inf,
                  segment.color = "grey", color = "black") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size=14),
        axis.text.y = element_text(face = "bold", size = 16),
        axis.line = element_line(size = 1.0),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold", margin = margin(r = 45)),
        plot.title = element_text(size = 16, face = "bold"),
        legend.text = element_text(size = 12, face = "bold"),
        legend.title = element_text(size = 14, face = "bold")) +
  ggtitle("DEG: TET-AB vs. Thym") +
  theme(plot.title = element_text(size = 16, face = "bold")) +
  scale_x_discrete(limits = celltype_order)
dev.off()

######TET_B## #
#Clustermarker juxtaposition celltype.tissue against adult counterpart###
Idents(th.total)<-"celltype.tissue"
TC_TET_Bvsadult_DEG       <- FindMarkers(th.total, ident.1= "TC_TET_B", ident.2="TC_adult_Thymus",assay = "RNA", min.pct = 0.1)
BC_TET_Bvsadult_DEG       <- FindMarkers(th.total, ident.1= "BC_TET_B", ident.2="BC_adult_Thymus",assay = "RNA", min.pct = 0.1)
PC_TET_Bvsadult_DEG       <- FindMarkers(th.total, ident.1= "PC_TET_B", ident.2="PC_adult_Thymus",assay = "RNA", min.pct = 0.1)
TEC_TET_Bvsadult_DEG       <- FindMarkers(th.total, ident.1= "TEC_TET_B", ident.2="TEC_adult_Thymus",assay = "RNA", min.pct = 0.1)
DC_TET_Bvsadult_DEG       <- FindMarkers(th.total, ident.1= "DC_TET_B", ident.2="DC_adult_Thymus",assay = "RNA", min.pct = 0.1)
MAC_Mono_TET_Bvsadult_DEG       <- FindMarkers(th.total, ident.1= "MAC/Mono_TET_B", ident.2="MAC/Mono_adult_Thymus",assay = "RNA", min.pct = 0.1)
EC_TET_Bvsadult_DEG       <- FindMarkers(th.total, ident.1= "EC_TET_B", ident.2="EC_adult_Thymus",assay = "RNA", min.pct = 0.1)
FB_TET_Bvsadult_DEG       <- FindMarkers(th.total, ident.1= "FB_TET_B", ident.2="FB_adult_Thymus",assay = "RNA", min.pct = 0.1)
VSMC_Peri_TET_Bvsadult_DEG       <- FindMarkers(th.total, ident.1= "VSMC/Peri_TET_B", ident.2="VSMC/Peri_adult_Thymus",assay = "RNA", min.pct = 0.1)

TC_TET_Bvsadult_DEG$celltype <- "TC"
BC_TET_Bvsadult_DEG$celltype <- "BC"
PC_TET_Bvsadult_DEG$celltype <- "PC"
TEC_TET_Bvsadult_DEG$celltype <- "TEC"
DC_TET_Bvsadult_DEG$celltype <- "DC"
MAC_Mono_TET_Bvsadult_DEG$celltype <- "MAC/Mono"
EC_TET_Bvsadult_DEG$celltype <- "EC"
FB_TET_Bvsadult_DEG$celltype <- "FB"
VSMC_Peri_TET_Bvsadult_DEG$celltype <- "VSMC/Peri"

write.xlsx(TC_TET_Bvsadult_DEG, "Lists/Stripplot_lists/TC_TET_Bvsadult_DEG.xlsx")
write.xlsx(BC_TET_Bvsadult_DEG, "Lists/Stripplot_lists/BC_TET_Bvsadult_DEG.xlsx")
write.xlsx(PC_TET_Bvsadult_DEG, "Lists/Stripplot_lists/PC_TET_Bvsadult_DEG.xlsx")
write.xlsx(TEC_TET_Bvsadult_DEG, "Lists/Stripplot_lists/TEC_TET_Bvsadult_DEG.xlsx")
write.xlsx(DC_TET_Bvsadult_DEG, "Lists/Stripplot_lists/DC_TET_Bvsadult_DEG.xlsx")
write.xlsx(MAC_Mono_TET_Bvsadult_DEG, "Lists/Stripplot_lists/MAC_Mono_TET_Bvsadult_DEG.xlsx")
write.xlsx(EC_TET_Bvsadult_DEG, "Lists/Stripplot_lists/EC_TET_Bvsadult_DEG.xlsx")
write.xlsx(FB_TET_Bvsadult_DEG, "Lists/Stripplot_lists/FB_TET_Bvsadult_DEG.xlsx")
write.xlsx(VSMC_Peri_TET_Bvsadult_DEG, "Lists/Stripplot_lists/VSMC_Peri_TET_Bvsadult_DEG.xlsx")


###### Strip plot for DEG###

#1 Give each dataframe an additional column with the gene names the rownames might change later this will help you keep track

TC_TET_Bvsadult_DEG$gene<- rownames(TC_TET_Bvsadult_DEG)
BC_TET_Bvsadult_DEG$gene<- rownames(BC_TET_Bvsadult_DEG)
PC_TET_Bvsadult_DEG$gene<- rownames(PC_TET_Bvsadult_DEG)
TEC_TET_Bvsadult_DEG$gene<- rownames(TEC_TET_Bvsadult_DEG)
DC_TET_Bvsadult_DEG$gene<- rownames(DC_TET_Bvsadult_DEG)
MAC_Mono_TET_Bvsadult_DEG$gene<- rownames(MAC_Mono_TET_Bvsadult_DEG)
EC_TET_Bvsadult_DEG$gene<- rownames(EC_TET_Bvsadult_DEG)
FB_TET_Bvsadult_DEG$gene<- rownames(FB_TET_Bvsadult_DEG)
VSMC_Peri_TET_Bvsadult_DEG$gene<- rownames(VSMC_Peri_TET_Bvsadult_DEG)


#2 Create a list of dataframes with differentially expressed genes that you want to combine
DEG_strip_list <- list(TC_TET_Bvsadult_DEG, BC_TET_Bvsadult_DEG, PC_TET_Bvsadult_DEG, TEC_TET_Bvsadult_DEG, DC_TET_Bvsadult_DEG, MAC_Mono_TET_Bvsadult_DEG, EC_TET_Bvsadult_DEG, FB_TET_Bvsadult_DEG, VSMC_Peri_TET_Bvsadult_DEG)



# combine the dataframes into a single dataframe with an added column indicating the cell population
DEG_strip_df <- bind_rows(DEG_strip_list, .id = "celltype.tissue")


DEG_strip_df <- DEG_strip_df %>% mutate(Threshold_met =
                                          ifelse(abs(avg_log2FC) > 0.58 & p_val_adj < 0.05, 'p-value_adj & log2FC',
                                                 ifelse(abs(avg_log2FC) < 0.58 & p_val_adj < 0.05, 'p-value_adj',
                                                        ifelse(abs(avg_log2FC) > 0.58 & p_val_adj > 0.05, 'log2FC', 'Not sig.'))))




celltype_order <- c("TC","BC","PC","TEC","DC", "MAC/Mono","EC","FB","VSMC/Peri")



pdf("Graphs/Main/DEG_TET_B_vs_adult_labeled.pdf", height= 5, width= 25)
ggplot(DEG_strip_df, aes(x = celltype, y = avg_log2FC, color = Threshold_met)) +
  geom_jitter(width = 0.25) +
  scale_color_manual(values = c('p-value_adj & log2FC'="#9c2162",
                                'p-value_adj' ="#0b9ed2",
                                'log2FC' ="#b2bd5b",
                                'Not sig.'="#A3A3A3")) +
  labs(x = "Cluster", y = "Average log2 fold change", color = "Threshold_met") +
  theme_classic() +
  geom_text_repel(aes(label = ifelse(p_val_adj < 0.05 & abs(avg_log2FC) >= 8, gene, "")),
                  size = 4, max.overlaps = 8, fontface = "bold", point.padding = 0, box.padding = 0.2, min.segment.length = Inf,
                  segment.color = "grey", color = "black")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size=14),
        axis.text.y = element_text(face = "bold", size = 16),
        axis.line = element_line(size = 1.0),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold", margin = margin(r = 45)),
        plot.title = element_text(size = 16, face = "bold"),
        legend.text = element_text(size = 12, face = "bold"),
        legend.title = element_text(size = 14, face = "bold")) +
  ggtitle("DEG: TET-B vs. Thym") +
  theme(plot.title = element_text(size = 16, face = "bold")) +
  scale_x_discrete(limits = celltype_order)
dev.off()

######TET_C## #
#Clustermarker juxtaposition celltype.tissue against adult counterpart###
Idents(th.total)<-"celltype.tissue"
TC_TET_Cvsadult_DEG       <- FindMarkers(th.total, ident.1= "TC_TET_C", ident.2="TC_adult_Thymus",assay = "RNA", min.pct = 0.1)
BC_TET_Cvsadult_DEG       <- FindMarkers(th.total, ident.1= "BC_TET_C", ident.2="BC_adult_Thymus",assay = "RNA", min.pct = 0.1)
PC_TET_Cvsadult_DEG       <- FindMarkers(th.total, ident.1= "PC_TET_C", ident.2="PC_adult_Thymus",assay = "RNA", min.pct = 0.1)
TEC_TET_Cvsadult_DEG       <- FindMarkers(th.total, ident.1= "TEC_TET_C", ident.2="TEC_adult_Thymus",assay = "RNA", min.pct = 0.1)
DC_TET_Cvsadult_DEG       <- FindMarkers(th.total, ident.1= "DC_TET_C", ident.2="DC_adult_Thymus",assay = "RNA", min.pct = 0.1)
MAC_Mono_TET_Cvsadult_DEG       <- FindMarkers(th.total, ident.1= "MAC/Mono_TET_C", ident.2="MAC/Mono_adult_Thymus",assay = "RNA", min.pct = 0.1)
EC_TET_Cvsadult_DEG       <- FindMarkers(th.total, ident.1= "EC_TET_C", ident.2="EC_adult_Thymus",assay = "RNA", min.pct = 0.1)
FB_TET_Cvsadult_DEG       <- FindMarkers(th.total, ident.1= "FB_TET_C", ident.2="FB_adult_Thymus",assay = "RNA", min.pct = 0.1)
VSMC_Peri_TET_Cvsadult_DEG       <- FindMarkers(th.total, ident.1= "VSMC/Peri_TET_C", ident.2="VSMC/Peri_adult_Thymus",assay = "RNA", min.pct = 0.1)

TC_TET_Cvsadult_DEG$celltype <- "TC"
BC_TET_Cvsadult_DEG$celltype <- "BC"
PC_TET_Cvsadult_DEG$celltype <- "PC"
TEC_TET_Cvsadult_DEG$celltype <- "TEC"
DC_TET_Cvsadult_DEG$celltype <- "DC"
MAC_Mono_TET_Cvsadult_DEG$celltype <- "MAC/Mono"
EC_TET_Cvsadult_DEG$celltype <- "EC"
FB_TET_Cvsadult_DEG$celltype <- "FB"
VSMC_Peri_TET_Cvsadult_DEG$celltype <- "VSMC/Peri"

write.xlsx(TC_TET_Cvsadult_DEG, "Lists/Stripplot_lists/TC_TET_Cvsadult_DEG.xlsx")
write.xlsx(BC_TET_Cvsadult_DEG, "Lists/Stripplot_lists/BC_TET_Cvsadult_DEG.xlsx")
write.xlsx(PC_TET_Cvsadult_DEG, "Lists/Stripplot_lists/PC_TET_Cvsadult_DEG.xlsx")
write.xlsx(TEC_TET_Cvsadult_DEG, "Lists/Stripplot_lists/TEC_TET_Cvsadult_DEG.xlsx")
write.xlsx(DC_TET_Cvsadult_DEG, "Lists/Stripplot_lists/DC_TET_Cvsadult_DEG.xlsx")
write.xlsx(MAC_Mono_TET_Cvsadult_DEG, "Lists/Stripplot_lists/MAC_Mono_TET_Cvsadult_DEG.xlsx")
write.xlsx(EC_TET_Cvsadult_DEG, "Lists/Stripplot_lists/EC_TET_Cvsadult_DEG.xlsx")
write.xlsx(FB_TET_Cvsadult_DEG, "Lists/Stripplot_lists/FB_TET_Cvsadult_DEG.xlsx")
write.xlsx(VSMC_Peri_TET_Cvsadult_DEG, "Lists/Stripplot_lists/VSMC_Peri_TET_Cvsadult_DEG.xlsx")


###### Strip plot for DEG###

#1 Give each dataframe an additional column with the gene names the rownames might change later this will help you keep track

TC_TET_Cvsadult_DEG$gene<- rownames(TC_TET_Cvsadult_DEG)
BC_TET_Cvsadult_DEG$gene<- rownames(BC_TET_Cvsadult_DEG)
PC_TET_Cvsadult_DEG$gene<- rownames(PC_TET_Cvsadult_DEG)
TEC_TET_Cvsadult_DEG$gene<- rownames(TEC_TET_Cvsadult_DEG)
DC_TET_Cvsadult_DEG$gene<- rownames(DC_TET_Cvsadult_DEG)
MAC_Mono_TET_Cvsadult_DEG$gene<- rownames(MAC_Mono_TET_Cvsadult_DEG)
EC_TET_Cvsadult_DEG$gene<- rownames(EC_TET_Cvsadult_DEG)
FB_TET_Cvsadult_DEG$gene<- rownames(FB_TET_Cvsadult_DEG)
VSMC_Peri_TET_Cvsadult_DEG$gene<- rownames(VSMC_Peri_TET_Cvsadult_DEG)


#2 Create a list of dataframes with differentially expressed genes that you want to combine
DEG_strip_list <- list(TC_TET_Cvsadult_DEG, BC_TET_Cvsadult_DEG, PC_TET_Cvsadult_DEG, TEC_TET_Cvsadult_DEG, DC_TET_Cvsadult_DEG, MAC_Mono_TET_Cvsadult_DEG, EC_TET_Cvsadult_DEG, FB_TET_Cvsadult_DEG, VSMC_Peri_TET_Cvsadult_DEG)



# combine the dataframes into a single dataframe with an added column indicating the cell population
DEG_strip_df <- bind_rows(DEG_strip_list, .id = "celltype.tissue")


DEG_strip_df <- DEG_strip_df %>% mutate(Threshold_met =
                                          ifelse(abs(avg_log2FC) > 0.58 & p_val_adj < 0.05, 'p-value_adj & log2FC',
                                                 ifelse(abs(avg_log2FC) < 0.58 & p_val_adj < 0.05, 'p-value_adj',
                                                        ifelse(abs(avg_log2FC) > 0.58 & p_val_adj > 0.05, 'log2FC', 'Not sig.'))))




celltype_order <- c("TC","BC","PC","TEC","DC", "MAC/Mono","EC","FB","VSMC/Peri")



pdf("Graphs/Main/DEG_TET_C_vs_adult_labeled.pdf", height= 5, width= 25)
ggplot(DEG_strip_df, aes(x = celltype, y = avg_log2FC, color = Threshold_met)) +
  geom_jitter(width = 0.25) +
  scale_color_manual(values = c('p-value_adj & log2FC'="#9c2162",
                                'p-value_adj' ="#0b9ed2",
                                'log2FC' ="#b2bd5b",
                                'Not sig.'="#A3A3A3")) +
  labs(x = "Cluster", y = "Average log2 fold change", color = "Threshold_met") +
  theme_classic() +
  geom_text_repel(aes(label = ifelse(p_val_adj < 0.05 & abs(avg_log2FC) >= 8, gene, "")),
                  size = 4, max.overlaps = 8, fontface = "bold", point.padding = 0, box.padding = 0.2, min.segment.length = Inf,
                  segment.color = "grey", color = "black")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size=14),
        axis.text.y = element_text(face = "bold", size = 16),
        axis.line = element_line(size = 1.0),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold", margin = margin(r = 45)),
        plot.title = element_text(size = 16, face = "bold"),
        legend.text = element_text(size = 12, face = "bold"),
        legend.title = element_text(size = 14, face = "bold")) +
  ggtitle("DEG: TET-C vs. Thym") +
  theme(plot.title = element_text(size = 16, face = "bold")) +
  scale_x_discrete(limits = celltype_order)
dev.off()

######TET_MNT## #
#Clustermarker juxtaposition celltype.tissue against adult counterpart###
Idents(th.total)<-"celltype.tissue"
TC_TET_MNTvsadult_DEG       <- FindMarkers(th.total, ident.1= "TC_TET_MNT", ident.2="TC_adult_Thymus",assay = "RNA", min.pct = 0.1)
BC_TET_MNTvsadult_DEG       <- FindMarkers(th.total, ident.1= "BC_TET_MNT", ident.2="BC_adult_Thymus",assay = "RNA", min.pct = 0.1)
PC_TET_MNTvsadult_DEG       <- FindMarkers(th.total, ident.1= "PC_TET_MNT", ident.2="PC_adult_Thymus",assay = "RNA", min.pct = 0.1)
TEC_TET_MNTvsadult_DEG       <- FindMarkers(th.total, ident.1= "TEC_TET_MNT", ident.2="TEC_adult_Thymus",assay = "RNA", min.pct = 0.1)
DC_TET_MNTvsadult_DEG       <- FindMarkers(th.total, ident.1= "DC_TET_MNT", ident.2="DC_adult_Thymus",assay = "RNA", min.pct = 0.1)
MAC_Mono_TET_MNTvsadult_DEG       <- FindMarkers(th.total, ident.1= "MAC/Mono_TET_MNT", ident.2="MAC/Mono_adult_Thymus",assay = "RNA", min.pct = 0.1)
EC_TET_MNTvsadult_DEG       <- FindMarkers(th.total, ident.1= "EC_TET_MNT", ident.2="EC_adult_Thymus",assay = "RNA", min.pct = 0.1)
FB_TET_MNTvsadult_DEG       <- FindMarkers(th.total, ident.1= "FB_TET_MNT", ident.2="FB_adult_Thymus",assay = "RNA", min.pct = 0.1)
VSMC_Peri_TET_MNTvsadult_DEG       <- FindMarkers(th.total, ident.1= "VSMC/Peri_TET_MNT", ident.2="VSMC/Peri_adult_Thymus",assay = "RNA", min.pct = 0.1)

TC_TET_MNTvsadult_DEG$celltype <- "TC"
BC_TET_MNTvsadult_DEG$celltype <- "BC"
PC_TET_MNTvsadult_DEG$celltype <- "PC"
TEC_TET_MNTvsadult_DEG$celltype <- "TEC"
DC_TET_MNTvsadult_DEG$celltype <- "DC"
MAC_Mono_TET_MNTvsadult_DEG$celltype <- "MAC/Mono"
EC_TET_MNTvsadult_DEG$celltype <- "EC"
FB_TET_MNTvsadult_DEG$celltype <- "FB"
VSMC_Peri_TET_MNTvsadult_DEG$celltype <- "VSMC/Peri"

write.xlsx(TC_TET_MNTvsadult_DEG, "Lists/Stripplot_lists/TC_TET_MNTvsadult_DEG.xlsx")
write.xlsx(BC_TET_MNTvsadult_DEG, "Lists/Stripplot_lists/BC_TET_MNTvsadult_DEG.xlsx")
write.xlsx(PC_TET_MNTvsadult_DEG, "Lists/Stripplot_lists/PC_TET_MNTvsadult_DEG.xlsx")
write.xlsx(TEC_TET_MNTvsadult_DEG, "Lists/Stripplot_lists/TEC_TET_MNTvsadult_DEG.xlsx")
write.xlsx(DC_TET_MNTvsadult_DEG, "Lists/Stripplot_lists/DC_TET_MNTvsadult_DEG.xlsx")
write.xlsx(MAC_Mono_TET_MNTvsadult_DEG, "Lists/Stripplot_lists/MAC_Mono_TET_MNTvsadult_DEG.xlsx")
write.xlsx(EC_TET_MNTvsadult_DEG, "Lists/Stripplot_lists/EC_TET_MNTvsadult_DEG.xlsx")
write.xlsx(FB_TET_MNTvsadult_DEG, "Lists/Stripplot_lists/FB_TET_MNTvsadult_DEG.xlsx")
write.xlsx(VSMC_Peri_TET_MNTvsadult_DEG, "Lists/Stripplot_lists/VSMC_Peri_TET_MNTvsadult_DEG.xlsx")


###### Strip plot for DEG###

#1 Give each dataframe an additional column with the gene names the rownames might change later this will help you keep track

TC_TET_MNTvsadult_DEG$gene<- rownames(TC_TET_MNTvsadult_DEG)
BC_TET_MNTvsadult_DEG$gene<- rownames(BC_TET_MNTvsadult_DEG)
PC_TET_MNTvsadult_DEG$gene<- rownames(PC_TET_MNTvsadult_DEG)
TEC_TET_MNTvsadult_DEG$gene<- rownames(TEC_TET_MNTvsadult_DEG)
DC_TET_MNTvsadult_DEG$gene<- rownames(DC_TET_MNTvsadult_DEG)
MAC_Mono_TET_MNTvsadult_DEG$gene<- rownames(MAC_Mono_TET_MNTvsadult_DEG)
EC_TET_MNTvsadult_DEG$gene<- rownames(EC_TET_MNTvsadult_DEG)
FB_TET_MNTvsadult_DEG$gene<- rownames(FB_TET_MNTvsadult_DEG)
VSMC_Peri_TET_MNTvsadult_DEG$gene<- rownames(VSMC_Peri_TET_MNTvsadult_DEG)


#2 Create a list of dataframes with differentially expressed genes that you want to combine
DEG_strip_list <- list(TC_TET_MNTvsadult_DEG, BC_TET_MNTvsadult_DEG, PC_TET_MNTvsadult_DEG, TEC_TET_MNTvsadult_DEG, DC_TET_MNTvsadult_DEG, MAC_Mono_TET_MNTvsadult_DEG, EC_TET_MNTvsadult_DEG, FB_TET_MNTvsadult_DEG, VSMC_Peri_TET_MNTvsadult_DEG)



# combine the dataframes into a single dataframe with an added column indicating the cell population
DEG_strip_df <- bind_rows(DEG_strip_list, .id = "celltype.tissue")


DEG_strip_df <- DEG_strip_df %>% mutate(Threshold_met =
                                          ifelse(abs(avg_log2FC) > 0.58 & p_val_adj < 0.05, 'p-value_adj & log2FC',
                                                 ifelse(abs(avg_log2FC) < 0.58 & p_val_adj < 0.05, 'p-value_adj',
                                                        ifelse(abs(avg_log2FC) > 0.58 & p_val_adj > 0.05, 'log2FC', 'Not sig.'))))




celltype_order <- c("TC","BC","PC","TEC","DC", "MAC/Mono","EC","FB","VSMC/Peri")



pdf("Graphs/Main/DEG_TET_MNT_vs_adult_labeled.pdf", height= 5, width= 25)
ggplot(DEG_strip_df, aes(x = celltype, y = avg_log2FC, color = Threshold_met)) +
  geom_jitter(width = 0.25) +
  scale_color_manual(values = c('p-value_adj & log2FC'="#9c2162",
                                'p-value_adj' ="#0b9ed2",
                                'log2FC' ="#b2bd5b",
                                'Not sig.'="#A3A3A3"))  +
  labs(x = "Cluster", y = "Average log2 fold change", color = "Threshold_met") +
  theme_classic() +
  geom_text_repel(aes(label = ifelse(p_val_adj < 0.05 & abs(avg_log2FC) >= 8, gene, "")),
                  size = 4, max.overlaps = 8, fontface = "bold", point.padding = 0, box.padding = 0.2, min.segment.length = Inf,
                  segment.color = "grey", color = "black")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size=14),
        axis.text.y = element_text(face = "bold", size = 16),
        axis.line = element_line(size = 1.0),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold", margin = margin(r = 45)),
        plot.title = element_text(size = 16, face = "bold"),
        legend.text = element_text(size = 12, face = "bold"),
        legend.title = element_text(size = 14, face = "bold")) +
  ggtitle("DEG: TET-MNT vs. Thym") +
  theme(plot.title = element_text(size = 16, face = "bold")) +
  scale_x_discrete(limits = celltype_order)
dev.off()


#(44)_Bulkseq GTF2I_Pattern####

library(TCGAbiolinks)
library(DESeq2)
library(biomaRt)
library(pheatmap)
library(EnhancedVolcano)
library(viridis)
library(plyr)
library(vsn)
library(pals)
library(tidyverse)
library(dplyr)
library(scales)
library(grid)
library(xlsx)
library(readxl)



AddGeneSetScore <- function(
    dds,
    features,
    pool = NULL,
    nbin = 24,
    ctrl = 100,
    name = 'Set',
    seed = 123
) {
  if (!is.null(x = seed)) {
    set.seed(seed = seed)
  }
  
  `%||%` <- function(a, b) if (!is.null(a)) a else b
  
  object <- counts(dds)
  object.old <- object
  object <- object %||% object.old
  
  features <- list(features)
  features.old <- features
  
  if (is.null(x = features)) {
    stop("Missing input feature list")
  }
  features <- lapply(
    X = features,
    FUN = function(x) {
      missing.features <- setdiff(x = x, y = rownames(x = object))
      if (length(x = missing.features) > 0) {
        warning(
          "The following features are not present in the object: ",
          paste(missing.features, collapse = ", ")
        ) 
        warning(
          paste0("\n ",
                 paste(missing.features, collapse = ", "),
                 " dropped for calculating the geneset score."
          )
        )
      }
      return(intersect(x = x, y = rownames(x = object)))
    }
  )
  
  geneset.length <- length(x = features)
  
  pool <- pool %||% rownames(x = object)
  data.avg <- Matrix::rowMeans(x = object[pool, ])
  data.avg <- data.avg[order(data.avg)]
  data.cut <- cut_number(x = data.avg + rnorm(n = length(data.avg))/1e30, n = nbin, labels = FALSE, right = FALSE)
  names(x = data.cut) <- names(x = data.avg)
  ctrl.use <- vector(mode = "list", length = geneset.length)
  for (i in 1:geneset.length) {
    features.use <- features[[i]]
    for (j in 1:length(x = features.use)) {
      ctrl.use[[i]] <- c(
        ctrl.use[[i]],
        names(x = sample(
          x = data.cut[which(x = data.cut == data.cut[features.use[j]])],
          size = ctrl,
          replace = FALSE
        ))
      )
    }
  }
  ctrl.use <- lapply(X = ctrl.use, FUN = unique)
  ctrl.scores <- matrix(
    data = numeric(length = 1L),
    nrow = length(x = ctrl.use),
    ncol = ncol(x = object)
  )
  for (i in 1:length(ctrl.use)) {
    features.use <- ctrl.use[[i]]
    ctrl.scores[i, ] <- Matrix::colMeans(x = object[features.use, ])
  }
  features.scores <- matrix(
    data = numeric(length = 1L),
    nrow = geneset.length,
    ncol = ncol(x = object)
  )
  for (i in 1:geneset.length) {
    features.use <- features[[i]]
    data.use <- object[features.use, , drop = FALSE]
    features.scores[i, ] <- Matrix::colMeans(x = data.use)
  }
  features.scores.use <- features.scores - ctrl.scores
  rownames(x = features.scores.use) <- paste0(name, 1:geneset.length)
  features.scores.use <- as.data.frame(x = t(x = features.scores.use))
  
  range01 <- lapply(
    X = features.scores.use,
    FUN = function(x) {
      range01 <- (x-min(x))/(max(x)-min(x))
    }
  )
  
  range01 <- as.data.frame(x = range01)
  rownames(x = range01) <- colnames(object)
  
  colData(dds) <- cbind(colData(dds), range01)
  
  return(dds)
}

query_exp_thym <- GDCquery(
  project = "TCGA-THYM",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification", 
  workflow.type = "STAR - Counts"
)

GDCdownload(query_exp_thym)
exp_thym <- GDCprepare(
  query = query_exp_thym
)

rowRanges(exp_thym)
data <- exp_thym
row.names(data) <- rowRanges(data)$gene_name

colnames(data[,duplicated(substr(colnames(data), 1, 12))])

data <- data[,!duplicated(substr(colnames(data), 1, 12))]
colnames(data) <- substr(colnames(data), 1, 12)

data <- data[,data$patient %in% c("TCGA-3G-AB0O","TCGA-3G-AB14","TCGA-3G-AB19","TCGA-3Q-A9WF","TCGA-3S-A8YW","TCGA-3S-AAYX","TCGA-3T-AA9L","TCGA-4V-A9QI","TCGA-4V-A9QJ","TCGA-4V-A9QL",
                                  "TCGA-4V-A9QM","TCGA-4V-A9QN","TCGA-4V-A9QR","TCGA-4V-A9QS","TCGA-4V-A9QT","TCGA-4V-A9QU","TCGA-4V-A9QW","TCGA-4V-A9QX","TCGA-4X-A9F9","TCGA-4X-A9FA",
                                  "TCGA-4X-A9FB","TCGA-4X-A9FC","TCGA-4X-A9FD","TCGA-5G-A9ZZ","TCGA-5K-AAAP","TCGA-5U-AB0D","TCGA-5U-AB0E","TCGA-5U-AB0F","TCGA-5V-A9RR","TCGA-X7-A8D6",
                                  "TCGA-X7-A8D7",
                                  "TCGA-X7-A8D8",
                                  "TCGA-X7-A8D9",
                                  "TCGA-X7-A8DB",
                                  "TCGA-X7-A8DD",
                                  "TCGA-X7-A8DE",
                                  "TCGA-X7-A8DF",
                                  "TCGA-X7-A8DG",
                                  "TCGA-X7-A8DJ",
                                  "TCGA-X7-A8M0",
                                  "TCGA-X7-A8M1",
                                  "TCGA-X7-A8M3",
                                  "TCGA-X7-A8M4",
                                  "TCGA-X7-A8M5",
                                  "TCGA-X7-A8M6",
                                  "TCGA-X7-A8M7",
                                  "TCGA-X7-A8M8",
                                  "TCGA-XH-A853",
                                  "TCGA-XM-A8R8",
                                  "TCGA-XM-A8R9",
                                  "TCGA-XM-A8RB",
                                  "TCGA-XM-A8RC",
                                  "TCGA-XM-A8RD",
                                  "TCGA-XM-A8RE",
                                  "TCGA-XM-A8RF",
                                  "TCGA-XM-A8RG",
                                  "TCGA-XM-A8RH",
                                  "TCGA-XM-A8RI",
                                  "TCGA-XM-A8RL",
                                  "TCGA-XM-AAZ1",
                                  "TCGA-XM-AAZ2",
                                  "TCGA-XM-AAZ3",
                                  "TCGA-XU-A92O",
                                  "TCGA-XU-A92Q",
                                  "TCGA-XU-A92T",
                                  "TCGA-XU-A92U",
                                  "TCGA-XU-A92V",
                                  "TCGA-XU-A92W",
                                  "TCGA-XU-A92X",
                                  "TCGA-XU-A92Y",
                                  "TCGA-XU-A92Z",
                                  "TCGA-XU-A930",
                                  "TCGA-XU-A931",
                                  "TCGA-XU-A932",
                                  "TCGA-XU-A933",
                                  "TCGA-XU-A936",
                                  "TCGA-XU-AAXW",
                                  "TCGA-XU-AAXX",
                                  "TCGA-XU-AAXY",
                                  "TCGA-XU-AAXZ",
                                  "TCGA-XU-AAY0",
                                  "TCGA-XU-AAY1",
                                  "TCGA-YT-A95D",
                                  "TCGA-YT-A95E",
                                  "TCGA-YT-A95F",
                                  "TCGA-YT-A95G",
                                  "TCGA-YT-A95H",
                                  "TCGA-ZB-A961",
                                  "TCGA-ZB-A962",
                                  "TCGA-ZB-A963",
                                  "TCGA-ZB-A964",
                                  "TCGA-ZB-A965",
                                  "TCGA-ZB-A966",
                                  "TCGA-ZB-A969",
                                  "TCGA-ZB-A96A",
                                  "TCGA-ZB-A96B",
                                  "TCGA-ZB-A96C",
                                  "TCGA-ZB-A96D",
                                  "TCGA-ZB-A96E",
                                  "TCGA-ZB-A96F",
                                  "TCGA-ZB-A96G",
                                  "TCGA-ZB-A96H",
                                  "TCGA-ZB-A96I",
                                  "TCGA-ZB-A96K",
                                  "TCGA-ZB-A96L",
                                  "TCGA-ZB-A96M",
                                  "TCGA-ZB-A96O",
                                  "TCGA-ZB-A96P",
                                  "TCGA-ZB-A96Q",
                                  "TCGA-ZB-A96R",
                                  "TCGA-ZB-A96V",
                                  "TCGA-ZC-AAA7",
                                  "TCGA-ZC-AAAA",
                                  "TCGA-ZC-AAAF",
                                  "TCGA-ZC-AAAH",
                                  "TCGA-ZL-A9V6",
                                  "TCGA-ZT-A8OM"
)]

exp_thym$patient
data$patient

#write.xlsx(assay(data), 'D:/Direder/Projekte/D003_10xThymom/D003.R2/Berechnung/D003R2_c1/Lists/201123_TCGATHYM_HTSeq_raocounts.xlsx')

#clinical <- GDCquery_clinic(project = "TCGA-THYM", type = "clinical")
#query_cli <- GDCquery(
#  project = "TCGA-THYM", 
#  data.category = "Clinical", data.format = "BCR XML")
#
#GDCdownload(query_cli)
#query_cli <- GDCprepare_clinic(query_cli, clinical.info = "patient")
clinical <- read_excel('D:/Direder/Projekte/D003_10xThymom/Revision/Berechnung/D003R_c1/TCGA_thymoma/clinical.xlsx')
head(clinical)
clinical <- clinical[!duplicated(clinical$Patient_Short_Barcode),]
rownames(clinical) <- clinical$Patient_Short_Barcode
head(clinical)
data$patient
patientorder<- data$patient
clinical <- clinical %>% arrange(factor(Patient_Short_Barcode, levels = patientorder))
#clinical <- clinical [patientorder,]
head(clinical)
rownames(clinical) <- clinical$Patient_Short_Barcode
rownames(clinical)

colData(data) <- cbind(colData(data), clinical)
data <- data[!is.na(row.names(data)),!is.na(data$`GTF2I_positivity`)]
data <- data[!duplicated(row.names(data)),]
data$GTF2I_positivity <- 
  factor(data$GTF2I_positivity, levels = c("NO", "YES"))

setwd(dir ="D:/Direder/Projekte/D003_10xThymom/D003.R2/Berechnung/D003R2_c1/TCGA_thymoma")

save(data,file="Bulkseq_TCGA_Thym.RData")



ddsSE <- DESeqDataSet(data, design = ~ GTF2I_positivity)
keep <- rowMeans(counts(ddsSE)) >= 5
ddsSE <- ddsSE[keep,]
ddsSE <- DESeq(ddsSE)

resultsNames(ddsSE)

res <- results(ddsSE, name = "GTF2I_positivity_YES_vs_NO", alpha = 0.01, lfcThreshold = 2)
resOrdered <- res[order(res$padj),]


dea <- as.data.frame(resOrdered)
write.xlsx(dea, file = "DESeq2.res.xlsx")
summary(res)
dea
dea[c("GFRA3","POTEF","GPM6A","SNED1","LEPR","BCAM","LINC00408","ZNF534"),]
plotMA(res, alpha = 0.1)

degs <- subset(resOrdered, padj < 0.01)
#degs <- degs[abs(degs$log2FoldChange) > 2,]
genes.up <- row.names(degs[degs$log2FoldChange >= 0,])
genes.down <- row.names(degs[degs$log2FoldChange < 0,])
write(genes.up, file = "DESeq2.padj0.1.lfc1.upgenes.txt")
write(genes.down, file = "DESeq2.padj0.1lfc1.downgenes.txt")

# pheatmap
vsd <- vst(ddsSE)
meanSdPlot(assay(vsd))
df <- as.data.frame(colData(ddsSE)[,c("GTF2I_positivity", "FINAL_diagnosis")])
row.names(df) <- colnames(ddsSE)
#levels(df$FINAL_diagnosis) <- c( "Thymoma; Type A", "Thymoma; Type AThymoma; Type AB", "Thymoma; Type AB", "Thymoma; Type B1", "Thymoma; Type B1Thymoma; Type B2",
#                                                          "Thymoma; Type B2", "Thymoma; Type B2Thymoma; Type B3", "Thymoma; Type B3", "Thymoma; Type C" )
#df$FINAL_diagnosis <- revalue(df$FINAL_diagnosis, c("Thymoma; Type A" = "A", "Thymoma; Type AThymoma; Type AB" = "A;AB", "Thymoma; Type AB" = "AB", "Thymoma; Type B1" = "B1", "Thymoma; Type B1Thymoma; Type B2" = "B1;B2",
#                                                                                                      "Thymoma; Type B2" = "B2", "Thymoma; Type B2Thymoma; Type B3" = "B2;B3", "Thymoma; Type B3" = "B3", "Thymoma; Type C" = "C"))

# Specify colors
# cpal <- viridis::plasma(9)
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
cpal <- gg_color_hue(9)
ann_colors = list(
  GTF2I_positivity = c(YES = "black", NO = "grey"),
  FINAL_diagnosis = c(A = "#e81123", AB="orange", B1="#00188f", B2="#00188f", B3="#00188f", CA="#009e49", "CA-LCNE"="#009e49", "CA-SQ"="#009e49", "CA-UN"="#009e49", "MN-T"="#bad80a")
)

zvsd <- (assay(vsd) - rowMeans(assay(vsd))) / rowSds(assay(vsd))

my_col_order <- c( "TCGA-3G-AB19","TCGA-4X-A9FA","TCGA-XM-A8RG","TCGA-XU-A92Q",
                   "TCGA-ZB-A962","TCGA-ZB-A964",
                   "TCGA-ZB-A965","TCGA-ZB-A969",
                   "TCGA-ZB-A96E","TCGA-ZL-A9V6",
                   "TCGA-4V-A9QL","TCGA-4V-A9QM",
                   "TCGA-4V-A9QU","TCGA-4V-A9QW",
                   "TCGA-4V-A9QX","TCGA-5K-AAAP",
                   "TCGA-5U-AB0F","TCGA-5V-A9RR",
                   "TCGA-X7-A8D8","TCGA-X7-A8D9",
                   "TCGA-X7-A8DF","TCGA-X7-A8DJ",
                   "TCGA-X7-A8M1","TCGA-XM-A8RF",
                   "TCGA-XM-AAZ2","TCGA-XM-AAZ3",
                   "TCGA-XU-A92O","TCGA-XU-A92T",
                   "TCGA-XU-A92U","TCGA-XU-A92V",
                   "TCGA-XU-AAY0","TCGA-YT-A95D",
                   "TCGA-YT-A95E","TCGA-YT-A95F",
                   "TCGA-ZB-A963","TCGA-ZB-A96A",
                   "TCGA-ZB-A96C","TCGA-ZB-A96D",
                   "TCGA-ZB-A96H","TCGA-ZB-A96K",
                   "TCGA-ZB-A96L","TCGA-ZB-A96M",
                   "TCGA-ZB-A96Q","TCGA-ZC-AAAF",
                   "TCGA-X7-A8M5","TCGA-X7-A8M0",
                   "TCGA-4X-A9FC","TCGA-X7-A8D6",
                   "TCGA-X7-A8DD","TCGA-X7-A8M3",
                   "TCGA-X7-A8M4","TCGA-X7-A8M8",
                   "TCGA-XH-A853","TCGA-XM-A8RE",
                   "TCGA-XM-A8RH","TCGA-XU-AAY1",
                   "TCGA-YT-A95H","TCGA-ZB-A96B",
                   "TCGA-ZB-A96G","TCGA-ZT-A8OM",
                   "TCGA-3G-AB0O","TCGA-3Q-A9WF",
                   "TCGA-3S-AAYX","TCGA-3T-AA9L",
                   "TCGA-4V-A9QT","TCGA-X7-A8D7",
                   "TCGA-X7-A8DB","TCGA-XM-A8RB",
                   "TCGA-XM-A8RL","TCGA-ZB-A96O",
                   "TCGA-ZB-A96P","TCGA-ZB-A96R",
                   "TCGA-4V-A9QI","TCGA-4V-A9QJ",
                   "TCGA-4V-A9QR","TCGA-4V-A9QS",
                   "TCGA-4X-A9F9","TCGA-4X-A9FB",
                   "TCGA-X7-A8DE","TCGA-X7-A8DG",
                   "TCGA-X7-A8M6","TCGA-X7-A8M7",
                   "TCGA-XM-A8RC","TCGA-XM-A8RD",
                   "TCGA-XM-AAZ1","TCGA-XU-A92W",
                   "TCGA-XU-A92X","TCGA-XU-A92Y",
                   "TCGA-XU-A932","TCGA-XU-AAXW",
                   "TCGA-XU-AAXX","TCGA-XU-AAXY",
                   "TCGA-XU-AAXZ","TCGA-YT-A95G",
                   "TCGA-ZB-A96F","TCGA-ZC-AAAA",
                   "TCGA-ZC-AAAH","TCGA-3G-AB14",
                   "TCGA-4V-A9QN","TCGA-4X-A9FD",
                   "TCGA-5U-AB0E","TCGA-XM-A8RI",
                   "TCGA-XU-A930","TCGA-XU-A931",
                   "TCGA-ZB-A961","TCGA-ZB-A96I",
                   "TCGA-5G-A9ZZ","TCGA-5U-AB0D",
                   "TCGA-3S-A8YW","TCGA-XM-A8R8",
                   "TCGA-XU-A92Z","TCGA-ZB-A966",
                   "TCGA-XU-A933","TCGA-XU-A936",
                   "TCGA-ZB-A96V","TCGA-ZC-AAA7",
                   "TCGA-XM-A8R9")


pdf("pheatmap.DESeq2.padj0.01.lfc2.pdf", width = 16, height = 16) 
pheatmap(assay(vsd)[row.names(degs),my_col_order], annotation_col=df, color = viridis(10), annotation_colors = ann_colors,cluster_cols= F , 
         clustering_distance_rows = "correlation", clustering_method = "ward.D2")
dev.off()


#convert to NCBI Symbol & remove duplicates
genes.up.modified <- lapply(genes.up, function(x){
  Signlist1=alias2SymbolUsingNCBI(x,paste0(c("D:/Direder/SC-EnsembID/Homo_sapiens.gene_info")),required.columns = c("GeneID","Symbol"))
  Signlist2=as.data.frame(cbind(Signlist1,x))
  colnames(Signlist2)[3]=c("rname")
  Signlist2$new = ifelse(is.na(Signlist2$Symbol), as.character(Signlist2$rname), Signlist2$Symbol)
  x= as.character(Signlist2$new)
})
genes.down.modified <- lapply(genes.down, function(x){
  Signlist1=alias2SymbolUsingNCBI(x,paste0(c("D:/Direder/SC-EnsembID/Homo_sapiens.gene_info")),required.columns = c("GeneID","Symbol"))
  Signlist2=as.data.frame(cbind(Signlist1,x))
  colnames(Signlist2)[3]=c("rname")
  Signlist2$new = ifelse(is.na(Signlist2$Symbol), as.character(Signlist2$rname), Signlist2$Symbol)
  x= as.character(Signlist2$new)
})

#keep only genes present in seurat object
#up
genelist.th.total<-as.list(rownames(th.total))
genelist.th.total<-as.character(genelist.th.total)
genes.up.modified<-as.character(genes.up.modified)
genes.up.newlist <- c(genes.up.modified,genelist.th.total)
sum(duplicated(genes.up.newlist))
genes.up.newlist<-genes.up.newlist[duplicated(genes.up.newlist) | duplicated(genes.up.newlist, fromLast=TRUE)]
genes.up.newlist <- sapply(genes.up.newlist, unique)
sum(duplicated(genes.up.newlist))
genes.up.newlist<- unique(genes.up.newlist)
write.xlsx(genes.up.newlist, "Lists/genes.up.newlist.xlsx")

#down
genelist.th.total<-as.list(rownames(th.total))
genelist.th.total<-as.character(genelist.th.total)
genes.down.modified<-as.character(genes.down.modified)
genes.down.newlist <- c(genes.down.modified,genelist.th.total)
sum(duplicated(genes.down.newlist))
genes.down.newlist<-genes.down.newlist[duplicated(genes.down.newlist) | duplicated(genes.down.newlist, fromLast=TRUE)]
genes.down.newlist <- sapply(genes.down.newlist, unique)
sum(duplicated(genes.down.newlist))
genes.down.newlist<- unique(genes.down.newlist)
write.xlsx(genes.down.newlist, "Lists/genes.down.newlist.xlsx")

th.total <- AddModuleScore(th.total, features = genes.up.newlist, name = "genes.up.newlist module score")
th.total <- AddModuleScore(th.total, features = genes.down.newlist, name = "genes.down.newlist module score")


DEG_bulkseq_signature <- c("genes.down.newlist module score1","genes.up.newlist module score1")

Idents(th.total)<-"condition"

pdf("Graphs/Main/Bulkseq_condition_comparison.pdf")
DotPlot(object = th.total, features = DEG_bulkseq_signature,assay="RNA") + coord_flip()+ scale_color_viridis(option = "B")+theme(axis.text.x = element_text(angle = 45)) + ggtitle("Main_modulescore") + scale_y_discrete(guide = guide_axis(angle = 90))
VlnPlot(th.total, features = DEG_bulkseq_signature, pt.size=0)
dev.off()

Idents(th.total)<-"tissue"

pdf("Graphs/Main/Bulkseq_tissue_comparison.pdf")
DotPlot(object = th.total, features = DEG_bulkseq_signature,assay="RNA") + coord_flip()+ scale_color_viridis(option = "B")+theme(axis.text.x = element_text(angle = 45)) + ggtitle("Main_modulescore") + scale_y_discrete(guide = guide_axis(angle = 90))
VlnPlot(th.total, features = DEG_bulkseq_signature, pt.size=0)
dev.off()

Idents(th.total)<-"sample"

pdf("Graphs/Main/Bulkseq_sample_comparison.pdf", width = 15, height=4)
DotPlot(object = th.total, features = DEG_bulkseq_signature,assay="RNA") + coord_flip()+ scale_color_viridis(option = "H")+theme(axis.text.x = element_text(angle = 45)) + ggtitle("Main_modulescore") + scale_y_discrete(guide = guide_axis(angle = 90))
VlnPlot(th.total, features = DEG_bulkseq_signature, pt.size=0)
dev.off()


pdf("Heatmap_adapted_Module_Genelist_GTF2I.pdf", height=10, width=10)
pheatmap(assay(vsd)[genes.up.newlist,], annotation_col=df, color = viridis(10), annotation_colors = ann_colors,
         clustering_distance_rows = "correlation", clustering_method = "ward.D2")

pheatmap(assay(vsd)[genes.down.newlist,], annotation_col=df, color = viridis(10), annotation_colors = ann_colors,
         clustering_distance_rows = "correlation", clustering_method = "ward.D2")
dev.off()

pdf("Heatmap_adapted_Module_Genelist_GTF2I_full.pdf", height=25, width=10)
pheatmap(assay(vsd)[genes.up.newlist,], annotation_col=df, color = viridis(10), annotation_colors = ann_colors,
         clustering_distance_rows = "correlation", clustering_method = "ward.D2")

pheatmap(assay(vsd)[genes.down.newlist,], annotation_col=df, color = viridis(10), annotation_colors = ann_colors,
         clustering_distance_rows = "correlation", clustering_method = "ward.D2")
dev.off()

# PCA
levels(df$FINAL_diagnosis) <- c("A","AB","B1","B2","B3","CA","CA-LCNE","CA-SQ","CA-UN","MN-T")
df$FINAL_diagnosis <- revalue(df$FINAL_diagnosis, c("A"="A","AB"="AB","B1"="B","B2"="B","B3"="B","CA"="B","CA-LCNE"="C","CA-SQ"="C","CA-UN"="C","MN-T"="MNT"))
levels(df$FINAL_diagnosis) <- c("A","AB","B","C","MNT")
#

levels(ddsSE$FINAL_diagnosis) <- c("A","AB","B1","B2","B3","CA","CA-LCNE","CA-SQ","CA-UN","MN-T")
ddsSE$FINAL_diagnosis <- revalue(ddsSE$FINAL_diagnosis, c("A"="A","AB"="AB","B1"="B","B2"="B","B3"="B","CA"="B","CA-LCNE"="C","CA-SQ"="C","CA-UN"="C","MN-T"="MNT"))
levels(ddsSE$FINAL_diagnosis) <- c("A","AB","B","C","MNT")

vsd <- vst(ddsSE)
meanSdPlot(assay(vsd))

pdf("GTF2I_PCA.pdf", width = 10, height = 5)
p1 <- plotPCA(vsd, intgroup=c("GTF2I_positivity")) + 
  stat_ellipse(type = "norm", level = 0.67, geom = "polygon", alpha = 0.1, aes(fill = GTF2I_positivity)) +
  ggtitle("GTF2I_positive") + scale_color_manual(values = c("grey","black"))+
  theme(plot.title = element_text(hjust = 0.5), legend.position = "bottom")+ theme_light()
p2 <- plotPCA(vsd, intgroup=c("FINAL_diagnosis")) +  
  ggtitle("Histological type") + scale_color_manual(values = c("#e81123","#ff8c00","#00188f","#009e49","#bad80a"))+
  theme(plot.title = element_text(hjust = 0.5), legend.position = "bottom")+ theme_light()
gridExtra::grid.arrange(p1, p2, nrow = 1)
dev.off()


# volcano plot
res.vol <- res[res$pvalue < 0.05,]

# all DEGs
pdf("GTF2I_vol.alls.pdf", width = 20, height = 12)
res1 <- lfcShrink(ddsSE, contrast = c("GTF2I_positivity", 'YES', 'NO'), 
                  res=res, type = 'normal')
EnhancedVolcano(res.vol,
                lab = rownames(res.vol),
                x = "log2FoldChange",
                y = "padj",
                pCutoff = 0.1,
                FCcutoff = 1,
                xlab = bquote(~Log[2]~ "fold change"),
                ylab = bquote(~-Log[10]~adjusted~italic(P)),
                drawConnectors = TRUE,
                widthConnectors = 0.2,
                cutoffLineWidth = 0,
                ylim = c(0,13),
                xlim = c(-20,10),
                col=c("black", "black", "black", "red3"),)
dev.off()

pdf("GTF2I_vol.alls.annot30.pdf", width = 6, height = 6)
EnhancedVolcano(res.vol,
                lab = rownames(res.vol),
                x = "log2FoldChange",
                y = "padj",
                pCutoff = 0.1,
                FCcutoff = 1,
                xlab = bquote(~Log[2]~ "fold change"),
                ylab = bquote(~-Log[10]~adjusted~italic(P)),
                drawConnectors = TRUE,
                widthConnectors = 0.2,
                cutoffLineWidth = 0,
                selectLab = row.names(resOrdered[1:30,]),
                ylim = c(0,13),
                xlim = c(-10,10),
                col=c("black", "black", "black", "red3"),)
dev.off()

# stacked bar plot
GTF2I.who <- df
GTF2I.who <- tibble::as_tibble(GTF2I.who)
# GTF2I.who$GTF2I_positivity
# GTF2I.who <- GTF2I.who %>% mutate(count = 1)
GTF2I.who <- GTF2I.who %>% 
  dplyr::group_by(FINAL_diagnosis, GTF2I_positivity) %>%
  dplyr::tally()???# %>%
# dplyr::filter(n > 5)

ggplot(GTF2I.who, aes(x = FINAL_diagnosis, y = n, fill = GTF2I_positivity)) + 
  geom_bar( stat='identity')+ scale_fill_manual(values = c("grey","black"))
ggsave("GTF2I.WHO.bar.pdf", width = 7, height = 2)



#Modulescore


ddsSE <- AddGeneSetScore(ddsSE, 
                         features = genes.up,
                         ctrl = 5,
                         name = 'GTF2I_genesUP'
)



ddsSE <- AddGeneSetScore(ddsSE, 
                         features = genes.down,
                         ctrl = 5,
                         name = 'GTF2I_genesdown'
)

pdf("Modulescore_genes_GTF2I_Bulkseq.pdf", width=10, height=8)
gg1 <- ggplot(as.data.frame(colData(ddsSE)), aes(x = FINAL_diagnosis, y = GTF2I_genesUP1, color = GTF2I_positivity)) + 
  geom_point() +
  theme_bw() + scale_color_manual(values = c("grey","black"))+ 
  ylab("Regulation of GTF2I_UP regulated genes") + 
  theme(text = element_text(size = 18))

gg2 <- ggplot(as.data.frame(colData(ddsSE)), aes(x = FINAL_diagnosis, y = GTF2I_genesdown1, color = GTF2I_positivity)) + 
  geom_point() +
  theme_bw() + scale_color_manual(values = c("grey","black"))+ 
  ylab("Regulation of GTF2I_DOWN regulated genes") + 
  theme(text = element_text(size = 18))
gg1+gg2
dev.off()

#(45)_Mutations_TCGA####

library(TCGAbiolinks)
library(SummarizedExperiment)
library(maftools)
library(TCGAmutations)
library(dplyr)
library(DT)

TCGAmutations::tcga_available()

tcga_thym_mc3 <- TCGAmutations::tcga_load(study = "THYM")
tcga_thym_mc3
getSampleSummary(x = tcga_thym_mc3)

getGeneSummary(x = tcga_thym_mc3)

getClinicalData(x = tcga_thym_mc3)[1:10, 1:10]

frame()
plotmafSummary(maf = tcga_thym_mc3, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)

plot.new()  
oncoplot(maf = tcga_thym_mc3) #, top = 5, fontSize = 12)
table(tcga_thym_mc3@clinical.data$histological_type)
sample_order = sample(tcga_thym_mc3@clinical.data$histological_type )
plot.new()  
pdf("D:/Direder/Projekte/D003_10xThymom/Revision/Berechnung/D003R_c1/TCGA_thymoma/Mutation_plot_thym.pdf", height=10, width=20)
oncoplot(maf = tcga_thym_mc3, clinicalFeatures = 'histological_type', sortByAnnotation = T)
dev.off()
thym.sig = oncodrive(maf = tcga_thym_mc3, AACol = 'Protein_Change', minMut = 5, pvalMethod = 'zscore')


#
#(46)_pathways####
#(47)_pw_FB_pathways####
net <- get_progeny(organism = 'human', top = 500)
net

data <- FB_sub
# Extract the normalized log-transformed counts
mat <- as.matrix(data@assays$RNA$data)

# Run mlm
acts <- run_mlm(mat=mat, net=net, .source='source', .target='target',
                .mor='weight', minsize = 5)
acts

#data$condition <- data$tissue
#levels(data$condition)

# Extract mlm and store it in pathwaysmlm in data
data[['pathwaysmlm']] <- acts %>%
  pivot_wider(id_cols = 'source', names_from = 'condition',
              values_from = 'score') %>%
  column_to_rownames('source') %>%
  Seurat::CreateAssayObject(.)

# Change assay
DefaultAssay(object = data) <- "pathwaysmlm"

# Scale the data
data <- ScaleData(data)
data@assays$pathwaysmlm@data <- data@assays$pathwaysmlm@scale.data

p1 <- DimPlot(data, reduction = "umap", label = TRUE, pt.size = 0.5) + 
  NoLegend() + ggtitle('Cell types')
p2 <- (FeaturePlot(data, features = c("Trail")) & 
         scale_colour_gradient2(low = 'blue', mid = 'white', high = 'red')) +
  ggtitle('Trail activity')
pdf("Graphs/Main/Pathways_Trailactivity_FB.pdf", width=15, height=10)
p1 | p2
dev.off()

# Extract activities from object as a long dataframe
df <- t(as.matrix(data@assays$pathwaysmlm@data)) %>%
  as.data.frame() %>%
  mutate(cluster = Idents(data)) %>%
  pivot_longer(cols = -cluster, names_to = "source", values_to = "score") %>%
  group_by(cluster, source) %>%
  summarise(mean = mean(score))

# Transform to wide matrix
top_acts_mat <- df %>%
  pivot_wider(id_cols = 'cluster', names_from = 'source',
              values_from = 'mean') %>%
  column_to_rownames('cluster') %>%
  as.matrix()

write.xlsx(top_acts_mat, "Lists/Pathways_FB.xlsx")


# Choose color palette
palette_length = 100
my_color = colorRampPalette(c("Darkblue", "white","red"))(palette_length)

my_breaks <- c(seq(-2, 0, length.out=ceiling(palette_length/2) + 1),
               seq(0.05, 2, length.out=floor(palette_length/2)))

# Plot
pdf("Graphs/Main/Pathways_Heatmap_FB.pdf")
pheatmap(top_acts_mat, border_color = NA, color=my_color, breaks = my_breaks) 
dev.off()



#(48)_pw_th.total_pathways####
#pw_th.total_pathways#
Idents(th.total) <- th.total$celltype.full
th.total$celltype.tissue <- paste(Idents(th.total), th.total$tissue, sep = "_")
Idents(th.total)<-th.total$celltype.tissue
levels(th.total)
Idents(th.total) <- th.total$celltype.tissue
th.total$celltype.tissue.donor <- paste(Idents(th.total), th.total$donor, sep = "_")
Idents(th.total)<-th.total$celltype.tissue.donor
levels(th.total)
Idents(th.total)<-th.total$tissue
levels(th.total)

#(49)_pw_prenatal_thymus####
prenatal_Thymus_pathway <- subset(th.total, idents = "prenatal_Thymus")
Idents(prenatal_Thymus_pathway)<- prenatal_Thymus_pathway$celltype.tissue.donor
levels(prenatal_Thymus_pathway)
table(prenatal_Thymus_pathway$celltype.tissue.donor)
net <- get_progeny(organism = 'human', top = 500)
net

data <- prenatal_Thymus_pathway
# Extract the normalized log-transformed counts
mat <- as.matrix(data@assays$RNA$data)

# Run mlm
acts <- run_mlm(mat=mat, net=net, .source='source', .target='target',
                .mor='weight', minsize = 5)
acts

# Extract mlm and store it in pathwaysmlm in data
data[['pathwaysmlm']] <- acts %>%
  pivot_wider(id_cols = 'source', names_from = 'condition',
              values_from = 'score') %>%
  column_to_rownames('source') %>%
  Seurat::CreateAssayObject(.)

# Change assay
DefaultAssay(object = data) <- "pathwaysmlm"

# Scale the data
data <- ScaleData(data)
data@assays$pathwaysmlm@data <- data@assays$pathwaysmlm@scale.data

p1 <- DimPlot(data, reduction = "umap", label = TRUE, pt.size = 0.5) + 
  NoLegend() + ggtitle('Cell types')
p2 <- (FeaturePlot(data, features = c("Trail")) & 
         scale_colour_gradient2(low = 'blue', mid = 'white', high = 'red')) +
  ggtitle('Trail activity')
pdf("Graphs/Main/Pathways_Trailactivity_prenatal_Thymus.pdf", width=15, height=10)
p1 | p2
dev.off()

# Extract activities from object as a long dataframe
df <- t(as.matrix(data@assays$pathwaysmlm@data)) %>%
  as.data.frame() %>%
  mutate(cluster = Idents(data)) %>%
  pivot_longer(cols = -cluster, names_to = "source", values_to = "score") %>%
  group_by(cluster, source) %>%
  summarise(mean = mean(score))

# Transform to wide matrix
top_acts_mat <- df %>%
  pivot_wider(id_cols = 'cluster', names_from = 'source',
              values_from = 'mean') %>%
  column_to_rownames('cluster') %>%
  as.matrix()

write.xlsx(top_acts_mat, "Lists/Pathways_prenatal_Thymus.xlsx")

#plot mean values
daten<-read_excel("Lists/Pathways_prenatal_Thymus_mean.xlsx")
rownames(daten) <- daten[[1]]
rn <- rownames(daten)
daten <- daten[-1]
top_acts_mat <- as.matrix(daten)
rownames(top_acts_mat) <- rn

# Choose color palette
palette_length = 100
my_color = colorRampPalette(c("Darkblue", "white","red"))(palette_length)

my_breaks <- c(seq(-2, 0, length.out=ceiling(palette_length/2) + 1),
               seq(0.05, 2, length.out=floor(palette_length/2)))

# Plot
pdf("Graphs/Main/Pathways_Heatmap_prenatal_Thymus.pdf")
pheatmap(top_acts_mat, border_color = NA, color=my_color, breaks = my_breaks) 
dev.off()





#(50)_pw_pediatric_Thymus####
pediatric_Thymus_pathway <- subset(th.total, idents = "pediatric_Thymus")
Idents(pediatric_Thymus_pathway)<- pediatric_Thymus_pathway$celltype.tissue.donor
table(pediatric_Thymus_pathway$celltype.tissue.donor)
net <- get_progeny(organism = 'human', top = 500)
net

data <- pediatric_Thymus_pathway
# Extract the normalized log-transformed counts
mat <- as.matrix(data@assays$RNA$data)

# Run mlm
acts <- run_mlm(mat=mat, net=net, .source='source', .target='target',
                .mor='weight', minsize = 5)
acts

# Extract mlm and store it in pathwaysmlm in data
data[['pathwaysmlm']] <- acts %>%
  pivot_wider(id_cols = 'source', names_from = 'condition',
              values_from = 'score') %>%
  column_to_rownames('source') %>%
  Seurat::CreateAssayObject(.)

# Change assay
DefaultAssay(object = data) <- "pathwaysmlm"

# Scale the data
data <- ScaleData(data)
data@assays$pathwaysmlm@data <- data@assays$pathwaysmlm@scale.data

p1 <- DimPlot(data, reduction = "umap", label = TRUE, pt.size = 0.5) + 
  NoLegend() + ggtitle('Cell types')
p2 <- (FeaturePlot(data, features = c("Trail")) & 
         scale_colour_gradient2(low = 'blue', mid = 'white', high = 'red')) +
  ggtitle('Trail activity')
pdf("Graphs/Main/Pathways_Trailactivity_pediatric_Thymus.pdf", width=15, height=10)
p1 | p2
dev.off()

# Extract activities from object as a long dataframe
df <- t(as.matrix(data@assays$pathwaysmlm@data)) %>%
  as.data.frame() %>%
  mutate(cluster = Idents(data)) %>%
  pivot_longer(cols = -cluster, names_to = "source", values_to = "score") %>%
  group_by(cluster, source) %>%
  summarise(mean = mean(score))

# Transform to wide matrix
top_acts_mat <- df %>%
  pivot_wider(id_cols = 'cluster', names_from = 'source',
              values_from = 'mean') %>%
  column_to_rownames('cluster') %>%
  as.matrix()

write.xlsx(top_acts_mat, "Lists/Pathways_pediatric_Thymus.xlsx")

#plot mean values
daten<-read_excel("Lists/Pathways_pediatric_Thymus_mean.xlsx")
rownames(daten) <- daten[[1]]
rn <- rownames(daten)
daten <- daten[-1]
top_acts_mat <- as.matrix(daten)
rownames(top_acts_mat) <- rn

# Choose color palette
palette_length = 100
my_color = colorRampPalette(c("Darkblue", "white","red"))(palette_length)

my_breaks <- c(seq(-2, 0, length.out=ceiling(palette_length/2) + 1),
               seq(0.05, 2, length.out=floor(palette_length/2)))

top_acts_mat_pediatric_Thymus <- top_acts_mat

# Plot
pdf("Graphs/Main/Pathways_Heatmap_pediatric_Thymus.pdf")
pheatmap(top_acts_mat, border_color = NA, color=my_color, breaks = my_breaks) 
dev.off()





#(51)_pw_adult_Thymus####
adult_Thymus_pathway <- subset(th.total, idents = "adult_Thymus")
Idents(adult_Thymus_pathway)<- adult_Thymus_pathway$celltype.tissue.donor
table(adult_Thymus_pathway$celltype.tissue.donor)
net <- get_progeny(organism = 'human', top = 500)
net

data <- adult_Thymus_pathway
# Extract the normalized log-transformed counts
mat <- as.matrix(data@assays$RNA$data)

# Run mlm
acts <- run_mlm(mat=mat, net=net, .source='source', .target='target',
                .mor='weight', minsize = 5)
acts

# Extract mlm and store it in pathwaysmlm in data
data[['pathwaysmlm']] <- acts %>%
  pivot_wider(id_cols = 'source', names_from = 'condition',
              values_from = 'score') %>%
  column_to_rownames('source') %>%
  Seurat::CreateAssayObject(.)

# Change assay
DefaultAssay(object = data) <- "pathwaysmlm"

# Scale the data
data <- ScaleData(data)
data@assays$pathwaysmlm@data <- data@assays$pathwaysmlm@scale.data

p1 <- DimPlot(data, reduction = "umap", label = TRUE, pt.size = 0.5) + 
  NoLegend() + ggtitle('Cell types')
p2 <- (FeaturePlot(data, features = c("Trail")) & 
         scale_colour_gradient2(low = 'blue', mid = 'white', high = 'red')) +
  ggtitle('Trail activity')
pdf("Graphs/Main/Pathways_Trailactivity_adult_Thymus.pdf", width=15, height=10)
p1 | p2
dev.off()

# Extract activities from object as a long dataframe
df <- t(as.matrix(data@assays$pathwaysmlm@data)) %>%
  as.data.frame() %>%
  mutate(cluster = Idents(data)) %>%
  pivot_longer(cols = -cluster, names_to = "source", values_to = "score") %>%
  group_by(cluster, source) %>%
  summarise(mean = mean(score))

# Transform to wide matrix
top_acts_mat <- df %>%
  pivot_wider(id_cols = 'cluster', names_from = 'source',
              values_from = 'mean') %>%
  column_to_rownames('cluster') %>%
  as.matrix()

write.xlsx(top_acts_mat, "Lists/Pathways_adult_Thymus.xlsx")

#plot mean values
daten<-read_excel("Lists/Lists/Pathways_adult_Thymus_mean.xlsx")
rownames(daten) <- daten[[1]]
rn <- rownames(daten)
daten <- daten[-1]
top_acts_mat <- as.matrix(daten)
rownames(top_acts_mat) <- rn

# Choose color palette
palette_length = 100
my_color = colorRampPalette(c("Darkblue", "white","red"))(palette_length)

my_breaks <- c(seq(-2, 0, length.out=ceiling(palette_length/2) + 1),
               seq(0.05, 2, length.out=floor(palette_length/2)))

top_acts_mat_adult_Thymus <- top_acts_mat

# Plot
pdf("Graphs/Main/Pathways_Heatmap_adult_Thymus.pdf")
pheatmap(top_acts_mat, border_color = NA, color=my_color, breaks = my_breaks) 
dev.off()





#(52)_pw_TTH####
TTH_pathway <- subset(th.total, idents = "TTH")
Idents(TTH_pathway)<- TTH_pathway$celltype.tissue.donor
table(TTH_pathway$celltype.tissue.donor)
net <- get_progeny(organism = 'human', top = 500)
net

data <- TTH_pathway
# Extract the normalized log-transformed counts
mat <- as.matrix(data@assays$RNA$data)

# Run mlm
acts <- run_mlm(mat=mat, net=net, .source='source', .target='target',
                .mor='weight', minsize = 5)
acts

# Extract mlm and store it in pathwaysmlm in data
data[['pathwaysmlm']] <- acts %>%
  pivot_wider(id_cols = 'source', names_from = 'condition',
              values_from = 'score') %>%
  column_to_rownames('source') %>%
  Seurat::CreateAssayObject(.)

# Change assay
DefaultAssay(object = data) <- "pathwaysmlm"

# Scale the data
data <- ScaleData(data)
data@assays$pathwaysmlm@data <- data@assays$pathwaysmlm@scale.data

p1 <- DimPlot(data, reduction = "umap", label = TRUE, pt.size = 0.5) + 
  NoLegend() + ggtitle('Cell types')
p2 <- (FeaturePlot(data, features = c("Trail")) & 
         scale_colour_gradient2(low = 'blue', mid = 'white', high = 'red')) +
  ggtitle('Trail activity')
pdf("Graphs/Main/Pathways_Trailactivity_TTH.pdf", width=15, height=10)
p1 | p2
dev.off()

# Extract activities from object as a long dataframe
df <- t(as.matrix(data@assays$pathwaysmlm@data)) %>%
  as.data.frame() %>%
  mutate(cluster = Idents(data)) %>%
  pivot_longer(cols = -cluster, names_to = "source", values_to = "score") %>%
  group_by(cluster, source) %>%
  summarise(mean = mean(score))

# Transform to wide matrix
top_acts_mat <- df %>%
  pivot_wider(id_cols = 'cluster', names_from = 'source',
              values_from = 'mean') %>%
  column_to_rownames('cluster') %>%
  as.matrix()

write.xlsx(top_acts_mat, "Lists/Pathways_TTH.xlsx")

#plot mean values
daten<-read_excel("Lists/Pathways_TTH_mean.xlsx")
rownames(daten) <- daten[[1]]
rn <- rownames(daten)
daten <- daten[-1]
top_acts_mat <- as.matrix(daten)
rownames(top_acts_mat) <- rn

# Choose color palette
palette_length = 100
my_color = colorRampPalette(c("Darkblue", "white","red"))(palette_length)

my_breaks <- c(seq(-2, 0, length.out=ceiling(palette_length/2) + 1),
               seq(0.05, 2, length.out=floor(palette_length/2)))

top_acts_mat_adult_TTH <- top_acts_mat

# Plot
pdf("Graphs/Main/Pathways_Heatmap_TTH.pdf")
pheatmap(top_acts_mat, border_color = NA, color=my_color, breaks = my_breaks) 
dev.off()





#(53)_pw_TET_A####
TET_A_pathway <- subset(th.total, idents = "TET_A")
Idents(TET_A_pathway)<- TET_A_pathway$celltype.tissue.donor
table(TET_A_pathway$celltype.tissue.donor)
net <- get_progeny(organism = 'human', top = 500)
net

data <- TET_A_pathway
# Extract the normalized log-transformed counts
mat <- as.matrix(data@assays$RNA$data)

# Run mlm
acts <- run_mlm(mat=mat, net=net, .source='source', .target='target',
                .mor='weight', minsize = 5)
acts

# Extract mlm and store it in pathwaysmlm in data
data[['pathwaysmlm']] <- acts %>%
  pivot_wider(id_cols = 'source', names_from = 'condition',
              values_from = 'score') %>%
  column_to_rownames('source') %>%
  Seurat::CreateAssayObject(.)

# Change assay
DefaultAssay(object = data) <- "pathwaysmlm"

# Scale the data
data <- ScaleData(data)
data@assays$pathwaysmlm@data <- data@assays$pathwaysmlm@scale.data

p1 <- DimPlot(data, reduction = "umap", label = TRUE, pt.size = 0.5) + 
  NoLegend() + ggtitle('Cell types')
p2 <- (FeaturePlot(data, features = c("Trail")) & 
         scale_colour_gradient2(low = 'blue', mid = 'white', high = 'red')) +
  ggtitle('Trail activity')
pdf("Graphs/Main/Pathways_Trailactivity_TET_A.pdf", width=15, height=10)
p1 | p2
dev.off()

# Extract activities from object as a long dataframe
df <- t(as.matrix(data@assays$pathwaysmlm@data)) %>%
  as.data.frame() %>%
  mutate(cluster = Idents(data)) %>%
  pivot_longer(cols = -cluster, names_to = "source", values_to = "score") %>%
  group_by(cluster, source) %>%
  summarise(mean = mean(score))

# Transform to wide matrix
top_acts_mat <- df %>%
  pivot_wider(id_cols = 'cluster', names_from = 'source',
              values_from = 'mean') %>%
  column_to_rownames('cluster') %>%
  as.matrix()

write.xlsx(top_acts_mat, "Lists/Pathways_TET_A.xlsx")

#plot mean values
daten<-read_excel("Lists/Pathways_TET_A_mean.xlsx")
rownames(daten) <- daten[[1]]
rn <- rownames(daten)
daten <- daten[-1]
top_acts_mat <- as.matrix(daten)
rownames(top_acts_mat) <- rn

# Choose color palette
palette_length = 100
my_color = colorRampPalette(c("Darkblue", "white","red"))(palette_length)

my_breaks <- c(seq(-2, 0, length.out=ceiling(palette_length/2) + 1),
               seq(0.05, 2, length.out=floor(palette_length/2)))

top_acts_mat_adult_TET_A <- top_acts_mat

# Plot
pdf("Graphs/Main/Pathways_Heatmap_TET_A.pdf")
pheatmap(top_acts_mat, border_color = NA, color=my_color, breaks = my_breaks) 
dev.off()





#(54)_pw_TET_AB####
TET_AB_pathway <- subset(th.total, idents = "TET_AB")
Idents(TET_AB_pathway)<- TET_AB_pathway$celltype.tissue.donor
table(TET_AB_pathway$celltype.tissue.donor)
net <- get_progeny(organism = 'human', top = 500)
net

data <- TET_AB_pathway
# Extract the normalized log-transformed counts
mat <- as.matrix(data@assays$RNA$data)

# Run mlm
acts <- run_mlm(mat=mat, net=net, .source='source', .target='target',
                .mor='weight', minsize = 5)
acts

# Extract mlm and store it in pathwaysmlm in data
data[['pathwaysmlm']] <- acts %>%
  pivot_wider(id_cols = 'source', names_from = 'condition',
              values_from = 'score') %>%
  column_to_rownames('source') %>%
  Seurat::CreateAssayObject(.)

# Change assay
DefaultAssay(object = data) <- "pathwaysmlm"

# Scale the data
data <- ScaleData(data)
data@assays$pathwaysmlm@data <- data@assays$pathwaysmlm@scale.data

p1 <- DimPlot(data, reduction = "umap", label = TRUE, pt.size = 0.5) + 
  NoLegend() + ggtitle('Cell types')
p2 <- (FeaturePlot(data, features = c("Trail")) & 
         scale_colour_gradient2(low = 'blue', mid = 'white', high = 'red')) +
  ggtitle('Trail activity')
pdf("Graphs/Main/Pathways_Trailactivity_TET_AB.pdf", width=15, height=10)
p1 | p2
dev.off()

# Extract activities from object as a long dataframe
df <- t(as.matrix(data@assays$pathwaysmlm@data)) %>%
  as.data.frame() %>%
  mutate(cluster = Idents(data)) %>%
  pivot_longer(cols = -cluster, names_to = "source", values_to = "score") %>%
  group_by(cluster, source) %>%
  summarise(mean = mean(score))

# Transform to wide matrix
top_acts_mat <- df %>%
  pivot_wider(id_cols = 'cluster', names_from = 'source',
              values_from = 'mean') %>%
  column_to_rownames('cluster') %>%
  as.matrix()

write.xlsx(top_acts_mat, "Lists/Pathways_TET_AB.xlsx")

#plot mean values
daten<-read_excel("Lists/Pathways_TET_AB_mean.xlsx")
rownames(daten) <- daten[[1]]
rn <- rownames(daten)
daten <- daten[-1]
top_acts_mat <- as.matrix(daten)
rownames(top_acts_mat) <- rn

# Choose color palette
palette_length = 100
my_color = colorRampPalette(c("Darkblue", "white","red"))(palette_length)

my_breaks <- c(seq(-2, 0, length.out=ceiling(palette_length/2) + 1),
               seq(0.05, 2, length.out=floor(palette_length/2)))

top_acts_mat_adult_TET_AB <- top_acts_mat

# Plot
pdf("Graphs/Main/Pathways_Heatmap_TET_AB.pdf")
pheatmap(top_acts_mat, border_color = NA, color=my_color, breaks = my_breaks) 
dev.off()





#(55)_pw_TET_B####
TET_B_pathway <- subset(th.total, idents = "TET_B")
Idents(TET_B_pathway)<- TET_B_pathway$celltype.tissue.donor
table(TET_B_pathway$celltype.tissue.donor)
net <- get_progeny(organism = 'human', top = 500)
net

data <- TET_B_pathway
# Extract the normalized log-transformed counts
mat <- as.matrix(data@assays$RNA$data)

# Run mlm
acts <- run_mlm(mat=mat, net=net, .source='source', .target='target',
                .mor='weight', minsize = 5)
acts

# Extract mlm and store it in pathwaysmlm in data
data[['pathwaysmlm']] <- acts %>%
  pivot_wider(id_cols = 'source', names_from = 'condition',
              values_from = 'score') %>%
  column_to_rownames('source') %>%
  Seurat::CreateAssayObject(.)

# Change assay
DefaultAssay(object = data) <- "pathwaysmlm"

# Scale the data
data <- ScaleData(data)
data@assays$pathwaysmlm@data <- data@assays$pathwaysmlm@scale.data

p1 <- DimPlot(data, reduction = "umap", label = TRUE, pt.size = 0.5) + 
  NoLegend() + ggtitle('Cell types')
p2 <- (FeaturePlot(data, features = c("Trail")) & 
         scale_colour_gradient2(low = 'blue', mid = 'white', high = 'red')) +
  ggtitle('Trail activity')
pdf("Graphs/Main/Pathways_Trailactivity_TET_B.pdf", width=15, height=10)
p1 | p2
dev.off()

# Extract activities from object as a long dataframe
df <- t(as.matrix(data@assays$pathwaysmlm@data)) %>%
  as.data.frame() %>%
  mutate(cluster = Idents(data)) %>%
  pivot_longer(cols = -cluster, names_to = "source", values_to = "score") %>%
  group_by(cluster, source) %>%
  summarise(mean = mean(score))

# Transform to wide matrix
top_acts_mat <- df %>%
  pivot_wider(id_cols = 'cluster', names_from = 'source',
              values_from = 'mean') %>%
  column_to_rownames('cluster') %>%
  as.matrix()

write.xlsx(top_acts_mat, "Lists/Pathways_TET_B.xlsx")

#plot mean values
daten<-read_excel("Lists/Pathways_TET_B_mean.xlsx")
rownames(daten) <- daten[[1]]
rn <- rownames(daten)
daten <- daten[-1]
top_acts_mat <- as.matrix(daten)
rownames(top_acts_mat) <- rn

# Choose color palette
palette_length = 100
my_color = colorRampPalette(c("Darkblue", "white","red"))(palette_length)

my_breaks <- c(seq(-2, 0, length.out=ceiling(palette_length/2) + 1),
               seq(0.05, 2, length.out=floor(palette_length/2)))

top_acts_mat_adult_TET_B <- top_acts_mat

# Plot
pdf("Graphs/Main/Pathways_Heatmap_TET_B.pdf")
pheatmap(top_acts_mat, border_color = NA, color=my_color, breaks = my_breaks) 
dev.off()





#(56)_pw_TET_C####
TET_C_pathway <- subset(th.total, idents = "TET_C")
Idents(TET_C_pathway)<- TET_C_pathway$celltype.tissue.donor
table(TET_C_pathway$celltype.tissue.donor)
net <- get_progeny(organism = 'human', top = 500)
net

data <- TET_C_pathway
# Extract the normalized log-transformed counts
mat <- as.matrix(data@assays$RNA$data)

# Run mlm
acts <- run_mlm(mat=mat, net=net, .source='source', .target='target',
                .mor='weight', minsize = 5)
acts

# Extract mlm and store it in pathwaysmlm in data
data[['pathwaysmlm']] <- acts %>%
  pivot_wider(id_cols = 'source', names_from = 'condition',
              values_from = 'score') %>%
  column_to_rownames('source') %>%
  Seurat::CreateAssayObject(.)

# Change assay
DefaultAssay(object = data) <- "pathwaysmlm"

# Scale the data
data <- ScaleData(data)
data@assays$pathwaysmlm@data <- data@assays$pathwaysmlm@scale.data

p1 <- DimPlot(data, reduction = "umap", label = TRUE, pt.size = 0.5) + 
  NoLegend() + ggtitle('Cell types')
p2 <- (FeaturePlot(data, features = c("Trail")) & 
         scale_colour_gradient2(low = 'blue', mid = 'white', high = 'red')) +
  ggtitle('Trail activity')
pdf("Graphs/Main/Pathways_Trailactivity_TET_C.pdf", width=15, height=10)
p1 | p2
dev.off()

# Extract activities from object as a long dataframe
df <- t(as.matrix(data@assays$pathwaysmlm@data)) %>%
  as.data.frame() %>%
  mutate(cluster = Idents(data)) %>%
  pivot_longer(cols = -cluster, names_to = "source", values_to = "score") %>%
  group_by(cluster, source) %>%
  summarise(mean = mean(score))

# Transform to wide matrix
top_acts_mat <- df %>%
  pivot_wider(id_cols = 'cluster', names_from = 'source',
              values_from = 'mean') %>%
  column_to_rownames('cluster') %>%
  as.matrix()

write.xlsx(top_acts_mat, "Lists/Pathways_TET_C.xlsx")

#plot mean values
daten<-read_excel("Lists/Pathways_TET_C_mean.xlsx")
rownames(daten) <- daten[[1]]
rn <- rownames(daten)
daten <- daten[-1]
top_acts_mat <- as.matrix(daten)
rownames(top_acts_mat) <- rn

# Choose color palette
palette_length = 100
my_color = colorRampPalette(c("Darkblue", "white","red"))(palette_length)

my_breaks <- c(seq(-2, 0, length.out=ceiling(palette_length/2) + 1),
               seq(0.05, 2, length.out=floor(palette_length/2)))

top_acts_mat_adult_TET_C <- top_acts_mat

# Plot
pdf("Graphs/Main/Pathways_Heatmap_TET_C.pdf")
pheatmap(top_acts_mat, border_color = NA, color=my_color, breaks = my_breaks) 
dev.off()





#(57)_pw_TET_MNT####
TET_MNT_pathway <- subset(th.total, idents = "TET_MNT")
Idents(TET_MNT_pathway)<- TET_MNT_pathway$celltype.tissue.donor
table(TET_MNT_pathway$celltype.tissue.donor)
net <- get_progeny(organism = 'human', top = 500)
net

data <- TET_MNT_pathway
# Extract the normalized log-transformed counts
mat <- as.matrix(data@assays$RNA$data)

# Run mlm
acts <- run_mlm(mat=mat, net=net, .source='source', .target='target',
                .mor='weight', minsize = 5)
acts

# Extract mlm and store it in pathwaysmlm in data
data[['pathwaysmlm']] <- acts %>%
  pivot_wider(id_cols = 'source', names_from = 'condition',
              values_from = 'score') %>%
  column_to_rownames('source') %>%
  Seurat::CreateAssayObject(.)

# Change assay
DefaultAssay(object = data) <- "pathwaysmlm"

# Scale the data
data <- ScaleData(data)
data@assays$pathwaysmlm@data <- data@assays$pathwaysmlm@scale.data

p1 <- DimPlot(data, reduction = "umap", label = TRUE, pt.size = 0.5) + 
  NoLegend() + ggtitle('Cell types')
p2 <- (FeaturePlot(data, features = c("Trail")) & 
         scale_colour_gradient2(low = 'blue', mid = 'white', high = 'red')) +
  ggtitle('Trail activity')
pdf("Graphs/Main/Pathways_Trailactivity_TET_MNT.pdf", width=15, height=10)
p1 | p2
dev.off()

# Extract activities from object as a long dataframe
df <- t(as.matrix(data@assays$pathwaysmlm@data)) %>%
  as.data.frame() %>%
  mutate(cluster = Idents(data)) %>%
  pivot_longer(cols = -cluster, names_to = "source", values_to = "score") %>%
  group_by(cluster, source) %>%
  summarise(mean = mean(score))

# Transform to wide matrix
top_acts_mat <- df %>%
  pivot_wider(id_cols = 'cluster', names_from = 'source',
              values_from = 'mean') %>%
  column_to_rownames('cluster') %>%
  as.matrix()

write.xlsx(top_acts_mat, "Lists/Pathways_TET_MNT.xlsx")

#plot mean values
daten<-read_excel("Lists/Pathways_TET_MNT_mean.xlsx")
rownames(daten) <- daten[[1]]
rn <- rownames(daten)
daten <- daten[-1]
top_acts_mat <- as.matrix(daten)
rownames(top_acts_mat) <- rn

# Choose color palette
palette_length = 100
my_color = colorRampPalette(c("Darkblue", "white","red"))(palette_length)

my_breaks <- c(seq(-2, 0, length.out=ceiling(palette_length/2) + 1),
               seq(0.05, 2, length.out=floor(palette_length/2)))

top_acts_mat_adult_TET_MNT <- top_acts_mat

# Plot
pdf("Graphs/Main/Pathways_Heatmap_TET_MNT.pdf")
pheatmap(top_acts_mat, border_color = NA, color=my_color, breaks = my_breaks) 
dev.off()










#(58)_APP-CD74 Interaction TET_A ####

Idents(TC_sub)<- TC_sub$tissue
Idents(BC_sub)<- BC_sub$tissue
Idents(PC_sub)<- PC_sub$tissue
Idents(TEC_sub)<- TEC_sub$tissue
Idents(DC_sub)<- DC_sub$tissue
Idents(MAC_Mono_sub)<- MAC_Mono_sub$tissue
Idents(EC_sub)<- EC_sub$tissue
Idents(FB_sub)<- FB_sub$tissue
Idents(VSMC_Peri_sub)<- VSMC_Peri_sub$tissue

TC_sub_TET_A <- subset(TC_sub, idents = "TET_A")
BC_sub_TET_A <- subset(BC_sub, idents = "TET_A")
PC_sub_TET_A <- subset(PC_sub, idents = "TET_A")
TEC_sub_TET_A <- subset(TEC_sub, idents = "TET_A")
DC_sub_TET_A <- subset(DC_sub, idents = "TET_A")
MAC_Mono_sub_TET_A <- subset(MAC_Mono_sub, idents = "TET_A")
EC_sub_TET_A <- subset(EC_sub, idents = "TET_A")
FB_sub_TET_A <- subset(FB_sub, idents = "TET_A")
VSMC_Peri_sub_TET_A <- subset(VSMC_Peri_sub, idents = "TET_A")

TET_A_sub <- merge( TC_sub_TET_A, y= c(BC_sub_TET_A,PC_sub_TET_A,TEC_sub_TET_A,DC_sub_TET_A,MAC_Mono_sub_TET_A,EC_sub_TET_A,FB_sub_TET_A,VSMC_Peri_sub_TET_A))

pdf("Graphs/Main/Dotplot_TET_A_APP_CD74.pdf", width = 10, height = 5)
DotPlot(TET_A_sub, features=c("APP","CD74"), group.by = "celltype.full",assay="RNA") + coord_flip()+ scale_color_viridis(option = "H")+theme(axis.text.x = element_text(angle = 45)) + ggtitle("TEC Cluster_DEG") + scale_y_discrete(guide = guide_axis(angle = 90))
dev.off()

#Heatmap#
## Define Gene lists that you want to plot as Heatmap

Idents(TET_A_sub)<- TET_A_sub$celltype.full
genes <- as.character(c("APP","CD74"))
avgexp <- TET_A_sub

## Plot Heatmaps
## Calculate Average gene expressions (pseudobulked) for each of the celltypes in each condition
levels(avgexp)
DefaultAssay(avgexp) <- "RNA"
avgexp <- NormalizeData(avgexp)
avgexp <- AggregateExpression(avgexp, return.seurat =T, assays = "RNA")


# Modify the heatmap plot
library(RColorBrewer)
heatmap_plot <- DoHeatmap(
  avgexp,
  features = genes,
  size = 6,
  group.bar.height = 0.02,
  draw.lines = FALSE,angle = 90,
  group.colors = color_levelTEC
) +
  scale_fill_gradientn(colours = rev(brewer.pal(8, name = "RdYlBu")))

# Modify the theme to make rownames cursive
heatmap_plot 
ggsave("Graphs/Main/Heatmap_TET_A_APP_CD74.pdf", width = 40, height = 20, units = c("cm"))


pdf("Graphs/Main/Featureblend_TET_A_APP_CD74.pdf", width=10, height=5)
FeaturePlot(TET_A_pathway,c("APP","CD74"), blend=T, order=T, label=T, cols= c("lightgrey","green","blue"), max.cutoff = 2)
dev.off()

save(TET_A_sub,file="TET_A_sub.RData")

