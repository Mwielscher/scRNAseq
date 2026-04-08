#this script accompanies the publication "The transcriptional profile of keloidal Schwann cells"
#code was written by Dr. Martin Direder
#martin.direder@meduniwien.ac.at

setwd(dir ="..")

#load packages and define color schemes
{
  library(Polychrome)
  library(Seurat)
  library(dplyr)
  library(magrittr)
  library(ggplot2)
  library(xlsx)
  library(patchwork)
  library(sctransform)
  library(monocle3)
  library(SeuratWrappers)
  library(ggrepel)
  library(tidyr)
  library(clustree)
  library(ggpubr)
  library(ggsignif)
  library(limma)
  library(Matrix.utils)
  library(rstatix)
  library(CellChat)
  library(enrichplot)
  library(org.Hs.eg.db)
  library(clusterProfiler)
  
  color_samples <- as.character(c("#ff8c00","#e81123","#68217a","#00188f","#ec008c","#00bcf2","#00b294","#009e49","#bad80a","#fff100"))
  single_col <- as.character("#02025f")
  color_DEG <- as.character(c("#54990F","#C51D34"))
    # sXL#009fdd, kXL#0f6bb3, sTT#891f05, kCCD#6eb440, kMD#f6911f
  color_datasets <- as.character(c("#891f05","#009fdd","#0f6bb3","#6eb440","#f6911f"))
  color_datasets_Keloid_only <- as.character(c("#f6911f","#6eb440","#0f6bb3","#009fdd","#000000"))
  color_datasets_order1 <- as.character(c("#009fdd","#0f6bb3","#891f05","#6eb440","#f6911f","#009fdd","#0f6bb3","#0f6bb3","#009fdd","#009fdd","#6eb440","#f6911f","#f6911f","#f6911f","#891f05","#f6911f"))
    # FB#990F26, SMC#B33E52, PC#CC7A88, KC#99600F,  EC#54990F, LEC#78B33E, TC#45b6fe, BC#3792cb, MAC/DC#296d98, DC#, SC#F9A602, MEL#967ACC, SGC#666666
  color_UMAP_all <- as.character(c("#990F26","#CC7A88","#CC7A88","#CC7A88","#99600F","#54990F","#78B33E","#45b6fe","#3792cb","#296d98","#296d98","#296d98","#F9A602","#967ACC","#666666"))
  color_UMAP_sTT <- as.character(c("#990F26","#B33E52","#CC7A88","#99600F","#54990F","#45b6fe","#296d98","#F9A602","#967ACC","#666666"))
  color_UMAP_sXL <- as.character(c("#990F26","#B33E52","#CC7A88","#99600F","#54990F","#78B33E","#45b6fe","#296d98","#F9A602","#666666"))
  color_UMAP_nCCD <- as.character(c("#990F26","#CC7A88","#99600F","#54990F","#78B33E","#45b6fe","#296d98"))
  color_UMAP_nMD <- as.character(c("#990F26","#CC7A88","#99600F","#54990F","#78B33E","#45b6fe","#296d98","#967ACC"))
  color_UMAP_kMD <- as.character(c("#990F26","#CC7A88","#99600F","#54990F","#78B33E","#45b6fe","#3792cb","#296d98","#F9A602"))
  color_UMAP_kCCD <- as.character(c("#990F26","#B33E52","#CC7A88","#99600F","#54990F","#78B33E","#296d98","#F9A602","#666666"))
  color_UMAP_kXL <- as.character(c("#990F26","#CC7A88","#99600F","#54990F","#78B33E","#45b6fe","#296d98","#F9A602","#666666"))
  # Myel+Nonmyel#338333, Myel#344d0e, Nonmyel1#68991c, Nonmyel2#8acc26, Keloid#FFBF00, Keloid2#fcae1e Prolif#1F78C8, EC#C8308C, FB#FF83FA
  color_SC_sTT <- as.character(c("#338333","#FF83FA"))
  color_SC_sXL <- as.character(c("#344d0e","#68991c","#FFBF00"))
  color_SC_kMD <- as.character(c("#338333","#FFBF00","#1F78C8","#C8308C","#FF83FA"))
  color_SC_kCCD <- as.character(c("#338333","#FFBF00"))
  color_SC_kXL <- as.character(c("#344d0e","#68991c","#8acc26","#FFBF00"))
  color_SC_all <- as.character(c("#344d0e","#68991c","#FFBF00","#1F78C8","#C8308C","#FF83FA"))
  color_SC_basal <- as.character(c("#344d0e","#338333","#68991c","#68991c","#8acc26","#FFBF00","#1F78C8","#C8308C","#FF83FA"))
  color_SC_all_celltype <- as.character(c("#344d0e", "#344d0e", "#338333", "#338333", "#338333", "#68991c", "#68991c", "#8acc26" ,"#FFBF00","#FFBF00","#FFBF00","#FFBF00","#1F78C8","#C8308C","#FF83FA","#FF83FA"))
  color_SC_all_celltype_order<-as.character(c("#338333","#FF83FA","#344d0e","#68991c","#FFBF00","#344d0e","#68991c","#8acc26","#FFBF00","#338333","#FFBF00","#338333","#FFBF00","#1F78C8","#C8308C","#FF83FA"))
  color_SC_all_donor <- as.character(c("#338333","#FF83FA","#344d0e","#68991c","#FFBF00","#344d0e","#68991c","#8acc26","#FFBF00","#338333","#FFBF00","#338333","#FFBF00","#1F78C8","#C8308C","#FF83FA"))
  color_SC_k.i.all_donor <- as.character(c("#344d0e","#68991c","#8acc26","#FFBF00","#fcae1e","#1F78C8","#C8308c","#FF83FA"))
  color_SC_repvsmat <- as.character(c("#338333","#338333","#FFBF00","#338333","#FFBF00","#338333","#FFBF00","#338333","#FFBF00"))
  color_cellchat_keloid <- as.character(c("#338333","#FFBF00","#1F78C8","#C8308C","#FF83FA","#990F26","#CC7A88","#99600F","#54990F","#78B33E","#45b6fe","#296d98","#666666"))
  color_cellchat_keloid_mix <- as.character(c("#990F26","#CC7A88","#99600F","#54990F","#78B33E","#45b6fe","#296d98","#FFBF00","#666666"))
  color_cellchat_branch <- as.character(c("#990F26","#CC7A88","#99600F","#54990F","#78B33E","#45b6fe","#296d98","#ff6600","#666666"))
  color_cellchat_skin <- as.character(c("#338333","#990F26","#CC7A88","#99600F","#54990F","#45b6fe","#296d98","#967ACC","#666666"))
  color_cellchat_skin_mix <- as.character(c("#990F26","#CC7A88","#99600F","#54990F","#45b6fe","#296d98","#FFBF00","#967ACC","#666666"))
  color_branch <- as.character(c("#7f7f7f","#ff6600","#7f7f7f","#7f7f7f","#7f7f7f","#7f7f7f","#7f7f7f","#7f7f7f"))
}

#prepare data

#############################Skin-Data wrangling############################
###-###-###-###-###-###-###-###-###-###-###-###-###

#load tabib_skin.data
tabib.data<-read.csv(file = "../Skin_6Control_rawUMI.csv")
rownames(tabib.data)<-tabib.data$X
tabib.data<-tabib.data[,-1]

tabib.metadata<-read.csv(file = "../Skin_6Control_Metadata.csv")  
rownames(tabib.metadata)<-tabib.metadata$X
tabib.metadata<-tabib.metadata[,-1]


#setup Tabib data
tabib.skin<-CreateSeuratObject(counts=tabib.data, assay="RNA", meta.data = tabib.metadata)

#subset samples
skin1.sct <- subset(tabib.skin, ident=c("SC18control"))
skin2.sct <- subset(tabib.skin, ident=c("SC1control"))
skin3.sct <- subset(tabib.skin, ident=c("SC32control"))
skin4.sct <- subset(tabib.skin, ident=c("SC33control"))
skin5.sct <- subset(tabib.skin, ident=c("SC34control"))
skin6.sct <- subset(tabib.skin, ident=c("SC4control"))

#Get original 10x Data, convert to NCBI Symbol, Create Seurat Object
#skin1
skin1.sct.data <- GetAssay(skin1.sct, assay="RNA")
skin1.sct.data<-skin1.sct.data@counts
Signlist1=alias2SymbolUsingNCBI(as.character(rownames(skin1.sct.data)),paste0(c("../Homo_sapiens.gene_info")),required.columns = c("GeneID","Symbol"))
Signlist2=as.data.frame(cbind(Signlist1,as.character(rownames(skin1.sct.data))))
colnames(Signlist2)[3]=c("rname")
Signlist2$new = ifelse(is.na(Signlist2$Symbol), as.character(Signlist2$rname), Signlist2$Symbol)
rownames(skin1.sct.data)= as.character(Signlist2$new)


#skin2
skin2.sct.data <- GetAssay(skin2.sct, assay="RNA")
skin2.sct.data<-skin2.sct.data@counts
Signlist1=alias2SymbolUsingNCBI(as.character(rownames(skin2.sct.data)),paste0(c("../Homo_sapiens.gene_info")),required.columns = c("GeneID","Symbol"))
Signlist2=as.data.frame(cbind(Signlist1,as.character(rownames(skin2.sct.data))))
colnames(Signlist2)[3]=c("rname")
Signlist2$new = ifelse(is.na(Signlist2$Symbol), as.character(Signlist2$rname), Signlist2$Symbol)
rownames(skin2.sct.data)= as.character(Signlist2$new)


#skin3
skin3.sct.data <- GetAssay(skin3.sct, assay="RNA")
skin3.sct.data<-skin3.sct.data@counts
Signlist1=alias2SymbolUsingNCBI(as.character(rownames(skin3.sct.data)),paste0(c("../Homo_sapiens.gene_info")),required.columns = c("GeneID","Symbol"))
Signlist2=as.data.frame(cbind(Signlist1,as.character(rownames(skin3.sct.data))))
colnames(Signlist2)[3]=c("rname")
Signlist2$new = ifelse(is.na(Signlist2$Symbol), as.character(Signlist2$rname), Signlist2$Symbol)
rownames(skin3.sct.data)= as.character(Signlist2$new)



#skin4
skin4.sct.data <- GetAssay(skin4.sct, assay="RNA")
skin4.sct.data<-skin4.sct.data@counts
Signlist1=alias2SymbolUsingNCBI(as.character(rownames(skin4.sct.data)),paste0(c("../Homo_sapiens.gene_info")),required.columns = c("GeneID","Symbol"))
Signlist2=as.data.frame(cbind(Signlist1,as.character(rownames(skin4.sct.data))))
colnames(Signlist2)[3]=c("rname")
Signlist2$new = ifelse(is.na(Signlist2$Symbol), as.character(Signlist2$rname), Signlist2$Symbol)
rownames(skin4.sct.data)= as.character(Signlist2$new)


#skin5
skin5.sct.data <- GetAssay(skin5.sct, assay="RNA")
skin5.sct.data<-skin5.sct.data@counts
Signlist1=alias2SymbolUsingNCBI(as.character(rownames(skin5.sct.data)),paste0(c("../Homo_sapiens.gene_info")),required.columns = c("GeneID","Symbol"))
Signlist2=as.data.frame(cbind(Signlist1,as.character(rownames(skin5.sct.data))))
colnames(Signlist2)[3]=c("rname")
Signlist2$new = ifelse(is.na(Signlist2$Symbol), as.character(Signlist2$rname), Signlist2$Symbol)
rownames(skin5.sct.data)= as.character(Signlist2$new)


#skin6
skin6.sct.data <- GetAssay(skin6.sct, assay="RNA")
skin6.sct.data<-skin6.sct.data@counts
Signlist1=alias2SymbolUsingNCBI(as.character(rownames(skin6.sct.data)),paste0(c("../Homo_sapiens.gene_info")),required.columns = c("GeneID","Symbol"))
Signlist2=as.data.frame(cbind(Signlist1,as.character(rownames(skin6.sct.data))))
colnames(Signlist2)[3]=c("rname")
Signlist2$new = ifelse(is.na(Signlist2$Symbol), as.character(Signlist2$rname), Signlist2$Symbol)
rownames(skin6.sct.data)= as.character(Signlist2$new)


# load close skin
skincl1.sct.data<-Read10X("../K007 - skin/filtered data")
skincl2.sct.data<-Read10X("../K009 - skin/filtered data")
skincl3.sct.data<-Read10X("../K013 - skin/filtered data")
skincl4.sct.data<-Read10X("../K012 - skin/filtered data")

#convert to NCBI Symbol
#skincl1
Signlist1=alias2SymbolUsingNCBI(as.character(rownames(skincl1.sct.data)),paste0(c("../Homo_sapiens.gene_info")),required.columns = c("GeneID","Symbol"))
Signlist2=as.data.frame(cbind(Signlist1,as.character(rownames(skincl1.sct.data))))
colnames(Signlist2)[3]=c("rname")
Signlist2$new = ifelse(is.na(Signlist2$Symbol), as.character(Signlist2$rname), Signlist2$Symbol)
rownames(skincl1.sct.data)= as.character(Signlist2$new)

#skincl2
Signlist1=alias2SymbolUsingNCBI(as.character(rownames(skincl2.sct.data)),paste0(c("..D/Homo_sapiens.gene_info")),required.columns = c("GeneID","Symbol"))
Signlist2=as.data.frame(cbind(Signlist1,as.character(rownames(skincl2.sct.data))))
colnames(Signlist2)[3]=c("rname")
Signlist2$new = ifelse(is.na(Signlist2$Symbol), as.character(Signlist2$rname), Signlist2$Symbol)
rownames(skincl2.sct.data)= as.character(Signlist2$new)

#skincl3
Signlist1=alias2SymbolUsingNCBI(as.character(rownames(skincl3.sct.data)),paste0(c("../Homo_sapiens.gene_info")),required.columns = c("GeneID","Symbol"))
Signlist2=as.data.frame(cbind(Signlist1,as.character(rownames(skincl3.sct.data))))
colnames(Signlist2)[3]=c("rname")
Signlist2$new = ifelse(is.na(Signlist2$Symbol), as.character(Signlist2$rname), Signlist2$Symbol)
rownames(skincl3.sct.data)= as.character(Signlist2$new)

#skincl4
Signlist1=alias2SymbolUsingNCBI(as.character(rownames(skincl4.sct.data)),paste0(c("../Homo_sapiens.gene_info")),required.columns = c("GeneID","Symbol"))
Signlist2=as.data.frame(cbind(Signlist1,as.character(rownames(skincl4.sct.data))))
colnames(Signlist2)[3]=c("rname")
Signlist2$new = ifelse(is.na(Signlist2$Symbol), as.character(Signlist2$rname), Signlist2$Symbol)
rownames(skincl4.sct.data)= as.character(Signlist2$new)

# load normscar
nscar1.sct.data<-Read10X("../scar_1/filtered feature bc matrix")
nscar2.sct.data<-Read10X("../scar_2/filtered feature bc matrix")
nscar3.sct.data<-Read10X("../scar_3/filtered feature bc matrix")
nscar4.sct.data<-Read10X("../10x_PMID34140509_normalscar/NF1_matrix")
nscar5.sct.data<-Read10X("../10x_PMID34140509_normalscar/NF2_matrix")
nscar6.sct.data<-Read10X("../10x_PMID34140509_normalscar/NF3_matrix")

#convert to NCBI Symbol
#nscar1
Signlist1=alias2SymbolUsingNCBI(as.character(rownames(nscar1.sct.data)),paste0(c("../Homo_sapiens.gene_info")),required.columns = c("GeneID","Symbol"))
Signlist2=as.data.frame(cbind(Signlist1,as.character(rownames(nscar1.sct.data))))
colnames(Signlist2)[3]=c("rname")
Signlist2$new = ifelse(is.na(Signlist2$Symbol), as.character(Signlist2$rname), Signlist2$Symbol)
rownames(nscar1.sct.data)= as.character(Signlist2$new)

#nscar2
Signlist1=alias2SymbolUsingNCBI(as.character(rownames(nscar2.sct.data)),paste0(c("../Homo_sapiens.gene_info")),required.columns = c("GeneID","Symbol"))
Signlist2=as.data.frame(cbind(Signlist1,as.character(rownames(nscar2.sct.data))))
colnames(Signlist2)[3]=c("rname")
Signlist2$new = ifelse(is.na(Signlist2$Symbol), as.character(Signlist2$rname), Signlist2$Symbol)
rownames(nscar2.sct.data)= as.character(Signlist2$new)

#nscar3
Signlist1=alias2SymbolUsingNCBI(as.character(rownames(nscar3.sct.data)),paste0(c("../Homo_sapiens.gene_info")),required.columns = c("GeneID","Symbol"))
Signlist2=as.data.frame(cbind(Signlist1,as.character(rownames(nscar3.sct.data))))
colnames(Signlist2)[3]=c("rname")
Signlist2$new = ifelse(is.na(Signlist2$Symbol), as.character(Signlist2$rname), Signlist2$Symbol)
rownames(nscar3.sct.data)= as.character(Signlist2$new)

#nscar4
Signlist1=alias2SymbolUsingNCBI(as.character(rownames(nscar4.sct.data)),paste0(c("../Homo_sapiens.gene_info")),required.columns = c("GeneID","Symbol"))
Signlist2=as.data.frame(cbind(Signlist1,as.character(rownames(nscar4.sct.data))))
colnames(Signlist2)[3]=c("rname")
Signlist2$new = ifelse(is.na(Signlist2$Symbol), as.character(Signlist2$rname), Signlist2$Symbol)
rownames(nscar4.sct.data)= as.character(Signlist2$new)

#nscar5
Signlist1=alias2SymbolUsingNCBI(as.character(rownames(nscar5.sct.data)),paste0(c("../Homo_sapiens.gene_info")),required.columns = c("GeneID","Symbol"))
Signlist2=as.data.frame(cbind(Signlist1,as.character(rownames(nscar5.sct.data))))
colnames(Signlist2)[3]=c("rname")
Signlist2$new = ifelse(is.na(Signlist2$Symbol), as.character(Signlist2$rname), Signlist2$Symbol)
rownames(nscar5.sct.data)= as.character(Signlist2$new)

#nscar6
Signlist1=alias2SymbolUsingNCBI(as.character(rownames(nscar6.sct.data)),paste0(c("../Homo_sapiens.gene_info")),required.columns = c("GeneID","Symbol"))
Signlist2=as.data.frame(cbind(Signlist1,as.character(rownames(nscar6.sct.data))))
colnames(Signlist2)[3]=c("rname")
Signlist2$new = ifelse(is.na(Signlist2$Symbol), as.character(Signlist2$rname), Signlist2$Symbol)
rownames(nscar6.sct.data)= as.character(Signlist2$new)


# load keloid
keloid1.sct.data<-Read10X("../keloid_1/filtered feature bc matrix")
keloid2.sct.data<-Read10X("../keloid_2/filtered feature bc matrix")
keloid3L.sct.data<-Read10X("../keloid_3L/filtered feature bc matrix")
keloid3R.sct.data<-Read10X("../keloid_3R/filtered feature bc matrix")
keloid4.sct.data<-Read10X("../10x_PMID34140509_Keloid/KF1_matrix")
keloid5.sct.data<-Read10X("../10x_PMID34140509_Keloid/KF2_matrix")
keloid6.sct.data<-Read10X("../10x_PMID34140509_Keloid/KF3_matrix")
keloid7.sct.data<-Read10X("../PMID34242659_Keloid/K007 - keloid/filtered data")
keloid8.sct.data<-Read10X("../PMID34242659_Keloid/K009 - keloid/filtered data")
keloid9.sct.data<-Read10X("../PMID34242659_Keloid/K013 - keloid/filtered data")
keloid10.sct.data<-Read10X("../PMID34242659_Keloid/K012 - keloid/filtered data")



#convert to NCBI Symbol
#keloid1
Signlist1=alias2SymbolUsingNCBI(as.character(rownames(keloid1.sct.data)),paste0(c("../Homo_sapiens.gene_info")),required.columns = c("GeneID","Symbol"))
Signlist2=as.data.frame(cbind(Signlist1,as.character(rownames(keloid1.sct.data))))
colnames(Signlist2)[3]=c("rname")
Signlist2$new = ifelse(is.na(Signlist2$Symbol), as.character(Signlist2$rname), Signlist2$Symbol)
rownames(keloid1.sct.data)= as.character(Signlist2$new)

#keloid2
Signlist1=alias2SymbolUsingNCBI(as.character(rownames(keloid2.sct.data)),paste0(c("../Homo_sapiens.gene_info")),required.columns = c("GeneID","Symbol"))
Signlist2=as.data.frame(cbind(Signlist1,as.character(rownames(keloid2.sct.data))))
colnames(Signlist2)[3]=c("rname")
Signlist2$new = ifelse(is.na(Signlist2$Symbol), as.character(Signlist2$rname), Signlist2$Symbol)
rownames(keloid2.sct.data)= as.character(Signlist2$new)

#keloid3L
Signlist1=alias2SymbolUsingNCBI(as.character(rownames(keloid3L.sct.data)),paste0(c("../Homo_sapiens.gene_info")),required.columns = c("GeneID","Symbol"))
Signlist2=as.data.frame(cbind(Signlist1,as.character(rownames(keloid3L.sct.data))))
colnames(Signlist2)[3]=c("rname")
Signlist2$new = ifelse(is.na(Signlist2$Symbol), as.character(Signlist2$rname), Signlist2$Symbol)
rownames(keloid3L.sct.data)= as.character(Signlist2$new)

#keloid3R
Signlist1=alias2SymbolUsingNCBI(as.character(rownames(keloid3R.sct.data)),paste0(c("../Homo_sapiens.gene_info")),required.columns = c("GeneID","Symbol"))
Signlist2=as.data.frame(cbind(Signlist1,as.character(rownames(keloid3R.sct.data))))
colnames(Signlist2)[3]=c("rname")
Signlist2$new = ifelse(is.na(Signlist2$Symbol), as.character(Signlist2$rname), Signlist2$Symbol)
rownames(keloid3R.sct.data)= as.character(Signlist2$new)

#keloid4
Signlist1=alias2SymbolUsingNCBI(as.character(rownames(keloid4.sct.data)),paste0(c("../Homo_sapiens.gene_info")),required.columns = c("GeneID","Symbol"))
Signlist2=as.data.frame(cbind(Signlist1,as.character(rownames(keloid4.sct.data))))
colnames(Signlist2)[3]=c("rname")
Signlist2$new = ifelse(is.na(Signlist2$Symbol), as.character(Signlist2$rname), Signlist2$Symbol)
rownames(keloid4.sct.data)= as.character(Signlist2$new)

#keloid5
Signlist1=alias2SymbolUsingNCBI(as.character(rownames(keloid5.sct.data)),paste0(c("../Homo_sapiens.gene_info")),required.columns = c("GeneID","Symbol"))
Signlist2=as.data.frame(cbind(Signlist1,as.character(rownames(keloid5.sct.data))))
colnames(Signlist2)[3]=c("rname")
Signlist2$new = ifelse(is.na(Signlist2$Symbol), as.character(Signlist2$rname), Signlist2$Symbol)
rownames(keloid5.sct.data)= as.character(Signlist2$new)

#keloid6
Signlist1=alias2SymbolUsingNCBI(as.character(rownames(keloid6.sct.data)),paste0(c("../Homo_sapiens.gene_info")),required.columns = c("GeneID","Symbol"))
Signlist2=as.data.frame(cbind(Signlist1,as.character(rownames(keloid6.sct.data))))
colnames(Signlist2)[3]=c("rname")
Signlist2$new = ifelse(is.na(Signlist2$Symbol), as.character(Signlist2$rname), Signlist2$Symbol)
rownames(keloid6.sct.data)= as.character(Signlist2$new)

#keloid7
Signlist1=alias2SymbolUsingNCBI(as.character(rownames(keloid7.sct.data)),paste0(c("../Homo_sapiens.gene_info")),required.columns = c("GeneID","Symbol"))
Signlist2=as.data.frame(cbind(Signlist1,as.character(rownames(keloid7.sct.data))))
colnames(Signlist2)[3]=c("rname")
Signlist2$new = ifelse(is.na(Signlist2$Symbol), as.character(Signlist2$rname), Signlist2$Symbol)
rownames(keloid7.sct.data)= as.character(Signlist2$new)

#keloid8
Signlist1=alias2SymbolUsingNCBI(as.character(rownames(keloid8.sct.data)),paste0(c("../Homo_sapiens.gene_info")),required.columns = c("GeneID","Symbol"))
Signlist2=as.data.frame(cbind(Signlist1,as.character(rownames(keloid8.sct.data))))
colnames(Signlist2)[3]=c("rname")
Signlist2$new = ifelse(is.na(Signlist2$Symbol), as.character(Signlist2$rname), Signlist2$Symbol)
rownames(keloid8.sct.data)= as.character(Signlist2$new)

#keloid9
Signlist1=alias2SymbolUsingNCBI(as.character(rownames(keloid9.sct.data)),paste0(c("../Homo_sapiens.gene_info")),required.columns = c("GeneID","Symbol"))
Signlist2=as.data.frame(cbind(Signlist1,as.character(rownames(keloid9.sct.data))))
colnames(Signlist2)[3]=c("rname")
Signlist2$new = ifelse(is.na(Signlist2$Symbol), as.character(Signlist2$rname), Signlist2$Symbol)
rownames(keloid9.sct.data)= as.character(Signlist2$new)

#keloid10
Signlist1=alias2SymbolUsingNCBI(as.character(rownames(keloid10.sct.data)),paste0(c("../Homo_sapiens.gene_info")),required.columns = c("GeneID","Symbol"))
Signlist2=as.data.frame(cbind(Signlist1,as.character(rownames(keloid10.sct.data))))
colnames(Signlist2)[3]=c("rname")
Signlist2$new = ifelse(is.na(Signlist2$Symbol), as.character(Signlist2$rname), Signlist2$Symbol)
rownames(keloid10.sct.data)= as.character(Signlist2$new)


# keep only Genes detected in all dataset Template
#dimension_check
dim(skin1.sct.data)
dim(keloid1.sct.data)
dim(skincl1.sct.data)
dim(keloid8.sct.data)
dim(keloid4.sct.data)

ROWLIST1 <- skin1.sct.data[!duplicated(rownames(skin1.sct.data)),]
ROWLIST2 <- keloid1.sct.data[!duplicated(rownames(keloid1.sct.data)),]
ROWLIST3 <- skincl1.sct.data[!duplicated(rownames(skincl1.sct.data)),]


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


#remove duplicates
skin1.sct.data <- skin1.sct.data[!duplicated(rownames(skin1.sct.data)),]
skin2.sct.data <- skin2.sct.data[!duplicated(rownames(skin2.sct.data)),]
skin3.sct.data <- skin3.sct.data[!duplicated(rownames(skin3.sct.data)),]
skin4.sct.data <- skin4.sct.data[!duplicated(rownames(skin4.sct.data)),]
skin5.sct.data <- skin5.sct.data[!duplicated(rownames(skin5.sct.data)),]
skin6.sct.data <- skin6.sct.data[!duplicated(rownames(skin6.sct.data)),]
skincl1.sct.data <- skincl1.sct.data[!duplicated(rownames(skincl1.sct.data)),]
skincl2.sct.data <- skincl2.sct.data[!duplicated(rownames(skincl2.sct.data)),]
skincl3.sct.data <- skincl3.sct.data[!duplicated(rownames(skincl3.sct.data)),]
skincl4.sct.data <- skincl4.sct.data[!duplicated(rownames(skincl4.sct.data)),]
nscar1.sct.data <- nscar1.sct.data[!duplicated(rownames(nscar1.sct.data)),]
nscar2.sct.data <- nscar2.sct.data[!duplicated(rownames(nscar2.sct.data)),]
nscar3.sct.data <- nscar3.sct.data[!duplicated(rownames(nscar3.sct.data)),]
nscar4.sct.data <- nscar4.sct.data[!duplicated(rownames(nscar4.sct.data)),]
nscar5.sct.data <- nscar5.sct.data[!duplicated(rownames(nscar5.sct.data)),]
nscar6.sct.data <- nscar6.sct.data[!duplicated(rownames(nscar6.sct.data)),]
keloid1.sct.data <- keloid1.sct.data[!duplicated(rownames(keloid1.sct.data)),]
keloid2.sct.data <- keloid2.sct.data[!duplicated(rownames(keloid2.sct.data)),]
keloid3L.sct.data <- keloid3L.sct.data[!duplicated(rownames(keloid3L.sct.data)),]
keloid3R.sct.data <- keloid3R.sct.data[!duplicated(rownames(keloid3R.sct.data)),]
keloid4.sct.data <- keloid4.sct.data[!duplicated(rownames(keloid4.sct.data)),]
keloid5.sct.data <- keloid5.sct.data[!duplicated(rownames(keloid5.sct.data)),]
keloid6.sct.data <- keloid6.sct.data[!duplicated(rownames(keloid6.sct.data)),]
keloid7.sct.data <- keloid7.sct.data[!duplicated(rownames(keloid7.sct.data)),]
keloid8.sct.data <- keloid8.sct.data[!duplicated(rownames(keloid8.sct.data)),]
keloid9.sct.data <- keloid9.sct.data[!duplicated(rownames(keloid9.sct.data)),]
keloid10.sct.data <- keloid10.sct.data[!duplicated(rownames(keloid10.sct.data)),]


# adapt all
Temp3<-join.Matrix(Temp2,skin1.sct.data, by.x = rownames(Temp2), by.y=rownames(skin1.sct.data), all.x = F, all.y = F)
skin1.sct.data <- Temp3[,3:ncol(Temp3)]

Temp3<-join.Matrix(Temp2,skin2.sct.data, by.x = rownames(Temp2), by.y=rownames(skin2.sct.data), all.x = F, all.y = F)
skin2.sct.data <- Temp3[,3:ncol(Temp3)]

Temp3<-join.Matrix(Temp2,skin3.sct.data, by.x = rownames(Temp2), by.y=rownames(skin3.sct.data), all.x = F, all.y = F)
skin3.sct.data <- Temp3[,3:ncol(Temp3)]

Temp3<-join.Matrix(Temp2,skin4.sct.data, by.x = rownames(Temp2), by.y=rownames(skin4.sct.data), all.x = F, all.y = F)
skin4.sct.data <- Temp3[,3:ncol(Temp3)]

Temp3<-join.Matrix(Temp2,skin5.sct.data, by.x = rownames(Temp2), by.y=rownames(skin5.sct.data), all.x = F, all.y = F)
skin5.sct.data <- Temp3[,3:ncol(Temp3)]

Temp3<-join.Matrix(Temp2,skin6.sct.data, by.x = rownames(Temp2), by.y=rownames(skin6.sct.data), all.x = F, all.y = F)
skin6.sct.data <- Temp3[,3:ncol(Temp3)]

Temp3<-join.Matrix(Temp2,skincl1.sct.data, by.x = rownames(Temp2), by.y=rownames(skincl1.sct.data), all.x = F, all.y = F)
skincl1.sct.data <- Temp3[,3:ncol(Temp3)]

Temp3<-join.Matrix(Temp2,skincl2.sct.data, by.x = rownames(Temp2), by.y=rownames(skincl2.sct.data), all.x = F, all.y = F)
skincl2.sct.data <- Temp3[,3:ncol(Temp3)]

Temp3<-join.Matrix(Temp2,skincl3.sct.data, by.x = rownames(Temp2), by.y=rownames(skincl3.sct.data), all.x = F, all.y = F)
skincl3.sct.data <- Temp3[,3:ncol(Temp3)]

Temp3<-join.Matrix(Temp2,skincl4.sct.data, by.x = rownames(Temp2), by.y=rownames(skincl4.sct.data), all.x = F, all.y = F)
skincl4.sct.data <- Temp3[,3:ncol(Temp3)]

Temp3<-join.Matrix(Temp2,nscar1.sct.data, by.x = rownames(Temp2), by.y=rownames(nscar1.sct.data), all.x = F, all.y = F)
nscar1.sct.data <- Temp3[,3:ncol(Temp3)]

Temp3<-join.Matrix(Temp2,nscar2.sct.data, by.x = rownames(Temp2), by.y=rownames(nscar2.sct.data), all.x = F, all.y = F)
nscar2.sct.data <- Temp3[,3:ncol(Temp3)]

Temp3<-join.Matrix(Temp2,nscar3.sct.data, by.x = rownames(Temp2), by.y=rownames(nscar3.sct.data), all.x = F, all.y = F)
nscar3.sct.data <- Temp3[,3:ncol(Temp3)]

Temp3<-join.Matrix(Temp2,nscar4.sct.data, by.x = rownames(Temp2), by.y=rownames(nscar4.sct.data), all.x = F, all.y = F)
nscar4.sct.data <- Temp3[,3:ncol(Temp3)]

Temp3<-join.Matrix(Temp2,nscar5.sct.data, by.x = rownames(Temp2), by.y=rownames(nscar5.sct.data), all.x = F, all.y = F)
nscar5.sct.data <- Temp3[,3:ncol(Temp3)]

Temp3<-join.Matrix(Temp2,nscar6.sct.data, by.x = rownames(Temp2), by.y=rownames(nscar6.sct.data), all.x = F, all.y = F)
nscar6.sct.data <- Temp3[,3:ncol(Temp3)]

Temp3<-join.Matrix(Temp2,keloid1.sct.data, by.x = rownames(Temp2), by.y=rownames(keloid1.sct.data), all.x = F, all.y = F)
keloid1.sct.data <- Temp3[,3:ncol(Temp3)]

Temp3<-join.Matrix(Temp2,keloid2.sct.data, by.x = rownames(Temp2), by.y=rownames(keloid2.sct.data), all.x = F, all.y = F)
keloid2.sct.data <- Temp3[,3:ncol(Temp3)]

Temp3<-join.Matrix(Temp2,keloid3L.sct.data, by.x = rownames(Temp2), by.y=rownames(keloid3L.sct.data), all.x = F, all.y = F)
keloid3L.sct.data <- Temp3[,3:ncol(Temp3)]

Temp3<-join.Matrix(Temp2,keloid3R.sct.data, by.x = rownames(Temp2), by.y=rownames(keloid3R.sct.data), all.x = F, all.y = F)
keloid3R.sct.data <- Temp3[,3:ncol(Temp3)]

Temp3<-join.Matrix(Temp2,keloid4.sct.data, by.x = rownames(Temp2), by.y=rownames(keloid4.sct.data), all.x = F, all.y = F)
keloid4.sct.data <- Temp3[,3:ncol(Temp3)]

Temp3<-join.Matrix(Temp2,keloid5.sct.data, by.x = rownames(Temp2), by.y=rownames(keloid5.sct.data), all.x = F, all.y = F)
keloid5.sct.data <- Temp3[,3:ncol(Temp3)]

Temp3<-join.Matrix(Temp2,keloid6.sct.data, by.x = rownames(Temp2), by.y=rownames(keloid6.sct.data), all.x = F, all.y = F)
keloid6.sct.data <- Temp3[,3:ncol(Temp3)]

Temp3<-join.Matrix(Temp2,keloid7.sct.data, by.x = rownames(Temp2), by.y=rownames(keloid7.sct.data), all.x = F, all.y = F)
keloid7.sct.data <- Temp3[,3:ncol(Temp3)]

Temp3<-join.Matrix(Temp2,keloid8.sct.data, by.x = rownames(Temp2), by.y=rownames(keloid8.sct.data), all.x = F, all.y = F)
keloid8.sct.data <- Temp3[,3:ncol(Temp3)]

Temp3<-join.Matrix(Temp2,keloid9.sct.data, by.x = rownames(Temp2), by.y=rownames(keloid9.sct.data), all.x = F, all.y = F)
keloid9.sct.data <- Temp3[,3:ncol(Temp3)]

Temp3<-join.Matrix(Temp2,keloid10.sct.data, by.x = rownames(Temp2), by.y=rownames(keloid10.sct.data), all.x = F, all.y = F)
keloid10.sct.data <- Temp3[,3:ncol(Temp3)]

#Create Sparcematrix
skin1.sct<-CreateSeuratObject(skin1.sct.data)
skin2.sct<-CreateSeuratObject(skin2.sct.data)
skin3.sct<-CreateSeuratObject(skin3.sct.data)
skin4.sct<-CreateSeuratObject(skin4.sct.data)
skin5.sct<-CreateSeuratObject(skin5.sct.data)
skin6.sct<-CreateSeuratObject(skin6.sct.data)
skincl1.sct<-CreateSeuratObject(skincl1.sct.data)
skincl2.sct<-CreateSeuratObject(skincl2.sct.data)
skincl3.sct<-CreateSeuratObject(skincl3.sct.data)
skincl4.sct<-CreateSeuratObject(skincl4.sct.data)
nscar1.sct<-CreateSeuratObject(nscar1.sct.data)
nscar2.sct<-CreateSeuratObject(nscar2.sct.data)
nscar3.sct<-CreateSeuratObject(nscar3.sct.data)
nscar4.sct<-CreateSeuratObject(nscar4.sct.data)
nscar5.sct<-CreateSeuratObject(nscar5.sct.data)
nscar6.sct<-CreateSeuratObject(nscar6.sct.data)
keloid1.sct<-CreateSeuratObject(keloid1.sct.data)
keloid2.sct<-CreateSeuratObject(keloid2.sct.data)
keloid3L.sct<-CreateSeuratObject(keloid3L.sct.data)
keloid3R.sct<-CreateSeuratObject(keloid3R.sct.data)
keloid4.sct<-CreateSeuratObject(keloid4.sct.data)
keloid5.sct<-CreateSeuratObject(keloid5.sct.data)
keloid6.sct<-CreateSeuratObject(keloid6.sct.data)
keloid7.sct<-CreateSeuratObject(keloid7.sct.data)
keloid8.sct<-CreateSeuratObject(keloid8.sct.data)
keloid9.sct<-CreateSeuratObject(keloid9.sct.data)
keloid10.sct<-CreateSeuratObject(keloid10.sct.data)

#flag-1_sample
skin1.sct<-AddMetaData(skin1.sct, "skin_sep1", col.name = "sample")
skin2.sct<-AddMetaData(skin2.sct, "skin_sep2", col.name = "sample")
skin3.sct<-AddMetaData(skin3.sct, "skin_sep3", col.name = "sample")
skin4.sct<-AddMetaData(skin4.sct, "skin_sep4", col.name = "sample")
skin5.sct<-AddMetaData(skin5.sct, "skin_sep5", col.name = "sample")
skin6.sct<-AddMetaData(skin6.sct, "skin_sep6", col.name = "sample")
skincl1.sct<-AddMetaData(skincl1.sct, "skin_adj1", col.name = "sample")
skincl2.sct<-AddMetaData(skincl2.sct, "skin_adj2", col.name = "sample")
skincl3.sct<-AddMetaData(skincl3.sct, "skin_adj3", col.name = "sample")
skincl4.sct<-AddMetaData(skincl4.sct, "skin_adj4", col.name = "sample")
nscar1.sct<-AddMetaData(nscar1.sct, "nscar_1", col.name = "sample")
nscar2.sct<-AddMetaData(nscar2.sct, "nscar_2", col.name = "sample")
nscar3.sct<-AddMetaData(nscar3.sct, "nscar_3", col.name = "sample")
nscar4.sct<-AddMetaData(nscar4.sct, "nscar_4", col.name = "sample")
nscar5.sct<-AddMetaData(nscar5.sct, "nscar_5", col.name = "sample")
nscar6.sct<-AddMetaData(nscar6.sct, "nscar_6", col.name = "sample")
keloid1.sct<-AddMetaData(keloid1.sct, "keloid_centr1", col.name = "sample")
keloid2.sct<-AddMetaData(keloid2.sct, "keloid_centr2", col.name = "sample")
keloid3L.sct<-AddMetaData(keloid3L.sct, "keloid_centr3", col.name = "sample")
keloid3R.sct<-AddMetaData(keloid3R.sct, "keloid_centr4", col.name = "sample")
keloid4.sct<-AddMetaData(keloid4.sct, "keloid_mix1", col.name = "sample")
keloid5.sct<-AddMetaData(keloid5.sct, "keloid_mix2", col.name = "sample")
keloid6.sct<-AddMetaData(keloid6.sct, "keloid_mix3", col.name = "sample")
keloid7.sct<-AddMetaData(keloid7.sct, "keloid_centr5", col.name = "sample")
keloid8.sct<-AddMetaData(keloid8.sct, "keloid_centr6", col.name = "sample")
keloid9.sct<-AddMetaData(keloid9.sct, "keloid_centr7", col.name = "sample")
keloid10.sct<-AddMetaData(keloid10.sct, "keloid_centr8", col.name = "sample")

#flag-2_tissue
skin1.sct<-AddMetaData(skin1.sct, "skin", col.name = "tissue")
skin2.sct<-AddMetaData(skin2.sct, "skin", col.name = "tissue")
skin3.sct<-AddMetaData(skin3.sct, "skin", col.name = "tissue")
skin4.sct<-AddMetaData(skin4.sct, "skin", col.name = "tissue")
skin5.sct<-AddMetaData(skin5.sct, "skin", col.name = "tissue")
skin6.sct<-AddMetaData(skin6.sct, "skin", col.name = "tissue")
skincl1.sct<-AddMetaData(skincl1.sct, "skin", col.name = "tissue")
skincl2.sct<-AddMetaData(skincl2.sct, "skin", col.name = "tissue")
skincl3.sct<-AddMetaData(skincl3.sct, "skin", col.name = "tissue")
skincl4.sct<-AddMetaData(skincl4.sct, "skin", col.name = "tissue")
nscar1.sct<-AddMetaData(nscar1.sct, "nscar", col.name = "tissue")
nscar2.sct<-AddMetaData(nscar2.sct, "nscar", col.name = "tissue")
nscar3.sct<-AddMetaData(nscar3.sct, "nscar", col.name = "tissue")
nscar4.sct<-AddMetaData(nscar4.sct, "nscar", col.name = "tissue")
nscar5.sct<-AddMetaData(nscar5.sct, "nscar", col.name = "tissue")
nscar6.sct<-AddMetaData(nscar6.sct, "nscar", col.name = "tissue")
keloid1.sct<-AddMetaData(keloid1.sct, "keloid", col.name = "tissue")
keloid2.sct<-AddMetaData(keloid2.sct, "keloid", col.name = "tissue")
keloid3L.sct<-AddMetaData(keloid3L.sct, "keloid", col.name = "tissue")
keloid3R.sct<-AddMetaData(keloid3R.sct, "keloid", col.name = "tissue")
keloid4.sct<-AddMetaData(keloid4.sct, "keloid", col.name = "tissue")
keloid5.sct<-AddMetaData(keloid5.sct, "keloid", col.name = "tissue")
keloid6.sct<-AddMetaData(keloid6.sct, "keloid", col.name = "tissue")
keloid7.sct<-AddMetaData(keloid7.sct, "keloid", col.name = "tissue")
keloid8.sct<-AddMetaData(keloid8.sct, "keloid", col.name = "tissue")
keloid9.sct<-AddMetaData(keloid9.sct, "keloid", col.name = "tissue")
keloid10.sct<-AddMetaData(keloid10.sct, "keloid", col.name = "tissue")

#flag-3_origin
skin1.sct<-AddMetaData(skin1.sct, "TT", col.name = "origin")
skin2.sct<-AddMetaData(skin2.sct, "TT", col.name = "origin")
skin3.sct<-AddMetaData(skin3.sct, "TT", col.name = "origin")
skin4.sct<-AddMetaData(skin4.sct, "TT", col.name = "origin")
skin5.sct<-AddMetaData(skin5.sct, "TT", col.name = "origin")
skin6.sct<-AddMetaData(skin6.sct, "TT", col.name = "origin")
skincl1.sct<-AddMetaData(skincl1.sct, "XL", col.name = "origin")
skincl2.sct<-AddMetaData(skincl2.sct, "XL", col.name = "origin")
skincl3.sct<-AddMetaData(skincl3.sct, "XL", col.name = "origin")
skincl4.sct<-AddMetaData(skincl4.sct, "XL", col.name = "origin")
nscar1.sct<-AddMetaData(nscar1.sct, "MD", col.name = "origin")
nscar2.sct<-AddMetaData(nscar2.sct, "MD", col.name = "origin")
nscar3.sct<-AddMetaData(nscar3.sct, "MD", col.name = "origin")
nscar4.sct<-AddMetaData(nscar4.sct, "CCD", col.name = "origin")
nscar5.sct<-AddMetaData(nscar5.sct, "CCD", col.name = "origin")
nscar6.sct<-AddMetaData(nscar6.sct, "CCD", col.name = "origin")
keloid1.sct<-AddMetaData(keloid1.sct, "MD", col.name = "origin")
keloid2.sct<-AddMetaData(keloid2.sct, "MD", col.name = "origin")
keloid3L.sct<-AddMetaData(keloid3L.sct, "MD", col.name = "origin")
keloid3R.sct<-AddMetaData(keloid3R.sct, "MD", col.name = "origin")
keloid4.sct<-AddMetaData(keloid4.sct, "CCD", col.name = "origin")
keloid5.sct<-AddMetaData(keloid5.sct, "CCD", col.name = "origin")
keloid6.sct<-AddMetaData(keloid6.sct, "CCD", col.name = "origin")
keloid7.sct<-AddMetaData(keloid7.sct, "XL", col.name = "origin")
keloid8.sct<-AddMetaData(keloid8.sct, "XL", col.name = "origin")
keloid9.sct<-AddMetaData(keloid9.sct, "XL", col.name = "origin")
keloid10.sct<-AddMetaData(keloid10.sct, "XL", col.name = "origin")

#flag-4_location
skin1.sct<-AddMetaData(skin1.sct, "skin_sep", col.name = "location")
skin2.sct<-AddMetaData(skin2.sct, "skin_sep", col.name = "location")
skin3.sct<-AddMetaData(skin3.sct, "skin_sep", col.name = "location")
skin4.sct<-AddMetaData(skin4.sct, "skin_sep", col.name = "location")
skin5.sct<-AddMetaData(skin5.sct, "skin_sep", col.name = "location")
skin6.sct<-AddMetaData(skin6.sct, "skin_sep", col.name = "location")
skincl1.sct<-AddMetaData(skincl1.sct, "skin_adj", col.name = "location")
skincl2.sct<-AddMetaData(skincl2.sct, "skin_adj", col.name = "location")
skincl3.sct<-AddMetaData(skincl3.sct, "skin_adj", col.name = "location")
skincl4.sct<-AddMetaData(skincl4.sct, "skin_adj", col.name = "location")
nscar1.sct<-AddMetaData(nscar1.sct, "nscar_centr", col.name = "location")
nscar2.sct<-AddMetaData(nscar2.sct, "nscar_centr", col.name = "location")
nscar3.sct<-AddMetaData(nscar3.sct, "nscar_centr", col.name = "location")
nscar4.sct<-AddMetaData(nscar4.sct, "nscar_centr", col.name = "location")
nscar5.sct<-AddMetaData(nscar5.sct, "nscar_centr", col.name = "location")
nscar6.sct<-AddMetaData(nscar6.sct, "nscar_centr", col.name = "location")
keloid1.sct<-AddMetaData(keloid1.sct, "keloid_centr", col.name = "location")
keloid2.sct<-AddMetaData(keloid2.sct, "keloid_centr", col.name = "location")
keloid3L.sct<-AddMetaData(keloid3L.sct, "keloid_centr", col.name = "location")
keloid3R.sct<-AddMetaData(keloid3R.sct, "keloid_centr", col.name = "location")
keloid4.sct<-AddMetaData(keloid4.sct, "keloid_mix", col.name = "location")
keloid5.sct<-AddMetaData(keloid5.sct, "keloid_mix", col.name = "location")
keloid6.sct<-AddMetaData(keloid6.sct, "keloid_mix", col.name = "location")
keloid7.sct<-AddMetaData(keloid7.sct, "keloid_centr", col.name = "location")
keloid8.sct<-AddMetaData(keloid8.sct, "keloid_centr", col.name = "location")
keloid9.sct<-AddMetaData(keloid9.sct, "keloid_centr", col.name = "location")
keloid10.sct<-AddMetaData(keloid10.sct, "keloid_centr", col.name = "location")

###skin_tabib###Quality ctrl, scale, normalize
#Quality ctrl, scale, normalize, skin1.sct
skin1.sct[["percent.mt"]] <- PercentageFeatureSet(skin1.sct, pattern = "^MT-")
skin1.sct[["percent.ERY"]] <- PercentageFeatureSet(skin1.sct, pattern = "HBB")
skin1.sct <- subset(skin1.sct, subset = percent.ERY < 5 & percent.mt < 5)

VlnPlot(skin1.sct, features = c("nFeature_RNA", "nCount_RNA","percent.mt"), ncol = 3)
FeatureScatter(skin1.sct, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
FeatureScatter(skin1.sct, feature1 = "nCount_RNA", feature2 = "percent.mt")

skin1.sct <- SCTransform(skin1.sct, method = "glmGamPoi",vars.to.regress = "percent.mt",verbose = F)

#Quality ctrl, scale, normalize, skin2.sct
skin2.sct[["percent.mt"]] <- PercentageFeatureSet(skin2.sct, pattern = "^MT-")
skin2.sct[["percent.ERY"]] <- PercentageFeatureSet(skin2.sct, pattern = "HBB")
skin2.sct <- subset(skin2.sct, subset = percent.ERY < 5 & percent.mt < 5)

VlnPlot(skin2.sct, features = c("nFeature_RNA", "nCount_RNA","percent.mt"), ncol = 3)
FeatureScatter(skin2.sct, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
FeatureScatter(skin2.sct, feature1 = "nCount_RNA", feature2 = "percent.mt")

skin2.sct <- SCTransform(skin2.sct, method = "glmGamPoi",vars.to.regress = "percent.mt",verbose = F)

#Quality ctrl, scale, normalize, skin3.sct
skin3.sct[["percent.mt"]] <- PercentageFeatureSet(skin3.sct, pattern = "^MT-")
skin3.sct[["percent.ERY"]] <- PercentageFeatureSet(skin3.sct, pattern = "HBB")
skin3.sct <- subset(skin3.sct, subset = percent.ERY < 5 & percent.mt < 5)

VlnPlot(skin3.sct, features = c("nFeature_RNA", "nCount_RNA","percent.mt"), ncol = 3)
FeatureScatter(skin3.sct, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
FeatureScatter(skin3.sct, feature1 = "nCount_RNA", feature2 = "percent.mt")

skin3.sct <- SCTransform(skin3.sct, method = "glmGamPoi",vars.to.regress = "percent.mt",verbose = F)

#Quality ctrl, scale, normalize, skin4.sct
skin4.sct[["percent.mt"]] <- PercentageFeatureSet(skin4.sct, pattern = "^MT-")
skin4.sct[["percent.ERY"]] <- PercentageFeatureSet(skin4.sct, pattern = "HBB")
skin4.sct <- subset(skin4.sct, subset = percent.ERY < 5 & percent.mt < 5)

VlnPlot(skin4.sct, features = c("nFeature_RNA", "nCount_RNA","percent.mt"), ncol = 3)
FeatureScatter(skin4.sct, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
FeatureScatter(skin4.sct, feature1 = "nCount_RNA", feature2 = "percent.mt")

skin4.sct <- SCTransform(skin4.sct, method = "glmGamPoi",vars.to.regress = "percent.mt",verbose = F)

#Quality ctrl, scale, normalize, skin5.sct
skin5.sct[["percent.mt"]] <- PercentageFeatureSet(skin5.sct, pattern = "^MT-")
skin5.sct[["percent.ERY"]] <- PercentageFeatureSet(skin5.sct, pattern = "HBB")
skin5.sct <- subset(skin5.sct, subset = percent.ERY < 5 & percent.mt < 5)

VlnPlot(skin5.sct, features = c("nFeature_RNA", "nCount_RNA","percent.mt"), ncol = 3)
FeatureScatter(skin5.sct, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
FeatureScatter(skin5.sct, feature1 = "nCount_RNA", feature2 = "percent.mt")

skin5.sct <- SCTransform(skin5.sct, method = "glmGamPoi",vars.to.regress = "percent.mt",verbose = F)

#Quality ctrl, scale, normalize, skin6.sct
skin6.sct[["percent.mt"]] <- PercentageFeatureSet(skin6.sct, pattern = "^MT-")
skin6.sct[["percent.ERY"]] <- PercentageFeatureSet(skin6.sct, pattern = "HBB")
skin6.sct <- subset(skin6.sct, subset = percent.ERY < 5 & percent.mt < 5)

VlnPlot(skin6.sct, features = c("nFeature_RNA", "nCount_RNA","percent.mt"), ncol = 3)
FeatureScatter(skin6.sct, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
FeatureScatter(skin6.sct, feature1 = "nCount_RNA", feature2 = "percent.mt")

skin6.sct <- SCTransform(skin6.sct, method = "glmGamPoi",vars.to.regress = "percent.mt",verbose = F)

###skincl_Direder###Quality ctrl, scale, normalize
#Quality ctrl, scale, normalize, skincl1.sct
skincl1.sct[["percent.mt"]] <- PercentageFeatureSet(skincl1.sct, pattern = "^MT-")
skincl1.sct[["percent.ERY"]] <- PercentageFeatureSet(skincl1.sct, pattern = "HBB")
skincl1.sct <- subset(skincl1.sct, subset = percent.ERY < 5 & percent.mt < 5)

VlnPlot(skincl1.sct, features = c("nFeature_RNA", "nCount_RNA","percent.mt"), ncol = 3)
FeatureScatter(skincl1.sct, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
FeatureScatter(skincl1.sct, feature1 = "nCount_RNA", feature2 = "percent.mt")

skincl1.sct <- SCTransform(skincl1.sct, method = "glmGamPoi",vars.to.regress = "percent.mt",verbose = F)

#Quality ctrl, scale, normalize, skincl2.sct
skincl2.sct[["percent.mt"]] <- PercentageFeatureSet(skincl2.sct, pattern = "^MT-")
skincl2.sct[["percent.ERY"]] <- PercentageFeatureSet(skincl2.sct, pattern = "HBB")
skincl2.sct <- subset(skincl2.sct, subset = percent.ERY < 5 & percent.mt < 5)

VlnPlot(skincl2.sct, features = c("nFeature_RNA", "nCount_RNA","percent.mt"), ncol = 3)
FeatureScatter(skincl2.sct, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
FeatureScatter(skincl2.sct, feature1 = "nCount_RNA", feature2 = "percent.mt")

skincl2.sct <- SCTransform(skincl2.sct, method = "glmGamPoi",vars.to.regress = "percent.mt",verbose = F)

#Quality ctrl, scale, normalize, skincl3.sct
skincl3.sct[["percent.mt"]] <- PercentageFeatureSet(skincl3.sct, pattern = "^MT-")
skincl3.sct[["percent.ERY"]] <- PercentageFeatureSet(skincl3.sct, pattern = "HBB")
skincl3.sct <- subset(skincl3.sct, subset = percent.ERY < 5 & percent.mt < 5)

VlnPlot(skincl3.sct, features = c("nFeature_RNA", "nCount_RNA","percent.mt"), ncol = 3)
FeatureScatter(skincl3.sct, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
FeatureScatter(skincl3.sct, feature1 = "nCount_RNA", feature2 = "percent.mt")

skincl3.sct <- SCTransform(skincl3.sct, method = "glmGamPoi",vars.to.regress = "percent.mt",verbose = F)

#Quality ctrl, scale, normalize, skincl4.sct
skincl4.sct[["percent.mt"]] <- PercentageFeatureSet(skincl4.sct, pattern = "^MT-")
skincl4.sct[["percent.ERY"]] <- PercentageFeatureSet(skincl4.sct, pattern = "HBB")
skincl4.sct <- subset(skincl4.sct, subset = percent.ERY < 5 & percent.mt < 5)

VlnPlot(skincl4.sct, features = c("nFeature_RNA", "nCount_RNA","percent.mt"), ncol = 3)
FeatureScatter(skincl4.sct, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
FeatureScatter(skincl4.sct, feature1 = "nCount_RNA", feature2 = "percent.mt")

skincl4.sct <- SCTransform(skincl4.sct, method = "glmGamPoi",vars.to.regress = "percent.mt",verbose = F)

###nscar_Direder###Quality ctrl, scale, normalize
#Quality ctrl, scale, normalize, nscar1.sct
nscar1.sct[["percent.mt"]] <- PercentageFeatureSet(nscar1.sct, pattern = "^MT-")
nscar1.sct[["percent.ERY"]] <- PercentageFeatureSet(nscar1.sct, pattern = "HBB")
nscar1.sct <- subset(nscar1.sct, subset = percent.ERY < 5 & percent.mt < 5)

VlnPlot(nscar1.sct, features = c("nFeature_RNA", "nCount_RNA","percent.mt"), ncol = 3)
FeatureScatter(nscar1.sct, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
FeatureScatter(nscar1.sct, feature1 = "nCount_RNA", feature2 = "percent.mt")

nscar1.sct <- SCTransform(nscar1.sct, method = "glmGamPoi",vars.to.regress = "percent.mt",verbose = F)

#Quality ctrl, scale, normalize, nscar2.sct
nscar2.sct[["percent.mt"]] <- PercentageFeatureSet(nscar2.sct, pattern = "^MT-")
nscar2.sct[["percent.ERY"]] <- PercentageFeatureSet(nscar2.sct, pattern = "HBB")
nscar2.sct <- subset(nscar2.sct, subset = percent.ERY < 5 & percent.mt < 5)

VlnPlot(nscar2.sct, features = c("nFeature_RNA", "nCount_RNA","percent.mt"), ncol = 3)
FeatureScatter(nscar2.sct, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
FeatureScatter(nscar2.sct, feature1 = "nCount_RNA", feature2 = "percent.mt")

nscar2.sct <- SCTransform(nscar2.sct, method = "glmGamPoi",vars.to.regress = "percent.mt",verbose = F)

#Quality ctrl, scale, normalize, nscar3.sct
nscar3.sct[["percent.mt"]] <- PercentageFeatureSet(nscar3.sct, pattern = "^MT-")
nscar3.sct[["percent.ERY"]] <- PercentageFeatureSet(nscar3.sct, pattern = "HBB")
nscar3.sct <- subset(nscar3.sct, subset = percent.ERY < 5 & percent.mt < 5)

VlnPlot(nscar3.sct, features = c("nFeature_RNA", "nCount_RNA","percent.mt"), ncol = 3)
FeatureScatter(nscar3.sct, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
FeatureScatter(nscar3.sct, feature1 = "nCount_RNA", feature2 = "percent.mt")

nscar3.sct <- SCTransform(nscar3.sct, method = "glmGamPoi",vars.to.regress = "percent.mt",verbose = F)

#Quality ctrl, scale, normalize, nscar4.sct
nscar4.sct[["percent.mt"]] <- PercentageFeatureSet(nscar4.sct, pattern = "^MT-")
nscar4.sct[["percent.ERY"]] <- PercentageFeatureSet(nscar4.sct, pattern = "HBB")
nscar4.sct <- subset(nscar4.sct, subset = percent.ERY < 5 & percent.mt < 5)

VlnPlot(nscar4.sct, features = c("nFeature_RNA", "nCount_RNA","percent.mt"), ncol = 3)
FeatureScatter(nscar4.sct, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
FeatureScatter(nscar4.sct, feature1 = "nCount_RNA", feature2 = "percent.mt")

nscar4.sct <- SCTransform(nscar4.sct, method = "glmGamPoi",vars.to.regress = "percent.mt",verbose = F)

#Quality ctrl, scale, normalize, nscar5.sct
nscar5.sct[["percent.mt"]] <- PercentageFeatureSet(nscar5.sct, pattern = "^MT-")
nscar5.sct[["percent.ERY"]] <- PercentageFeatureSet(nscar5.sct, pattern = "HBB")
nscar5.sct <- subset(nscar5.sct, subset = percent.ERY < 5 & percent.mt < 5)

VlnPlot(nscar5.sct, features = c("nFeature_RNA", "nCount_RNA","percent.mt"), ncol = 3)
FeatureScatter(nscar5.sct, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
FeatureScatter(nscar5.sct, feature1 = "nCount_RNA", feature2 = "percent.mt")

nscar5.sct <- SCTransform(nscar5.sct, method = "glmGamPoi",vars.to.regress = "percent.mt",verbose = F)

#Quality ctrl, scale, normalize, nscar6.sct
nscar6.sct[["percent.mt"]] <- PercentageFeatureSet(nscar6.sct, pattern = "^MT-")
nscar6.sct[["percent.ERY"]] <- PercentageFeatureSet(nscar6.sct, pattern = "HBB")
nscar6.sct <- subset(nscar6.sct, subset = percent.ERY < 5 & percent.mt < 5)

VlnPlot(nscar6.sct, features = c("nFeature_RNA", "nCount_RNA","percent.mt"), ncol = 3)
FeatureScatter(nscar6.sct, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
FeatureScatter(nscar6.sct, feature1 = "nCount_RNA", feature2 = "percent.mt")

nscar6.sct <- SCTransform(nscar6.sct, method = "glmGamPoi",vars.to.regress = "percent.mt",verbose = F)

###keloid###Quality ctrl, scale, normalize
#Quality ctrl, scale, normalize, keloid1.sct
keloid1.sct[["percent.mt"]] <- PercentageFeatureSet(keloid1.sct, pattern = "^MT-")
keloid1.sct[["percent.ERY"]] <- PercentageFeatureSet(keloid1.sct, pattern = "HBB")
keloid1.sct <- subset(keloid1.sct, subset = percent.ERY < 5 & percent.mt < 5)

VlnPlot(keloid1.sct, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
FeatureScatter(keloid1.sct, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
FeatureScatter(keloid1.sct, feature1 = "nCount_RNA", feature2 = "percent.mt")

keloid1.sct <- SCTransform(keloid1.sct,method = "glmGamPoi", vars.to.regress = "percent.mt" ,verbose = F)


#Quality ctrl, scale, normalize, keloid2.sct
keloid2.sct[["percent.mt"]] <- PercentageFeatureSet(keloid2.sct, pattern = "^MT-")
keloid2.sct[["percent.ERY"]] <- PercentageFeatureSet(keloid2.sct, pattern = "HBB")
keloid2.sct <- subset(keloid2.sct, subset = percent.ERY < 5 & percent.mt < 5)

VlnPlot(keloid2.sct, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
FeatureScatter(keloid2.sct, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
FeatureScatter(keloid2.sct, feature1 = "nCount_RNA", feature2 = "percent.mt")

keloid2.sct <- SCTransform(keloid2.sct,method = "glmGamPoi", vars.to.regress = "percent.mt",verbose = F)


#Quality ctrl, scale, normalize, keloid3L.sct
keloid3L.sct[["percent.mt"]] <- PercentageFeatureSet(keloid3L.sct, pattern = "^MT-")
keloid3L.sct[["percent.ERY"]] <- PercentageFeatureSet(keloid3L.sct, pattern = "HBB")
keloid3L.sct <- subset(keloid3L.sct, subset = percent.ERY < 5 & percent.mt < 5)

VlnPlot(keloid3L.sct, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
FeatureScatter(keloid3L.sct, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
FeatureScatter(keloid3L.sct, feature1 = "nCount_RNA", feature2 = "percent.mt")

keloid3L.sct <- SCTransform(keloid3L.sct, method = "glmGamPoi",vars.to.regress = "percent.mt",verbose = F)


#Quality ctrl, scale, normalize, keloid3R.sct
keloid3R.sct[["percent.mt"]] <- PercentageFeatureSet(keloid3R.sct, pattern = "^MT-")
keloid3R.sct[["percent.ERY"]] <- PercentageFeatureSet(keloid3R.sct, pattern = "HBB")
keloid3R.sct <- subset(keloid3R.sct, subset = percent.ERY < 5 & percent.mt < 5)

VlnPlot(keloid3R.sct, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
FeatureScatter(keloid3R.sct, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
FeatureScatter(keloid3R.sct, feature1 = "nCount_RNA", feature2 = "percent.mt")

keloid3R.sct <- SCTransform(keloid3R.sct, method = "glmGamPoi",vars.to.regress = "percent.mt",verbose = F)


#Quality ctrl, scale, normalize, keloid4.sct
keloid4.sct[["percent.mt"]] <- PercentageFeatureSet(keloid4.sct, pattern = "^MT-")
keloid4.sct[["percent.ERY"]] <- PercentageFeatureSet(keloid4.sct, pattern = "HBB")
keloid4.sct <- subset(keloid4.sct, subset = percent.ERY < 5 & percent.mt < 5)

VlnPlot(keloid4.sct, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
FeatureScatter(keloid4.sct, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
FeatureScatter(keloid4.sct, feature1 = "nCount_RNA", feature2 = "percent.mt")

keloid4.sct <- SCTransform(keloid4.sct, method = "glmGamPoi",vars.to.regress = "percent.mt",verbose = F)


#Quality ctrl, scale, normalize, keloid5.sct
keloid5.sct[["percent.mt"]] <- PercentageFeatureSet(keloid5.sct, pattern = "^MT-")
keloid5.sct[["percent.ERY"]] <- PercentageFeatureSet(keloid5.sct, pattern = "HBB")
keloid5.sct <- subset(keloid5.sct, subset = percent.ERY < 5 & percent.mt < 5)

VlnPlot(keloid5.sct, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
FeatureScatter(keloid5.sct, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
FeatureScatter(keloid5.sct, feature1 = "nCount_RNA", feature2 = "percent.mt")

keloid5.sct <- SCTransform(keloid5.sct, method = "glmGamPoi",vars.to.regress = "percent.mt",verbose = F)


#Quality ctrl, scale, normalize, keloid6.sct
keloid6.sct[["percent.mt"]] <- PercentageFeatureSet(keloid6.sct, pattern = "^MT-")
keloid6.sct[["percent.ERY"]] <- PercentageFeatureSet(keloid6.sct, pattern = "HBB")
keloid6.sct <- subset(keloid6.sct, subset = percent.ERY < 5 & percent.mt < 5)

VlnPlot(keloid6.sct, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
FeatureScatter(keloid6.sct, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
FeatureScatter(keloid6.sct, feature1 = "nCount_RNA", feature2 = "percent.mt")

keloid6.sct <- SCTransform(keloid6.sct, method = "glmGamPoi",vars.to.regress = "percent.mt",verbose = F)


#Quality ctrl, scale, normalize, keloid7.sct
keloid7.sct[["percent.mt"]] <- PercentageFeatureSet(keloid7.sct, pattern = "^MT-")
keloid7.sct[["percent.ERY"]] <- PercentageFeatureSet(keloid7.sct, pattern = "HBB")
keloid7.sct <- subset(keloid7.sct, subset = percent.ERY < 5 & percent.mt < 5)

VlnPlot(keloid7.sct, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
FeatureScatter(keloid7.sct, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
FeatureScatter(keloid7.sct, feature1 = "nCount_RNA", feature2 = "percent.mt")

keloid7.sct <- SCTransform(keloid7.sct, method = "glmGamPoi",vars.to.regress = "percent.mt",verbose = F)

#Quality ctrl, scale, normalize, keloid8.sct
keloid8.sct[["percent.mt"]] <- PercentageFeatureSet(keloid8.sct, pattern = "^MT-")
keloid8.sct[["percent.ERY"]] <- PercentageFeatureSet(keloid8.sct, pattern = "HBB")
keloid8.sct <- subset(keloid8.sct, subset = percent.ERY < 5 & percent.mt < 5)

VlnPlot(keloid8.sct, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
FeatureScatter(keloid8.sct, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
FeatureScatter(keloid8.sct, feature1 = "nCount_RNA", feature2 = "percent.mt")

keloid8.sct <- SCTransform(keloid8.sct, method = "glmGamPoi",vars.to.regress = "percent.mt",verbose = F)

#Quality ctrl, scale, normalize, keloid9.sct
keloid9.sct[["percent.mt"]] <- PercentageFeatureSet(keloid9.sct, pattern = "^MT-")
keloid9.sct[["percent.ERY"]] <- PercentageFeatureSet(keloid9.sct, pattern = "HBB")
keloid9.sct <- subset(keloid9.sct, subset = percent.ERY < 5 & percent.mt < 5)

VlnPlot(keloid9.sct, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
FeatureScatter(keloid9.sct, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
FeatureScatter(keloid9.sct, feature1 = "nCount_RNA", feature2 = "percent.mt")

keloid9.sct <- SCTransform(keloid9.sct, method = "glmGamPoi",vars.to.regress = "percent.mt",verbose = F)

#Quality ctrl, scale, normalize, keloid10.sct
keloid10.sct[["percent.mt"]] <- PercentageFeatureSet(keloid10.sct, pattern = "^MT-")
keloid10.sct[["percent.ERY"]] <- PercentageFeatureSet(keloid10.sct, pattern = "HBB")
keloid10.sct <- subset(keloid10.sct, subset = percent.ERY < 5 & percent.mt < 5)

VlnPlot(keloid10.sct, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
FeatureScatter(keloid10.sct, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
FeatureScatter(keloid10.sct, feature1 = "nCount_RNA", feature2 = "percent.mt")

keloid10.sct <- SCTransform(keloid10.sct, method = "glmGamPoi",vars.to.regress = "percent.mt",verbose = F)

#integrate data



##### integrate skin TT####
s.TT.list<-list(skin1.sct, skin2.sct, skin3.sct, skin4.sct, skin5.sct, skin6.sct)
s.TT.features <- SelectIntegrationFeatures(object.list = s.TT.list, nfeatures = 500)
s.TT.list <- PrepSCTIntegration(s.TT.list, anchor.features = s.TT.features)
s.TT.list <- lapply(s.TT.list, RunPCA, verbose = F, features= s.TT.features)
s.TT.anchors<-FindIntegrationAnchors(s.TT.list,normalization.method = "SCT", anchor.features = s.TT.features, reduction = "rpca")
s.TT.x <- IntegrateData(anchorset=s.TT.anchors, normalization.method = "SCT", dims = 1:50)
s.TT.y <- s.TT.x
s.TT.y <- RunPCA(s.TT.y, npcs = 60)
ElbowPlot(s.TT.y, ndims = 30)
s.TT.y <- RunUMAP(s.TT.y, dims = 1:30)
s.TT.y <- FindNeighbors(s.TT.y, dims = 1:30)
s.TT.y <- FindClusters(s.TT.y, resolution = 0.6)
UMAPPlot(s.TT.y, label=T)
s.TT<-s.TT.y
DefaultAssay(s.TT)<- "RNA"
s.TT <- NormalizeData(s.TT)

#Cluster sum
  DotPlot(s.TT, features=c("PDGFRA","LUM","COL1A1", "DCN","FBLN1","ACTA2","RGS5","KRT10","KRT1","KRT14","KRT5","SELE","PECAM1","VWF","LYVE1","CD3D","CD2","CXCR4","CD68","AIF1","FCER1A","S100B","NGFR","PMEL","MLANA","MUCL1","SCGB1B2P","SCGB1D2"), assay = "RNA") + coord_flip()+scale_color_gradient(low="grey", "high"=single_col)+theme(axis.text.x = element_text(angle = 90)) + ggtitle("Basic Cluster identification")
  
  #Marker Celltype specific
  DefaultAssay(s.TT)<- "RNA"
  #Name Cluster
  Idents(s.TT)<-s.TT$seurat_clusters
  s.TT<-RenameIdents(s.TT,
                       `0`="FB", `1`="KC", `2`="FB", `3`="EC", `4`="KC", 
                       `5`="FB", `6`="SMC", `7`="KC", `8`="MAC/DC", `9`="TC", `10`="PC",
                       `11`="FB", `12`="EC", `13`="KC", `14`="FB", `15`="SMC", `16`="SGC", `17`="SC", `18`="SMC", `19`="MEL", `20`="FB")
  s.TT$celltype<-Idents(s.TT)
  
  
  #Define cluster levels
  clusters_ordered<-c(
    "FB", 
    "SMC","PC",
    "KC",
    "EC",  
    "TC", "MAC/DC", 
    "SC","MEL","SGC")
  s.TT$celltype<- factor(s.TT$celltype, levels = clusters_ordered)
  Idents(s.TT)<-s.TT$celltype
  
  #clustermarker recheck
  DotPlot(s.TT, features=c("PDGFRA","LUM","COL1A1", "DCN","FBLN1","ACTA2","RGS5","KRT10","KRT1","KRT14","KRT5","SELE","PECAM1","VWF","LYVE1","CD3D","CD2","CXCR4","CD68","AIF1","FCER1A","S100B","NGFR","PMEL","MLANA","MUCL1","SCGB1B2P","SCGB1D2"), assay = "RNA") + coord_flip()+scale_color_gradient(low="grey", "high"=single_col)+theme(axis.text.x = element_text(angle = 90)) + ggtitle("Basic Cluster identification")

  
  #UMAPPlots, named clusters
  DefaultAssay(s.TT)<-"integrated"
  Idents(s.TT)<-"celltype"
  UMAPPlot(s.TT, label=T, label.color = "black",raster=F)+scale_color_manual(values=color_UMAP_sTT)+ ggtitle("Basic-UMAP")

  #Pieplot/cellnumber /tissue/cluster
  s.TT.cluster.frequencies<-as.data.frame(table(s.TT$celltype))
  s.TT.cluster.frequencies$tissue="Skin_TT"
  s.TT.cluster.frequencies<-mutate(s.TT.cluster.frequencies,Percentage= s.TT.cluster.frequencies$Freq/sum(s.TT.cluster.frequencies$Freq)*100)
 
  pie(s.TT.cluster.frequencies$Freq, labels = F , main="Skin_TT",  clockwise = T, col = c(color_UMAP_sTT))

  
#subset SC
Idents(s.TT)<- s.TT$celltype
SC_s.TT <- subset(s.TT,idents = "SC")
DefaultAssay(SC_s.TT)<-"RNA"
SC_s.TT[["percent.mt"]] <- PercentageFeatureSet(SC_s.TT, pattern = "^MT-")
SC_s.TT <- SCTransform(SC_s.TT, method = "glmGamPoi",vars.to.regress = "percent.mt",verbose = F)
SC_s.TT <- RunPCA(SC_s.TT, npcs = 50)
ElbowPlot(SC_s.TT, ndims = 50)
SC_s.TT <- RunUMAP(SC_s.TT, dims = 1:30)
SC_s.TT <- FindNeighbors(SC_s.TT, dims = 1:30)
SC_s.TT <- FindClusters(SC_s.TT, resolution = 0.3)
UMAPPlot(SC_s.TT, label=T)
DefaultAssay(SC_s.TT)<- "RNA"
SC_s.TT <- NormalizeData(SC_s.TT)
FeaturePlot(SC_s.TT,c("S100B","TOP2A","SELE","DCN","SCN7A"))

#Name Cluster
DotPlot(SC_s.TT,features=c("MBP","MPZ","PLP1","PMP22","NCAM1","L1CAM","SCN7A","S100B","NGFR","NES","IGFBP5","CCN3","TOP2A","MKI67","ICAM1","SELE","DCN","LUM"),assay="RNA") +scale_color_gradient(low="grey", "high"=single_col)+ ggtitle("SC-characterisation")

#Idents(SC_s.TT)<-SC_s.TT$seurat_clusters
SC_s.TT<-RenameIdents(SC_s.TT,
                        `0`="SC-Myel+Nonmyel_skin_TT", `1`="SC-FB_skin_TT")
SC_s.TT$celltype<-Idents(SC_s.TT)

#Define cluster levels
clusters_ordered<-c("SC-Myel+Nonmyel_skin_TT", "SC-FB_skin_TT")
SC_s.TT$celltype<- factor(SC_s.TT$celltype, levels = clusters_ordered)
Idents(SC_s.TT)<-SC_s.TT$celltype

#clustermarker recheck
DotPlot(SC_s.TT,features=c("MBP","MPZ","PLP1","PMP22","NCAM1","L1CAM","SCN7A","S100B","NGFR","NES","IGFBP5","CCN3","TOP2A","MKI67","ICAM1","SELE","DCN","LUM"),assay="RNA") +scale_color_gradient(low="grey", "high"=single_col)+ ggtitle("SC-characterisation")+coord_flip()


#UMAPPlots, named clusters
DefaultAssay(SC_s.TT)<-"integrated"
Idents(SC_s.TT)<-"celltype"

UMAPPlot(SC_s.TT, label=T, label.color = "black",raster=F)+scale_color_manual(values=color_SC_sTT)+ ggtitle("SC-UMAP")

#Clustermarker
SC_s.TT.clustermarker<-FindAllMarkers(SC_s.TT,  assay = "RNA", min.pct = 0.25, logfc.threshold = 0.25)
SC_s.TT.clustermarker$Foldchange_UP <- 2^(SC_s.TT.clustermarker$avg_log2FC)
SC_s.TT.clustermarker$Foldchange_DOWN <- 2^(-SC_s.TT.clustermarker$avg_log2FC)
SC_s.TT.clustermarker$Ratio_pct1_pct2 <- (SC_s.TT.clustermarker$pct.1)/(SC_s.TT.clustermarker$pct.2)
write.xlsx(SC_s.TT.clustermarker, "../SC_s.TT.clustermarker.xlsx")

#Barplot Freq up and downregulated Genes (FC>2)
SC_s.TT.clustermarker_subset <- subset(SC_s.TT.clustermarker, Foldchange_UP >=2 | Foldchange_DOWN >=2)
SC_s.TT.clustermarker_subset$Direction <- ifelse (SC_s.TT.clustermarker_subset$Foldchange_UP>=2,"UP","DOWN")
SC_s.TT.clustermarker_subset<- select(SC_s.TT.clustermarker_subset, Direction, cluster)
row.names(SC_s.TT.clustermarker_subset)<-NULL
SC_s.TT.clustermarker_subset$cluster<- as.character(SC_s.TT.clustermarker_subset$cluster)
SC_s.TT.clustermarker_subset<- dplyr::count(SC_s.TT.clustermarker_subset, cluster,Direction) %>% ungroup()
SC_s.TT.clustermarker_subset$Direction <- factor(x = SC_s.TT.clustermarker_subset$Direction, levels = c("UP","DOWN"))
SC_s.TT.clustermarker_subset$cluster <- factor(x = SC_s.TT.clustermarker_subset$cluster, levels = c("SC-Myel+Nonmyel_skin_TT", "SC-FB_skin_TT"))

ggplot(data=SC_s.TT.clustermarker_subset, aes(x=cluster, y=n, fill=Direction)) + theme_classic() + geom_bar(stat="identity", position=position_dodge())+ ggtitle("DEG-SC") +scale_fill_manual(values =color_DEG) + xlab("SC-type")+ylab("total number of DEG")+ scale_x_discrete(guide = guide_axis(angle = 45))


#subset samples
SC_1.subset<-subset(SC_s.TT, sample=="skin_sep1")
SC_2.subset<-subset(SC_s.TT, sample=="skin_sep2")
SC_3.subset<-subset(SC_s.TT, sample=="skin_sep3")
SC_4.subset<-subset(SC_s.TT, sample=="skin_sep4")
SC_5.subset<-subset(SC_s.TT, sample=="skin_sep5")
SC_6.subset<-subset(SC_s.TT, sample=="skin_sep6")

#Data frame frequencies
SC_1.cluster.frequencies<-as.data.frame(table(SC_1.subset$celltype))
SC_2.cluster.frequencies<-as.data.frame(table(SC_2.subset$celltype))
SC_3.cluster.frequencies<-as.data.frame(table(SC_3.subset$celltype))
SC_4.cluster.frequencies<-as.data.frame(table(SC_4.subset$celltype))
SC_5.cluster.frequencies<-as.data.frame(table(SC_5.subset$celltype))
SC_6.cluster.frequencies<-as.data.frame(table(SC_6.subset$celltype))

#Barplot Cluster percentage
SC_1.cluster.frequencies$sample="skin_sep1"
SC_2.cluster.frequencies$sample="skin_sep2"
SC_3.cluster.frequencies$sample="skin_sep3"
SC_4.cluster.frequencies$sample="skin_sep4"
SC_5.cluster.frequencies$sample="skin_sep5"
SC_6.cluster.frequencies$sample="skin_sep6"

SC_1.cluster.frequencies<-mutate(SC_1.cluster.frequencies,Percentage= SC_1.cluster.frequencies$Freq/sum(SC_1.cluster.frequencies$Freq)*100)
SC_2.cluster.frequencies<-mutate(SC_2.cluster.frequencies,Percentage= SC_2.cluster.frequencies$Freq/sum(SC_2.cluster.frequencies$Freq)*100)
SC_3.cluster.frequencies<-mutate(SC_3.cluster.frequencies,Percentage= SC_3.cluster.frequencies$Freq/sum(SC_3.cluster.frequencies$Freq)*100)
SC_4.cluster.frequencies<-mutate(SC_4.cluster.frequencies,Percentage= SC_4.cluster.frequencies$Freq/sum(SC_4.cluster.frequencies$Freq)*100)
SC_5.cluster.frequencies<-mutate(SC_5.cluster.frequencies,Percentage= SC_5.cluster.frequencies$Freq/sum(SC_5.cluster.frequencies$Freq)*100)
SC_6.cluster.frequencies<-mutate(SC_6.cluster.frequencies,Percentage= SC_6.cluster.frequencies$Freq/sum(SC_6.cluster.frequencies$Freq)*100)

SC_samp.cluster.df <- rbind(SC_1.cluster.frequencies,SC_2.cluster.frequencies,SC_3.cluster.frequencies,SC_4.cluster.frequencies,SC_5.cluster.frequencies,SC_6.cluster.frequencies)
SC_samp.cluster.df$sample <- factor(x = SC_samp.cluster.df$sample, levels = c("skin_sep1","skin_sep2","skin_sep3","skin_sep4","skin_sep5","skin_sep6"))


##### integrate skin XL####
s.XL.list<-list(skincl1.sct, skincl2.sct, skincl3.sct, skincl4.sct)
s.XL.features <- SelectIntegrationFeatures(object.list = s.XL.list, nfeatures = 500)
s.XL.list <- PrepSCTIntegration(s.XL.list, anchor.features = s.XL.features)
s.XL.list <- lapply(s.XL.list, RunPCA, verbose = F, features= s.XL.features)
s.XL.anchors<-FindIntegrationAnchors(s.XL.list,normalization.method = "SCT", anchor.features = s.XL.features, reduction = "rpca")
s.XL.x <- IntegrateData(anchorset=s.XL.anchors, normalization.method = "SCT", dims = 1:50)
s.XL.y <- s.XL.x
s.XL.y <- RunPCA(s.XL.y, npcs = 60)
ElbowPlot(s.XL.y, ndims = 30)
s.XL.y <- RunUMAP(s.XL.y, dims = 1:30)
s.XL.y <- FindNeighbors(s.XL.y, dims = 1:30)
s.XL.y <- FindClusters(s.XL.y, resolution = 0.6)
UMAPPlot(s.XL.y, label=T)
s.XL<-s.XL.y
DefaultAssay(s.XL)<- "RNA"
s.XL <- NormalizeData(s.XL)

#Cluster sum

  DotPlot(s.XL, features=c("PDGFRA","LUM","COL1A1", "DCN","FBLN1","ACTA2","RGS5","KRT10","KRT1","KRT14","KRT5","SELE","PECAM1","VWF","LYVE1","CD3D","CD2","CXCR4","CD68","AIF1","FCER1A","S100B","NGFR","PMEL","MLANA","MUCL1","SCGB1B2P","SCGB1D2"), assay = "RNA") + coord_flip()+scale_color_gradient(low="grey", "high"=single_col)+theme(axis.text.x = element_text(angle = 90)) + ggtitle("Basic Cluster identification")
  
  #Marker Celltype specific
  DefaultAssay(s.XL)<- "RNA"
  #Name Cluster
  Idents(s.XL)<-s.XL$seurat_clusters
  s.XL<-RenameIdents(s.XL,
                     `0`="KC", `1`="KC", `2`="KC", `3`="TC", `4`="EC", 
                     `5`="FB", `6`="EC", `7`="PC", `8`="FB", `9`="KC", `10`="SC",
                     `11`="EC", `12`="KC", `13`="SMC", `14`="MAC/DC", `15`="LEC", `16`="SGC", `17`="SGC", `18`="MAC/DC", `19`="SMC")
  s.XL$celltype<-Idents(s.XL)
  
  
  #Define cluster levels
  clusters_ordered<-c(
    "FB", 
    "SMC","PC",
    "KC",
    "EC", "LEC", 
    "TC", "MAC/DC", 
    "SC","SGC")
  s.XL$celltype<- factor(s.XL$celltype, levels = clusters_ordered)
  Idents(s.XL)<-s.XL$celltype
  
 # clustermarker recheck 
  DotPlot(s.XL, features=c("PDGFRA","LUM","COL1A1", "DCN","FBLN1","ACTA2","RGS5","KRT10","KRT1","KRT14","KRT5","SELE","PECAM1","VWF","LYVE1","CD3D","CD2","CXCR4","CD68","AIF1","FCER1A","S100B","NGFR","PMEL","MLANA","MUCL1","SCGB1B2P","SCGB1D2"), assay = "RNA") + coord_flip()+scale_color_gradient(low="grey", "high"=single_col)+theme(axis.text.x = element_text(angle = 90)) + ggtitle("Basic Cluster identification")

  
  #UMAPPlots, named clusters
  DefaultAssay(s.XL)<-"integrated"
  Idents(s.XL)<-"celltype"
  
  UMAPPlot(s.XL, label=T, label.color = "black",raster=F)+scale_color_manual(values=color_UMAP_sXL)+ ggtitle("Basic-UMAP")
  
  #Pieplot/Barplot cellnumber /tissue/cluster
  s.XL.cluster.frequencies<-as.data.frame(table(s.XL$celltype))
  s.XL.cluster.frequencies$tissue="Skin_XL"
  s.XL.cluster.frequencies<-mutate(s.XL.cluster.frequencies,Percentage= s.XL.cluster.frequencies$Freq/sum(s.XL.cluster.frequencies$Freq)*100)
  
  pie(s.XL.cluster.frequencies$Freq, labels = F, main="Skin_XL",  clockwise = T, col = c(color_UMAP_sXL))
  
#subset SC
Idents(s.XL)<- s.XL$celltype
SC_s.XL <- subset(s.XL,idents = "SC")
DefaultAssay(SC_s.XL)<-"RNA"
SC_s.XL[["percent.mt"]] <- PercentageFeatureSet(SC_s.XL, pattern = "^MT-")
SC_s.XL <- SCTransform(SC_s.XL, method = "glmGamPoi",vars.to.regress = "percent.mt",verbose = F)
SC_s.XL <- RunPCA(SC_s.XL, npcs = 50)
ElbowPlot(SC_s.XL, ndims = 50)
SC_s.XL <- RunUMAP(SC_s.XL, dims = 1:30)
SC_s.XL <- FindNeighbors(SC_s.XL, dims = 1:30)
#to exclude MEL more precise and prevent from loosing SCs 
# a resolution =0.4 was chosen, 
# cluster identification was performed afterwards as if it was computed with a resolution of 0.3
SC_s.XL <- FindClusters(SC_s.XL, resolution = 0.3)
UMAPPlot(SC_s.XL, label=T)
SC_s.XL <- FindClusters(SC_s.XL, resolution = 0.4)
UMAPPlot(SC_s.XL, label=T)
DefaultAssay(SC_s.XL)<- "RNA"
SC_s.XL <- NormalizeData(SC_s.XL)
DotPlot(SC_s.XL,features=c("MBP","MPZ","PLP1","PMP22","NCAM1","L1CAM","SCN7A","S100B","NGFR","NES","IGFBP5","CCN3","TOP2A","MKI67","ICAM1","SELE","DCN","LUM","MEL","PMEL"),assay="RNA") +scale_color_gradient(low="grey", "high"=single_col)+ ggtitle("SC-characterisation")
SC_s.XL <- subset(SC_s.XL,idents = c("4","3","2","0"))
DotPlot(SC_s.XL,features=c("MBP","MPZ","PLP1","PMP22","NCAM1","L1CAM","SCN7A","S100B","NGFR","NES","IGFBP5","CCN3","TOP2A","MKI67","ICAM1","SELE","DCN","LUM","MEL","PMEL"),assay="RNA") +scale_color_gradient(low="grey", "high"=single_col)+ ggtitle("SC-characterisation")

Idents(SC_s.XL)<-SC_s.XL$seurat_clusters
SC_s.XL<-RenameIdents(SC_s.XL,
                      `0`="SC-Nonmyel_skin_XL", `2`="SC-Nonmyel_skin_XL", `3`="SC-Keloid_skin_XL", `4`="SC-Myel_skin_XL")
SC_s.XL$celltype<-Idents(SC_s.XL)

#Define cluster levels
clusters_ordered<-c("SC-Myel_skin_XL","SC-Nonmyel_skin_XL","SC-Keloid_skin_XL")
SC_s.XL$celltype<- factor(SC_s.XL$celltype, levels = clusters_ordered)
Idents(SC_s.XL)<-SC_s.XL$celltype

DotPlot(SC_s.XL,features=c("MBP","MPZ","PLP1","PMP22","NCAM1","L1CAM","SCN7A","S100B","NGFR","NES","IGFBP5","CCN3","TOP2A","MKI67","ICAM1","SELE","DCN","LUM","PMEL","MLANA"),assay="RNA") +scale_color_gradient(low="grey", "high"=single_col)+ ggtitle("SC-characterisation")+coord_flip()


#UMAPPlots, named clusters
DefaultAssay(SC_s.XL)<-"integrated"
Idents(SC_s.XL)<-"celltype"

UMAPPlot(SC_s.XL, label=T, label.color = "black",raster=F)+scale_color_manual(values=color_SC_sXL)+ ggtitle("SC-UMAP")

#Clustermarker
SC_s.XL.clustermarker<-FindAllMarkers(SC_s.XL,  assay = "RNA", min.pct = 0.25, logfc.threshold = 0.25)
SC_s.XL.clustermarker$Foldchange_UP <- 2^(SC_s.XL.clustermarker$avg_log2FC)
SC_s.XL.clustermarker$Foldchange_DOWN <- 2^(-SC_s.XL.clustermarker$avg_log2FC)
SC_s.XL.clustermarker$Ratio_pct1_pct2 <- (SC_s.XL.clustermarker$pct.1)/(SC_s.XL.clustermarker$pct.2)
write.xlsx(SC_s.XL.clustermarker, "Lists/SC_s.XL.clustermarker.xlsx")

#Barplot Freq up and downregulated Genes (FC>2)
SC_s.XL.clustermarker_subset <- subset(SC_s.XL.clustermarker, Foldchange_UP >=2 | Foldchange_DOWN >=2)
SC_s.XL.clustermarker_subset$Direction <- ifelse (SC_s.XL.clustermarker_subset$Foldchange_UP>=2,"UP","DOWN")
SC_s.XL.clustermarker_subset<- select(SC_s.XL.clustermarker_subset, Direction, cluster)
row.names(SC_s.XL.clustermarker_subset)<-NULL
SC_s.XL.clustermarker_subset$cluster<- as.character(SC_s.XL.clustermarker_subset$cluster)
SC_s.XL.clustermarker_subset<- dplyr::count(SC_s.XL.clustermarker_subset, cluster,Direction) %>% ungroup()
SC_s.XL.clustermarker_subset$Direction <- factor(x = SC_s.XL.clustermarker_subset$Direction, levels = c("UP","DOWN"))
SC_s.XL.clustermarker_subset$cluster <- factor(x = SC_s.XL.clustermarker_subset$cluster, levels = c("SC-Myel_skin_XL","SC-Nonmyel_skin_XL","SC-Keloid_skin_XL"))

ggplot(data=SC_s.XL.clustermarker_subset, aes(x=cluster, y=n, fill=Direction)) + theme_classic() + geom_bar(stat="identity", position=position_dodge())+ ggtitle("DEG-SC") +scale_fill_manual(values =color_DEG) + xlab("SC-type")+ylab("total number of DEG")+ scale_x_discrete(guide = guide_axis(angle = 45))

#subset samples
SC_1.subset<-subset(SC_s.XL, sample=="skin_adj1")
SC_2.subset<-subset(SC_s.XL, sample=="skin_adj2")
SC_3.subset<-subset(SC_s.XL, sample=="skin_adj3")
SC_4.subset<-subset(SC_s.XL, sample=="skin_adj4")

#Data frame frequencies
SC_1.cluster.frequencies<-as.data.frame(table(SC_1.subset$celltype))
SC_2.cluster.frequencies<-as.data.frame(table(SC_2.subset$celltype))
SC_3.cluster.frequencies<-as.data.frame(table(SC_3.subset$celltype))
SC_4.cluster.frequencies<-as.data.frame(table(SC_4.subset$celltype))

#Barplot Cluster percentage
SC_1.cluster.frequencies$sample="skin_adj1"
SC_2.cluster.frequencies$sample="skin_adj2"
SC_3.cluster.frequencies$sample="skin_adj3"
SC_4.cluster.frequencies$sample="skin_adj4"

SC_1.cluster.frequencies<-mutate(SC_1.cluster.frequencies,Percentage= SC_1.cluster.frequencies$Freq/sum(SC_1.cluster.frequencies$Freq)*100)
SC_2.cluster.frequencies<-mutate(SC_2.cluster.frequencies,Percentage= SC_2.cluster.frequencies$Freq/sum(SC_2.cluster.frequencies$Freq)*100)
SC_3.cluster.frequencies<-mutate(SC_3.cluster.frequencies,Percentage= SC_3.cluster.frequencies$Freq/sum(SC_3.cluster.frequencies$Freq)*100)
SC_4.cluster.frequencies<-mutate(SC_4.cluster.frequencies,Percentage= SC_4.cluster.frequencies$Freq/sum(SC_4.cluster.frequencies$Freq)*100)

SC_samp.cluster.df <- rbind(SC_1.cluster.frequencies,SC_2.cluster.frequencies,SC_3.cluster.frequencies,SC_4.cluster.frequencies)
SC_samp.cluster.df$sample <- factor(x = SC_samp.cluster.df$sample, levels = c("skin_adj1","skin_adj2","skin_adj3","skin_adj4"))


##### integrate normal scar CCD####
ns.CCD.list<-list(nscar4.sct, nscar5.sct, nscar6.sct)
ns.CCD.features <- SelectIntegrationFeatures(object.list = ns.CCD.list, nfeatures = 500)
ns.CCD.list <- PrepSCTIntegration(ns.CCD.list, anchor.features = ns.CCD.features)
ns.CCD.list <- lapply(ns.CCD.list, RunPCA, verbose = F, features= ns.CCD.features)
ns.CCD.anchors<-FindIntegrationAnchors(ns.CCD.list,normalization.method = "SCT", anchor.features = ns.CCD.features, reduction = "rpca")
ns.CCD.x <- IntegrateData(anchorset=ns.CCD.anchors, normalization.method = "SCT", dims = 1:50)
ns.CCD.y <- ns.CCD.x
ns.CCD.y <- RunPCA(ns.CCD.y, npcs = 60)
ElbowPlot(ns.CCD.y, ndims = 30)
ns.CCD.y <- RunUMAP(ns.CCD.y, dims = 1:30)
ns.CCD.y <- FindNeighbors(ns.CCD.y, dims = 1:30)
ns.CCD.y <- FindClusters(ns.CCD.y, resolution = 0.6)
UMAPPlot(ns.CCD.y, label=T)
ns.CCD<-ns.CCD.y
DefaultAssay(ns.CCD)<- "RNA"
ns.CCD <- NormalizeData(ns.CCD)

#Cluster sum

  DotPlot(ns.CCD, features=c("PDGFRA","LUM","COL1A1", "DCN","FBLN1","ACTA2","RGS5","KRT10","KRT1","KRT14","KRT5","SELE","PECAM1","VWF","LYVE1","CD3D","CD2","CXCR4","CD68","AIF1","FCER1A","S100B","NGFR","PMEL","MLANA","MUCL1","SCGB1B2P","SCGB1D2"), assay = "RNA") + coord_flip()+scale_color_gradient(low="grey", "high"=single_col)+theme(axis.text.x = element_text(angle = 90)) + ggtitle("Basic Cluster identification")
  
  #Marker Celltype specific
  DefaultAssay(ns.CCD)<- "RNA"
  #Name Cluster
  Idents(ns.CCD)<-ns.CCD$seurat_clusters
  ns.CCD<-RenameIdents(ns.CCD,
                      `0`="FB", `1`="FB", `2`="EC", `3`="KC", `4`="SMC/PC", 
                      `5`="EC", `6`="EC", `7`="FB", `8`="KC", `9`="FB", `10`="KC",
                      `11`="KC", `12`="KC", `13`="LEC", `14`="MAC/DC", `15`="FB", `16`="TC")
  ns.CCD$celltype<-Idents(ns.CCD)
  
  
  #Define cluster levels
  clusters_ordered<-c(
    "FB", 
    "SMC/PC",
    "KC",
    "EC", "LEC", 
    "TC", "MAC/DC")
  ns.CCD$celltype<- factor(ns.CCD$celltype, levels = clusters_ordered)
  Idents(ns.CCD)<-ns.CCD$celltype
  
  #Clustermarker recheck
  DotPlot(ns.CCD, features=c("PDGFRA","LUM","COL1A1", "DCN","FBLN1","ACTA2","RGS5","KRT10","KRT1","KRT14","KRT5","SELE","PECAM1","VWF","LYVE1","CD3D","CD2","CXCR4","CD68","AIF1","FCER1A","S100B","NGFR","PMEL","MLANA","MUCL1","SCGB1B2P","SCGB1D2"), assay = "RNA") + coord_flip()+scale_color_gradient(low="grey", "high"=single_col)+theme(axis.text.x = element_text(angle = 90)) + ggtitle("Basic Cluster identification")
  
  #UMAPPlots, named clusters
  DefaultAssay(ns.CCD)<-"integrated"
  Idents(ns.CCD)<-"celltype"
  
  UMAPPlot(ns.CCD, label=T, label.color = "black",raster=F)+scale_color_manual(values=color_UMAP_nCCD)+ ggtitle("Basic-UMAP")
  
  #Pieplot/Barplot cellnumber /tissue/cluster
  ns.CCD.cluster.frequencies<-as.data.frame(table(ns.CCD$celltype))
  ns.CCD.cluster.frequencies$tissue="Normalscar_CCD"
  ns.CCD.cluster.frequencies<-mutate(ns.CCD.cluster.frequencies,Percentage= ns.CCD.cluster.frequencies$Freq/sum(ns.CCD.cluster.frequencies$Freq)*100)
  
  pie(ns.CCD.cluster.frequencies$Freq, labels = F, main="Normalscar_CCD",  clockwise = T, col = c(color_UMAP_nCCD))
  
#subset SC -> no SC identified


##### integrate normal scar MD####
ns.MD.list<-list(nscar1.sct, nscar2.sct, nscar3.sct)
ns.MD.features <- SelectIntegrationFeatures(object.list = ns.MD.list, nfeatures = 500)
ns.MD.list <- PrepSCTIntegration(ns.MD.list, anchor.features = ns.MD.features)
ns.MD.list <- lapply(ns.MD.list, RunPCA, verbose = F, features= ns.MD.features)
ns.MD.anchors<-FindIntegrationAnchors(ns.MD.list,normalization.method = "SCT", anchor.features = ns.MD.features, reduction = "rpca")
ns.MD.x <- IntegrateData(anchorset=ns.MD.anchors, normalization.method = "SCT", dims = 1:50)
ns.MD.y <- ns.MD.x
ns.MD.y <- RunPCA(ns.MD.y, npcs = 60)
ElbowPlot(ns.MD.y, ndims = 30)
ns.MD.y <- RunUMAP(ns.MD.y, dims = 1:30)
ns.MD.y <- FindNeighbors(ns.MD.y, dims = 1:30)
ns.MD.y <- FindClusters(ns.MD.y, resolution = 0.6)
UMAPPlot(ns.MD.y, label=T)
ns.MD<-ns.MD.y
DefaultAssay(ns.MD)<- "RNA"
ns.MD <- NormalizeData(ns.MD)

#Cluster sum

  DotPlot(ns.MD, features=c("PDGFRA","LUM","COL1A1", "DCN","FBLN1","ACTA2","RGS5","KRT10","KRT1","KRT14","KRT5","SELE","PECAM1","VWF","LYVE1","CD3D","CD2","CXCR4","CD68","AIF1","FCER1A","S100B","NGFR","PMEL","MLANA","MUCL1","SCGB1B2P","SCGB1D2"), assay = "RNA") + coord_flip()+scale_color_gradient(low="grey", "high"=single_col)+theme(axis.text.x = element_text(angle = 90)) + ggtitle("Basic Cluster identification")
  
  #Marker Celltype specific
  DefaultAssay(ns.MD)<- "RNA"
  #Name Cluster
  Idents(ns.MD)<-ns.MD$seurat_clusters
  ns.MD<-RenameIdents(ns.MD,
                     `0`="FB", `1`="FB", `2`="EC", `3`="FB", `4`="KC", 
                     `5`="MAC/DC", `6`="SMC/PC", `7`="EC", `8`="KC", `9`="KC", `10`="TC",
                     `11`="LEC", `12`="MEL", `13`="FB")
  ns.MD$celltype<-Idents(ns.MD)
  
  
  #Define cluster levels
  clusters_ordered<-c(
    "FB", 
    "SMC/PC",
    "KC",
    "EC", "LEC", 
    "TC", "MAC/DC", 
    "MEL")
  ns.MD$celltype<- factor(ns.MD$celltype, levels = clusters_ordered)
  Idents(ns.MD)<-ns.MD$celltype
  
  #Clustermarker reckeck
  DotPlot(ns.MD, features=c("PDGFRA","LUM","COL1A1", "DCN","FBLN1","ACTA2","RGS5","KRT10","KRT1","KRT14","KRT5","SELE","PECAM1","VWF","LYVE1","CD3D","CD2","CXCR4","CD68","AIF1","FCER1A","S100B","NGFR","PMEL","MLANA","MUCL1","SCGB1B2P","SCGB1D2"), assay = "RNA") + coord_flip()+scale_color_gradient(low="grey", "high"=single_col)+theme(axis.text.x = element_text(angle = 90)) + ggtitle("Basic Cluster identification")
  
  #UMAPPlots, named clusters
  DefaultAssay(ns.MD)<-"integrated"
  Idents(ns.MD)<-"celltype"
  
  UMAPPlot(ns.MD, label=T, label.color = "black",raster=F)+scale_color_manual(values=color_UMAP_nMD)+ ggtitle("Basic-UMAP")
  
  #Pieplot/Barplot cellnumber /tissue/cluster
  ns.MD.cluster.frequencies<-as.data.frame(table(ns.MD$celltype))
  ns.MD.cluster.frequencies$tissue="Normalscar_MD"
  ns.MD.cluster.frequencies<-mutate(ns.MD.cluster.frequencies,Percentage= ns.MD.cluster.frequencies$Freq/sum(ns.MD.cluster.frequencies$Freq)*100)
  
  pie(ns.MD.cluster.frequencies$Freq, labels = F, main="Normalscar_MD",  clockwise = T, col = c(color_UMAP_nMD))
 
#subset SC -> no SC identified


##### integrate keloid MD####
k.MD.list<-list(keloid1.sct,keloid2.sct, keloid3L.sct,keloid3R.sct)
k.MD.features <- SelectIntegrationFeatures(object.list = k.MD.list, nfeatures = 500)
k.MD.list <- PrepSCTIntegration(k.MD.list, anchor.features = k.MD.features)
k.MD.list <- lapply(k.MD.list, RunPCA, verbose = F, features= k.MD.features)
k.MD.anchors<-FindIntegrationAnchors(k.MD.list,normalization.method = "SCT", anchor.features = k.MD.features, reduction = "rpca")
k.MD.x <- IntegrateData(anchorset=k.MD.anchors, normalization.method = "SCT", dims = 1:50)
k.MD.y <- k.MD.x
k.MD.y <- RunPCA(k.MD.y, npcs = 60)
ElbowPlot(k.MD.y, ndims = 30)
k.MD.y <- RunUMAP(k.MD.y, dims = 1:30)
k.MD.y <- FindNeighbors(k.MD.y, dims = 1:30)
k.MD.y <- FindClusters(k.MD.y, resolution = 0.6)
UMAPPlot(k.MD.y, label=T)
k.MD<-k.MD.y
DefaultAssay(k.MD)<- "RNA"
k.MD <- NormalizeData(k.MD)

#Cluster sum
  DotPlot(k.MD, features=c("PDGFRA","LUM","COL1A1", "DCN","FBLN1","ACTA2","RGS5","KRT10","KRT1","KRT14","KRT5","SELE","PECAM1","VWF","LYVE1","CD3D","CD2","CXCR4","CD68","AIF1","FCER1A","S100B","NGFR","PMEL","MLANA","MUCL1","SCGB1B2P","SCGB1D2"), assay = "RNA") + coord_flip()+scale_color_gradient(low="grey", "high"=single_col)+theme(axis.text.x = element_text(angle = 90)) + ggtitle("Basic Cluster identification")
  
  #Marker Celltype specific
  DefaultAssay(k.MD)<- "RNA"
  #Name Cluster
  Idents(k.MD)<-k.MD$seurat_clusters
  k.MD<-RenameIdents(k.MD,
                          `0`="FB", `1`="EC", `2`="FB", `3`="FB", `4`="SMC/PC", 
                          `5`="EC", `6`="FB", `7`="FB", `8`="KC", `9`="LEC", `10`="FB",
                          `11`="DC", `12`="FB", `13`="TC", `14`="SC", `15`="MAC", `16`="EC", `17`="KC", `18`="FB")
  k.MD$celltype<-Idents(k.MD)
  
  
  #Define cluster levels
  clusters_ordered<-c(
    "FB", 
    "SMC/PC",
    "KC",
    "EC", "LEC", 
    "TC", "MAC","DC", 
    "SC")
  k.MD$celltype<- factor(k.MD$celltype, levels = clusters_ordered)
  Idents(k.MD)<-k.MD$celltype
  
  #Clustermarker recheck
  DotPlot(k.MD, features=c("PDGFRA","LUM","COL1A1", "DCN","FBLN1","ACTA2","RGS5","KRT10","KRT1","KRT14","KRT5","SELE","PECAM1","VWF","LYVE1","CD3D","CD2","CXCR4","CD68","AIF1","FCER1A","S100B","NGFR","PMEL","MLANA","MUCL1","SCGB1B2P","SCGB1D2"), assay = "RNA") + coord_flip()+scale_color_gradient(low="grey", "high"=single_col)+theme(axis.text.x = element_text(angle = 90)) + ggtitle("Basic Cluster identification")
  
  #UMAPPlots, named clusters
  DefaultAssay(k.MD)<-"integrated"
  Idents(k.MD)<-"celltype"
  
  UMAPPlot(k.MD, label=T, label.color = "black",raster=F)+scale_color_manual(values=color_UMAP_kMD)+ ggtitle("Basic-UMAP")
  
  #Pieplot/Barplot cellnumber /tissue/cluster
  k.MD.cluster.frequencies<-as.data.frame(table(k.MD$celltype))
  k.MD.cluster.frequencies$tissue="Keloid_MD"
  k.MD.cluster.frequencies<-mutate(k.MD.cluster.frequencies,Percentage= k.MD.cluster.frequencies$Freq/sum(k.MD.cluster.frequencies$Freq)*100)
  
  pie(k.MD.cluster.frequencies$Freq, labels = F, main="Keloid_MD",  clockwise = T, col = c(color_UMAP_kMD))

#subset SC
Idents(k.MD)<- k.MD$celltype
SC_k.MD <- subset(k.MD,idents = "SC")
DefaultAssay(SC_k.MD)<-"RNA"
SC_k.MD[["percent.mt"]] <- PercentageFeatureSet(SC_k.MD, pattern = "^MT-")
SC_k.MD <- SCTransform(SC_k.MD, method = "glmGamPoi",vars.to.regress = "percent.mt",verbose = F)
SC_k.MD <- RunPCA(SC_k.MD, npcs = 50)
ElbowPlot(SC_k.MD, ndims = 50)
SC_k.MD <- RunUMAP(SC_k.MD, dims = 1:30)
SC_k.MD <- FindNeighbors(SC_k.MD, dims = 1:30)
SC_k.MD <- FindClusters(SC_k.MD, resolution = 0.3)
UMAPPlot(SC_k.MD, label=T)
DefaultAssay(SC_k.MD)<- "RNA"
SC_k.MD <- NormalizeData(SC_k.MD)

#Name Cluster
DotPlot(SC_k.MD,features=c("MBP","MPZ","PLP1","PMP22","NCAM1","L1CAM","SCN7A","S100B","NGFR","NES","IGFBP5","CCN3","TOP2A","MKI67","ICAM1","SELE","DCN","LUM"),assay="RNA") +scale_color_gradient(low="grey", "high"=single_col)+ ggtitle("SC-characterisation")

Idents(SC_k.MD)<-SC_k.MD$seurat_clusters
SC_k.MD<-RenameIdents(SC_k.MD,
                          `0`="SC-Keloid_MD", `1`="SC-FB_MD", `2`="SC-Prolif_MD", `3`="SC-EC_MD",`4`="SC-Myel+Nonmyel_MD")
SC_k.MD$celltype<-Idents(SC_k.MD)

#Define cluster levels
clusters_ordered<-c("SC-Myel+Nonmyel_MD","SC-Keloid_MD","SC-Prolif_MD","SC-EC_MD","SC-FB_MD")
SC_k.MD$celltype<- factor(SC_k.MD$celltype, levels = clusters_ordered)
Idents(SC_k.MD)<-SC_k.MD$celltype

DotPlot(SC_k.MD,features=c("MBP","MPZ","PLP1","PMP22","NCAM1","L1CAM","SCN7A","S100B","NGFR","NES","IGFBP5","CCN3","TOP2A","MKI67","ICAM1","SELE","DCN","LUM"),assay="RNA") +scale_color_gradient(low="grey", "high"=single_col)+ ggtitle("SC-characterisation")+coord_flip()

#UMAPPlots, named clusters
DefaultAssay(SC_k.MD)<-"integrated"
Idents(SC_k.MD)<-"celltype"

UMAPPlot(SC_k.MD, label=T, label.color = "black",raster=F)+scale_color_manual(values=color_SC_kMD)+ ggtitle("SC-UMAP")

#Clustermarker
SC_k.MD.clustermarker<-FindAllMarkers(SC_k.MD,  assay = "RNA", min.pct = 0.25, logfc.threshold = 0.25)
SC_k.MD.clustermarker$Foldchange_UP <- 2^(SC_k.MD.clustermarker$avg_log2FC)
SC_k.MD.clustermarker$Foldchange_DOWN <- 2^(-SC_k.MD.clustermarker$avg_log2FC)
SC_k.MD.clustermarker$Ratio_pct1_pct2 <- (SC_k.MD.clustermarker$pct.1)/(SC_k.MD.clustermarker$pct.2)
write.xlsx(SC_k.MD.clustermarker, "Lists/SC_k.MD.clustermarker.xlsx")

#Barplot Freq up and downregulated Genes (FC>2)
SC_k.MD.clustermarker_subset <- subset(SC_k.MD.clustermarker, Foldchange_UP >=2 | Foldchange_DOWN >=2)
SC_k.MD.clustermarker_subset$Direction <- ifelse (SC_k.MD.clustermarker_subset$Foldchange_UP>=2,"UP","DOWN")
SC_k.MD.clustermarker_subset<- select(SC_k.MD.clustermarker_subset, Direction, cluster)
row.names(SC_k.MD.clustermarker_subset)<-NULL
SC_k.MD.clustermarker_subset$cluster<- as.character(SC_k.MD.clustermarker_subset$cluster)
SC_k.MD.clustermarker_subset<- dplyr::count(SC_k.MD.clustermarker_subset, cluster,Direction) %>% ungroup()
SC_k.MD.clustermarker_subset$Direction <- factor(x = SC_k.MD.clustermarker_subset$Direction, levels = c("UP","DOWN"))
SC_k.MD.clustermarker_subset$cluster <- factor(x = SC_k.MD.clustermarker_subset$cluster, levels = c("SC-Myel+Nonmyel_MD","SC-Keloid_MD","SC-Prolif_MD","SC-EC_MD","SC-FB_MD"))

ggplot(data=SC_k.MD.clustermarker_subset, aes(x=cluster, y=n, fill=Direction)) + theme_classic() + geom_bar(stat="identity", position=position_dodge())+ ggtitle("DEG-SC") +scale_fill_manual(values =color_DEG) + xlab("SC-type")+ylab("total number of DEG")+ scale_x_discrete(guide = guide_axis(angle = 45))


#subset samples
SC_1.subset<-subset(SC_k.MD, sample=="keloid_centr1")
SC_2.subset<-subset(SC_k.MD, sample=="keloid_centr2")
SC_3.subset<-subset(SC_k.MD, sample=="keloid_centr3")
SC_4.subset<-subset(SC_k.MD, sample=="keloid_centr4")

#Data frame frequencies
SC_1.cluster.frequencies<-as.data.frame(table(SC_1.subset$celltype))
SC_2.cluster.frequencies<-as.data.frame(table(SC_2.subset$celltype))
SC_3.cluster.frequencies<-as.data.frame(table(SC_3.subset$celltype))
SC_4.cluster.frequencies<-as.data.frame(table(SC_4.subset$celltype))

#Barplot Cluster percentage
SC_1.cluster.frequencies$sample="keloid_centr1"
SC_2.cluster.frequencies$sample="keloid_centr2"
SC_3.cluster.frequencies$sample="keloid_centr3"
SC_4.cluster.frequencies$sample="keloid_centr4"

SC_1.cluster.frequencies<-mutate(SC_1.cluster.frequencies,Percentage= SC_1.cluster.frequencies$Freq/sum(SC_1.cluster.frequencies$Freq)*100)
SC_2.cluster.frequencies<-mutate(SC_2.cluster.frequencies,Percentage= SC_2.cluster.frequencies$Freq/sum(SC_2.cluster.frequencies$Freq)*100)
SC_3.cluster.frequencies<-mutate(SC_3.cluster.frequencies,Percentage= SC_3.cluster.frequencies$Freq/sum(SC_3.cluster.frequencies$Freq)*100)
SC_4.cluster.frequencies<-mutate(SC_4.cluster.frequencies,Percentage= SC_4.cluster.frequencies$Freq/sum(SC_4.cluster.frequencies$Freq)*100)

SC_samp.cluster.df <- rbind(SC_1.cluster.frequencies,SC_2.cluster.frequencies,SC_3.cluster.frequencies,SC_4.cluster.frequencies)
SC_samp.cluster.df$sample <- factor(x = SC_samp.cluster.df$sample, levels = c("keloid_centr1","keloid_centr2","keloid_centr3","keloid_centr4"))

ggplot(data=SC_samp.cluster.df, aes(x=sample, y=Freq, fill=Var1)) +
  geom_bar(stat="identity", position="stack") +  theme_classic() +
  ggtitle("Cell source_schwann cell cluster") + xlab("condition")+ylab("total cell count") +scale_fill_manual(values =color_SC_kMD)


##### integrate keloid CCD####
k.CCD.list<-list(keloid4.sct,keloid5.sct,keloid6.sct)
k.CCD.features <- SelectIntegrationFeatures(object.list = k.CCD.list, nfeatures = 500)
k.CCD.list <- PrepSCTIntegration(k.CCD.list, anchor.features = k.CCD.features)
k.CCD.list <- lapply(k.CCD.list, RunPCA, verbose = F, features= k.CCD.features)
k.CCD.anchors<-FindIntegrationAnchors(k.CCD.list,normalization.method = "SCT", anchor.features = k.CCD.features, reduction = "rpca")
k.CCD.x <- IntegrateData(anchorset=k.CCD.anchors, normalization.method = "SCT", dims = 1:50)
k.CCD.y <- k.CCD.x
k.CCD.y <- RunPCA(k.CCD.y, npcs = 60)
ElbowPlot(k.CCD.y, ndims = 30)
k.CCD.y <- RunUMAP(k.CCD.y, dims = 1:30)
k.CCD.y <- FindNeighbors(k.CCD.y, dims = 1:30)
k.CCD.y <- FindClusters(k.CCD.y, resolution = 0.6)
UMAPPlot(k.CCD.y, label=T)
k.CCD<-k.CCD.y
DefaultAssay(k.CCD)<- "RNA"
k.CCD <- NormalizeData(k.CCD)

#Cluster sum
  DotPlot(k.CCD, features=c("PDGFRA","LUM","COL1A1", "DCN","FBLN1","ACTA2","RGS5","KRT10","KRT1","KRT14","KRT5","SELE","PECAM1","VWF","LYVE1","CD3D","CD2","CXCR4","CD68","AIF1","FCER1A","S100B","NGFR","PMEL","MLANA","MUCL1","SCGB1B2P","SCGB1D2"), assay = "RNA") + coord_flip()+scale_color_gradient(low="grey", "high"=single_col)+theme(axis.text.x = element_text(angle = 90)) + ggtitle("Basic Cluster identification")
  
  #Marker Celltype specific
  DefaultAssay(k.CCD)<- "RNA"
  #Name Cluster
  Idents(k.CCD)<-k.CCD$seurat_clusters
  k.CCD<-RenameIdents(k.CCD,
                     `0`="EC", `1`="PC", `2`="EC", `3`="EC", `4`="FB", 
                     `5`="FB", `6`="FB", `7`="FB", `8`="KC", `9`="SGC", `10`="LEC",
                     `11`="FB", `12`="MAC/DC", `13`="SMC", `14`="FB", `15`="SC", `16`="FB", `17`="SGC")
  k.CCD$celltype<-Idents(k.CCD)
  
  
  #Define cluster levels
  clusters_ordered<-c(
    "FB", 
    "SMC","PC",
    "KC",
    "EC", "LEC", 
    "MAC/DC", 
    "SC","SGC")
  k.CCD$celltype<- factor(k.CCD$celltype, levels = clusters_ordered)
  Idents(k.CCD)<-k.CCD$celltype
  
  #Clustermarker recheck
  DotPlot(k.CCD, features=c("PDGFRA","LUM","COL1A1", "DCN","FBLN1","ACTA2","RGS5","KRT10","KRT1","KRT14","KRT5","SELE","PECAM1","VWF","LYVE1","CD3D","CD2","CXCR4","CD68","AIF1","FCER1A","S100B","NGFR","PMEL","MLANA","MUCL1","SCGB1B2P","SCGB1D2"), assay = "RNA") + coord_flip()+scale_color_gradient(low="grey", "high"=single_col)+theme(axis.text.x = element_text(angle = 90)) + ggtitle("Basic Cluster identification")
  
  #UMAPPlots, named clusters
  DefaultAssay(k.CCD)<-"integrated"
  Idents(k.CCD)<-"celltype"
  
  UMAPPlot(k.CCD, label=T, label.color = "black",raster=F)+scale_color_manual(values=color_UMAP_kCCD)+ ggtitle("Basic-UMAP")
  
  #Pieplot/Barplot cellnumber /tissue/cluster
  k.CCD.cluster.frequencies<-as.data.frame(table(k.CCD$celltype))
  k.CCD.cluster.frequencies$tissue="Keloid_CCD"
  k.CCD.cluster.frequencies<-mutate(k.CCD.cluster.frequencies,Percentage= k.CCD.cluster.frequencies$Freq/sum(k.CCD.cluster.frequencies$Freq)*100)
 
  pie(k.CCD.cluster.frequencies$Freq, labels = F, main="Keloid_CCD",  clockwise = T, col = c(color_UMAP_kCCD))
  
#subset SC
Idents(k.CCD)<- k.CCD$celltype
SC_k.CCD <- subset(k.CCD,idents = "SC")
DefaultAssay(SC_k.CCD)<-"RNA"
SC_k.CCD[["percent.mt"]] <- PercentageFeatureSet(SC_k.CCD, pattern = "^MT-")
SC_k.CCD <- SCTransform(SC_k.CCD, method = "glmGamPoi",vars.to.regress = "percent.mt",verbose = F)
SC_k.CCD <- RunPCA(SC_k.CCD, npcs = 50)
ElbowPlot(SC_k.CCD, ndims = 50)
SC_k.CCD <- RunUMAP(SC_k.CCD, dims = 1:30)
SC_k.CCD <- FindNeighbors(SC_k.CCD, dims = 1:30)
SC_k.CCD <- FindClusters(SC_k.CCD, resolution = 0.3)
UMAPPlot(SC_k.CCD, label=T)
DefaultAssay(SC_k.CCD)<- "RNA"
SC_k.CCD <- NormalizeData(SC_k.CCD)

#Name Cluster
DotPlot(SC_k.CCD,features=c("MBP","MPZ","PLP1","PMP22","NCAM1","L1CAM","SCN7A","S100B","NGFR","NES","IGFBP5","CCN3","TOP2A","MKI67","ICAM1","SELE","DCN","LUM"),assay="RNA") +scale_color_gradient(low="grey", "high"=single_col)+ ggtitle("SC-characterisation")+coord_flip()

Idents(SC_k.CCD)<-SC_k.CCD$seurat_clusters
SC_k.CCD<-RenameIdents(SC_k.CCD,
                      `0`="SC-Myel+Nonmyel_CCD", `1`="SC-Keloid_CCD")
SC_k.CCD$celltype<-Idents(SC_k.CCD)

#Define cluster levels
clusters_ordered<-c("SC-Myel+Nonmyel_CCD","SC-Keloid_CCD")
SC_k.CCD$celltype<- factor(SC_k.CCD$celltype, levels = clusters_ordered)
Idents(SC_k.CCD)<-SC_k.CCD$celltype

#Clustermarker recheck
DotPlot(SC_k.CCD,features=c("MBP","MPZ","PLP1","PMP22","NCAM1","L1CAM","SCN7A","S100B","NGFR","NES","IGFBP5","CCN3","TOP2A","MKI67","ICAM1","SELE","DCN","LUM"),assay="RNA") +scale_color_gradient(low="grey", "high"=single_col)+ ggtitle("SC-characterisation")+coord_flip()

#UMAPPlots, named clusters
DefaultAssay(SC_k.CCD)<-"integrated"
Idents(SC_k.CCD)<-"celltype"

UMAPPlot(SC_k.CCD, label=T, label.color = "black",raster=F)+scale_color_manual(values=color_SC_kCCD)+ ggtitle("SC-UMAP")

#Clustermarker
SC_k.CCD.clustermarker<-FindAllMarkers(SC_k.CCD,  assay = "RNA", min.pct = 0.25, logfc.threshold = 0.25)
SC_k.CCD.clustermarker$Foldchange_UP <- 2^(SC_k.CCD.clustermarker$avg_log2FC)
SC_k.CCD.clustermarker$Foldchange_DOWN <- 2^(-SC_k.CCD.clustermarker$avg_log2FC)
SC_k.CCD.clustermarker$Ratio_pct1_pct2 <- (SC_k.CCD.clustermarker$pct.1)/(SC_k.CCD.clustermarker$pct.2)
write.xlsx(SC_k.CCD.clustermarker, "Lists/SC_k.CCD.clustermarker.xlsx")

#Barplot Freq up and downregulated Genes (FC>2)
SC_k.CCD.clustermarker_subset <- subset(SC_k.CCD.clustermarker, Foldchange_UP >=2 | Foldchange_DOWN >=2)
SC_k.CCD.clustermarker_subset$Direction <- ifelse (SC_k.CCD.clustermarker_subset$Foldchange_UP>=2,"UP","DOWN")
SC_k.CCD.clustermarker_subset<- select(SC_k.CCD.clustermarker_subset, Direction, cluster)
row.names(SC_k.CCD.clustermarker_subset)<-NULL
SC_k.CCD.clustermarker_subset$cluster<- as.character(SC_k.CCD.clustermarker_subset$cluster)
SC_k.CCD.clustermarker_subset<- dplyr::count(SC_k.CCD.clustermarker_subset, cluster,Direction) %>% ungroup()
SC_k.CCD.clustermarker_subset$Direction <- factor(x = SC_k.CCD.clustermarker_subset$Direction, levels = c("UP","DOWN"))
SC_k.CCD.clustermarker_subset$cluster <- factor(x = SC_k.CCD.clustermarker_subset$cluster, levels = c("SC-Myel+Nonmyel_CCD","SC-Keloid_CCD"))

ggplot(data=SC_k.CCD.clustermarker_subset, aes(x=cluster, y=n, fill=Direction)) + theme_classic() + geom_bar(stat="identity", position=position_dodge())+ ggtitle("DEG-SC") +scale_fill_manual(values =color_DEG) + xlab("SC-type")+ylab("total number of DEG")+  scale_x_discrete(guide = guide_axis(angle = 45))


#subset samples
SC_1.subset<-subset(SC_k.CCD, sample=="keloid_mix1")
SC_2.subset<-subset(SC_k.CCD, sample=="keloid_mix2")
SC_3.subset<-subset(SC_k.CCD, sample=="keloid_mix3")

#Data frame frequencies
SC_1.cluster.frequencies<-as.data.frame(table(SC_1.subset$celltype))
SC_2.cluster.frequencies<-as.data.frame(table(SC_2.subset$celltype))
SC_3.cluster.frequencies<-as.data.frame(table(SC_3.subset$celltype))

#Barplot Cluster percentage
SC_1.cluster.frequencies$sample="keloid_mix1"
SC_2.cluster.frequencies$sample="keloid_mix2"
SC_3.cluster.frequencies$sample="keloid_mix3"

SC_1.cluster.frequencies<-mutate(SC_1.cluster.frequencies,Percentage= SC_1.cluster.frequencies$Freq/sum(SC_1.cluster.frequencies$Freq)*100)
SC_2.cluster.frequencies<-mutate(SC_2.cluster.frequencies,Percentage= SC_2.cluster.frequencies$Freq/sum(SC_2.cluster.frequencies$Freq)*100)
SC_3.cluster.frequencies<-mutate(SC_3.cluster.frequencies,Percentage= SC_3.cluster.frequencies$Freq/sum(SC_3.cluster.frequencies$Freq)*100)

SC_samp.cluster.df <- rbind(SC_1.cluster.frequencies,SC_2.cluster.frequencies,SC_3.cluster.frequencies)
SC_samp.cluster.df$sample <- factor(x = SC_samp.cluster.df$sample, levels = c("keloid_mix1","keloid_mix2","keloid_mix3"))

ggplot(data=SC_samp.cluster.df, aes(x=sample, y=Freq, fill=Var1)) +
  geom_bar(stat="identity", position="stack") +  theme_classic() +
  ggtitle("Cell source_schwann cell cluster") + xlab("condition")+ylab("total cell count") +scale_fill_manual(values =color_SC_kCCD)


##### integrate keloid XL####
k.XL.list<-list(keloid7.sct,keloid8.sct,keloid9.sct,keloid10.sct)
k.XL.features <- SelectIntegrationFeatures(object.list = k.XL.list, nfeatures = 500)
k.XL.list <- PrepSCTIntegration(k.XL.list, anchor.features = k.XL.features)
k.XL.list <- lapply(k.XL.list, RunPCA, verbose = F, features= k.XL.features)
k.XL.anchors<-FindIntegrationAnchors(k.XL.list,normalization.method = "SCT", anchor.features = k.XL.features, reduction = "rpca")
k.XL.x <- IntegrateData(anchorset=k.XL.anchors, normalization.method = "SCT", dims = 1:50)
k.XL.y <- k.XL.x
k.XL.y <- RunPCA(k.XL.y, npcs = 60)
ElbowPlot(k.XL.y, ndims = 30)
k.XL.y <- RunUMAP(k.XL.y, dims = 1:30)
k.XL.y <- FindNeighbors(k.XL.y, dims = 1:30)
k.XL.y <- FindClusters(k.XL.y, resolution = 0.6)
UMAPPlot(k.XL.y, label=T)
k.XL<-k.XL.y
DefaultAssay(k.XL)<- "RNA"
k.XL <- NormalizeData(k.XL)

#Cluster sum
  DotPlot(k.XL, features=c("PDGFRA","LUM","COL1A1", "DCN","FBLN1","ACTA2","RGS5","KRT10","KRT1","KRT14","KRT5","SELE","PECAM1","VWF","LYVE1","CD3D","CD2","CXCR4","CD68","AIF1","FCER1A","S100B","NGFR","PMEL","MLANA","MUCL1","SCGB1B2P","SCGB1D2"), assay = "RNA") + coord_flip()+scale_color_gradient(low="grey", "high"=single_col)+theme(axis.text.x = element_text(angle = 90)) + ggtitle("Basic Cluster identification")
  
  #Marker Celltype specific
  DefaultAssay(k.XL)<- "RNA"
  #Name Cluster
  Idents(k.XL)<-k.XL$seurat_clusters
  k.XL<-RenameIdents(k.XL,
                     `0`="EC", `1`="KC", `2`="EC", `3`="FB", `4`="KC", 
                     `5`="FB", `6`="SMC/PC", `7`="TC", `8`="EC", `9`="EC", `10`="KC",
                     `11`="SC", `12`="EC", `13`="SGC", `14`="KC", `15`="SGC", `16`="LEC", `17`="MAC/DC")
  k.XL$celltype<-Idents(k.XL)
  
  
  #Define cluster levels
  clusters_ordered<-c(
    "FB", 
    "SMC/PC",
    "KC",
    "EC", "LEC", 
    "TC", "MAC/DC", 
    "SC","SGC")
  k.XL$celltype<- factor(k.XL$celltype, levels = clusters_ordered)
  Idents(k.XL)<-k.XL$celltype
  
  #clustermarker recheck
  DotPlot(k.XL, features=c("PDGFRA","LUM","COL1A1", "DCN","FBLN1","ACTA2","RGS5","KRT10","KRT1","KRT14","KRT5","SELE","PECAM1","VWF","LYVE1","CD3D","CD2","CXCR4","CD68","AIF1","FCER1A","S100B","NGFR","PMEL","MLANA","MUCL1","SCGB1B2P","SCGB1D2"), assay = "RNA") + coord_flip()+scale_color_gradient(low="grey", "high"=single_col)+theme(axis.text.x = element_text(angle = 90)) + ggtitle("Basic Cluster identification")
  dev.off()
  
  #UMAPPlots, named clusters
  DefaultAssay(k.XL)<-"integrated"
  Idents(k.XL)<-"celltype"
  
  UMAPPlot(k.XL, label=T, label.color = "black",raster=F)+scale_color_manual(values=color_UMAP_kXL)+ ggtitle("Basic-UMAP")
  
  #Pieplot/Barplot cellnumber /tissue/cluster
  k.XL.cluster.frequencies<-as.data.frame(table(k.XL$celltype))
  k.XL.cluster.frequencies$tissue="Keloid_XL"
  k.XL.cluster.frequencies<-mutate(k.XL.cluster.frequencies,Percentage= k.XL.cluster.frequencies$Freq/sum(k.XL.cluster.frequencies$Freq)*100)
  
  pie(k.XL.cluster.frequencies$Freq, labels = F, main="Keloid_XL",  clockwise = T, col = c(color_UMAP_kXL))
  
#subset SC
Idents(k.XL)<- k.XL$celltype
SC_k.XL <- subset(k.XL,idents = "SC")
DefaultAssay(SC_k.XL)<-"RNA"
SC_k.XL[["percent.mt"]] <- PercentageFeatureSet(SC_k.XL, pattern = "^MT-")
SC_k.XL <- SCTransform(SC_k.XL, method = "glmGamPoi",vars.to.regress = "percent.mt",verbose = F)
SC_k.XL <- RunPCA(SC_k.XL, npcs = 50)
ElbowPlot(SC_k.XL, ndims = 50)
SC_k.XL <- RunUMAP(SC_k.XL, dims = 1:30)
SC_k.XL <- FindNeighbors(SC_k.XL, dims = 1:30)
SC_k.XL <- FindClusters(SC_k.XL, resolution = 0.3)
UMAPPlot(SC_k.XL, label=T)
DefaultAssay(SC_k.XL)<- "RNA"
SC_k.XL <- NormalizeData(SC_k.XL)

#Name Cluster
DotPlot(SC_k.XL,features=c("MBP","MPZ","PLP1","PMP22","NCAM1","L1CAM","SCN7A","S100B","NGFR","NES","IGFBP5","CCN3","TOP2A","MKI67","ICAM1","SELE","DCN","LUM"),assay="RNA") +scale_color_gradient(low="grey", "high"=single_col)+ ggtitle("SC-characterisation")

Idents(SC_k.XL)<-SC_k.XL$seurat_clusters
SC_k.XL<-RenameIdents(SC_k.XL,
                      `0`="SC-Nonmyel1_XL", `1`="SC-Keloid_XL", `2`="SC-Nonmyel2_XL", `3`="SC-Myel_XL")
SC_k.XL$celltype<-Idents(SC_k.XL)

#Define cluster levels
clusters_ordered<-c("SC-Myel_XL","SC-Nonmyel1_XL","SC-Nonmyel2_XL","SC-Keloid_XL")
SC_k.XL$celltype<- factor(SC_k.XL$celltype, levels = clusters_ordered)
Idents(SC_k.XL)<-SC_k.XL$celltype

DotPlot(SC_k.XL,features=c("MBP","MPZ","PLP1","PMP22","NCAM1","L1CAM","SCN7A","S100B","NGFR","NES","IGFBP5","CCN3","TOP2A","MKI67","ICAM1","SELE","DCN","LUM"),assay="RNA") +scale_color_gradient(low="grey", "high"=single_col)+ ggtitle("SC-characterisation")+coord_flip()

#UMAPPlots, named clusters
DefaultAssay(SC_k.XL)<-"integrated"
Idents(SC_k.XL)<-"celltype"

UMAPPlot(SC_k.XL, label=T, label.color = "black",raster=F)+scale_color_manual(values=color_SC_kXL)+ ggtitle("SC-UMAP")

#Clustermarker
SC_k.XL.clustermarker<-FindAllMarkers(SC_k.XL,  assay = "RNA", min.pct = 0.25, logfc.threshold = 0.25)
SC_k.XL.clustermarker$Foldchange_UP <- 2^(SC_k.XL.clustermarker$avg_log2FC)
SC_k.XL.clustermarker$Foldchange_DOWN <- 2^(-SC_k.XL.clustermarker$avg_log2FC)
SC_k.XL.clustermarker$Ratio_pct1_pct2 <- (SC_k.XL.clustermarker$pct.1)/(SC_k.XL.clustermarker$pct.2)
write.xlsx(SC_k.XL.clustermarker, "Lists/SC_k.XL.clustermarker.xlsx")

#Barplot Freq up and downregulated Genes (FC>2)
SC_k.XL.clustermarker_subset <- subset(SC_k.XL.clustermarker, Foldchange_UP >=2 | Foldchange_DOWN >=2)
SC_k.XL.clustermarker_subset$Direction <- ifelse (SC_k.XL.clustermarker_subset$Foldchange_UP>=2,"UP","DOWN")
SC_k.XL.clustermarker_subset<- select(SC_k.XL.clustermarker_subset, Direction, cluster)
row.names(SC_k.XL.clustermarker_subset)<-NULL
SC_k.XL.clustermarker_subset$cluster<- as.character(SC_k.XL.clustermarker_subset$cluster)
SC_k.XL.clustermarker_subset<- dplyr::count(SC_k.XL.clustermarker_subset, cluster,Direction) %>% ungroup()
SC_k.XL.clustermarker_subset$Direction <- factor(x = SC_k.XL.clustermarker_subset$Direction, levels = c("UP","DOWN"))
SC_k.XL.clustermarker_subset$cluster <- factor(x = SC_k.XL.clustermarker_subset$cluster, levels = c("SC-Myel_XL","SC-Nonmyel1_XL","SC-Nonmyel2_XL","SC-Keloid_XL"))

ggplot(data=SC_k.XL.clustermarker_subset, aes(x=cluster, y=n, fill=Direction)) + theme_classic() + geom_bar(stat="identity", position=position_dodge())+ ggtitle("DEG-SC") +scale_fill_manual(values =color_DEG) + xlab("SC-type")+ylab("total number of DEG")+  scale_x_discrete(guide = guide_axis(angle = 45))

#subset samples
SC_1.subset<-subset(SC_k.XL, sample=="keloid_centr5")
SC_2.subset<-subset(SC_k.XL, sample=="keloid_centr6")
SC_3.subset<-subset(SC_k.XL, sample=="keloid_centr7")
SC_4.subset<-subset(SC_k.XL, sample=="keloid_centr8")

#Data frame frequencies
SC_1.cluster.frequencies<-as.data.frame(table(SC_1.subset$celltype))
SC_2.cluster.frequencies<-as.data.frame(table(SC_2.subset$celltype))
SC_3.cluster.frequencies<-as.data.frame(table(SC_3.subset$celltype))
SC_4.cluster.frequencies<-as.data.frame(table(SC_4.subset$celltype))

#Barplot Cluster percentage
SC_1.cluster.frequencies$sample="keloid_centr5"
SC_2.cluster.frequencies$sample="keloid_centr6"
SC_3.cluster.frequencies$sample="keloid_centr7"
SC_4.cluster.frequencies$sample="keloid_centr8"

SC_1.cluster.frequencies<-mutate(SC_1.cluster.frequencies,Percentage= SC_1.cluster.frequencies$Freq/sum(SC_1.cluster.frequencies$Freq)*100)
SC_2.cluster.frequencies<-mutate(SC_2.cluster.frequencies,Percentage= SC_2.cluster.frequencies$Freq/sum(SC_2.cluster.frequencies$Freq)*100)
SC_3.cluster.frequencies<-mutate(SC_3.cluster.frequencies,Percentage= SC_3.cluster.frequencies$Freq/sum(SC_3.cluster.frequencies$Freq)*100)
SC_4.cluster.frequencies<-mutate(SC_4.cluster.frequencies,Percentage= SC_4.cluster.frequencies$Freq/sum(SC_4.cluster.frequencies$Freq)*100)

SC_samp.cluster.df <- rbind(SC_1.cluster.frequencies,SC_2.cluster.frequencies,SC_3.cluster.frequencies,SC_4.cluster.frequencies)
SC_samp.cluster.df$sample <- factor(x = SC_samp.cluster.df$sample, levels = c("keloid_centr5","keloid_centr6","keloid_centr7","keloid_centr8"))

ggplot(data=SC_samp.cluster.df, aes(x=sample, y=Freq, fill=Var1)) +
  geom_bar(stat="identity", position="stack") +  theme_classic() +
  ggtitle("Cell source_schwann cell cluster") + xlab("condition")+ylab("total cell count") +scale_fill_manual(values =color_SC_kXL)


###### All celltypes
All_cells.cluster.df <- rbind(s.TT.cluster.frequencies, s.XL.cluster.frequencies, ns.CCD.cluster.frequencies, ns.MD.cluster.frequencies, k.XL.cluster.frequencies, k.CCD.cluster.frequencies, k.MD.cluster.frequencies)
All_cells.cluster.df$tissue <- factor(x = All_cells.cluster.df$tissue, levels = c("Skin_TT","Skin_XL","Normalscar_CCD","Normalscar_MD", "Keloid_XL", "Keloid_CCD", "Keloid_MD"))
All_cells.cluster.df$Var1 <- factor(x = All_cells.cluster.df$Var1, levels = c("FB","SMC/PC","SMC","PC","KC","EC","LEC","TC","BC","MAC/DC","MAC","DC","SC","MEL","SGC"))

#####Featureblends

FeaturePlot(SC_s.TT,features = c("S100B","LUM"), blend=T,blend.threshold = 0.5, cols = c("lightgrey","#54990F","#990F26"), pt.size = 1, order=T, min.cutoff = 0, max.cutoff = 0.5)
FeaturePlot(SC_s.XL,features = c("S100B","LUM"), blend=T,blend.threshold = 0.5, cols = c("lightgrey","#54990F","#990F26"), pt.size = 1, order=T, min.cutoff = 0, max.cutoff = 0.5)
FeaturePlot(SC_k.XL,features = c("S100B","LUM"), blend=T,blend.threshold = 0.5, cols = c("lightgrey","#54990F","#990F26"), pt.size = 1, order=T, min.cutoff = 0, max.cutoff = 0.5)
FeaturePlot(SC_k.CCD,features = c("S100B","LUM"), blend=T,blend.threshold = 0.5, cols = c("lightgrey","#54990F","#990F26"), pt.size = 1, order=T, min.cutoff = 0, max.cutoff = 3) #max.cutoff changed to 3 as 0.5 results in error
FeaturePlot(SC_k.MD,features = c("S100B","LUM"), blend=T,blend.threshold = 0.5, cols = c("lightgrey","#54990F","#990F26"), pt.size = 1, order=T, min.cutoff = 0, max.cutoff = 0.5)


##### integrate all SCs #####
#integrate all SC#
SC_s.i.TT<-SC_s.TT
Idents(SC_s.i.TT)<- SC_s.i.TT$celltype
DefaultAssay(SC_s.i.TT)<-"RNA"
SC_s.i.TT[["percent.mt"]] <- PercentageFeatureSet(SC_s.i.TT, pattern = "^MT-")
SC_s.i.TT <- SCTransform(SC_s.i.TT, method = "glmGamPoi",vars.to.regress = "percent.mt",verbose = F)
SC_s.i.TT[['integrated']] <- NULL
Idents(SC_s.i.TT)<- SC_s.i.TT$celltype

SC_s.i.XL<-SC_s.XL
Idents(SC_s.i.XL)<- SC_s.i.XL$celltype
DefaultAssay(SC_s.i.XL)<-"RNA"
SC_s.i.XL[["percent.mt"]] <- PercentageFeatureSet(SC_s.i.XL, pattern = "^MT-")
SC_s.i.XL <- SCTransform(SC_s.i.XL, method = "glmGamPoi",vars.to.regress = "percent.mt",verbose = F)
SC_s.i.XL[['integrated']] <- NULL
Idents(SC_s.i.XL)<- SC_s.i.XL$celltype

SC_k.i.MD<-SC_k.MD
Idents(SC_k.i.MD)<- SC_k.i.MD$celltype
DefaultAssay(SC_k.i.MD)<-"RNA"
SC_k.i.MD[["percent.mt"]] <- PercentageFeatureSet(SC_k.i.MD, pattern = "^MT-")
SC_k.i.MD <- SCTransform(SC_k.i.MD, method = "glmGamPoi",vars.to.regress = "percent.mt",verbose = F)
SC_k.i.MD[['integrated']] <- NULL

SC_k.i.CCD<-SC_k.CCD
Idents(SC_k.i.CCD)<- SC_k.i.CCD$celltype
DefaultAssay(SC_k.i.CCD)<-"RNA"
SC_k.i.CCD[["percent.mt"]] <- PercentageFeatureSet(SC_k.i.CCD, pattern = "^MT-")
SC_k.i.CCD <- SCTransform(SC_k.i.CCD, method = "glmGamPoi",vars.to.regress = "percent.mt",verbose = F)
SC_k.i.CCD[['integrated']] <- NULL
Idents(SC_k.i.CCD)<- SC_k.i.CCD$celltype

SC_k.i.XL<-SC_k.XL
Idents(SC_k.i.XL)<- SC_k.i.XL$celltype
DefaultAssay(SC_k.i.XL)<-"RNA"
SC_k.i.XL[["percent.mt"]] <- PercentageFeatureSet(SC_k.i.XL, pattern = "^MT-")
SC_k.i.XL <- SCTransform(SC_k.i.XL, method = "glmGamPoi",vars.to.regress = "percent.mt",verbose = F)
SC_k.i.XL[['integrated']] <- NULL
Idents(SC_k.i.XL)<- SC_k.i.XL$celltype

SC_k.i.all.list<-list(SC_k.i.MD,SC_k.i.CCD,SC_k.i.XL,SC_s.i.TT,SC_s.i.XL)
SC_k.i.all.features <- SelectIntegrationFeatures(object.list = SC_k.i.all.list, nfeatures = 3000)
SC_k.i.all.list <- PrepSCTIntegration(SC_k.i.all.list, anchor.features = SC_k.i.all.features)
SC_k.i.all.list <- lapply(SC_k.i.all.list, RunPCA, verbose = F, features= SC_k.i.all.features)
SC_k.i.all.anchors<-FindIntegrationAnchors(SC_k.i.all.list,normalization.method = "SCT", anchor.features = SC_k.i.all.features, reduction = "rpca")
SC_k.i.all<-IntegrateData(anchorset=SC_k.i.all.anchors, normalization.method = "SCT", k.weight = 50)
SC_k.i.all <- RunPCA(SC_k.i.all)
ElbowPlot(SC_k.i.all, ndims =30)
SC_k.i.all <- RunUMAP(SC_k.i.all, dims = 1:18)
SC_k.i.all <- FindNeighbors(SC_k.i.all, dims = 1:18)
SC_k.i.all <- FindClusters(SC_k.i.all, resolution = 0.5)
UMAPPlot(SC_k.i.all, label=F, group.by="celltype")
UMAPPlot(SC_k.i.all, label=F, group.by="seurat_clusters")
DefaultAssay(SC_k.i.all)<- "RNA"
SC_k.i.all <- NormalizeData(SC_k.i.all)


#UMAPPlots, named clusters
DefaultAssay(SC_k.i.all)<-"integrated"
Idents(SC_k.i.all)<-"celltype"

Idents(SC_k.i.all) <- SC_k.i.all$tissue
SC_k.i.all$tissue_origin <- paste(Idents(SC_k.i.all), SC_k.i.all$origin, sep = "_")
Idents(SC_k.i.all)<-SC_k.i.all$tissue_origin
clusters_tissue_ordered<-c("skin_TT","skin_XL","keloid_XL","keloid_CCD","keloid_MD")
SC_k.i.all$tissue_origin <- factor(SC_k.i.all$tissue_origin, levels = clusters_tissue_ordered)
Idents(SC_k.i.all)<-SC_k.i.all$tissue_origin

DefaultAssay(SC_k.i.all)<-"RNA"
SC_k.i.all$celltype <- factor(x = SC_k.i.all$celltype, levels = c("SC-Myel_skin_XL","SC-Myel_XL", "SC-Myel+Nonmyel_skin_TT", "SC-Myel+Nonmyel_CCD", "SC-Myel+Nonmyel_MD","SC-Nonmyel_skin_XL","SC-Nonmyel1_XL","SC-Nonmyel2_XL","SC-Keloid_skin_XL","SC-Keloid_XL", "SC-Keloid_CCD", "SC-Keloid_MD",  "SC-Prolif_MD", "SC-EC_MD",  "SC-FB_skin_TT","SC-FB_MD"))
Idents(SC_k.i.all)<-SC_k.i.all$celltype

DefaultAssay(SC_k.i.all)<- "RNA"
#Name Cluster
Idents(SC_k.i.all)<-SC_k.i.all$celltype
SC_k.i.all<-RenameIdents(SC_k.i.all,
                         "SC-Myel_skin_XL"="SC-Myel", "SC-Myel_XL"="SC-Myel",  "SC-Myel+Nonmyel_skin_TT"="SC-Myel+Nonmyel", "SC-Myel+Nonmyel_CCD"="SC-Myel+Nonmyel", "SC-Myel+Nonmyel_MD"="SC-Myel+Nonmyel", "SC-Nonmyel_skin_XL"="SC-Nonmyel", 
                         "SC-Nonmyel1_XL"="SC-Nonmyel1", "SC-Nonmyel2_XL"="SC-Nonmyel2", "SC-Keloid_skin_XL"="SC-Keloid", "SC-Keloid_XL"="SC-Keloid", "SC-Keloid_CCD"="SC-Keloid", "SC-Keloid_MD"="SC-Keloid", "SC-Prolif_MD"="SC-Prolif", "SC-EC_MD"="SC-EC", "SC-FB_skin_TT"="SC-FB", "SC-FB_MD"="SC-FB")
SC_k.i.all$celltype_basal<-Idents(SC_k.i.all)

#Define cluster levels
clusters_ordered<-c(
  "SC-Myel","SC-Myel+Nonmyel","SC-Nonmyel","SC-Nonmyel1","SC-Nonmyel2","SC-Keloid","SC-Prolif","SC-EC","SC-FB")
SC_k.i.all$celltype_basal<- factor(SC_k.i.all$celltype_basal, levels = clusters_ordered)
Idents(SC_k.i.all)<-SC_k.i.all$celltype_basal

SC_k.i.all$sample <- factor(x = SC_k.i.all$sample, levels = c("skin_sep1","skin_sep2","skin_sep3", "skin_sep4", "skin_sep5", "skin_sep6", "skin_adj1","skin_adj2","skin_adj3", "skin_adj4", "nscar_1","nscar_2","nscar_3","nscar_4", "nscar_5", "nscar_6", "keloid_mix1", "keloid_mix2", "keloid_mix3","keloid_centr1","keloid_centr2","keloid_centr3","keloid_centr4","keloid_centr5","keloid_centr6","keloid_centr7","keloid_centr8"))


p1<-UMAPPlot(SC_k.i.all, label=F, group.by="tissue_origin",  raster=F)+ ggtitle("SC dataset combination")+ scale_color_manual(values=color_datasets)
p2<-UMAPPlot(SC_k.i.all, label=F, group.by="celltype_basal",  raster=F)+ ggtitle("SC-celltype distribution")+ scale_color_manual(values=color_SC_basal)
wrap_plots(p1,p2)

#Percentage total
###Pieplot/Barplot cellnumber /sample/cluster
#subset conditions
SC_skin1.subset<-subset(SC_k.i.all, sample=="skin_sep1")
SC_skin2.subset<-subset(SC_k.i.all, sample=="skin_sep2")
SC_skin4.subset<-subset(SC_k.i.all, sample=="skin_sep4")
SC_skin3.subset<-subset(SC_k.i.all, sample=="skin_sep3")
SC_skin5.subset<-subset(SC_k.i.all, sample=="skin_sep5")
SC_skin6.subset<-subset(SC_k.i.all, sample=="skin_sep6")
SC_skina1.subset<-subset(SC_k.i.all, sample=="skin_adj1")
SC_skina2.subset<-subset(SC_k.i.all, sample=="skin_adj2")
SC_skina3.subset<-subset(SC_k.i.all, sample=="skin_adj3")
SC_skina4.subset<-subset(SC_k.i.all, sample=="skin_adj4")
SC_nscar1.subset<-subset(SC_k.i.all, sample=="nscar_1")
SC_nscar2.subset<-subset(SC_k.i.all, sample=="nscar_2")
SC_nscar3.subset<-subset(SC_k.i.all, sample=="nscar_3")
SC_nscar4.subset<-subset(SC_k.i.all, sample=="nscar_4")
SC_nscar5.subset<-subset(SC_k.i.all, sample=="nscar_5")
SC_nscar6.subset<-subset(SC_k.i.all, sample=="nscar_6")
SC_keloidm1.subset<-subset(SC_k.i.all, sample=="keloid_mix1")
SC_keloidm2.subset<-subset(SC_k.i.all, sample=="keloid_mix2")
SC_keloidm3.subset<-subset(SC_k.i.all, sample=="keloid_mix3")
SC_keloidc1.subset<-subset(SC_k.i.all, sample=="keloid_centr1")
SC_keloidc2.subset<-subset(SC_k.i.all, sample=="keloid_centr2")
SC_keloidc3.subset<-subset(SC_k.i.all, sample=="keloid_centr3")
SC_keloidc4.subset<-subset(SC_k.i.all, sample=="keloid_centr4")
SC_keloidc5.subset<-subset(SC_k.i.all, sample=="keloid_centr5")
SC_keloidc6.subset<-subset(SC_k.i.all, sample=="keloid_centr6")
SC_keloidc7.subset<-subset(SC_k.i.all, sample=="keloid_centr7")
SC_keloidc8.subset<-subset(SC_k.i.all, sample=="keloid_centr8")


#Data frame frequencies
SC_skins1.cluster.frequencies<-as.data.frame(table(SC_skin1.subset$celltype_basal))
SC_skins2.cluster.frequencies<-as.data.frame(table(SC_skin2.subset$celltype_basal))
SC_skins3.cluster.frequencies<-as.data.frame(table(SC_skin3.subset$celltype_basal))
SC_skins4.cluster.frequencies<-as.data.frame(table(SC_skin4.subset$celltype_basal))
SC_skins5.cluster.frequencies<-as.data.frame(table(SC_skin5.subset$celltype_basal))
SC_skins6.cluster.frequencies<-as.data.frame(table(SC_skin6.subset$celltype_basal))
SC_skina1.cluster.frequencies<-as.data.frame(table(SC_skina1.subset$celltype_basal))
SC_skina2.cluster.frequencies<-as.data.frame(table(SC_skina2.subset$celltype_basal))
SC_skina3.cluster.frequencies<-as.data.frame(table(SC_skina3.subset$celltype_basal))
SC_skina4.cluster.frequencies<-as.data.frame(table(SC_skina4.subset$celltype_basal))
SC_nscar1.cluster.frequencies<-as.data.frame(table(SC_nscar1.subset$celltype_basal))
SC_nscar2.cluster.frequencies<-as.data.frame(table(SC_nscar2.subset$celltype_basal))
SC_nscar3.cluster.frequencies<-as.data.frame(table(SC_nscar3.subset$celltype_basal))
SC_nscar4.cluster.frequencies<-as.data.frame(table(SC_nscar4.subset$celltype_basal))
SC_nscar5.cluster.frequencies<-as.data.frame(table(SC_nscar5.subset$celltype_basal))
SC_nscar6.cluster.frequencies<-as.data.frame(table(SC_nscar6.subset$celltype_basal))
SC_keloidm1.cluster.frequencies<-as.data.frame(table(SC_keloidm1.subset$celltype_basal))
SC_keloidm2.cluster.frequencies<-as.data.frame(table(SC_keloidm2.subset$celltype_basal))
SC_keloidm3.cluster.frequencies<-as.data.frame(table(SC_keloidm3.subset$celltype_basal))
SC_keloidc1.cluster.frequencies<-as.data.frame(table(SC_keloidc1.subset$celltype_basal))
SC_keloidc2.cluster.frequencies<-as.data.frame(table(SC_keloidc2.subset$celltype_basal))
SC_keloidc3.cluster.frequencies<-as.data.frame(table(SC_keloidc3.subset$celltype_basal))
SC_keloidc4.cluster.frequencies<-as.data.frame(table(SC_keloidc4.subset$celltype_basal))
SC_keloidc5.cluster.frequencies<-as.data.frame(table(SC_keloidc5.subset$celltype_basal))
SC_keloidc6.cluster.frequencies<-as.data.frame(table(SC_keloidc6.subset$celltype_basal))
SC_keloidc7.cluster.frequencies<-as.data.frame(table(SC_keloidc7.subset$celltype_basal))
SC_keloidc8.cluster.frequencies<-as.data.frame(table(SC_keloidc8.subset$celltype_basal))

#Barplot Cluster percentage
SC_skins1.cluster.frequencies$sample="skin_sep1"
SC_skins2.cluster.frequencies$sample="skin_sep2"
SC_skins3.cluster.frequencies$sample="skin_sep3"
SC_skins4.cluster.frequencies$sample="skin_sep4"
SC_skins5.cluster.frequencies$sample="skin_sep5"
SC_skins6.cluster.frequencies$sample="skin_sep6"
SC_skina1.cluster.frequencies$sample="skin_adj1"
SC_skina2.cluster.frequencies$sample="skin_adj2"
SC_skina3.cluster.frequencies$sample="skin_adj3"
SC_skina4.cluster.frequencies$sample="skin_adj4"
SC_nscar1.cluster.frequencies$sample="nscar_1"
SC_nscar2.cluster.frequencies$sample="nscar_2"
SC_nscar3.cluster.frequencies$sample="nscar_3"
SC_nscar4.cluster.frequencies$sample="nscar_4"
SC_nscar5.cluster.frequencies$sample="nscar_5"
SC_nscar6.cluster.frequencies$sample="nscar_6"
SC_keloidm1.cluster.frequencies$sample="keloid_mix1"
SC_keloidm2.cluster.frequencies$sample="keloid_mix2"
SC_keloidm3.cluster.frequencies$sample="keloid_mix3"
SC_keloidc1.cluster.frequencies$sample="keloid_centr1"
SC_keloidc2.cluster.frequencies$sample="keloid_centr2"
SC_keloidc3.cluster.frequencies$sample="keloid_centr3"
SC_keloidc4.cluster.frequencies$sample="keloid_centr4"
SC_keloidc5.cluster.frequencies$sample="keloid_centr5"
SC_keloidc6.cluster.frequencies$sample="keloid_centr6"
SC_keloidc7.cluster.frequencies$sample="keloid_centr7"
SC_keloidc8.cluster.frequencies$sample="keloid_centr8"

SC_skins1.cluster.frequencies<-mutate(SC_skins1.cluster.frequencies,Percentage= SC_skins1.cluster.frequencies$Freq/sum(SC_skins1.cluster.frequencies$Freq)*100)
SC_skins2.cluster.frequencies<-mutate(SC_skins2.cluster.frequencies,Percentage= SC_skins2.cluster.frequencies$Freq/sum(SC_skins2.cluster.frequencies$Freq)*100)
SC_skins3.cluster.frequencies<-mutate(SC_skins3.cluster.frequencies,Percentage= SC_skins3.cluster.frequencies$Freq/sum(SC_skins3.cluster.frequencies$Freq)*100)
SC_skins4.cluster.frequencies<-mutate(SC_skins4.cluster.frequencies,Percentage= SC_skins4.cluster.frequencies$Freq/sum(SC_skins4.cluster.frequencies$Freq)*100)
SC_skins5.cluster.frequencies<-mutate(SC_skins5.cluster.frequencies,Percentage= SC_skins5.cluster.frequencies$Freq/sum(SC_skins5.cluster.frequencies$Freq)*100)
SC_skins6.cluster.frequencies<-mutate(SC_skins6.cluster.frequencies,Percentage= SC_skins6.cluster.frequencies$Freq/sum(SC_skins6.cluster.frequencies$Freq)*100)
SC_skina1.cluster.frequencies<-mutate(SC_skina1.cluster.frequencies,Percentage= SC_skina1.cluster.frequencies$Freq/sum(SC_skina1.cluster.frequencies$Freq)*100)
SC_skina2.cluster.frequencies<-mutate(SC_skina2.cluster.frequencies,Percentage= SC_skina2.cluster.frequencies$Freq/sum(SC_skina2.cluster.frequencies$Freq)*100)
SC_skina3.cluster.frequencies<-mutate(SC_skina3.cluster.frequencies,Percentage= SC_skina3.cluster.frequencies$Freq/sum(SC_skina3.cluster.frequencies$Freq)*100)
SC_skina4.cluster.frequencies<-mutate(SC_skina4.cluster.frequencies,Percentage= SC_skina4.cluster.frequencies$Freq/sum(SC_skina4.cluster.frequencies$Freq)*100)
SC_nscar1.cluster.frequencies<-mutate(SC_nscar1.cluster.frequencies,Percentage= SC_nscar1.cluster.frequencies$Freq/sum(SC_nscar1.cluster.frequencies$Freq)*100)
SC_nscar2.cluster.frequencies<-mutate(SC_nscar2.cluster.frequencies,Percentage= SC_nscar2.cluster.frequencies$Freq/sum(SC_nscar2.cluster.frequencies$Freq)*100)
SC_nscar3.cluster.frequencies<-mutate(SC_nscar3.cluster.frequencies,Percentage= SC_nscar3.cluster.frequencies$Freq/sum(SC_nscar3.cluster.frequencies$Freq)*100)
SC_nscar4.cluster.frequencies<-mutate(SC_nscar4.cluster.frequencies,Percentage= SC_nscar4.cluster.frequencies$Freq/sum(SC_nscar4.cluster.frequencies$Freq)*100)
SC_nscar5.cluster.frequencies<-mutate(SC_nscar5.cluster.frequencies,Percentage= SC_nscar5.cluster.frequencies$Freq/sum(SC_nscar5.cluster.frequencies$Freq)*100)
SC_nscar6.cluster.frequencies<-mutate(SC_nscar6.cluster.frequencies,Percentage= SC_nscar6.cluster.frequencies$Freq/sum(SC_nscar6.cluster.frequencies$Freq)*100)
SC_keloidm1.cluster.frequencies<-mutate(SC_keloidm1.cluster.frequencies,Percentage= SC_keloidm1.cluster.frequencies$Freq/sum(SC_keloidm1.cluster.frequencies$Freq)*100)
SC_keloidm2.cluster.frequencies<-mutate(SC_keloidm2.cluster.frequencies,Percentage= SC_keloidm2.cluster.frequencies$Freq/sum(SC_keloidm2.cluster.frequencies$Freq)*100)
SC_keloidm3.cluster.frequencies<-mutate(SC_keloidm3.cluster.frequencies,Percentage= SC_keloidm3.cluster.frequencies$Freq/sum(SC_keloidm3.cluster.frequencies$Freq)*100)
SC_keloidc1.cluster.frequencies<-mutate(SC_keloidc1.cluster.frequencies,Percentage= SC_keloidc1.cluster.frequencies$Freq/sum(SC_keloidc1.cluster.frequencies$Freq)*100)
SC_keloidc2.cluster.frequencies<-mutate(SC_keloidc2.cluster.frequencies,Percentage= SC_keloidc2.cluster.frequencies$Freq/sum(SC_keloidc2.cluster.frequencies$Freq)*100)
SC_keloidc3.cluster.frequencies<-mutate(SC_keloidc3.cluster.frequencies,Percentage= SC_keloidc3.cluster.frequencies$Freq/sum(SC_keloidc3.cluster.frequencies$Freq)*100)
SC_keloidc4.cluster.frequencies<-mutate(SC_keloidc4.cluster.frequencies,Percentage= SC_keloidc4.cluster.frequencies$Freq/sum(SC_keloidc4.cluster.frequencies$Freq)*100)
SC_keloidc5.cluster.frequencies<-mutate(SC_keloidc5.cluster.frequencies,Percentage= SC_keloidc5.cluster.frequencies$Freq/sum(SC_keloidc5.cluster.frequencies$Freq)*100)
SC_keloidc6.cluster.frequencies<-mutate(SC_keloidc6.cluster.frequencies,Percentage= SC_keloidc6.cluster.frequencies$Freq/sum(SC_keloidc6.cluster.frequencies$Freq)*100)
SC_keloidc7.cluster.frequencies<-mutate(SC_keloidc7.cluster.frequencies,Percentage= SC_keloidc7.cluster.frequencies$Freq/sum(SC_keloidc7.cluster.frequencies$Freq)*100)
SC_keloidc8.cluster.frequencies<-mutate(SC_keloidc8.cluster.frequencies,Percentage= SC_keloidc8.cluster.frequencies$Freq/sum(SC_keloidc8.cluster.frequencies$Freq)*100)


SC_skin_kel_samp.cluster.df <- rbind(SC_skins1.cluster.frequencies,SC_skins2.cluster.frequencies,SC_skins3.cluster.frequencies,SC_skins4.cluster.frequencies,SC_skins5.cluster.frequencies,SC_skins6.cluster.frequencies,SC_skina1.cluster.frequencies,SC_skina2.cluster.frequencies,SC_skina3.cluster.frequencies,SC_skina4.cluster.frequencies,SC_keloidm1.cluster.frequencies,SC_keloidm2.cluster.frequencies,SC_keloidm3.cluster.frequencies,SC_keloidc1.cluster.frequencies,SC_keloidc2.cluster.frequencies,SC_keloidc3.cluster.frequencies,SC_keloidc4.cluster.frequencies,SC_keloidc5.cluster.frequencies,SC_keloidc6.cluster.frequencies,SC_keloidc7.cluster.frequencies,SC_keloidc8.cluster.frequencies)
SC_skin_kel_samp.cluster.df$sample <- factor(x = SC_skin_kel_samp.cluster.df$sample, levels = c("skin_sep1","skin_sep2","skin_sep3", "skin_sep4", "skin_sep5", "skin_sep6", "skin_adj1","skin_adj2","skin_adj3", "skin_adj4", "keloid_mix1", "keloid_mix2", "keloid_mix3","keloid_centr1","keloid_centr2","keloid_centr3","keloid_centr4","keloid_centr5","keloid_centr6","keloid_centr7","keloid_centr8"))

ggplot(data=SC_skin_kel_samp.cluster.df, aes(x=sample, y=Freq, fill=Var1)) +
  geom_bar(stat="identity", position="stack") +  theme_classic() +
  ggtitle("Cell source_schwann cell cluster") + xlab("condition")+ylab("total cell count") +scale_fill_manual(values =color_SC_basal)

#common genes in top 100
DefaultAssay(SC_k.i.all)<-"RNA"
SC_k.i.all$celltype <- factor(x = SC_k.i.all$celltype, levels = c("SC-Myel+Nonmyel_skin_TT", "SC-FB_skin_TT", "SC-Myel_skin_XL", "SC-Nonmyel_skin_XL", "SC-Keloid_skin_XL","SC-Myel_XL","SC-Nonmyel1_XL","SC-Nonmyel2_XL", "SC-Keloid_XL",  "SC-Myel+Nonmyel_CCD", "SC-Keloid_CCD", "SC-Myel+Nonmyel_MD",   "SC-Keloid_MD", "SC-Prolif_MD", "SC-EC_MD", "SC-FB_MD"))
Idents(SC_k.i.all)<-SC_k.i.all$celltype
SC_k.i.all<-RenameIdents(SC_k.i.all,
                         `SC-Myel+Nonmyel_skin_TT`="SC_skin_TT", `SC-Myel_skin_XL`="SC_skin_XL", `SC-Nonmyel_skin_XL`="SC_skin_XL", `SC-Myel_XL`="SC_keloid_XL",`SC-Nonmyel1_XL`="SC_keloid_XL",`SC-Nonmyel2_XL`="SC_keloid_XL",`SC-Myel+Nonmyel_CCD`="SC_keloid_CCD", `SC-Myel+Nonmyel_MD`="SC_keloid_MD", `SC-Keloid_skin_XL`="SC_skin_XL", `SC-Keloid_XL`="SC_keloid_XL", `SC-Keloid_CCD`="SC_keloid_CCD",`SC-Keloid_MD`="SC_keloid_MD", `SC-Prolif_MD`="SC_keloid_MD", `SC-FB_MD`="SC_keloid_MD", `SC-FB_skin_TT`="SC_skin_TT", `SC-EC_MD`="SC_keloid_MD")
SC_k.i.all$ct.t<-Idents(SC_k.i.all)

#Define cluster levels
clusters_ordered<-c("SC_skin_TT","SC_skin_XL","SC_keloid_XL","SC_keloid_CCD","SC_keloid_MD")
SC_k.i.all$ct.t<- factor(SC_k.i.all$ct.t, levels = clusters_ordered)
Idents(SC_k.i.all)<-SC_k.i.all$ct.t

DotPlot(SC_k.i.all,features=c("MBP","MPZ","PLP1","PMP22","NCAM1","L1CAM","SCN7A","S100B","NGFR","NES","IGFBP5","CCN3","TOP2A","MKI67","DCN","LUM","ICAM1","SELE"),assay="RNA") +scale_color_gradient(low="grey", "high"=single_col)+ ggtitle("SC-characterisation")+ coord_flip()+ theme(axis.text.x =  element_text(angle = 45, hjust = 1, colour = color_SC_all_donor))

DefaultAssay(SC_k.i.all)<-"RNA"
Idents(SC_k.i.all)<-SC_k.i.all$celltype
##subset SC_Keloid, SC myel , SC nonmyel and categorize rep vs myel/nonmyel
SC_rep_mature <- subset(SC_k.i.all,idents = c("SC-Myel_skin_XL","SC-Myel_XL","SC-Myel+Nonmyel_skin_TT","SC-Myel+Nonmyel_CCD","SC-Myel+Nonmyel_MD","SC-Nonmyel_skin_XL","SC-Nonmyel1_XL","SC-Nonmyel2_XL","SC-Keloid_skin_XL","SC-Keloid_XL","SC-Keloid_CCD","SC-Keloid_MD"))

Idents(SC_rep_mature)<-SC_rep_mature$celltype
SC_rep_mature<-RenameIdents(SC_rep_mature,
                            "SC-Myel_skin_XL"="SC-Myel+Nonmyel_skin_XL","SC-Myel_XL"="SC-Myel+Nonmyel_XL","SC-Myel+Nonmyel_skin_TT"="SC-Myel+Nonmyel_skin_TT","SC-Myel+Nonmyel_CCD"="SC-Myel+Nonmyel_CCD","SC-Myel+Nonmyel_MD"="SC-Myel+Nonmyel_MD","SC-Nonmyel_skin_XL"="SC-Myel+Nonmyel_skin_XL","SC-Nonmyel1_XL"="SC-Myel+Nonmyel_XL","SC-Nonmyel2_XL"="SC-Myel+Nonmyel_XL","SC-Keloid_skin_XL"="SC-Keloid_skin_XL","SC-Keloid_XL"="SC-Keloid_XL","SC-Keloid_CCD"="SC-Keloid_CCD","SC-Keloid_MD"="SC-Keloid_MD")
SC_rep_mature$repvsmat<-Idents(SC_rep_mature)

#Define cluster levels
clusters_ordered<-c("SC-Myel+Nonmyel_skin_TT",
                    "SC-Myel+Nonmyel_skin_XL",   "SC-Keloid_skin_XL",
                    "SC-Myel+Nonmyel_XL",        "SC-Keloid_XL",
                    "SC-Myel+Nonmyel_CCD",       "SC-Keloid_CCD",
                    "SC-Myel+Nonmyel_MD",        "SC-Keloid_MD")
SC_rep_mature$repvsmat<- factor(SC_rep_mature$repvsmat, levels = clusters_ordered)
Idents(SC_rep_mature)<-SC_rep_mature$repvsmat



top21 <- c("CALB2","CCN3","COL1A1","COL7A1","CSRP2","ELN","ENC1","IGFBP3","IGFBP5","ITGB1","LOXL2","LYPD1","MMP15","NES","PPP1R14B","S100A16","SH3BGRL3","SPARCL1","TAGLN","TGFBI","TPM2")

VlnPlot(SC_rep_mature,features = top21, pt.size = 0, cols=color_SC_repvsmat, ncol =3) & stat_summary(fun = "mean", geom = "crossbar")

SC_matrisome_shortlist <- c("CCN3","COL1A1","COL3A1","COL4A1","COL4A2","COL5A1","COL5A2","COL5A3","COL6A1","COL6A2","COL7A1","COL8A1","COL9A3","COL12A1","COL16A1","COL18A1","ELN","IGFBP3","IGFBP5","TNC","TGFBI")
Idents(SC_k.i.all)<-"celltype"
VlnPlot(SC_k.i.all,features = SC_matrisome_shortlist, pt.size = 0, cols=color_SC_all_celltype_order, ncol = 3) & stat_summary(fun = "mean", geom = "crossbar")

##### Pseudotime ####

library(monocle3)

#### Create a Monocle CDS Object
# Project PC dimensions to whole data set
SC_convert <- ProjectDim(SC_k.i.all, reduction = "pca")

# Create an expression matrix
expression_matrix <- SC_convert@assays$RNA@counts

# Get cell metadata
cell_metadata <- SC_convert@meta.data
if (all.equal(colnames(expression_matrix), rownames(cell_metadata))) {
  print(sprintf("Cell identifiers match"))
} else {
  print(sprintf("Cell identifier mismatch - %i cells in expression matrix, %i cells in metadata",
                ncol(expression_matrix), nrow(cell_metadata)))
  print("If the counts are equal, sort differences will throw this error")
}

# get gene annotations
gene_annotation <- data.frame(gene_short_name = rownames(SC_convert@assays$RNA), row.names = rownames(SC_convert@assays$RNA))
if (all.equal(rownames(expression_matrix), rownames(gene_annotation))) {
  print(sprintf("Gene identifiers all match"))
} else {
  print(sprintf("Gene identifier mismatch - %i genes in expression matrix, %i gene in gene annotation",
                nrow(expression_matrix), nrow(gene_annotation)))
  print("If the counts are equal, sort differences will throw this error")
}

# Seurat-derived CDS
SC_pseudo.cds <- new_cell_data_set(expression_matrix,
                                     cell_metadata = cell_metadata,
                                     gene_metadata = gene_annotation)

# Transfer Seurat embeddings
# Note that these may be calculated on the Integrated object, not the counts
#   and thus will involve fewer genes
reducedDim(SC_pseudo.cds, type = "PCA") <- SC_k.i.all@reductions$pca@cell.embeddings 
SC_pseudo.cds@preprocess_aux$prop_var_expl <- SC_k.i.all@reductions$pca@stdev
plot_pc_variance_explained(SC_pseudo.cds)

# Transfer Seurat UMAP embeddings
SC_pseudo.cds@int_colData@listData$reducedDims$UMAP <- SC_k.i.all@reductions$umap@cell.embeddings
#    plot_cells(SC_pseudo.cds)

# Copy cluster info from Seurat
SC_pseudo.cds@clusters$UMAP_so$clusters <- SC_k.i.all@meta.data$gt_tp_cell_type_integrated_.0.9

SC_pseudo.cds <- cluster_cells(SC_pseudo.cds, reduction_method = "UMAP")

# Fix from https://gitmemory.com/cole-trapnell-lab
rownames(SC_pseudo.cds@principal_graph_aux$UMAP$dp_mst) <- NULL
colnames(SC_pseudo.cds@int_colData@listData$reducedDims$UMAP) <- NULL

DimPlot(SC_k.i.all, reduction = "umap")

SC_pseudo.cds = cluster_cells(SC_pseudo.cds, resolution = 0.00000000000000000000000001)

p1 <- plot_cells(SC_pseudo.cds,group_label_size = 3.5, cell_size = 2) + ggtitle("SC_k.i.all_ _ cluster")
p2 <- plot_cells(SC_pseudo.cds, color_cells_by="ct", group_cells_by="partition",group_label_size = 3.5, cell_size = 2) + ggtitle("SC_k.i.all_partition")
wrap_plots(p1, p2)

colData(SC_pseudo.cds)$monocle_cluster <- as.character(clusters(SC_pseudo.cds))

SC_pseudo.cds <- learn_graph(SC_pseudo.cds,learn_graph_control = list(prune_graph=F, ncenter=7), use_partition = F)
plot_cells(SC_pseudo.cds, label_groups_by_cluster = T,color_cells_by="celltype", label_leaves = FALSE, label_branch_points = FALSE, cell_size=2, label_cell_groups = F)+ scale_color_manual(values=color_SC_all_celltype) + ggtitle("Principal graph")

SC_pseudo.cds<- order_cells(SC_pseudo.cds)
plot_cells(SC_pseudo.cds, label_groups_by_cluster = T,color_cells_by="pseudotime", label_leaves = FALSE, label_branch_points = FALSE, cell_size=2, label_cell_groups = F) + ggtitle("Pseudotime_ SC-Skin as root")

p1<-plot_cells(SC_pseudo.cds, label_groups_by_cluster = T,color_cells_by="tissue_origin", label_leaves = FALSE, label_branch_points = FALSE, cell_size=1, label_cell_groups = F)+ scale_color_manual(values=color_datasets) + ggtitle("SC dataset combination")
p2<-plot_cells(SC_pseudo.cds, label_groups_by_cluster = T,color_cells_by="celltype", label_leaves = FALSE, label_branch_points = FALSE, cell_size=1, label_cell_groups = F)+ scale_color_manual(values=color_SC_all_celltype) + ggtitle("SC celltype distribution")
p3<-plot_cells(SC_pseudo.cds, label_groups_by_cluster = T,color_cells_by="pseudotime", label_leaves = FALSE, label_branch_points = FALSE, cell_size=1, label_cell_groups = F) + ggtitle("Pseudotime_ SC-Skin as root")
plot_cells(SC_pseudo.cds, label_groups_by_cluster = T,color_cells_by="seurat_clusters", label_leaves = FALSE, label_branch_points = FALSE, cell_size=1, label_cell_groups = F) + ggtitle("Pseudotime_ SC-Skin as root")
wrap_plots(p1,p2,p3)

#UMAP branchingpoint SCs for cellchat + principal graph
plot_cells(SC_pseudo.cds, label_groups_by_cluster = T,color_cells_by="seurat_clusters", label_leaves = FALSE, label_branch_points = FALSE, cell_size=1, label_cell_groups = F)+ scale_color_manual(values=color_branch) + ggtitle("SC celltype distribution")



#genes high in branching point
Idents(SC_k.i.all)<-"seurat_clusters"
branchingpoint.clustermarker<-FindMarkers(SC_k.i.all,  ident.1= "1", ident.2= c("0","2","3","4","5","6","7"), assay = "RNA", min.pct = 0.01, logfc.threshold = 0.01)
branchingpoint.clustermarker$Foldchange_UP <- 2^(branchingpoint.clustermarker$avg_log2FC)
branchingpoint.clustermarker$Foldchange_DOWN <- 2^(-branchingpoint.clustermarker$avg_log2FC)
branchingpoint.clustermarker$Ratio_pct1_pct2 <- (branchingpoint.clustermarker$pct.1)/(branchingpoint.clustermarker$pct.2)
write.xlsx(branchingpoint.clustermarker, "Lists/DEG_branchingpoint.xlsx")


p1<-plot_cells(SC_pseudo.cds,
               genes=c("JUNB"),           
               label_cell_groups=F,
               show_trajectory_graph=T, cell_size=1, trajectory_graph_color = "black", label_branch_points = F, label_roots = F, label_leaves = F)+scale_color_viridis(option="inferno")

p2<-plot_cells(SC_pseudo.cds,
                genes=c("JUND"),           
                label_cell_groups=F,
                show_trajectory_graph=T, cell_size=1, trajectory_graph_color = "black", label_branch_points = F, label_roots = F, label_leaves = F)+scale_color_viridis(option="inferno")

p3<-plot_cells(SC_pseudo.cds,
                genes=c("FOS"),           
                label_cell_groups=F,
                show_trajectory_graph=T, cell_size=1, trajectory_graph_color = "black", label_branch_points = F, label_roots = F, label_leaves = F)+scale_color_viridis(option="inferno")

p4<-plot_cells(SC_pseudo.cds,
               genes=c("FOSB"),           
               label_cell_groups=F,
               show_trajectory_graph=T, cell_size=1, trajectory_graph_color = "black", label_branch_points = F, label_roots = F, label_leaves = F)+scale_color_viridis(option="inferno")

p5<-plot_cells(SC_pseudo.cds,
               genes=c("KLF2"),           
               label_cell_groups=F,
               show_trajectory_graph=T, cell_size=1, trajectory_graph_color = "black", label_branch_points = F, label_roots = F, label_leaves = F)+scale_color_viridis(option="inferno")

p6<-plot_cells(SC_pseudo.cds,
               genes=c("KLF4"),           
               label_cell_groups=F,
               show_trajectory_graph=T, cell_size=1, trajectory_graph_color = "black", label_branch_points = F, label_roots = F, label_leaves = F)+scale_color_viridis(option="inferno")

p7<-plot_cells(SC_pseudo.cds,
               genes=c("KLF6"),           
               label_cell_groups=F,
               show_trajectory_graph=T, cell_size=1, trajectory_graph_color = "black", label_branch_points = F, label_roots = F, label_leaves = F)+scale_color_viridis(option="inferno")

p8<-plot_cells(SC_pseudo.cds,
                genes=c("KLF10"),           
                label_cell_groups=F,
                show_trajectory_graph=T, cell_size=1, trajectory_graph_color = "black", label_branch_points = F, label_roots = F, label_leaves = F)+scale_color_viridis(option="inferno")

p9<-plot_cells(SC_pseudo.cds,
               genes=c("IER2"),           
               label_cell_groups=F,
               show_trajectory_graph=T, cell_size=1, trajectory_graph_color = "black", label_branch_points = F, label_roots = F, label_leaves = F)+scale_color_viridis(option="inferno")

p10<-plot_cells(SC_pseudo.cds,
                genes=c("IER3"),           
                label_cell_groups=F,
                show_trajectory_graph=T, cell_size=1, trajectory_graph_color = "black", label_branch_points = F, label_roots = F, label_leaves = F)+scale_color_viridis(option="inferno")

p11<-plot_cells(SC_pseudo.cds,
               genes=c("IER5"),           
               label_cell_groups=F,
               show_trajectory_graph=T, cell_size=1, trajectory_graph_color = "black", label_branch_points = F, label_roots = F, label_leaves = F)+scale_color_viridis(option="inferno")

p12<-plot_cells(SC_pseudo.cds,
               genes=c("IER5L"),           
               label_cell_groups=F,
               show_trajectory_graph=T, cell_size=1, trajectory_graph_color = "black", label_branch_points = F, label_roots = F, label_leaves = F)+scale_color_viridis(option="inferno")

p13<-plot_cells(SC_pseudo.cds,
                genes=c("ATF3"),           
                label_cell_groups=F,
                show_trajectory_graph=T, cell_size=1, trajectory_graph_color = "black", label_branch_points = F, label_roots = F, label_leaves = F)+scale_color_viridis(option="inferno")

p14<-plot_cells(SC_pseudo.cds,
                genes=c("EGR1"),           
                label_cell_groups=F,
                show_trajectory_graph=T, cell_size=1, trajectory_graph_color = "black", label_branch_points = F, label_roots = F, label_leaves = F)+scale_color_viridis(option="inferno")

p15<-plot_cells(SC_pseudo.cds,
                genes=c("HES1"),           
                label_cell_groups=F,
                show_trajectory_graph=T, cell_size=1, trajectory_graph_color = "black", label_branch_points = F, label_roots = F, label_leaves = F)+scale_color_viridis(option="inferno")

p16<-plot_cells(SC_pseudo.cds,
                genes=c("ZFP36"),           
                label_cell_groups=F,
                show_trajectory_graph=T, cell_size=1, trajectory_graph_color = "black", label_branch_points = F, label_roots = F, label_leaves = F)+scale_color_viridis(option="inferno")
wrap_plots(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16)


#known repair SC factors
p1<-plot_cells(SC_pseudo.cds,
               genes=c("JUN"),           
               label_cell_groups=F,
               show_trajectory_graph=T, cell_size=1, trajectory_graph_color = "black", label_branch_points = F, label_roots = F, label_leaves = F)+scale_color_viridis(option="inferno")

p2<-plot_cells(SC_pseudo.cds,
               genes=c("STAT3"),           
               label_cell_groups=F,
               show_trajectory_graph=T, cell_size=1, trajectory_graph_color = "black", label_branch_points = F, label_roots = F, label_leaves = F)+scale_color_viridis(option="inferno")

p3<-plot_cells(SC_pseudo.cds,
               genes=c("ARTN"),           
               label_cell_groups=F,
               show_trajectory_graph=T, cell_size=1, trajectory_graph_color = "black", label_branch_points = F, label_roots = F, label_leaves = F)+scale_color_viridis(option="inferno")
p4<-plot_cells(SC_pseudo.cds,
               genes=c("BDNF"),           
               label_cell_groups=F,
               show_trajectory_graph=T, cell_size=1, trajectory_graph_color = "black", label_branch_points = F, label_roots = F, label_leaves = F)+scale_color_viridis(option="inferno")

p5<-plot_cells(SC_pseudo.cds,
               genes=c("GDNF"),           
               label_cell_groups=F,
               show_trajectory_graph=T, cell_size=1, trajectory_graph_color = "black", label_branch_points = F, label_roots = F, label_leaves = F)+scale_color_viridis(option="inferno")

p6<-plot_cells(SC_pseudo.cds,
               genes=c("SHH"),           
               label_cell_groups=F,
               show_trajectory_graph=T, cell_size=1, trajectory_graph_color = "black", label_branch_points = F, label_roots = F, label_leaves = F)+scale_color_viridis(option="inferno")

p7<-plot_cells(SC_pseudo.cds,
               genes=c("OLIG1"),           
               label_cell_groups=F,
               show_trajectory_graph=T, cell_size=1, trajectory_graph_color = "black", label_branch_points = F, label_roots = F, label_leaves = F)+scale_color_viridis(option="inferno")

p8<-plot_cells(SC_pseudo.cds,
           genes=c("IGFBP2"),           
           label_cell_groups=F,
           show_trajectory_graph=T, cell_size=1, trajectory_graph_color = "black", label_branch_points = F, label_roots = F, label_leaves = F)+scale_color_viridis(option="inferno")
wrap_plots(p1,p2,p3,p4,p5,p6,p7,p8, ncol=4)


##### Cellchat ####
# load necessary libraries
options(stringsAsFactors = FALSE)

#Keloid_ cellchat
Idents(SC_k.i.all)="tissue"
SC.cluster.k <- subset(SC_k.i.all, idents = "keloid")
Idents(SC.cluster.k)<- SC.cluster.k$celltype
SC.cluster.k <- subset(SC.cluster.k, idents = c("SC-Myel_XL","SC-Myel+Nonmyel_CCD","SC-Myel+Nonmyel_MD","SC-Nonmyel1_XL","SC-Nonmyel2_XL","SC-Keloid_XL","SC-Keloid_CCD","SC-Keloid_MD"))
Idents(SC.cluster.k)<- SC.cluster.k$celltype
DefaultAssay(SC.cluster.k)<-"RNA"
SC.cluster.k[["percent.mt"]] <- PercentageFeatureSet(SC.cluster.k, pattern = "^MT-")
SC.cluster.k <- SCTransform(SC.cluster.k, method = "glmGamPoi",vars.to.regress = "percent.mt",verbose = F)
SC.cluster.k[['integrated']] <- NULL
Idents(SC.cluster.k)<- SC.cluster.k$celltype


k.i.MD<-k.MD
Idents(k.i.MD)<- k.i.MD$celltype
DefaultAssay(k.i.MD)<-"RNA"
k.i.MD[["percent.mt"]] <- PercentageFeatureSet(k.i.MD, pattern = "^MT-")
k.i.MD <- SCTransform(k.i.MD, method = "glmGamPoi",vars.to.regress = "percent.mt",verbose = F)
k.i.MD[['integrated']] <- NULL

k.i.CCD<-k.CCD
Idents(k.i.CCD)<- k.i.CCD$celltype
DefaultAssay(k.i.CCD)<-"RNA"
k.i.CCD[["percent.mt"]] <- PercentageFeatureSet(k.i.CCD, pattern = "^MT-")
k.i.CCD <- SCTransform(k.i.CCD, method = "glmGamPoi",vars.to.regress = "percent.mt",verbose = F)
k.i.CCD[['integrated']] <- NULL
Idents(k.i.CCD)<- k.i.CCD$celltype

k.i.XL<-k.XL
Idents(k.i.XL)<- k.i.XL$celltype
DefaultAssay(k.i.XL)<-"RNA"
k.i.XL[["percent.mt"]] <- PercentageFeatureSet(k.i.XL, pattern = "^MT-")
k.i.XL <- SCTransform(k.i.XL, method = "glmGamPoi",vars.to.regress = "percent.mt",verbose = F)
k.i.XL[['integrated']] <- NULL
Idents(k.i.XL)<- k.i.XL$celltype

total_k.all.list<-list(k.i.MD,k.i.CCD,k.i.XL)#,SC.cluster.k
total_k.all.features <- SelectIntegrationFeatures(object.list = total_k.all.list, nfeatures = 3000)
total_k.all.list <- PrepSCTIntegration(total_k.all.list, anchor.features = total_k.all.features)
total_k.all.list <- lapply(total_k.all.list, RunPCA, verbose = F, features= total_k.all.features)
total_k.all.anchors<-FindIntegrationAnchors(total_k.all.list,normalization.method = "SCT", anchor.features = total_k.all.features, reduction = "rpca")
total_k.all<-IntegrateData(anchorset=total_k.all.anchors, normalization.method = "SCT", k.weight = 50)
total_k.all <- RunPCA(total_k.all)
ElbowPlot(total_k.all, ndims =30)
total_k.all <- RunUMAP(total_k.all, dims = 1:18)
total_k.all <- FindNeighbors(total_k.all, dims = 1:18)
total_k.all <- FindClusters(total_k.all, resolution = 0.5)
UMAPPlot(total_k.all)
DefaultAssay(total_k.all)<- "RNA"
total_k.all <- NormalizeData(total_k.all)

table(total_k.all$celltype)
Idents(total_k.all)<-total_k.all$celltype
total_k.all<-RenameIdents(total_k.all,
                          `DC`="MAC/DC",`EC`="EC",`FB`="FB",`KC`="KC",`LEC`="LEC",`MAC`="MAC/DC",`MAC/DC`="MAC/DC",`PC`="SMC/PC",
                          `SC`="SC",`SGC`="SGC",`SMC`="SMC/PC",`SMC/PC`="SMC/PC",`TC`="TC")
total_k.all$celltype<-Idents(total_k.all)

#Define cluster levels
clusters_ordered<-c(
  
  "FB", 
  "SMC/PC",
  "KC",
  "EC", "LEC", 
  "TC", "MAC/DC", "SC",
  "SGC")
total_k.all$celltype<- factor(total_k.all$celltype, levels = clusters_ordered)
Idents(total_k.all)<-total_k.all$celltype


###CellChat###
# create test seurat object, set identity/check levels
levels(total_k.all)
Idents(total_k.all)<- "celltype"
keloid_cellchat <- total_k.all
Idents(keloid_cellchat)<- keloid_cellchat$celltype

# Create cell Chat object, use RNA assay as negative values can be found in the integrated data
keloid_cellchat.chat <- createCellChat(keloid_cellchat, group.by = "ident",assay = "RNA")

# load CellChat Database is a manually curated database of literature-supported ligand-receptor interactions
#in both human and mouse. CellChatDB in mouse contains 2,021 validated molecular interactions, including 60%
#of secrete autocrine/paracrine signaling interactions, 21% of extracellular matrix (ECM)-receptor interactions
#and 19% of cell-cell contact interactions. CellChatDB in human contains 1,939 validated molecular interactions,
#including 61.8% of paracrine/autocrine signaling interactions, 21.7% of extracellular matrix (ECM)-receptor interactions
#and 16.5% of cell-cell contact interactions.

#If needed one can manually add ligand receptor pairs to the database, for tutorial see to
#https://htmlpreview.github.io/?https://github.com/sqjin/CellChat/blob/master/tutorial/CellChat-vignette.html#load-the-required-libraries

CellChatDB <- CellChatDB.human
showDatabaseCategory(CellChatDB)
showDatabaseCategory(CellChatDB)

CellChatDB.use <- CellChatDB
keloid_cellchat.chat@DB<- CellChatDB.use

keloid_cellchat.chat<- subsetData(keloid_cellchat.chat)


# identify overexpressed genes
keloid_cellchat.chat <- identifyOverExpressedGenes(keloid_cellchat.chat)
keloid_cellchat.chat <- identifyOverExpressedInteractions(keloid_cellchat.chat)

# compute communication probability
keloid_cellchat.chat <- computeCommunProb(keloid_cellchat.chat)
keloid_cellchat.chat <- filterCommunication(keloid_cellchat.chat, min.cells = 35)

keloid_cellchat.chat <- computeCommunProbPathway(keloid_cellchat.chat)
keloid_cellchat.chat <- aggregateNet(keloid_cellchat.chat)

df.net <- subsetCommunication(keloid_cellchat.chat)
df.netP <- subsetCommunication(keloid_cellchat.chat,slot.name = "netP")
write.xlsx(df.net, "Lists/df_net_keloidcells_SCmix.xlsx")
write.xlsx(df.netP, "Lists/df_netP_keloidcells_SCmix.xlsx")

groupSize <- as.numeric(table(keloid_cellchat.chat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(keloid_cellchat.chat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions", color.use = color_cellchat_keloid_mix)
netVisual_circle(keloid_cellchat.chat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength", color.use = color_cellchat_keloid_mix)

netVisual_aggregate(keloid_cellchat.chat, signaling =c("COLLAGEN"), layout = "chord")

par(mfrow=c(1,1))
netVisual_aggregate(keloid_cellchat.chat, signaling = c("COLLAGEN"), layout = "chord", color.use = color_cellchat_keloid_mix)
par(mfrow=c(1,1))
netVisual_aggregate(keloid_cellchat.chat, signaling = c("TNF"), layout = "chord", color.use = color_cellchat_keloid_mix)

mat <- keloid_cellchat.chat@net$weight
par(mfrow = c(3,6), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i], color.use = color_cellchat_keloid_mix)
}

vertex.receiver = seq(1,4) # a numeric vector. 
netVisual_aggregate(keloid_cellchat.chat, signaling = c("COLLAGEN"),  vertex.receiver = vertex.receiver)


#Keloid_subtype cellchat
Idents(SC_k.i.all)="tissue"
SC.cluster.k <- subset(SC_k.i.all, idents = "keloid")
Idents(SC.cluster.k)<- SC.cluster.k$celltype
DefaultAssay(SC.cluster.k)<-"RNA"
SC.cluster.k[["percent.mt"]] <- PercentageFeatureSet(SC.cluster.k, pattern = "^MT-")
SC.cluster.k <- SCTransform(SC.cluster.k, method = "glmGamPoi",vars.to.regress = "percent.mt",verbose = F)
SC.cluster.k[['integrated']] <- NULL
Idents(SC.cluster.k)<- SC.cluster.k$celltype


k.i.MD<-k.MD
Idents(k.i.MD)<- k.i.MD$celltype
k.i.MD <- subset(k.i.MD, idents = c("FB","SMC/PC","KC","EC","LEC","TC","MAC","DC"))
DefaultAssay(k.i.MD)<-"RNA"
k.i.MD[["percent.mt"]] <- PercentageFeatureSet(k.i.MD, pattern = "^MT-")
k.i.MD <- SCTransform(k.i.MD, method = "glmGamPoi",vars.to.regress = "percent.mt",verbose = F)
k.i.MD[['integrated']] <- NULL

k.i.CCD<-k.CCD
Idents(k.i.CCD)<- k.i.CCD$celltype
k.i.CCD <- subset(k.i.CCD, idents = c("FB","SMC","PC","KC","EC","LEC","MAC/DC","SGC"))
DefaultAssay(k.i.CCD)<-"RNA"
k.i.CCD[["percent.mt"]] <- PercentageFeatureSet(k.i.CCD, pattern = "^MT-")
k.i.CCD <- SCTransform(k.i.CCD, method = "glmGamPoi",vars.to.regress = "percent.mt",verbose = F)
k.i.CCD[['integrated']] <- NULL
Idents(k.i.CCD)<- k.i.CCD$celltype

k.i.XL<-k.XL
Idents(k.i.XL)<- k.i.XL$celltype
k.i.XL <- subset(k.i.XL, idents = c("FB","SMC/PC","KC","EC","LEC","TC","MAC/DC","SGC" ))
DefaultAssay(k.i.XL)<-"RNA"
k.i.XL[["percent.mt"]] <- PercentageFeatureSet(k.i.XL, pattern = "^MT-")
k.i.XL <- SCTransform(k.i.XL, method = "glmGamPoi",vars.to.regress = "percent.mt",verbose = F)
k.i.XL[['integrated']] <- NULL
Idents(k.i.XL)<- k.i.XL$celltype

total_k.all.list<-list(k.i.MD,k.i.CCD,k.i.XL,SC.cluster.k)
total_k.all.features <- SelectIntegrationFeatures(object.list = total_k.all.list, nfeatures = 3000)
total_k.all.list <- PrepSCTIntegration(total_k.all.list, anchor.features = total_k.all.features)
total_k.all.list <- lapply(total_k.all.list, RunPCA, verbose = F, features= total_k.all.features)
total_k.all.anchors<-FindIntegrationAnchors(total_k.all.list,normalization.method = "SCT", anchor.features = total_k.all.features, reduction = "rpca")
total_k.all<-IntegrateData(anchorset=total_k.all.anchors, normalization.method = "SCT", k.weight = 50)
total_k.all <- RunPCA(total_k.all)
ElbowPlot(total_k.all, ndims =30)
total_k.all <- RunUMAP(total_k.all, dims = 1:18)
total_k.all <- FindNeighbors(total_k.all, dims = 1:18)
total_k.all <- FindClusters(total_k.all, resolution = 0.5)
UMAPPlot(total_k.all)
DefaultAssay(total_k.all)<- "RNA"
total_k.all <- NormalizeData(total_k.all)

table(total_k.all$celltype)
backup1<-total_k.all
Idents(total_k.all)<-total_k.all$celltype
total_k.all<-RenameIdents(total_k.all,
                          `DC`="MAC/DC",`EC`="EC",`FB`="FB",`KC`="KC",`LEC`="LEC",`MAC`="MAC/DC",`MAC/DC`="MAC/DC",`PC`="SMC/PC",
                          `SC-EC_MD`="SC-EC",`SC-FB_MD`="SC-FB",`SC-Myel_XL`="SC-Myel+Nonmyel",`SC-Myel+Nonmyel_CCD`="SC-Myel+Nonmyel",`SC-Myel+Nonmyel_MD`="SC-Myel+Nonmyel",`SC-Nonmyel1_XL`="SC-Myel+Nonmyel",`SC-Nonmyel2_XL`="SC-Myel+Nonmyel",`SC-Prolif_MD`="SC-Prolif",
                          `SC-Keloid_CCD`="SC-Keloid",`SC-Keloid_MD`="SC-Keloid",`SC-Keloid_XL`="SC-Keloid",`SGC`="SGC",`SMC`="SMC/PC",`SMC/PC`="SMC/PC",`TC`="TC")
total_k.all$celltype<-Idents(total_k.all)

#Define cluster levels
clusters_ordered<-c(
  "SC-Myel+Nonmyel", "SC-Keloid","SC-Prolif", "SC-EC","SC-FB",
  "FB", 
  "SMC/PC",
  "KC",
  "EC", "LEC", 
  "TC", "MAC/DC", 
  "SGC")
total_k.all$celltype<- factor(total_k.all$celltype, levels = clusters_ordered)
Idents(total_k.all)<-total_k.all$celltype


###CellChat###
# create test seurat object, set identity/check levels
levels(total_k.all)
Idents(total_k.all)<- "celltype"
keloid_cellchat <- total_k.all
Idents(keloid_cellchat)<- keloid_cellchat$celltype

# Create cell Chat object, use RNA assay as negative values can be found in the integrated data
keloid_cellchat.chat <- createCellChat(keloid_cellchat, group.by = "ident",assay = "RNA")

# load CellChat Database is a manually curated database of literature-supported ligand-receptor interactions
#in both human and mouse. CellChatDB in mouse contains 2,021 validated molecular interactions, including 60%
#of secrete autocrine/paracrine signaling interactions, 21% of extracellular matrix (ECM)-receptor interactions
#and 19% of cell-cell contact interactions. CellChatDB in human contains 1,939 validated molecular interactions,
#including 61.8% of paracrine/autocrine signaling interactions, 21.7% of extracellular matrix (ECM)-receptor interactions
#and 16.5% of cell-cell contact interactions.

#If needed one can manually add ligand receptor pairs to the database, for tutorial see to
#https://htmlpreview.github.io/?https://github.com/sqjin/CellChat/blob/master/tutorial/CellChat-vignette.html#load-the-required-libraries

CellChatDB <- CellChatDB.human
showDatabaseCategory(CellChatDB)
showDatabaseCategory(CellChatDB)

CellChatDB.use <- CellChatDB
keloid_cellchat.chat@DB<- CellChatDB.use

keloid_cellchat.chat<- subsetData(keloid_cellchat.chat)


# identify overexpressed genes
keloid_cellchat.chat <- identifyOverExpressedGenes(keloid_cellchat.chat)
keloid_cellchat.chat <- identifyOverExpressedInteractions(keloid_cellchat.chat)

# compute communication probability
keloid_cellchat.chat <- computeCommunProb(keloid_cellchat.chat)
keloid_cellchat.chat <- filterCommunication(keloid_cellchat.chat, min.cells = 35)

keloid_cellchat.chat <- computeCommunProbPathway(keloid_cellchat.chat)
keloid_cellchat.chat <- aggregateNet(keloid_cellchat.chat)

df.net <- subsetCommunication(keloid_cellchat.chat)
df.netP <- subsetCommunication(keloid_cellchat.chat,slot.name = "netP")
write.xlsx(df.net, "Lists/df_net_keloidcells.xlsx")
write.xlsx(df.netP, "Lists/df_netP_keloidcells.xlsx")

groupSize <- as.numeric(table(keloid_cellchat.chat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(keloid_cellchat.chat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions", color.use = color_cellchat_keloid)
netVisual_circle(keloid_cellchat.chat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength", color.use = color_cellchat_keloid)

netVisual_aggregate(keloid_cellchat.chat, signaling =c("COLLAGEN"), layout = "chord")

par(mfrow=c(1,1))
netVisual_aggregate(keloid_cellchat.chat, signaling = c("COLLAGEN"), layout = "chord", color.use = color_cellchat_keloid)
par(mfrow=c(1,1))
netVisual_aggregate(keloid_cellchat.chat, signaling = c("TNF"), layout = "chord", color.use = color_cellchat_keloid)

mat <- keloid_cellchat.chat@net$weight
par(mfrow = c(3,6), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i], color.use = color_cellchat_keloid)
}

vertex.receiver = seq(1,4) # a numeric vector. 
netVisual_aggregate(keloid_cellchat.chat, signaling = c("COLLAGEN"),  vertex.receiver = vertex.receiver)
# Circle plot
par(mfrow=c(1,1))
netVisual_aggregate(keloid_cellchat.chat, signaling = c("COLLAGEN"), layout = "circle")

par(mfrow=c(1,1))
netVisual_heatmap(keloid_cellchat.chat, signaling = c("AGRN"), color.heatmap = "Reds", color.use = color_cellchat_keloid)

par(mfrow=c(1,1))
netVisual_heatmap(keloid_cellchat.chat, signaling = c("BSP"), color.heatmap = "Reds", color.use = color_cellchat_keloid)

par(mfrow=c(1,1))
netVisual_heatmap(keloid_cellchat.chat, signaling = c("CADM"), color.heatmap = "Reds", color.use = color_cellchat_keloid)

par(mfrow=c(1,1))
netVisual_heatmap(keloid_cellchat.chat, signaling = c("CDH"), color.heatmap = "Reds", color.use = color_cellchat_keloid)

par(mfrow=c(1,1))
netVisual_heatmap(keloid_cellchat.chat, signaling = c("COLLAGEN"), color.heatmap = "Reds", color.use = color_cellchat_keloid)

par(mfrow=c(1,1))
netVisual_heatmap(keloid_cellchat.chat, signaling = c("EGF"), color.heatmap = "Reds", color.use = color_cellchat_keloid)

par(mfrow=c(1,1))
netVisual_heatmap(keloid_cellchat.chat, signaling = c("FN1"), color.heatmap = "Reds", color.use = color_cellchat_keloid)

par(mfrow=c(1,1))
netVisual_heatmap(keloid_cellchat.chat, signaling = c("PDGF"), color.heatmap = "Reds", color.use = color_cellchat_keloid)

par(mfrow=c(1,1))
netVisual_heatmap(keloid_cellchat.chat, signaling = c("PERIOSTIN"), color.heatmap = "Reds", color.use = color_cellchat_keloid)

par(mfrow=c(1,1))
netVisual_heatmap(keloid_cellchat.chat, signaling = c("PTN"), color.heatmap = "Reds", color.use = color_cellchat_keloid)

par(mfrow=c(1,1))
netVisual_heatmap(keloid_cellchat.chat, signaling = c("SEMA3"), color.heatmap = "Reds", color.use = color_cellchat_keloid)

par(mfrow=c(1,1))
netVisual_heatmap(keloid_cellchat.chat, signaling = c("TENASCIN"), color.heatmap = "Reds", color.use = color_cellchat_keloid)

par(mfrow=c(1,2))
netVisual_heatmap(keloid_cellchat.chat, signaling = c("LAMININ"), color.heatmap = "Reds", color.use = color_cellchat_keloid)


par(mfrow=c(1,2))
netVisual_heatmap(keloid_cellchat.chat, signaling = c("MIF"), color.heatmap = "Reds", color.use = color_cellchat_keloid)

netVisual_bubble(keloid_cellchat.chat, sources.use = 1, targets.use = c(1:13), color.heatmap = "Spectral" ,remove.isolate = FALSE)

netVisual_bubble(keloid_cellchat.chat, sources.use = 2, targets.use = c(1:13),remove.isolate = FALSE)

netVisual_bubble(keloid_cellchat.chat, sources.use = 3, targets.use = c(1:13),remove.isolate = FALSE)

netVisual_bubble(keloid_cellchat.chat, sources.use = 4, targets.use = c(1:13),remove.isolate = FALSE)

netVisual_bubble(keloid_cellchat.chat, sources.use = 5, targets.use = c(1:13),remove.isolate = FALSE)

netVisual_bubble(keloid_cellchat.chat, sources.use = 6, targets.use = c(1:5), signaling= c("ANGPTL","COLLAGEN","MK","PTN"), color.heatmap = "Spectral", remove.isolate = FALSE)

netVisual_bubble(keloid_cellchat.chat, sources.use = 7, targets.use = c(1:5), signaling= c("COLLAGEN"), remove.isolate = FALSE)


#cellchat_healthy skin mix
total_s.all<-s.TT
table(total_s.all$celltype)
Idents(total_s.all)<-total_s.all$celltype
total_s.all<-RenameIdents(total_s.all,
                          `EC`="EC",`FB`="FB",`KC`="KC",`MAC/DC`="MAC/DC",`MEL`="MEL",
                          `PC`="SMC/PC",`SC`="SC",`SGC`="SGC",`SMC`="SMC/PC",`TC`="TC")
total_s.all$celltype<-Idents(total_s.all)

#Define cluster levels
clusters_ordered<-c(
  "FB", 
  "SMC/PC",
  "KC",
  "EC", 
  "TC", "MAC/DC","SC","MEL", 
  "SGC")
total_s.all$celltype<- factor(total_s.all$celltype, levels = clusters_ordered)
Idents(total_s.all)<-total_s.all$celltype


###CellChat###
# set identity/check levels
levels(total_s.all)
Idents(total_s.all)<- "celltype"
skin_cellchat <- total_s.all
Idents(skin_cellchat)<- skin_cellchat$celltype

# Create cell Chat object, use RNA assay as negative values can be found in the integrated data
skin_cellchat.chat <- createCellChat(skin_cellchat, group.by = "ident",assay = "RNA")

# load CellChat Database is a manually curated database of literature-supported ligand-receptor interactions
#in both human and mouse. CellChatDB in mouse contains 2,021 validated molecular interactions, including 60%
#of secrete autocrine/paracrine signaling interactions, 21% of extracellular matrix (ECM)-receptor interactions
#and 19% of cell-cell contact interactions. CellChatDB in human contains 1,939 validated molecular interactions,
#including 61.8% of paracrine/autocrine signaling interactions, 21.7% of extracellular matrix (ECM)-receptor interactions
#and 16.5% of cell-cell contact interactions.

#If needed one can manually add ligand receptor pairs to the database, for tutorial see to
#https://htmlpreview.github.io/?https://github.com/sqjin/CellChat/blob/master/tutorial/CellChat-vignette.html#load-the-required-libraries

CellChatDB <- CellChatDB.human
showDatabaseCategory(CellChatDB)

CellChatDB.use <- CellChatDB
skin_cellchat.chat@DB<- CellChatDB.use

skin_cellchat.chat<- subsetData(skin_cellchat.chat)


# identify overexpressed genes
skin_cellchat.chat <- identifyOverExpressedGenes(skin_cellchat.chat)
skin_cellchat.chat <- identifyOverExpressedInteractions(skin_cellchat.chat)

# compute communication probability
skin_cellchat.chat <- computeCommunProb(skin_cellchat.chat)
skin_cellchat.chat <- filterCommunication(skin_cellchat.chat, min.cells = 35)

skin_cellchat.chat <- computeCommunProbPathway(skin_cellchat.chat)
skin_cellchat.chat <- aggregateNet(skin_cellchat.chat)

df.net <- subsetCommunication(skin_cellchat.chat)
df.netP <- subsetCommunication(skin_cellchat.chat,slot.name = "netP")
write.xlsx(df.net, "Lists/df_net_skincells_mix.xlsx")
write.xlsx(df.netP, "Lists/df_netP_skincells_mix.xlsx")


groupSize <- as.numeric(table(skin_cellchat.chat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(skin_cellchat.chat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions", color.use = color_cellchat_skin_mix)
netVisual_circle(skin_cellchat.chat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength", color.use = color_cellchat_skin_mix)

netVisual_aggregate(skin_cellchat.chat, signaling =c("COLLAGEN"), layout = "chord")

par(mfrow=c(1,1))
netVisual_aggregate(skin_cellchat.chat, signaling = c("COLLAGEN"), layout = "chord", color.use = color_cellchat_skin_mix)
par(mfrow=c(1,1))
netVisual_aggregate(skin_cellchat.chat, signaling = c("TNF"), layout = "chord", color.use = color_cellchat_skin_mix)

mat <- skin_cellchat.chat@net$weight
par(mfrow = c(3,6), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i], color.use = color_cellchat_skin_mix)
}




#healthy skin cellchat
Idents(SC_k.i.all)="origin"
SC.cluster.s <- subset(SC_k.i.all, idents = "TT")
Idents(SC.cluster.s)<- SC.cluster.s$celltype
DefaultAssay(SC.cluster.s)<-"RNA"
SC.cluster.s[["percent.mt"]] <- PercentageFeatureSet(SC.cluster.s, pattern = "^MT-")
SC.cluster.s <- SCTransform(SC.cluster.s, method = "glmGamPoi",vars.to.regress = "percent.mt",verbose = F)
SC.cluster.s[['integrated']] <- NULL
Idents(SC.cluster.s)<- SC.cluster.s$celltype


s.i.TT<-s.TT
Idents(s.i.TT)<- s.i.TT$celltype
s.i.TT <- subset(s.i.TT, idents = c("FB","SMC","PC","KC","EC","TC","MAC/DC","MEL","SGC"))
DefaultAssay(s.i.TT)<-"RNA"
s.i.TT[["percent.mt"]] <- PercentageFeatureSet(s.i.TT, pattern = "^MT-")
s.i.TT <- SCTransform(s.i.TT, method = "glmGamPoi",vars.to.regress = "percent.mt",verbose = F)
s.i.TT[['integrated']] <- NULL

total_s.all.list<-list(s.i.TT,SC.cluster.s)
total_s.all.features <- SelectIntegrationFeatures(object.list = total_s.all.list, nfeatures = 3000)
total_s.all.list <- PrepSCTIntegration(total_s.all.list, anchor.features = total_s.all.features)
total_s.all.list <- lapply(total_s.all.list, RunPCA, verbose = F, features= total_s.all.features)
total_s.all.anchors<-FindIntegrationAnchors(total_s.all.list,normalization.method = "SCT", anchor.features = total_s.all.features, reduction = "rpca")
total_s.all<-IntegrateData(anchorset=total_s.all.anchors, normalization.method = "SCT", k.weight = 50)
total_s.all <- RunPCA(total_s.all)
ElbowPlot(total_s.all, ndims =30)
total_s.all <- RunUMAP(total_s.all, dims = 1:18)
total_s.all <- FindNeighbors(total_s.all, dims = 1:18)
total_s.all <- FindClusters(total_s.all, resolution = 0.5)
UMAPPlot(total_s.all)
DefaultAssay(total_s.all)<- "RNA"
total_s.all <- NormalizeData(total_s.all)

table(total_s.all$celltype)
Idents(total_s.all)<-total_s.all$celltype
#as SC-FB_skin_TT were supposed to be FBs, in the cellchat analysis they were assigned as FBs
total_s.all<-RenameIdents(total_s.all,
                          `EC`="EC",`FB`="FB",`KC`="KC",`MAC/DC`="MAC/DC",`MEL`="MEL",
                          `PC`="SMC/PC",`SC-FB_skin_TT`="FB",`SC-Myel+Nonmyel_skin_TT`="SC-Myel+Nonmyel",`SGC`="SGC",`SMC`="SMC/PC",`TC`="TC")
total_s.all$celltype<-Idents(total_s.all)

#Define cluster levels
clusters_ordered<-c(
  "SC-Myel+Nonmyel",
  "FB", 
  "SMC/PC",
  "KC",
  "EC", 
  "TC", "MAC/DC","MEL", 
  "SGC")
total_s.all$celltype<- factor(total_s.all$celltype, levels = clusters_ordered)
Idents(total_s.all)<-total_s.all$celltype


###CellChat###
# set identity/check levels
levels(total_s.all)
Idents(total_s.all)<- "celltype"
skin_cellchat <- total_s.all
Idents(skin_cellchat)<- skin_cellchat$celltype

# Create cell Chat object, use RNA assay as negative values can be found in the integrated data
skin_cellchat.chat <- createCellChat(skin_cellchat, group.by = "ident",assay = "RNA")

# load CellChat Database is a manually curated database of literature-supported ligand-receptor interactions
#in both human and mouse. CellChatDB in mouse contains 2,021 validated molecular interactions, including 60%
#of secrete autocrine/paracrine signaling interactions, 21% of extracellular matrix (ECM)-receptor interactions
#and 19% of cell-cell contact interactions. CellChatDB in human contains 1,939 validated molecular interactions,
#including 61.8% of paracrine/autocrine signaling interactions, 21.7% of extracellular matrix (ECM)-receptor interactions
#and 16.5% of cell-cell contact interactions.

#If needed one can manually add ligand receptor pairs to the database, for tutorial see to
#https://htmlpreview.github.io/?https://github.com/sqjin/CellChat/blob/master/tutorial/CellChat-vignette.html#load-the-required-libraries

CellChatDB <- CellChatDB.human
showDatabaseCategory(CellChatDB)

CellChatDB.use <- CellChatDB
skin_cellchat.chat@DB<- CellChatDB.use

skin_cellchat.chat<- subsetData(skin_cellchat.chat)


# identify overexpressed genes
skin_cellchat.chat <- identifyOverExpressedGenes(skin_cellchat.chat)
skin_cellchat.chat <- identifyOverExpressedInteractions(skin_cellchat.chat)

# compute communication probability
skin_cellchat.chat <- computeCommunProb(skin_cellchat.chat)
skin_cellchat.chat <- filterCommunication(skin_cellchat.chat, min.cells = 35)

skin_cellchat.chat <- computeCommunProbPathway(skin_cellchat.chat)
skin_cellchat.chat <- aggregateNet(skin_cellchat.chat)

df.net <- subsetCommunication(skin_cellchat.chat)
df.netP <- subsetCommunication(skin_cellchat.chat,slot.name = "netP")
write.xlsx(df.net, "Lists/df_net_skincells.xlsx")
write.xlsx(df.netP, "Lists/df_netP_skincells.xlsx")


groupSize <- as.numeric(table(skin_cellchat.chat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(skin_cellchat.chat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions", color.use = color_cellchat_skin)
netVisual_circle(skin_cellchat.chat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength", color.use = color_cellchat_skin)

netVisual_aggregate(skin_cellchat.chat, signaling =c("COLLAGEN"), layout = "chord")

par(mfrow=c(1,1))
netVisual_aggregate(skin_cellchat.chat, signaling = c("COLLAGEN"), layout = "chord", color.use = color_cellchat_skin)
par(mfrow=c(1,1))
netVisual_aggregate(skin_cellchat.chat, signaling = c("TNF"), layout = "chord", color.use = color_cellchat_skin)

mat <- skin_cellchat.chat@net$weight
par(mfrow = c(3,6), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i], color.use = color_cellchat_skin)
}

vertex.receiver = seq(1,4) # a numeric vector. 
netVisual_aggregate(skin_cellchat.chat, signaling = c("COLLAGEN"),  vertex.receiver = vertex.receiver)
# Circle plot
par(mfrow=c(1,1))
netVisual_aggregate(skin_cellchat.chat, signaling = c("COLLAGEN"), layout = "circle")

par(mfrow=c(1,1))
netVisual_heatmap(skin_cellchat.chat, signaling = c("AGRN"), color.heatmap = "Reds", color.use = color_cellchat_skin)

par(mfrow=c(1,1))
netVisual_heatmap(skin_cellchat.chat, signaling = c("BSP"), color.heatmap = "Reds", color.use = color_cellchat_skin)

par(mfrow=c(1,1))
netVisual_heatmap(skin_cellchat.chat, signaling = c("CADM"), color.heatmap = "Reds", color.use = color_cellchat_skin) 

par(mfrow=c(1,1))
netVisual_heatmap(skin_cellchat.chat, signaling = c("CDH"), color.heatmap = "Reds", color.use = color_cellchat_skin)

par(mfrow=c(1,1))
netVisual_heatmap(skin_cellchat.chat, signaling = c("COLLAGEN"), color.heatmap = "Reds", color.use = color_cellchat_skin)

par(mfrow=c(1,1))
netVisual_heatmap(skin_cellchat.chat, signaling = c("EGF"), color.heatmap = "Reds", color.use = color_cellchat_skin)

par(mfrow=c(1,1))
netVisual_heatmap(skin_cellchat.chat, signaling = c("FN1"), color.heatmap = "Reds", color.use = color_cellchat_skin)

par(mfrow=c(1,1))
netVisual_heatmap(skin_cellchat.chat, signaling = c("PDGF"), color.heatmap = "Reds", color.use = color_cellchat_skin)

par(mfrow=c(1,1))
netVisual_heatmap(skin_cellchat.chat, signaling = c("PERIOSTIN"), color.heatmap = "Reds", color.use = color_cellchat_skin)

par(mfrow=c(1,1))
netVisual_heatmap(skin_cellchat.chat, signaling = c("PTN"), color.heatmap = "Reds", color.use = color_cellchat_skin)

par(mfrow=c(1,1))
netVisual_heatmap(skin_cellchat.chat, signaling = c("SEMA3"), color.heatmap = "Reds", color.use = color_cellchat_skin)

par(mfrow=c(1,1))
netVisual_heatmap(skin_cellchat.chat, signaling = c("TENASCIN"), color.heatmap = "Reds", color.use = color_cellchat_skin)

par(mfrow=c(1,2))
netVisual_heatmap(skin_cellchat.chat, signaling = c("LAMININ"), color.heatmap = "Reds", color.use = color_cellchat_skin)

par(mfrow=c(1,2))
netVisual_heatmap(skin_cellchat.chat, signaling = c("MIF"), color.heatmap = "Reds", color.use = color_cellchat_skin)


par(mfrow=c(1,2))
netVisual_heatmap(skin_cellchat.chat, signaling = c("L1CAM"), color.heatmap = "Reds", color.use = color_cellchat_skin)

netVisual_bubble(skin_cellchat.chat, sources.use = 1, targets.use = c(1:9), color.heatmap = "Spectral" ,remove.isolate = FALSE)

netVisual_bubble(skin_cellchat.chat, sources.use = 2, targets.use = c(1),remove.isolate = FALSE)

netVisual_bubble(skin_cellchat.chat, sources.use = 3, targets.use = c(1),remove.isolate = FALSE)



#SC-branching point_subtype cellchat
Idents(SC_k.i.all)="tissue"
SC.cluster.k <- subset(SC_k.i.all, idents = "keloid")
Idents(SC.cluster.k)="seurat_clusters"
SC.cluster.k <- subset(SC.cluster.k, idents = "1")
Idents(SC.cluster.k)<- SC.cluster.k$celltype
DefaultAssay(SC.cluster.k)<-"RNA"
SC.cluster.k[["percent.mt"]] <- PercentageFeatureSet(SC.cluster.k, pattern = "^MT-")
SC.cluster.k <- SCTransform(SC.cluster.k, method = "glmGamPoi",vars.to.regress = "percent.mt",verbose = F)
SC.cluster.k[['integrated']] <- NULL
Idents(SC.cluster.k)<- SC.cluster.k$celltype


k.i.MD<-k.MD
Idents(k.i.MD)<- k.i.MD$celltype
k.i.MD <- subset(k.i.MD, idents = c("FB","SMC/PC","KC","EC","LEC","TC","MAC","DC"))
DefaultAssay(k.i.MD)<-"RNA"
k.i.MD[["percent.mt"]] <- PercentageFeatureSet(k.i.MD, pattern = "^MT-")
k.i.MD <- SCTransform(k.i.MD, method = "glmGamPoi",vars.to.regress = "percent.mt",verbose = F)
k.i.MD[['integrated']] <- NULL

k.i.CCD<-k.CCD
Idents(k.i.CCD)<- k.i.CCD$celltype
k.i.CCD <- subset(k.i.CCD, idents = c("FB","SMC","PC","KC","EC","LEC","MAC/DC","SGC"))
DefaultAssay(k.i.CCD)<-"RNA"
k.i.CCD[["percent.mt"]] <- PercentageFeatureSet(k.i.CCD, pattern = "^MT-")
k.i.CCD <- SCTransform(k.i.CCD, method = "glmGamPoi",vars.to.regress = "percent.mt",verbose = F)
k.i.CCD[['integrated']] <- NULL
Idents(k.i.CCD)<- k.i.CCD$celltype

k.i.XL<-k.XL
Idents(k.i.XL)<- k.i.XL$celltype
k.i.XL <- subset(k.i.XL, idents = c("FB","SMC/PC","KC","EC","LEC","TC","MAC/DC","SGC" ))
DefaultAssay(k.i.XL)<-"RNA"
k.i.XL[["percent.mt"]] <- PercentageFeatureSet(k.i.XL, pattern = "^MT-")
k.i.XL <- SCTransform(k.i.XL, method = "glmGamPoi",vars.to.regress = "percent.mt",verbose = F)
k.i.XL[['integrated']] <- NULL
Idents(k.i.XL)<- k.i.XL$celltype

total_kb.all.list<-list(k.i.MD,k.i.CCD,k.i.XL,SC.cluster.k)
total_kb.all.features <- SelectIntegrationFeatures(object.list = total_kb.all.list, nfeatures = 3000)
total_kb.all.list <- PrepSCTIntegration(total_kb.all.list, anchor.features = total_kb.all.features)
total_kb.all.list <- lapply(total_kb.all.list, RunPCA, verbose = F, features= total_kb.all.features)
total_kb.all.anchors<-FindIntegrationAnchors(total_kb.all.list,normalization.method = "SCT", anchor.features = total_kb.all.features, reduction = "rpca")
total_kb.all<-IntegrateData(anchorset=total_kb.all.anchors, normalization.method = "SCT", k.weight = 50)
total_kb.all <- RunPCA(total_kb.all)
ElbowPlot(total_kb.all, ndims =30)
total_kb.all <- RunUMAP(total_kb.all, dims = 1:18)
total_kb.all <- FindNeighbors(total_kb.all, dims = 1:18)
total_kb.all <- FindClusters(total_kb.all, resolution = 0.5)
UMAPPlot(total_kb.all)
DefaultAssay(total_kb.all)<- "RNA"
total_kb.all <- NormalizeData(total_kb.all)

table(total_kb.all$celltype)
Idents(total_kb.all)<-total_kb.all$celltype
total_kb.all<-RenameIdents(total_kb.all,
                           `DC`="MAC/DC",`EC`="EC",`FB`="FB",`KC`="KC",`LEC`="LEC",`MAC`="MAC/DC",`MAC/DC`="MAC/DC",`PC`="SMC/PC",
                           `SC-EC_MD`="SC-branch", `SC-FB_MD`="SC-branch", `SC-Myel+Nonmyel_CCD`="SC-branch", `SC-Myel+Nonmyel_MD`="SC-branch", `SC-Nonmyel1_XL`="SC-branch", 
                           `SC-Nonmyel2_XL`="SC-branch", `SC-Prolif_MD`="SC-branch", `SC-Keloid_MD`="SC-branch", `SC-Keloid_XL`="SC-branch",
                           `SGC`="SGC",`SMC`="SMC/PC",`SMC/PC`="SMC/PC",`TC`="TC")
total_kb.all$celltype<-Idents(total_kb.all)

#Define cluster levels
clusters_ordered<-c(
  "FB", 
  "SMC/PC",
  "KC",
  "EC", "LEC", 
  "TC", "MAC/DC","SC-branch", 
  "SGC")
total_kb.all$celltype<- factor(total_kb.all$celltype, levels = clusters_ordered)
Idents(total_kb.all)<-total_kb.all$celltype


###CellChat###
# create test seurat object, set identity/check levels
levels(total_kb.all)
Idents(total_kb.all)<- "celltype"
keloid_cellchat <- total_kb.all
Idents(keloid_cellchat)<- keloid_cellchat$celltype

# Create cell Chat object, use RNA assay as negative values can be found in the integrated data
keloid_cellchat.chat <- createCellChat(keloid_cellchat, group.by = "ident",assay = "RNA")

# load CellChat Database is a manually curated database of literature-supported ligand-receptor interactions
#in both human and mouse. CellChatDB in mouse contains 2,021 validated molecular interactions, including 60%
#of secrete autocrine/paracrine signaling interactions, 21% of extracellular matrix (ECM)-receptor interactions
#and 19% of cell-cell contact interactions. CellChatDB in human contains 1,939 validated molecular interactions,
#including 61.8% of paracrine/autocrine signaling interactions, 21.7% of extracellular matrix (ECM)-receptor interactions
#and 16.5% of cell-cell contact interactions.

#If needed one can manually add ligand receptor pairs to the database, for tutorial see to
#https://htmlpreview.github.io/?https://github.com/sqjin/CellChat/blob/master/tutorial/CellChat-vignette.html#load-the-required-libraries

CellChatDB <- CellChatDB.human
showDatabaseCategory(CellChatDB)
showDatabaseCategory(CellChatDB)

CellChatDB.use <- CellChatDB
keloid_cellchat.chat@DB<- CellChatDB.use

keloid_cellchat.chat<- subsetData(keloid_cellchat.chat)


# identify overexpressed genes
keloid_cellchat.chat <- identifyOverExpressedGenes(keloid_cellchat.chat)
keloid_cellchat.chat <- identifyOverExpressedInteractions(keloid_cellchat.chat)

# compute communication probability
keloid_cellchat.chat <- computeCommunProb(keloid_cellchat.chat)
keloid_cellchat.chat <- filterCommunication(keloid_cellchat.chat, min.cells = 35)

keloid_cellchat.chat <- computeCommunProbPathway(keloid_cellchat.chat)
keloid_cellchat.chat <- aggregateNet(keloid_cellchat.chat)

df.net <- subsetCommunication(keloid_cellchat.chat)
df.netP <- subsetCommunication(keloid_cellchat.chat,slot.name = "netP")
write.xlsx(df.net, "Lists/df_net_keloidcells_branch.xlsx")
write.xlsx(df.netP, "Lists/df_netP_keloidcells_branch.xlsx")

groupSize <- as.numeric(table(keloid_cellchat.chat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(keloid_cellchat.chat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions", color.use = color_cellchat_branch)
netVisual_circle(keloid_cellchat.chat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength", color.use = color_cellchat_branch)

netVisual_aggregate(keloid_cellchat.chat, signaling =c("COLLAGEN"), layout = "chord")

par(mfrow=c(1,1))
netVisual_aggregate(keloid_cellchat.chat, signaling = c("COLLAGEN"), layout = "chord", color.use = color_cellchat_branch)
par(mfrow=c(1,1))
netVisual_aggregate(keloid_cellchat.chat, signaling = c("TNF"), layout = "chord", color.use = color_cellchat_branch)

mat <- keloid_cellchat.chat@net$weight
par(mfrow = c(3,6), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i], color.use = color_cellchat_branch)
}


netVisual_bubble(keloid_cellchat.chat, sources.use = c(1:9), targets.use = c(8),remove.isolate = FALSE)+coord_flip()


#compute DEG Keloid vs. myel/nonmyel
Idents(SC_s.XL)<-"celltype"
s.XL.clustermarker<-FindMarkers(SC_s.XL,  ident.1= "SC-Keloid_skin_XL", ident.2= c("SC-Myel_skin_XL", "SC-Nonmyel_skin_XL"), assay = "RNA", min.pct = 0.01, logfc.threshold = 0.01)
s.XL.clustermarker$Foldchange_UP <- 2^(s.XL.clustermarker$avg_log2FC)
s.XL.clustermarker$Foldchange_DOWN <- 2^(-s.XL.clustermarker$avg_log2FC)
s.XL.clustermarker$Ratio_pct1_pct2 <- (s.XL.clustermarker$pct.1)/(s.XL.clustermarker$pct.2)
write.xlsx(s.XL.clustermarker, "Lists/DEG_repvsmyelnonmyel_s_XL.xlsx")

Idents(SC_k.XL)<-"celltype"
k.XL.clustermarker<-FindMarkers(SC_k.XL,  ident.1= "SC-Keloid_XL", ident.2= c("SC-Myel_XL","SC-Nonmyel1_XL","SC-Nonmyel2_XL"), assay = "RNA", min.pct = 0.01, logfc.threshold = 0.01)
k.XL.clustermarker$Foldchange_UP <- 2^(k.XL.clustermarker$avg_log2FC)
k.XL.clustermarker$Foldchange_DOWN <- 2^(-k.XL.clustermarker$avg_log2FC)
k.XL.clustermarker$Ratio_pct1_pct2 <- (k.XL.clustermarker$pct.1)/(k.XL.clustermarker$pct.2)
write.xlsx(k.XL.clustermarker, "Lists/DEG_repvsmyelnonmyel_k_XL.xlsx")

Idents(SC_k.CCD)<-"celltype"
k.CCD.clustermarker<-FindMarkers(SC_k.CCD,  ident.1= "SC-Keloid_CCD", ident.2= c("SC-Myel+Nonmyel_CCD"), assay = "RNA", min.pct = 0.01, logfc.threshold = 0.01)
k.CCD.clustermarker$Foldchange_UP <- 2^(k.CCD.clustermarker$avg_log2FC)
k.CCD.clustermarker$Foldchange_DOWN <- 2^(-k.CCD.clustermarker$avg_log2FC)
k.CCD.clustermarker$Ratio_pct1_pct2 <- (k.CCD.clustermarker$pct.1)/(k.CCD.clustermarker$pct.2)
write.xlsx(k.CCD.clustermarker, "Lists/DEG_repvsmyelnonmyel_k_CCD.xlsx")

Idents(SC_k.MD)<-"celltype"
k.MD.clustermarker<-FindMarkers(SC_k.MD,  ident.1= "SC-Keloid_MD", ident.2= c("SC-Myel+Nonmyel_MD"), assay = "RNA", min.pct = 0.01, logfc.threshold = 0.01)
k.MD.clustermarker$Foldchange_UP <- 2^(k.MD.clustermarker$avg_log2FC)
k.MD.clustermarker$Foldchange_DOWN <- 2^(-k.MD.clustermarker$avg_log2FC)
k.MD.clustermarker$Ratio_pct1_pct2 <- (k.MD.clustermarker$pct.1)/(k.MD.clustermarker$pct.2)
write.xlsx(k.MD.clustermarker, "Lists/DEG_repvsmyelnonmyel_k_MD.xlsx")
