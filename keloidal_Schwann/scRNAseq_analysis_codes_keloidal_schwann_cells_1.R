#this script accompanies the publication "The transcriptional profile of keloidal Schwann cells"
#code was written by Dr. Martin Direder
#martin.direder@meduniwien.ac.at

setwd(dir ="..")


#load packages and define color schemes
{ library(Seurat)
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
  library(ggpubr)
  library(ggsignif)
  library(limma)
  library(Matrix.utils)
  library(rstatix)
  single_col <- as.character("#02025f")
  color_UMAP_short <- as.character(c("#990F26", "#B33E52","#CC7A88", 
                                     "#99600F",   
                                     "#54990F", "#78B33E",
                                     "#45b6fe", "#296d98", 
                                     "#F9A602", "#967ACC", 
                                     "#666666", "#ff6600", "#333333", "#999999", "#CCCCCC","#3D0F99"))
  color_SC <- as.character(c("#344d0e","#68991c","#FFBF00","#1F78C8","#C8308C","#FF83FA","#fcae1e"))
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


##### integrate all data####
s.k.total.list<-list(keloid1.sct,keloid2.sct, keloid3L.sct,keloid3R.sct,keloid4.sct,keloid5.sct,keloid6.sct,keloid7.sct,keloid8.sct,keloid9.sct,keloid10.sct,skincl1.sct, skincl2.sct, skincl3.sct, skincl4.sct, skin1.sct, skin2.sct, skin3.sct, skin4.sct, skin5.sct, skin6.sct, nscar1.sct, nscar2.sct, nscar3.sct, nscar4.sct, nscar5.sct, nscar6.sct)
s.k.total.features <- SelectIntegrationFeatures(object.list = s.k.total.list, nfeatures = 500)
s.k.total.list <- PrepSCTIntegration(s.k.total.list, anchor.features = s.k.total.features)
s.k.total.list <- lapply(s.k.total.list, RunPCA, verbose = F, features= s.k.total.features)
s.k.total.anchors<-FindIntegrationAnchors(s.k.total.list,normalization.method = "SCT", anchor.features = s.k.total.features, reduction = "rpca")
s.k.total.x <- IntegrateData(anchorset=s.k.total.anchors, normalization.method = "SCT", dims = 1:50)
s.k.total.y <- s.k.total.x
s.k.total.y <- RunPCA(s.k.total.y, npcs = 60)
ElbowPlot(s.k.total.y, ndims = 30)
s.k.total.y <- RunUMAP(s.k.total.y, dims = 1:30)
s.k.total.y <- FindNeighbors(s.k.total.y, dims = 1:30)
s.k.total.y <- FindClusters(s.k.total.y, resolution = 0.6)
UMAPPlot(s.k.total.y, label=T)
s.k.total<-s.k.total.y
DefaultAssay(s.k.total)<- "RNA"
s.k.total <- NormalizeData(s.k.total)

#order conditions
{ s.k.total$tissue <- factor(x = s.k.total$tissue, levels = c("skin", "nscar", "keloid"))
  s.k.total$sample <- factor(x = s.k.total$sample, levels = c("skin_sep1","skin_sep2","skin_sep3", "skin_sep4", "skin_sep5", "skin_sep6", "skin_adj1","skin_adj2","skin_adj3", "skin_adj4", "nscar_1","nscar_2","nscar_3","nscar_4", "nscar_5", "nscar_6", "keloid_mix1", "keloid_mix2", "keloid_mix3","keloid_centr1","keloid_centr2","keloid_centr3","keloid_centr4","keloid_centr5","keloid_centr6","keloid_centr7","keloid_centr8"))
  s.k.total$origin <- factor(x = s.k.total$origin, levels = c("MD", "XL","TT","CCD"))
  s.k.total$location <- factor(x = s.k.total$location, levels = c("skin_sep", "skin_adj", "nscar_centr","keloid_mix","keloid_centr"))}

#Name Cluster
Idents(s.k.total)<-s.k.total$location
s.k.total<-RenameIdents(s.k.total,
                        `skin_sep`="SKIN", `skin_adj`="adjacent SKIN", `nscar_centr`="NORMAL SCAR",   `keloid_mix`="KELOID", `keloid_centr`="KELOID")
s.k.total$area<-Idents(s.k.total)


#Define cluster levels
clusters_ordered<-c("SKIN","adjacent SKIN","NORMAL SCAR","KELOID")
s.k.total$area<- factor(s.k.total$area, levels = clusters_ordered)
Idents(s.k.total)<-s.k.total$celltype

#UMAPPlots, unnamed clusters
UMAPPlot(s.k.total, label=T)
UMAPPlot(s.k.total, label=F, split.by="tissue")
UMAPPlot(s.k.total, label=F, split.by="sample", ncol = 3)
UMAPPlot(s.k.total, label=F, split.by="origin", ncol = 3)
UMAPPlot(s.k.total, label=F, split.by="location", ncol = 3)
UMAPPlot(s.k.total, label=F, split.by="area", ncol = 3)


#Clusteridentification
{
  DotPlot(s.k.total, features=c("PDGFRA","LUM","COL1A1", "DCN","FBLN1","ACTA2","RGS5","KRT10","KRT1","KRT14","KRT5","SELE","VWF","LYVE1","CD3D","CD2","CXCR4","CD79A","MS4A1","CD68","AIF1","FCER1A","S100B","NGFR","PMEL","MLANA","HBB","HBA1","TOP2A","MKI67", "CD20","CD19"), group.by = "seurat_clusters", assay = "RNA") + coord_flip()+scale_color_gradient(low="grey", "high"=single_col)+theme(axis.text.x = element_text(angle = 90)) + ggtitle("Basic Cluster identification")
  DefaultAssay(s.k.total)<- "RNA"
  
  #Name Cluster
  Idents(s.k.total)<-s.k.total$seurat_clusters
  s.k.total<-RenameIdents(s.k.total,
                          `0`="FB", `1`="EC", `2`="FB", `3`="FB", `4`="PC", 
                          `5`="KC", `6`="EC", `7`="KC", `8`="KC", `9`="MAC/DC", `10`="TC",
                          `11`="KC", `12`="EC", `13`="KC", `14`="LEC", `15`="SC/MEL", `16`="FB", `17`="EC", `18`="SMC", `19`="KC")
  s.k.total$celltype<-Idents(s.k.total)
  
  
  #Define cluster levels
  clusters_ordered<-c(
    "FB", 
    "SMC","PC",
    "KC",
    "EC", "LEC", 
    "TC", "MAC/DC", 
    "SC/MEL")
  s.k.total$celltype<- factor(s.k.total$celltype, levels = clusters_ordered)
  Idents(s.k.total)<-s.k.total$celltype
}

#Marker Celltype specific
DotPlot(s.k.total, features=c("PDGFRA","LUM","COL1A1", "DCN","FBLN1","ACTA2","RGS5","KRT10","KRT1","KRT14","KRT5","SELE","VWF","LYVE1","CD3D","CD2","CXCR4","CD68","AIF1","FCER1A","S100B","NGFR","PMEL","MLANA"), group.by = "celltype", assay = "RNA") + coord_flip()+scale_color_gradient(low="grey", "high"=single_col)+theme(axis.text.x = element_text(angle = 90)) + ggtitle("Basic Cluster identification")

#UMAPPlots, named clusters
DefaultAssay(s.k.total)<-"integrated"
Idents(s.k.total)<-"celltype"

UMAPPlot(s.k.total, label=T, label.color = "black",raster=F)+scale_color_manual(values=color_UMAP_short)+ ggtitle("Basic-UMAP")
UMAPPlot(s.k.total, label=F, split.by="area", ncol = 4, raster=F)+ scale_color_manual(values=color_UMAP_short)


#Clustermarker total celltype
s.k.total.clustermarker<-FindAllMarkers(s.k.total, assay = "RNA", min.pct = 0.25, logfc.threshold = 0.25)
s.k.total.clustermarker$Foldchange_UP <- 2^(s.k.total.clustermarker$avg_log2FC)
s.k.total.clustermarker$Foldchange_DOWN <- 2^(-s.k.total.clustermarker$avg_log2FC)
s.k.total.clustermarker$Ratio_pct1_pct2 <- (s.k.total.clustermarker$pct.1)/(s.k.total.clustermarker$pct.2)
write.xlsx(s.k.total.clustermarker, "../Clustermarker_total.xlsx")

#Percentage total
###Pieplot cellnumber /area/cluster
#subset conditions
skin.subset<-subset(s.k.total, area=="SKIN")
adjskin.subset<-subset(s.k.total, area=="adjacent SKIN")
nscar.subset<-subset(s.k.total, area=="NORMAL SCAR")
keloid.subset<-subset(s.k.total, area=="KELOID")

#Data frame frequencies
skin.cluster.frequencies<-as.data.frame(table(skin.subset$celltype))
adjskin.cluster.frequencies<-as.data.frame(table(adjskin.subset$celltype))
nscar.cluster.frequencies<-as.data.frame(table(nscar.subset$celltype))
keloid.cluster.frequencies<-as.data.frame(table(keloid.subset$celltype))

#Barplot Cluster percentage
skin.cluster.frequencies$area="SKIN"
adjskin.cluster.frequencies$area="adjacent SKIN"
nscar.cluster.frequencies$area="NORMAL SCAR"
keloid.cluster.frequencies$area="KELOID"

skin.cluster.frequencies<-mutate(skin.cluster.frequencies,Percentage= skin.cluster.frequencies$Freq/sum(skin.cluster.frequencies$Freq)*100)
adjskin.cluster.frequencies<-mutate(adjskin.cluster.frequencies,Percentage= adjskin.cluster.frequencies$Freq/sum(adjskin.cluster.frequencies$Freq)*100)
nscar.cluster.frequencies<-mutate(nscar.cluster.frequencies,Percentage= nscar.cluster.frequencies$Freq/sum(nscar.cluster.frequencies$Freq)*100)
keloid.cluster.frequencies<-mutate(keloid.cluster.frequencies,Percentage= keloid.cluster.frequencies$Freq/sum(keloid.cluster.frequencies$Freq)*100)

skin_kel.cluster.df <- rbind(skin.cluster.frequencies,adjskin.cluster.frequencies,nscar.cluster.frequencies,keloid.cluster.frequencies)
skin_kel.cluster.df$area <- factor(x = skin_kel.cluster.df$area, levels = c("SKIN","adjacent SKIN", "NORMAL SCAR", "KELOID"))

pdf("Graphs/Main/PiePlot_area_skin_kel_cellfrequencie.pdf", width = 10, height = 10)
pie(skin.cluster.frequencies$Freq, labels = skin.cluster.frequencies$Var1, main="SKIN",  clockwise = T, col = c(color_UMAP_short))
pie(adjskin.cluster.frequencies$Freq, labels = adjskin.cluster.frequencies$Var1, main="adjacent SKIN",  clockwise = T, col = c(color_UMAP_short))
pie(nscar.cluster.frequencies$Freq, labels = nscar.cluster.frequencies$Var1, main="NORMAL SCAR",  clockwise = T, col = c(color_UMAP_short))
pie(keloid.cluster.frequencies$Freq, labels = keloid.cluster.frequencies$Var1, main="KELOID",  clockwise = T, col = c(color_UMAP_short))
dev.off()

#Percentage total
#subset conditions
skins1.subset<-subset(s.k.total, sample=="skin_sep1")
skins2.subset<-subset(s.k.total, sample=="skin_sep2")
skins3.subset<-subset(s.k.total, sample=="skin_sep3")
skins4.subset<-subset(s.k.total, sample=="skin_sep4")
skins5.subset<-subset(s.k.total, sample=="skin_sep5")
skins6.subset<-subset(s.k.total, sample=="skin_sep6")
skina1.subset<-subset(s.k.total, sample=="skin_adj1")
skina2.subset<-subset(s.k.total, sample=="skin_adj2")
skina3.subset<-subset(s.k.total, sample=="skin_adj3")
skina4.subset<-subset(s.k.total, sample=="skin_adj4")
nscar1.subset<-subset(s.k.total, sample=="nscar_1")
nscar2.subset<-subset(s.k.total, sample=="nscar_2")
nscar3.subset<-subset(s.k.total, sample=="nscar_3")
nscar4.subset<-subset(s.k.total, sample=="nscar_4")
nscar5.subset<-subset(s.k.total, sample=="nscar_5")
nscar6.subset<-subset(s.k.total, sample=="nscar_6")
keloidm1.subset<-subset(s.k.total, sample=="keloid_mix1")
keloidm2.subset<-subset(s.k.total, sample=="keloid_mix2")
keloidm3.subset<-subset(s.k.total, sample=="keloid_mix3")
keloidc1.subset<-subset(s.k.total, sample=="keloid_centr1")
keloidc2.subset<-subset(s.k.total, sample=="keloid_centr2")
keloidc3.subset<-subset(s.k.total, sample=="keloid_centr3")
keloidc4.subset<-subset(s.k.total, sample=="keloid_centr4")
keloidc5.subset<-subset(s.k.total, sample=="keloid_centr5")
keloidc6.subset<-subset(s.k.total, sample=="keloid_centr6")
keloidc7.subset<-subset(s.k.total, sample=="keloid_centr7")
keloidc8.subset<-subset(s.k.total, sample=="keloid_centr8")

#Data frame frequencies
skins1.cluster.frequencies<-as.data.frame(table(skins1.subset$celltype))
skins2.cluster.frequencies<-as.data.frame(table(skins2.subset$celltype))
skins3.cluster.frequencies<-as.data.frame(table(skins3.subset$celltype))
skins4.cluster.frequencies<-as.data.frame(table(skins4.subset$celltype))
skins5.cluster.frequencies<-as.data.frame(table(skins5.subset$celltype))
skins6.cluster.frequencies<-as.data.frame(table(skins6.subset$celltype))
skina1.cluster.frequencies<-as.data.frame(table(skina1.subset$celltype))
skina2.cluster.frequencies<-as.data.frame(table(skina2.subset$celltype))
skina3.cluster.frequencies<-as.data.frame(table(skina3.subset$celltype))
skina4.cluster.frequencies<-as.data.frame(table(skina4.subset$celltype))
nscar1.cluster.frequencies<-as.data.frame(table(nscar1.subset$celltype))
nscar2.cluster.frequencies<-as.data.frame(table(nscar2.subset$celltype))
nscar3.cluster.frequencies<-as.data.frame(table(nscar3.subset$celltype))
nscar4.cluster.frequencies<-as.data.frame(table(nscar4.subset$celltype))
nscar5.cluster.frequencies<-as.data.frame(table(nscar5.subset$celltype))
nscar6.cluster.frequencies<-as.data.frame(table(nscar6.subset$celltype))
keloidm1.cluster.frequencies<-as.data.frame(table(keloidm1.subset$celltype))
keloidm2.cluster.frequencies<-as.data.frame(table(keloidm2.subset$celltype))
keloidm3.cluster.frequencies<-as.data.frame(table(keloidm3.subset$celltype))
keloidc1.cluster.frequencies<-as.data.frame(table(keloidc1.subset$celltype))
keloidc2.cluster.frequencies<-as.data.frame(table(keloidc2.subset$celltype))
keloidc3.cluster.frequencies<-as.data.frame(table(keloidc3.subset$celltype))
keloidc4.cluster.frequencies<-as.data.frame(table(keloidc4.subset$celltype))
keloidc5.cluster.frequencies<-as.data.frame(table(keloidc5.subset$celltype))
keloidc6.cluster.frequencies<-as.data.frame(table(keloidc6.subset$celltype))
keloidc7.cluster.frequencies<-as.data.frame(table(keloidc7.subset$celltype))
keloidc8.cluster.frequencies<-as.data.frame(table(keloidc8.subset$celltype))


#Barplot Cluster percentage
skins1.cluster.frequencies$sample="skin_sep1"
skins2.cluster.frequencies$sample="skin_sep2"
skins3.cluster.frequencies$sample="skin_sep3"
skins4.cluster.frequencies$sample="skin_sep4"
skins5.cluster.frequencies$sample="skin_sep5"
skins6.cluster.frequencies$sample="skin_sep6"
skina1.cluster.frequencies$sample="skin_adj1"
skina2.cluster.frequencies$sample="skin_adj2"
skina3.cluster.frequencies$sample="skin_adj3"
skina4.cluster.frequencies$sample="skin_adj4"
nscar1.cluster.frequencies$sample="nscar_1"
nscar2.cluster.frequencies$sample="nscar_2"
nscar3.cluster.frequencies$sample="nscar_3"
nscar4.cluster.frequencies$sample="nscar_4"
nscar5.cluster.frequencies$sample="nscar_5"
nscar6.cluster.frequencies$sample="nscar_6"
keloidm1.cluster.frequencies$sample="keloid_mix1"
keloidm2.cluster.frequencies$sample="keloid_mix2"
keloidm3.cluster.frequencies$sample="keloid_mix3"
keloidc1.cluster.frequencies$sample="keloid_centr1"
keloidc2.cluster.frequencies$sample="keloid_centr2"
keloidc3.cluster.frequencies$sample="keloid_centr3"
keloidc4.cluster.frequencies$sample="keloid_centr4"
keloidc5.cluster.frequencies$sample="keloid_centr5"
keloidc6.cluster.frequencies$sample="keloid_centr6"
keloidc7.cluster.frequencies$sample="keloid_centr7"
keloidc8.cluster.frequencies$sample="keloid_centr8"


#Barplot Cluster percentage
skins1.cluster.frequencies$area="SKIN"
skins2.cluster.frequencies$area="SKIN"
skins3.cluster.frequencies$area="SKIN"
skins4.cluster.frequencies$area="SKIN"
skins5.cluster.frequencies$area="SKIN"
skins6.cluster.frequencies$area="SKIN"
skina1.cluster.frequencies$area="adjacent SKIN"
skina2.cluster.frequencies$area="adjacent SKIN"
skina3.cluster.frequencies$area="adjacent SKIN"
skina4.cluster.frequencies$area="adjacent SKIN"
nscar1.cluster.frequencies$area="NORMAL SCAR"
nscar2.cluster.frequencies$area="NORMAL SCAR"
nscar3.cluster.frequencies$area="NORMAL SCAR"
nscar4.cluster.frequencies$area="NORMAL SCAR"
nscar5.cluster.frequencies$area="NORMAL SCAR"
nscar6.cluster.frequencies$area="NORMAL SCAR"
keloidm1.cluster.frequencies$area="KELOID"
keloidm2.cluster.frequencies$area="KELOID"
keloidm3.cluster.frequencies$area="KELOID"
keloidc1.cluster.frequencies$area="KELOID"
keloidc2.cluster.frequencies$area="KELOID"
keloidc3.cluster.frequencies$area="KELOID"
keloidc4.cluster.frequencies$area="KELOID"
keloidc5.cluster.frequencies$area="KELOID"
keloidc6.cluster.frequencies$area="KELOID"
keloidc7.cluster.frequencies$area="KELOID"
keloidc8.cluster.frequencies$area="KELOID"

skins1.cluster.frequencies<-mutate(skins1.cluster.frequencies,Percentage= skins1.cluster.frequencies$Freq/sum(skins1.cluster.frequencies$Freq)*100)
skins2.cluster.frequencies<-mutate(skins2.cluster.frequencies,Percentage= skins2.cluster.frequencies$Freq/sum(skins2.cluster.frequencies$Freq)*100)
skins3.cluster.frequencies<-mutate(skins3.cluster.frequencies,Percentage= skins3.cluster.frequencies$Freq/sum(skins3.cluster.frequencies$Freq)*100)
skins4.cluster.frequencies<-mutate(skins4.cluster.frequencies,Percentage= skins4.cluster.frequencies$Freq/sum(skins4.cluster.frequencies$Freq)*100)
skins5.cluster.frequencies<-mutate(skins5.cluster.frequencies,Percentage= skins5.cluster.frequencies$Freq/sum(skins5.cluster.frequencies$Freq)*100)
skins6.cluster.frequencies<-mutate(skins6.cluster.frequencies,Percentage= skins6.cluster.frequencies$Freq/sum(skins6.cluster.frequencies$Freq)*100)
skina1.cluster.frequencies<-mutate(skina1.cluster.frequencies,Percentage= skina1.cluster.frequencies$Freq/sum(skina1.cluster.frequencies$Freq)*100)
skina2.cluster.frequencies<-mutate(skina2.cluster.frequencies,Percentage= skina2.cluster.frequencies$Freq/sum(skina2.cluster.frequencies$Freq)*100)
skina3.cluster.frequencies<-mutate(skina3.cluster.frequencies,Percentage= skina3.cluster.frequencies$Freq/sum(skina3.cluster.frequencies$Freq)*100)
skina4.cluster.frequencies<-mutate(skina4.cluster.frequencies,Percentage= skina4.cluster.frequencies$Freq/sum(skina4.cluster.frequencies$Freq)*100)
nscar1.cluster.frequencies<-mutate(nscar1.cluster.frequencies,Percentage= nscar1.cluster.frequencies$Freq/sum(nscar1.cluster.frequencies$Freq)*100)
nscar2.cluster.frequencies<-mutate(nscar2.cluster.frequencies,Percentage= nscar2.cluster.frequencies$Freq/sum(nscar2.cluster.frequencies$Freq)*100)
nscar3.cluster.frequencies<-mutate(nscar3.cluster.frequencies,Percentage= nscar3.cluster.frequencies$Freq/sum(nscar3.cluster.frequencies$Freq)*100)
nscar4.cluster.frequencies<-mutate(nscar4.cluster.frequencies,Percentage= nscar4.cluster.frequencies$Freq/sum(nscar4.cluster.frequencies$Freq)*100)
nscar5.cluster.frequencies<-mutate(nscar5.cluster.frequencies,Percentage= nscar5.cluster.frequencies$Freq/sum(nscar5.cluster.frequencies$Freq)*100)
nscar6.cluster.frequencies<-mutate(nscar6.cluster.frequencies,Percentage= nscar6.cluster.frequencies$Freq/sum(nscar6.cluster.frequencies$Freq)*100)
keloidm1.cluster.frequencies<-mutate(keloidm1.cluster.frequencies,Percentage= keloidm1.cluster.frequencies$Freq/sum(keloidm1.cluster.frequencies$Freq)*100)
keloidm2.cluster.frequencies<-mutate(keloidm2.cluster.frequencies,Percentage= keloidm2.cluster.frequencies$Freq/sum(keloidm2.cluster.frequencies$Freq)*100)
keloidm3.cluster.frequencies<-mutate(keloidm3.cluster.frequencies,Percentage= keloidm3.cluster.frequencies$Freq/sum(keloidm3.cluster.frequencies$Freq)*100)
keloidc1.cluster.frequencies<-mutate(keloidc1.cluster.frequencies,Percentage= keloidc1.cluster.frequencies$Freq/sum(keloidc1.cluster.frequencies$Freq)*100)
keloidc2.cluster.frequencies<-mutate(keloidc2.cluster.frequencies,Percentage= keloidc2.cluster.frequencies$Freq/sum(keloidc2.cluster.frequencies$Freq)*100)
keloidc3.cluster.frequencies<-mutate(keloidc3.cluster.frequencies,Percentage= keloidc3.cluster.frequencies$Freq/sum(keloidc3.cluster.frequencies$Freq)*100)
keloidc4.cluster.frequencies<-mutate(keloidc4.cluster.frequencies,Percentage= keloidc4.cluster.frequencies$Freq/sum(keloidc4.cluster.frequencies$Freq)*100)
keloidc5.cluster.frequencies<-mutate(keloidc5.cluster.frequencies,Percentage= keloidc5.cluster.frequencies$Freq/sum(keloidc5.cluster.frequencies$Freq)*100)
keloidc6.cluster.frequencies<-mutate(keloidc6.cluster.frequencies,Percentage= keloidc6.cluster.frequencies$Freq/sum(keloidc6.cluster.frequencies$Freq)*100)
keloidc7.cluster.frequencies<-mutate(keloidc7.cluster.frequencies,Percentage= keloidc7.cluster.frequencies$Freq/sum(keloidc7.cluster.frequencies$Freq)*100)
keloidc8.cluster.frequencies<-mutate(keloidc8.cluster.frequencies,Percentage= keloidc8.cluster.frequencies$Freq/sum(keloidc8.cluster.frequencies$Freq)*100)


skin_kel_samp.cluster.df <- rbind(skins1.cluster.frequencies,skins2.cluster.frequencies,skins3.cluster.frequencies,skins4.cluster.frequencies,skins5.cluster.frequencies,skins6.cluster.frequencies,skina1.cluster.frequencies,skina2.cluster.frequencies,skina3.cluster.frequencies,skina4.cluster.frequencies,nscar1.cluster.frequencies,nscar2.cluster.frequencies,nscar3.cluster.frequencies,nscar4.cluster.frequencies,nscar5.cluster.frequencies,nscar6.cluster.frequencies,keloidm1.cluster.frequencies,keloidm2.cluster.frequencies,keloidm3.cluster.frequencies,keloidc1.cluster.frequencies,keloidc2.cluster.frequencies,keloidc3.cluster.frequencies,keloidc4.cluster.frequencies,keloidc5.cluster.frequencies,keloidc6.cluster.frequencies,keloidc7.cluster.frequencies,keloidc8.cluster.frequencies)
skin_kel_samp.cluster.df$sample <- factor(x = skin_kel_samp.cluster.df$sample, levels = c("skin_sep1","skin_sep2","skin_sep3","skin_sep4","skin_sep5","skin_sep6","skin_adj1","skin_adj2","skin_adj3","skin_adj4","nscar_1","nscar_2","nscar_3","nscar_4","nscar_5","nscar_6","keloid_mix1","keloid_mix2", "keloid_mix3", "keloid_centr1", "keloid_centr2", "keloid_centr3", "keloid_centr4", "keloid_centr5", "keloid_centr6", "keloid_centr7", "keloid_centr8"))
write.xlsx(skin_kel_samp.cluster.df, "../celltype_distribution_total.xlsx")




####Subset Schwann Cells####
{
  #subset SC/MEL
  Idents(s.k.total)<- s.k.total$celltype
  SC_11s6n11k.x <- subset(s.k.total,idents = "SC/MEL")
  DefaultAssay(SC_11s6n11k.x)<-"RNA"
  SC_11s6n11k.x[["percent.mt"]] <- PercentageFeatureSet(SC_11s6n11k.x, pattern = "^MT-")
  SC_11s6n11k.x <- SCTransform(SC_11s6n11k.x, method = "glmGamPoi",vars.to.regress = "percent.mt",verbose = F)
  SC_11s6n11k.x <- RunPCA(SC_11s6n11k.x, npcs = 50)
  ElbowPlot(SC_11s6n11k.x, ndims = 50)
  SC_11s6n11k.y<-SC_11s6n11k.x
  SC_11s6n11k.y <- RunUMAP(SC_11s6n11k.y, dims = 1:35)
  SC_11s6n11k.y <- FindNeighbors(SC_11s6n11k.y, dims = 1:35)
  SC_11s6n11k.y <- FindClusters(SC_11s6n11k.y, resolution = 0.2)
  UMAPPlot(SC_11s6n11k.y, label=T)
  DefaultAssay(SC_11s6n11k.y)<- "RNA"
  SC_11s6n11k.y <- NormalizeData(SC_11s6n11k.y)
  #Identification of melanocyte subcluster
  DotPlot(SC_11s6n11k.y,features=c("S100B","PMEL","MLANA","DCN","NGFR","CD1A"))
  Idents(SC_11s6n11k.y)<- SC_11s6n11k.y$seurat_clusters
  #subset of Schwann cell cluster without mel
  SC_11s6n11k.x <- subset(SC_11s6n11k.y,idents = c("0","1","3","4","5","6","7"))
  DefaultAssay(SC_11s6n11k.x)<-"RNA"
  SC_11s6n11k.x[["percent.mt"]] <- PercentageFeatureSet(SC_11s6n11k.x, pattern = "^MT-")
  SC_11s6n11k.x <- SCTransform(SC_11s6n11k.x, method = "glmGamPoi",vars.to.regress = "percent.mt",verbose = F)
  SC_11s6n11k.x <- RunPCA(SC_11s6n11k.x, npcs = 50)
  ElbowPlot(SC_11s6n11k.x, ndims = 50)
  SC_11s6n11k.y<-SC_11s6n11k.x
  SC_11s6n11k.y <- RunUMAP(SC_11s6n11k.y, dims = 1:25)
  SC_11s6n11k.y <- FindNeighbors(SC_11s6n11k.y, dims = 1:25)
  SC_11s6n11k.y <- FindClusters(SC_11s6n11k.y, resolution = 0.1)
  UMAPPlot(SC_11s6n11k.y, label=T)
  SC_11s6n11k<-SC_11s6n11k.y
  DefaultAssay(SC_11s6n11k)<- "RNA"
  SC_11s6n11k <- NormalizeData(SC_11s6n11k)
}




#SC identification

DotPlot(SC_11s6n11k,features=c("MBP","MPZ","PLP1","PMP22","NCAM1","L1CAM","SCN7A","S100B","NGFR","NES","IGFBP5","CCN3","TOP2A","MKI67","ICAM1","SELE","DCN","LUM"))

#Name Cluster
Idents(SC_11s6n11k)<-SC_11s6n11k$seurat_clusters
SC_11s6n11k<-RenameIdents(SC_11s6n11k,
                          `0`="SC-Keloid", `1`="SC-Nonmyel", `2`="SC-FB", `3`="SC-Myel",`4`="SC-EC",`5`="SC-Prolif")
SC_11s6n11k$celltype<-Idents(SC_11s6n11k)


#Define cluster levels
clusters_ordered<-c("SC-Myel","SC-Nonmyel","SC-Keloid","SC-Prolif",
                    "SC-EC","SC-FB")
SC_11s6n11k$celltype<- factor(SC_11s6n11k$celltype, levels = clusters_ordered)
Idents(SC_11s6n11k)<-SC_11s6n11k$celltype

#Name Cluster
Idents(SC_11s6n11k)<-SC_11s6n11k$location
SC_11s6n11k<-RenameIdents(SC_11s6n11k,
                          `skin_sep`="SKIN", `skin_adj`="adjacent SKIN", `nscar_centr`="NORMAL SCAR",   `keloid_mix`="KELOID", `keloid_centr`="KELOID")
SC_11s6n11k$area<-Idents(SC_11s6n11k)

FeaturePlot(SC_11s6n11k,features = c("S100B","LUM"), blend=T, split.by="area",blend.threshold = 0.5, cols = c("lightgrey","#54990F","#990F26"), pt.size = 1, order=T, min.cutoff = 0, max.cutoff = 4)

#Define cluster levels
clusters_ordered<-c("SKIN","adjacent SKIN","NORMAL SCAR","KELOID")
SC_11s6n11k$area<- factor(SC_11s6n11k$area, levels = clusters_ordered)
Idents(SC_11s6n11k)<-SC_11s6n11k$celltype

#Clustermarker-recheck
DotPlot(SC_11s6n11k,features=c("MBP","MPZ","PLP1","PMP22","NCAM1","L1CAM","SCN7A","S100B","NGFR","NES","IGFBP5","CCN3","TOP2A","MKI67","ICAM1","SELE","DCN","LUM"),assay="RNA")+coord_flip() +scale_color_gradient(low="grey", "high"=single_col)+ ggtitle("SC subset_total characterisation") + theme(axis.text.x = element_text(angle = 90))

#Clustermarker
  SC_11s6n11k.clustermarker<-FindAllMarkers(SC_11s6n11k,  assay = "RNA", min.pct = 0.25, logfc.threshold = 0.25)
  SC_11s6n11k.clustermarker$Foldchange_UP <- 2^(SC_11s6n11k.clustermarker$avg_log2FC)
  SC_11s6n11k.clustermarker$Foldchange_DOWN <- 2^(-SC_11s6n11k.clustermarker$avg_log2FC)
  SC_11s6n11k.clustermarker$Ratio_pct1_pct2 <- (SC_11s6n11k.clustermarker$pct.1)/(SC_11s6n11k.clustermarker$pct.2)
  write.xlsx(SC_11s6n11k.clustermarker, "../Clustermarker_SC_11s6n11k.xlsx")

#compute DEG Keloid vs. myel/nonmyel
  Idents(SC_11s6n11k)<-"celltype"
  SC_11s6n11k.kvsmnm.clustermarker<-FindMarkers(SC_11s6n11k,  ident.1= "SC-Repair", ident.2= c("SC-Myel","SC-Nonmyel"), assay = "RNA", min.pct = 0.01, logfc.threshold = 0.01)
  SC_11s6n11k.kvsmnm.clustermarker$Foldchange_UP <- 2^(SC_11s6n11k.kvsmnm.clustermarker$avg_log2FC)
  SC_11s6n11k.kvsmnm.clustermarker$Foldchange_DOWN <- 2^(-SC_11s6n11k.kvsmnm.clustermarker$avg_log2FC)
  SC_11s6n11k.kvsmnm.clustermarker$Ratio_pct1_pct2 <- (SC_11s6n11k.kvsmnm.clustermarker$pct.1)/(SC_11s6n11k.kvsmnm.clustermarker$pct.2)
  write.xlsx(SC_11s6n11k.kvsmnm.clustermarker, "Lists/DEG_SCrepaire_vs_SCmyel_nonmyel_total.xlsx")

#UMAPPlots, named clusters

  UMAPPlot(SC_11s6n11k, label=F)+ scale_color_manual(values=color_SC)+ ggtitle("SC-UMAP")


#Percentage total
###Pieplot/Barplot cellnumber /sample/cluster
#subset conditions
SC_skin1.subset<-subset(SC_11s6n11k, sample=="skin_sep1")
SC_skin2.subset<-subset(SC_11s6n11k, sample=="skin_sep2")
SC_skin4.subset<-subset(SC_11s6n11k, sample=="skin_sep4")
SC_skin3.subset<-subset(SC_11s6n11k, sample=="skin_sep3")
SC_skin5.subset<-subset(SC_11s6n11k, sample=="skin_sep5")
SC_skin6.subset<-subset(SC_11s6n11k, sample=="skin_sep6")
SC_skina1.subset<-subset(SC_11s6n11k, sample=="skin_adj1")
SC_skina2.subset<-subset(SC_11s6n11k, sample=="skin_adj2")
SC_skina3.subset<-subset(SC_11s6n11k, sample=="skin_adj3")
SC_skina4.subset<-subset(SC_11s6n11k, sample=="skin_adj4")
SC_nscar1.subset<-subset(SC_11s6n11k, sample=="nscar_1")
SC_nscar2.subset<-subset(SC_11s6n11k, sample=="nscar_2")
SC_nscar3.subset<-subset(SC_11s6n11k, sample=="nscar_3")
SC_nscar4.subset<-subset(SC_11s6n11k, sample=="nscar_4")
SC_nscar5.subset<-subset(SC_11s6n11k, sample=="nscar_5")
SC_nscar6.subset<-subset(SC_11s6n11k, sample=="nscar_6")
SC_keloidm1.subset<-subset(SC_11s6n11k, sample=="keloid_mix1")
SC_keloidm2.subset<-subset(SC_11s6n11k, sample=="keloid_mix2")
SC_keloidm3.subset<-subset(SC_11s6n11k, sample=="keloid_mix3")
SC_keloidc1.subset<-subset(SC_11s6n11k, sample=="keloid_centr1")
SC_keloidc2.subset<-subset(SC_11s6n11k, sample=="keloid_centr2")
SC_keloidc3.subset<-subset(SC_11s6n11k, sample=="keloid_centr3")
SC_keloidc4.subset<-subset(SC_11s6n11k, sample=="keloid_centr4")
SC_keloidc5.subset<-subset(SC_11s6n11k, sample=="keloid_centr5")
SC_keloidc6.subset<-subset(SC_11s6n11k, sample=="keloid_centr6")
SC_keloidc7.subset<-subset(SC_11s6n11k, sample=="keloid_centr7")
SC_keloidc8.subset<-subset(SC_11s6n11k, sample=="keloid_centr8")


#Data frame frequencies
SC_skins1.cluster.frequencies<-as.data.frame(table(SC_skin1.subset$celltype))
SC_skins2.cluster.frequencies<-as.data.frame(table(SC_skin2.subset$celltype))
SC_skins3.cluster.frequencies<-as.data.frame(table(SC_skin3.subset$celltype))
SC_skins4.cluster.frequencies<-as.data.frame(table(SC_skin4.subset$celltype))
SC_skins5.cluster.frequencies<-as.data.frame(table(SC_skin5.subset$celltype))
SC_skins6.cluster.frequencies<-as.data.frame(table(SC_skin6.subset$celltype))
SC_skina1.cluster.frequencies<-as.data.frame(table(SC_skina1.subset$celltype))
SC_skina2.cluster.frequencies<-as.data.frame(table(SC_skina2.subset$celltype))
SC_skina3.cluster.frequencies<-as.data.frame(table(SC_skina3.subset$celltype))
SC_skina4.cluster.frequencies<-as.data.frame(table(SC_skina4.subset$celltype))
SC_nscar1.cluster.frequencies<-as.data.frame(table(SC_nscar1.subset$celltype))
SC_nscar2.cluster.frequencies<-as.data.frame(table(SC_nscar2.subset$celltype))
SC_nscar3.cluster.frequencies<-as.data.frame(table(SC_nscar3.subset$celltype))
SC_nscar4.cluster.frequencies<-as.data.frame(table(SC_nscar4.subset$celltype))
SC_nscar5.cluster.frequencies<-as.data.frame(table(SC_nscar5.subset$celltype))
SC_nscar6.cluster.frequencies<-as.data.frame(table(SC_nscar6.subset$celltype))
SC_keloidm1.cluster.frequencies<-as.data.frame(table(SC_keloidm1.subset$celltype))
SC_keloidm2.cluster.frequencies<-as.data.frame(table(SC_keloidm2.subset$celltype))
SC_keloidm3.cluster.frequencies<-as.data.frame(table(SC_keloidm3.subset$celltype))
SC_keloidc1.cluster.frequencies<-as.data.frame(table(SC_keloidc1.subset$celltype))
SC_keloidc2.cluster.frequencies<-as.data.frame(table(SC_keloidc2.subset$celltype))
SC_keloidc3.cluster.frequencies<-as.data.frame(table(SC_keloidc3.subset$celltype))
SC_keloidc4.cluster.frequencies<-as.data.frame(table(SC_keloidc4.subset$celltype))
SC_keloidc5.cluster.frequencies<-as.data.frame(table(SC_keloidc5.subset$celltype))
SC_keloidc6.cluster.frequencies<-as.data.frame(table(SC_keloidc6.subset$celltype))
SC_keloidc7.cluster.frequencies<-as.data.frame(table(SC_keloidc7.subset$celltype))
SC_keloidc8.cluster.frequencies<-as.data.frame(table(SC_keloidc8.subset$celltype))

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


SC_skin_kel_samp.cluster.df <- rbind(SC_skins1.cluster.frequencies,SC_skins2.cluster.frequencies,SC_skins3.cluster.frequencies,SC_skins4.cluster.frequencies,SC_skins5.cluster.frequencies,SC_skins6.cluster.frequencies,SC_skina1.cluster.frequencies,SC_skina2.cluster.frequencies,SC_skina3.cluster.frequencies,SC_skina4.cluster.frequencies,SC_nscar1.cluster.frequencies,SC_nscar2.cluster.frequencies,SC_nscar3.cluster.frequencies,SC_nscar4.cluster.frequencies,SC_nscar5.cluster.frequencies,SC_nscar6.cluster.frequencies,SC_keloidm1.cluster.frequencies,SC_keloidm2.cluster.frequencies,SC_keloidm3.cluster.frequencies,SC_keloidc1.cluster.frequencies,SC_keloidc2.cluster.frequencies,SC_keloidc3.cluster.frequencies,SC_keloidc4.cluster.frequencies,SC_keloidc5.cluster.frequencies,SC_keloidc6.cluster.frequencies,SC_keloidc7.cluster.frequencies,SC_keloidc8.cluster.frequencies)
SC_skin_kel_samp.cluster.df$sample <- factor(x = SC_skin_kel_samp.cluster.df$sample, levels = c("skin_sep1","skin_sep2","skin_sep3", "skin_sep4", "skin_sep5", "skin_sep6", "skin_adj1","skin_adj2","skin_adj3", "skin_adj4", "nscar_1","nscar_2","nscar_3","nscar_4","nscar_5","nscar_6","keloid_mix1", "keloid_mix2", "keloid_mix3","keloid_centr1","keloid_centr2","keloid_centr3","keloid_centr4","keloid_centr5","keloid_centr6","keloid_centr7","keloid_centr8"))
write.xlsx(SC_skin_kel_samp.cluster.df, "../celltype_distribution_SC.xlsx")

ggplot(data=SC_skin_kel_samp.cluster.df, aes(x=sample, y=Freq, fill=Var1)) +
  geom_bar(stat="identity", position="stack") +  theme_classic() +
  ggtitle("Cell source_schwann cell cluster") + xlab("condition")+ylab("total cell count") +scale_fill_manual(values =color_SC)
