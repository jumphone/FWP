
TAGS=c('1142243F','1160920F','CID4290','CID4465','CID44971','CID4535')

library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)


TAG='CID4465'
print(TAG)

this.mat=Read10X(paste0('data/spatial/raw_count_matrices/',TAG,'_raw_feature_bc_matrix/'),gene.column=2)
this.image=Read10X_Image(paste0('data/spatial/spatial/',TAG,'_spatial/'))
DefaultAssay(object = this.image) <- 'Spatial'
##########
this.obj=CreateSeuratObject(counts =this.mat, project = 'this', assay = 'Spatial',min.cells=1, min.features=0)
this.obj[['image']] <- this.image
########
this.meta=read.csv(paste0('data/spatial/metadata/',TAG,'_metadata.csv'),row.names=1)
this.obj=AddMetaData(this.obj,this.meta)
this.obj=subset(this.obj, cells=colnames(this.obj)[which(!is.na(this.obj$Classification))])
table(this.obj$Classification)
this.obj$tag=rep(0,ncol(this.obj))
this.obj$tag[which(stringr::str_detect(this.obj$Classification,'cancer'))]=1



this.obj <- NormalizeData(object = this.obj, normalization.method = "LogNormalize", scale.factor = 10000)

#SpatialDimPlot(this.obj,group.by='Classification', label = TRUE, label.size = 3)
this.obj$rtag=1-this.obj$tag

pdf('plot/p14_spatial_type.pdf',width=4,height=3.5)
SpatialDimPlot(this.obj,group.by='rtag', label = FALSE ,pt.size.factor=2)
dev.off()

this.data=as.matrix(this.obj[['Spatial']]@data)

source('fwp.R')
FW=readRDS('feature_weight/FW_bulk_TCGA_CESC.rds')
out.fwp=fwp(this.data,FW)
.evaluate(out.fwp,this.obj$tag)



this.obj$fwp.CESC=out.fwp
pdf('plot/p15_spatial_fwp_CESC.pdf',width=4,height=3.5)
SpatialFeaturePlot(this.obj, features = c("fwp.CESC"),alpha=c(1,1),pt.size.factor=2)
dev.off()



pdf('plot/p16_spatial_typeOrigin.pdf',width=4,height=3.5)
SpatialDimPlot(this.obj,group.by='Classification', label = FALSE ,pt.size.factor=2)
dev.off()


source('fwp.R')
FW=readRDS('feature_weight/FW_bulk_TCGA_BRCA.rds')
out.fwp=fwp(this.data,FW)
.evaluate(out.fwp,this.obj$tag)


this.obj$fwp.BRCA=out.fwp
pdf('plot/p17_spatial_fwp_BRCA.pdf',width=4,height=3.5)
SpatialFeaturePlot(this.obj, features = c("fwp.BRCA"),alpha=c(1,1),pt.size.factor=2)
dev.off()










this.obj$fwp=out.fwp

SpatialFeaturePlot(this.obj, features = c("fwp"),alpha=c(1,1),pt.size.factor=3)
#SpatialFeaturePlot(this.obj, features = c("fwp"),alpha=c(0.1,0.5),pt.size.factor=3)

plot(this.obj$fwp,this.obj$nFeature_RNA)

SpatialFeaturePlot(this.obj, features = c("nFeature_RNA"),alpha=c(1,1),pt.size.factor=3)


plot(this.obj$fwp,this.obj$nFeature_RNA)

.evaluate(out.fwp,this.obj$tag)





out.fwo=fwo(this.data,FW)
.evaluate(out.fwo,this.obj$tag)


this.obj$fwp=out.fwp








