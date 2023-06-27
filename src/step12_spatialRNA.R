
TAGS=c('1142243F','1160920F','CID4290','CID4465','CID44971','CID4535')

library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)

i=1
while(i<=length(TAGS)){

TAG=TAGS[i]
print(TAG)

this.mat=Read10X(paste0('data/spatial/raw_count_matrices/',TAG,'_raw_feature_bc_matrix/'),gene.column=2)
this.image=Read10X_Image(paste0('data/spatial/spatial/',TAG,'_spatial/'))
DefaultAssay(object = this.image) <- 'Spatial'
##########
this.obj=CreateSeuratObject(counts =this.mat, project = 'this', assay = 'Spatial',min.cells=1, min.features=0)
this.obj[['image']] <- this.image
########
this.obj <- NormalizeData(object = this.obj, normalization.method = "LogNormalize", scale.factor = 10000)
###########
this.meta=read.csv(paste0('data/spatial/metadata/',TAG,'_metadata.csv'),row.names=1)
this.obj=AddMetaData(this.obj,this.meta)
this.obj=subset(this.obj, cells=colnames(this.obj)[which(!is.na(this.obj$Classification))])
#SpatialDimPlot(this.obj,group.by='Classification', label = TRUE, label.size = 3)
#SpatialFeaturePlot(this.obj, features = c("CD44"))
#boxplot(this.obj$nFeature_Spatial~this.obj$Classification)
table(this.obj$Classification)
this.obj$tag=rep(0,ncol(this.obj))
this.obj$tag[which(stringr::str_detect(this.obj$Classification,'cancer'))]=1
this.data=as.matrix(this.obj[['Spatial']]@data)

ddd=apply(this.data,2,round,3)
ttt=this.obj$tag

write.table(ttt,file=paste0('./data/spatial_processed/',TAG,'/CorrectDP_train.csv'),sep=',',row.names=T,col.names=F,quote=F)
write.table(ddd,paste0('./data/spatial_processed/',TAG,'/mat_train.tsv'),sep='\t',row.names=T,col.names=T,quote=F)

i=i+1
}








source('fwp.R')
FW=readRDS('feature_weight/FW_bulk_TCGA_BRCA.rds')
out.fwp=fwp(this.data,FW)
out.fwo=fwo(this.data,FW)

.evaluate(out.fwp,this.obj$tag)
.evaluate(out.fwo,this.obj$tag)


this.obj$fwp=out.fwp








