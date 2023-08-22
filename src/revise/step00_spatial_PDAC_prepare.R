
library(Seurat)

##############
pbmc=readRDS('PDAC_DonorA_FFPE-probes.rds')
SpatialDimPlot(pbmc,group.by='integrated_snn_res.0.2')


###################
MAT=Read10X_h5('GSM6505133_DonorA_FFPE-probes_filtered_feature_bc_matrix.h5')
PN=colnames(pbmc)
MN=colnames(MAT)

################
TMP=t(matrix(unlist(stringr::str_split(PN,'_')),nrow=6))[,6]
CI=match(TMP,MN)
MAT=MAT[,CI]

####################
ppp= CreateSeuratObject(counts = MAT, project = "pbmc3k", min.cells = 0, min.features = 0)
ppp <- NormalizeData(ppp, normalization.method = "LogNormalize", scale.factor = 10000)
MAT=ppp[['RNA']]@data
MAT=as.matrix(MAT)
MAT=apply(MAT,2,round,3)
############################

CLST=as.character(pbmc@meta.data$integrated_snn_res.0.2)
TYPE=rep('0',length(CLST))
TYPE[which(CLST %in% c('4'))]='1'

pbmc$type=TYPE
SpatialDimPlot(pbmc,group.by='type')

############################################
OTYPE=cbind(colnames(MAT),TYPE)
write.table(MAT, 'PDAC_DA_FFPE_PR/mat_test.tsv',sep='\t',row.names=T,col.names=T,quote=F)
write.table(OTYPE, 'PDAC_DA_FFPE_PR/CorrectDP_test.csv',sep=',',row.names=F,col.names=F,quote=F)









############################
#PDAC_DonorB_FFPE-probes.rds

pbmc=readRDS('PDAC_DonorB_FFPE-probes.rds')
SpatialDimPlot(pbmc,group.by='SCT_snn_res.0.3')


#################
MAT=Read10X_h5('GSM6505134_DonorB_FFPE-probes_filtered_feature_bc_matrix.h5')
PN=colnames(pbmc)
MN=colnames(MAT)

################
#TMP=t(matrix(unlist(stringr::str_split(PN,'_')),nrow=6))[,6]
TMP=PN
CI=match(TMP,MN)
MAT=MAT[,CI]

####################
ppp= CreateSeuratObject(counts = MAT, project = "pbmc3k", min.cells = 0, min.features = 0)
ppp <- NormalizeData(ppp, normalization.method = "LogNormalize", scale.factor = 10000)
MAT=ppp[['RNA']]@data
MAT=as.matrix(MAT)
MAT=apply(MAT,2,round,3)
############################

CLST=as.character(pbmc@meta.data$SCT_snn_res.0.3)
TYPE=rep('0',length(CLST))
TYPE[which(CLST %in% c('2'))]='1'

pbmc$type=TYPE
SpatialDimPlot(pbmc,group.by='type')

OTYPE=cbind(colnames(MAT),TYPE)
write.table(MAT, 'PDAC_DB_FFPE_PR/mat_test.tsv',sep='\t',row.names=T,col.names=T,quote=F)
write.table(OTYPE, 'PDAC_DB_FFPE_PR/CorrectDP_test.csv',sep=',',row.names=F,col.names=F,quote=F)



#########################
#PDAC_DonorC_FFPE-probes.rds

pbmc=readRDS('PDAC_DonorC_FFPE-probes.rds')
SpatialDimPlot(pbmc,group.by='SCT_snn_res.0.6')


#################
MAT=Read10X_h5('GSM6505135_DonorC_FFPE-probes_filtered_feature_bc_matrix.h5')
PN=colnames(pbmc)
MN=colnames(MAT)

################
#TMP=t(matrix(unlist(stringr::str_split(PN,'_')),nrow=6))[,6]
TMP=PN
CI=match(TMP,MN)
MAT=MAT[,CI]

####################
ppp= CreateSeuratObject(counts = MAT, project = "pbmc3k", min.cells = 0, min.features = 0)
ppp <- NormalizeData(ppp, normalization.method = "LogNormalize", scale.factor = 10000)
MAT=ppp[['RNA']]@data
MAT=as.matrix(MAT)
MAT=apply(MAT,2,round,3)
############################


CLST=as.character(pbmc@meta.data$SCT_snn_res.0.6)
TYPE=rep('0',length(CLST))
TYPE[which(CLST %in% c('1'))]='1'

pbmc$type=TYPE
SpatialDimPlot(pbmc,group.by='type')

OTYPE=cbind(colnames(MAT),TYPE)
write.table(MAT, 'PDAC_DC_FFPE_PR/mat_test.tsv',sep='\t',row.names=T,col.names=T,quote=F)
write.table(OTYPE, 'PDAC_DC_FFPE_PR/CorrectDP_test.csv',sep=',',row.names=F,col.names=F,quote=F)
















