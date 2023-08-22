


a=read.table('GSE70630_OG_processed_data_v2.txt',header=T)

library(Seurat)

#a[which(is.na(a))]=0
RS=rowSums(a)
MAT=a[which(RS>0),]
CS=colSums(MAT)
saveRDS(MAT, file='MAT.rds')






library(Seurat)
pbmc=CreateSeuratObject(counts = MAT, project = "pbmc3k", min.cells = 0, min.features = 0)
BATCH=as.character(pbmc$orig.ident)
table(BATCH)

BATCH[which(BATCH=='X93')]='MGH93'
BATCH[which(BATCH=='X97')]='MGH97'


DATA=MAT

source('/home/toolkit/src/BEER.R')

mybeer=BEER(DATA, BATCH, GNUM=30, PCNUM=50, ROUND=1, GN=2000, SEED=1, COMBAT=TRUE )


pbmc <- mybeer$seurat
PCUSE <- mybeer$select
pbmc <- RunUMAP(object = pbmc, reduction='pca',dims = PCUSE, check_duplicates=FALSE)

DimPlot(pbmc, reduction='umap', group.by='batch', pt.size=0.1) 


VEC=pbmc@reductions$umap@cell.embeddings
set.seed(123)

KM=kmeans(VEC,center=10)

TYPE=rep('1',ncol(pbmc))
TYPE[which(KM$cluster %in% c(3,4))]=0

plot(VEC, col=as.factor(TYPE),pch=as.character(TYPE))

MAT=DATA


OTYPE=cbind(colnames(MAT),TYPE)
write.table(MAT, 'OLIG/mat_test.tsv',sep='\t',row.names=T,col.names=T,quote=F)
write.table(OTYPE, 'OLIG/CorrectDP_test.csv',sep=',',row.names=F,col.names=F,quote=F)

saveRDS(pbmc, file='pbmc.rds')






























