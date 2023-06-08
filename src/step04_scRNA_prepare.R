

library(Seurat)


DATA=Read10X('/home/database/data/fitCancer/data/Wu_etal_2021_BRCA_scRNASeq/',gene.column =1)
META=read.table('/home/database/data/fitCancer/data/Wu_etal_2021_BRCA_scRNASeq/metadata.csv',sep=',',header=TRUE,row.names=1)

pbmc=CreateSeuratObject(counts = DATA, min.cells = 0, min.features = 0, project = "ALL",meta.data=META)

saveRDS(pbmc, '/home/database/data/fitCancer/data/Wu_etal_2021_BRCA_scRNASeq/pbmc.rds')

set.seed(123)
pbmc10k=subset(pbmc, cells=sample(colnames(pbmc),10000))
saveRDS(pbmc10k, '/home/database/data/fitCancer/data/Wu_etal_2021_BRCA_scRNASeq/pbmc10k.rds')


library(Seurat)
pbmc10k=readRDS( '/home/database/data/fitTest/data/Wu_etal_2021_BRCA_scRNASeq/pbmc10k.rds')

set.seed(123)
pbmc5k=subset(pbmc10k, cells=sample(colnames(pbmc10k),5000))
saveRDS(pbmc5k, '/home/database/data/fitTest/data/Wu_etal_2021_BRCA_scRNASeq/pbmc5k.rds')



library(Seurat)
pbmc=readRDS( '/home/database/data/fitTest/data/Wu_etal_2021_BRCA_scRNASeq/pbmc5k.rds')
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
DATA=as.matrix(pbmc[['RNA']]@data)
DATA=apply(DATA,2,round,3)
CT=pbmc$celltype_major
TAG=rep(0,length(CT))
TAG[which(CT=='Cancer Epithelial')]=1
names(TAG)=names(CT)
write.table(TAG,file='./data/scRNA/scBRCA/CorrectDP_test.csv',
            sep=',',row.names=T,col.names=F,quote=F)

write.table(DATA,'./data/scRNA/scBRCA/mat_test.tsv',sep='\t',row.names=T,col.names=T,quote=F)


########################



library(Seurat)
pbmc=readRDS( '/home/database/data/fitTest/data/Wu_etal_2021_BRCA_scRNASeq/pbmc5k.rds')

source('BEER.R')
DATA=pbmc[['RNA']]@counts
BATCH=pbmc$orig.ident

mybeer=BEER(DATA, BATCH, GNUM=30, PCNUM=50, ROUND=1, GN=2000, SEED=1, COMBAT=TRUE )
saveRDS(mybeer,'/home/database/data/fitTest/data/Wu_etal_2021_BRCA_scRNASeq/mybeer_pbmc5k.rds')

# Check selected PCs
PCUSE=mybeer$select
COL=rep('black',length(mybeer$cor))
COL[PCUSE]='red'
plot(mybeer$cor,mybeer$lcor,pch=16,col=COL,
    xlab='Rank Correlation',ylab='Linear Correlation',xlim=c(0,1),ylim=c(0,1))

pb <- mybeer$seurat
PCUSE <- mybeer$select
pb <- RunUMAP(object = pb, reduction='pca',dims = PCUSE, check_duplicates=FALSE,min.dist=0.5)
DimPlot(pb, reduction='umap', group.by='batch', pt.size=0.1) +NoLegend()
pb$celltype_major=pbmc$celltype_major



pb$type=pb$celltype_major
pb$type[which(pb$celltype_major=='Endothelial')]='End'
pb$type[which(pb$celltype_major=='Normal Epithelial')]='NEP'
pb$type[which(pb$celltype_major=='Plasmablasts')]='PLB'
pb$type[which(pb$celltype_major=='Cancer Epithelial')]='Cancer'


saveRDS(pb,'/home/database/data/fitTest/data/Wu_etal_2021_BRCA_scRNASeq/pb_mybeer_pbmc5k.rds')

pdf('plot/p08_umap_label.pdf',width=4,height=4)
DimPlot(pb, reduction='umap', group.by='type', pt.size=0.1, label=TRUE) +NoLegend()
DimPlot(pb, reduction='umap', group.by='type', pt.size=0.1, label=FALSE) +NoLegend()
dev.off()

library(Seurat)
source('fwp.R')
source('fwp_source.R')

pb=readRDS('/home/database/data/fitTest/data/Wu_etal_2021_BRCA_scRNASeq/pb_mybeer_pbmc5k.rds')

DIR='data/first/'
TRAIN.TYPE='BRCA'
TRAIN.DATA=.loadFileNoGap(paste0(DIR,TRAIN.TYPE,'/mat_train.tsv'))
TRAIN.TAG=read.csv(paste0(DIR,TRAIN.TYPE,'/CorrectDP_train.csv'),header=F)

BRCA.FW=.calFW(TRAIN.DATA,TRAIN.TAG[,2])

#source('fwp_bac4.R')
source('fwp.R')
source('fwp_source.R')

DATA=as.matrix(pb[['RNA']]@data)
out.fwp=fwp(DATA, BRCA.FW,npcs=2)

pb$fwp=out.fwp
saveRDS(pb,'/home/database/data/fitTest/data/Wu_etal_2021_BRCA_scRNASeq/pb_mybeer_pbmc5k.rds')

source('source.R')

.qcolY<-function(V){
    COL=visa.vcol(V,
              c(quantile(V,0),quantile(V,0.5),quantile(V,0.85),
                quantile(V,0.90),quantile(V,1)),
              c('blue1','royalblue1','grey90','indianred1','red1')
             )
    return(COL)
    }


UMAP=pb@reductions$umap@cell.embeddings

pdf('plot/p09_umap_fwp.pdf',width=3.5,height=4)
plot(UMAP,pch=16,col=.qcolY(pb$fwp),cex=0.5)
V=c(-35:224)/1000
plot(x=V,y=rep(1,length(V)),col=.qcolY(V),type='h',lwd=3,ylim=c(0,1))
dev.off()






##########
library(Seurat)
source('fwp.R')
source('fwp_source.R')

pb=readRDS('/home/database/data/fitTest/data/Wu_etal_2021_BRCA_scRNASeq/pb_mybeer_pbmc5k.rds')

DIR='data/second/'
TRAIN.TYPE='CESC'
TRAIN.DATA=.loadFileNoGap(paste0(DIR,TRAIN.TYPE,'/mat_train.tsv'))
TRAIN.TAG=read.csv(paste0(DIR,TRAIN.TYPE,'/CorrectDP_train.csv'),header=F)

BRCA.FW=.calFW(TRAIN.DATA,TRAIN.TAG[,2])

#source('fwp_bac4.R')
source('fwp.R')
source('fwp_source.R')

DATA=as.matrix(pb[['RNA']]@data)
out.fwp=fwp(DATA, BRCA.FW,npcs=2)

pb$fwp=out.fwp
#saveRDS(pb,'/home/database/data/fitTest/data/Wu_etal_2021_BRCA_scRNASeq/pb_mybeer_pbmc5k.rds')


source('source.R')

.qcolY<-function(V){
    COL=visa.vcol(V,
              c(quantile(V,0),quantile(V,0.5),quantile(V,0.85),
                quantile(V,0.9),quantile(V,1)),
              c('blue1','royalblue1','grey90','indianred1','red1')
             )
    return(COL)
    }

MIN=min(pb$fwp)
MAX=max(pb$fwp)

UMAP=pb@reductions$umap@cell.embeddings

pdf('plot/p09_umap_fwp_CESC.pdf',width=3.5,height=4)
plot(UMAP,pch=16,col=.qcolY(pb$fwp),cex=0.5)
V=c(-(abs(MIN)*1000):(MAX*1000))/1000
plot(x=V,y=rep(1,length(V)),col=.qcolY(V),type='h',lwd=3,ylim=c(0,1))
dev.off()






