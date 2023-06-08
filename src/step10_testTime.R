

library(Seurat)

pbmc=readRDS( '/home/database/data/fitTest/data/Wu_etal_2021_BRCA_scRNASeq/pbmc5k.rds')
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

CT=pbmc$celltype_major
TAG=rep(0,length(CT))
TAG[which(CT=='Cancer Epithelial')]=1
names(TAG)=names(CT)
########################

pbmc$tag=TAG

source('fwp.R')
source('fwp_source.R')
DIR='data/first/'
TRAIN.TYPE='BRCA'
TRAIN.DATA=.loadFileNoGap(paste0(DIR,TRAIN.TYPE,'/mat_train.tsv'))
TRAIN.TAG=read.csv(paste0(DIR,TRAIN.TYPE,'/CorrectDP_train.csv'),header=F)
BRCA.FW=.calFW(TRAIN.DATA,TRAIN.TAG[,2])


DATA=as.matrix(pbmc[['RNA']]@data)
out.fwp=fwp(DATA, BRCA.FW,npcs=2)
.evaluate(out.fwp,pbmc$tag)
#Number: 5000
#AUC: 0.958
#Time: 24s

this_pbmc=subset(pbmc, cells=sample(colnames(pbmc),4000))
DATA=as.matrix(this_pbmc[['RNA']]@data)
out.fwp=fwp(DATA, BRCA.FW,npcs=2)
.evaluate(out.fwp,this_pbmc$tag)
#Number: 4000
#AUC: 0.956
#Time: 17s

this_pbmc=subset(pbmc, cells=sample(colnames(pbmc),3000))
DATA=as.matrix(this_pbmc[['RNA']]@data)
out.fwp=fwp(DATA, BRCA.FW,npcs=2)
.evaluate(out.fwp,this_pbmc$tag)
#Number: 3000
#AUC: 0.957
#Time: 12s

this_pbmc=subset(pbmc, cells=sample(colnames(pbmc),2000))
DATA=as.matrix(this_pbmc[['RNA']]@data)
out.fwp=fwp(DATA, BRCA.FW,npcs=2)
.evaluate(out.fwp,this_pbmc$tag)
#Number: 2000
#AUC: 0.955
#Time: 7s

this_pbmc=subset(pbmc, cells=sample(colnames(pbmc),1000))
DATA=as.matrix(this_pbmc[['RNA']]@data)
out.fwp=fwp(DATA, BRCA.FW,npcs=2)
.evaluate(out.fwp,this_pbmc$tag)
#Number: 1000
#AUC: 0.956
#Time: 3s

this_pbmc=subset(pbmc, cells=sample(colnames(pbmc),100))
DATA=as.matrix(this_pbmc[['RNA']]@data)
out.fwp=fwp(DATA, BRCA.FW,npcs=2)
.evaluate(out.fwp,this_pbmc$tag)
#Number: 100
#AUC: 0.912
#Time: 1s





AUC=c(0.958,0.956,0.957,0.955,0.956,0.912)
TIME=c(24,17,12,7,3,1)
N=c(5000,4000,3000,2000,1000,100)

AUC=AUC[order(N)]
TIME=TIME[order(N)]
N=N[order(N)]

pdf('plot/p11_testTime.pdf',width=5,height=5)
barplot(AUC,col='indianred2')
abline(h=0.5,lty=2)

plot(N,TIME,pch='+',cex=2,col='indianred2',type='both',ylim=c(0,25))
dev.off()









