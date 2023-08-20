########################
library(pROC)
library(MuSiC)
library(Biobase)
library(Seurat)

source('fwp.R')
source('fwp_source.R')

F_TYPE=c('BRCA','KIRP','LIHC','PCPG')
F_DIR=rep('./data/first/',length(F_TYPE))

S_TYPE=c('BLCA','CESC','CHOL','COAD','ESCA','HNSC','KICH','LUSC','PAAD','THYM','UCEC')
S_DIR=rep('./data/second/',length(S_TYPE))

TYPE=c(F_TYPE, S_TYPE)
DIR=c(F_DIR, S_DIR)

TRAIN_TAB=c()
TEST_TAB=c()

i=1
while(i<=length(TYPE)){
this_dir=DIR[i]
this_type=TYPE[i]
############
print(this_dir)
print(this_type)
############
TRAIN.DATA=.loadFileNoGap(paste0(this_dir,'/',this_type,'/mat_train.tsv'))
TRAIN.TAG=read.csv(paste0(this_dir,'/',this_type,'/CorrectDP_train.csv'),header=F)
TEST.DATA=.loadFileNoGap(paste0(this_dir,'/',this_type,'/mat_test.tsv'))
TEST.TAG=read.csv(paste0(this_dir,'/',this_type,'/CorrectDP_test.csv'),header=F)
this_train_tab=table(TRAIN.TAG[,2])
this_test_tab=table(TEST.TAG[,2])
TRAIN_TAB=rbind(TRAIN_TAB,this_train_tab)
TEST_TAB=rbind(TEST_TAB,this_test_tab)
i=i+1
}

rownames(TRAIN_TAB)=TYPE
rownames(TEST_TAB)=TYPE

barplot(t(TRAIN_TAB))
barplot(t(TEST_TAB))

TAG=c(rep(1,length(F_TYPE)), rep(2,length(S_TYPE)))
SUM1=rowSums(TRAIN_TAB[,c(1:2)])
SUM2=rowSums(TEST_TAB[,c(1:2)])

TRAIN_TAB=cbind(TRAIN_TAB,SUM1,TAG)
TEST_TAB=cbind(TEST_TAB,SUM2,TAG)
#######################################
saveRDS(TRAIN_TAB,file='./data/TRAIN_TAB.rds')
saveRDS(TEST_TAB,file='./data/TEST_TAB.rds')




source('fwp.R')
source('fwp_source.R')

F_TYPE=c('BRCA','KIRP','LIHC','PCPG')
F_DIR=rep('./data/first/',length(F_TYPE))

S_TYPE=c('BLCA','CESC','CHOL','COAD','ESCA','HNSC','KICH','LUSC','PAAD','THYM','UCEC')
S_DIR=rep('./data/second/',length(S_TYPE))

TYPE=c(F_TYPE, S_TYPE)
DIR=c(F_DIR, S_DIR)

TRAIN_TAB=readRDS(file='./data/TRAIN_TAB.rds')
TEST_TAB=readRDS(file='./data/TEST_TAB.rds')


TMP_TRAIN_TAB=TRAIN_TAB[which(TRAIN_TAB[,4]==1),]
TMP_TEST_TAB=TEST_TAB[which(TEST_TAB[,4]==1),]
TMP_ORDER=order(-TRAIN_TAB[which(TRAIN_TAB[,4]==1),3])

TMP_TRAIN_TAB=TMP_TRAIN_TAB[TMP_ORDER,]
TMP_TEST_TAB=TMP_TEST_TAB[TMP_ORDER,]

pdf('plot/p01_barplot_first.pdf',width=3,height=4)
barplot(t(TMP_TRAIN_TAB[,c(1:2)]),col=c('red','pink1'),las=2)
barplot(t(TMP_TEST_TAB[,c(1:2)]),col=c('red','pink1'),las=2)
dev.off()



TMP_TRAIN_TAB=TRAIN_TAB[which(TRAIN_TAB[,4]==2),]
TMP_TEST_TAB=TEST_TAB[which(TEST_TAB[,4]==2),]
TMP_ORDER=order(-TRAIN_TAB[which(TRAIN_TAB[,4]==2),3])

TMP_TRAIN_TAB=TMP_TRAIN_TAB[TMP_ORDER,]
TMP_TEST_TAB=TMP_TEST_TAB[TMP_ORDER,]

pdf('plot/p02_barplot_second.pdf',width=6,height=4)
barplot(t(TMP_TRAIN_TAB[,c(1:2)]),col=c('royalblue3','skyblue1'),las=2)
barplot(t(TMP_TEST_TAB[,c(1:2)]),col=c('royalblue3','skyblue1'),las=2)
dev.off()





########################
library(pROC)
library(MuSiC)
library(Biobase)
library(Seurat)

source('fwp.R')
source('fwp_source.R')

TYPE=c('BRCA','KIRP','LIHC','PCPG')
DIR='./data/first/'

OUT=.testAll(TYPE, DIR, .test_fwp, npcs=2,mode='g')
saveRDS(OUT, paste0(DIR,'/','eval_fwp.g.rds'))

OUT=.testAll(TYPE, DIR, .test_fwp, mode='l')
saveRDS(OUT, paste0(DIR,'/','eval_fwp.l.rds'))

Og=readRDS(paste0(DIR,'/','eval_fwp.g.rds'))
Ol=readRDS(paste0(DIR,'/','eval_fwp.l.rds'))
O1=readRDS(paste0(DIR,'/','eval_fwp.rds'))
O2=readRDS(paste0(DIR,'/','eval_fwo.rds'))

OOO=cbind(O2[,4],Og[,4],Ol[,4],O1[,4])
colnames(OOO)=c('fwo','global','local','fwp')

source('source.R')
library('ComplexHeatmap')
library('circlize')
train_tag=O2[,1]
test_tag=O2[,2]
ha= rowAnnotation(train=train_tag,test=test_tag,
     col = list(train = c("BRCA" = "indianred2","KIRP" = "gold1",
                          'LIHC'='seagreen3','PCPG'='royalblue1'),
                test  = c("BRCA" = "indianred2","KIRP" = "gold1",
                          'LIHC'='seagreen3','PCPG'='royalblue1')
                ),
     simple_anno_size = unit(0.7, "cm"), gap = unit(1, "mm") )

ht=.drawHeatmap(OOO,show_column_names= T, cluster_columns=F,right_annotation=ha)

pdf('plot/p03_first_gl_heatmap.pdf',width=4,height=3)
print(ht)
dev.off()


X=OOO[,2]
shapiro.test(X)

wilcox.test(OOO[,2],OOO[,1],paired=T,alternative='greater')
wilcox.test(OOO[,3],OOO[,1],paired=T,alternative='greater')
wilcox.test(OOO[,4],OOO[,1],paired=T,alternative='greater')
0.1897
0.001188
0.0005035


wilcox.test(OOO[,1],OOO[,2],paired=T,alternative='greater')
0.8269
wilcox.test(OOO[,3],OOO[,2],paired=T,alternative='greater')
0.08768
wilcox.test(OOO[,4],OOO[,2],paired=T,alternative='greater')
0.05229


wilcox.test(OOO[,1],OOO[,3],paired=T,alternative='greater')
0.999
wilcox.test(OOO[,2],OOO[,3],paired=T,alternative='greater')
0.9205
wilcox.test(OOO[,4],OOO[,3],paired=T,alternative='greater')
0.3991


wilcox.test(OOO[,1],OOO[,4],paired=T,alternative='greater')
0.9996
wilcox.test(OOO[,2],OOO[,4],paired=T,alternative='greater')
0.9533
wilcox.test(OOO[,3],OOO[,4],paired=T,alternative='greater')
0.6226












