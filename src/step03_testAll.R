########################
library(pROC)
library(MuSiC)
library(Biobase)
library(Seurat)

source('fwp.R')
source('fwp_source.R')

TYPE=c('BRCA','KIRP','LIHC','PCPG')
DIR='./data/first/'

OUT=.testAll(TYPE, DIR, .test_fwp, npcs=1)
saveRDS(OUT, paste0(DIR,'/','eval_fwp.1.rds'))
OUT=.testAll(TYPE, DIR, .test_fwp, npcs=2)
saveRDS(OUT, paste0(DIR,'/','eval_fwp.rds'))
OUT=.testAll(TYPE, DIR, .test_fwp, npcs=3)
saveRDS(OUT, paste0(DIR,'/','eval_fwp.3.rds'))
OUT=.testAll(TYPE, DIR, .test_fwp, npcs=4)
saveRDS(OUT, paste0(DIR,'/','eval_fwp.4.rds'))
OUT=.testAll(TYPE, DIR, .test_fwp, npcs=5)
saveRDS(OUT, paste0(DIR,'/','eval_fwp.5.rds'))
OUT=.testAll(TYPE, DIR, .test_fwo)
saveRDS(OUT, paste0(DIR,'/','eval_fwo.rds'))
OUT=.testAll(TYPE, DIR, .test_music)
saveRDS(OUT, paste0(DIR,'/','eval_music.rds'))
OUT=.testAll(TYPE, DIR, .test_ikarus)
saveRDS(OUT, paste0(DIR,'/','eval_ikarus.rds'))

DIR='./data/first/'
O1.1=readRDS(paste0(DIR,'/','eval_fwp.1.rds'))
O1=readRDS(paste0(DIR,'/','eval_fwp.rds'))
O1.3=readRDS(paste0(DIR,'/','eval_fwp.3.rds'))
O1.4=readRDS(paste0(DIR,'/','eval_fwp.4.rds'))
O1.5=readRDS(paste0(DIR,'/','eval_fwp.5.rds'))
O2=readRDS(paste0(DIR,'/','eval_fwo.rds'))
O3=readRDS(paste0(DIR,'/','eval_music.rds'))
O4=readRDS(paste0(DIR,'/','eval_ikarus.rds'))

OOO=cbind(O1.1[,4],O1[,4],O1.3[,4],O1.4[,4],O1.5[,4],O2[,4],O3[,4],O4[,4])

pdf('plot/p04_boxplot_first.pdf',width=4,height=3.5)
COL=c(rep('indianred2',5),'royalblue1','grey70','grey70')
boxplot(OOO, pch='+',col=COL,ylim=c(0,1))
abline(h=0.5,lty=2)
dev.off()

apply(OOO,2,mean)
[1] 0.9181632 0.9262461 0.9169794 0.9181807 0.9180131 0.8920425 0.8813760
[8] 0.85379

DIR='./data/first/'
O1=readRDS(paste0(DIR,'/','eval_fwp.rds'))
O2=readRDS(paste0(DIR,'/','eval_fwo.rds'))
plot(O1[,4],O2[,4])
abline(a=0,b=1)

mean(O1[,4])
mean(O2[,4])
median(O1[,4])
median(O2[,4])

t.test(O1[,4],O2[,4],paired=T)
wilcox.test(O1[,4],O2[,4],paired=T)



plot(TMP,type='h',lwd=10,col='grey50')


O3=readRDS(paste0(DIR,'/','eval_music.rds'))
plot(O1[,4],O3[,4])
abline(a=0,b=1)

O4=readRDS(paste0(DIR,'/','eval_ikarus.rds'))
plot(O1[,4],O4[,4])
abline(a=0,b=1)




source('fwp.R')
source('fwp_source.R')

TYPE=c('BLCA','CESC','CHOL','COAD','ESCA','HNSC','KICH','LUSC','PAAD','THYM','UCEC')
DIR='./data/second/'

OUT=.testAll(TYPE, DIR, .test_fwp,npcs=2)
saveRDS(OUT, paste0(DIR,'/','eval_fwp.rds'))
OUT=.testAll(TYPE, DIR, .test_fwo)
saveRDS(OUT, paste0(DIR,'/','eval_fwo.rds'))
OUT=.testAll(TYPE, DIR, .test_music)
saveRDS(OUT, paste0(DIR,'/','eval_music.rds'))
OUT=.testAll(TYPE, DIR, .test_ikarus)
saveRDS(OUT, paste0(DIR,'/','eval_ikarus.rds'))

DIR='./data/second/'
O1=readRDS(paste0(DIR,'/','eval_fwp.rds'))
O2=readRDS(paste0(DIR,'/','eval_fwo.rds'))
O3=readRDS(paste0(DIR,'/','eval_music.rds'))
O4=readRDS(paste0(DIR,'/','eval_ikarus.rds'))

OOO=cbind(O1[,4],O2[,4],O3[,4],O4[,4])

pdf('plot/p05_boxplot_second.pdf',width=3,height=3.5)
COL=c(rep('indianred2',1),'royalblue1','grey70','grey70')
boxplot(OOO, pch='+',col=COL,ylim=c(0,1))
abline(h=0.5,lty=2)
dev.off()

apply(OOO,2,mean)
#[1] 0.8817769 0.8660364 0.8505413 0.7008754 


C.INDEX=which(O1[,1]==O1[,2])
I.INDEX=which(O1[,1]!=O1[,2])
pdf('plot/p12_boxplot_second_consistent.pdf',width=3,height=3.5)
COL=c(rep('indianred2',1),'royalblue1','grey70','grey70')
boxplot(OOO,pch='+',col=COL,ylim=c(0,1))
abline(h=0.5,lty=2)
boxplot(OOO[C.INDEX,],pch='+',col=COL,ylim=c(0,1))
abline(h=0.5,lty=2)
boxplot(OOO[I.INDEX,],pch='+',col=COL,ylim=c(0,1))
abline(h=0.5,lty=2)
dev.off()

apply(OOO[C.INDEX,],2,mean)
# 0.9705992 0.9697952 0.9463702 0.8586975
apply(OOO[I.INDEX,],2,mean)
# 0.8728947 0.8556605 0.8409584 0.6850931





DIR='./data/second/'
O1=readRDS(paste0(DIR,'/','eval_fwp.rds'))
O2=readRDS(paste0(DIR,'/','eval_fwo.rds'))

DIFF=O1[,4]-O2[,4]
names(DIFF)=paste0(O1[,1],'_',O1[,2])
O.DIFF=DIFF[order(abs(DIFF))]
O.COL=rep('indianred2',length(O.DIFF))
O.COL[which(O.DIFF<0)]='royalblue1'
O.COL[which(O.DIFF==0)]='grey70'

pdf('plot/p06_diff_second.pdf',width=4,height=4)
plot(O.DIFF,col=O.COL,type='h',lwd=2,ylim=c(-0.3,0.3))
abline(h=0,lty=2)
dev.off()


length(O.DIFF[which(O.DIFF>0.05 )])
length(O.DIFF[which(O.DIFF< -0.05)])

table(O.COL)/length(O.COL)


# grey70 indianred2 royalblue1
# 0.2479339  0.5123967  0.2396694


#plot(O1[,4],O2[,4])
#abline(a=0,b=1)

mean(O1[,4])
mean(O2[,4])

t.test(O1[,4],O2[,4],paired=T,alternative='greater')
#0.001005
wilcox.test(O1[,4],O2[,4],paired=T,alternative='greater')
#0.0001199






