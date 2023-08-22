########################
library(pROC)
library(MuSiC)
library(Biobase)
library(Seurat)

source('fwp.R')
source('fwp_source.R')





########################
library(pROC)
library(MuSiC)
library(Biobase)
library(Seurat)

source('fwp.R')
source('fwp_source.R')

TYPE=c('BRCA','KIRP','LIHC','PCPG')
DIR='./data/first/'


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

ht=.drawHeatmap(OOO,show_column_names= T, cluster_columns=F,right_annotation=ha,
    cell_fun = function(j, i, x, y, width, height, fill) {
        grid.text(sprintf("%.3f", OOO[i, j]), x, y, gp = gpar(fontsize = 10))
        })

pdf('revisePlot/p02_first_gl_heatmap_withValue.pdf',width=4.5,height=5)
print(ht)
dev.off()












