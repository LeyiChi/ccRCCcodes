
dir.create('results')
suppressPackageStartupMessages({
    library(tidyverse)
    library(MOVICS)
    library(data.table)
})
mo.data   <- readRDS('../01.bulk.RNAseq/TCGA_multiomics/results/mo.data.rds')
#tcga.cli  <- read.table('../01.bulk.RNAseq/TCGA/results/tcga.cli.txt') %>% select(fustat=OS,futime=OS.time)
tcga.cli  <- read.table('../01.bulk.RNAseq/TCGA_multiomics/results/tcga.cli.txt')
colnames(tcga.cli)[1:2] = c("fustat","futime") 
tcga.cli <- tcga.cli[,c("fustat","futime")]

### This function provides several methods to help selecting elites from input features, which aims to reduce data dimension for multi-omics integrative clustering analysis.
mo.data$mRMA <- mo.data$mRNA.exp %>% getElites(method= "mad",na.action = "rm",elite.num = 1000) %>% pluck('elite.dat') %>% 
   getElites(method= "cox",surv.info = tcga.cli,p.cutoff = 0.05) %>% pluck('elite.dat')

mo.data$lncRNA <- mo.data$lncRNA.exp %>% getElites(method= "mad",na.action = "rm",elite.num = 1000) %>% pluck('elite.dat') %>% 
   getElites(method= "cox",surv.info = tcga.cli,p.cutoff = 0.05) %>% pluck('elite.dat')

mo.data$miRMA <- mo.data$miRMA.exp %>% getElites(method= "mad",na.action = "rm",elite.num = 1000) %>% pluck('elite.dat') %>% 
   getElites(method= "cox",surv.info = tcga.cli,p.cutoff = 0.05) %>% pluck('elite.dat')

mo.data$meth <- mo.data$meth.beta %>% getElites(method= "mad",na.action = "rm",elite.num = 1000) %>% pluck('elite.dat') %>% 
   getElites(method= "cox",surv.info = tcga.cli,p.cutoff = 0.05) %>% pluck('elite.dat')

mo.data$mut <- mo.data$mut.status %>% as.matrix %>% getElites(method= "freq",na.action = "rm",elite.pct = 0.05) %>% pluck('elite.dat')

mo.data <- mo.data[c("mRMA","lncRNA","miRMA","meth","mut")]

load("../2.MOVICS/all data.RData")
optk <- getClustNum(data = mo.data,
                         is.binary = c(F,F,F,F,T), 
                         try.N.clust = 2:8, 
                         fig.path = 'results',
                         fig.name = "Fig.S2")

# optk <- getClustNum(data = mo.data,
#                     is.binary = c(F,F,F,F,T), 
#                     try.N.clust = 2:8, 
#                     fig.path = 'results',
#                     fig.name = "Fig.S2")

# methodslist = list("SNF", "CIMLR", "PINSPlus", "NEMO", "COCA", "MoCluster",
#                    "LRAcluster", "ConsensusClustering", "IntNMF", "iClusterBayes")

moic.res.list <- getMOIC(data = mo.data, 
                         N.clust = optk$N.clust,
                         type = c("gaussian", "gaussian", "gaussian", "gaussian" ,"binomial"))

cmoic <- getConsensusMOIC(moic.res.list = moic.res.list,
                               distance = "euclidean", 
                               linkage  = "average",
                               fig.path = 'results',
                               fig.name = "Fig.2C")
saveRDS(moic.res.list,"results/moic.res.list.rds",compress =F)
saveRDS(cmoic,"results/cmoic.rds",compress =F)

moic.res.list<-read_rds('results/moic.res.list.rds')
cmoic<-read_rds('results/cmoic.rds')

####### Fig.S3-SILHOUETTE##############################
getSilhouette(sil = cmoic$sil,
              fig.path = 'results',
              fig.name = "Fig.S3")

##### Fig.2B-10 CLUSTERING METHODS################################
cmoic$clust.res %>% select(clust) %>% write.table("results/tcga.clust.txt",sep="\t",quote = F)

clust.res <- cmoic$clust.res$clust

for (moic.res in moic.res.list){
     clust.res <- cbind(clust.res,moic.res$clust.res$clust)
}
colnames(clust.res) <- c('Subtype',names(moic.res.list))

col = c("#2EC4B6", "#E71D36", "#FF9F1C", "#BDD5EA", "#FFA5AB", "#011627","#023E8A", "#9D4EDD")
names(col) <- 1:length(col)

clust.res = as.data.frame(clust.res) %>% arrange(Subtype)
anno = ComplexHeatmap::HeatmapAnnotation(df=clust.res,show_legend = F,col=list(Subtype= col))
library(ComplexHeatmap)
pdf(file = 'results/Fig.2B.pdf',width = 5,height = 3)
Heatmap(matrix(nrow = 0, ncol = length(clust.res$Subtype)),top_annotation =anno)
dev.off()

#### Fig.2D-KAPLAN-MEIER CURVE OF CONSENSUSMOIC#################################
surv <- compSurv(moic.res         = cmoic,
                      surv.info        = tcga.cli,
                      convt.time       = "m",
                      surv.median.line = "h", 
                      xyrs.est         = c(5,10),
                      fig.path      = 'results',
                      fig.name         = "Fig.2D")



