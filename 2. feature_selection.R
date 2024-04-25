
rm(list = ls())
setwd("F:/007-research/002-papers/[2024-01]/224.肾癌MOVICS")

library(dplyr)
library(survival)
library(Boruta)
library(future.apply)
library(survivalROC)

repstress.gene<-read.csv(file = "code/marker.up.txt", header = FALSE)
names(repstress.gene) <- "Gene"
dfs <- read.csv(file = "code/tcga_cli.txt", header = TRUE, sep = "\t")
names(dfs) <- c("Sample ID", "Age", "Gender", "T.Stage", "N.Stage", "M.Stage", 
                "Stage", "DFS_Status", "DFS", "Age1")

load("./code/Bulkseq.RData")

# tcgaexp_dfs <- read.csv(file = "01.bulk.RNAseq/Bulkdatpre/1.TCGA/results/tcga_dat_T.txt", sep = "\t")

tcgaexp_dfs <- TCGA

# tcgaexp_dfs <- tcgaexp_dfs %>% t() %>% 
#   as.data.frame() %>%
#   tibble::rownames_to_column("Sample ID") 

tcgaexp_dfs$`Sample ID` <- row.names(tcgaexp_dfs)

repstress.gene <- unique(repstress.gene$Gene)
intersect(colnames(tcgaexp_dfs),repstress.gene)

tcga_rep <- dfs %>%
  inner_join(tcgaexp_dfs[,c("Sample ID",intersect(colnames(tcgaexp_dfs),repstress.gene))],by= "Sample ID")


##  univariate cox regression--------------------------------------------------------------------
res <- data.frame()
data <- tcga_rep
genes <- colnames(data)[11:ncol(data)]

for (i in 1:length(genes)) {
  print(genes[i])
  surv =as.formula(paste('Surv(DFS, DFS_Status)~', genes[i]))
  x = coxph(surv, data = data)
  x = summary(x)
  coef=x$coef[1]
  p.value=signif(x$wald["pvalue"], digits=2)
  HR =signif(x$coef[2], digits=2);#exp(beta)
  HR.confint.lower = signif(x$conf.int[,"lower .95"], 2)
  HR.confint.upper = signif(x$conf.int[,"upper .95"],2)
  CI <- paste0("(", 
               HR.confint.lower, "-", HR.confint.upper, ")")
  res[i,1] = genes[i]
  res[i,2] = coef
  res[i,3] = HR
  res[i,4] = CI
  res[i,5] = p.value
}

names(res) <- c("gene","coef","HR","95% CI","p.value")
res.sig <- filter(res, p.value <0.01)

### Bootstrap
outTab <- NULL
surv <- tcga_rep[,c("DFS", "DFS_Status",res.sig$gene)]
rownames(surv) <- surv$`Patient ID`
for(i in 3:ncol(surv)){ # survival information (OS in this case)
  
  # customized function
  display.progress = function (index, totalN, breakN=20) {
    if ( index %% ceiling(totalN/breakN)  ==0  ) {
      cat(paste(round(index*100/totalN), "% ", sep=""))
    }
  }    
  
  display.progress(index = i, totalN = ncol(surv)) # show running progression
  gene <- colnames(surv)[i]
  Mboot <- future_replicate(1000, expr = { # bootstrap for 1,000 times
    indices <- sample(rownames(surv), size = nrow(surv) * 0.8, replace = F) # extract 80% samples at each bootsratp
    data <- surv[indices,]
    fmla1 <- as.formula(Surv(data[,"DFS"],data[,"DFS_Status"]) ~ data[,gene])
    mycox <- coxph(fmla1,data = data)
    coxResult <- summary(mycox)
    P <- coxResult$coefficients[,"Pr(>|z|)"]
  }
  )
  times <- length(Mboot[which(Mboot < 0.01)])
  outTab <- rbind(outTab,
                  cbind(gene = gene,
                        times = times))
}

outTab <- as.data.frame(outTab)

bootGene <- outTab[as.numeric(as.character(outTab$times)) > 800,] 


# Boruta ------------------------------------------------------------------
# boruta =Boruta (DFS_Status ~ ., data = surv[,c("DFS_Status",share)], doTrace=2, ntree=1000, maxRuns = 1000)
boruta =Boruta (DFS_Status ~ ., data = surv[,c("DFS_Status", bootGene$gene)], doTrace=2, ntree=1000, maxRuns = 1000)

table(boruta$finalDecision)

pdf("./code/boruta_plotImpHistory.pdf",height = 4,width = 5)
Boruta::plotImpHistory(boruta)
dev.off()
plot(boruta, whichShadow = c(TRUE, TRUE, TRUE))


help(package="Boruta")
tt <- as.data.frame(boruta$finalDecision)
names(tt) <- "Var"
selected.vars <- c("shadowMax", "shadowMean", "shadowMin", row.names(tt)[tt$Var == "Confirmed"])

ttt <- as.data.frame(boruta$ImpHistory)[, selected.vars]



library(tidyr)
library(dplyr)
library(ggplot2)

# 使用 gather() 函数将宽表转换为长表
boruta.data <- gather(ttt, key = "Gene", value = "Importance")
boruta.data$Group <- c(rep("shadowMax", 999), rep("shadowMean", 999), rep("shadowMin", 999), rep("Confirmed", 999 * 43))

result <- boruta.data %>%
  group_by(Gene) %>%
  summarize(Median = median(Importance))

result <- arrange(result, Median)



# 使用 ggplot2 创建箱线图
ggplot(boruta.data, aes(x = Gene, y = Importance, fill = Group)) +
  geom_boxplot(outlier.colour = "red", outlier.size = 2) +
  theme_bw() +
  labs(title = "", x = "", y = "Importance") +
  scale_x_discrete(limits = result$Gene) +
  theme(title =  element_text(size = rel(1.1), colour = "black"),
        axis.text.x = element_text(size = rel(1.1), colour = "black", angle = 90),
        axis.text.y = element_text(size = rel(1.1), colour = "black"),
        axis.title.x = element_text(size = rel(1.1), colour = "black"),
        axis.title.y = element_text(size = rel(1.1), colour = "black"),
        plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = c(1, 0), legend.justification = c(1, 0),
        panel.grid.major = element_blank(),  # 去除主要网格线
        panel.grid.minor = element_blank()   # 去除次要网格线
  ) +
  theme(
    axis.line.x = element_line(color = "black"),  # 设置 x 轴的颜色
    axis.line.y = element_line(color = "black"),  # 设置 y 轴的颜色
    panel.border = element_blank()  # 隐藏面板的边框
  ) +
  guides(fill = guide_legend(title="", nrow=4,byrow=TRUE)) +
  # guides(fill = guide_legend(reverse = T)) +
  # scale_fill_discrete(breaks=c("NYHA I", "NYHA II", "NYHA III", "NYHA IV"))
  scale_fill_manual(values = c("green", "#F3D660", "darkblue", "#AC60F3"))






