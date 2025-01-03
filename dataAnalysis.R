# Data analysis-----------------------------------
data <- read_excel(path = "D:/desktop/1_ToDoList/10_wt论文_图片修改/2023_11_20_wt原始数据/3.1.CircNovel/3.3.CircDiff/normalized_expression/count_TPM.xls")
data <- as.data.frame(data)
rownames(data) <- data$circRNA.readcount
data <- data[,-1]
count <- data[,1:6]
TPM <- data[,7:12]
colnames(count) <- c('Control_1','Control_2','Control_3','IR_1','IR_2','IR_3')
colnames(TPM) <- c('Control_1','Control_2','Control_3','IR_1','IR_2','IR_3')
colData <- data.frame(row.names = colnames(count),
                      name_list = colnames(count),
                      group = factor(rep(c("Control","IR"), each = 3))
)


group <- colData$group

exp <- "IR"
ctr <- "Control"

library(DESeq2)
if(T){
  dds <- DESeqDataSetFromMatrix(countData = count,
                                colData = colData,
                                design = ~ group)
}

dds$group <- relevel(dds$group, ref = exp)

keep <- rowSums(counts(dds)) >= 1.5*ncol(count)

table(keep)
# FALSE  TRUE 
# 4765  9708
dds <- dds[keep,]
dds <- DESeq(dds,quiet = F) 
res <- results(dds, contrast=c("group", exp, ctr))
resOrdered <- res[order(res$padj),]  
DEG <- as.data.frame(resOrdered)

DEG$sig <- as.factor(ifelse(DEG$pvalue < 0.05 & abs(DEG$log2FoldChange) > 1,
                            ifelse(DEG$log2FoldChange > 1 ,'Up','Down'),'None'))

# levels(DEG$sig) <- c('Down','Up','None')

table(DEG$sig)
# Down None Up 
# 82 9566   44

DEG$gene <- rownames(DEG)
DEG <- DEG %>%
  group_by(sig) %>%
  arrange(pvalue, .by_group = TRUE)



save(DEG, file="DEG.RData")
setwd('D:/desktop/1_ToDoList/10_wt论文_图片修改/1_热图火山图')
write.csv(DEG, file = 'DEG.csv')

# Figure supplementary 3A
library(EnhancedVolcano)
p <- EnhancedVolcano(DEG,
                     lab = DEG$gene,
                     # selectLab = c('mmu_circ_0000823','mmu_circ_0000850','mmu_circ_0001111','mmu_circ_0000505','mmu_circ_0000239'),
                     selectLab = c('mmu_circ_0000823'),
                     x = 'log2FoldChange',
                     y = 'pvalue',
                     pCutoff = 0.05,
                     FCcutoff = 1,
                     pointSize = 2, 
                     labSize = 4, 
                     title = 'IR vs. Control',
                     legendPosition = 'right',
                     legendLabSize = 12,
                     legendIconSize = 4.0,
                     drawConnectors = TRUE,
                     widthConnectors = 0.25,
                     col = c("#bebbb8","#66c2a5","#3288bd","#f46d43"))+
  theme(
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank()   
  )+ggtitle("pCutoff = 0.05 & log2FCcutoff = 1")





circRNA <- c('mmu_circ_0000823','mmu_circ_0000850','mmu_circ_0001111',
             'mmu_circ_0000505','mmu_circ_0000239','mmu_circ_0001382',
             'mmu_circ_0000095','mmu_circ_0000597','mmu_circ_0001406',
             'mmu_circ_0001183')
TPM.hp <- TPM[circRNA,]
TPM.hp[1,3] <- 42524.53629
write.csv(TPM.hp, file = 'TPM.hp.csv')

# Figure 2B
pheatmap(TPM.hp,
         scale = "row", 
         border="white", 
         cluster_cols = F,
         cluster_rows = F,
         angle_col = 45,
         display_numbers = F,
         fontsize_col = 8,
         fontsize_row = 8,
         color = colorRampPalette(c("#0da9ce", "white", "#e74a32"))(100))

