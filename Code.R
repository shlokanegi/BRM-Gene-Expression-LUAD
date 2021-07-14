## Importing required Packages and Libraries
#---------------------------------------------
install.packages("ggpubr")
install.packages("ggplot2")
install.packages("reshape2")
install.packages("pheatmap")
library(ggplot2)
library(reshape2)
library(pheatmap)
library(ggpubr)

#-------------------------------------------------------------
## Importing LUAD TCGA dataset, downloaded from GDAC Firehose
#-------------------------------------------------------------
x <- LUAD.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data <- read.delim("/Users/ShlokaNegi_1/Desktop/RECORDS/LUAD.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt", row.names=1)
x <- x[-1,] # 20531 genes
y <- x


#-------------------
## Quality Controls 
#----------------- FOR GENES -------------#
"Decided to keep only the genes which have >= 50% of its samples with 
non-zero expression values & samples with > 20% of genes showing non-zero
expression values"
# Introducing NA values, for all cells with 0
for (i in colnames(y)){
  c = 0
  for (j in as.numeric(y[[i]])){
    c = c + 1
    if(j!= 0){
      y[[i]][c] <- NA
    }
  }
}
View(y)
write.table(y,file = "ytable(NA)")

y <- as.matrix(y)
zeros_genes <- apply(y, 1, function(x) length(which(!is.na(x)))) #counting number of 0s per gene
zeros_genes <- as.data.frame(zeros_genes)
zeros_genes <- cbind(rownames(zeros_genes), zeros_genes)
rownames(zeros_genes) <- NULL
colnames(zeros_genes)[1] <- "GeneID"
colnames(zeros_genes)[2] <- "#0s"
head(zeros_genes)

mean_samples <- ncol(x)/2 
exp_genes <- zeros_genes[zeros_genes$`#0s`< mean_samples,] 
head(exp_genes) #Choosing genes with >=50% of its samples having non-zero expression values

x <- cbind(rownames(x), x)
new_data <- x[x$`rownames(x)` %in% exp_genes$GeneID,]
new_data <- new_data[,-1]
head(new_data) # has 17754 genes after QC (removed 2777 genes)

#--------------- FOR SAMPLES ---------------#
zeros_samples <- apply(y, 2, function(x) length(which(!is.na(x)))) # Counting number of 0s per sample
zeros_samples <- as.data.frame(zeros_samples)
colnames(zeros_samples)[1] <- "#0s"
zeros_samples <- cbind(rownames(zeros_samples), zeros_samples)
rownames(zeros_samples) <- NULL
colnames(zeros_samples)[1] <- "Samples"
head(zeros_samples)

mean_genes <- 0.2*nrow(x)
View(new_data)
exp_samples <- zeros_samples[zeros_samples$`#0s`< mean_genes,]
head(exp_samples) ##Choosing samples with >20% of its genes having non-zero expression values

final_data <- new_data[,match(exp_samples$Samples, colnames(new_data))]
View(final_data) # Final dataset - # 17754 genes and 575 samples after QC


#-------------------------------------------
## Splitting Datasets into Normal and Tumor
#-------------------------------------------
dat1 <- final_data
dat1 <- t(dat1) # now in matrix form, with samples as rows
control <- dat1[grep("11A", rownames(dat1), fixed = TRUE),]
control_indices <- grep("11A", rownames(dat1), fixed = TRUE, value = FALSE)
tumor <- dat1[-control_indices,]
View(control) # 58 normal samples
View(tumor) # 517 tumor samples


#------------------------------------
## Studying SMARCA2 gene expression
#------------------------------------
SMARCA2_tumor <- as.data.frame(tumor[,grep("SMARCA2", colnames(tumor), fixed = TRUE)])
colnames(SMARCA2_tumor)[1] <- "value"
SMARCA2_tumor <- cbind(rownames(SMARCA2_tumor), SMARCA2_tumor)
colnames(SMARCA2_tumor)[1] <- "Tumor_samples"
SMARCA2_tumor$value <- as.numeric(SMARCA2_tumor$value)
head(SMARCA2_tumor)

# Histogram of SMARCA2 in tumor samples
vals1 <- as.numeric(SMARCA2_tumor$value)
h_t <- hist(vals1, breaks = 20, col = "grey", labels = TRUE,
            ylim = c(0,200),
            xlab = "SMARCA2 Expression Values", 
            ylab = "Number of Samples", 
            main = "Histogram of SMARCA2 gene in Tumor Samples", 
            xlim = c(0,8000))
xfit<- seq(min(vals1), max(vals1), length = 517)
yfit <- dnorm(xfit, mean = mean(vals1), sd = sd(vals1))
yfit <- yfit*diff(h_t$mids[1:2])*length(vals1)
lines(xfit, yfit, col="red", lwd=3) 
abline(v = mean(vals1), col = "blue", lwd = 3)
abline(v = median(vals1), col = "dark green", lwd = 3)

legend(4000, 180, legend = c("Normal Curve", "Mean = 1622.098", 
       "Median = 1489.089"), col = c("red", "blue", "dark green"), 
       lwd = 3)
mean(vals1) # 1622.098
median(vals1) # 1489.089

# Histogram of SMARCA2 in normal samples
SMARCA2_control <- as.data.frame(control[,grep("SMARCA2", colnames(control), fixed = TRUE)])
colnames(SMARCA2_control)[1] <- "value"
SMARCA2_control <- cbind(rownames(SMARCA2_control), SMARCA2_control)
colnames(SMARCA2_control)[1] <- "Control_samples"
SMARCA2_control$value <- as.numeric(SMARCA2_control$value)
View(SMARCA2_control)

vals2 <- as.numeric(SMARCA2_control$value)
h_c <- hist(vals2, breaks = 20, col = "grey", labels = TRUE,
            ylim = c(0,15),
            xlab = "SMARCA2 Expression Values", 
            ylab = "Number of Samples", 
            main = "Histogram of SMARCA2 gene in Control Samples", 
            xlim = c(1000,6000))
xfit<- seq(min(vals2), max(vals2), length = 58)
yfit <- dnorm(xfit, mean = mean(vals2), sd = sd(vals2))
yfit <- yfit*diff(h_c$mids[1:2])*length(vals2)
lines(xfit, yfit, col="red", lwd=3) 
abline(v = mean(vals2), col = "blue", lwd = 3)
abline(v = median(vals2), col = "dark green", lwd = 3)
legend(4000, 15, legend = c("Normal Curve", "Mean = 3069.014", 
       "Median = 3195.834"), col = c("red", "blue", "dark green"), 
       lwd = 3)
mean(vals2) # 3069.014
median(vals2) # 3195.834


#---------------------------------------------
## Categorizing samples into BRM-H and BRM-L
#---------------------------------------------
'The "1st quartile" is the median of the first half of the data set 
and marks the point at which 25% of the data values are lower and 
75% are higher. 
The "3rd quartile" is the median of the second half of the data set 
and marks the point at which 25% of the data values are higher and 
75% lower.'

quantile(vals1, probs = seq(0, 1, 0.25), names = TRUE, type = 8)

h_s <- hist(vals1, breaks = 20, col = "grey", labels = TRUE,
            ylim = c(0,200),
            xlab = "SMARCA2 Expression Values", 
            ylab = "Number of Samples", 
            main = "Tumor samples divided into quartiles", 
            xlim = c(0,8000))
abline(v = 1088.2563, col = "blue", lwd = 3)
abline(v = 1489.0892, col = "black", lwd = 3)
abline(v = 2020.9780, col = "dark green", lwd = 3)
legend(3000, 200, legend = c("Normal Curve", "1st quartile 25%", 
                          "Median = 50%", "3rd quartile 75%"), 
       col = c("red", "blue", "black", "dark green"), 
       lwd = 3)
                          
q1 <- SMARCA2_tumor[SMARCA2_tumor$value<=1088.2563,]
q2 <- SMARCA2_tumor[SMARCA2_tumor$value>1088.2563 & SMARCA2_tumor$value<2020.9780,]
q3 <- SMARCA2_tumor[SMARCA2_tumor$value>=2020.9780,]
View(q1)  # 129 samples -> BRM-Low category
View(q2) # 259 samples
View(q3) # 129 samples -> BRM-High category

q1$BRM_level <- factor(rep("low", 129))
q3$BRM_level <- factor(rep("high", 129))
SMARCA2_new <- rbind(q1, q3)
rownames(SMARCA2_new) <- NULL

#---------------------------------------------------------
## Segregating tumor data based on BRM-L and BRM-H groups
#---------------------------------------------------------
# For BRM-Low Category
temp <- as.data.frame(tumor[match(q1$Tumor_samples, rownames(tumor)),])
g <- grep("SMARCA2", colnames(temp), fixed = TRUE, value = FALSE)
temp <- temp[,-g] #remove SMARCA2 gene data
temp <- cbind(q1$BRM_level, temp)
colnames(temp)[1] = "BRM_level"
View(temp)

# For BRM-High Category
temp2 <- as.data.frame(tumor[match(q3$Tumor_samples, rownames(tumor)),])
temp2 <- temp2[,-g]
temp2 <- cbind(q3$BRM_level, temp2)
colnames(temp2)[1] = "BRM_level"
View(temp2)

data1 <- rbind(temp, temp2)
View(data1) # Contains BRM-High (n = 129) and BRM-low samples (n = 129)


#---------------------------------------------------------
## Identifying genes associated with BRM expression
#---------------------------------------------------------
"Performing non-parametric test (Mann Whitney Wilcoxon Test) for 
identifying associations between a gene and BRM-levels"

pvals <- vector(mode = "numeric", length = 17753)
for(i in 1:17753){
        r <- wilcox.test(as.numeric(data1[,i+1]) ~ BRM_level, data = data1)
        pvals[i] <- r$p.value
}
pvals <- as.data.frame(pvals)
genes <- colnames(data1)[-1]
pvals <- cbind(genes, pvals)
View(pvals)

# Identifying significant genes, associated with BRM-levels
"According to The Bonferroni Correction (correction for multiple 
testing), we can set the threshold as, 0.05/17753 = 2.816425e-06
0 --> non-significant; p >= 2.816425e-06
1 --> significant; p < 2.816425e-06"

signif <- ifelse(pvals$pvals >= 2.816425e-06, 0, 1)
View(signif)
ones <- 0  # 6034 significant genes
zeros <- 0 # 11719 non-significant genes
for(i in signif){
  if(i == 1){
    ones = ones+1
  }else{
    zeros = zeros + 1
  }
}
print(paste("significant genes", ones))
print(paste("Non-significant genes", zeros))

pvals <- cbind(pvals, signif)
pvals<- pvals[order(pvals$pvals),] #Ranking the genes, based on their p-values
View(pvals) #contains ranked genes

signif_genes <- vector(mode = "character", length = length(ones))
for(i in 1:nrow(pvals)){
        if(pvals[i,3] == 1){
            signif_genes= c(signif_genes, pvals[i,1])
        }
}
signif_genes <- signif_genes[-1] #remove the first empty string
signif_genes <- as.data.frame(signif_genes)
data1[,-1] <- lapply(data1[,-1], as.numeric) # change data type of each column values from character to numeric
signif_index <- match(signif_genes$signif_genes, colnames(data1))
signif_genes <- data1[,signif_index] 
signif_genes <- cbind(data1$BRM_level, signif_genes) # now signif_genes contains only ranked significant genes data 
colnames(signif_genes)[1] <- "BRM_level"
View(signif_genes)


#-------------------------------------------------------
## Plotting top 15 significant genes (till p-value E-29)
#-------------------------------------------------------
# DIABLO
compare_means(`DIABLO|56616` ~ BRM_level, data = signif_genes)
p <- ggboxplot(data = signif_genes, x = "BRM_level", 
               y = "DIABLO|56616", color = "BRM_level", 
               palette = "lancet", add = "jitter", 
               ylab = "DIABLO gene expression", 
               main = "Boxplot of DIABLO")
p + stat_compare_means()

# PIK3R1
compare_means(`PIK3R1|5295` ~ BRM_level, data = signif_genes)
p <- ggboxplot(data = signif_genes, x = "BRM_level", 
               y = "PIK3R1|5295", color = "BRM_level", 
               palette = "lancet", add = "jitter", 
               ylab = "PIK3R1 gene expression", 
               main = "Boxplot of PIK3R1")
p + stat_compare_means()

# SMUG1
compare_means(`SMUG1|23583` ~ BRM_level, data = signif_genes)
p <- ggboxplot(data = signif_genes, x = "BRM_level", 
               y = "SMUG1|23583", color = "BRM_level", 
               palette = "lancet", add = "jitter", 
               ylab = "SMUG1 gene expression", 
               main = "Boxplot of SMUG1")
p + stat_compare_means()

# NFIX
compare_means(`NFIX|4784` ~ BRM_level, data = signif_genes)
p <- ggboxplot(data = signif_genes, x = "BRM_level", 
               y = "NFIX|4784", color = "BRM_level", 
               palette = "lancet", add = "jitter", 
               ylab = "NFIX gene expression", 
               main = "Boxplot of NFIX")
p + stat_compare_means()

# MRPL12
compare_means(`MRPL12|6182` ~ BRM_level, data = signif_genes)
p <- ggboxplot(data = signif_genes, x = "BRM_level", 
               y = "MRPL12|6182", color = "BRM_level", 
               palette = "lancet", add = "jitter", 
               ylab = "MRPL12 gene expression", 
               main = "Boxplot of MRPL12")
p + stat_compare_means()

# SNRPF
compare_means(`SNRPF|6636` ~ BRM_level, data = signif_genes)
p <- ggboxplot(data = signif_genes, x = "BRM_level", 
               y = "SNRPF|6636", color = "BRM_level", 
               palette = "lancet", add = "jitter", 
               ylab = "SNRPF gene expression", 
               main = "Boxplot of SNRPF")
p + stat_compare_means()

# TUBA1C
compare_means(`TUBA1C|84790` ~ BRM_level, data = signif_genes)
p <- ggboxplot(data = signif_genes, x = "BRM_level", 
               y = "TUBA1C|84790", color = "BRM_level", 
               palette = "lancet", add = "jitter", 
               ylab = "TUBA1C gene expression", 
               main = "Boxplot of TUBA1C")
p + stat_compare_means()

# RASSF5
compare_means(`RASSF5|83593` ~ BRM_level, data = signif_genes)
p <- ggboxplot(data = signif_genes, x = "BRM_level", 
               y = "RASSF5|83593", color = "BRM_level", 
               palette = "lancet", add = "jitter", 
               ylab = "RASSF5 gene expression", 
               main = "Boxplot of RASSF5")
p + stat_compare_means()

# C12orf73
compare_means(`C12orf73|728568` ~ BRM_level, data = signif_genes)
p <- ggboxplot(data = signif_genes, x = "BRM_level", 
               y = "C12orf73|728568", color = "BRM_level", 
               palette = "lancet", add = "jitter", 
               ylab = "C12orf73 gene expression", 
               main = "Boxplot of C12orf73")
p + stat_compare_means()

# BOLA3
compare_means(`BOLA3|388962` ~ BRM_level, data = signif_genes)
p <- ggboxplot(data = signif_genes, x = "BRM_level", 
               y = "BOLA3|388962", color = "BRM_level", 
               palette = "lancet", add = "jitter", 
               ylab = "BOLA3 gene expression", 
               main = "Boxplot of BOLA3")
p + stat_compare_means()

# COLEC12
compare_means(`COLEC12|81035` ~ BRM_level, data = signif_genes)
p <- ggboxplot(data = signif_genes, x = "BRM_level", 
               y = "COLEC12|81035", color = "BRM_level", 
               palette = "lancet", add = "jitter", 
               ylab = "COLEC12 gene expression", 
               main = "Boxplot of COLEC12")
p + stat_compare_means()

# CKS1B
compare_means(`CKS1B|1163` ~ BRM_level, data = signif_genes)
p <- ggboxplot(data = signif_genes, x = "BRM_level", 
               y = "CKS1B|1163", color = "BRM_level", 
               palette = "lancet", add = "jitter", 
               ylab = "CKS1B gene expression", 
               main = "Boxplot of CKS1B")
p + stat_compare_means()

# DTYMK|1841
compare_means(`DTYMK|1841` ~ BRM_level, data = signif_genes)
p <- ggboxplot(data = signif_genes, x = "BRM_level", 
               y = "DTYMK|1841", color = "BRM_level", 
               palette = "lancet", add = "jitter", 
               ylab = "DTYMK gene expression", 
               main = "Boxplot of DTYMK")
p + stat_compare_means()

# AKAP13
compare_means(`AKAP13|11214` ~ BRM_level, data = signif_genes)
p <- ggboxplot(data = signif_genes, x = "BRM_level", 
               y = "AKAP13|11214", color = "BRM_level", 
               palette = "lancet", add = "jitter", 
               ylab = "AKAP13 gene expression", 
               main = "Boxplot of AKAP13")
p + stat_compare_means()

# MCC
compare_means(`MCC|4163` ~ BRM_level, data = signif_genes)
p <- ggboxplot(data = signif_genes, x = "BRM_level", 
               y = "MCC|4163", color = "BRM_level", 
               palette = "lancet", add = "jitter", 
               ylab = "MCC gene expression", 
               main = "Boxplot of MCC")
p + stat_compare_means()


#-------------------------------------------------------
## Perform Functional Annotation Clustering using DAVID 
#-------------------------------------------------------
geneIDs <- ""
for (i in colnames(signif_genes)[-1]){
        s <- strsplit(i,"|")
        s <- unlist(s)
        m <- match("|", s)
        m <- m+1
        E_ID <- paste(s[m:length(s)], collapse = "")
        geneIDs <- c(geneIDs, E_ID)
}
geneIDs <- geneIDs[-1]
write(geneIDs, file = "signif_geneIDs.txt")
length(geneIDs) # 6034 significant genes

write(geneIDs[1:3000], file = "signif_geneIDs_3000.txt")


#--------------------------------------------------------------------------
## Comparing Tumor v/s Normal gene expression for top 15 significant genes
#--------------------------------------------------------------------------
tumor <- as.data.frame(tumor)
control <- as.data.frame(control)

mean(as.numeric(tumor$`RASSF5|83593`)) #1046.543
mean(as.numeric(control$`RASSF5|83593`)) #1619.791

mean(as.numeric(tumor$`C12orf73|728568`)) #136.8939
mean(as.numeric(control$`C12orf73|728568`)) #92.67397

mean(as.numeric(tumor$`COLEC12|81035`)) #827.1525
mean(as.numeric(control$`COLEC12|81035`)) #2548.703

mean(as.numeric(tumor$`AKAP13|11214`)) #4178.481
mean(as.numeric(control$`AKAP13|11214`)) #8544.608

mean(as.numeric(tumor$`MCC|4163`)) #265.4271
mean(as.numeric(control$`MCC|4163`)) #642.501


