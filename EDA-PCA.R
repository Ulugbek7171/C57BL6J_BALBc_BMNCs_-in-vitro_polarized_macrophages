#Data loading
getwd()
setwd("/Users/surgeon7171/Library/CloudStorage/OneDrive-OUMail(OsakaUniversity)/Research/TAC/RNA-seq_analyses/Data_all/Raw counts_all")

#Loading both data sets (F9808 and F10115);
counts1 <- read.csv("F10115_vx0191_raw_count.csv", row.names = 1)
counts2 <- read.csv("F9808_vx0169_raw_count.csv", row.names = 1)

#Check gene overlap
genes_common <- intersect(rownames(counts1), rownames(counts2))
length(genes_common)

#Subset both to shared genes only:
counts1 <- counts1[genes_common,]
counts2 <- counts2[genes_common,]

#Compare distribution of raw counts
library_size1 <- colSums(counts1)
library_size2 <- colSums(counts2)

#Calculate global correlation between data sets
avg1 <- rowMeans(counts1) #Not necessary for PCA
avg2 <- rowMeans(counts2)

#Next step is merging the data sets and running PCA pipeline again
#Combine them column-wise (assuming genes are aligned)
combined_counts <- cbind(counts1, counts2)
dim(combined_counts)
head(colnames(combined_counts))

#Create a combined group vector
group <- factor(c("BALBc_BM", "BALBc_BM", "BALBc_BM", 
                  "BALBc_B", "BALBc_B", "BALBc_B", 
                  "C57BL6J_BM", "C57BL6J_BM", "C57BL6J_BM",
                  "C57BL6J_B", "C57BL6J_B", "C57BL6J_B")) #THIS GROUPING ALLOWS TO HAVE FOUR GROUPS: 2 MACROPHAGE AND 2 BONE MARROW CELLS, 4 COLORS ARE FOR 4 GROUPS;
group <- factor(c("BALBc_BM_1", "BALBc_BM_2", "BALBc_BM_3", 
                  "BALBc_B_1", "BALBc_B_2", "BALBc_B_3", 
                  "C57BL6J_BM_1", "C57BL6J_BM_2", "C57BL6J_BM_3",
                  "C57BL6J_B_1", "C57BL6J_B_2", "C57BL6J_B_3")) # THIS GROUPING BRINGS ALL SAMPLES SEPARATELY IN THE PLOT WITH 12 DIFFERENT COLORS; 

#Installation of edgeR package
install.packages("BiocManager")
BiocManager::install("edgeR")
library(edgeR) #after package installation, library installation is compulsory

#Create a DGEList object
dge <- DGEList(counts = combined_counts, group = group) #DON'T REPLACE FIRST g with the capital letter

#Filter lowly expressed genes
keep <- filterByExpr(dge)
dge <- dge[keep,, keep.lib.size = FALSE]

#Normalization (TMM)
dge <- calcNormFactors(dge)

#How many genes are kept
sum(keep)

#Log-CPM transformation
logCPM <- cpm(dge, log = TRUE, prior.count = 1)

#Run PCA
pca <- prcomp(t(logCPM), scale. = TRUE)

#Visualization of PCA
pca_df <- as.data.frame(pca$x)
pca_df$group <- group

#Installation of ggplot package
install.packages("ggplot2")

#Library installation of ggplot2
library(ggplot2)

#Creating a plot using ggplot
ggplot(pca_df, aes(PC1, PC2, color = group, shape = group, size = group)) + 
  geom_point() +
  
  # Axes labels with explained variance
  xlab(paste0("PC1: ", round(summary(pca)$importance[2,1]*100, 1), "%")) +
  ylab(paste0("PC2: ", round(summary(pca)$importance[2,2]*100, 1), "%")) +
  
  # Minimal theme with font/size tweaks
  theme_minimal() +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.1), #hjust is horizontal justification, left to right from 0 to 1
    text = element_text(size = 16),
    axis.title.x = element_text(size = 16, face = "plain", color="black", hjust = 0.5, margin = margin(t=6, b=5)), #t is top, b is bottom, but these features works in the context of title box, but not the whole plot;
    axis.title.y = element_text(size = 16, face = "plain", color = "black", margin = margin(r=0, l=5)), #r is right, l is left side distance;
    axis.text.x = element_text(size = 16, color = "black"),
    axis.text.y = element_text(size = 16, color = "black"),
    legend.title = element_text(size = 16, face = "plain"),
    legend.text = element_text(size = 16, margin = margin(r=25, l=0))
  ) +
  
  # Capitalize legend title
  labs(color = "Subgroups", shape = "Subgroups", size= "Subgroups") + 
  scale_shape_manual(values = c(
    "BALBc_B" = 17,  # filled circle
    "BALBc_BM" = 16,  # filled triangle
    "C57BL6J_B" = 17,  # filled circle
    "C57BL6J_BM" = 16  # filled triangle
  )) +
  
  # Customize sizes (make diamond bigger, others normal)
  scale_size_manual(values = c(
    "BALBc_B" = 5.5,
    "BALBc_BM" = 7,
    "C57BL6J_B" = 5.5,
    "C57BL6J_BM" = 7 
  )) +
  
  # Customize colors (set per group)
  scale_color_manual(values = c(
    "BALBc_B" = "red",
    "BALBc_BM" = "red",
    "C57BL6J_B" = "blue",
    "C57BL6J_BM" = "blue"   
  ))

#The end!
