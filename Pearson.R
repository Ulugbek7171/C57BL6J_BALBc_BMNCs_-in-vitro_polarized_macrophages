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
avg1 <- rowMeans(counts1)
avg2 <- rowMeans(counts2)

#Shrink space around x/y labels:
par(mgp = c(2, 1, 0))
par(mar = c(5, 5, 4, 2))
par(family = "Helvetica") #This is default font family in MacOS

#Installing packages and loading libraries
install.packages("scales")
library(scales)

#Color genes by dominant strain
point_colors <- ifelse(
  avg1 > avg2,
  alpha("red", 0.1),   # Higher in BALB/c
  alpha("blue", 0.1)     # Higher in C57BL/6J
)
col = point_colors

#Creating a plot
plot(log2(avg1 + 1), 
     log2(avg2 + 1),
     xlab = "", ylab = "",  
     main = "",
     cex.main = 1.6,
     font.main = 1, 
     cex.axis = 1.6,
     cex = 1, pch = 16, #controls the size of text and symbols; pch = means symbols, i.e. pch = 1 is open circle, pch = 16 is filled circle, pch = 17 is filled triangle;
     col = point_colors) #alpha value means transparency, i.e. 0.3 means 30%

#Adding a legend
legend("topleft",
       legend = c("C57BL/6J", "BALB/c"),
       col = c("blue", "red"),
       pch = 16,
       bty = "n",
       cex = 1.6)

#For changing x and y axes titles
mtext("Log2 avg expression (BALB/c, n=6)", 
      side = 1, line = 3, cex = 1.6, font = 1) #manual customization of x axis, line = places the text n lines away from the plot edge;

mtext("Log2 avg expression (C57BL/6J, n=6)", 
      side = 2, line = 3, cex = 1.6, font = 1) #manual customization of y axis, side = means bottom (1), left (2), top (3), right (4);

#Calculate correlation
r <- cor(log2(avg1 +1), 
         log2(avg2 + 1), 
         method = "pearson")

#Add the result to the plot
text(x = 16, y = 4, 
     labels = paste0("r = ", 
                     round(r, 2)), pos = 1, 
     cex = 1.6, font = 1, col = "black") #col = color, font = 2 - bold, cex = 
text(
  x = 16, y = 4,
  labels = bquote(italic(r) == .(round(r, 2))),
  pos = 1,
  cex = 1.6,
  col = "black"
)

#The end!