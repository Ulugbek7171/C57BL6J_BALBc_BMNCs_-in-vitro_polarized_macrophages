#Creating bar graph highlighting read numbers
#Install and load required packages
install.packages(c("readxl", "ggplot2", "dplyr")) 
library(readxl) #reads excel
library(ggplot2) #creates plots
library(dplyr) #allows to manipulate data

#Check the WD
getwd()
#Read the data
data <- read_excel("Results_all.xlsx", sheet = "Metadata_RNA-seq") #reading excel file, sheet = reads the defined worksheet in excel file;
View(data) #view excel file
colnames(data) #view column names

#Extract columns from loaded excel file
df <- data[1:12, c("Sample ID", "Read Number")] #data extraction from the uploaded excel file; 1:30 means only rows from 1 to 30; 
df$Group <- rep(1:4, each = 3)  # 12 samples → 4 groups, 3 samples in each
View(df) #View extracted data


# Create group-of-three color variable
df <- df %>%
  mutate(Condition = ifelse(Group %in% 1:2, "C57BL/6J", "BALB/c"))  # 12 samples → 4 groups, 3 samples in each
#Define custom colors for each subgroup 
group_colors <- c("blue", "blue", "red",
                   "red")

# Force Sample ID order to match Excel order
df$`Sample ID` <- factor(df$`Sample ID`, levels = df$`Sample ID`)

# Make sure Group is a factor with levels in correct order
df$Group <- factor(df$Group, levels = unique(df$Group))

#Plot creation
ggplot(df, aes(x = `Sample ID`, y = `Read Number`, fill = factor(Group))) +
  geom_bar(stat = "identity", width = 0.8) +  # bar width customization
  
  #Subgroup colors 
  scale_fill_manual(
    values = group_colors,
    name = "Subgroups", # legend title
    labels = c("BMNC_M", "BMNC",
               "BMNC_M", "BMNC") #Subgroup titles matching the sample order
  ) +
  # Add Condition legend using dummy points
  geom_point(aes(y = 0, color = Condition), size = 5, alpha = 0) +
  scale_color_manual(
  name = "Mouse Strain",
  values = c("C57BL/6J" = "blue", "BALB/c" = "red")
  ) +
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 5))) +
  
  theme_minimal() +
  labs(title = "", # plot title
       x = "Sample ID", # x axis title
       y = "Read Number") + # y axis title
  theme( # for customization of plot
    axis.text.x = element_text(angle = 90, hjust = 1, size = 16, face = "plain", color = "black"), # x axis label customization
    axis.text.y = element_text(size = 16, face = "plain", color = "black"), # y axis label customization
    axis.title.x = element_text(size = 16, face = "plain", vjust = 0.1), # x axis title customization, vjust = vertical justification;
    axis.title.y = element_text(size = 16, face = "plain"), # y axis title customization
    plot.title   = element_text(size = 16, face = "plain", hjust = 0.5), # plot title customization, hjust = horizontal justification;
    legend.position = "right", # legend position in the plot
    legend.box = "vertical", # legend box position itself
    legend.title = element_text(size = 16, face = "plain"), # legend title customization, face = bold, italic, or bold.italic;
    legend.text = element_text(size = 16, face = "plain"), # legend text customization
    plot.margin = margin(t = 10, r = 10, b = 30, l = 10) # plot margins customization
  )

# The end
