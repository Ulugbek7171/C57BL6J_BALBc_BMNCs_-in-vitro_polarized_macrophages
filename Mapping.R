library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)

#Check the WD
getwd()
setwd("/Users/surgeon7171/Library/CloudStorage/OneDrive-OUMail(OsakaUniversity)/Research/TAC/RNA-seq_analyses/R projects/Scientific_Data_Read_Numbers/Data")
#Read the data
data <- read_excel("F10115_F9808_Mapping_Statistics.xlsx", sheet = "F10115_statistics") #reading excel file, sheet = reads the defined worksheet in excel file;
View(data) #view excel file
colnames(data) #view column names

#Extract columns from loaded excel file
df <- data[1:12, c("LIB_ID", "BAM_file", "Total_reads", "Mapped_reads", "Unmapped_reads", "Primary_mapped_reads")] #data extraction from the uploaded excel file; 1:30 means only rows from 1 to 30; 
df$Group <- rep(1:4, each = 3)  # 12 samples → 4 groups, 3 samples in each
View(df) #View extracted data
#Step 1: Prepare data for stacked-percent bars
df_long <- df %>%
  mutate(
    secondary_mapped_reads = Mapped_reads - Primary_mapped_reads,
    LineType = ifelse(row_number() <= 6, "BALB/c", "C57BL/6J")
  ) %>%
  select(
    BAM_file,
    Total_reads,
    Primary_mapped_reads,
    secondary_mapped_reads,
    Unmapped_reads,
    LineType
  ) %>%
  pivot_longer(
    cols = c(Primary_mapped_reads,
             secondary_mapped_reads,
             Unmapped_reads),
    names_to = "ReadType",
    values_to = "Count"
  ) %>%
  group_by(BAM_file) %>%
  mutate(Percent = Count / sum(Count))
#Step 2: Define factor levels (controls order & legends)
df_long$ReadType <- factor(
  df_long$ReadType,
  levels = c("Primary_mapped_reads",
             "secondary_mapped_reads",
             "Unmapped_reads"),
  labels = c("Primary mapped", "Secondary mapped", "Unmapped")
)

df_long$LineType <- factor(
  df_long$LineType,
  levels = c("BALB/c", "C57BL/6J")
)
#Step 3: Define colors & linetypes
read_colors <- c(
  "Primary mapped"   = "salmon2",
  "Secondary mapped" = "orange",
  "Unmapped"         = "purple"
)

line_types <- c(
  "BALB/c"    = "dashed",
  "C57BL/6J"  = "solid"
)
#Step 4: Make the plot (this meets ALL 4 requirements)
ggplot(df_long,
       aes(x = BAM_file,
           y = Percent,
           fill = ReadType,
           linetype = LineType)) +
  geom_col(
    color = "black",
    linewidth = 0.6
  ) +
  scale_fill_manual(
    values = read_colors,
    name = "Read category"
  ) +
  scale_linetype_manual(
    values = line_types,
    name = "Strain"
  ) +
  scale_y_continuous(
    labels = percent_format()
  ) +
  labs(
    title = "",
    x = "Sample ID",
    y = "Percentage of total reads"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    text = element_text(color = "black"),
    axis.text.x = element_text(angle = 90, hjust = 1, color = "black"),
    axis.text.y = element_text(color = "black"),
    legend.position = "left",
    legend.box = "vertical",
    legend.key = element_rect(fill = NA, colour = NA)   # ← remove Read category borders
  ) +
  guides(
    linetype = guide_legend(
      override.aes = list(
        fill = "grey",    # ← Strain legend box color
        colour = "black"
      )
    )
  )
