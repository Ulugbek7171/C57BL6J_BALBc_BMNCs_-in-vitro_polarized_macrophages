library(ggplot2)

rin <- c(9.7, 9.8, 9.7, 9.7, 9.6, 8.3,
         9.2, 9.5, 9.6, 8.7, 9.7, 9.9)

# Create grouping variable
group_color <- c(rep("Group1", 6), rep("Group2", 6))

df <- data.frame(
  Group = "Samples",
  RIN = rin,
  ColorGroup = group_color
)

p <- ggplot(df, aes(x = Group, y = RIN)) +
  
  # Violin
  geom_violin(fill = "skyblue", color = NA, width = 0.6, alpha = 0.6) +
  
  # Boxplot
  geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white") +
  
  # Jitter points with color mapping
  geom_jitter(aes(color = ColorGroup),
              width = 0.08,
              size = 4,
              alpha = 0.9,
              shape = 16) +
  
  # Manual colors
  scale_color_manual(values = c("Group1" = "blue",
                                "Group2" = "red")) +
  
  # Mean point
  stat_summary(fun = mean,
               geom = "point",
               shape = 23,
               size = 5,
               fill = "orange") +
  
  labs(title = "RNA Integrity Number (RIN) Distribution",
       y = "RIN Score",
       x = NULL,
       color = NULL) +
  
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

p
