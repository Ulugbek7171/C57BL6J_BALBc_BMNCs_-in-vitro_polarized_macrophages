# ============================================================
# Figure 1c — Flow cytometry quantification
# C57BL/6J vs BALB/c
# n = 3 biological replicates per strain
# ============================================================

# Load package
library(tidyverse)

# ------------------------------------------------------------
# 1. Import data
# ------------------------------------------------------------

data <- read.csv(
  "/Users/surgeon7171/Desktop/Book1.csv",
  stringsAsFactors = FALSE
)

# Check imported data
print(data)


# ------------------------------------------------------------
# 2. Convert from wide to long format
# ------------------------------------------------------------

data_long <- data %>%
  pivot_longer(
    cols = c(CD11b, F4_80, CD206, F4_80_CD206),
    names_to = "Marker",
    values_to = "Percentage"
  )


# ------------------------------------------------------------
# 3. Clean strain names
# ------------------------------------------------------------

data_long <- data_long %>%
  mutate(
    Strain = case_when(
      Strain == "C57BL6J" ~ "C57BL/6J",
      Strain == "BALBc"   ~ "BALB/c",
      TRUE ~ Strain
    )
  )


# ------------------------------------------------------------
# 4. Set order of markers and strains
# ------------------------------------------------------------

data_long$Strain <- factor(
  data_long$Strain,
  levels = c("C57BL/6J", "BALB/c")
)

data_long$Marker <- factor(
  data_long$Marker,
  levels = c(
    "CD11b",
    "F4_80",
    "CD206",
    "F4_80_CD206"
  ),
  labels = c(
    "CD11b",
    "F4/80",
    "CD206",
    "F4/80+CD206+"
  )
)


# ------------------------------------------------------------
# 5. Check that there are 3 replicates per group
# ------------------------------------------------------------

print(
  data_long %>%
    count(Marker, Strain)
)


# ------------------------------------------------------------
# 6. Calculate mean and SD
# ------------------------------------------------------------

summary_data <- data_long %>%
  group_by(Marker, Strain) %>%
  summarise(
    Mean = mean(Percentage, na.rm = TRUE),
    SD = sd(Percentage, na.rm = TRUE),
    n = sum(!is.na(Percentage)),
    .groups = "drop"
  )

print(summary_data)


# ------------------------------------------------------------
# 7. Calculate Welch's two-sided t-test
# ------------------------------------------------------------

p_values <- data_long %>%
  group_by(Marker) %>%
  summarise(
    p_value = t.test(
      Percentage ~ Strain,
      data = cur_data(),
      var.equal = FALSE,
      alternative = "two.sided"
    )$p.value,
    .groups = "drop"
  )

print(p_values)


# ------------------------------------------------------------
# 8. Format P values
# ------------------------------------------------------------

p_values <- p_values %>%
  mutate(
    p_label = case_when(
      p_value < 0.0001 ~ "P < 0.0001",
      p_value < 0.001  ~ "P < 0.001",
      p_value < 0.01   ~ "P < 0.01",
      p_value < 0.05   ~ "P < 0.05",
      TRUE ~ paste0("P = ", sprintf("%.3f", p_value))
    )
  )

print(p_values)


# ------------------------------------------------------------
# 9. Determine position of P values
# ------------------------------------------------------------

p_positions <- summary_data %>%
  group_by(Marker) %>%
  summarise(
    y.position = max(Mean + SD) + 3,
    .groups = "drop"
  ) %>%
  left_join(
    p_values,
    by = "Marker"
  )


# ------------------------------------------------------------
# 10. Create Figure 1c
# ------------------------------------------------------------

figure_1c <- ggplot(
  summary_data,
  aes(
    x = Strain,
    y = Mean,
    fill = Strain
  )
) +
  
  # Mean bars
  geom_col(
    width = 0.65,
    color = "black"
  ) +
  
  # SD error bars
  geom_errorbar(
    aes(
      ymin = Mean - SD,
      ymax = Mean + SD
    ),
    width = 0.15,
    linewidth = 0.6
  ) +
  
  # Individual biological replicates
  geom_jitter(
    data = data_long,
    aes(
      x = Strain,
      y = Percentage
    ),
    width = 0.08,
    size = 2.5,
    shape = 21,
    fill = "white",
    color = "black"
  ) +
  
  # P values
  geom_text(
    data = p_positions,
    aes(
      x = 1.5,
      y = y.position,
      label = p_label
    ),
    inherit.aes = FALSE,
    size = 3.5
  ) +
  
  # Four panels
  facet_wrap(
    ~ Marker,
    nrow = 1
  ) +
  scale_fill_manual(
    values = c(
      "C57BL/6J" = "blue",
      "BALB/c" = "red"
    )
  ) +
  # Axis labels
  labs(
    x = NULL,
    y = "Marker-positive cells (%)"
  ) +
  
  # Clean publication-style appearance
  theme_classic(
    base_size = 12
  ) +
  
  theme(
    legend.position = "none",
    
    strip.background = element_blank(),
    
    strip.text = element_text(
      face = "bold",
      size = 12
    ),
    
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      color = "black"
    ),
    
    axis.text.y = element_text(
      color = "black"
    ),
    
    axis.title.y = element_text(
      size = 12
    )
  )


# ------------------------------------------------------------
# 11. Display figure
# ------------------------------------------------------------

print(figure_1c)


# ------------------------------------------------------------
# 12. Save figure
# ------------------------------------------------------------

ggsave(
  filename = "/Users/surgeon7171/Desktop/Figure_1c.pdf",
  plot = figure_1c,
  width = 8,
  height = 3.8,
  units = "in"
)

ggsave(
  filename = "/Users/surgeon7171/Desktop/Figure_1c.png",
  plot = figure_1c,
  width = 9,
  height = 3,
  units = "in",
  dpi = 600
)
