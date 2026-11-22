rin <- c(9.7, 9.8, 9.7, 9.7, 9.6, 8.3, 9.2, 9.5, 9.6, 8.7, 9.7, 9.9)
library(ggplot2)

df <- data.frame(RIN = rin)

ggplot(df, aes(x = RIN)) +
  geom_histogram(aes(y = after_stat(density)),
                 bins = 6,
                 fill = "orange",
                 color = "black") +
  geom_density(bw = 0.15,
               color = "blue",
               linewidth = 1) +
  geom_vline(aes(xintercept = mean(RIN)),
             color = "red",
             linetype = "dashed",
             linewidth = 1) +
  geom_rug() +
  labs(title = "Density Plot of RIN Values",
       x = "RIN",
       y = "Density") +
  theme_minimal()
