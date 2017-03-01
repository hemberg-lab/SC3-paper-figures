library(cowplot)

d <- read.csv("data/S12a.csv")

d$Genotype <- factor(d$Genotype, levels = c("WT", "Tet2"))

p <- ggplot(d, aes(Genotype, CoV)) +
    geom_boxplot(width = 0.3) +
    theme_classic(base_size = 9) +
    labs(x = "", y = "Coefficient of variation")

plot_grid(p, NULL, ncol = 2, labels = c("a", "b"))

ggsave("pdf/S12.pdf", w = 9, h = 6)
ggsave("jpeg/S12.jpeg", w = 9, h = 6)
