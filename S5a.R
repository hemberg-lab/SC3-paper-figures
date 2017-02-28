library(cowplot)

d <- read.csv("data/S5a.csv")

p <- ggplot(d, aes(variable, value)) +
    geom_boxplot() +
    labs(x = "", y = "ARI") +
    ylim(0, 1) +
    theme_classic(base_size = 8)

ggsave("pdf/S5a.pdf", w = 4, h = 6)
ggsave("jpeg/S5a.jpeg", w = 4, h = 6)
