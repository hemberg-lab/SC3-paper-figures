library(cowplot)
library(dplyr)

options(stringsAsFactors = FALSE)

d <- read.csv("data/S1.csv")

d_med <- d %>%
    group_by(Dataset, gene_filter, Distance, Transformation, d) %>%
    dplyr::summarise(ARI = median(ARI))

d$Group <- paste(d$Transformation, d$d)
d$Dataset <- factor(d$Dataset, levels = c("Biase", "Yan", "Goolam"))
d_med$Dataset <- factor(d_med$Dataset, levels = c("Biase", "Yan", "Goolam"))

p <- ggplot(d, aes(d, ARI, group = Group, color = Transformation)) +
    geom_boxplot(outlier.size = 0.3, width = 1) +
    facet_grid(Dataset + Distance ~ gene_filter, scales = "free_x") +
    labs(x = "d as % of N", y = "ARI") +
    geom_vline(xintercept = c(4, 7)) +
    geom_line(data = d_med, aes(group = Transformation), size = 0.5) +
    theme_classic(base_size = 8) +
    theme(axis.line=element_blank(),
          strip.background = element_rect(colour = "white")) +
    annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, color = "black")+
    annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, color = "black")

ggsave("jpeg/S1.jpeg", w = 9, h = 6)
