library(cowplot)
library(dplyr)

options(stringsAsFactors = FALSE)

d <- read.csv("data/1c.csv")

d$Dataset <- factor(
    d$Dataset, 
    levels = c("Biase", "Yan", "Goolam", "Deng", "Pollen", "Kolodziejczyk")
)

cols <- c("Biase" = "#bc80bd", 
          "Yan" = "#ccebc5", 
          "Goolam" = "#ffed6f", 
          "Deng" = "#bebada", 
          "Pollen" = "#fb8072", 
          "Kolodziejczyk" = "#bf812d")

p <- ggplot(d, aes(d, fill = Dataset)) +
    geom_histogram(bins = 22) +
    scale_fill_manual(values = cols) +
    geom_vline(xintercept = c(3.5,7.5)) +
    labs(x = "d as % of N", y = "# of solutions with ARI > 95% of max.") +
    theme_classic(base_size = 14)

ggsave("pdf/1c.pdf", w = 6, h = 4.5)
ggsave("jpeg/1c.jpeg", w = 6, h = 4.5)
