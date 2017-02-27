library(cowplot)
library(dplyr)

options(stringsAsFactors = FALSE)

# Download raw data from https://s3.eu-west-2.amazonaws.com/sc3-paper-figures/fig1-raw.csv
# then:
d <- read.csv("fig1-raw.csv")
d <- d[d$Dataset %in% c("Biase", "Yan", "Goolam", "Deng", "Pollen", "Kolodziejczyk"), ]
d <- d[d$gene_filter == "With Gene Filter", ]
d <- d[, c(1:4, 6)]

d_med <- d %>%
    group_by(Dataset, Distance, Transformation, d) %>%
    dplyr::summarise(ARI = median(ARI))

tmp <- d_med %>%
    group_by(Dataset, Distance, Transformation) %>%
    dplyr::summarise(peak = max(ARI))

res <- merge(d_med, tmp)
res <- res[res$ARI >= 0.95*res$peak, ]
res <- res[,c(1,4)]

# save data for online publication
write.csv(res, "1c.csv", quote = FALSE, row.names = FALSE)

res$Dataset <- factor(
    res$Dataset, 
    levels = c("Biase", "Yan", "Goolam", "Deng", "Pollen", "Kolodziejczyk")
)

cols <- c("Biase" = "#bc80bd", 
          "Yan" = "#ccebc5", 
          "Goolam" = "#ffed6f", 
          "Deng" = "#bebada", 
          "Pollen" = "#fb8072", 
          "Kolodziejczyk" = "#bf812d")

p <- ggplot(res, aes(d, fill = Dataset)) +
    geom_histogram(bins = 22) +
    scale_fill_manual(values = cols) +
    geom_vline(xintercept = c(3.5,7.5)) +
    labs(x = "d as % of N", y = "# of solutions with ARI > 95% of max.") +
    theme_classic(base_size = 14)

ggsave("1c.pdf", w = 6, h = 4.5)
