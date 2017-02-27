library(ggplot2)
library(dplyr)
library(RCurl)

options(stringsAsFactors = FALSE)

# Download raw data from https://s3.eu-west-2.amazonaws.com/sc3-paper-figures/fig1-raw.csv
# then:
d <- read.csv("fig1-raw.csv")
d <- d[, c(1:2,8:9,11,13)]
d <- d[d$Dataset %in% c("Treutlein", "Ting", "Patel", "Usoskin", "Klein", "Zeisel"), ]
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
write.csv(res, "S3b.csv", quote = FALSE, row.names = FALSE)

res$Dataset <- factor(
    res$Dataset, 
    levels = c("Treutlein", "Ting", "Patel", "Usoskin", "Klein", "Zeisel")
)

cols <- c("Treutlein" = "#8dd3c7", 
          "Ting" = "#ffffb3", 
          "Patel" = "#80b1d3", 
          "Usoskin" = "#fdb462", 
          "Klein" = "#b3de69", 
          "Zeisel" = "#fccde5")

p <- ggplot(res, aes(d, fill = Dataset)) +
    geom_histogram(bins = 24) +
    scale_fill_manual(values = cols) +
    geom_vline(xintercept = c(3.5,7.6)) +
    labs(x = "d as % of N", y = "# of solutions with ARI > 95% of max.") +
    theme_classic(base_size = 14)

ggsave("S3b.pdf", w = 6, h = 4.5)
