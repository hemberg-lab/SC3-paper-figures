library(pheatmap)
library(cowplot)

options(stringsAsFactors = FALSE)
font_size <- 8

d <- readRDS("data/3.rds")

dat <- d$exprs
write.csv(dat[, d$hc$order], file = "data/3.csv", quote = FALSE, row.names = FALSE)

anno_colors <- list(Cluster = c("#e41a1c", "#984ea3", "#40E0D0", "#FFFF33"))
names(anno_colors$Cluster) <- levels(d$ann$Cluster)

p <- pheatmap(d$exprs,
         cluster_cols = d$hc,
         cluster_rows = F,
         cutree_cols = 3,
         gaps_row = c(10, 20),
         annotation_col = d$ann,
         annotation_colors = anno_colors,
         show_colnames = F,
         silent = T)

plot_grid(p$gtable, nrow = 1)

ggsave("jpeg/3.jpeg", w = 9, h = 6)
ggsave("pdf/3.pdf", w = 9, h = 6)
