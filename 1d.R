library(cowplot)
library(dplyr)

options(stringsAsFactors = FALSE)

gtable_select <- function (x, ...)
{
    matches <- c(...)
    x$layout <- x$layout[matches, , drop = FALSE]
    x$grobs <- x$grobs[matches]
    x
}

gtable_stack <- function(g1, g2){
    g1$grobs <- c(g1$grobs, g2$grobs)
    g1$layout <- transform(g1$layout, z= z-max(z), name="g2")
    g1$layout <- rbind(g1$layout, g2$layout)
    g1
}

d <- read.csv("data/1d.csv")

d[d$Dataset == "Kolodziejczyk", ]$Dataset <- "Kolodz."

cols <- c("Biase" = "#bc80bd", "Treutlein" = "#8dd3c7", "Ting" = "#ffffb3", 
          "Yan" = "#ccebc5", "Goolam" = "#ffed6f", "Deng" = "#bebada",
          "Pollen1" = "#fb8072", "Pollen1-TopHat" = "#fb8072",
          "Pollen2-TopHat" = "#fb8072", "Pollen2" = "#fb8072",
          "Patel" = "#80b1d3", "Usoskin1" = "#fdb462", "Usoskin2" = "#fdb462",
          "Usoskin3" = "#fdb462", "Petropoulos" = "#8dd3c7", "Kolodz." = "#bf812d",
          "Klein" = "#b3de69", "Zeisel" = "#fccde5")
d$Dataset <- factor(
    d$Dataset, 
    levels = c(
        "Biase", "Yan", "Goolam", "Deng", "Pollen1", "Pollen2", "Kolodz.",
        "Treutlein", "Ting", "Patel", "Usoskin1", "Usoskin2", "Usoskin3",
        "Klein", "Zeisel"
    )
)

colnames(d)[2] <- "Clustering"
d$Clustering <- factor(d$Clustering, levels = c("Individual", "Consensus"))

cols.clust <- c("Individual" = "#999999", "Consensus" = "#e41a1c")

d1 <- d %>%
    group_by(Clustering, Dataset) %>%
    dplyr::summarise(Median = median(ARI))

p <- ggplot(d, aes(x = 1, y = ARI, fill = Clustering, group = Clustering)) +
    geom_bar(data = d1, aes(y = Median), position="dodge", stat="identity", width = 0.1) +
    geom_point(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.1), size = 0.4) +
    facet_wrap(ncol = 5, ~ Dataset) +
    scale_fill_manual(values = cols.clust) +
    geom_hline(yintercept = 0.8) +
    labs(x = "") +
    theme_classic(base_size = 14) +
    theme(axis.ticks.x = element_blank(), axis.text.x=element_blank(),
          axis.title.x=element_blank(), axis.line=element_blank()) +
    annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, color = "black")+
    annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, color = "black")

dummy <- ggplot(d, aes(x = 1, y = ARI, fill = Clustering)) +
    facet_wrap(ncol = 5, ~ Dataset) +
    geom_rect(aes(fill = Dataset), xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +
    scale_fill_manual(values = cols) +
    theme_minimal()

g1 <- ggplotGrob(p)
g2 <- ggplotGrob(dummy)

panels <- grepl(pattern="panel", g2$layout$name)
strips <- grepl(pattern="strip-t", g2$layout$name)
g2$layout$t[panels] <- g2$layout$t[panels] - 1
g2$layout$b[panels] <- g2$layout$b[panels] - 1

new_strips <- gtable_select(g2, panels | strips)

## ideally you'd remove the old strips, for now they're just covered
new_plot <- gtable_stack(g1, new_strips)

plot_grid(new_plot)

ggsave("pdf/1d.pdf", w = 6.3, h = 4.5)
ggsave("jpeg/1d.jpeg", w = 6.3, h = 4.5)
