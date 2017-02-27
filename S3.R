library(cowplot)
library(dplyr)

options(stringsAsFactors = FALSE)
font_size <- 6

get_a <- function() {
    d <- read.csv("data/S3a.csv")
    
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
        geom_histogram(bins = 21) +
        scale_fill_manual(values = cols) +
        geom_vline(xintercept = c(3.4,7.7)) +
        labs(x = "d as % of N", y = "# of runs with ARI > 95% of max") +
        theme_classic(base_size = font_size)
    return(p)
}

get_b <- function() {
    d <- read.csv("data/S3b.csv")
    
    d$Dataset <- factor(
        d$Dataset, 
        levels = c("Treutlein", "Ting", "Patel", "Usoskin", "Klein", "Zeisel")
    )
    
    cols <- c("Treutlein" = "#8dd3c7", 
              "Ting" = "#ffffb3", 
              "Patel" = "#80b1d3", 
              "Usoskin" = "#fdb462", 
              "Klein" = "#b3de69", 
              "Zeisel" = "#fccde5")
    
    p <- ggplot(d, aes(d, fill = Dataset)) +
        geom_histogram(bins = 24) +
        scale_fill_manual(values = cols) +
        geom_vline(xintercept = c(3.5,7.6)) +
        labs(x = "d as % of N", y = "# of runs with ARI > 95% of max") +
        theme_classic(base_size = font_size)
    return(p)
}

get_c <- function() {
    d <- read.csv("data/S3c.csv")
    
    d$Dataset <- factor(
        d$Dataset, 
        levels = c("Biase", "Yan", "Goolam", "Deng", "Pollen1", "Pollen2", "Kolodziejczyk")
    )
    
    d_med <- d %>%
        group_by(Dataset, x, Param) %>%
        dplyr::summarise(Median = median(ARI))
    
    p1 <- ggplot(d[d$Param == "Gene Filter Fraction",], 
                  aes(x, ARI, fill = Param)) +
        geom_bar(data = d_med[d_med$Param == "Gene Filter Fraction",], 
                 aes(y = Median), stat = "identity", position = "dodge") +
        facet_grid(Param ~ Dataset, scales = "free") +
        geom_point(position = 
                       position_jitterdodge(jitter.width = 0.7, dodge.width = 0.9), 
                   size = 0.1) +
        theme_classic(base_size = font_size) +
        theme(axis.line=element_blank(),
              strip.background = element_rect(colour = "white")) +
        annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, color = "black")+
        annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, color = "black")+
        guides(fill = FALSE) +
        labs(x = "% of N") +
        scale_fill_manual(values = c("#e41a1c"))
    
    p2 <- ggplot(d[d$Param != "Gene Filter Fraction",], aes(x, ARI, fill = "red")) +
        geom_bar(data = d_med[d_med$Param != "Gene Filter Fraction",], 
                 aes(y = Median), stat = "identity", position = "dodge") +
        facet_grid(Param ~ Dataset, scales = "free") +
        geom_point(position = 
                       position_jitterdodge(jitter.width = 0.7, dodge.width = 0.9), 
                   size = 0.1) +
        theme_classic(base_size = font_size) +
        theme(axis.line=element_blank(),
              strip.background = element_rect(colour = "white")) +
        annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, color = "black")+
        annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, color = "black")+
        labs(x = "Min. expression value") +
        guides(fill = FALSE) +
        scale_fill_manual(values = c("#e41a1c"))

    p <- plot_grid(p1, p2, ncol = 1, align = 'v')
    return(p)
}

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

get_d <- function() {
    set.seed(1234567)
    cols <- c("Biase" = "#bc80bd", "Treutlein" = "#8dd3c7", "Ting" = "#ffffb3", 
              "Yan" = "#ccebc5", "Goolam" = "#ffed6f", "Deng" = "#bebada",
              "Pollen1" = "#fb8072", "Pollen2" = "#fb8072",
              "Patel" = "#80b1d3", "Usoskin1" = "#fdb462", "Usoskin2" = "#fdb462",
              "Usoskin3" = "#fdb462", "Kolodziejczyk" = "#bf812d",
              "Klein" = "#b3de69", "Zeisel" = "#fccde5")
    
    d <- read.csv("data/S3d.csv")
    
    d$Hierarchy <- factor(
        d$Hierarchy, 
        levels = c("Biase", "Yan", "Goolam", "Deng", "Pollen1", "Pollen2", 
                   "Kolodziejczyk", "Treutlein", "Ting", "Patel", "Usoskin1", 
                   "Usoskin2", "Usoskin3", "Klein", "Zeisel")
    )
    
    d$Dropouts <- factor(d$Dropouts, levels = c("Without", "With"))
    
    cols.clust <- c("Without" = "#999999", "With" = "#e41a1c")
    
    d1 <- d %>%
        group_by(Dropouts, Hierarchy) %>%
        dplyr::summarise(Median = median(ARI))
    
    p <- ggplot(d, aes(x = 1, y = ARI, fill = Dropouts, group = Dropouts)) +
        geom_bar(data = d1, aes(y = Median), position="dodge", stat="identity", width = 0.6) +
        geom_point(position = position_jitterdodge(jitter.width = 0.6, dodge.width = 0.6), size = 0.1) +
        facet_wrap(ncol = 5, ~ Hierarchy) +
        scale_fill_manual(values = cols.clust) +
        geom_hline(yintercept = 0.8) +
        labs(x = "") +
        theme_classic(base_size = font_size) +
        theme(axis.ticks.x = element_blank(), axis.text.x=element_blank(),
              axis.title.x=element_blank(), axis.line=element_blank()) +
        annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, color = "black")+
        annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, color = "black")
    
    dummy <- ggplot(d, aes(x = 1, y = ARI, fill = Dropouts)) +
        facet_wrap(ncol = 5, ~ Hierarchy) +
        geom_rect(aes(fill = Hierarchy), xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +
        scale_fill_manual(values = cols) +
        theme_minimal(base_size = font_size)
    
    g1 <- ggplotGrob(p)
    g2 <- ggplotGrob(dummy)
    
    panels <- grepl(pattern="panel", g2$layout$name)
    strips <- grepl(pattern="strip-t", g2$layout$name)
    g2$layout$t[panels] <- g2$layout$t[panels] - 1
    g2$layout$b[panels] <- g2$layout$b[panels] - 1
    
    new_strips <- gtable_select(g2, panels | strips)
    
    ## ideally you'd remove the old strips, for now they're just covered
    new_plot <- gtable_stack(g1, new_strips)
    return(new_plot)
}

plot_grid(get_a(), get_b(), get_c(), get_d(), ncol = 2, labels = c("a", "b", "c", "d"), rel_heights = c(1, 1.8))

ggsave("jpeg/S3.jpeg", w = 9, h = 6)
