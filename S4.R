library(cowplot)
library(dplyr)

options(stringsAsFactors = FALSE)
font_size <- 7

get_a <- function() {
    d <- read.csv("data/S4a.csv")
    d$group <- paste(d$nstart, d$Method)
    d$Method <- factor(
        d$Method, 
        levels = c("SC3", "tSNE+kmeans", "pcaReduce", "SNN-Cliq", "SINCERA", "SEURAT")
    )
    p <- ggplot(d, aes(x = cells, y = time, color = Method, group = group)) +
        geom_point() +
        geom_line() +
        scale_fill_brewer(palette="Set1") +
        scale_colour_brewer(palette="Set1") +
        labs(x = "Number of cells", y = "Time, min") +
        theme_classic(base_size = font_size)
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

get_b <- function() {
    set.seed(1234567)
    cols <- c("Biase" = "#bc80bd", "Treutlein" = "#8dd3c7", "Ting" = "#ffffb3", 
              "Yan" = "#ccebc5", "Goolam" = "#ffed6f", "Deng" = "#bebada",
              "Pollen1" = "#fb8072", "Pollen2" = "#fb8072",
              "Patel" = "#80b1d3", "Usoskin1" = "#fdb462", "Usoskin2" = "#fdb462",
              "Usoskin3" = "#fdb462", "Kolodziejczyk" = "#bf812d",
              "Klein" = "#b3de69", "Zeisel" = "#fccde5")
    
    d <- read.csv("data/S4b.csv")

    d$Dataset <- factor(
        d$Dataset, 
        levels = c("Biase", "Yan", "Goolam", "Deng", "Pollen1", "Pollen2", 
                   "Kolodziejczyk", "Treutlein", "Ting", "Patel", "Usoskin1", 
                   "Usoskin2", "Usoskin3", "Klein", "Zeisel")
    )
    
    d$nstart <- factor(d$nstart, levels = c("50", "1000"))
    cols.clust <- c("50" = "#999999", "1000" = "#e41a1c")

    d1 <- d %>%
        group_by(nstart, Dataset) %>%
        dplyr::summarise(Median = median(ARI))
    
    p <- ggplot(d, aes(x = 1, y = ARI, fill = nstart, group = nstart)) +
        geom_bar(data = d1, aes(y = Median), position="dodge", stat="identity", width = 0.6) +
        geom_point(position = position_jitterdodge(jitter.width = 0.6, dodge.width = 0.6), size = 0.1) +
        facet_wrap(ncol = 5, ~ Dataset) +
        scale_fill_manual(values = cols.clust) +
        geom_hline(yintercept = 0.8) +
        labs(x = "") +
        theme_classic(base_size = font_size) +
        theme(axis.ticks.x = element_blank(), axis.text.x=element_blank(),
              axis.title.x=element_blank(), axis.line=element_blank()) +
        annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, color = "black")+
        annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, color = "black")
    
    dummy <- ggplot(d, aes(x = 1, y = ARI, fill = nstart)) +
        facet_wrap(ncol = 5, ~ Dataset) +
        geom_rect(aes(fill = Dataset), xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +
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

get_c <- function() {
    cols <- c("Treutlein" = "#8dd3c7", "Ting" = "#ffffb3", "Deng" = "#bebada",
              "Pollen2" = "#fb8072", "Patel" = "#80b1d3", 
              "Kolodziejczyk" = "#bf812d", "Usoskin3" = "#fdb462",
              "Klein" = "#40E0D0", "Zeisel" = "#fccde5", "Macosko" = "#d9d9d9")
    
    d <- read.csv("data/S4c.csv")
    
    d$Dataset <- factor(
        d$Dataset,
        levels = c(
            "Deng",
            "Pollen2",
            "Kolodziejczyk",
            "Patel",
            "Usoskin3",
            "Klein",
            "Zeisel",
            "Macosko"
        )
    )
    
    d$Fraction <- factor(
        d$Fraction,
        levels = sort(unique(as.numeric(d$Fraction)))
    )
    
    p <- ggplot(d, aes(x = 1, ARI, fill = Dataset, color = Dataset)) +
        geom_boxplot(position = position_dodge(width = 1.5)) +
        geom_hline(yintercept = 0.8) +
        labs(x = "# of training cells as % of N", y = "ARI") +
        scale_fill_manual(values = cols) +
        scale_colour_manual(values = cols) +
        facet_grid(. ~ Fraction) +
        # coord_cartesian(ylim = c(0, 1)) +
        theme_classic(base_size = font_size) +
        theme(axis.ticks.x = element_blank(), axis.text.x=element_blank(),
              axis.title.x=element_blank(), axis.line=element_blank(),
              strip.background = element_rect(colour = "white"),
              legend.key.size = unit(0.4, "cm")) +
        ylim(0,1) +
        annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, color = "black")+
        annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, color = "black")
    p <- ggdraw(p) + 
        draw_label("% of total # of cells\nin a training set", 
                   fontface = "bold",
                   size = font_size-3,
                   x = 0.87, y = 0.93)
    return(p)
}

get_d <- function() {
    cols <- c("Pollen2" = "#fb8072",
              "Kolodziejczyk" = "#bf812d")
    
    d <- read.csv("data/S4d.csv")
    d$Dataset <- factor(d$Dataset, levels <- c("Pollen2", "Kolodziejczyk"))
    
    limits <- aes(ymax = mean + se, ymin = mean - se)
    
    p <- ggplot(d, aes(cells, mean, color = Dataset, group = cells)) +
        geom_point(aes(group = Dataset)) +
        geom_errorbar(limits, width=0.2) +
        geom_line(aes(group = Dataset)) +
        scale_colour_manual(values = cols) +
        labs(x = "% of N", y = "% of runs") +
        theme_classic(base_size = font_size)
    return(p)
}

get_e <- function() {
    cols <- c("Pollen2" = "#fb8072",
              "Kolodziejczyk" = "#bf812d")
    
    d <- read.csv("data/S4e.csv")
    d$Dataset <- factor(d$Dataset, levels <- c("Pollen2", "Kolodziejczyk"))
    
    limits <- aes(ymax = mean + se, ymin = mean - se)
    
    p <- ggplot(d, aes(cells, mean, color = Dataset, group = cells)) +
        geom_point(aes(group = Dataset)) +
        geom_errorbar(limits, width=0.1) +
        geom_line(aes(group = Dataset)) +
        scale_colour_manual(values = cols) +
        labs(x = "% of N", y = "% of runs") +
        theme_classic(base_size = font_size)
    return(p)
}


first_row <- plot_grid(get_a(), get_b(), ncol = 2, labels = c("a", "b"))
second_row <- plot_grid(get_c(), ncol = 1, labels = c("c"))
third_row <- plot_grid(get_d(), get_e(), ncol = 2, labels = c("d", "e"))

plot_grid(first_row, second_row, third_row, nrow = 3)

ggsave("jpeg/S4.jpeg", w = 9, h = 6)
