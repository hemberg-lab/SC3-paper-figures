library(cowplot)
library(dplyr)
library(pheatmap)

options(stringsAsFactors = FALSE)
font_size <- 7

get_a <- function() {
    de.genes = rbind(
        c(1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
        c(0,0,0,0,0,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0),
        c(0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,0,0,0,0,0),
        c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1),
        c(1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0),
        c(1,1,1,1,1,0,0,0,0,0,1,1,1,1,1,0,0,0,0,0),
        c(1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1),
        c(0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0),
        c(0,0,0,0,0,1,1,1,1,1,0,0,0,0,0,1,1,1,1,1),
        c(0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1),
        c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0),
        c(1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,1,1,1,1,1),
        c(1,1,1,1,1,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1),
        c(0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)
    )
    
    colnames(de.genes) = paste("Cell ", 1:20, sep = "")
    rownames(de.genes) = paste("Gene ", 1:14, sep = "")
    
    p <- pheatmap(de.genes,
             cluster_row = FALSE,
             cluster_cols = FALSE,
             gaps_col = c(5, 10, 15),
             gaps_row = 4,
             legend = F,
             fontsize = font_size,
             silent = TRUE)
    return(p$gtable)
}

get_b <- function() {
    cols <- c("Biase" = "#bc80bd", "Treutlein" = "#8dd3c7", "Ting" = "#ffffb3", 
              "Yan" = "#ccebc5", "Goolam" = "#ffed6f", "Deng" = "#bebada",
              "Pollen1" = "#fb8072", "Pollen1-TopHat" = "#fb8072",
              "Pollen2-TopHat" = "#fb8072", "Pollen2" = "#fb8072",
              "Patel" = "#80b1d3", "Usoskin1" = "#fdb462", "Usoskin2" = "#fdb462",
              "Usoskin3" = "#fdb462", "Kolodziejczyk" = "#bf812d",
              "Klein" = "#b3de69", "Zeisel" = "#fccde5", "Macosko" = "#d9d9d9")
    
    d <- read.csv("data/S6c.csv")
    
    d$Dataset <- factor(d$Dataset, levels = c("Biase", "Yan", "Goolam", "Deng",
                                              "Pollen1", "Pollen2", "Kolodziejczyk",
                                              "Treutlein", "Ting", "Patel",
                                              "Usoskin1", "Usoskin2", "Usoskin3",
                                              "Klein", "Zeisel", "Macosko"))
    
    p <- ggplot(d, aes(x = AUROC, colour = Dataset, fill = Dataset)) +
        geom_density(alpha = 0.6) +
        scale_fill_manual(values = cols) +
        scale_colour_manual(values = cols) +
        labs(y = "Density") +
        theme_classic(base_size = font_size) +
        theme(legend.key.size = unit(0.4, "cm"))
    
    return(p)
}

get_c <- function() {
    d <- read.csv("data/S6b.csv")
    d$Stage <- factor(d$Stage, levels = unique(d$Stage))
    cols <- c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", 
              "#FF00FF", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A")
    p <- ggplot(d, aes(x = Cells, y = outl, fill = Stage, color = Stage)) +
        geom_bar(stat = "identity") +
        scale_fill_manual(values = cols) +
        scale_color_manual(values = cols) +
        labs(y = "Outlier score") +
        theme_classic(base_size = font_size) +
        theme(legend.key.size = unit(0.4, "cm"))
    return(p)
}

first_row <- plot_grid(NULL, get_a(), NULL, get_b(), ncol = 4, labels = c("a", "", "", "b"), rel_widths = c(0.2,0.6,0.2,1))
second_row <- plot_grid(get_c(), ncol = 1, labels = c("c"))

plot_grid(first_row, second_row, nrow = 2, rel_heights = c(1.5, 1))

ggsave("jpeg/S6.jpeg", w = 9, h = 6)
ggsave("pdf/S6.pdf", w = 9, h = 6)

