library(cowplot)

options(stringsAsFactors = FALSE)
font_size <- 9

d1 <- read.csv("data/S7a1.csv")
d2 <- read.csv("data/S7a2.csv")
ercc.fraction1 <- read.csv("data/S7b1.csv")
ercc.fraction2 <- read.csv("data/S7b2.csv")

p1 <- ggplot(d1, aes(x)) +
    geom_histogram(bins = 130) + 
    geom_vline(xintercept = 420, color = "red") +
    labs(x = "# of expressed genes", y = "# of cells", title = "Patient 1") +
    theme_classic(base_size = font_size)

p2 <- ggplot(d2, aes(x)) +
    geom_histogram(bins = 130) + 
    geom_vline(xintercept = 900, color = "red") +
    labs(x = "# of expressed genes", y = "# of cells", title = "Patient 2") +
    theme_classic(base_size = font_size)

p3 <- ggplot(ercc.fraction1, aes(x)) +
    geom_histogram(bins = 100) + 
    labs(x = "(# of ERCC reads)/(# endogenous reads)", y = "# of cells") +
    xlim(0, 0.15) +
    theme_classic(base_size = font_size)
    
p4 <- ggplot(ercc.fraction2, aes(x)) +
    geom_histogram(bins = 100) + 
    labs(x = "(# of ERCC reads)/(# endogenous reads)", y = "# of cells") +
    xlim(0, 0.15) +
    geom_vline(xintercept = 0.05, color = "red") +
    theme_classic(base_size = font_size)

plot_grid(p1, p2, p3, p4, ncol = 2, labels = c("a", "", "b", ""))

ggsave("jpeg/S8.jpeg", w = 9, h = 6)
