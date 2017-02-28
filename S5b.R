library(googleVis)
library(dplyr)
library(reshape2)

d <- read.csv("data/S5b.csv")

reference <- d$SEURAT
clusters <- d$SC3

colors <- c('#FF0000',
            '#FFA500',
            '#008000', '#008000', '#008000', '#008000', '#008000', '#008000', '#008000', '#008000', '#008000', '#008000', '#008000', '#008000', '#008000', '#008000', '#008000', '#008000', '#008000', '#008000', '#008000', '#008000', '#008000',
            '#57A5D6', '#57A5D6',
            '#0000FF', '#0000FF', '#0000FF', '#0000FF', '#0000FF', '#0000FF', '#0000FF', '#0000FF',
            '#800080', '#800080', '#800080', '#800080', '#800080', '#800080')

res.all <- NULL
for(j in 1:39){
    res <- NULL
    for(i in 1:39) {
        tmp <- length(intersect(which(clusters == i), which(reference == j)))
        res <- c(res, tmp)
    }
    res.all <- rbind(res.all, res)
}

colnames(res.all) <- 1:39
rownames(res.all) <- 1:39

res <- melt(res.all)
res <- res[res$value != 0, ]

maxs <- res %>%
    group_by(Var1) %>%
    dplyr::summarise(max = max(value))

res <- merge(res, maxs)
maxs <- res[res$value == res$max, ]
res <- res[res$value != res$max, ]

res <- rbind(maxs, res)
res <- res[,1:3]

# remove cycles from the data
res[, 1] <- paste0(res[, 1], " ")
res[, 2] <- paste0(" ", res[, 2])

colnames(res) <- c("From", "To", "Weight")

if(!is.null(colors)) {
    colors <- paste(colors, collapse = "', '")
    colors <- paste0("['", colors, "']")
}

Sankey <- gvisSankey(
    res,
    from="From",
    to="To",
    weight="Weight",
    options = list(
        width = 700,
        height = 1000,
        sankey = paste0("{
            node:{
                label:{
                    fontName:'Arial',
                    fontSize:11,color:
                    '#000000',
                    bold:true,
                    italic:false
                },
                colors:'#FFFFFF',
                nodePadding:12
            },", if(!is.null(colors)) {
                paste0("link:{
                colorMode: 'source',
                colors: ", colors, "
            },")},
                "iterations:0
        }"
        ))
)

plot(Sankey)

