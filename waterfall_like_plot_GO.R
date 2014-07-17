library(ggplot2)

df <- read.table("/Volumes/xi.embl.de/g/bork8/mkuhn/data/se_protein/go_per_se.tsv", header=F, sep="\t", quote="", comment.char="")

colnames(df) <- c("se", "protein", "qv", "go_id", "name")

df$name <- sub(" ", "\n", df$name)

df <- subset(df, name != "?")

df$overall_rank <- rank(df$qv)

ggplot(df, aes(overall_rank, qv, fill=name, color=name)) + geom_point() + facet_grid(name~.) + scale_shape_manual(values=c(19,1)) + xlab("Rank") + ylab("q-value") + scale_y_log10()  + theme(legend.position="bottom")
