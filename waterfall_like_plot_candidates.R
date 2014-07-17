library(ggplot2)

groups <- data.frame(status=c("literature support", "unknown", "related phenotype", "same phenotype", "anecdotal evidence", "hypersensitivity reaction", "wrong",              "false positive", "opposite phenotype", "not verified"), 
                     group =c("known",              "unknown", "known"            , "known"         , "known"             , "erroneous association" , "erroneous association", "erroneous association", "known"      , "not verified"        )
                     )

# df <- read.table("/Volumes/xi.embl.de/g/bork8/mkuhn/data/se_protein/pred.mouse.0.01.tsv", header=T, sep="\t", quote="", comment.char="")
df <- read.table("/Volumes/xi.embl.de/g/bork8/mkuhn/data/se_protein/pred.full_merged.0.01.tsv", header=T, sep="\t", quote="", comment.char="")

df <- subset(df, only.metabolizing. == "")

colnames(df)[1] <- "qv"
# df <- ddply(df, .(UMLS.code.of.SE), function(df) df[ order(df$qv)[1], ])
df <- merge(df, groups, all=T)

df$candidate <- "no"
df[df$side.effect == "hyperesthesia",]$candidate <- "yes"
df[df$side.effect == "hyperesthesia",]$group <- "unknown"

df <- df[order(df$qv),]
df$overall_rank <- 1:nrow(df)

ggplot(df, aes(overall_rank, qv, fill=status, color=status, size=candidate)) + geom_point() + facet_grid(group~.) + scale_size_manual(values=c(2,6)) + scale_shape_manual(values=c(19,1)) + xlab("Rank") + ylab("log10 of q-value")
