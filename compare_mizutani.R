library(ggplot2)
library(reshape2)

df1 <- read.table("/Volumes/xi.embl.de/g/bork8/mkuhn/data/se_protein/mizutani/filtered_protein_se_pv.tsv", header=T, sep="\t", quote="", comment.char="")
df2 <- read.table("/Volumes/xi.embl.de/g/bork8/mkuhn/data/se_protein/mizutani/protein_se_weights.tsv", header=T, sep="\t", quote="", comment.char="")

shared_proteins <- intersect(df1$protein, df2$protein)
shared_se <- intersect(df1$se, df2$se)

df1.common <- subset(df1, protein %in% shared_proteins & se %in% shared_se)
df2.common <- subset(df2, protein %in% shared_proteins & se %in% shared_se)

df <- merge(df1.common, df2.common, all=T)

ggplot(df, aes(q, weight)) + geom_point() + scale_x_log10() + xlab("q-value") + ylab("SCCA weight") 

cor(df$q, df$weight, use="complete")

df.mouse <- read.table("/Volumes/xi.embl.de/g/bork8/mkuhn/data/se_protein/mouse_se.tsv", header=F, sep="\t", quote="", comment.char="")[,1:2]
colnames(df.mouse) <- c("se", "protein")
df.mouse$mouse <- T

df.cardiovascular <- read.table("/Volumes/xi.embl.de/g/bork8/mkuhn/data/se_protein/whitebread.tsv", header=F, sep="\t", quote="", comment.char="")
colnames(df.cardiovascular) <- c("se", "protein")
df.cardiovascular$cardiovascular <- T

df.literature <- read.table("/Volumes/xi.embl.de/g/bork8/mkuhn/data/se_protein/caused.tsv", header=F, sep="\t", quote="", comment.char="")
colnames(df.literature) <- c("se", "protein")
df.literature$literature <- T

dfX <- merge(df, df.mouse, all=T)
dfX <- merge(dfX, df.cardiovascular, all=T)
dfX <- merge(dfX, df.literature, all=T)

dfX$p <- NULL
dfX$protein_name <- NULL
dfX$se_name <- NULL

dfX[ is.na(dfX$mouse), "mouse" ] <- F
dfX[ is.na(dfX$cardiovascular), "cardiovascular" ] <- F
dfX[ is.na(dfX$literature), "literature" ] <- F

dfX$weight <- sqrt(dfX$weight)

dfXcomplete <- dfX

dfXshared <- subset(dfXcomplete, !is.na(q) & !is.na(weight))

dfX <- melt(dfX, id.vars=c("se", "protein", "q", "weight"))
colnames(dfX)[5:6] <- c("dataset", "verified")
dfX <- melt(dfX, id.vars=c("se", "protein", "dataset", "verified"))

dfXshared <- melt(dfXshared, id.vars=c("se", "protein", "q", "weight"))
colnames(dfXshared)[5:6] <- c("dataset", "verified")
dfXshared <- melt(dfXshared, id.vars=c("se", "protein", "dataset", "verified"))

pv <- ddply(dfX, .(dataset, variable), function(df) ks.test(subset(df, verified==T)$value, subset(df, verified==F)$value)$p.value)
pv.shared <- ddply(dfXshared, .(dataset, variable), function(df) ks.test(subset(df, verified==T)$value, subset(df, verified==F)$value)$p.value)

dataset_labels <- list(
    "q" = "q-value",
    "weight" = "SCCA weight"
)

variable_labels <- list(
    "mouse" = "mouse\nphenotypes",
    "cardiovascular" = "cardiovascular\nside effects",
    "literature" = "literature-derived\nside effects"
)

facet_labeller <- function(variable,value){
    if (variable == "variable") {
        return(dataset_labels[value])
    }
    return(variable_labels[value])
}

p <- ggplot(dfX, aes(value, color=verified)) + geom_density() + facet_grid(dataset~variable, scales="free", labeller=facet_labeller)
p.shared <- ggplot(dfXshared, aes(value, color=verified)) + geom_density() + facet_grid(dataset~variable, scales="free", labeller=facet_labeller)

p <- p + scale_color_discrete(breaks=c(F,T), labels=c("no", "yes"))
p.shared <- p.shared + scale_color_discrete(breaks=c(F,T), labels=c("no", "yes"))

pv$label <- sprintf("p = %.2g", pv$V1)
pv.shared$label <- sprintf("p = %.2g", pv.shared$V1)

p + geom_text(data=pv, aes(0.1,0.5,label=label), color="black", hjust=0)
p.shared + geom_text(data=pv.shared, aes(0.1,5,label=label), color="black", hjust=0)



df.best_pred <- read.table("/Volumes/xi.embl.de/g/bork8/mkuhn/data/se_protein/go_per_se.tsv", header=F, sep="\t", quote="", comment.char="")
colnames(df.best_pred) <- c("se", "protein", "qv", "go_id", "name")

df.best_pred$kind <- "predicted"
df.mouse$kind <- "mouse phenotypes"
df.cardiovascular$kind <- "cardiovascular SE"
df.caused$kind <- "literature SE"

cols <- c("protein", "se", "kind")

df.merged <- rbind(df.best_pred[,cols], df.mouse[,cols], df.cardiovascular[,cols], df.caused[,cols])

df.stats <- ddply(df.merged, .(kind), summarise, n_proteins_per_kind = length(unique(protein)), n_se_per_kind = length(unique(se)))

abundance <- ddply(df.merged, .(kind, protein), summarise, n_se=length(se))
abundance.combined <- ddply(abundance, .(kind, n_se), summarise, n_protein=length(protein))

abundance.combined <- merge(abundance.combined, df.stats)

ggplot(abundance.combined, aes(n_se, n_protein, color=kind)) + geom_point() + scale_y_log10(breaks=c(1,5,10,25,50)) + scale_x_log10(breaks=c(1,5,10,25,50)) + geom_smooth(method="lm", se=F) + xlab("Number of side effect per protein") + ylab("Count of proteins")
ggplot(abundance.combined, aes(n_se, n_protein/n_proteins_per_kind, color=kind, shape=kind)) + geom_point() + scale_y_log10() + scale_x_log10(breaks=c(1,5,10,25,50)) + geom_smooth(method="lm", se=F) + xlab("Number of side effect per protein") + ylab("Fraction of proteins")


