library(plyr)
# library(reshape2)
library(hash)
library(ggplot2)

# source("http://bioconductor.org/biocLite.R")
# biocLite("qvalue")

library(qvalue)

setwd("/Volumes/xi.embl.de/g/bork8/mkuhn/data/se_protein")
setwd("/Volumes/xi.embl.de/g/bork8/mkuhn/data/se_protein/se_protein.cutoff2")

drug_protein <- read.table("interactions.tsv", as.is=TRUE, sep="\t", quote="")
names(drug_protein) <- c("drug", "protein")

protein_freq <- count(drug_protein$protein)
colnames(protein_freq) <- c("protein", "protein_freq")

drug_se <- read.table("side_effects.tsv", as.is=TRUE, sep="\t", quote="")
names(drug_se) <- c("drug", "se", "drug_name", "se_name")

se_freq <- count(drug_se$se)
colnames(se_freq) <- c("se", "se_freq")

counts_per_protein <- ddply(drug_protein, .(protein), function(df) nrow(df))
counts_per_protein <- hash( counts_per_protein$protein, counts_per_protein$V1 )

n_drugs <- length(unique(drug_protein$drug))

process_se_protein <- function(df, n_drugs_with_se)
{
    n_drugs_with_se_with_protein <- nrow(df)
    n_drugs_with_se_without_protein <- n_drugs_with_se - n_drugs_with_se_with_protein

    n_drugs_without_se_with_protein <- counts_per_protein[[ df[1,"protein"] ]] - n_drugs_with_se_with_protein
    n_drugs_without_se_without_protein <- n_drugs - n_drugs_without_se_with_protein - n_drugs_with_se
    
    m <- matrix(c(n_drugs_with_se_with_protein, n_drugs_with_se_without_protein, n_drugs_without_se_with_protein, n_drugs_without_se_without_protein), nrow=2 )
    
    fisher.test(m, alternative="g")$p.value
}

process_se <- function(df)
{
    n_drugs_with_se <- nrow(df)
    
    se_drug_proteins <- merge(df, drug_protein)
    
    ddply( se_drug_proteins, .(protein), process_se_protein, n_drugs_with_se = n_drugs_with_se)
}

protein_se_pv <- ddply( drug_se, .(se), process_se)
colnames(protein_se_pv)[3] <- "p"

q <- qvalue(protein_se_pv$p, gui=T)
protein_se_pv$q <- q$qvalues

write.table(protein_se_pv, file="protein_se_pv.tsv", quote=F, sep="\t", row.names=F)

###

protein_se_pv <- read.table("protein_se_pv.tsv", quote="", sep="\t", header=T)

protein_se_pv <- merge(protein_se_pv, se_freq)
protein_se_pv <- merge(protein_se_pv, protein_freq)

min_pv_per_se <- ddply(protein_se_pv, .(se), summarise, min_p = min(p), min_q = min(q), se_freq = min(se_freq))
min_pv_per_se.se <- ddply(min_pv_per_se, .(se_freq, min_q), summarise, point_count = length(se))

ggplot(min_pv_per_se.se, aes(se_freq, min_q, size=point_count)) + geom_point() + scale_y_log10(breaks=c(1e-15, 1e-10, 1e-5, 0.01, 1)) + scale_x_log10(breaks=c(2,5,20,100,500)) + xlab("Number of drugs per side effect") + ylab("Best q-value per side effect") + scale_size_area(max_size=20, breaks=c(1,10,50,100,200), name="Number of side effects\nwith same best q-value")  # + theme(legend.justification=c(1,0), legend.position=c(1,0))

min_pv_per_protein <- ddply(protein_se_pv, .(protein), summarise, min_p = min(p), min_q = min(q), protein_freq = min(protein_freq))
min_pv_per_protein.protein <- ddply(min_pv_per_protein, .(protein_freq, min_q), summarise, point_count = length(protein))

ggplot(min_pv_per_protein.protein, aes(protein_freq, min_q, size=point_count)) + geom_point()+ scale_y_log10(breaks=c(1e-15, 1e-10, 1e-5, 0.01, 1)) + scale_x_log10(breaks=c(2,5,20,100))+ xlab("Number of drugs per target") + ylab("Best q-value per target") + scale_size_area(max_size=20, breaks=c(1,10,50,100,200), name="Number of targets\nwith same best q-value") # + theme(legend.justification=c(1,0), legend.position=c(1,0))
