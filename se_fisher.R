library(plyr)
# library(reshape2)
library(hash)

# source("http://bioconductor.org/biocLite.R")
# biocLite("qvalue")

library(qvalue)

setwd("/Volumes/xi.embl.de/g/bork8/mkuhn/data/se_protein")

drug_protein <- read.table("interactions.tsv", as.is=TRUE, sep="\t", quote="")
names(drug_protein) <- c("drug", "protein")

protein_freq <- count(drug_protein$protein)
colnames(protein_freq) <- c("protein", "protein_freq")

drug_se <- read.table("side_effects.tsv", as.is=TRUE, sep="\t", quote="")
names(drug_se) <- c("drug", "se", "drug_name", "se_name")

se_freq <- count(drug_se$se)
colnames(se_freq) <- c("se", "protein_freq")

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

# protein_se_pv <- merge(protein_se_pv, se_freq)
# protein_se_pv <- merge(protein_se_pv, protein_freq)
# 
# min_pv_per_se <- ddply(protein_se_pv, .(se), summarise, min_p = min(p), se_freq = min(se_freq))
# ggplot(min_pv_per_se, aes(se_freq, min_p)) + geom_point() + scale_y_log10() + scale_x_log10()
# min_pv_per_protein <- ddply(protein_se_pv, .(protein), summarise, min_p = min(p), protein_freq = min(protein_freq))
# ggplot(min_pv_per_protein, aes(protein_freq, min_p)) + geom_point() + scale_y_log10() + scale_x_log10()
