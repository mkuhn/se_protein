library(ggplot2)
library(plyr)

s <- "caused"
filename <- paste("/Volumes/xi.embl.de/g/bork8/mkuhn/data/se_protein/plot.",s,".se_protein_known.all.txt", sep="")
d1 <- read.table(filename, header=TRUE, sep='\t')
d1$dataset <- s

s <- "whitebread"
filename <- paste("/Volumes/xi.embl.de/g/bork8/mkuhn/data/se_protein/plot.",s,".se_protein_known.all.txt", sep="")
d2 <- read.table(filename, header=TRUE, sep='\t')
d2$dataset <- s

d <- rbind(d1, d2)

d <- subset(d, known!="yes")

cols <- c("no" = "black","best" = "orange") 
scale <- scale_colour_manual(values = cols, breaks = c("no", "best"), labels = c("background", "best action")) 

p2 <- ggplot(d, aes(x=pvalue, group=known)) + geom_density(fill = NA, position="dodge", aes(colour=known)) + scale + xlab("q-value of protein-side effect association")+theme(aspect.ratio=2) + facet_grid(.~dataset)
#+opts(aspect.ratio=2, title="Unbiased set of side effects")

print(p2)

# for supplement:
print(p2 + scale_x_log10())

ddply(d, .(dataset), function(df) ks.test(subset(df,known=="best")$pvalue, subset(df,known=="no")$pvalue, alternative="g" )$p.value)
#print( ks.test(best$fraction, not_known$fraction, alternative="l" ) )
