library(trena)
targetGene <- "INPP5D"
chromosome <- "chr2"

if(!exists("tbl.enhancers")){
   source("~/s/data/public/GeneHancer/v4.7/getAllEnhancers.R")
   tbl.enhancers <- getAllEnhancers(targetGene)
   colnames(tbl.enhancers) <- c("chr", "start", "end", "type", "score")
   dim(tbl.enhancers)
   gr.enhancers <- GRanges(tbl.enhancers)
   }

extra.loc <- 5000
min.loc <- min(tbl.enhancers$start) - extra.loc
max.loc <- max(tbl.enhancers$end) + extra.loc
max.loc - min.loc

hdf <- HumanDHSFilter("hg38", "wgEncodeRegDnaseClustered", pwmMatchPercentageThreshold=0,
                      geneInfoDatabase.uri="bogus", regions=data.frame(), pfms=list())

tbl.dhs <- getRegulatoryRegions(hdf, "wgEncodeRegDnaseClustered", chromosome, min.loc, max.loc)
dim(tbl.dhs)
head(tbl.dhs)
tbl.dhs <- tbl.dhs[, c("chrom", "chromStart", "chromEnd", "score")]
colnames(tbl.dhs) <- c("chr", "start", "end", "value")
head(tbl.dhs)
save(tbl.dhs, file="../shinyApp/data/tbl.dhs.RData")
