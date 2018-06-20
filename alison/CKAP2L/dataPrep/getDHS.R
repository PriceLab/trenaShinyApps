library(trena)
hdf <- HumanDHSFilter("hg38", "wgEncodeRegDnaseClustered", pwmMatchPercentageThreshold=0,
                      geneInfoDatabase.uri="bogus", regions=data.frame(), pfms=list())
tableNames <- getEncodeRegulatoryTableNames(hdf)
tbl.proximalEnhancerRegion <- data.frame(chrom="chr2", start=112719115, end=112905047)


chrom <- "chr2"
start <- 112719115
end <- 112905047

tbl.dhs <- getRegulatoryRegions(hdf, "wgEncodeRegDnaseClustered", chrom, start, end)
dim(tbl.dhs)
colnames(tbl.dhs) <- c("chrom", "start", "end", "count", "score")
head(tbl.dhs)
save(tbl.dhs, file="../data/tbl.dhs.RData")
