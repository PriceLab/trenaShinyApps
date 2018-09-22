library(trena)
source("constants.R")
hdf <- HumanDHSFilter("hg38", "wgEncodeRegDnaseClustered", pwmMatchPercentageThreshold=0,
                       geneInfoDatabase.uri="bogus", regions=data.frame(), pfms=list())
chrom <- "chr1"
start <- snp.min.loc - 5000
end <-   snp.max.loc + 5000
tbl.dhs <- getRegulatoryRegions(hdf, "wgEncodeRegDnaseClustered", chrom, start, end)
dim(tbl.dhs)
save(tbl.dhs, file="tbl.dhs.RData")
