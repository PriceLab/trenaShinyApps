# minerva's snp: rs9357347
# contained in enhancer region: chr6:41182853
# what binds there?
#------------------------------------------------------------------------------------------------------------------------
library(trena)
library(motifStack)
library(igvR)
#------------------------------------------------------------------------------------------------------------------------
pfms <- query(MotifDb, "sapiens", c("jaspar2018", "hocomoco"))
mm <- MotifMatcher("hg38", as.list(pfms), quiet=FALSE)
chrom <- "chr6"
snp.loc <- 41182853
start <- snp.loc - 20
end   <- snp.loc + 20
tbl.regions <- data.frame(chrom=chrom, start=start, end=end, stringsAsFactors=FALSE)
tbl.hits <- findMatchesByChromosomalRegion(mm, tbl.regions, 80, variants = NA_character_)
tbl.hits <- associateTranscriptionFactors(MotifDb, tbl.hits, source=c("motifdb", "tfclass"), expand.rows=TRUE)

goi <- "SPI1"
tbl.goi <- subset(tbl.hits, geneSymbol==goi)
dim(tbl.goi)

goi.associated.motifs <- unique(tbl.goi$shortMotif)
goi.pfms <- query(MotifDb, c("sapiens", "jaspar2018"), goi.associated.motifs)
if(length(goi.pfms) == 1)
   goi.pfms <- rep(goi.pfms, 2)

motifStack(lapply(names(goi.pfms), function(mName) new("pfm", goi.pfms[[mName]], name=mName)))

tbl.snp <- data.frame(chrom="chr6", start=snp.loc, end=snp.loc, stringsAsFactors=FALSE)

igv <- igvR()
setGenome(igv, "hg38")
showGenomicRegion(igv, sprintf("%s:%d-%d", chrom, start, end))
track <- DataFrameAnnotationTrack("snp", tbl.snp, trackHeight=30, color="red")
displayTrack(igv, track)

track <- DataFrameAnnotationTrack(goi, tbl.goi[, c("chrom", "motifStart", "motifEnd", "motifName")], displayMode="EXPANDED")
displayTrack(igv, track)
