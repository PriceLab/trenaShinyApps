f <- "../variants/fam_based/EXT/UM0463F.results.all.txt"
tbl <- read.table(f, header=TRUE, nrows=-1, stringsAsFactors=FALSE)
dim(tbl)
head(tbl$rsID)
as.data.frame(t(tbl[3,]))
163785631 - 163754868
dim(subset(tbl, CHR.x == 1))   # all chromosome 1
head(tbl[, c("CHR.x", "BP", "BP")])
hg19.locs <- with(tbl, sprintf("chr%s:%d-%d", CHR.x, BP, BP))
write(hg19.locs, file="hg19.locs.txt", sep="\n")
tbl.hg38 <- read.table("hglft_genome_199d_9ae190.bed", stringsAsFactors=FALSE)
hg38.loc <- tbl.hg38$V1
hg38.loc <- sub("chr1:", "", hg38.loc)
tokens <- strsplit(hg38.loc, "-")
x <- as.numeric(unlist(lapply(tokens, "[", 1)))
tbl$bp.hg38 <- x
save(tbl, file="tbl.RData")

