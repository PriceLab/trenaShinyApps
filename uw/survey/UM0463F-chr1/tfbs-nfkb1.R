library(tfBindingSites)
library(BSgenome.Hsapiens.UCSC.hg38)
library(MotifDb)
source("~/github/trenaShinyApps/utils/plotMotifs.R")
motifs <- query(MotifDb, c("NFKB1", "sapiens"))
length(motifs)
names(motifs)
plotMotifs(motifs)
# choose among those displayed, get their names
names.chosen <- c("Hsapiens-jaspar2018-NFKB1-MA0105.2", "Hsapiens-jaspar2018-NFKB1-MA0105.3", "Hsapiens-jaspar2018-NFKB1-MA0105.4")
motifs.selected <- MotifDb[names.chosen]
export(motifs.selected, con="~/github/fimoService/pfms/nfkb1.human.meme", format="meme")
# buffer sh-fimo-server
# add to ~/github/fimoService/server/makefile.pshannon
#   MOTIFS="../pfms/nfkb2.human.meme"
# make -f makefile.pshannon

tf <- "NFKB1"
tfbs <- tfBindingSitesCtor(tf,
                           genome=BSgenome.Hsapiens.UCSC.hg38,
                           chipseq.db.file="~/s/data/public/human/remap-2018/remap-all.sqlite",
                           fimoHost="localHost",
                           fimoPort=5558,
                           quiet=TRUE)
chromosome <- "chr1"
start <- 163243750
end <- 170894880
tbl.wholeRegion <- data.frame(chr=chromosome, start=start, end=end)

   #----------------------------------------
   # get and display the ChIP-seq hits
   #----------------------------------------

tbl.cs.nfkb1 <- getChipSeqHits(tfbs, chromosome, start, end)
dim(tbl.cs.nfkb1)
head(tbl.cs.nfkb1)
track <- DataFrameAnnotationTrack("NFKB1-tfbs", tbl.cs.nfkb1[, c(1:3, 5)], color="purple", trackHeight=25)
displayTrack(igv, track)

   #--------------------------------------------
   # get and display the fimo sequence matches
   #--------------------------------------------

tbl.fimo.cs <- getFimoHits(tfbs, tbl.cs.nfkb1, 10e-3)
dim(tbl.fimo.cs)
tbl.wig <- tbl.fimo[, c("chr", "start", "stop", "pvalScore")]
colnames(tbl.wig) <- c("chr", "start", "end", "value")
track <- DataFrameQuantitativeTrack("nfkb1-fimo-cs", tbl.wig, "darkRed", trackHeight=25, autoscale=TRUE)
displayTrack(igv, track)

tbl.fimo.all <- getFimoHits(tfbs, tbl.wholeRegion, 10e-3)
dim(tbl.fimo.all)
tbl.wig <- tbl.fimo.all[, c("chr", "start", "stop", "pvalScore")]
colnames(tbl.wig) <- c("chr", "start", "end", "value")
tbl.wig$chr <- as.character(tbl.wig$chr)
track <- DataFrameQuantitativeTrack("nfkb1-fimo-al", tbl.wig, "darkRed", trackHeight=25, autoscale=TRUE)
displayTrack(igv, track)



   #--------------------------------------------
   # get and display the bioc sequence matches
   #--------------------------------------------

tbl.bioc <- getMotifMatcherHits(tfbs, tbl.cs.nfkb1, motifs.selected, matchThreshold=0.65)
dim(tbl.bioc)
tbl.wig <- tbl.bioc[, c("chrom", "motifStart", "motifEnd", "motifRelativeScore")]
colnames(tbl.wig) <- c("chr", "start", "end", "value")
tbl.wig$value <- as.integer(100 * tbl.wig$value)
track <- DataFrameQuantitativeTrack("nfkb1-bioc", tbl.wig, "purple", trackHeight=25, autoscale=TRUE)
displayTrack(igv, track)

load("tbl.RData")
tbl.snp <- data.frame(chr=rep("chr1", nrow(tbl)), start=tbl$bp.hg38, end=tbl$bp.hg38, value=-log10(tbl$PVAL), stringsAsFactors=FALSE)
dim(tbl.snp)

tbl.ov <- as.data.frame(findOverlaps(GRanges(tbl.snp), GRanges(tbl.wig)))
dim(tbl.ov)
colnames(tbl.ov) <- c("snp", "bioc-motif")
length(unique(tbl.ov$snp))  # 33
tbl.snp.hits <- tbl.snp[unique(tbl.ov$snp),]
track <- DataFrameAnnotationTrack("snp.in.nfkb1.bs", tbl.snp.hits, color="darkBrown")
displayTrack(igv, track)




