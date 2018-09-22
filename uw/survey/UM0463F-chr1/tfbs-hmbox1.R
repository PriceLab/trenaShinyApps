library(tfBindingSites)
library(BSgenome.Hsapiens.UCSC.hg38)
library(MotifDb)
source("~/github/trenaShinyApps/utils/plotMotifs.R")
motifs <- query(MotifDb, c("HMBOX1", "sapiens"))
length(motifs)
names(motifs)
plotMotifs(motifs)
motifs.selected <- MotifDb["Hsapiens-jaspar2018-HMBOX1-MA0895.1"]
export(motifs.selected, con="~/github/fimoService/pfms/hmbox1.human.meme", format="meme")
# buffer sh-fimo-server
# add to ~/github/fimoService/server/makefile.pshannon
#   MOTIFS="../pfms/nfkb2.human.meme"
# make -f makefile.pshannon

library(igvR)
igv <- igvR()
setGenome(igv, "hg38")

targetGene <- "PBX1"
chromosome <- "chr1"
snp.region.start <- 163448747
snp.region.end   <- 170558466
shoulder <- 10000
big.roi <- sprintf("%s:%d-%d", chromosome, snp.region.start - shoulder, snp.region.end + shoulder)

showGenomicRegion(igv, big.roi)

print(load("tbl.RData"))
tbl.snps <- tbl

tbl.wig <- data.frame(chr=rep("chr1", nrow(tbl.snps)),
                      start=tbl.snps$bp.hg38,
                      end=tbl.snps$bp.hg38,
                      value=-log10(tbl.snps$PVAL),
                      stringsAsFactors=FALSE)
label <- "-log10(PVAL)"
track <- DataFrameQuantitativeTrack(label, tbl.wig, color="darkRed", autoscale=TRUE)
displayTrack(igv, track)





print(load("tbl.enhancers.umo463f.RData"))
tbl.enhancers.pbx1 <- subset(tbl.enhancers.umo463f, geneSymbol==targetGene)
dim(tbl.enhancers.pbx1)   # 61 rows

tbl.wig <- tbl.enhancers.pbx1[, c("chr", "enhancer_start", "enhancer_end", "combined_score")]
colnames(tbl.wig) <- c("chr", "start", "end", "value")
label <- "pbx1 enhancers"
track <- DataFrameQuantitativeTrack(label, tbl.wig, color="black", autoscale=FALSE, min=0, max=30, trackHeight=25)
displayTrack(igv, track)


load("tbl.dhs.RData")
dim(tbl.dhs)
colnames(tbl.dhs)
tbl.wig <- tbl.dhs[, c("chrom", "chromStart", "chromEnd", "score")]
colnames(tbl.wig) <- c("chr", "start", "end", "value")
dim(tbl.wig)
label <- "dhs"
track <- DataFrameQuantitativeTrack(label, tbl.wig, color="darkGreen", autoscale=TRUE, trackHeight=25)
displayTrack(igv, track)



tf <- "HMBOX1"
tfbs <- tfBindingSitesCtor(tf,
                           genome=BSgenome.Hsapiens.UCSC.hg38,
                           chipseq.db.file="~/s/data/public/human/remap-2018/remap-all.sqlite",
                           fimoHost="localHost",
                           fimoPort=5558,
                           quiet=FALSE)

   #----------------------------------------
   # get and display the ChIP-seq hits
   #----------------------------------------

tbl.cs.hmbox1 <- getChipSeqHits(tfbs, chromosome, snp.region.start, snp.region.end)
dim(tbl.cs.hmbox1)
head(tbl.cs.hmbox1)
track <- DataFrameAnnotationTrack("HMBOX1-ChIP", tbl.cs.hmbox1[, c(1:3, 5)], color="purple", trackHeight=25)
displayTrack(igv, track)

   #-----------------------------------------------
   # get and display fimo binding sites for HMBOX1
   #-----------------------------------------------

tbl.pbx1 <- data.frame(chr=chromosome, start=164482061, end=164869549, stringsAsFactors=FALSE)
#tbl.fimo.cs <- getFimoHits(tfbs, tbl.cs.hmbox1, 10e-3)
tbl.fimo.hmbox1 <- getFimoHits(tfbs, tbl.pbx1, 10e-3)
dim(tbl.fimo.hmbox1)
tbl.wig <- tbl.fimo.hmbox1[, c("chr", "start", "stop", "pvalScore")]
colnames(tbl.wig) <- c("chr", "start", "end", "value")
track <- DataFrameQuantitativeTrack("hmbox-fimo", tbl.wig, "darkRed", trackHeight=25, autoscale=TRUE)
displayTrack(igv, track)


tbl.bioc.hmbox1 <- getMotifMatcherHits(tfbs, tbl.pbx1, motifs.selected, matchThreshold=0.65)
dim(tbl.bioc.hmbox1)
tbl.wig <- tbl.bioc.hmbox1[, c("chrom", "motifStart", "motifEnd", "motifRelativeScore")]
colnames(tbl.wig) <- c("chr", "start", "end", "value")
tbl.wig$value <- as.integer(100 * tbl.wig$value)
track <- DataFrameQuantitativeTrack("hmbox1-bioc", tbl.wig, "purple", trackHeight=25, autoscale=TRUE)
displayTrack(igv, track)

tbl.snp <- data.frame(chr=rep("chr1", nrow(tbl)), start=tbl.snps$bp.hg38, end=tbl.snps$bp.hg38,
                      value=-log10(tbl.snps$PVAL), stringsAsFactors=FALSE)
dim(tbl.snp)

tbl.ov <- as.data.frame(findOverlaps(GRanges(tbl.snp), GRanges(tbl.wig)))
dim(tbl.ov)  # 9 2
colnames(tbl.ov) <- c("snp", "hmbox.bs")

length(unique(tbl.ov$snp))  # 33
tbl.snp.hits <- tbl.snp[unique(tbl.ov$snp),]
track <- DataFrameAnnotationTrack("snp.in.hmbox1.bs", tbl.snp.hits, color="darkBrown", trackHeight=25)
displayTrack(igv, track)

tbl.bs.hits <- tbl.wig[unique(tbl.ov$hmbox.bs),]
track <- DataFrameAnnotationTrack("hmbox1.bs.snp", tbl.bs.hits, color="darkBrown", trackHeight=25)
displayTrack(igv, track)


