library(igvR)
library(RColorBrewer)
#------------------------------------------------------------------------------------------------------------------------
load("tbl.RData")
igv <- igvR()
setGenome(igv, "hg38")
showGenomicRegion(igv, "chr1")
tbl.wig <- data.frame(chr=rep("chr1", nrow(tbl)), start=tbl$bp.hg38, end=tbl$bp.hg38, value=-log10(tbl$PVAL), stringsAsFactors=FALSE)
label <- "-log10(PVAL)"
track <- DataFrameQuantitativeTrack(label, tbl.wig, color="darkRed", autoscale=TRUE)
displayTrack(igv, track)

tbl.wig <- data.frame(chr=rep("chr1", nrow(tbl)), start=tbl$bp.hg38, end=tbl$bp.hg38, value=-log10(tbl$PVAL.POP), stringsAsFactors=FALSE)
label <- "-log10(PVAL.POP)"
track <- DataFrameQuantitativeTrack(label, tbl.wig, color="darkOrange", autoscale=TRUE)
displayTrack(igv, track)


tbl.wig <- data.frame(chr=rep("chr1", nrow(tbl)), start=tbl$bp.hg38, end=tbl$bp.hg38, value=tbl$MAF, stringsAsFactors=FALSE)
label <- "MAF"
track <- DataFrameQuantitativeTrack(label, tbl.wig, color="darkGreen", autoscale=TRUE)
displayTrack(igv, track)

tbl.wig <- data.frame(chr=rep("chr1", nrow(tbl)), start=tbl$bp.hg38, end=tbl$bp.hg38, value=tbl$MAF.white, stringsAsFactors=FALSE)
label <- "MAF.white"
track <- DataFrameQuantitativeTrack(label, tbl.wig, color="blue", autoscale=TRUE)
displayTrack(igv, track)

print(load("tbl.enhancers.umo463f.RData"))

tbl.wig <- tbl.enhancers.umo463f[, c("chr", "enhancer_start", "enhancer_end", "combined_score")]
colnames(tbl.wig) <- c("chr", "start", "end", "value")
label <- "enhancers"
track <- DataFrameQuantitativeTrack(label, tbl.wig, color="black", autoscale=TRUE)
displayTrack(igv, track)

load("tbl.enhancers.chr1.RData")

tbl.wig <- tbl.enhancers.chr1[, c("chr", "enhancer_start", "enhancer_end", "combined_score")]
colnames(tbl.wig) <- c("chr", "start", "end", "value")
dim(tbl.wig)
label <- "enhancers.all"
track <- DataFrameQuantitativeTrack(label, tbl.wig, color="darkGreen", autoscale=TRUE)
displayTrack(igv, track)

load("tbl.dhs.RData")
dim(tbl.dhs)
colnames(tbl.dhs)
tbl.wig <- tbl.dhs[, c("chrom", "chromStart", "chromEnd", "score")]
colnames(tbl.wig) <- c("chr", "start", "end", "value")
dim(tbl.wig)
label <- "dhs"
track <- DataFrameQuantitativeTrack(label, tbl.wig, color="darkGreen", autoscale=TRUE)
displayTrack(igv, track)

showGenomicRegion(igv, "chr1:162,671,057-172,522,660")


tbl.xtab.enhancers <- as.data.frame(table(tbl.enhancers.umo463f$geneSymbol))
tbl.xtab.enhancers <- tbl.xtab.enhancers[order(tbl.xtab.enhancers$Freq, decreasing=TRUE),]
head(tbl.xtab.enhancers)
top.targets <- as.character(tbl.xtab.enhancers$Var1) # [1:5])
length(top.targets)

top.targets <- "PBX1"

colors <- brewer.pal(8, "Dark2")
totalColorCount <- length(colors)
colorNumber <- 0

for(target in top.targets){
   tbl.sub <- subset(tbl.enhancers.umo463f, geneSymbol == target)[, c("chr", "enhancer_start", "enhancer_end", "combined_score")]
   colnames(tbl.sub) <- c("chr", "start", "end", "value")
   colorNumber <- (colorNumber %% totalColorCount) + 1
   color <- colors[colorNumber]
   track <- DataFrameQuantitativeTrack(target, tbl.sub, color=color, autoscale=FALSE, min=0, max=20, trackHeight=25)
   displayTrack(igv, track)
   printf("%20s: %4d enhancer regions", target, nrow(tbl.sub))
   }


print(load("tbl.cs.mpzl1.RData"))
dim(tbl.cs);
head(tbl.cs)

cs.tfs <- unique(tbl.cs$tf)
length(cs.tfs)
for(cs.tf in cs.tfs){
   tbl.sub <- subset(tbl.cs, tf==cs.tf)
   tbl.sub <- tbl.sub[, c("chr", "start", "end", "name")]
   colorNumber <- (colorNumber %% totalColorCount) + 1
   color <- colors[colorNumber]
   track <- DataFrameAnnotationTrack(cs.tf, tbl.sub, color, trackHeight=30)
   displayTrack(igv, track)
   }


load("tbl.bindingSites.crem.RData")
track <- DataFrameQuantitativeTrack("crem2.bs", subset(tbl.wig, value > 3), color="green", trackHeight=30, autoscale=TRUE)
displayTrack(igv, track)


library(RSQLite)
db <- dbConnect(dbDriver("SQLite"), "~/s/data/public/human/remap-2018/remap-all.sqlite")
tbl.crem2 <- dbGetQuery(db, "select * from chipseq where tf='CREM2'")


print(load("/Users/paul/github/trenaSGM/inst/batch/cory/AD/tbl.models.all.RData"))
tbl.model.pbx1 <- head(subset(tbl.models, targetGene.hgnc=="PBX1"), n=10)
tbl.cs.nfkb1 <- dbGetQuery(db, "select * from chipseq where tf='NFKB1' and chr='chr1' and start > 161829985 and end < 172094633")
dim(tbl.cs.nfkb1)

track <- DataFrameAnnotationTrack("NKKB1 ChIP", tbl.cs.nfkb1, color="purple", trackHeight=30)
displayTrack(igv, track)

tbl.cs <- dbGetQuery(db, "select * from chipseq where tf='HMBOX1' and chr='chr1' and start > 161829985 and end < 172094633")
dim(tbl.cs)
track <- DataFrameAnnotationTrack("HMBOX1 ChIP", tbl.cs, color="purple", trackHeight=30)
displayTrack(igv, track)
