# From Genbank/TAIR data we know that there are 8 transcripts for the gene, with 3 different
# transcription start sites that are not far from each other - but coding just two different proteins:
#
#    https://www.ncbi.nlm.nih.gov/gene/818427
#    https://apps.araport.org/thalemine/portal.do?externalids=AT2G38460
#
# For now I will only be doing RT-PCRs
# with primers that do not distinguish between the 8 transcripts to see if there is any general
# upregulation of IREG1 transcripts. On your end you might have to have 3 models.
# https://www.ncbi.nlm.nih.gov/gene/818427
#------------------------------------------------------------------------------------------------------------------------
library(trena)
library(trenaSGM)
library(igvR)
library(RColorBrewer)
packageVersion("igvR")
library(GenomicRanges)
library(MotifDb)
library(FimoClient)
library(RUnit)
genomeName <- "BSgenome.Athaliana.TAIR.TAIR9"
library(genomeName, character.only=TRUE)
#------------------------------------------------------------------------------------------------------------------------
source("../../../utils/roundNumericColumnsInDataframe.R")
#------------------------------------------------------------------------------------------------------------------------
geneSymbol <- "IREG1"
orf <- "AT2G38460"

tss.1 <- 16103164
tss.2 <- 16103603
targetGene <- orf
chromosome <- "chr2"
loc.start <- 16103164
loc.end <- 16106177
chrom.loc.string <- sprintf("%s:%d-%d", chromosome, loc.start - 1000, loc.end + 1000)
#------------------------------------------------------------------------------------------------------------------------
if(!exists("tbl.bud.dhs")){
    extra.bases <- 8000
    f <- "../../arabidopsisData/Ath_buds_DNase.bw"
    file.exists(f)
    buds.gr <- import(f, which=seqinfo(BigWigFile(f))["Chr2"])
    class(buds.gr)
    tbl.buds.dhs <- as.data.frame(buds.gr)
    colnames(tbl.buds.dhs) <- c("chr", "start", "end", "width", "strand", "value")
    tbl.buds.dhs <- tbl.buds.dhs[, c("chr", "start", "end", "value")]
    fivenum(tbl.buds.dhs$value)  # [1] 0.000 0.410 0.670 1.009 8.610
    tbl.buds.dhs <- subset(tbl.buds.dhs, start > (loc.start - extra.bases) & end < (loc.end + extra.bases))
    tbl.buds.dhs$chr <- as.character(tbl.buds.dhs$chr)
    dim(tbl.buds.dhs)
    fivenum(tbl.buds.dhs$value)  # [1] 0.000 0.410 0.670 1.009 8.610

    f <- "../../arabidopsisData/Ath_leaf_DNase.bw"
    file.exists(f)
    leaf.gr <- import(f, which=seqinfo(BigWigFile(f))["Chr2"])
    class(leaf.gr)
    tbl.leaf.dhs <- as.data.frame(leaf.gr)
    dim(tbl.leaf.dhs)
    fivenum(tbl.leaf.dhs$score)
    colnames(tbl.leaf.dhs) <- c("chr", "start", "end", "width", "strand", "value")
    tbl.leaf.dhs <- tbl.leaf.dhs[, c("chr", "start", "end", "value")]

    fivenum(tbl.leaf.dhs$value)  # 0.000 0.500 0.760 1.060 7.492
    tbl.leaf.dhs <- subset(tbl.leaf.dhs, start > (loc.start - extra.bases) & end < (loc.end + extra.bases))
    tbl.leaf.dhs$chr <- as.character(tbl.leaf.dhs$chr)
    dim(tbl.leaf.dhs)
    fivenum(tbl.leaf.dhs$value)  # [1] 0.000 0.410 0.670 1.009 8.610
    save(tbl.leaf.dhs, tbl.buds.dhs, file="../shinyApp/data/tbls.dhs.RData")
    }

print(load("../shinyApp/data/tbls.dhs.RData"))

if(!exists("igv")){
   igv <- igvR()
   setBrowserWindowTitle(igv, geneSymbol)
   setGenome(igv, "tair10")
   showGenomicRegion(igv, chrom.loc.string)
   }
#------------------------------------------------------------------------------------------------------------------------
if(!exists("dhsTracks")){
   load("../shinyApp/data/tbls.dhs.RData")
   exists("tbl.buds.dhs"); exists("tbl.leaf.dhs")
   fivenum(tbl.buds.dhs$value)
   fivenum(tbl.leaf.dhs$value)
   track.buds <- DataFrameQuantitativeTrack("DHS buds", autoscale=FALSE, min=0, max=10, tbl.buds.dhs)
   displayTrack(igv, track.buds)
   fivenum(tbl.leaf.dhs$value)
   track.leaves <- DataFrameQuantitativeTrack("DHS leaves", autoscale=FALSE, min=0, max=10, tbl.leaf.dhs)
   displayTrack(igv, track.leaves)
   }
#------------------------------------------------------------------------------------------------------------------------
if(!exists("mtx")){
   print(load("../shinyApp/data/mtx.zinc.22810x42.RData")) # mtx
   print(load("../shinyApp/data/tbl.xref.RData"))          # "tbl.upOrf" "tbl.orfSym" "tbl.upOldgene" "tbl.xref"
   }
#------------------------------------------------------------------------------------------------------------------------
meme.file <- "~/github/fimoService/pfms/athaliana.all.meme"

if(!file.exists(meme.file)){
   # create this file
   library(MotifDb)
   pfms <- query(MotifDb, "athaliana", "jaspar2018")
   export(pfms, meme.file, 'meme')
   # restart fimo server
   # cd ~/github/fimoService/server
   # in makefile.pshannon:
   #    PORT=5558
   #    MOTIFS="../pfms/athaliana.all.meme"
   #  make -f makefile.pshannon
   #    python -i runFimoServer.py 5558 "/Users/paul/meme/bin/fimo" "../pfms/athaliana.all.meme"
   }
if(!exists("fc")){
   FIMO_HOST <- "localhost"
   FIMO_PORT <- 5558
   fc <- FimoClient(FIMO_HOST, FIMO_PORT, quiet=FALSE)
   }
#------------------------------------------------------------------------------------------------------------------------
if(!exists("candidate.tfs.1")){
   start <- tss.1 - 1999
   end   <- tss.1 + 500
   sequence <- as.character(getSeq(BSgenome.Athaliana.TAIR.TAIR9, "Chr2", start, end))
   nchar(sequence)  # 2500
   fimo.pval.threshold <- 10e-5
   tbl.fimo <- requestMatch(fc, list(frd3=sequence), pvalThreshold=fimo.pval.threshold)
   dim(tbl.fimo)
   tbl.fimo$start <- tbl.fimo$start + start
   tbl.fimo$stop <- tbl.fimo$stop + start
   tbl.fimo$chrom <- chromosome
   tbl.bed <- tbl.fimo[, c("chrom", "start", "stop", "score", "motif")]
   dim(tbl.bed)
   tbl.xtab <- as.data.frame(table(tbl.fimo$motif))
   fivenum(tbl.xtab$Freq)
   candidate.motifs <- unique(tbl.fimo$motif)
   length(candidate.motifs) # 198
   tfs.1 <- as.data.frame(mcols(MotifDb[candidate.motifs]))$geneSymbol
   tfs.orfs.1 <- tbl.xref$orfs[match(tfs.1, tbl.xref$geneSymbol)]
   tfs.orfs.2 <- tbl.xref$orfs[match(tfs.1, tbl.xref$SYMBOL)]
   length(tfs.orfs.1)
   length(which(is.na(tfs.orfs.1)))
   length(tfs.orfs.2)
   length(which(is.na(tfs.orfs.2)))
   tfs.orfs <- unique(c(tfs.orfs.1, tfs.orfs.2))
   deleters <- which(is.na(tfs.orfs))
   printf("%d/%d tf.orfs are NA", length(deleters), length(tfs.orfs))
   if(length(deleters) > 0)
      tfs.orfs <- tfs.orfs[-deleters]
   candidate.tfs.1 <- tfs.orfs
   length(candidate.tfs.1) # 168
   } # tfs


if(!exists("tbl.model")) {
   length(intersect(candidate.tfs.1, rownames(mtx))) # 155
   solver <- EnsembleSolver(mtx, targetGene, candidate.tfs.1)
   tbl.model <- run(solver)
   tbl.model$geneSymbol <- tbl.xref[match(tbl.model$gene, tbl.xref$orfs), "geneSymbol"]
   tbl.model$SYMBOL <- tbl.xref[match(tbl.model$gene, tbl.xref$orfs), "SYMBOL"]
   tbl.model <- tbl.model[order(tbl.model$pcaMax, decreasing=TRUE),]
   tbl.model <- roundNumericColumnsInDataframe(tbl.model, 3, "lassoPValue")
   tbl.model <- subset(tbl.model, abs(pearsonCoeff) > 0.4)
   if(nrow(tbl.model) > 10)
      tbl.model <- tbl.model[1:10,]
   dim(tbl.model)
   save(tbl.model, file="../shinyApp/data/tbl.model.RData")
   }  # tbl.model

add.orfs.geneSymbol.to.fimoResults <- function(tbl.bed){
   motifs <- tbl.bed$motif
   motifDB.geneSymbols <- unlist(lapply(motifs, function(motif) mcols(MotifDb[motif])$geneSymbol))
   length(motifs)
   length(motifDB.geneSymbols)
   length(which(is.na(motifDB.geneSymbols)))
   xref.indices <- match(motifDB.geneSymbols, tbl.xref$geneSymbol)
   orfs <- tbl.xref$orfs[xref.indices]
   geneSymbols <- tbl.xref$geneSymbol[xref.indices]
   tbl.bed$orf <- orfs
   tbl.bed$geneSymbol <- geneSymbols
   tbl.bed
   }

if(!exists("tbl.bindingSites")){
   colors <- brewer.pal(8, "Dark2")
   totalColorCount <- length(colors)
   colorNumber <- 0
   tbl.bindingSites <- add.orfs.geneSymbol.to.fimoResults(tbl.bed)
   tbl.bindingSites <-  subset(tbl.bindingSites, orf %in% tbl.model$gene)
   dim(tbl.bindingSites)  # 19 7
   tfs.unique <- sort(unique(tbl.bindingSites$geneSymbol))
   length(tfs.unique) # 10
   for(tf in tfs.unique){
      tbl.tfbs <- subset(tbl.bindingSites, geneSymbol==tf)[, c("chrom", "start", "stop", "score", "geneSymbol")]
      colnames(tbl.tfbs) <- c("chr", "start", "end", "value", "tf")
      colorNumber <- (colorNumber %% totalColorCount) + 1
      color <- colors[colorNumber]
      track <- DataFrameQuantitativeTrack(tf, tbl.tfbs, color=color, trackHeight=30, autoscale=FALSE, min=0, max=20)
      displayTrack(igv, track)
      } # for tf
   save(tbl.bindingSites, file="../shinyApp/data/tbl.bindingSites.RData")
   } # tbl.bindingSites


if(!exists("tbl.bindingSites.2")){
   tbl.bindingSites.2 <- add.orfs.geneSymbol.to.fimoResults(tbl.bed.2)
   tbl.bindingSites.2 <-  subset(tbl.bindingSites.2, orf %in% tbl.model.2$gene)
   dim(tbl.bindingSites.2)
   tfs.unique <- sort(unique(tbl.bindingSites.2$geneSymbol))
   for(tf in tfs.unique){
      tbl.tfbs <- subset(tbl.bindingSites.2, geneSymbol==tf)
      track <- DataFrameAnnotationTrack(sprintf("%s.2", tf), tbl.tfbs, color="red", trackHeight=22)
      displayTrack(igv, track)
      } # for tf
   save(tbl.bindingSites.2, file="../shinyApp/data/tbl.bindingSites.frd3-transcript2.RData")
   } # tbl.bindingSites.2



#------------------------------------------------------------------------------------------------------------------------
findRegionsAboveThreshold <- function(vec, threshold, minSpan) {
   vec[vec < threshold] <- NA
   vec[vec >= threshold] <- 1
   vec.rle <- rle(vec)
   rle.segments.above.threshold <- which(vec.rle$lengths >= minSpan)
       # make all runs less than desired.minSpan are NA'd
   all.subMinimumSpans <- setdiff(1:length(vec.rle$values), rle.segments.above.threshold)
   vec.rle$values[all.subMinimumSpans] <- NA    # vec.rle$values[18] <- NA

   rle.region.count <- length(vec.rle$values)
   actual.index <- 1

   peak.starts <- vector("numeric", length(vec))
   peak.ends <- vector("numeric", length(vec))
   peak.count <- 0

   for(i in 1:rle.region.count){
      size <- vec.rle$length[i]
      value <- vec.rle$values[i]
      if(!is.na(value)){
         peak.count <- peak.count + 1
         region.start <- actual.index
        region.end   <- actual.index + size - 1
        printf("peak found:  %d-%d", region.start, region.end)
        peak.starts[peak.count] <- region.start
        peak.ends[peak.count] <- region.end
        } # !is.na
    actual.index <- actual.index + size
   } # for i
   list(starts=peak.starts[1:peak.count], ends=peak.ends[1:peak.count])
} # findRegionsAboveThreshold
#------------------------------------------------------------------------------------------------------------------------
if(!exists("call-peaks")){
   x <- findRegionsAboveThreshold(tbl.buds.dhs$value, 5, 10)
   tbl.peaks <- with(x, data.frame(chr=rep("chr2", length(starts)), start=starts, end=ends), stringsAsFactors=FALSE)
   tbl.peaks$start <- tbl.peaks$start + tbl.buds.dhs$start[1]
   tbl.peaks$end <- tbl.peaks$end + tbl.buds.dhs$start[1]
   track <- DataFrameAnnotationTrack("bud peaks", tbl.peaks, color="darkGreen", trackHeight=25)
   displayTrack(igv, track)

}

