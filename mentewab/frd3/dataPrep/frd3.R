library(trena)
library(trenaSGM)
library(igvR)
library(GenomicRanges)
library(MotifDb)
library(FimoClient)
library(RUnit)
genomeName <- "BSgenome.Athaliana.TAIR.TAIR9"
library(genomeName, character.only=TRUE)
#------------------------------------------------------------------------------------------------------------------------
source("../../../utils/roundNumericColumnsInDataframe.R")
#------------------------------------------------------------------------------------------------------------------------
load("data/tbls.frd3.dhs.budAndLeaf.RData")   # "tbl.frd3buds" "tbl.frd3leaf"
#------------------------------------------------------------------------------------------------------------------------
geneSymbol <- "FRD3"
orf <- "AT3G08040"
tss.orf.1 <- 2569396
tss.orf.2 <- 2572151
targetGene <- orf
chromosome <- "chr3"
loc.start <- 2566277
loc.end <- 2572151
chrom.loc.string <- sprintf("%s:%d-%d", chromosome, loc.start - 1000, loc.end + 1000)
#------------------------------------------------------------------------------------------------------------------------
if(!exists("igv")){
   igv <- igvR()
   Sys.sleep(3)
   setBrowserWindowTitle(igv, "FRD3")
   setGenome(igv, "tair10")
   Sys.sleep(3)
   showGenomicRegion(igv, chrom.loc.string)
   Sys.sleep(5)
   }
#------------------------------------------------------------------------------------------------------------------------
if(!exists("tbl.frd3buds")){
   load("../shinyApp/data/tbls.frd3.dhs.budAndLeaf.RData")
   stopifnot(exists("tbl.frd3buds"))
   stopifnot(exists("tbl.frd3leaf"))
   track.buds <- DataFrameQuantitativeTrack("DHS buds", tbl.frd3buds)
   displayTrack(igv, track.buds)
   track.leaves <- DataFrameQuantitativeTrack("DHS leaves", tbl.frd3leaf)
   displayTrack(igv, track.leaves)
   }
#------------------------------------------------------------------------------------------------------------------------
if(!exists("mtx")){
   load("../shinyApp/data/mtx.zinc.22810x42.RData")
   load("../shinyApp/data/tbl.xref.RData")
   }
#------------------------------------------------------------------------------------------------------------------------
if(!exists("fc")){
   FIMO_HOST <- "localhost"
   FIMO_PORT <- 5558
   fc <- FimoClient(FIMO_HOST, FIMO_PORT, quiet=FALSE)
   }

fimo.pval.threshold <- 5 * 10e-5

if(!exists("candidate.tfs.1")){
     # assemble and export pfms for fimo, if not already loaded and running
     # pfms <- query(MotifDb, c("Athaliana", "jaspar2018"))
     # export(pfms, "~/github/fimoService/pfms/athaliana.all.meme", "meme")
   start <- tss.orf.1 - 499
   end   <- tss.orf.1 + 2000
   sequence <- as.character(getSeq(BSgenome.Athaliana.TAIR.TAIR9, "Chr3", start, end))
   nchar(sequence)
   tbl.fimo <- requestMatch(fc, list(frd3=sequence), pvalThreshold=fimo.pval.threshold)
   dim(tbl.fimo)
   tbl.fimo$start <- tbl.fimo$start + start
   tbl.fimo$stop <- tbl.fimo$stop + start
   tbl.fimo$chrom <- chromosome
   tbl.bed.1 <- tbl.fimo[, c("chrom", "start", "stop", "score", "motif")]
   dim(tbl.bed.1)
   tbl.xtab <- as.data.frame(table(tbl.fimo$motif))
   fivenum(tbl.xtab$Freq)
   candidate.motifs <- unique(tbl.fimo$motif)
   length(candidate.motifs) # 373
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
   length(candidate.tfs.1)
   } # tfs

if(!exists("candidate.tfs.2")){
     # assemble and export pfms for fimo, if not already loaded and running
     # pfms <- query(MotifDb, c("Athaliana", "jaspar2018"))
     # export(pfms, "~/github/fimoService/pfms/athaliana.all.meme", "meme")
   start <- tss.orf.2 - 499
   end   <- tss.orf.2 + 2000
   sequence <- as.character(getSeq(BSgenome.Athaliana.TAIR.TAIR9, "Chr3", start, end))
   nchar(sequence)
   tbl.fimo <- requestMatch(fc, list(frd3=sequence), pvalThreshold=fimo.pval.threshold)
   dim(tbl.fimo)
   tbl.fimo$start <- tbl.fimo$start + start
   tbl.fimo$stop <- tbl.fimo$stop + start
   tbl.fimo$chrom <- chromosome
   tbl.bed.2 <- tbl.fimo[, c("chrom", "start", "stop", "score", "motif")]
   dim(tbl.bed.2)
   tbl.xtab <- as.data.frame(table(tbl.fimo$motif))
   fivenum(tbl.xtab$Freq)
   candidate.motifs <- unique(tbl.fimo$motif)
   tfs.2 <- as.data.frame(mcols(MotifDb[candidate.motifs]))$geneSymbol
   tfs.orfs <- tbl.xref$orfs[match(tfs.2, tbl.xref$geneSymbol)]
   deleters <- which(is.na(tfs.orfs))
   printf("%d/%d tf.orfs are NA", length(deleters), length(tfs.orfs))
   if(length(deleters) > 0)
      tfs.orfs <- tfs.orfs[-deleters]
   candidate.tfs.2 <- tfs.orfs
   length(candidate.tfs.2)
   } # tfs

if(!exists("tbl.model.1")) {
   length(intersect(candidate.tfs.1, rownames(mtx))) # 294
   solver <- EnsembleSolver(mtx, targetGene, candidate.tfs.1)
   tbl.model <- run(solver)
   tbl.model$geneSymbol <- tbl.xref[match(tbl.model$gene, tbl.xref$orfs), "geneSymbol"]
   tbl.model$SYMBOL <- tbl.xref[match(tbl.model$gene, tbl.xref$orfs), "SYMBOL"]
   tbl.model <- tbl.model[order(tbl.model$pcaMax, decreasing=TRUE),]
   tbl.model <- roundNumericColumnsInDataframe(tbl.model, 3, "lassoPValue")
   tbl.model <- subset(tbl.model, abs(pearsonCoeff) > 0.4)
   if(nrow(tbl.model) > 10)
      tbl.model <- tbl.model[1:10,]
   tbl.model.1 <- tbl.model
   save(tbl.model.1, file="../shinyApp/data/tbl.model.frd3-transcript1.RData")
   }  # tbl.model

if(!exists("tbl.model.2")) {
   solver <- EnsembleSolver(mtx, targetGene, candidate.tfs.2)
   tbl.model <- run(solver)
   tbl.model$symbol <- tbl.xref[match(tbl.model$gene, tbl.xref$orfs), "geneSymbol"]
   tbl.model <- tbl.model[order(tbl.model$pcaMax, decreasing=TRUE),]
   tbl.model <- roundNumericColumnsInDataframe(tbl.model, 3, "lassoPValue")
   tbl.model <- subset(tbl.model, abs(pearsonCoeff) > 0.4)
   if(nrow(tbl.model) > 10)
      tbl.model <- tbl.model[1:10,]
   tbl.model.2 <- tbl.model
   save(tbl.model.2, file="../shinyApp/data/tbl.model.frd3-transcript2.RData")
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

if(!exists("tbl.bindingSites.1")){
   tbl.bindingSites.1 <- add.orfs.geneSymbol.to.fimoResults(tbl.bed.1)
   tbl.bindingSites.1 <-  subset(tbl.bindingSites.1, orf %in% tbl.model.1$gene)
   dim(tbl.bindingSites.1)
   tfs.unique <- sort(unique(tbl.bindingSites.1$geneSymbol))
   for(tf in tfs.unique){
      tbl.tfbs <- subset(tbl.bindingSites.1, geneSymbol==tf)
      track <- DataFrameAnnotationTrack(tf, tbl.tfbs, color="red", trackHeight=25)
      displayTrack(igv, track)
      } # for tf
   save(tbl.bindingSites.1, file="../shinyApp/data/tbl.bindingSites.frd3-transcript1.RData")
   } # tbl.bindingSites.1


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

