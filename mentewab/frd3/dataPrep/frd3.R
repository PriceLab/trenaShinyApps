library(trena)
library(trenaSGM)
library(igvR)
library(GenomicRanges)
library(MotifDb)
library(FimoClient)
library(RUnit)
genomeName <- "BSgenome.Athaliana.TAIR.TAIR9"
library(BSgenome.Athaliana.TAIR.TAIR9)
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
if(!exists("tfs")){
   pfms <- query(MotifDb, c("Athaliana", "jaspar2018"))
   export(pfms, "~/github/fimoService/pfms/athaliana.all.meme", "meme")
   start <- tss.orf.2 - 500
   end   <- tss.orf.2 + 2000
   sequence <- as.character(getSeq(BSgenome.Athaliana.TAIR.TAIR9, "Chr3", start, end))
   nchar(sequence)
   tbl.fimo <- requestMatch(fc, list(frd3=sequence), pvalThreshold=0.0001)
   dim(tbl.fimo)
   tbl.fimo$start <- tbl.fimo$start + start
   tbl.fimo$stop <- tbl.fimo$stop + start
   tbl.fimo$chrom <- chromosome
   tbl.bed <- tbl.fimo[, c("chrom", "start", "stop", "score", "motif")]
   tbl.xtab <- as.data.frame(table(tbl.fimo$motif))
   fivenum(tbl.xtab$Freq)
   candidate.motifs <- unique(tbl.fimo$motif)
   tfs <- as.data.frame(mcols(MotifDb[candidate.motifs]))$geneSymbol
   tfs.orfs <- tbl.xref$orfs[match(tfs, tbl.xref$geneSymbol)]
   deleters <- which(is.na(tfs.orfs))
   if(length(deleters) > 0)
      tfs.orfs <- tfs.orfs[-deleters]
   candidate.tfs <- tfs.orfs
   length(candidate.tfs)
   } # tfs

if(!exists("tbl.model")) {
   solver <- EnsembleSolver(mtx, targetGene, candidate.tfs)
   tbl.model <- run(solver)
   tbl.model <- tbl.model[order(tbl.model$rfScore, decreasing=TRUE),]
   tbl.model$symbol <- tbl.xref[match(tbl.model$gene, tbl.xref$orfs), "geneSymbol"]
   }  # tbl.model

if(!exists("tbl.bindingSites")){
   tbl.bindingSites <- subset(tbl.bed, motif %in% tbl.xref[match(tbl.model$gene, tbl.xref$orfs), "motif"])
   track <- DataFrameAnnotationTrack("model.bs", tbl.bindingSites, color="red")
   displayTrack(igv, track)
   } # tbl.bindingSites


#------------------------------------------------------------------------------------------------------------------------

