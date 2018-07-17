# findAllBindingSites.R: all 15 curated pfms, over the 318kb region
#------------------------------------------------------------------------------------------------------------------------
library(FimoClient)
library(trena)   # for motifMatcher, which obtains sequence
library(RUnit)
library(igvR)
#------------------------------------------------------------------------------------------------------------------------
FIMO_HOST <- "localhost"
FIMO_PORT <- 5558
if(!exists("fimo")){
   fimo <- FimoClient(FIMO_HOST, FIMO_PORT, quiet=FALSE)
   }
#------------------------------------------------------------------------------------------------------------------------
if(!exists("tbl.enhancers")){
   load("../shinyApp/data/tbl.enhancers.RData")
   gr.enhancers <- GRanges(tbl.enhancers)
   }

if(!exists("tbl.model")){
   load("../shinyApp/data/tbl.model.RData")
   tfs.all <- tbl.model$tf.hgnc[1:10]
   }
#------------------------------------------------------------------------------------------------------------------------
targetGene <- "INPP5D"
chromosome <- tbl.enhancers$chr[1]
loc.min <- min(tbl.enhancers$start) - 2000
loc.max <-  max(tbl.enhancers$end) + 2000

pfm.names <- c(sp100="Hsapiens-SwissRegulon-SP100.SwissRegulon",
               spi1="Hsapiens-jaspar2018-SPI1-MA0080.4",
               tal1="Hsapiens-HOCOMOCOv10-TAL1_HUMAN.H10MO.S",
               nfatc2="Mmusculus;Rnorvegicus;Hsapiens-jaspar2018-NFATC2-MA0152.1",
               fli1="Hsapiens-jaspar2018-FLI1-MA0475.2")
pfm.map <- toupper(names(pfm.names))
names(pfm.map) <- as.character(pfm.names)
pfm.map <- as.list(pfm.map)

#------------------------------------------------------------------------------------------------------------------------
if(!exists("igv")){
   igv <- igvR()
   Sys.sleep(5)
   setGenome(igv, "hg38")
   Sys.sleep(5)
   showGenomicRegion(igv, targetGene)
   showGenomicRegion(igv, sprintf("%s:%d-%d", chromosome, loc.min, loc.max))
   }

#------------------------------------------------------------------------------------------------------------------------
showRegion <- function(tag)
{
   showGenomicRegion(igv, switch(tag,
                       "enhancers" = "chr2:232,850,784-233,511,033",
                       "gene" = "INPP5D",
                       "snps" = "chr2:232,999,213-233,163,733"
                        ))
} # showRegion
#------------------------------------------------------------------------------------------------------------------------
displayEnhancers <- function()
{
   track <- DataFrameAnnotationTrack("enhancers", tbl.enhancers, color="blue", trackHeight=30)
   displayTrack(igv, track)

} # displayEnhancers
#------------------------------------------------------------------------------------------------------------------------
displaySnps <- function()
{
   track <- DataFrameAnnotationTrack("rsid", tbl.enhancer.snps[, 1:4], color="red", trackHeight=20)
   displayTrack(igv, track)

   track <- DataFrameQuantitativeTrack("snps", tbl.enhancer.snps[, c(1:3, 8)], color="red", trackHeight=60)
   displayTrack(igv, track)

   stopifnot(exists("tbl.breaks"))
   snps.in.model <- rownames(tbl.breaks)
   deleters <- grep("\\.[1-9]$", snps.in.model)
   if(length(deleters) > 0)
      snps.in.model <- snps.in.model[-deleters]

   tbl.snps.in.model <- subset(tbl.snp, rsid %in% snps.in.model)

   track <- DataFrameAnnotationTrack("snp in model", tbl.snps.in.model, color="maroon", trackHeight=20)
   displayTrack(igv, track)

} # displaySnps
#------------------------------------------------------------------------------------------------------------------------
run <- function()
{
   mm <- MotifMatcher("hg38", list(), quiet=TRUE)
   tbl.regionMax <- data.frame(chrom=chromosome, start=loc.min, end=loc.max, stringsAsFactors=FALSE)
   tbl.seq <- getSequence(mm, tbl.regionMax)
   sequence <- tbl.seq$seq  # 318k
   nchar(sequence)  # 582k
   stopifnot(nchar(sequence) == (1 + tbl.regionMax$end - tbl.regionMax$start))
   tbl.fimo.raw <- requestMatch(fimo, list(big=sequence))
   dim(tbl.fimo.raw) # [1] 325k x 9
   save(tbl.fimo.raw, file="tbl.fimo.raw.325k.regionRData")
   tbl.fimo <- tbl.fimo.raw
   tbl.fimo$start <- tbl.fimo$start + tbl.regionMax$start
   tbl.fimo$stop <-  tbl.fimo$stop + tbl.regionMax$start
   tbl.fimo$pScore <- -log10(tbl.fimo$p.value)
   dim(tbl.fimo)
   dim(subset(tbl.fimo, pScore >= 3))
   tbl.fimo <- subset(tbl.fimo, pScore >= 3)
   tbl.fimo$chrom <- chromosome

   tfs <- as.character(pfm.map[tbl.fimo$motif])
   checkEquals(length(tfs), nrow(tbl.fimo))
   tbl.fimo$tf <- tfs
   dim(tbl.fimo)   # 6089 x 12
   tbl.fimo <- tbl.fimo[, c("chrom", "start", "stop", "strand", "motif", "tf", "score", "pScore", "p.value", "q.value", "matched.sequence")]
   colnames(tbl.fimo) <- c("chrom", "motifStart", "motifEnd", "strand", "motifName", "tf", "motifScore",
                           "motif.pScore", "motif.pVal", "motif.qVal", "sequence")
   tbl.fimo <- tbl.fimo[order(tbl.fimo$motifStart, decreasing=FALSE),]
   save(tbl.fimo, file="../shinyApp/data/tbl.fimo.5tfs.6089bindingSites.RData")

   tfs <- unique(tbl.fimo$tf)
   length(tfs)
   table(tbl.fimo$tf)
   fivenum(tbl.fimo$motif.pScore)  # [1] 0.000000 0.000224 0.000486 0.000721 0.001000

   for(current.tf in tfs){
      tbl.tf <- subset(tbl.fimo, tf==current.tf & motif.pScore >= 4)
      dim(tbl.tf)
      tbl.tf <- tbl.tf[, c("chrom", "motifStart", "motifEnd", "motif.pScore")]
      colnames(tbl.tf) <- c("chrom", "start", "end", "score")
      fivenum(tbl.tf$score)  # [1] 3.000000 3.131063 3.336299 3.671620 6.290730
      track <- DataFrameQuantitativeTrack(current.tf, tbl.tf, color="darkBlue", min=2.5, max=7)
      displayTrack(igv, track)
      } # for current.tf

} # run
#------------------------------------------------------------------------------------------------------------------------
