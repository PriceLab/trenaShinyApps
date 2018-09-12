# findAllBindingSites.R: all 15 curated pfms, over the 318kb region
#------------------------------------------------------------------------------------------------------------------------
library(FimoClient)
library(trena)   # for motifMatcher, which obtains sequence
library(RUnit)
library(igvR)
library(GenomicRanges)
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
   tbl.fimo <- subset(tbl.fimo,  pScore >= -log10(0.01))
   dim(tbl.fimo)   # 145k
   tbl.fimo$chrom <- chromosome

   tfs <- as.character(pfm.map[tbl.fimo$motif])
   checkEquals(length(tfs), nrow(tbl.fimo))
   tbl.fimo$tf <- tfs
   dim(tbl.fimo)   # 145k x 12
   tbl.fimo <- tbl.fimo[, c("chrom", "start", "stop", "strand", "motif", "tf", "score", "pScore", "p.value", "q.value", "matched.sequence")]
   tbl.fimo <- tbl.fimo[order(tbl.fimo$start, decreasing=FALSE),]
   gr.fimo <- GRanges(tbl.fimo)   # save this for below
   colnames(tbl.fimo) <- c("chrom", "motifStart", "motifEnd", "strand", "motifName", "tf", "motifScore",
                           "motif.pScore", "motif.pVal", "motif.qVal", "sequence")
   tbl.ov <- as.data.frame(findOverlaps(gr.enhancers, gr.fimo))
   colnames(tbl.ov) <- c("enhancers", "tfbs")
   tbl.fimoE <- tbl.fimo[unique(tbl.ov$tfbs),]
   dim(tbl.fimoE)
   tbl.fimo <- tbl.fimoE
   save(tbl.fimo, file="../shinyApp/data/tbl.fimo.5tfs.bindingSites.RData")

   tfs <- unique(tbl.fimo$tf)
   length(tfs)
   table(tbl.fimo$tf)
   fivenum(tbl.fimo$motif.pScore)  # [1] 0.000000 0.000224 0.000486 0.000721 0.001000

   track <- DataFrameAnnotationTrack("enhancers", tbl.enhancers[, c(1:3, 4)], color="purple", trackHeight=50)
   displayTrack(igv, track)

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
# chr2:233,132,931-233,133,659
# nfatc2 motif breaking snp rs10933435
# intersects a FLI1 binding site, and an overlapping somewhat larger (upstream) SPI1 site as well
# but no binding site for nfatc2 shows up.
corysBug <- function()
{

   snp.loc <- 233133436
   subset(tbl.fimo, motifStart <= snp.loc & motifEnd >= snp.loc)
   tbl.x <- subset(tbl.fimo, start <= snp.loc & stop >= snp.loc)

   # nfatc2 has pScore of 2.89, just missed being included in tbl.fimo used by the shinyapp.
   # remedy:  add slider on pScore, be more inclusive in the delivered tbl.fimo

} # corysBug
#------------------------------------------------------------------------------------------------------------------------
