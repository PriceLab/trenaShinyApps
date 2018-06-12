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
if(!exists("igv")){
   igv <- igvR()
   Sys.sleep(5)
   setGenome(igv, genome)
   Sys.sleep(5)
   showGenomicRegion(igv, targetGene)
   showGenomicRegion(igv, "chr6:41,049,342-41,367,926")
   }
#------------------------------------------------------------------------------------------------------------------------

load("../shinyApp/trem2-models.RData")
tbl.summary <- all.models$summary$model
tfs.all <- subset(tbl.summary, observed > 1 & rank.sum < 50)$TF
length(tfs.all) # 17
#  "PLEK"   "IKZF1"  "LYL1"   "SPI1"   "VAV1"   "IRF5"   "IRF8"   "TFEC"   "BCL6"   "CEBPA"  "ZBTB18"
#  "NFATC2" "EGR4"   "FOXP1"  "STAT1"  "ELK3"   "TAL1"
#  setdiff(tfs, tfs.with.single.motifs)  #  "PLEK" "LYL1" "VAV1"

mm <- MotifMatcher("hg38", list(), quiet=TRUE)

tbl.regionMax <- data.frame(chrom="chr6", start=41049342, end=41367926)

tbl.seq <- getSequence(mm, tbl.regionMax)
sequence <- tbl.seq$seq  # 318k
stopifnot(nchar(sequence) == (1 + tbl.regionMax$end - tbl.regionMax$start))
tbl.fimo.raw <- requestMatch(fimo, list(big=sequence))
dim(tbl.fimo.raw) # [1] 1043063       9

save(tbl.fimo.raw, file="tbl.fimo.raw.318k.regionRData")
tbl.fimo <- tbl.fimo.raw
tbl.fimo$start <- tbl.fimo$start + tbl.regionMax$start
tbl.fimo$stop <-  tbl.fimo$stop + tbl.regionMax$start
dim(tbl.fimo)
tbl.fimo <- subset(tbl.fimo, p.value <= 0.001)
dim(tbl.fimo)
tbl.fimo$chrom <- "chr6"
   # add the tf column

motif.names <- c("Hsapiens-HOCOMOCOv10-STAT1_HUMAN.H10MO.A",
                 "Hsapiens-jaspar2018-SPI1-MA0080.1",
                 "Hsapiens-HOCOMOCOv10-BCL6B_HUMAN.H10MO.D",
                 "Hsapiens-HOCOMOCOv10-CEBPA_HUMAN.H10MO.A",
                 "Hsapiens-HOCOMOCOv10-EGR4_HUMAN.H10MO.D",
                 "Hsapiens-SwissRegulon-ELK3.SwissRegulon",
                 "Hsapiens-jaspar2018-FOXP1-MA0481.2",
                 "Hsapiens-HOCOMOCOv10-IRF5_HUMAN.H10MO.D",
                 "Hsapiens-HOCOMOCOv10-IRF8_HUMAN.H10MO.D",
                 "Hsapiens-HOCOMOCOv10-TAL1_HUMAN.H10MO.S",
                 "Hsapiens-SwissRegulon-BCL6.SwissRegulon",
                 "Hsapiens-SwissRegulon-IKZF1.SwissRegulon",
                 "Hsapiens-SwissRegulon-NFATC2.SwissRegulon",
                 "Hsapiens-SwissRegulon-TFEC.SwissRegulon",
                 "Hsapiens-SwissRegulon-ZBTB18.SwissRegulon")

pfms.chosen <- MotifDb[motif.names]

motif2tf.map <- mcols(pfms.chosen)$geneSymbol
names(motif2tf.map) <- motif.names
save(motif2tf.map, file="motif2tf.map.RData")
tfs <- as.character(motif2tf.map[tbl.fimo$motif])
checkEquals(length(tfs), nrow(tbl.fimo))
tbl.fimo$tf <- tfs
dim(tbl.fimo)   # 16435 x 11
tbl.fimo <- tbl.fimo[, c("chrom", "start", "stop", "strand", "motif", "tf", "score", "p.value", "q.value", "matched.sequence")]
colnames(tbl.fimo) <- c("chrom", "motifStart", "motifEnd", "strand", "motifName", "tf", "motifScore", "motif.pVal", "motif.qVal", "sequence")
tbl.fimo <- tbl.fimo[order(tbl.fimo$motifStart, decreasing=FALSE),]
save(tbl.fimo, file="tbl.fimo.15tfs.16435bindingSites.RData")

tfs <- unique(tbl.fimo$tf)
length(tfs)
table(tbl.fimo$tf)
fivenum(tbl.fimo$motif.pVal)  # [1] 0.000000 0.000224 0.000486 0.000721 0.001000

for(current.tf in tfs[1]){
   tbl.tf <- subset(tbl.fimo, tf==current.tf & motif.pVal <= 0.0001)
   dim(tbl.tf)
   tbl.tf <- tbl.tf[, c("chrom", "motifStart", "motifEnd", "motif.pVal")]
   colnames(tbl.tf) <- c("chrom", "start", "end", "score")
   tbl.tf$score <- -log10(tbl.tf$score)
   fivenum(tbl.tf$score)  # [1] 3.000000 3.131063 3.336299 3.671620 6.290730
   track <- DataFrameQuantitativeTrack(current.tf, tbl.tf, color="darkBlue", min=2.5, max=7)
   displayTrack(igv, track)
   } # for current.tf
