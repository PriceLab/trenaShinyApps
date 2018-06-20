# findAllBindingSites.R: all ~50 curated pfms, over the near and far enhancer regions (gaps filled)
#------------------------------------------------------------------------------------------------------------------------
proximal.enhancer.region <- "chr2:112,719,115-112,905,047"  # 185 kb
distant.enhancer.region  <- "chr2:112,022,627-112,067,553"  #  44 kb
load("../data/enhancers.RData")
#------------------------------------------------------------------------------------------------------------------------
library(FimoClient)
library(trenaSGM)   # for motifMatcher, which obtains sequence
library(RUnit)
library(igvR)
#------------------------------------------------------------------------------------------------------------------------
FIMO_HOST <- "localhost"
FIMO_PORT <- 5558
if(!exists("fimo")){
   fimo <- FimoClient(FIMO_HOST, FIMO_PORT, quiet=FALSE)
   }
#------------------------------------------------------------------------------------------------------------------------
genome <- "hg38"
targetGene <- "CKAP2L"
chromosome <- "chr2"
tss <- 112764686
snp.rsid <- "rs7594852"
snp.bp <- 112764177
tbl.snp <- data.frame(chrom=chromosome, start=snp.bp, end=snp.bp, rsid=snp.rsid, stringsAsFactors=FALSE)
save(tbl.snp, file="../data/tbl.snp.RData")
if(!exists("sgm")){
   sgm <- trenaSGM(genome, targetGene)
   }
#------------------------------------------------------------------------------------------------------------------------
if(!exists("igv")){
   igv <- igvR()
   Sys.sleep(5)
   setGenome(igv, genome)
   Sys.sleep(5)
   showGenomicRegion(igv, targetGene)
   showGenomicRegion(igv, proximal.enhancer.region)
   track <- DataFrameAnnotationTrack("enhancers", tbl.enhancers, "green")
   displayTrack(igv, track)
   track <- DataFrameAnnotationTrack("snp", tbl.snp, "red")
   displayTrack(igv, track)
   }
#------------------------------------------------------------------------------------------------------------------------
load("../data/models.RData")
tbl.summary <- summarizeModels(sgm, "rfScore", 10, models)
tfs.all <- rownames(subset(tbl.summary, observed > 1 & rank.sum < 50))
length(tfs.all) # 17
# "NFE2L2" "FOXM1"  "MXI1"   "ETV6"   "ETV5"   "CEBPA"  "MEF2D"  "SP2"    "ZFHX3"  "MYBL2"  "POU3F3" "POU3F2"
# "NR6A1"  "STAT4"  "ETS1"   "GRHL1"  "HBP1"   "RFX7"

#  "PLEK"   "IKZF1"  "LYL1"   "SPI1"   "VAV1"   "IRF5"   "IRF8"   "TFEC"   "BCL6"   "CEBPA"  "ZBTB18"
#  "NFATC2" "EGR4"   "FOXP1"  "STAT1"  "ELK3"   "TAL1"
#  setdiff(tfs, tfs.with.single.motifs)  #  "PLEK" "LYL1" "VAV1"

mm <- MotifMatcher("hg38", list(), quiet=TRUE)

tbl.proximalEnhancerRegion <- data.frame(chrom="chr2", start=112719115, end=112905047)
tbl.region <- tbl.proximalEnhancerRegion

tbl.seq <- getSequence(mm, tbl.region)
sequence <- tbl.seq$seq  # 318k
stopifnot(nchar(sequence) == (1 + tbl.region$end - tbl.region$start))
tbl.fimo.raw <- requestMatch(fimo, list(big=sequence))
dim(tbl.fimo.raw) # [1] 2228048       9
save(tbl.fimo.raw, file="tbl.fimo.raw.proximalEnhancerRegion.Data")
tbl.fimo <- tbl.fimo.raw
tbl.fimo$start <- tbl.fimo$start + tbl.region$start
tbl.fimo$stop <-  tbl.fimo$stop + tbl.region$start
dim(tbl.fimo)
tbl.fimo <- subset(tbl.fimo, p.value <= 0.0001)
dim(tbl.fimo)   # 29406 9
tbl.fimo$chrom <- "chr2"
   # add the tf column
tfs <- mcols(MotifDb[tbl.fimo$motif])$geneSymbol
tbl.fimo$tf <- tfs

dim(tbl.fimo)   # 29406 x 11
tbl.fimo <- tbl.fimo[, c("chrom", "start", "stop", "strand", "motif", "tf", "score", "p.value", "q.value", "matched.sequence")]
colnames(tbl.fimo) <- c("chrom", "motifStart", "motifEnd", "strand", "motifName", "tf", "motifScore", "motif.pVal", "motif.qVal", "sequence")
tbl.fimo <- tbl.fimo[order(tbl.fimo$motifStart, decreasing=FALSE),]
save(tbl.fimo, file="../data/tbl.fimo.55tfs.29406.bindingSites.RData")

tfs <- unique(tbl.fimo$tf)
length(tfs)
table(tbl.fimo$tf)
fivenum(tbl.fimo$motif.pVal)  # [1] 0.000000 0.000224 0.000486 0.000721 0.001000

top.tfs <- c("NFE2L2","POU3F3","FOXM1","ETV6","CEBPA","SP2","POU3F2","MEF2D","ZFHX3","NR6A1","RFX7",
             "ETS1","TFCP2L1","E2F1","HBP1","SOX12","GLIS2","ELF3","SREBF1","CDC5L")

tbl.fimo <- subset(tbl.fimo, tf %in% top.tfs)
dim(tbl.fimo)
save(tbl.fimo, file="../data/tbl.fimo.20tfs.1115.bindingSites.RData")
library(RColorBrewer)
#display.brewer.all()
colors <- brewer.pal(12, "Paired")
colors <- rep(colors,3)
stopifnot(length(colors) > length(top.tfs))

i <- 0
for(current.tf in top.tfs){
   i <- i + 1
   color <- colors[1]
   tbl.tf <- subset(tbl.fimo, tf==current.tf & motif.pVal <= 0.0001)
   dim(tbl.tf)
   tbl.tf <- tbl.tf[, c("chrom", "motifStart", "motifEnd", "motif.pVal")]
   colnames(tbl.tf) <- c("chrom", "start", "end", "score")
   tbl.tf$score <- -log10(tbl.tf$score)
   fivenum(tbl.tf$score)  # [1] 3.000000 3.131063 3.336299 3.671620 6.290730
   track <- DataFrameQuantitativeTrack(current.tf, tbl.tf, color=color, min=2.5, max=7)
   displayTrack(igv, track)
   } # for current.tf
