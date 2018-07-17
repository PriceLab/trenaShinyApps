library(GenomicRanges)
library(MotifDb)
library(FimoClient)
library(trena)
library(motifStack)
library(motifbreakR)
library(SNPlocs.Hsapiens.dbSNP150.GRCh38)
library(BSgenome.Hsapiens.UCSC.hg38)
library(igvR)
library(GenomicRanges)
source("../../../utils/roundNumericColumnsInDataframe.R")
#------------------------------------------------------------------------------------------------------------------------
# FIMO_HOST <- "localhost"
# FIMO_PORT <- 5558
# fc <- FimoClient(FIMO_HOST, FIMO_PORT, quiet=FALSE)
#------------------------------------------------------------------------------------------------------------------------
targetGene <- "INPP5D"

# Inositol polyphosphate-5-phosphatase (INPP5D) was reported to be
# associated with Alzheimer's disease (AD) through modulating the
# inflammatory process and immune response.
#------------------------------------------------------------------------------------------------------------------------
if(!exists("tbl.enhancers")){
   source("~/s/data/public/GeneHancer/v4.7/getAllEnhancers.R")
   tbl.enhancers <- getAllEnhancers(targetGene)
   colnames(tbl.enhancers) <- c("chr", "start", "end", "type", "score")
   gr.enhancers <- GRanges(tbl.enhancers)
   }

#------------------------------------------------------------------------------------------------------------------------
loc.min <- min(tbl.enhancers$start) - 2000
loc.max <-  max(tbl.enhancers$end) + 2000
chromosome <- tbl.enhancers$chr[1]
#------------------------------------------------------------------------------------------------------------------------
if(!exists("tbl.model")){
   load("../../survey/tbl.models.all.RData")
   tbl.model <- subset(tbl.models, targetGene.hgnc==targetGene)
   tbl.model <- tbl.model[order(tbl.model$rfScore, decreasing=TRUE),][1:10,]
   tbl.model <- tbl.model[, c(2,5:11,12,13)]
   tbl.model <- roundNumericColumnsInDataframe(tbl.model, 3, pvalColumnNames="lassoPValue")
   save(tbl.model, file="../shinyApp/data/tbl.model.RData")
   }
#------------------------------------------------------------------------------------------------------------------------
if(!exists("tbl.gwas")){
   load("~/s/work/priceLab/AD/tbl.gwas.level_1.RData") # tbl.gwas
   tbl.snp <- tbl.gwas[, c("CHR", "BP", "BP", "SNP", "Effect_allele", "Non_Effect_allele", "P")]
   colnames(tbl.snp) <- c("chrom", "start", "end", "rsid", "variant", "reference", "pval")
   tbl.snp$pScore <- -log10(tbl.snp$pval)
   gr.snps <- GRanges(tbl.snp)
   }
#------------------------------------------------------------------------------------------------------------------------
if(!exists("enhancer.snps")){
   tbl.ov <- as.data.frame(findOverlaps(gr.snps, gr.enhancers))
   colnames(tbl.ov) <- c("snp", "enhancer")
   tbl.enhancer.snps <- cbind(tbl.snp[tbl.ov$snp,], tbl.enhancers[tbl.ov$enhancer,])
   dim(tbl.enhancer.snps) # 65 13
   fivenum(tbl.enhancer.snps$pScore) # 1.302422 2.210139 3.193956 4.349514 8.471598
   fivenum(tbl.gwas$P) # [1] 1.000e-150  1.005e-02  2.270e-02  3.620e-02  5.000e-02
   enhancer.snps <- unique(subset(tbl.enhancer.snps, pScore >= 1.30103)$rsid)  # this is everything...
   length(enhancer.snps)
   }
#------------------------------------------------------------------------------------------------------------------------
pfm.names <- c(sp100="Hsapiens-SwissRegulon-SP100.SwissRegulon",
               spi1="Hsapiens-jaspar2018-SPI1-MA0080.4",
               tal1="Hsapiens-HOCOMOCOv10-TAL1_HUMAN.H10MO.S",
               nfatc2="Mmusculus;Rnorvegicus;Hsapiens-jaspar2018-NFATC2-MA0152.1",
               fli1="Hsapiens-jaspar2018-FLI1-MA0475.2")
               #lyl1="Hsapiens-jaspar2018-NHLH1-MA0048.2")   # inferred loosely from TFclass

#------------------------------------------------------------------------------------------------------------------------
if(!exists("igv")){
   igv <- igvR()
   Sys.sleep(3)
   setGenome(igv, "hg38")
   Sys.sleep(3)
   showGenomicRegion(igv, targetGene)
   Sys.sleep(5)

   }
#------------------------------------------------------------------------------------------------------------------------
selectMotif <- function()
{
   pfms <- query(MotifDb, "LYL1")
   length(pfms)
   motifs <- geneToMotif(MotifDb, "LYL1", "tfclass")$motif
   pfms <- query(MotifDb, "sapiens", motifs)
   length(pfms)
   motifStack(lapply(names(pfms[11:15]), function(mName) new("pfm", pfms[[mName]], name=mName)))

} # selectMotif
#------------------------------------------------------------------------------------------------------------------------
showRegion <- function(tag)
{
   showGenomicRegion(igv, switch(tag,
                       "enhancers" = "chr2:232,911,286-233,511,033",
                       "gene" = "INPP5D",
                       "snps" = "chr2:232,999,213-233,163,733"
                        ))

} # showRegion
#------------------------------------------------------------------------------------------------------------------------
displayEnhancers <- function()
{
   track <- DataFrameAnnotationTrack("enhancers", tbl.enhancers, color="blue", trackHeight=30)
   displayTrack(igv, track)

} # displayPromoters
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

} # displayHotSnps
#------------------------------------------------------------------------------------------------------------------------
mbreak <- function()
{
   snps.gr <- snps.from.rsid(rsid=enhancer.snps,
                             dbSNP=SNPlocs.Hsapiens.dbSNP150.GRCh38,
                             search.genome=BSgenome.Hsapiens.UCSC.hg38)
   motifs <- MotifDb[as.character(pfm.names)]
   length(snps.gr)
   length(motifs)
   results <- motifbreakR(snpList = snps.gr,
                          filterp = TRUE,
                          threshold = 0.85,
                          pwmList = motifs,
                          show.neutral=TRUE,
                          method = c("ic", "log", "notrans")[2],
                          bkg = c(A=0.25, C=0.25, G=0.25, T=0.25),
                          BPPARAM = BiocParallel::bpparam(),
                          verbose=TRUE)
   results.filtered <- results[results$pctRef > 0.8]
   results.filtered <- results.filtered[results.filtered$effect != "neut"]
   results.filtered.pvals <- calculatePvalue(results.filtered)
   tbl.breaks <- as.data.frame(mcols(results.filtered.pvals))
   tbl.breaksBed <- tbl.breaks[, c("snpPos", "snpPos")]
   colnames(tbl.breaksBed) <- c("start", "end")
   tbl.breaksBed$chrom <- chromosome
   tbl.breaksBed$rsid <- rownames(tbl.breaks)
   rownames(tbl.breaksBed) <- NULL
   coi <- c("chrom", "start", "end", "rsid")
   tbl.breaksBed <- tbl.breaksBed[, coi]
   rsid.with.prefix <- sprintf("tfbs-snp:%s", tbl.breaksBed$rsid)
   tbl.breaksBed$rsid <- rsid.with.prefix
   save(tbl.breaksBed, file="../shinyApp/data/tbl.breaksBed.RData")
   tbl.breaks$pctRef <- round(tbl.breaks$pctRef, 2)
   tbl.breaks$pctAlt <- round(tbl.breaks$pctAlt, 2)
   save(tbl.breaks, file="../shinyApp/data/tbl.breaks.RData")

   save(results, results.filtered.pvals, tbl.breaks, file="motifbreakR.results.RData")
   rsid <-  "rs13425849"
   for(rsid in rownames(tbl.breaks[1,])){
      filename <- sprintf("../shinyApp/png/%s.png", rsid)
      png(filename, width=800)
      plotMB(results.filtered.pvals, rsid)
      dev.off()
      printf("created %s", filename)
      }

} # mbreak
#------------------------------------------------------------------------------------------------------------------------
