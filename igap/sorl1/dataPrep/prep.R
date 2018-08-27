# ~/github/trenaShinyApps/igap/sorl1/prep.R
# extract gene model using top 15 tfs
# load gwas variants
# use sorl1-specific fimo to identify all possible binding sites
# intersect those binding sites with the gwas variants, invoke motifbreakR
#------------------------------------------------------------------------------------------------------------------------
source("../../../utils/roundNumericColumnsInDataframe.R")
#------------------------------------------------------------------------------------------------------------------------
library(GenomicRanges)
library(FimoClient)
library(BSgenome.Hsapiens.UCSC.hg38)
library(SNPlocs.Hsapiens.dbSNP150.GRCh38)
library(igvR)
library(motifbreakR)
#------------------------------------------------------------------------------------------------------------------------
genome <- "hg38"
targetGene <- "SORL1"
chromosome <- "chr11"
tss <- 121452363
#------------------------------------------------------------------------------------------------------------------------
FIMO_HOST <- "localhost"
FIMO_PORT <- 5558
if(!exists("fimo")){
   fimo <- FimoClient(FIMO_HOST, FIMO_PORT, quiet=FALSE)
   }
#------------------------------------------------------------------------------------------------------------------------
findFimoHits <- function(tbl.regions, pval=10e-4)
{
   seqs <- as.list(as.character(with(tbl.regions, getSeq(BSgenome.Hsapiens.UCSC.hg38, chrom, start, end))))
   x <- requestMatch(fimo, seqs, pval)
   invisible(x)

} # findFimoHits
#------------------------------------------------------------------------------------------------------------------------
test_findFimoHits <- function()
{
   tbl.bigHit <- subset(tbl.regions, combinedScore > 500)
   tbl.x <- findFimoHits(tbl.bigHit, pval=10e-6)
   tbl.moreHits <- subset(tbl.regions, combinedScore > 30)
   tbl.xfindFimoHits(tbl.bigHit)

} # test_findFimoHits
#------------------------------------------------------------------------------------------------------------------------
if(!exists(tbl.enhancers)){
   load("~/github/genomescale_scripts/geneHancer/v4.7/tbl.enhancers.allGenes.RData")
   dim(tbl.enhancers) #  716477      6
   tbl.enhancers <- subset(tbl.enhancers, geneSymbol=="SORL1")
   save(tbl.enhancers, file="tbl.enhancers.RData")
   region.start <- min(tbl.enhancers$start) - 1000
   region.end   <- max(tbl.enhancers$start) + 1000
   }

if(!exists("tbl.model")){
   load("~/github/trenaShinyApps/igap/survey/tbl.models.all.RData")
   tbl.model <- subset(tbl.models, targetGene.hgnc==targetGene) # & rfScore > 1.5 & pcaMax > 0.75)
   rf.threshold <- 1.5       # fivenum(tbl.model$rfScore)[4]
   pcaMax.threshold <- 0.75  # mean(fivenum(tbl.model$pcaMax)[4], fivenum(tbl.model$pcaMax)[5])
   tbl.model <- subset(tbl.model, rfScore >= rf.threshold & pcaMax >= pcaMax.threshold)
   tbl.model <- tbl.model [, c(2,5:13)]
   tbl.model <- roundNumericColumnsInDataframe(tbl.model, 3, "lassoPValue")
   colnames(tbl.model)[1] <- "TF"
   dim(tbl.model)
   save(tbl.model, file="../shinyApp/data/tbl.model.RData")
   }

if(!exists("igv")){
   igv <- igvR()
   setGenome(igv, genome)
   showGenomicRegion(igv, targetGene)
   tbl.regions <- subset(tbl.enhancers, geneSymbol == "SORL1")
   track <- DataFrameQuantitativeTrack("geneHancer", tbl.regions[, c("chrom", "start", "end", "combinedScore")])
   displayTrack(igv, track)
   big <- sprintf("chr11:%d-%d", min(tbl.regions$start)-5000, max(tbl.regions$end)+5000)
   showGenomicRegion(igv, big)
   }

if(!exists("tbl.gwas")){
   load("~/s/work/priceLab/AD/tbl.gwas.level_1.RData") # tbl.gwas
   tbl.gwas.fixed <- tbl.gwas[, c("CHR", "BP", "BP", "SNP", "P")]
   colnames(tbl.gwas.fixed) <- c("chrom", "start", "end", "rsid", "score")
   tbl.gwas.thisGene <- subset(tbl.gwas.fixed, chrom==chromosome & start >= region.start & end <= region.end)
   tbl.gwas.thisGene$score <- -log10(tbl.gwas.thisGene$score)
   dim(tbl.gwas.thisGene)  # 122 x 5
   save(tbl.gwas.thisGene, file="../shinyApp/data/tbl.gwas.thisGene.RData")
   gr.snp <- GRanges(tbl.gwas.thisGene)
   track <- DataFrameAnnotationTrack("snp", tbl.gwas.thisGene, color="red")
   displayTrack(igv, track)
   tbl.gwas.scored <- tbl.gwas[, c("CHR", "BP", "BP", "P", "SNP")]
   colnames(tbl.gwas.scored) <- c("chrom", "start", "end", "score", "rsid")
   tbl.gwas.scored$score <- -log10(tbl.gwas.scored$score)
   tbl.gwas.scored.thisGene <- subset(tbl.gwas.scored, chrom==chromosome &
                                      start >= min(tbl.regions$start)-5000 &
                                      end <=   max(tbl.regions$end)+5000)
   track <- DataFrameQuantitativeTrack("score", tbl.gwas.scored.thisGene, color="red")
   displayTrack(igv, track)
   }

if(!exists("tbl.bindingSites")){
   tbl.region <- data.frame(chrom=chromosome, start=region.start, end=region.end, stringsAsFactors=FALSE)
   tbl.bindingSites <- findFimoHits(tbl.region, pval=10e-4)
   tbl.bindingSites$start <- tbl.bindingSites$start + region.start
   tbl.bindingSites$stop <- tbl.bindingSites$stop + region.start
   tbl.bindingSites$chrom <- chromosome
   #  with(tbl.bindingSites, plot(-log10(p.value), score))
   with(tbl.bindingSites, plot(-log10(p.value), score))
   tbl.bs <- subset(tbl.bindingSites, score >= 5)
   dim(tbl.bs)   # 2154 x 11
   tbl.bs <- tbl.bs[, c("chrom", "start", "stop", "score", "motif")]
   track <- DataFrameQuantitativeTrack("bs", tbl.bs, color="darkgreen")
   displayTrack(igv, track)
   gr.bs <- GRanges(tbl.bs)
   tbl.ov <- as.data.frame(findOverlaps(gr.snp, gr.bs))
   dim(tbl.ov)
   colnames(tbl.ov) <- c("snp", "bs")
   tbl.snp.bs <- tbl.bs[tbl.ov$bs, c("chrom", "start", "stop", "motif")]
   track <- DataFrameAnnotationTrack("snp/bs", tbl.snp.bs, color="black")
   displayTrack(igv, track)
   save(tbl.snp.bs, file="tbl.snp.bs.RData")
   }

if(!exists("tbl.breaks")){
   motif.names <- unique(tbl.bindingSites$motif)   # 17
   motifs <- MotifDb[motif.names]
   snps.gr <- snps.from.rsid(rsid = tbl.gwas.thisGene$rsid,
                             dbSNP=SNPlocs.Hsapiens.dbSNP150.GRCh38,
                             search.genome=BSgenome.Hsapiens.UCSC.hg38)

   results <- motifbreakR(snpList = snps.gr,
                          filterp = TRUE,
                          pwmList = motifs,
                          show.neutral=TRUE,
                          method = c("ic", "log", "notrans")[2],
                          bkg = c(A=0.25, C=0.25, G=0.25, T=0.25),
                          BPPARAM = BiocParallel::bpparam(),
                          verbose=TRUE)
   results$delta <- results$pctRef - results$pctAlt
   results.filtered <- results[(results$pctRef > 0.70 & results$delta > 0.15) |
                                (results$pctRef > 0.80 & results$delta > 0.10)] # 44
   length(results.filtered)

   range <- seq_len(length(results.filtered))
   range <- 1:3

   for(i in range){
      x <- results.filtered[i]
      rsid <- names(x)
      tf <- mcols(x)$geneSymbol
      filename <- sprintf("../shinyApp/data/png/%s.%s.png", tf, rsid)
      printf("--- writing %s", filename)
      png(filename)
      plotMB(x, rsid)
      dev.off()
      #Sys.sleep(2)
      }
