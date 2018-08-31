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
   printf("finding fimo hits, pval = %f", pval)
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

if(!exists("tbl.fp")){
   genome.db.uri  <- "postgres://khaleesi/hg38"              # has gtf and motifsgenes tables
   fpdbs <- c("brain_hint_20", "brain_hint_16", "brain_wellington_20", "brain_wellington_16")
   tbl.fp.all <- data.frame()
   for(fpdb in fpdbs){
      project.db.uri <- sprintf("postgres://khaleesi/%s", fpdb)
      fpf <- FootprintFinder(genome.db.uri, project.db.uri, quiet=TRUE)
      tbl.fp <- getFootprintsInRegion(fpf, chromosome, region.start, region.end) # tss-200, tss+200)
      dim(tbl.fp)
      tbl.fp.reduced <- unique(tbl.fp[, c("chrom", "start", "endpos")])
      colnames(tbl.fp.reduced) <- c("chrom", "start", "end")
      dim(tbl.fp.reduced)
      tbl.fp.all <- rbind(tbl.fp.all, tbl.fp.reduced)
      }
   dim(tbl.fp.all)
   tbl.fp <- as.data.frame(reduce(GRanges(tbl.fp.all)))
   tbl.fp <- tbl.fp[, c("seqnames", "start", "end")]
   colnames(tbl.fp) <- c("chr", "start", "end")
   tbl.fp$chr <- as.character(tbl.fp$chr)
   dim(tbl.fp)
   track <- DataFrameAnnotationTrack("bh20", tbl.fp.merged, color="darkgreen")
   displayTrack(igv, track)
   save(tbl.fp, file="../shinyApp/data/tbl.fp.RData")
   }

if(!exists("tbl.dhs")){
   hdf <- HumanDHSFilter("hg38", "wgEncodeRegDnaseClustered", pwmMatchPercentageThreshold=0,
                         geneInfoDatabase.uri="bogus", regions=data.frame(), pfms=list())
   tableNames <- getEncodeRegulatoryTableNames(hdf)
   tbl.dhs <- getRegulatoryRegions(hdf, "wgEncodeRegDnaseClustered", chromosome, region.start, region.end)
   colnames(tbl.dhs) <- c("chr", "start", "end", "count", "value")
   tbl.dhs <- tbl.dhs[, c("chr", "start", "end", "value")]
   dim(tbl.dhs)
   track <- DataFrameQuantitativeTrack("dhs", tbl.dhs, color="maroon")
   displayTrack(igv, track)
   save(tbl.dhs, file="../shinyApp/data/tbl.dhs.RData")
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
   tbl.bindingSites <- findFimoHits(tbl.region, pval=10e-2)
   dim(tbl.bindingSites)
   tbl.bindingSites$start <- tbl.bindingSites$start + region.start
   tbl.bindingSites$stop <- tbl.bindingSites$stop + region.start
      # overwrite the native fimo score with -log10(fimo.p.value)
   tbl.bindingSites$score <- -log10(tbl.bindingSites$p.value)
   tbl.bindingSites$chrom <- chromosome
   #  with(tbl.bindingSites, plot(-log10(p.value), score))

   tbl.bs <- subset(tbl.bindingSites, score >= 2)
   dim(tbl.bs)   # 151k for score >= 2
   tbl.bs <- tbl.bs[, c("chrom", "start", "stop", "score", "motif")]
   colnames(tbl.bs) <- c("chr", "start", "end", "score", "motif")
   id.map <- mcols(MotifDb[unique(tbl.bs$motif)])$geneSymbol
   names(id.map) <- unique(tbl.bs$motif)
   tfs <- as.character(id.map[tbl.bs$motif])
   tbl.bs$tf <- toupper(tfs)    #  humanize any mouse tf geneSymbols (as in Id2 for SORL1)
      # this binding site must be found.  it overlaps with motifbreakR's rs17125473 for CEBPA
      # with a score of 2.3
   checkEquals(nrow(subset(tbl.bs, tf == "CEBPA" & start > 121574082 & end < 121574123)), 1)

   dim(tbl.bs)
   save(tbl.bs, file="../shinyApp/data/tbl.bs.RData")
   track <- DataFrameQuantitativeTrack("bs", tbl.bs, color="darkgreen")
   displayTrack(igv, track)
   gr.bs <- GRanges(tbl.bs)
   tbl.ov <- as.data.frame(findOverlaps(gr.snp, gr.bs))
   dim(tbl.ov)
   colnames(tbl.ov) <- c("snp", "bs")
   tbl.bs.snp <- cbind(tbl.bs[tbl.ov$bs,], tbl.gwas.thisGene[tbl.ov$snp,])
   dim(tbl.bs.snp)
   colnames(tbl.bs.snp) <- c("chrom", "start", "end", "motifScore", "motif", "tf", "chrom.snp", "loc", "junk", "rsid", "snpScore")
   dim(tbl.bs.snp)
   tbl.snp.bs <- tbl.bs.snp[,  c("chrom", "start", "end", "motifScore", "tf", "motif", "snpScore", "rsid", "loc")]
   track <- DataFrameAnnotationTrack("snp/bs", tbl.snp.bs, color="black")
   displayTrack(igv, track)
   save(tbl.snp.bs, file="../shinyApp/data/tbl.snp.bs.RData")
   }

if(!exists("tbl.breaks")){
   motif.names <- unique(tbl.bs$motif)   # 17
   length(motif.names)
   motifs <- MotifDb[motif.names]
   snps.gr <- snps.from.rsid(rsid = tbl.gwas.thisGene$rsid,
                             dbSNP=SNPlocs.Hsapiens.dbSNP150.GRCh38,
                             search.genome=BSgenome.Hsapiens.UCSC.hg38)

   results <- motifbreakR(snpList = snps.gr,
                          filterp = TRUE,
                          pwmList = motifs,
                          show.neutral=FALSE,
                          method = c("ic", "log", "notrans")[2],
                          bkg = c(A=0.25, C=0.25, G=0.25, T=0.25),
                          BPPARAM = BiocParallel::bpparam(),
                          verbose=TRUE)
   results$delta <- results$pctRef - results$pctAlt
   results.filtered <- results[(results$pctRef > 0.70 & results$delta > 0.15) |
                                (results$pctRef > 0.80 & results$delta > 0.10)] # 44
   length(results.filtered)
   tbl.broken <- as.data.frame(mcols(results.filtered))
   tbl.broken$geneSymbol <-  toupper(tbl.broken$geneSymbol)
   all(tbl.broken$geneSymbol %in% tbl.model$TF)
   setdiff(tbl.broken$geneSymbol, tbl.model$TF)
   tbl.broken$chr <- chromosome
   tbl.broken$rsid <- rownames(tbl.broken)
   rownames(tbl.broken) <- NULL

   save(tbl.broken, file="../shinyApp/data/tbl.broken.RData")

   range <- seq_len(length(results.filtered))

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

   } # tbl.breaks
