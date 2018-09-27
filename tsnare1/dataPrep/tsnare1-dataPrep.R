#   The variant is rs4976977, described in the 2014 scz gwas paper as a “credible” (imputed?) snp
#   related to the index snp rs4129585, and with a crediblePVal of 3.043e-11.
# chr8:142,212,080-142,403,291(GRCh38/hg38)
# Size:191,212 basesOrientation:Minus strand
# chr8:143,293,441-143,484,601(GRCh37/hg19)
#------------------------------------------------------------------------------------------------------------------------
library(igvR)
library(trenaSGM)
library(RSQLite)
library(MotifDb)
source("~/github/trenaShinyApps/utils/plotMotifs.R")
source("~/s/data/public/GeneHancer/v4.7/getAllEnhancers.R")
#------------------------------------------------------------------------------------------------------------------------
genome <- "hg38"
targetGene <- "TSNARE1"
chromosome <- "chr8"
tss <- 142403291
snp.loc <- 142232793
snp <- "rs4976977"

if(!exists("igv")){
   igv <- igvR()
   setGenome(igv, "hg38")
   showGenomicRegion(igv, "TSNARE1")
   }

if(!exists("tbl.enhancers")){
   tbl.enhancers <- getAllEnhancers(targetGene)
   colnames(tbl.enhancers) <- c("chr", "start", "end", "type", "score")
   tbl.wig <- tbl.enhancers[, c("chr", "start", "end", "score")]
   gr.enhancers <- GRanges(tbl.enhancers)
   track <- DataFrameQuantitativeTrack("enhancers", tbl.wig, color="black", trackHeight=25, autoscale=TRUE)
   displayTrack(igv, track)
   track <- DataFrameAnnotationTrack(snp, data.frame(chr=chromosome, start=snp.loc, end=snp.loc, stringsAsFactors=FALSE),
                                     color="red", trackHeight=25)
   displayTrack(igv, track)
   }

shoulder <- 50000
region.min <- min(tbl.enhancers$start) - shoulder
region.max <- max(tbl.enhancers$end) + shoulder

if(!exists("tbl.dhs")){
   hdf <- HumanDHSFilter("hg38", "wgEncodeRegDnaseClustered", pwmMatchPercentageThreshold=0,
                         geneInfoDatabase.uri="bogus", regions=data.frame(), pfms=list())
   tableNames <- getEncodeRegulatoryTableNames(hdf)
   tbl.dhs <- getRegulatoryRegions(hdf, "wgEncodeRegDnaseClustered", chromosome, region.min, region.max)
   colnames(tbl.dhs) <- c("chr", "start", "end", "count", "value")
   tbl.wig <- tbl.dhs[, c("chr", "start", "end", "value")]
   track <- DataFrameQuantitativeTrack("DHS", tbl.wig, color="gray", trackHeight=25, autoscale=TRUE)
   displayTrack(igv, track)

   }

if(!exists("tbl.gwas")){
   load("~/s/work/priceLab/AD/tbl.gwas.level_1.RData") # tbl.gwas
   tbl.snps <- tbl.gwas[, c("CHR", "BP", "BP", "SNP", "Effect_allele", "Non_Effect_allele", "P")]
   colnames(tbl.snps) <- c("chrom", "start", "end", "rsid", "variant", "reference", "pval")
   tbl.snps$pScore <- -log10(tbl.snps$pval)
   gr.snps <- GRanges(tbl.snps)
   tbl.wig <- tbl.snps[, c("chrom", "start", "end", "pScore")]
   tbl.wig <- subset(tbl.wig, chrom==chromosome & start >= region.min & end <= region.max)
   track <- DataFrameQuantitativeTrack("gwas", tbl.wig, color="red", trackHeight=25, autoscale=TRUE)
   displayTrack(igv, track)
   }

if(!exists("tbl.nfix.chipSeq")){
   db <- dbConnect(dbDriver("SQLite"), "~/s/data/public/human/remap-2018/remap-all.sqlite")
   dbGetQuery(db, "select * from chipseq limit 5")
   dbGetQuery(db, "select * from chipseq where tf='NFIX' limit 5")   # no results
   }

if(!exists("tbl.model")){
   load("/Users/paul/github/trena/inst/extdata/tbl.models.all.RData")
   }

if(!exists("mtx.cer"))
   print(load("mayo.rnaSeq.cer.and.tcx.matrices.RData"))   # mtx.cer, mtx.tcx

if(!exists("sgm"))
   sgm <- trenaSGM(genome, targetGene, quiet=FALSE)



#------------------------------------------------------------------------------------------------------------------------
viz.setup <- function()
{
   showGenomicRegion(igv, "TSNARE1")
   enhancers.track <- DataFrameAnnotationTrack("enhancers", tbl.enhancers, color="black", displayMode="SQUISHED")
   displayTrack(igv, enhancers.track)
   snp.track <- DataFrameAnnotationTrack("snp",
                                         data.frame(chrom=chromosome, start=snp.loc, end=snp.loc, stringsAsFactors=FALSE),
                                         color="red")
   displayTrack(igv, snp.track)

} # viz.setup
#------------------------------------------------------------------------------------------------------------------------
fp.promoter.enhancers <- function()
{
   tbl.promoter <- data.frame(chrom=chromosome, start=tss-200, end=tss+2000, stringsAsFactors=FALSE)
   tbl.snpRegion <- data.frame(chrom=chromosome, start=142232639, end=142233021, stringsAsFactors=FALSE)
   #tbl.regions <- rbind(rbind(tbl.enhancers, tbl.promoter), tbl.snpRegion)
   #tbl.regions <- tbl.regions[order(tbl.regions$start, decreasing=FALSE),]
   tbl.regions <- tbl.promoter

   spec <- list(title="fp.2kup.200down.02",
                type="footprint.database",
                regions=tbl.regions,
                tss=tss,
                matrix=mtx.cer,
                db.host="khaleesi.systemsbiology.net",
                databases=list("brain_hint_20", "brain_hint_16", "brain_wellington_20", "brain_wellington_16"),
                motifDiscovery="builtinFimo",
                tfMapping=c("MotifDB"), #, "TFClass"),
                tfPrefilterCorrelation=0.2,
                orderModelByColumn="rfScore",
                solverNames=c("lasso", "lassopv", "pearson", "randomForest", "ridge", "spearman"))

   return(spec)


} # fp.promoter.enhancers
#------------------------------------------------------------------------------------------------------------------------
mm.promoter <- function()
{
   tbl.promoter <- data.frame(chrom=chromosome, start=tss-200, end=tss+2000, stringsAsFactors=FALSE)
   tbl.snpRegion <- data.frame(chrom=chromosome, start=142232639, end=142233021, stringsAsFactors=FALSE)
   #tbl.regions <- rbind(rbind(tbl.enhancers, tbl.promoter), tbl.snpRegion)
   #tbl.regions <- tbl.regions[order(tbl.regions$start, decreasing=FALSE),]
   tbl.regions <- tbl.promoter

   spec <- list(title="mm.promoter",
                type="regions.motifMatching",
                tss=tss,
                regions=tbl.regions,
                matrix=mtx.cer,
                pfms=query(MotifDb, "sapiens", c("jaspar2018", "hocomoco")),
                matchThreshold=80,
                motifDiscovery="matchPWM",
                tfMapping=c("MotifDb"), # , "TFClass"),
                tfPrefilterCorrelation=0.2,
                orderModelByColumn="rfScore",
                solverNames=c("lasso", "lassopv", "pearson", "randomForest", "ridge", "spearman"))

   return(spec)


} # mm.promoter
#------------------------------------------------------------------------------------------------------------------------
nfix.bindingSites <- function()
{
   tbl.regions <- data.frame(chrom="chr8", start=142147780, end=142567603, stringsAsFactors=FALSE)
   pfms <- query(MotifDb, c("sapiens", "jaspar2018", "NFIX"))
   mm <- MotifMatcher("hg38", as.list(pfms), quiet=TRUE)
   tbl.hits <- findMatchesByChromosomalRegion(mm, tbl.regions, 90, variants = NA_character_)
   track <- DataFrameQuantitativeTrack("MA0671.1", tbl.hits[, c("chrom", "motifStart", "motifEnd", "motifRelativeScore")], color="green")
   displayTrack(igv, track)
   # assessSnp(trena, as.list(pfms), "rs4976977", 10, 80)
   #                           motifName status assessed motifRelativeScore delta                                     signature chrom motifStart  motifEnd strand     match   variant
   # 1 Hsapiens-jaspar2018-NFIX-MA0671.1     wt  in.both          0.9675002     0 Hsapiens-jaspar2018-NFIX-MA0671.1;142232792;-  chr8  142232792 142232800      - GTTGGCAGA rs4976977
   # 2 Hsapiens-jaspar2018-NFIX-MA0671.1    mut  in.both          0.8618954     0 Hsapiens-jaspar2018-NFIX-MA0671.1;142232792;-  chr8  142232792 142232800      - GATGGCAGA rs4976977


} # nfix.bindingSites
#------------------------------------------------------------------------------------------------------------------------
run <- function()
{
   models.01 <- calculate(sgm, list(fp.promoter.enhancers=fp.promoter.enhancers()))
   models.02 <- calculate(sgm, list(mm.promoter=mm.promoter()))   # NFIX is tf rank 3, 37 binding sites
   tbl.nfixBindingSites <- subset(tbl.reg, geneSymbol=="NFIX")[, c("chrom", "motifStart", "motifEnd", "motifRelativeScore")]
   track.nfix <- DataFrameAnnotationTrack("NFIX", tbl.nfixBindingSites, color="green")
   displayTrack(igv, track.nfix)

} # run
#------------------------------------------------------------------------------------------------------------------------
