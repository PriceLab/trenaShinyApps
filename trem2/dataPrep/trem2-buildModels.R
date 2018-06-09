library(trenaSGM)
library(igvR)
library(FimoClient)
library(GenomicRanges)
library(RUnit)
#------------------------------------------------------------------------------------------------------------------------
source("~/github/trenaShinyApps/utils/getEnhancers.R")
#------------------------------------------------------------------------------------------------------------------------
genome <- "hg38"
targetGene <- "TREM2"
sgm <- trenaSGM(genome, targetGene, quiet=FALSE)
chromosome <- "chr6"
tss <- 41163186
rs9357347.loc <- "chr6:41182853-41182853"
rs9357347.bp <- 41182853
#------------------------------------------------------------------------------------------------------------------------
if(!exists("igv")){
   igv <- igvR()
   Sys.sleep(5)
   setGenome(igv, genome)
   Sys.sleep(5)
   showGenomicRegion(igv, targetGene)
   }

if(!all(exists(c("mtx.cer", "mtx.ros", "mtx.tcx"))))
   load("mtx.withDimers.cer.ros.tcx.RData")

#------------------------------------------------------------------------------------------------------------------------
FIMO_HOST <- "localhost"
FIMO_PORT <- 5558
if(!exists("fimo")){
   fimo <- FimoClient(FIMO_HOST, FIMO_PORT, quiet=FALSE)
   }
#------------------------------------------------------------------------------------------------------------------------
# minerva 2016:
# We restricted our analysis to variants located within 100kb of any coding TREM family gene at the
# chromosome 6p21.1 TREM gene cluster (Fig. 1). Variants were further selected based on the
# statistical significance of their AD-risk association in the IGAP stage 1 meta- analysis [12]
# (Supplementary Methods), where only those variants with p-values  0.0015 were kept. This p-value
# cut-off was arbitrarily chosen to select those variants that existed in both the IGAP stage 1 AD
# GWAS and our discovery eQTL cohort, Mayo Clinic Whole Genome-DASL dataset, and that could be
# genotyped, if needed, in the replication eQTL cohorts, using cost- Regulome scores were obtained
# from the Regulome database, which annotates variants with regulatory information from 962 different
# datasets and a variety of sources, including ENCODE [16]. Regulome scores are on a scale from 1 to
# 6, and these numerical categories are sub- classified with letters based on the number of lines of
# evidence of functional consequence. A value of 1a is assigned to the variant with the most evidence
# of regulatory potential, while a score of 6 has the least [16].

# grep("^TREM", as.data.frame(org.Hs.egSYMBOL)$symbol, v=TRUE)
#  "TREM2"   "TREM1"   "TREML2"  "TREML5P" "TREML4"  "TREML1"  "TREML3P"

if(!exists("tbl.tremGeneFamilyGenes")){
   library(RPostgreSQL)
   database.host <- "bddsrds.globusgenomics.org"
   db.hg38 <- dbConnect(PostgreSQL(), user= "trena", password="trena", dbname="hg38", host=database.host)
   db.query <- "select * from gtf where moleculetype='gene' and gene_biotype='protein_coding' and gene_name like 'TREM%'"
   tbl.tremGeneFamilyGenes <- dbGetQuery(db.hg38, db.query)[, c("chr", "start", "endpos", "gene_name")]
   tbl.regionMax <- with(tbl.tremGeneFamilyGenes, data.frame(chrom=unique(chr), start=min(start)-100000, end=max(start)+100000))
   # 318,585 bases
   regionsMax <- with(tbl.regionMax, sprintf("%s:%d-%d", chrom, start, end))
   showGenomicRegion(igv, regionsMax)
   }

if(!exists("tbl.snp")){
    load("~/s/work/priceLab/AD/tbl.gwas.level_1.RData") # tbl.gwas
    tbl.snp <- subset(tbl.gwas, BP >= tbl.regionMax$start & BP <= tbl.regionMax$end & CHR=="chr6")
    tbl.snp <- tbl.snp[, c("CHR", "BP", "BP", "SNP", "Effect_allele", "Non_Effect_allele", "Beta", "SE", "P")]
    colnames(tbl.snp) <- c("chrom", "start", "end", "rsid", "wt", "mut", "beta", "se", "snpPval")
    tbl.snp <- subset(tbl.snp, snpPval <= 0.0015)   # 28 rows
    track <- DataFrameAnnotationTrack("gwas", tbl.snp[, c("chrom", "start", "end", "rsid")], color="darkRed")
    displayTrack(igv, track)
    }

if(!exists("tbl.eqtl")){
   tbl.eqtl <- read.table("eqtl-hg38.tsv", sep="\t", as.is=TRUE, header=TRUE)
   track <- DataFrameAnnotationTrack("eqtl", tbl.eqtl[, 1:4], "blue")
   displayTrack(igv, track)
   }

if(!exists("sgm"))
   sgm <- trenaSGM(genome, targetGene)

if(!exists("tbl.enhancers.trem2")){
   tbl.enhancers.trem2 <- getEnhancerRegions(targetGene)
   track <- DataFrameAnnotationTrack("trem2 enhancers", tbl.enhancers.trem2, color="darkgrey")
   displayTrack(igv, track)
   tbl.enhancers.treml1 <- getEnhancerRegions("TREML1")
   track <- DataFrameAnnotationTrack("treml1 enhancers", tbl.enhancers.treml1, color="darkgrey")
   displayTrack(igv, track)
   track <- DataFrameAnnotationTrack("rs9357347", data.frame(chrom="chr6", start=rs9357347.bp, end=rs9357347.bp),
                                     color="red")
   displayTrack(igv, track)
   }

if(!exists("tbl.dhs")){
   hdf <- HumanDHSFilter("hg38", "wgEncodeRegDnaseClustered", pwmMatchPercentageThreshold=0,
                         geneInfoDatabase.uri="bogus", regions=data.frame(), pfms=list())
   tbl.dhs <- with(tbl.regionMax, getRegulatoryRegions(hdf, "wgEncodeRegDnaseClustered", chrom, start, end))
   save(tbl.dhs, file="../shinyApp/dhs.RData")
   track <- DataFrameAnnotationTrack("DHS", tbl.dhs, color="darkgray")
   displayTrack(igv, track)
   }

#------------------------------------------------------------------------------------------------------------------------
find.lyl1.spec <- function(mtx)
{
   tbl.regions <- data.frame(chrom=chromosome, start=tss-200, end=tss+2000, stringsAsFactors=FALSE)
   track <-  DataFrameAnnotationTrack("e2", tbl.regions, color="magenta")
   displayTrack(igv, track)

   print(load(system.file(package="trenaSGM", "extdata", "enhancers.TREM2.RData")))

   track <- DataFrameAnnotationTrack("e1", tbl.enhancers, color="orange")
   displayTrack(igv, track)

   spec <- list(title="find.lyl1",
                type="footprint.database",
                regions=tbl.enhancers,
                tss=tss,
                matrix=mtx,
                db.host="khaleesi.systemsbiology.net",
                databases=list("brain_hint_20", "brain_hint_16", "brain_wellington_20", "brain_wellington_16"),
                motifDiscovery="builtinFimo",
                tfMapping=c("MotifDb"),
                tfPrefilterCorrelation=0.2,
                orderModelByColumn="rfScore",
                solverNames=c("lasso", "lassopv", "pearson", "randomForest", "ridge", "spearman"))

   return(spec)

} # find.lyl1.spec
#------------------------------------------------------------------------------------------------------------------------
fp.2000up.200down.cor02 <- function(mtx)
{
   tbl.regions <- data.frame(chrom=chromosome, start=tss-200, end=tss+2000, stringsAsFactors=FALSE)

   spec <- list(title="fp.2000up.200down.02",
                type="footprint.database",
                regions=tbl.regions,
                tss=tss,
                matrix=mtx,
                db.host="khaleesi.systemsbiology.net",
                databases=list("brain_hint_20", "brain_hint_16", "brain_wellington_20", "brain_wellington_16"),
                motifDiscovery="builtinFimo",
                tfMapping=c("MotifDb"),
                tfPrefilterCorrelation=0.2,
                orderModelByColumn="rfScore",
                solverNames=c("lasso", "lassopv", "pearson", "randomForest", "ridge", "spearman"))

   return(spec)

} # fp.2000up.200down.cor02
#------------------------------------------------------------------------------------------------------------------------
fp.5kup.5kdown.cor02 <- function(mtx)
{
   tbl.regions <- data.frame(chrom=chromosome, start=tss-5000, end=tss+5000, stringsAsFactors=FALSE)

   spec <- list(title="fp.5kup.5kdown.02",
                type="footprint.database",
                regions=tbl.regions,
                tss=tss,
                matrix=mtx,
                db.host="khaleesi.systemsbiology.net",
                databases=list("brain_hint_20", "brain_hint_16", "brain_wellington_20", "brain_wellington_16"),
                motifDiscovery="builtinFimo",
                tfMapping=c("MotifDb"),
                tfPrefilterCorrelation=0.2,
                orderModelByColumn="rfScore",
                solverNames=c("lasso", "lassopv", "pearson", "randomForest", "ridge", "spearman"))

  #  build.spec <- list(title="fp.enhancers",
#                       type="footprint.database",
#                       regions=data.frame(chrom="chr6", start=tss-5000, end=tss+5000),
#                       tss=tss,
#                       matrix=mtx,
#                       db.host="khaleesi.systemsbiology.net",
#                       databases=c("brain_hint_20",  "brain_hint_16", "brain_wellington_20", "brain_wellington_16"),
#                       motifDiscovery="builtinFimo",
#                       tfMapping=c("MotifDb"),
#                       tfPrefilterCorrelation=0.4,
#                       orderModelByColumn="rfScore",
#                       solverNames=c("lasso", "lassopv", "pearson", "randomForest", "ridge", "spearman"))



   return(spec)

} # fp.5kup.5kdown.cor02
#------------------------------------------------------------------------------------------------------------------------
fp.enhancers.cor02 <- function(tbl.enhancers, mtx)
{
   spec <- list(title="fp.enhancers.02",
                type="footprint.database",
                regions=tbl.enhancers,
                tss=tss,
                matrix=mtx,
                db.host="khaleesi.systemsbiology.net",
                databases=list("brain_hint_20", "brain_hint_16", "brain_wellington_20", "brain_wellington_16"),
                motifDiscovery="builtinFimo",
                tfMapping=c("MotifDb"),
                tfPrefilterCorrelation=0.2,
                orderModelByColumn="rfScore",
                solverNames=c("lasso", "lassopv", "pearson", "randomForest", "ridge", "spearman"))

   return(spec)

} # fp.enhancers.cor02
#------------------------------------------------------------------------------------------------------------------------
fp.enhancers.cor02.5kup.5kdown <- function()
{

   tbl.regions <- rbind(getEnhancerRegions(targetGene),
                        data.frame(chrom=chromosome, start=tss-5000, end=tss+5000, stringsAsFactors=FALSE))
   tbl.regions <- tbl.regions[order(tbl.regions$start, decreasing=FALSE),]

   spec <- list(title="fp.enhancers.02",
                type="footprint.database",
                regions=tbl.regions,
                tss=tss,
                matrix=mtx,
                db.host="khaleesi.systemsbiology.net",
                databases=list("brain_hint_20", "brain_hint_16", "brain_wellington_20", "brain_wellington_16"),
                motifDiscovery="builtinFimo",
                tfMapping=c("MotifDb"),
                tfPrefilterCorrelation=0.2,
                orderModelByColumn="rfScore",
                solverNames=c("lasso", "lassopv", "pearson", "randomForest", "ridge", "spearman"))

   return(spec)


}  # fp.enhancers.cor02.5kup.5kdown
#------------------------------------------------------------------------------------------------------------------------
regions.motifMatching.5kup.5kdown.pwm80.cor02 <- function()
{
   tbl.regions <- data.frame(chrom=chromosome, start=tss-5000, end=tss+5000, stringsAsFactors=FALSE)

   spec <- list(title="rmm.5kup.5kdown",
                type="regions.motifMatching",
                tss=tss,
                regions=tbl.regions,
                matrix=mtx,
                pfms=query(MotifDb, "sapiens", "jaspar2018"),
                matchThreshold=80,
                motifDiscovery="matchPWM",
                tfMapping=c("MotifDb"),
                tfPrefilterCorrelation=0.2,
                orderModelByColumn="rfScore",
                solverNames=c("lasso", "lassopv", "pearson", "randomForest", "ridge", "spearman"))


   return(spec)

} # regions.motifMatching.5kup.5kdown.pwm80.cor02
#------------------------------------------------------------------------------------------------------------------------
regions.motifMatching.2kup.200down.pwm80.cor02 <- function()
{
   tbl.regions <- data.frame(chrom=chromosome, start=tss-200, end=tss+2000, stringsAsFactors=FALSE)

   spec <- list(title="rmm.5kup.5kdown",
                type="regions.motifMatching",
                tss=tss,
                regions=tbl.regions,
                matrix=mtx,
                pfms=query(MotifDb, "sapiens", "jaspar2018"),
                matchThreshold=80,
                motifDiscovery="matchPWM",
                tfMapping=c("MotifDb"),
                tfPrefilterCorrelation=0.2,
                orderModelByColumn="rfScore",
                solverNames=c("lasso", "lassopv", "pearson", "randomForest", "ridge", "spearman"))


   return(spec)

} # regions.motifMatching.2kup.200down.pwm80.cor02
#------------------------------------------------------------------------------------------------------------------------
# * cory, minerva describe rs9357347 possibly affecting TREM2 regulation  (5 apr 2016)
#     hg38: chr6:41182853-41182853
#    --> falls right between 2 low-scoring lyl1-related motif matches, for Hsapiens-jaspar2018-NHLH1-MA0048.2
#
# rs2093395 intersects lyl1-related 233  chr6 41187287 41187304 Hsapiens-jaspar2018-GATA1::TAL1-MA0140.2
# --- fimo finds this
#                                     motif     start     stop strand   score  p.value q.value   matched.sequence
#  Hsapiens-jaspar2018-GATA1::TAL1-MA0140.2  41187287 41187304      - 7.01124 0.000294    0.95 CTTATCTATCTACGACTT
# ---- pwmMatch finds this
#                                motifName chrom motifStart motifEnd strand motifScore motifRelativeScore              match status shortMotif
#  Hsapiens-jaspar2018-GATA1::TAL1-MA0140.2  chr6   41187286 41187303      -   8.546317          0.7681438 AAGTCGTAGATAGATAAG     wt   MA0140.2
#
find.lyl1 <- function()
{
   track <- DataFrameAnnotationTrack("rs9357347", data.frame(chrom="chr6", start=41182853, end=41182853), color="darkGreen")
   displayTrack(igv, track)

   tbl.bigRegion <- with(tbl.enhancers, data.frame(chrom=unique(chrom), start=min(start), end=max(end), stringsAsFactors=FALSE))
   bigRegion <- with(tbl.bigRegion, sprintf("%s:%d-%d", unique(chrom), min(start)-2000, max(end)+2000))
   showGenomicRegion(igv, bigRegion)
   track <- DataFrameAnnotationTrack("enhancers", tbl.enhancers.trem2, color="black")
   displayTrack(igv, track)

      #--------------------------------------------------------------------------------
      # first try MotifMatcher on rs2093395
      #--------------------------------------------------------------------------------

   pfms.lyl1 <- query(MotifDb, c("jaspar2018", "sapiens"), geneToMotif(MotifDb, "LYL1", c("MotifDb", "TFClass"))$motif)
   mm <- MotifMatcher("hg38", as.list(pfms.lyl1), quiet=TRUE)
   tbl.matches.lyl1 <- findMatchesByChromosomalRegion(mm, tbl.bigRegion, pwmMatchMinimumAsPercentage=70)
   tbl.tmp <- tbl.matches.lyl1[, c("chrom", "motifStart", "motifEnd", "motifName")]
   track <- DataFrameAnnotationTrack("mm lyl1", tbl.tmp, color="blue")
   displayTrack(igv, track)

   track <- DataFrameAnnotationTrack("eqtl", tbl.eqtl, color="red")
   displayTrack(igv, track)

   track <- DataFrameAnnotationTrack("gwas", tbl.snp, color="red")
   displayTrack(igv, track)

      #--------------------------------------------------------------------------------
      # now try fimo on rs2093395
      #--------------------------------------------------------------------------------

   export(pfms.lyl1, "~/github/fimoService/pfms/lyl1-related.meme", 'meme')
     # restart fimo: cd ~/github/fimoService/server;
   tbl.sequence <- getSequence(mm, tbl.bigRegion)
   tbl.fimo.lyl1 <- requestMatch(fimo, list(big=tbl.sequence$seq))
   tbl.fimo.lyl1 <- subset(tbl.fimo.lyl1, p.value < 0.001)  # 7356 lines
   tbl.fimo.lyl1$start <- tbl.fimo.lyl1$start + tbl.bigRegion$start
   tbl.fimo.lyl1$stop <- tbl.fimo.lyl1$stop + tbl.bigRegion$start
      tbl.tmp <- tbl.fimo.lyl1[, c("start", "stop", "motif")]
   tbl.tmp$chrom <- "chr6"
   tbl.tmp <- tbl.tmp[, c("chrom", "start", "stop", "motif")]
   track.lyl1 <- DataFrameAnnotationTrack("fimo", tbl.tmp, color="darkred")
   displayTrack(igv, track.lyl1)


     #--------------------------------------------------------------------------------
     # and the footprint databases, via trenaSGM?
     #--------------------------------------------------------------------------------

   tbl.regions <- tbl.enhancers.trem2 #data.frame(chrom="chr6", start=41186160, end=41187527)

   spec <- list(title="fp.trem2.ehnancers",
                type="footprint.database",
                regions=tbl.regions,
                tss=tss,
                matrix=mtx,
                db.host="khaleesi.systemsbiology.net",
                databases=list("brain_hint_20", "brain_hint_16", "brain_wellington_20", "brain_wellington_16"),
                motifDiscovery="builtinFimo",
                tfMapping=c("MotifDb"),
                tfPrefilterCorrelation=0.2,
                orderModelByColumn="rfScore",
                solverNames=c("lasso", "lassopv", "pearson", "randomForest", "ridge", "spearman"))

   x <- calculate(sgm, list(fp.368=spec))
   checkTrue("LYL1" %in% head(x[[1]]$model$gene))

} # find.lyl1
#------------------------------------------------------------------------------------------------------------------------
regions.motifMatching.2kup.200down.pwm90.cor02 <- function()
{
   tbl.regions <- data.frame(chrom=chromosome, start=tss-200, end=tss+2000, stringsAsFactors=FALSE)

   spec <- list(title="rmm.5kup.5kdown",
                type="regions.motifMatching",
                tss=tss,
                regions=tbl.regions,
                matrix=mtx,
                pfms=query(MotifDb, "sapiens", "jaspar2018"),
                matchThreshold=90,
                motifDiscovery="matchPWM",
                tfMapping=c("MotifDb"),
                tfPrefilterCorrelation=0.2,
                orderModelByColumn="rfScore",
                solverNames=c("lasso", "lassopv", "pearson", "randomForest", "ridge", "spearman"))


   return(spec)

} # regions.motifMatching.2kup.200down.pwm90.cor02
#------------------------------------------------------------------------------------------------------------------------
regions.motifMatching.2kup.200down.pwm95.cor02 <- function()
{
   tbl.regions <- data.frame(chrom=chromosome, start=tss-200, end=tss+2000, stringsAsFactors=FALSE)

   spec <- list(title="rmm.5kup.5kdown",
                type="regions.motifMatching",
                tss=tss,
                regions=tbl.regions,
                matrix=mtx,
                pfms=query(MotifDb, "sapiens", "jaspar2018"),
                matchThreshold=95,
                motifDiscovery="matchPWM",
                tfMapping=c("MotifDb"),
                tfPrefilterCorrelation=0.2,
                orderModelByColumn="rfScore",
                solverNames=c("lasso", "lassopv", "pearson", "randomForest", "ridge", "spearman"))


   return(spec)

} # regions.motifMatching.2kup.200down.pwm95.cor02
#------------------------------------------------------------------------------------------------------------------------
regions.motifMatching.enhancers.5kbup.5kbdown.pwm80.cor04 <- function()
{
   tbl.regions <- rbind(getEnhancerRegions(targetGene),
                        data.frame(chrom=chromosome, start=tss-5000, end=tss+5000, stringsAsFactors=FALSE))

   spec <- list(title="rmm.enhancers,5kup.5kdown",
                type="regions.motifMatching",
                tss=tss,
                regions=tbl.regions,
                matrix=mtx,
                pfms=query(MotifDb, "sapiens", "jaspar2018"),
                matchThreshold=80,
                motifDiscovery="matchPWM",
                tfMapping=c("MotifDb"),
                tfPrefilterCorrelation=0.4,
                orderModelByColumn="rfScore",
                solverNames=c("lasso", "lassopv", "pearson", "randomForest", "ridge", "spearman"))


   return(spec)

} # regions.motifMatching.enhancers.5kbup.5kbdown.pwm80.cor04
#------------------------------------------------------------------------------------------------------------------------
regions.motifMatching.snpEnhancerOnly <- function()
{
   tbl.regions <- data.frame(chrom=chromosome, start=41182688, end=41183027, stringsAsFactors=FALSE)

   spec <- list(title="rmm.snpEnhancerOnly",
                type="regions.motifMatching",
                tss=tss,
                regions=tbl.regions,
                matrix=mtx,
                pfms=query(MotifDb, "sapiens", "jaspar2018"),
                matchThreshold=85,
                motifDiscovery="matchPWM",
                tfMapping=c("MotifDb"),
                tfPrefilterCorrelation=0.4,
                orderModelByColumn="rfScore",
                solverNames=c("lasso", "lassopv", "pearson", "randomForest", "ridge", "spearman"))

   return(spec)

} # regions.motifMatching.snpEnhancerOnly
#------------------------------------------------------------------------------------------------------------------------
noDNA.allTFs.cor04 <- function(mtx)
{
   spec <- list(title="noDNA.allTFs.cor04",
                type="noDNA.tfsSupplied",
                matrix=mtx,
                tfs=allKnownTFs(),
                tfPrefilterCorrelation=0.4,
                orderModelByColumn="rfScore",
                solverNames=c("lasso", "lassopv", "pearson", "randomForest", "ridge", "spearman"))

   return(spec)

} # noDNA.allTFs.cor04
#------------------------------------------------------------------------------------------------------------------------
identifyMotifs <- function()
{

  load("../shinyApp/trem2-models.RData")
  tbl.summary <- all.models$summary$model
  tfs <- subset(tbl.summary, observed > 1 & rank.sum < 50)$TF
  length(tfs)

  pfms.human.jhs <- query(MotifDb, "sapiens", c("jaspar2018", "hocomoco", "swissregulon"))
  #pfms.human.jhs <- query(MotifDb, "sapiens", c("jaspar2018", "hocomoco"))
  length(pfms.human.jhs)
  pfms.tfs <- query(pfms.human.jhs, "", tfs)
  length(pfms.tfs)
  dimers <-  grep("::", names(pfms.tfs))
  if(length(dimers) > 0)
     pfms.tfs <- pfms.tfs[-dimers]
  length(pfms.tfs)

  tfs.with.motifs <- mcols(pfms.tfs)$geneSymbol
  setdiff(tfs, tfs.with.motifs)   # "PLEK" "LYL1" "VAV1"
  tbl.xtab <- as.data.frame(table(tfs.with.motifs))
  tbl.xtab <- tbl.xtab[order(tbl.xtab$Freq, decreasing=TRUE),]

  indices <- which(tfs.with.motifs == "STAT1")
  indices  # 34 35 36 37 38 39
  plotMotifs(pfms.tfs[indices])  # chose the first one, #34

  indices <- which(tfs.with.motifs == "SPI1")
  indices  # 3 4 5 6 7
  plotMotifs(pfms.tfs[indices])  # chose the long third one, # 5

  indices <- which(tfs.with.motifs == "BCL6B")
  indices  # 16 18 20
  plotMotifs(pfms.tfs[indices])  # chose the long first one, # 16

  indices <- which(tfs.with.motifs == "CEBPA")
  indices   # 11 22 23
  plotMotifs(pfms.tfs[indices])  # chose the long first one, # 21

  indices <- which(tfs.with.motifs == "EGR4")
  indices   #  28 29 30
  plotMotifs(pfms.tfs[indices])  # chose 28

  indices <- which(tfs.with.motifs == "ELK3")
  indices   # 40 41 42
  plotMotifs(pfms.tfs[indices])  # chose the first one, 41

  indices <- which(tfs.with.motifs == "FOXP1")
  indices   #  31 32 33
  plotMotifs(pfms.tfs[indices])  # chose 33

  indices <- which(tfs.with.motifs == "IRF5")
  indices   #  8 9 10
  plotMotifs(pfms.tfs[indices])  # chose 8, but 9 hocomoco more strict, interesting

  indices <- which(tfs.with.motifs == "IRF8")
  indices   # 11 12 13
  plotMotifs(pfms.tfs[indices])  # chose   11

  indices <- which(tfs.with.motifs == "TAL1")
  indices   #  43 44 45
  plotMotifs(pfms.tfs[indices])  # chose  44

  indices <- which(tfs.with.motifs == "BCL6")
  indices   # 17 19
  plotMotifs(pfms.tfs[indices])  # chose 19

  indices <- which(tfs.with.motifs == "IKZF1")
  indices   # 1 2
  plotMotifs(pfms.tfs[indices])  # chose 2

  indices <- which(tfs.with.motifs == "NFATC2")
  indices   # 26 27
  plotMotifs(pfms.tfs[indices])  # chose 26

  indices <- which(tfs.with.motifs == "TFEC")
  indices   #  14 15
  plotMotifs(pfms.tfs[indices])  # chose 14

  indices <- which(tfs.with.motifs == "ZBTB18")
  indices   # 24 25
  plotMotifs(pfms.tfs[indices])  # chose 24

  chosen <- c(34, 5, 16, 21, 28, 41, 33, 8, 11, 44, 19, 2, 26, 14, 24)
  pfms.chosen <- pfms.tfs[chosen]
  tfs.with.single.motifs <- mcols(pfms.chosen)$geneSymbol
  setdiff(tfs, tfs.with.single.motifs)  #  "PLEK" "LYL1" "VAV1"

  export(pfms.chosen, "~/github/fimoService/pfms/trem2.curated.15.pfms.meme", 'meme')

} # identifyMotifs
#------------------------------------------------------------------------------------------------------------------------
calculateDisruptions <- function()
{
    # extract tfs from the models, then set up fimo with their motifs
#  load("../shinyApp/trem2-models.RData")
#  tbl.summary <- all.models$summary$model
#  tfs <- subset(tbl.summary, observed > 1 & rank.sum < 50)$TF
#  length(tfs)
#
#  pfms.human.jhs <- query(MotifDb, "sapiens", c("jaspar2018", "hocomoco", "swissregulon"))
#  #pfms.human.jhs <- query(MotifDb, "sapiens", c("jaspar2018", "hocomoco"))
#  length(pfms.human.jhs)
#  pfms.tfs <- query(pfms.human.jhs, "", tfs)
#  length(pfms.tfs)
#  dimers <-  grep("::", names(pfms.tfs))
#  if(length(dimers) > 0)
#     pfms.tfs <- pfms.tfs[-dimers]
#  length(pfms.tfs)
#  pfm.names <- rev(sort(names(pfms.tfs)))
#  names(pfm.names) <- pfm.names
#  pfm.names <- sub("\\.[1-9]$", "", pfm.names)
#  dups <- which(duplicated(as.character(pfm.names)))
#  length(dups)
#  if(length(dups) > 0)
#      pfm.names <- pfm.names[-dups]
#  pfm.names.unique <- names(pfm.names)
#  pfms.tfs <- MotifDb[names(pfm.names)]
#  length(pfms.tfs)
#   # take special steps to get tfclass LYL1 motifs
#  #pfms.lyl1 <- query(pfms.human.jhs, "", geneToMotif(MotifDb, "LYL1", c("MotifDb", "TFClass"))$motif)
#  #length(pfms.lyl1)
#  #dups <- intersect(names(pfms.tfs), names(pfms.lyl1))
#  #if(length(dups) > 0){
#  #   deleters <- match(dups, names(pfms.lyl1))
#  #   pfms.lyl1 <- pfms.lyl1[-deleters]
#  #   }
#  #pfms.oi <- c(pfms.tfs, pfms.lyl1)
#  pfms.oi <- pfms.tfs
   load("pfms.chosen.RData")
  pfms.oi <- pfms.chosen
  length(pfms.oi)

  #export(pfms.oi, "~/github/fimoService/pfms/trem2.pfms.meme", 'meme')

     # restart fimo service:
     #   cd ~/github/fimoService/server
     #   add this to makefile.pshannon:
     # MOTIFS="../pfms/trem2.pfms.meme"
     # make -f makefile.pshannon

     # find all fimo matches in the 318kb region defined by minerva: trem family genes +/- 100kb

   mm <- MotifMatcher(genome, as.list(pfms.oi), quiet=TRUE)
   tbl.seq <- getSequence(mm, tbl.regionMax)
   sequence <- tbl.seq$seq  # 318k
   stopifnot(nchar(sequence) == (1 + tbl.regionMax$end - tbl.regionMax$start))
   tbl.fimo.raw <- requestMatch(fimo, list(big=sequence))
   save(tbl.fimo.raw, file="tbl.fimo.raw.RData")
   tbl.fimo <- tbl.fimo.raw
   tbl.fimo$start <- tbl.fimo$start + tbl.regionMax$start
   tbl.fimo$stop <-  tbl.fimo$stop + tbl.regionMax$start
   dim(tbl.fimo)
   tbl.fimo <- subset(tbl.fimo, p.value <= 0.001)
   dim(tbl.fimo)
   tbl.fimo$chrom <- "chr6"
      # add the tf column
   motif.names <- names(pfms.chosen)
   motif2tf.map <- mcols(pfms.chosen)$geneSymbol
   names(motif2tf.map) <- motif.names
   save(motif2tf.map, file="motif2tf.map.RData")
   tfs <- as.character(motif2tf.map[tbl.fimo$motif])
   checkEquals(length(tfs), nrow(tbl.fimo))
   tbl.fimo$tf <- tfs
   dim(tbl.fimo)
   tbl.fimo <- tbl.fimo[, c("chrom", "start", "stop", "strand", "motif", "tf", "score", "p.value", "q.value", "matched.sequence")]
   colnames(tbl.fimo) <- c("chrom", "motifStart", "motifEnd", "strand", "motifName", "tf", "motifScore", "motif.pVal", "motif.qVal", "sequence")
   tbl.fimo <- tbl.fimo[order(tbl.fimo$motifStart, decreasing=FALSE),]
   #motif.names <- unique(tbl.fimo$motifName)
   #length(motif.names)
   #tbl.fimo <- associateTranscriptionFactors(MotifDb, tbl.fimo, source="MotifDb", expand.rows=TRUE)
   #tbl.fimo <- tbl.fimo[, c(2:5, 1, 6:13)]
   #dim(tbl.fimo)
      # find the union of gwas (tbl.snp) and eqtl(tbl.eqtl) snps
   tbl.eqtl.slim <- tbl.eqtl[, c("chrom", "start", "end", "rsid", "ref", "mut", "snpPval")]
   #tbl.eqtl.slim$source <- "eqtl"
   colnames(tbl.eqtl.slim)[7] <- "eqtlPval"
   rownames(tbl.eqtl.slim) <- NULL
   tbl.snp.slim  <-  tbl.snp[, c("chrom", "start", "end", "rsid", "wt", "mut", "snpPval")]
   colnames(tbl.snp.slim)[5] <- "ref"
   colnames(tbl.snp.slim)[7] <- "gwasPval"
   rownames(tbl.snp.slim) <- NULL
   #tbl.snp.slim$source <- "gwas"

   tbl.snpsAll <- merge(tbl.snp.slim, tbl.eqtl.slim, by=c("rsid", "chrom", "start", "end"), all=TRUE)

   tbl.snpsAll <- tbl.snpsAll[, c(2,3,4,1,5:ncol(tbl.snpsAll))]
   colnames(tbl.snpsAll)[c(5,6,8,9)] <- c("ref.gwas", "mut.gwas", "ref.eqtl", "ref.eqtl")
   coi <- c("chrom", "start", "end", "rsid", "ref.gwas", "mut.gwas", "gwasPval", "ref.eqtl", "ref.eqtl", "eqtlPval")
   tbl.snpsAll <- tbl.snpsAll[, coi]
   track <- DataFrameAnnotationTrack("snps", tbl.snpsAll[, c("chrom", "start", "end", "rsid")], color="red")
   displayTrack(igv, track)
   gr.snpsAll <- GRanges(tbl.snpsAll[, 1:3])
   tbl.fimo.locsOnly <- tbl.fimo[, 1:3]
   #colnames(tbl.fimo.locsOnly) <- c("chrom", "start", "end")
   gr.fimo <- GRanges(tbl.fimo.locsOnly)
   tbl.ov <- as.data.frame(findOverlaps(gr.snpsAll, gr.fimo))
   nrow(tbl.ov)
   colnames(tbl.ov) <- c("snp", "fimo")

   tbl.bs <- cbind(tbl.fimo[tbl.ov$fimo,], tbl.snpsAll[tbl.ov$snp,])
     # remove duplicate chrom column, collapse start end down to snpLoc
   tbl.bs <- tbl.bs[, -c(11,13)]   # get rid of extra chrom.1 column, and snp end
   colnames(tbl.bs)[11] <- "snpLoc"
   colnames(tbl.bs)[c(2,3)] <- c("start", "end")
   dim(tbl.bs)   # 18, 18
     # track <- DataFrameAnnotationTrack("ov", tbl.bs[, c("chrom", "start", "stop", "motif")], color="green")
     # displayTrack(igv, track)
   save(tbl.bs, file="tbl.bs.inProcess.RData")

   sequence.mut <- vector(mode="character", length=nrow(tbl.bs))
     # not sure where it creeps in, but our chrom locs are +1 compared to hg38 reference
     # evident, for example, if you blat
     # 17244 chr14 92274808 92274829    RREB1 GGCACGGGGCTGGCTCTGGGGC

     # now add the actual, full mutant sequence for each region

   rc <- function(seq) as.character(reverseComplement(DNAString(seq)))
   tbl.regions <- tbl.bs[, c("chrom", "start", "end")]
   tbl.regions$start <- tbl.regions$start - 1
   tbl.regions$end <- tbl.regions$end - 1

   for(i in seq_len(nrow(tbl.bs))){
      strand <- tbl.bs$strand[i]
      seq.wt <-  getSequence(mm, tbl.regions[i,])$seq
      seq.mut <- tryCatch({getSequence(mm, tbl.regions[i,], tbl.bs$rsid[i])}$seq, error=function(e) return("N"))
      if(strand == "-"){
        seq.wt <- rc(seq.wt)
        seq.mut <- rc(seq.mut)
        }
      stopifnot(seq.wt == tbl.bs$matched.sequence[i])
      printf("%2d.  wt: %s   strand: %s   rsid: %s", i, tbl.bs$sequence[i], tbl.bs$strand[i], tbl.bs$rsid[i])
      printf("    mut: %s", seq.mut)
      #browser()
      sequence.mut[i] <- seq.mut
      } # for i

   tbl.bs$sequence.mut <- sequence.mut
   save(tbl.bs, file="tbl.bindingSitesWithSnps.RData")


} # calculateDisrutpions
#------------------------------------------------------------------------------------------------------------------------
addFimoScores <- function()
{
   load("tbl.bindingSitesWithSnps.RData")
   deltas <- vector(mode="numeric", length=nrow(tbl.bs))
   deltas <- rep(NA_real_, nrow(tbl.bs))

   for(i in seq_len(nrow(tbl.bs))){
      if(tbl.bs$sequence.mut[i] == "N")
         next
      seqs <- list(wt=tbl.bs$sequence[i],  mut=tbl.bs$sequence.mut[i])
      tbl.match <- requestMatch(fimo, seqs)
      tbl.moi <- subset(tbl.match, motif==tbl.bs$motifName[i])  # motif of interest
      scores <- -log10(tbl.moi$p.value)
      delta <- scores[1] - scores[2]
      deltas[i] <- delta
      print(tbl.bs[i, c("chrom", "start")])
      printf("   delta: %5.2f", delta)
      } # for i

   tbl.bs$wtMutMatchDelta <- deltas
   save(tbl.bs, file="../shinyApp/tbl.bindingSitesWithSnpsAndDeltaScores.RData")
   browser()
   xyz <- "end of addFimoScores"

} # addFimoScores

#------------------------------------------------------------------------------------------------------------------------
findFimoMatchesInEnhancers <- function()
{
   load("motif2tf.map.RData")
   mm <- MotifMatcher(genome, as.list(pfms.oi), quiet=TRUE)
   tbl.fimo <- data.frame()

   for(r in 1:nrow(tbl.enhancers.trem2)){
     tbl.seq <- getSequence(mm, tbl.enhancers.trem2[r,])
     sequence <- tbl.seq$seq  # 318k
     tbl.fimo.raw <- requestMatch(fimo, list(big=sequence))
     tf <- as.character(motif2tf.map[tbl.fimo.raw$motif])
     tbl.fimo.raw$tf <- tf
     tbl.fimo.raw$start <- tbl.fimo.raw$start + tbl.enhancers.trem2[r, "start"]
     tbl.fimo.raw$stop <- tbl.fimo.raw$stop + tbl.enhancers.trem2[r, "start"]
     tbl.fimo <- rbind(tbl.fimo, tbl.fimo.raw)
     } # for r

   coi <- c("motif", "start", "stop", "strand", "score", "p.value",  "q.value", "matched.sequence", "tf")
   tbl.fimo <- tbl.fimo[, coi]
   colnames(tbl.fimo)[3] <- "end"
   dim(tbl.fimo)                             # 53832 9
   dim(subset(tbl.fimo, p.value <= 1e-2))    #  7259
   dim(subset(tbl.fimo, p.value <= 1e-3))    #  970
   dim(subset(tbl.fimo, p.value <= 1e-4))    #  117
   tbl.fimoHitsInEnhancers <- tbl.fimo
   save(tbl.fimoHitsInEnhancers, file="../shinyApp/tbl.fimoHitsInEnhancers.RData")

} # findFimoMatchesInEnhancers
#------------------------------------------------------------------------------------------------------------------------
run <- function()
{
   models.ly11 <- calculate(sgm, list(find.lyl1=find.lyl1.spec(mtx.cer)))

   build.spec <- list(title="fp.enhancers",
                      type="footprint.database",
                      regions=data.frame(chrom="chr6", start=tss-5000, end=tss+5000),
                      tss=tss,
                      matrix=mtx,
                      db.host="khaleesi.systemsbiology.net",
                      databases=c("brain_hint_20"),  #"brain_hint_16", "brain_wellington_20", "brain_wellington_16"),
                      motifDiscovery="builtinFimo",
                      tfMapping=c("MotifDb"),
                      tfPrefilterCorrelation=0.4,
                      orderModelByColumn="rfScore",
                      solverNames=c("lasso", "lassopv", "pearson", "randomForest", "ridge", "spearman"))

   fpBuilder <- FootprintDatabaseModelBuilder("hg38", "TREM2", build.spec, quiet=FALSE)
   x <- build(fpBuilder)

   tbl.reg <- x$regulatoryRegions
   tbl.model <- x$model


   model.lyl1 <-  calculate(sgm, list(lyl1.10k=build.spec))

   models.01.cer <- calculate(sgm, list(fp.2000up.200down.cor02=fp.2000up.200down.cor02(mtx.cer)))
   models.01.tcx <- calculate(sgm, list(fp.2000up.200down.cor02=fp.2000up.200down.cor02(mtx.tcx)))

   models.02.cer <- calculate(sgm, list(fp.5kup.5kdown.cor02=fp.5kup.5kdown.cor02(mtx.cer)))
   models.02.tcx <- calculate(sgm, list(fp.5kup.5kdown.cor02=fp.5kup.5kdown.cor02(mtx.tcx)))

   models.03.cer <- calculate(sgm, list(fp.enhancers.cor02=fp.enhancers.cor02(tbl.enhancers.trem2, mtx.cer)))
   models.03.tcx <- calculate(sgm, list(fp.enhancers.cor02=fp.enhancers.cor02(tbl.enhancers.trem2, mtx.tcx)))

   #models.04.cer <- calculate(sgm, list(fp.enhancers.cor02=fp.enhancers.cor02(tbl.enhancers.trem2, mtx.cer)))
   #models.04.tcx <- calculate(sgm, list(fp.enhancers.cor02=fp.enhancers.cor02(tbl.enhancers.trem2, mtx.tcx)))

   models.05.cer <- calculate(sgm, list(noDNA.allTFs.cor04=noDNA.allTFs.cor04(mtx.cer)))
   models.05.tcx <- calculate(sgm, list(noDNA.allTFs.cor04=noDNA.allTFs.cor04(mtx.tcx)))


   all.models <- list(fp.2000up.200down.cer=models.01.cer[[1]],
                      fp.2000up.200down.tcx=models.01.tcx[[1]],
                      fp.5kup.5kdown.cer=models.03.cer[[1]],
                      fp.5kup.5kdown.tcx=models.03.tcx[[1]],
                      #fp.enhancers.cer=models.04.cer[[1]],
                      #fp.enhancers.txc=models.04.tcx[[1]],
                      allTFs.cer=models.05.cer[[1]],
                      allTFs.tcx=models.05.tcx[[1]]
                      )


   for(i in seq_len(length(all.models))){
      tbl <- all.models[[i]]$model
      printf("%s: %d %d", names(all.models)[i], nrow(tbl), ncol(tbl))
      tbl2 <- roundNumericColumns(tbl, 3, "lassoPValue")
      all.models[[i]]$model <- tbl2
      print(head(tbl2, n=3))
      }

   tbl.summary <- summarizeModels(sgm, orderBy="rfScore", maxTFpredictors=10, all.models)
   tbl.summary
   tbl.summary$TF <- rownames(tbl.summary)
   tbl.summary <- tbl.summary[, c("TF", colnames(tbl.summary)[1:ncol(tbl.summary)-1])]
   rownames(tbl.summary) <- NULL
   all.models[["summary"]] <- list(model=tbl.summary, regulatoryRegions=data.frame())

   save(all.models, file="../shinyApp/trem2-models.RData")
   save(tbl.enhancers.trem2, tbl.enhancers.treml1, file="../shinyApp/enhancers.RData")
   save(tbl.snp, tbl.eqtl, file="../shinyApp/snps.RData")

} # run
#------------------------------------------------------------------------------------------------------------------------
