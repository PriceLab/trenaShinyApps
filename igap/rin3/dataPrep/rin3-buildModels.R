library(trenaSGM)
library(igvR)
library(motifStack)
library(motifbreakR)
library(SNPlocs.Hsapiens.dbSNP150.GRCh38)
library(BSgenome.Hsapiens.UCSC.hg38)

  #   pcm <- new("pcm", mat=pfms[[1]], name=names(pfms)); plot(pcm)
  #   motifStack(lapply(names(x), function(mName) new("pfm", x[[mName]], name=mName)))

source("~/github/trenaShinyApps/utils/getEnhancers.R")
library(FimoClient)
library(GenomicRanges)
#------------------------------------------------------------------------------------------------------------------------
FIMO_HOST <- "localhost"
FIMO_PORT <- 5558
if(!exists("fimo")){
   fimo <- FimoClient(FIMO_HOST, FIMO_PORT, quiet=FALSE)
   }
#------------------------------------------------------------------------------------------------------------------------
# https://escholarship.umassmed.edu/cgi/viewcontent.cgi?article=1938&context=gsbs_diss: rreb1 in neuronal regenertion
#------------------------------------------------------------------------------------------------------------------------
genome <- "hg38"
targetGene <- "RIN3"
targetGeneID <- "79890"
chromosome <- "chr14"
tss <- 92513774
#------------------------------------------------------------------------------------------------------------------------
if(!exists("mtx"))
   load(system.file(package="trenaSGM", "extdata", "mayo.tcx.RData"))

if(!exists("sgm"))
   sgm <- trenaSGM(genome, targetGene, quiet=FALSE)

if(!exists("tbl.enhancers"))
   tbl.enhancers <- getEnhancerRegions(targetGene)

if(!exists("igv")){
   igv <- igvR()
   Sys.sleep(5)
   setGenome(igv, genome)
   Sys.sleep(5)
   showGenomicRegion(igv, targetGene)
   }
#------------------------------------------------------------------------------------------------------------------------
proximal.promoter.fp <- function()
{
   tbl.regions <- tbl.enhancers[42,]  # by direct inspection.  10119 bp straddling TSS

   spec <- list(title="proximal.promoter.fp",
                type="footprint.database",
                regions=tbl.regions,
                tss=tss,
                matrix=mtx,
                db.host="khaleesi.systemsbiology.net",
                databases=list("brain_hint_20", "brain_hint_16", "brain_wellington_20", "brain_wellington_16"),
                motifDiscovery="builtinFimo",
                tfMapping="MotifDB",
                tfPrefilterCorrelation=0.4,
                orderModelByColumn="rfScore",
                solverNames=c("lasso", "lassopv", "pearson", "randomForest", "ridge", "spearman"))

   return(spec)

} # proximal.promoter.fp
#------------------------------------------------------------------------------------------------------------------------
all.enhancers.fp <- function()
{
   tbl.regions <- tbl.enhancers  # by direct inspection.  10119 bp straddling TSS

   spec <- list(title="all.enhancers.fp",
                type="footprint.database",
                regions=tbl.regions,
                tss=tss,
                matrix=mtx,
                db.host="khaleesi.systemsbiology.net",
                databases=list("brain_hint_20", "brain_hint_16", "brain_wellington_20", "brain_wellington_16"),
                motifDiscovery="builtinFimo",
                tfMapping="MotifDB",
                tfPrefilterCorrelation=0.4,
                orderModelByColumn="rfScore",
                solverNames=c("lasso", "lassopv", "pearson", "randomForest", "ridge", "spearman"))

   return(spec)

} # all.enhancers.fp
#------------------------------------------------------------------------------------------------------------------------
proximal.promoter.pwmMatch <- function()
{
   tbl.regions <- tbl.enhancers[42,]  # by direct inspection.  10119 bp straddling TSS

   spec <- list(title="proximal.promoter.pwmMatch",
                type="regions.motifMatching",
                tss=tss,
                regions=tbl.regions,
                matrix=mtx,
                pfms=query(MotifDb, "sapiens", c("jaspar2018", "hocomoco")),
                matchThreshold=80,
                motifDiscovery="matchPWM",
                tfMapping=c("MotifDb", "TFClass"),
                tfPrefilterCorrelation=0.4,
                orderModelByColumn="rfScore",
                solverNames=c("lasso", "lassopv", "pearson", "randomForest", "ridge", "spearman"))

   return(spec)

} # proximal.promoter.pwmMatch
#------------------------------------------------------------------------------------------------------------------------
# rreb1 is the top tf using footprints, but not found by pwmMatch.  why not?
# see ~/github/fimoService/client-R/FimoClient/inst/unitTests/test_FimoClient.R and ~/github/notes/log entry
# titled * igap survey: rin3 & and the RREB1 transcription factor (27 may 2018)
study.rreb1.pwmatch <- function()
{
  pfms=query(MotifDb, c("sapiens", "MA0073"), c("jaspar2018", "hocomoco"))
  mm <- MotifMatcher(genome, as.list(pfms), quiet=TRUE)
  tbl.motifs <- findMatchesByChromosomalRegion(mm, tbl.enhancers[42,], pwmMatchMinimumAsPercentage=75)

  tbl.rreb1.fp <- unique(subset(models.01[[1]]$regulatoryRegions, geneSymbol=="RREB1")[, c("chrom", "fp_start", "fp_end")])
  tbl.rreb1.pwmMatch <- unique(tbl.motifs[, c("chrom", "motifStart", "motifEnd")])

  track.fp <- DataFrameAnnotationTrack("rreb1.fp", tbl.rreb1.fp, color="blue")
  displayTrack(igv, track.fp)
  track.pwm <- DataFrameAnnotationTrack("rreb1.pwm", tbl.rreb1.pwmMatch, color="green")
  displayTrack(igv, track.pwm)

  # 33 bp highlights the difference: chr14:92,513,603-92,513,635
  # getSequence(mm, data.frame(chrom="chr14", start=92513603, end=92513635))
  #   chrom    start      end                               seq status
  #   chr14 92513603 92513635   CTTGGCCCCAGCACCCCCCGCCCCGAGGCCCGG     wt
  #   consensusString(x[[4]])   ] "CCCCAAACCACCCCC??CC?"
  #                                        ||
  x <- query(MotifDb, c("sapiens", "MA0073", "jaspar"))
  pcm <- new("pcm", mat=x[[3]], name=names(x)[3])
  plot(pcm)

} # study.rreb1.pwmatch
#------------------------------------------------------------------------------------------------------------------------
find_RREB1_bindingSites <- function()
{
    pfms <- query(MotifDb, c("sapiens", "RREB1"))
    mm <- MotifMatcher(genome, as.list(pfms), quiet=TRUE)
    tbl.bigRegion <- data.frame(chrom="chr14",
                                start=min(tbl.enhancers$start),
                                end=max(tbl.enhancers$end),
                                stringsAsFactors=FALSE)
    tbl.motifs <- findMatchesByChromosomalRegion(mm, tbl.bigRegion, pwmMatchMinimumAsPercentage=75)

    fivenum(tbl.motifs$motifScore)
    fivenum(tbl.motifs$motifRelativeScore)

    tbl.motifs.oi <- subset(tbl.motifs, motifScore > 10)
    dim(tbl.motifs.oi)
    tbl.motifs.oi <- tbl.motifs.oi[, c("chrom", "motifStart", "motifEnd", "motifScore")]
    tbl.motifs.oi <- tbl.motifs.oi[order(tbl.motifs.oi$motifScore, decreasing=TRUE),]
    dups <- which(duplicated(tbl.motifs.oi[, c("chrom", "motifStart", "motifEnd")]))
    if(length(dups) > 0)
       tbl.motifs.oi <- tbl.motifs.oi[-dups, ]
    dim(tbl.motifs.oi)
    track.motifs.rreb1 <- DataFrameAnnotationTrack("mm RREB1", tbl.motifs.oi, color="green")
    displayTrack(igv, track.motifs.rreb1)


    tbl.seq <- getSequence(mm, tbl.bigRegion)
    sequence <- list(big=tbl.seq[1,"seq"])
    stopifnot(nchar(sequence) > 1356530)
    tbl.fimo <- requestMatch(fimo, sequence)
    fivenum(tbl.fimo$score)

    tbl.fimo.oi <- subset(tbl.fimo, score > 7)[, c("start", "stop", "score")]
    tbl.fimo.oi$chrom <- "chr14"
    tbl.fimo.oi$start <- tbl.fimo.oi$start + tbl.bigRegion$start[1]
    tbl.fimo.oi$stop <- tbl.fimo.oi$stop + tbl.bigRegion$start[1]

    tbl.fimo.oi <- tbl.fimo.oi[, c("chrom", "start", "stop", "score")]
    colnames(tbl.fimo.oi)[3] <- "end"
    tbl.fimo.oi <- tbl.fimo.oi[order(tbl.fimo.oi$score, decreasing=FALSE),]
    dups <- which(duplicated(tbl.fimo.oi[, c("chrom", "start", "end")]))
    printf("tbl.fimo.oi dups: %d", length(dups))
    if(length(dups) > 0)
       tbl.fimo.oi <- tbl.fimo.oi[-dups,]

    track.fimo.rreb1 <- DataFrameAnnotationTrack("fimo RREB1", tbl.fimo.oi, color="magenta")
    displayTrack(igv, track.fimo.rreb1)

    track.enhancers <- DataFrameAnnotationTrack("enhancers", tbl.enhancers, color="black")
    displayTrack(igv, track.enhancers)

    print(load("~/s/work/priceLab/AD/tbl.gwas.level_1.RData")) # tbl.gwas
    tbl.snp <- subset(tbl.gwas, BP >= tbl.bigRegion$start & BP <= tbl.bigRegion$end & CHR=="chr14")
    tbl.snp <- tbl.snp[, c("CHR", "BP", "BP", "SNP")]
    colnames(tbl.snp) <- c("chrom", "start", "end", "rsid")

    track.snp <- DataFrameAnnotationTrack("snp", tbl.snp, color="red")
    displayTrack(igv, track.snp)

    shoulder <- 10

    tbl.fimo.oi.shoulder <- tbl.fimo.oi
    tbl.fimo.oi.shoulder$start <- tbl.fimo.oi.shoulder$start - shoulder
    tbl.fimo.oi.shoulder$end   <- tbl.fimo.oi.shoulder$end   + shoulder

    tbl.motifs.oi.shoulder <- tbl.motifs.oi
    tbl.motifs.oi.shoulder$motifStart <- tbl.motifs.oi.shoulder$motifStart - shoulder
    tbl.motifs.oi.shoulder$motifEnd   <- tbl.motifs.oi.shoulder$motifEnd   + shoulder
    colnames(tbl.motifs.oi.shoulder) <- c("chrom", "start", "end", "score")

    gr.snps <- GRanges(tbl.snp)
    gr.fimo <- GRanges(tbl.fimo.oi.shoulder)
    gr.mm   <- GRanges(tbl.motifs.oi.shoulder)
    tbl.ov.fimo <- as.data.frame(findOverlaps(gr.snps, gr.fimo))
    colnames(tbl.ov.fimo) <- c("snp", "fimo")
    tbl.ov.mm   <- as.data.frame(findOverlaps(gr.snps, gr.mm))
    colnames(tbl.ov.mm) <- c("snp", "mm")

    snp.hits <- sort(unique(c(tbl.ov.fimo$snp, tbl.ov.mm$snp)))
    tbl.snp.hit <- tbl.snp[snp.hits,]
    track.snp.hit <- DataFrameAnnotationTrack("snpHit", tbl.snp.hit, color="red")
    displayTrack(igv, track.snp.hit)

} # find_RREB1_bindingSites
#------------------------------------------------------------------------------------------------------------------------
fp.5kup.5kdown.cor02 <- function()
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
                tfMapping="MotifDB",
                tfPrefilterCorrelation=0.2,
                orderModelByColumn="rfScore",
                solverNames=c("lasso", "lassopv", "pearson", "randomForest", "ridge", "spearman"))

   return(spec)

} # fp.5kup.5kdown.cor02
#------------------------------------------------------------------------------------------------------------------------
fp.enhancers.cor02 <- function()
{
   tbl.regions <- getEnhancerRegions(targetGene)
   spec <- list(title="fp.enhancers.02",
                type="footprint.database",
                regions=tbl.regions,
                tss=tss,
                matrix=mtx,
                db.host="khaleesi.systemsbiology.net",
                databases=list("brain_hint_20", "brain_hint_16", "brain_wellington_20", "brain_wellington_16"),
                motifDiscovery="builtinFimo",
                tfMapping="MotifDB",
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
                tfMapping="MotifDB",
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
                tfMapping=c("MotifDb", "TFClass"),
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
                tfMapping=c("MotifDb", "TFClass"),
                tfPrefilterCorrelation=0.2,
                orderModelByColumn="rfScore",
                solverNames=c("lasso", "lassopv", "pearson", "randomForest", "ridge", "spearman"))


   return(spec)

} # regions.motifMatching.2kup.200down.pwm80.cor02
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
                tfMapping=c("MotifDb", "TFClass"),
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
                tfMapping=c("MotifDb", "TFClass"),
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
                tfMapping=c("MotifDb", "TFClass"),
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
                tfMapping=c("MotifDb", "TFClass"),
                tfPrefilterCorrelation=0.4,
                orderModelByColumn="rfScore",
                solverNames=c("lasso", "lassopv", "pearson", "randomForest", "ridge", "spearman"))

   return(spec)

} # regions.motifMatching.snpEnhancerOnly
#------------------------------------------------------------------------------------------------------------------------
noDNA.allTFs.cor04 <- function()
{
   tbl.regions <- data.frame(chrom=chromosome, start=tss-200, end=tss+2000, stringsAsFactors=FALSE)

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
run <- function()
{
   models.01 <- calculate(sgm, list(proximal.promoter.fp=proximal.promoter.fp()))
   tbl.fpInModel <- subset(models.01[[1]]$regulatoryRegions, geneSymbol %in% models.01[[1]]$model$gene[1:5])
   tbl.fpInModel <- tbl.fpInModel[, c("chrom", "fp_start", "fp_end", "geneSymbol")]
   colnames(tbl.fpInModel) <- c("chrom", "start", "end", "geneSymbol")

   models.03 <- calculate(sgm, list(all.enhancers.fp=all.enhancers.fp()))
   tbl.fpInModel <- subset(models.01[[1]]$regulatoryRegions, geneSymbol %in% models.01[[1]]$model$gene[1:5])
   tbl.fpInModel <- tbl.fpInModel[, c("chrom", "fp_start", "fp_end", "geneSymbol")]
   colnames(tbl.fpInModel) <- c("chrom", "start", "end", "geneSymbol")

   model.04 <- calculate(sgm, list(allPossibleTFs=noDNA.allTFs.cor04()))

   tbl.seq <- getSequence(mm, tbl.bigRegion)

   sequence <- list(big=tbl.seq[1,"seq"])
   stopifnot(nchar(sequence) > 1356530)
      # restart fimo server with pfms from the top 10 tfs
   tfs <- models.03[[1]]$model$gene[1:10]
   pfms <- query(MotifDb, c("sapiens"), c(tfs))
   pfms <- query(pfms, "",  c("jaspar2018", "hocomoco"))
   export(pfms, "~/github/fimoService/pfms/model03.pfms.meme", 'meme')
    # restart fimo service, load these 22 matrices
   tbl.fimo <- requestMatch(fimo, sequence)  # 12525 9
   tbl.fimo$start <- tbl.fimo$start + tbl.bigRegion$start
   tbl.fimo$stop <- tbl.fimo$stop + tbl.bigRegion$start
   tbl.fimo$chrom <- "chr14"
   tbl.fimo.hits <- tbl.fimo[, c("chrom", "start", "stop", "motif", "q.value", "score")]
   colnames(tbl.fimo.hits)[3] <- "end"
   tbl.fimo.hits <- subset(tbl.fimo.hits, q.value < 0.1)
   gr.fimo <- GRanges(tbl.fimo.hits)
   print(load("~/s/work/priceLab/AD/tbl.gwas.level_1.RData")) # tbl.gwas
   tbl.snp <- subset(tbl.gwas, BP >= tbl.bigRegion$start & BP <= tbl.bigRegion$end & CHR=="chr14")
   tbl.snp <- tbl.snp[, c("CHR", "BP", "BP", "SNP")]
   colnames(tbl.snp) <- c("chrom", "start", "end", "rsid")
   gr.snps <- GRanges(tbl.snp)

   tbl.ov.fimo <- as.data.frame(findOverlaps(gr.snps, gr.fimo))
   colnames(tbl.ov.fimo) <- c("snp", "fimo")
   tbl.snp.hits <- tbl.snp[tbl.ov.fimo$snp,]
   tbl.fimo.tf.hits  <- tbl.fimo.hits[tbl.ov.fimo$fimo,]
   dim(tbl.snp.hits)
   dim(tbl.fimo.tf.hits)
   save(tbl.snp.hits, tbl.fimo, tbl.fimo.tf.hits, file="forDisruptionStudy.RData")


   track.snp.hit <- DataFrameAnnotationTrack("snpHit", tbl.snp.hits, color="red")
   displayTrack(igv, track.snp.hit)

   track.fimo.hit <- DataFrameAnnotationTrack("fimoHit", tbl.fimo.tf.hits, color="purple")
   displayTrack(igv, track.fimo.hit)


   xyz <- 99
#
#    print(load("~/s/work/priceLab/AD/tbl.gwas.level_1.RData")) # tbl.gwas
#    tbl.snp <- subset(tbl.gwas, BP >= tbl.bigRegion$start & BP <= tbl.bigRegion$end & CHR=="chr14")
#    tbl.snp <- tbl.snp[, c("CHR", "BP", "BP", "SNP")]
#    colnames(tbl.snp) <- c("chrom", "start", "end", "rsid")
#
#    gr.snps <- GRanges(tbl.snp)
#    shoulder <- 10
#    tbl.fpInModelShoulder <- tbl.fpInModel
#    tbl.fpInModelShoulder$start <- tbl.fpInModelShoulder$start - shoulder
#    tbl.fpInModelShoulder$end <- tbl.fpInModelShoulder$end + shoulder
#    gr.fpInModel <- GRanges(tbl.fpInModelShoulder)
#
#    tbl.ov.fp <- as.data.frame(findOverlaps(gr.snps, gr.fpInModel))
#    dim(tbl.ov.fp)
#    if(nrow(tbl.ov.fp) > 0){
#       colnames(tbl.ov.fimo) <- c("snp", "fimo")
#       tbl.ov.mm   <- as.data.frame(findOverlaps(gr.snps, gr.mm))
#       colnames(tbl.ov.mm) <- c("snp", "mm")
#       snp.hits <- sort(unique(c(tbl.ov.fimo$snp, tbl.ov.mm$snp)))
#       }
#
#
#    track.fp01 <- DataFrameAnnotationTrack("fp01", tbl.fpInModel(
#    models.02 <- calculate(sgm, list(proximal.promoter.pwmMatch=proximal.promoter.pwmMatch()))
#      # models.02[[1]]$model finds no RREB1, MA0073.1
#    all.models <- c(models.01, models.02)
#
#
#
#    model.small <- calculate(sgm, list(regions.motifMatching.snpEnhancerOnly=regions.motifMatching.snpEnhancerOnly()))
#
#    tbl.summary <- summarizeModels(sgm, orderBy="rfScore", maxTFpredictors=10, all.models)
#    top.tfs <- rownames(tbl.summary)
#
#    for(i in seq_len(length(all.models))){
#       tbl <- all.models[[i]]$model
#       printf("%s: %d %d", names(all.models)[i], nrow(tbl), ncol(tbl))
#       tbl2 <- roundNumericColumns(tbl, 3, "lassoPValue")
#       all.models[[i]]$model <- tbl2
#       print(head(tbl2, n=3))
#       }
#
#    save(all.models, tbl.summary, file="../shinyApp/trem2-models.RData")

} # run
#------------------------------------------------------------------------------------------------------------------------
#  head(models.01[[1]]$model)
#    gene  betaLasso  lassoPValue pearsonCoeff   rfScore  betaRidge spearmanCoeff concordance    pcaMax bindingSites
#   RREB1 0.33758371 1.234729e-48    0.7476279 48.319651 0.09419397     0.7432798   0.5418992 5.2227741           90
#    SPI1 0.35870087 1.160131e-49    0.7536468 44.466933 0.11484289     0.7541331   0.5368280 5.5780258          120
#  NFATC2 0.00000000 8.233839e-01    0.6689543 32.781795 0.04383253     0.6943073   0.4476320 1.1427128           44
#    ELK3 0.00000000 8.568806e-01    0.6086289 12.171360 0.02014866     0.6347150   0.3663466 0.8766895           23
#    IRF5 0.04410938 7.679872e-07    0.6046606  9.216642 0.09226013     0.6042614   0.3861056 0.7775808            6
#   STAT3 0.00000000 6.942070e-01    0.6081331  9.120325 0.01604489     0.6348048   0.3582188 0.8675912           55
disruptionStudy <- function()
{
   if(!all(exists(c("tbl.snp.hits", "tbl.fimo.tf.hits", "tbl.fimo"))))
      load("forDisruptionStudy.RData")

   if(!exists("tbl.gwas"))
      load("~/s/work/priceLab/AD/tbl.gwas.level_1.RData")

   rsids <-tbl.snp.hits$rsid
   tbl.gwas.oi <- subset(tbl.gwas, SNP %in% rsids)
   tbl.gwas.oi$score <- -log10(tbl.gwas.oi$P)
   tbl.gwas.oi <- tbl.gwas.oi[order(tbl.gwas.oi$score, decreasing=TRUE),]
   rsid.oi <- "rs10130041"
   tbl.fimo.tf.hits[grep(rsid.oi, tbl.snp.hits$rsid),] # [1] 60
   x <- MotifDb["Hsapiens-jaspar2018-RREB1-MA0073.1"]
   plot(new("pcm", mat=x[[1]], name=names(x)))
   consensusString(x[[1]])                                                           # CCCCAAACCACCCCC??CC?
   getSequence(mm, data.frame(chrom="chr14", start=92465743, end=92465762))          # CCCCTCCCGACCCCCGAGAG
   getSequence(mm, data.frame(chrom="chr14", start=92465743, end=92465762), rsid.oi) # CCCCTCCCGACCCCCAAGAG
   requestMatch(fimo, list(wt="CCCCTCCCGACCCCCGAGAG", mut="CCCCTCCCGACCCCCAAGAG", con="CCCCAAACCACCCCCCCCCC"))
    # just one match
    # 1 Hsapiens-jaspar2018-RREB1-MA0073.1           con     1   20      + 30.7865       0       0 CCCCAAACCACCCCCCCCCC

   requestMatch(fimo, list(wt="GGGGAGCTGGTGTTTTCTGTGG", mut="GGGGAGCTGATGTTTTCTGTGG", "good?"="ACAGAAAACAACAGCTCCCC"))
    #                                motif sequence.name start stop strand     score  p.value  q.value     matched.sequence
    #   Hsapiens-jaspar2018-RREB1-MA0073.1            wt     1   20      -  9.707870 5.22e-06 5.22e-05 ACAGAAAACACCAGCTCCCC
    #   Hsapiens-jaspar2018-RREB1-MA0073.1         good?     1   20      +  8.617980 7.46e-06 5.22e-05 ACAGAAAACAACAGCTCCCC
    #   Hsapiens-jaspar2018-RREB1-MA0073.1           mut     1   20      - -0.640449 8.68e-05 4.05e-04 ACAGAAAACATCAGCTCCCC

   mm <- MotifMatcher("hg38", as.list(MotifDb["Hsapiens-jaspar2018-RREB1-MA0073.1"]))
   tbl.regions <- data.frame(chrom="chr14", start=92432732, end=92432753)

   findMatchesByChromosomalRegion(mm, tbl.regions, 70) # , variants="rs2402140"))
   #                        motifName chrom motifStart motifEnd strand motifScore motifRelativeScore                match chromStart chromEnd                  seq status shortMotif
   #Hsapiens-jaspar2018-RREB1-MA0073.1 chr14   92432732 92432751      -   10.90909          0.7692308 GGGGAGCTGGTGTTTTCTGT   92432732 92432753 GGGGAGCTGGTGTTTTC...     wt   MA0073.1
   #                                                                                                 rc ACAGAAAACACCAGCTCCCC
   #  Hsapiens-jaspar2018-RREB1-MA0073.1 chr14   92432732 92432751      -   10.27273          0.724359  GGGGAGCTGATGTTTTCTGT   92432732 92432753 GGGGAGCTGATGTTTTC...    mut   MA0073.1
   #                                                                                                 rc ACAGAAAACATCAGCTCCCC

} # disruptionStudy
#------------------------------------------------------------------------------------------------------------------------
dhs <- function()
{
   hdf <- HumanDHSFilter("hg38", "wgEncodeRegDnaseClustered", pwmMatchPercentageThreshold=0,
                         geneInfoDatabase.uri="bogus", regions=data.frame(), pfms=list())
   chrom <- "chr14"
   start <- 91654248
   end <-   93010778
   tbl.dhs <- getRegulatoryRegions(hdf, "wgEncodeRegDnaseClustered", chrom, start, end)
   colnames(tbl.dhs) <- c("chrom", "start", "end", "count", "score")
   save(tbl.dhs, file="../shinyApp/tbl.dhs.RData")

} # dhs
#------------------------------------------------------------------------------------------------------------------------
findAllBindingSitesForAllTFsInSummaryTable <- function()
{
   tfs <- tbl.summary$gene
   mdb.base <- query(MotifDb, "sapiens", c("jaspar2018", "hocomoco", "swissregulon"))
   mdb.tfs <- query(mdb.base, "", tfs)
   setdiff(tfs, mcols(mdb.tfs)$geneSymbol)    # 3 expression-only tfs not found: "AKNA"   "KLHL6"  "SP140L"
   export(mdb.tfs, "~/github/fimoService/pfms/rin3.trena.pfms.meme", 'meme')

   chrom <- "chr14"
   start <- 91654248
   end <-   93010778
   tbl.sequence <- getSequence(mm, data.frame(chrom=chrom, start=start, end=end))
   checkEquals(nchar(tbl.sequence$seq[1]), 1 + end - start)
   tbl.matches <- requestMatch(fimo, list(big=tbl.sequence$seq[1]))
   tf <- mcols(mdb.tfs[tbl.matches$motif])$geneSymbol
   tbl.matches$tf <- tf    # 19832 x 10
   coi <- c("chrom", "start", "stop", "strand", "motif", "score", "p.value", "q.value", "tf", "matched.sequence")
   tbl.matches <- tbl.matches[, coi]
   colnames(tbl.matches)[3] <- "end"
   colnames(tbl.matches)[5] <- "name"
   save(tbl.matches, file="../shinyApp/sequenceMatchesForTfsInModel.RData")

} # findAllBindingSitesForAllTFsInSummaryTable
#------------------------------------------------------------------------------------------------------------------------
test.motifBreaker <- function()
{
   snp <- "rs10130041"
   snps.gr <- snps.from.rsid(rsid = snp,
                            dbSNP=SNPlocs.Hsapiens.dbSNP150.GRCh38,
                            search.genome=BSgenome.Hsapiens.UCSC.hg38)
   results <- motifbreakR(snpList = snps.gr, filterp = TRUE,
                          pwmList = MotifDb["Hsapiens-jaspar2018-RREB1-MA0073.1"],
                          threshold = 0,
                          method = "ic",
                          bkg = c(A=0.25, C=0.25, G=0.25, T=0.25),
                          BPPARAM = BiocParallel::bpparam())

    #  rs2402140
    # chr14:92,432,733-92,432,754
    seq <- getSequence(mm, data.frame(chrom="chr14", start=92432730, end=92432755))
    seq.mut <- getSequence(mm, data.frame(chrom="chr14", start=92432730, end=92432755), variants="rs2402140")
      # seq.mut$seq:    "CTGGGGAGCTGATGTTTTCTGTGGCT"
      # seq$seq:        "CTGGGGAGCTGGTGTTTTCTGTGGCT"

   requestMatch(fimo, list(wt=seq$seq, mut=seq.mut$seq))
     # 1 Hsapiens-HOCOMOCOv10-RREB1_HUMAN.H10MO.D            wt     2   23      + 16.977800 9.17e-07 1.83e-05 TGGGGAGCTGGTGTTTTCTGTG
     # 2 Hsapiens-SwissRegulon-RREB1.SwissRegulon            wt     2   23      + 16.733300 1.01e-06 2.02e-05 TGGGGAGCTGGTGTTTTCTGTG
     # 3       Hsapiens-jaspar2018-RREB1-MA0073.1            wt     3   22      -  9.707870 5.22e-06 1.46e-04   ACAGAAAACACCAGCTCCCC
     # 4       Hsapiens-jaspar2018-RREB1-MA0073.1           mut     3   22      - -0.640449 8.68e-05 1.22e-03   ACAGAAAACATCAGCTCCCC



} # test.motifBreaker
#------------------------------------------------------------------------------------------------------------------------
