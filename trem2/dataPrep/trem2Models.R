ibrary(trenaSGM)
library(motifStack)
#------------------------------------------------------------------------------------------------------------------------
genome <- "hg38"
targetGene <- "TREM2"
sgm <- trenaSGM(genome, targetGene, quiet=FALSE)
chromosome <- "chr6"
tss <- 41163186
#------------------------------------------------------------------------------------------------------------------------
if(!exists("mtx"))
   load(system.file(package="trenaSGM", "extdata", "mayo.tcx.RData"))

if(!exists("sgm"))
   sgm <- trenaSGM(genome, targetGene)

#------------------------------------------------------------------------------------------------------------------------
getEnhancerRegions <- function(targetGene)
{
   base.dir <- "~/s/data/public/GeneHancer/v4.7"
   path <- file.path(base.dir, "enhancer_elite_ids.txt")
   tbl.elite <- read.table(path, sep="\t", as.is=TRUE, header=TRUE)            #  243281      6
   path <- file.path(base.dir, "enhancer_gene_scores.txt")
   tbl.geneScoresAll <- read.table(path, sep="\t", as.is=TRUE, header=TRUE)  #  934287      9
      # dim(tbl.elite); dim(tbl.geneScoresAll)
      # [1]  250733      7
      # [1] 1086552      9

   tbl.scores <- subset(tbl.geneScoresAll, symbol==targetGene)
   clusters <- tbl.scores$cluster_id
   tbl.enhancers <- subset(tbl.elite, cluster_id %in% clusters)[, c("chr", "enhancer_start", "enhancer_end")]
   colnames(tbl.enhancers) <- c("chrom", "start", "end")
   tbl.enhancers$chrom <- paste("chr", tbl.enhancers$chrom, sep="")
   rownames(tbl.enhancers) <- NULL

   tbl.enhancers

} # getEnhancerRegions
#------------------------------------------------------------------------------------------------------------------------
fp.2000up.200down.cor02 <- function()
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
                tfMapping="MotifDB",
                tfPrefilterCorrelation=0.2,
                orderModelByColumn="rfScore",
                solverNames=c("lasso", "lassopv", "pearson", "randomForest", "ridge", "spearman"))

   return(spec)

} # fp.2000up.200down.cor02
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
   models.01 <- calculate(sgm, list(fp.2000up.200down.cor02=fp.2000up.200down.cor02()))
   models.02 <- calculate(sgm, list(fp.5kup.5kdown.cor02=fp.5kup.5kdown.cor02()))
   models.03 <- calculate(sgm, list(fp.enhancers.cor02=fp.enhancers.cor02()))
   models.04 <- calculate(sgm, list(regions.motifMatching.5kup.5kdown.pwm80.cor02=regions.motifMatching.5kup.5kdown.pwm80.cor02()))
   models.05 <- calculate(sgm, list(regions.motifMatching.2kup.200down.pwm80.cor02=regions.motifMatching.2kup.200down.pwm80.cor02()))
   models.06 <- calculate(sgm, list(regions.motifMatching.2kup.200down.pwm90.cor02=regions.motifMatching.2kup.200down.pwm90.cor02()))
   models.07 <- calculate(sgm, list(regions.motifMatching.2kup.200down.pwm95.cor02=regions.motifMatching.2kup.200down.pwm95.cor02()))
   models.08 <- calculate(sgm, list(noDNA.allTFs.cor04=noDNA.allTFs.cor04()))
   models.09 <- calculate(sgm, list(fp.enhancers.cor02.5kup.5kdown=fp.enhancers.cor02.5kup.5kdown()))
   models.10 <- calculate(sgm, list(regions.motifMatching.enhancers.5kbup.5kbdown.pwm80.cor04=
                                       regions.motifMatching.enhancers.5kbup.5kbdown.pwm80.cor04()))
   all.models <- c(models.01, models.02, models.03, models.04, models.05, models.06, models.07, models.08, models.09)

   model.small <- calculate(sgm, list(regions.motifMatching.snpEnhancerOnly=regions.motifMatching.snpEnhancerOnly()))

   tbl.summary <- summarizeModels(sgm, orderBy="rfScore", maxTFpredictors=10, all.models)

   for(i in seq_len(length(all.models))){
      tbl <- all.models[[i]]$model
      printf("%s: %d %d", names(all.models)[i], nrow(tbl), ncol(tbl))
      tbl2 <- roundNumericColumns(tbl, 3, "lassoPValue")
      all.models[[i]]$model <- tbl2
      print(head(tbl2, n=3))
      }

   save(all.models, tbl.summary, file="../shinyApp/trem2-models.RData")

} # run
#------------------------------------------------------------------------------------------------------------------------
