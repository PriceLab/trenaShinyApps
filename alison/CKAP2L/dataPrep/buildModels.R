library(trenaSGM)
library(igvR)
library(FimoClient)
library(GenomicRanges)
library(RUnit)
library (RColorBrewer)
#------------------------------------------------------------------------------------------------------------------------
source("~/github/trenaShinyApps/utils/getEnhancers.R")
#------------------------------------------------------------------------------------------------------------------------
genome <- "hg38"
targetGene <- "CKAP2L"
#sgm <- trenaSGM(genome, targetGene, quiet=FALSE)
chromosome <- "chr2"
tss <- 112764686
snp <- "rs7594852"
snp.bp <- 112764177
snp.loc <- sprintf("%s:%d-%d", chromosome, snp.bp, snp.bp)
#------------------------------------------------------------------------------------------------------------------------
FIMO_HOST <- "localhost"
FIMO_PORT <- 5558
if(!exists("fimo")){
   fimo <- FimoClient(FIMO_HOST, FIMO_PORT, quiet=FALSE)
   }
#------------------------------------------------------------------------------------------------------------------------
load("../data/CPMNorm_5252018.RData")
mtx.cpm <- as.matrix(GeneData_CPM)
load("../data/TMMNorm_572018.RData")
mtx.tmm <- as.matrix(GeneData_TMM)
save(mtx.cpm, mtx.tmm, file="../data/matrices.RData")
load("../data/tbl.pheno.RData")
stopifnot(all(colnames(mtx.cpm) %in% tbl.pheno$SequencingID))
stopifnot(all(colnames(mtx.tmm) %in% tbl.pheno$SequencingID))
#------------------------------------------------------------------------------------------------------------------------
if(!exists(tbl.enhancers)){
   tbl.enhancers <- getEnhancerRegions(targetGene)
   save(tbl.enhancers, file="../data/enhancers.RData")
   }
#------------------------------------------------------------------------------------------------------------------------
initial.tf.roundup <- function()
{
   tbl.promoter <- data.frame(chrom=chromosome, start=tss-5000, end=tss+5000, stringsAsFactors=FALSE)
   tbl.enhancers <- getEnhancerRegions(targetGene)
   tbl.regions <- rbind(tbl.promoter, tbl.enhancers)

   build.spec <- list(title="ckap2l.cpm.5kup.5kdown",
                      type="regions.motifMatching",
                      tss=tss,
                      regions=tbl.regions,
                      matrix=mtx.cpm,
                      pfms=query(MotifDb, "sapiens", c("jaspar2018", "swissregulon", "hocomoco")),
                      matchThreshold=75,
                      motifDiscovery="matchPWM",
                      tfMapping=c("MotifDB", "TFClass"),
                      tfPrefilterCorrelation=0.4,
                      orderModelByColumn="rfScore",
                      solverNames=c("lasso", "lassopv", "pearson", "randomForest", "ridge", "spearman"))

   builder <- RegionsMotifMatchingModelBuilder(genome, targetGene,  build.spec, quiet=FALSE)
   x.cpm <- build(builder)
      # get candidate tfs from this permissive model
      # so we can create a curated fimo database
   tbl.basic.cpm <- subset(x.cpm$model, gene %in% unique(subset(x.cpm$regulatoryRegions, motifRelativeScore > 0.85)$geneSymbol) & pcaMax > 0.5)
   tfs.cpm <- tbl.basic$gene

   build.spec$matrix <- mtx.tmm
   build.spec$tfPrefilterCorrelation <- 0.2
   builder <- RegionsMotifMatchingModelBuilder(genome, targetGene,  build.spec, quiet=FALSE)
   x.tmm <- build(builder)
   tbl.basic.tmm <- subset(x.tmm$model,
                           gene %in% unique(subset(x.tmm$regulatoryRegions, motifRelativeScore > 0.85)$geneSymbol) & pcaMax > 0.5)
   tfs.tmm <- tbl.basic.tmm$gene
   tfs.oi <- unique(c(tfs.cpm, tfs.tmm))
   length(tfs.oi)

   save(x.cpm, x.tmm, tbl.basic, tfs.cpm, tbl.basic.tmm, tfs.tmm, tfs.oi, file="third-tf-roundup.RData")

} # initial.tf.roundup
#------------------------------------------------------------------------------------------------------------------------
build.models <- function()
{
   tbl.regions <- data.frame(chrom=chromosome, start=tss-5000, end=tss+5000, stringsAsFactors=FALSE)

   build.spec <- list(title="CKAP2l.fimo.5kup.5kdown",
                      type="regions.motifMatching",
                      tss=tss,
                      regions=tbl.regions,
                      matrix=mtx.cpm,
                      pfms=list(),           # already in fimo server, about 55 curated pfms
                      matchThreshold=1e-4,
                      motifDiscovery="fimo",
                      tfMapping="MotifDB",
                      tfPrefilterCorrelation=0.2,
                      orderModelByColumn="pcaMax",
                      solverNames=c("lasso", "lassopv", "pearson", "randomForest", "ridge", "spearman"))

   builder.10k.cpm <- RegionsFimoModelBuilder(genome, targetGene,  build.spec, fimo, quiet=TRUE)
   x.10k.cpm <- build(builder.10k.cpm)

   build.spec$matrix <- mtx.tmm
   builder.10k.tmm <- RegionsFimoModelBuilder(genome, targetGene,  build.spec, fimo, quiet=TRUE)
   x.10k.tmm <- build(builder.10k.tmm)

   build.spec$matrix <- mtx.cpm
   build.spec$regions <- tbl.enhancers
   builder.enhancers <- RegionsFimoModelBuilder(genome, targetGene,  build.spec, fimo, quiet=TRUE)
   x.enhancers.cpm <-build(builder.enhancers)

   build.spec$matrix <- mtx.tmm
   build.spec$regions <- tbl.enhancers
   builder.enhancers <- RegionsFimoModelBuilder(genome, targetGene,  build.spec, fimo, quiet=TRUE)
   x.enhancers.tmm <-build(builder.enhancers)


   tbl.x1 <- roundNumericColumns(x.10k.cpm$model, 2, "lassoPValue")
   tbl.x2 <- roundNumericColumns(x.10k.tmm$model, 2, "lassoPValue")
   tbl.x3 <- roundNumericColumns(x.enhancers.cpm$model, 2, "lassoPValue")
   tbl.x4 <- roundNumericColumns(x.enhancers.cpm$model, 2, "lassoPValue")

   x.10k.cpm$model <- tbl.x1
   x.10k.tmm$model <- tbl.x2
   x.enhancers.cpm$model <- tbl.x3
   x.enhancers.tmm$model <- tbl.x4

   models <- list(promoter.10k.cpm=x.10k.cpm,
                  promoter.10k.tmm=x.10k.tmm,
                  enhancers.cpm=x.enhancers.cpm,
                  enhancers.tmm=x.enhancers.tmm)
   save(models, file="../data/models.RData")


} # build.models
#------------------------------------------------------------------------------------------------------------------------
plotMotifs <- function(motifs)
{
   stopifnot(length(motifs) > 0)

  if(length(motifs) == 1){
    pcm <- new("pcm", mat=motifs[[1]], name=names(motifs))
    plot(pcm)
    return()
    }

  motifStack(lapply(names(motifs), function(mName) new("pfm", motifs[[mName]], name=mName)))

} # plotMotifs
#------------------------------------------------------------------------------------------------------------------------
identify.motifs <- function()
{
   mdb.base <- query(MotifDb, "sapiens", c("jaspar2018", "hocomoco", "swissregulon"))  # 1846
   mdb.tf <- query(mdb.base, "", tfs.oi)   # 140
   library(motifStack)
   pwm.oi <- list()
   for(tf in tfs.oi[35:57]){
      motifs <- query(mdb.tf, tf)
      if(length(motifs) == 0)
         next
      printf('pwm.oi["%s"] <- "%s"', tf, names(motifs))
      plotMotifs(motifs)
      browser()
      }
   pwm.oi["SP2"] <- "Hsapiens-jaspar2018-SP2-MA0516.1"
   pwm.oi["NFE2L2"] <- "Hsapiens-jaspar2018-NFE2L2-MA0150.1"
   pwm.oi["ZFHX3"] <- "Hsapiens-HOCOMOCOv10-ZFHX3_HUMAN.H10MO.D"
   pwm.oi["HBP1"] <- "Hsapiens-HOCOMOCOv10-HBP1_HUMAN.H10MO.D"
   pwm.oi["ETS1"] <- "Hsapiens-jaspar2018-ETS1-MA0098.3"
   pwm.oi["ETV6"] <- "Hsapiens-jaspar2018-ETV6-MA0645.1"
   pwm.oi["TFCP2L1"] <- "Hsapiens-SwissRegulon-TFCP2L1.SwissRegulon"
   pwm.oi["SOX12"] <- "Hsapiens-SwissRegulon-SOX12.SwissRegulon"
   pwm.oi["ETV5"] <- "Hsapiens-jaspar2018-ETV5-MA0765.1"
   pwm.oi["MEF2D"] <- "Hsapiens-jaspar2018-MEF2D-MA0773.1"
   pwm.oi["RFX7"] <- "Hsapiens-SwissRegulon-RFX7.SwissRegulon"
   pwm.oi["FOXM1"] <- "Hsapiens-HOCOMOCOv10-FOXM1_HUMAN.H10MO.D"
   pwm.oi["SREBF1"] <- "Hsapiens-jaspar2018-SREBF1-MA0595.1"
   pwm.oi["MYBL2"] <- "Hsapiens-jaspar2018-MYBL2-MA0777.1"
   pwm.oi["E2F1"] <- "Hsapiens-jaspar2018-E2F1-MA0024.3"
   pwm.oi["FOXO3"] <- "Hsapiens-jaspar2018-FOXO3-MA0157.2"
   pwm.oi["POU3F2"] <- "Hsapiens-jaspar2018-POU3F2-MA0787.1"
   pwm.oi["CEBPG"] <- "Hsapiens-jaspar2018-CEBPG-MA0838.1"
   pwm.oi["ZNF263"] <- "Hsapiens-jaspar2018-ZNF263-MA0528.1"
   pwm.oi["POU3F3"] <- "Hsapiens-jaspar2018-POU3F3-MA0788.1"
   pwm.oi["IRX2"] <- "Hsapiens-HOCOMOCOv10-IRX2_HUMAN.H10MO.D"
   pwm.oi["NR2F6"] <- "Hsapiens-HOCOMOCOv10-NR2F6_HUMAN.H10MO.D"
   pwm.oi["SMAD4"] <- "Hsapiens-HOCOMOCOv10-SMAD4_HUMAN.H10MO.C"
   pwm.oi["IRF6"] <- "Hsapiens-SwissRegulon-IRF6.SwissRegulon"
   pwm.oi["PITX2"] <- "Hsapiens-HOCOMOCOv10-PITX2_HUMAN.H10MO.D"
   pwm.oi["TGIF1"] <- "Hsapiens-jaspar2018-TGIF1-MA0796.1"
   pwm.oi["CEBPZ"] <- "Hsapiens-HOCOMOCOv10-CEBPZ_HUMAN.H10MO.D"
   pwm.oi["ATF4"] <- "Hsapiens-jaspar2018-ATF4-MA0833.1"
   pwm.oi["CEBPB"] <- "Hsapiens-jaspar2018-CEBPB-MA0466.2"
   pwm.oi["NFIL3"] <- "Hsapiens-jaspar2018-NFIL3-MA0025.1"
   pwm.oi["MXI1"] <- "Hsapiens-jaspar2018-MXI1-MA1108.1"
   pwm.oi["CEBPA"] <- "Hsapiens-jaspar2018-CEBPA-MA0102.3"
   pwm.oi["STAT4"] <- "Hsapiens-HOCOMOCOv10-STAT4_HUMAN.H10MO.D"
   pwm.oi["JUND"] <- "Hsapiens-jaspar2018-JUND-MA0491.1"
   pwm.oi["HSF5"] <- "Hsapiens-jaspar2016-HSF4-MA0771.1"
   pwm.oi["ELF3"] <- "Hsapiens-jaspar2018-ELF3-MA0640.1"
   pwm.oi["E4F1"] <- "Hsapiens-HOCOMOCOv10-E4F1_HUMAN.H10MO.D"
   pwm.oi["CBFB"] <- "Hsapiens-SwissRegulon-CBFB.SwissRegulon"
   pwm.oi["NR6A1"] <- "Hsapiens-HOCOMOCOv10-NR6A1_HUMAN.H10MO.B"
   pwm.oi["NR6A1"] <- "Hsapiens-HOCOMOCOv10-NR6A1_HUMAN.H10MO.B"
   pwm.oi["GRHL1"] <- "Hsapiens-HOCOMOCOv10-GRHL1_HUMAN.H10MO.D"
   pwm.oi["DLX6"] <- "Hsapiens-jaspar2018-DLX6-MA0882.1"
   pwm.oi["SP3"] <- "Hsapiens-jaspar2018-SP3-MA0746.1"
   pwm.oi["BACH1"] <- "Hsapiens-HOCOMOCOv10-BACH1_HUMAN.H10MO.A"
   pwm.oi["CDC5L"] <- "Hsapiens-HOCOMOCOv10-CDC5L_HUMAN.H10MO.D"
   pwm.oi["GABPA"] <- "Hsapiens-jaspar2018-GABPA-MA0062.1"
   pwm.oi["ELF2"] <- "Hsapiens-HOCOMOCOv10-ELF2_HUMAN.H10MO.C"
   pwm.oi["TFAP2C"] <- "Hsapiens-jaspar2018-TFAP2C-MA0524.2"
   pwm.oi["GCM1"] <- "Hsapiens-jaspar2018-GCM1-MA0646.1"
   pwm.oi["PRDM4"] <- "Hsapiens-HOCOMOCOv10-PRDM4_HUMAN.H10MO.D"
   pwm.oi["ZNF410"] <- "Hsapiens-jaspar2018-ZNF410-MA0752.1"
   pwm.oi["RXRB"] <- "Hsapiens-jaspar2018-RXRB-MA0855.1"
   pwm.oi["THAP1"] <- "Hsapiens-jaspar2018-THAP1-MA0597.1"
   pwm.oi["MAFK"] <- "Hsapiens-jaspar2018-MAFK-MA0496.2"
   pwm.oi["FOXI3"] <- "Hsapiens-jaspar2018-FOXD1-MA0031.1"
   pwm.oi["THAP4"] <- "Hsapiens-jaspar2018-THAP1-MA0597.1"
   pwm.oi["GLIS2"] <- "Hsapiens-jaspar2018-GLIS2-MA0736.1"
   pwm.oi["THAP8"] <- "Hsapiens-jaspar2018-THAP1-MA0597.1"

   x <- MotifDb[unlist(pwm.oi, use.names=FALSE)]
   export(x, "~/github/fimoService/pfms/CKAP2L-57tfs.meme", format="meme")

} # identify.motif
#------------------------------------------------------------------------------------------------------------------------




