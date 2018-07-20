library(trenaSGM)
library(igvR)
library(GenomicRanges)
library(RUnit)
library (RColorBrewer)
#------------------------------------------------------------------------------------------------------------------------
source("~/github/trenaShinyApps/utils/getEnhancers.R")
print(load("~/s/data/public/GeneHancer/tbl.enhancers.RData"))
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
print(load("CPMNorm_5252018.RData"))   # GeneData_CPM
print(load("TMMNorm_572018.RData"))    # GeneData_TMM  Data_TMM
print(load("~/github/projects/priceLab/alison/trenaDB/placenta.draft.RData"));
tbl.models <- placenta.rank   # 116256     12
tbl.ckap2l <- subset(tbl.models, target.gene == targetGene)  # [1] 10 12
tfs <- tbl.ckap2l$gene  # some very weak

#------------------------------------------------------------------------------------------------------------------------
if(!exists("igv")){
   igv <- igvR()
   Sys.sleep(5)
   setGenome(igv, genome)
   Sys.sleep(5)
   showGenomicRegion(igv, targetGene)
   track <- DataFrameAnnotationTrack(snp, data.frame(chrom=chromosome, start=snp.bp, end=snp.bp), color="red", trackHeight=30)
   displayTrack(igv, track)
   #tbl.enhancers <- getEnhancerRegions("IL1B") # targetGene)
   #track <- DataFrameAnnotationTrack("IL1B enhancers", tbl.enhancers, color="black", trackHeight=30)
   #displayTrack(igv, track)
   }
#------------------------------------------------------------------------------------------------------------------------
tbl.snpHits <- subset(tbl.enhancers, chrom==chromosome & start <= snp.bp & end >= snp.bp)
goi <- c(targetGene, tbl.snpHits$geneSymbol)
for(gene in goi){
   tbl.enhancers.1gene <- getEnhancerRegions(gene) # targetGene)
   title <- sprintf("%s enhancers", gene)
   track <- DataFrameAnnotationTrack(title, tbl.enhancers.1gene, color="black", trackHeight=30)
   displayTrack(igv, track)
   #Sys.sleep(3)
   }



findTfBindingSites <- function()
{
   tfs <- subset(tbl.models, target.gene == targetGene)$gene
   pfms.oi <- query(MotifDb, "sapiens", tfs)
   export(pfms.oi, "~/github/fimoService/pfms/ckap2l.meme", 'meme')

   mcols(pfms.oi)$geneSymbol
   tbl.motifs <- as.data.frame(table(mcols(pfms.oi)$geneSymbol))

   mm <- MotifMatcher(genome, as.list(pfms.oi), quiet=TRUE)
   tbl.region <- with(getGenomicRegion(igv), data.frame(chrom=chrom, start=start, end=end))
   tbl.seq <- getSequence(mm, tbl.region)
   sequence <- tbl.seq$seq  # 318k
   stopifnot(nchar(sequence) == (1 + tbl.region$end - tbl.region$start))
   tbl.fimo.raw <- requestMatch(fimo, list(try1=sequence))
   tbl.fimo <- tbl.fimo.raw
   tbl.fimo$start <- tbl.fimo$start + tbl.region$start
   tbl.fimo$stop <-  tbl.fimo$stop + tbl.region$start
   dim(tbl.fimo)
   tbl.fimo <- subset(tbl.fimo, p.value <= 0.001)
   dim(tbl.fimo)
   tbl.fimo$chrom <- chromosome

   tbl.foi <- tbl.fimo[grep("FOXM1", tbl.fimo$motif),]
   tbl.hoi <- subset(tbl.foi, -log10(p.value) > 6)

   track <- DataFrameAnnotationTrack("FOXM1", tbl.hoi[, c("chrom", "start", "stop")], color="blue", trackHeight=200)
   displayTrack(igv, track)

} # findTfBindingSites
#------------------------------------------------------------------------------------------------------------------------
findAllBindingSites <- function()
{

   pfms.oi <- query(MotifDb, "sapiens", c("jaspar2018", "hocomoco", "swissregulon"))
   export(pfms.oi, "~/github/fimoService/pfms/hsapiens1846.meme", 'meme')

   mcols(pfms.oi)$geneSymbol
   mm <- MotifMatcher(genome, as.list(pfms.oi), quiet=TRUE)
   tbl.region <- with(getGenomicRegion(igv), data.frame(chrom=chrom, start=start, end=end))
   tbl.seq <- getSequence(mm, tbl.region)
   sequence <- tbl.seq$seq  # 318k
   stopifnot(nchar(sequence) == (1 + tbl.region$end - tbl.region$start))
   tbl.fimo.raw <- requestMatch(fimo, list(try1=sequence))
   tbl.fimo <- tbl.fimo.raw
   tbl.fimo$start <- tbl.fimo$start + tbl.region$start
   tbl.fimo$stop <-  tbl.fimo$stop + tbl.region$start
   dim(tbl.fimo)
   tbl.fimo <- subset(tbl.fimo, p.value <= 0.001)
   dim(tbl.fimo)
   tbl.fimo$chrom <- chromosome
   dim(tbl.fimo)
   hist(-log10(tbl.fimo$q.value))
   tbl.foi <- subset(tbl.fimo, -log10(q.value) > 2)
   dim(tbl.foi)
   tbl.foi <- tbl.foi[order(tbl.foi$p.value, decreasing=FALSE),]
   track <- DataFrameAnnotationTrack("1846", tbl.foi[, c("chrom", "start", "stop", "motif")], color="blue", trackHeight=200)
   displayTrack(igv, track)




} # findAllBindingSites
#------------------------------------------------------------------------------------------------------------------------
nearbyGenes <- function()
{
   library(RPostgreSQL)
   database.host <- "bddsrds.globusgenomics.org"
   db.hg38 <- dbConnect(PostgreSQL(), user= "trena", password="trena", dbname="hg38", host=database.host)
     # dbListTables(db.hg38)   # [1] "motifsgenes" "gtf"

   base.query <- "select * from gtf where moleculetype='gene' and gene_biotype='protein_coding' and"
   loc.query <- sprintf("chr='chr2' and start > (%d-500000) and endpos < (%d+500000)", snp.bp, snp.bp)
   query <- sprintf("%s %s", base.query, loc.query)
   tbl <- dbGetQuery(db.hg38, query) # [, c(1:3,5,10)]
   geneSymbols <- tbl$gene_name
   save(geneSymbols, file="geneSymbols.1MB.RData")


   for(gene in geneSymbols){
      printf(gene)
      tbl.enhancers <- getEnhancerRegions(gene) # targetGene)
      track <- DataFrameAnnotationTrack(sprintf("%s enhancers", gene), tbl.enhancers, color="black", trackHeight=30)
      displayTrack(igv, track)
      }

   # build noDNA models for these 18 genes

} # nearbyGenes
#------------------------------------------------------------------------------------------------------------------------
newExpressionMatrix <- function()
{

  print(load("CPMNorm_5252018.RData"))   # GeneData_CPM
  print(load("TMMNorm_572018.RData"))    # GeneData_TMM  Data_TMM

  mtx.cpm <- as.matrix(GeneData_CPM)
  mtx.tmm <- as.matrix(GeneData_TMM)
  dim(mtx.cpm); dim(mtx.tmm)

  tfs <- allKnownTFs()
  length(tfs)
  print(load("geneSymbols.1MB.RData"))
  target.genes <- geneSymbols
  goi <- c(tfs, "CKAP2L")
  length(goi)
  goi <- intersect(rownames(mtx.cpm), goi)
  length(goi)
  mtx.oi <- t(mtx.cpm[goi,])
  mtx.cor <- cor(mtx.oi)
  dim(mtx.cor)
  which(mtx.cor["CKAP2L",] > 0.75)
   #  AFF1    AFF4  CKAP2L   FOXM1  NFE2L2   PRDM5   TERF1   ZC3H8  ZNF37A  ZNF480 ZNF518A  ZNF566
   #     6       9      77     169     399     488     657     759     918     958     970     998


} # newExpressionMatrix
#------------------------------------------------------------------------------------------------------------------------
buildMotifMatchingModel <- function()
{
   genome <- "hg38"
   targetGene <- "CKAP2L"
   chromosome <- "chr2"
      # strand-aware start and end: ckap2l is on the minus strand
   tbl.regions <- data.frame(chrom=chromosome, start=tss-200, end=tss+2000, stringsAsFactors=FALSE)

   build.spec <- list(title="ckap2l.rmm.2000up.200down",
                      type="regions.motifMatching",
                      tss=tss,
                      regions=tbl.regions,
                      matrix=mtx.cpm,
                      pfms=query(MotifDb, "sapiens", "jaspar2018"),
                      matchThreshold=90,
                      motifDiscovery="matchPWM",
                      tfMapping="MotifDB",
                      tfPrefilterCorrelation=0.4,
                      orderModelByColumn="rfScore",
                      solverNames=c("lasso", "lassopv", "pearson", "randomForest", "ridge", "spearman"))

   builder <- RegionsMotifMatchingModelBuilder(genome, targetGene,  build.spec, quiet=TRUE)
   x <- build(builder)

   tbl.regions2 <- data.frame(chrom=chromosome, start=tss-5000, end=tss+5000, stringsAsFactors=FALSE)
   build.spec2 <- build.spec
   build.spec2$title <- "ckap2l.rmm.5kup.5kdown"
   build.spec2$regions <- tbl.regions2
   build.spec2$pfms <- query(MotifDb, "sapiens", c("jaspar2018", "hocomoco", "swissregulon"))
   build.spec2$tfMapping=c("MotifDb", "TFClass")

   build.spec2$matchThreshold=75


   builder.2 <- RegionsMotifMatchingModelBuilder(genome, targetGene,  build.spec2, quiet=TRUE)
   x2 <- build(builder.2)
   save(x2, file="x2.regionsMotifMatchingModel.RData")

   color.palette = 'Set2'
   colors <- brewer.pal(length(model.tfs), color.palette)
   model.tfs <- x2$model$gene[1:6]

   for(tf in model.tfs){
      tbl.tf <- subset(x2$regulatoryRegions, geneSymbol==tf)[, c("chrom", "motifStart", "motifEnd")]
      tbl.tf <- unique(tbl.tf[order(tbl.tf$motifStart, decreasing=FALSE),])
      track <- DataFrameAnnotationTrack(tf, tbl.tf, color=colors[grep(tf, model.tfs)])
      displayTrack(igv, track)
      }

} # buildMotifMatchingModel
#------------------------------------------------------------------------------------------------------------------------
buildNoDnaModel <- function()
{
   build.spec <- list(title="ckap2l.noDNA.allTFs",
                      type="noDNA.tfsSupplied",
                      matrix=mtx.cpm,
                      tfs=allKnownTFs(),
                      tfPrefilterCorrelation=0.4,
                      orderModelByColumn="rfScore",
                      solverNames=c("lasso", "lassopv", "pearson", "randomForest", "ridge", "spearman"))

   builder <- NoDnaModelBuilder(genome, targetGene,  build.spec, quiet=TRUE)
   x <- build(builder)

} # buildNoDnaModel
#------------------------------------------------------------------------------------------------------------------------
