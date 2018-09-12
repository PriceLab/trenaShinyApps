library(RSQLite)
targetGene <- "INPP5D"
chromosome <- "chr2"


if(!exist("tbl.model")){
   load("~/github/trena/inst/extdata/tbl.models.all.RData")
   tbl.model <- subset(tbl.models, targetGene.hgnc == "INPP5D" & rfScore >= 5)
   tfs <- tbl.model$tf.hgnc
   }

if(!exists("tbl.enhancers")){
   source("~/s/data/public/GeneHancer/v4.7/getAllEnhancers.R")
   tbl.enhancers <- getAllEnhancers(targetGene)
   colnames(tbl.enhancers) <- c("chr", "start", "end", "type", "score")
   dim(tbl.enhancers)
   gr.enhancers <- GRanges(tbl.enhancers)
   }

extra.loc <- 5000
min.loc <- min(tbl.enhancers$start) - extra.loc
max.loc <- max(tbl.enhancers$end) + extra.loc
max.loc - min.loc

db.lite <- dbConnect(dbDriver("SQLite"), "~/s/data/public/human/remap-2018/remap-all.sqlite")

# learn the columns
# dbGetQuery(db.lite, "select * from chipseq limit 5")
#    chr   start     end   tf                  name peakStart peakEnd
# 1 chr1 1002062 1002357 ADNP ENCSR440VKE.ADNP.K562   1002179 1002180
# 2 chr1 1124275 1124570 ADNP ENCSR440VKE.ADNP.K562   1124417 1124418

query <- sprintf("select * from chipseq where start >= %d and end <= %d and chr=='%s'", min.loc, max.loc, chromosome)
tbl.cs <- dbGetQuery(db.lite, query)
dim(tbl.cs)

#gr.cs <- GRanges(tbl.cs)
#tbl.ov <- as.data.frame(findOverlaps(gr.cs, gr.enhancers))
#colnames(tbl.ov) <- c("cs", "enhancers")
#dim(tbl.ov)
tfs.with.cs <- intersect(tfs, tbl.cs$tf)   # [1] "SPI1" "TAL1" "FLI1" "LYL1"
tbl.chipSeqInModel <- subset(tbl.cs, tf %in% tfs)
dim(tbl.chipSeqInModel) # 1308 7
save(tbl.csInModel, file="../shinyApp/data/tbl.chipSeqInModel.RData")
