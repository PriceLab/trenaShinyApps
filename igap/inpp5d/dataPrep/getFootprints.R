library(trena)

targetGene <- "INPP5D"

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
chrom <- tbl.enhancers$chr[1]

genome.db.uri    <- "postgres://khaleesi/hg38"                       # has gtf and motifsgenes tables
footprint.db.uris <- c("postgres://khaleesi/brain_hint_20",
                       "postgres://khaleesi/brain_hint_16",
                       "postgres://khaleesi/brain_wellington_20",
                       "postgres://khaleesi/brain_wellington_16")


tbl.fp.raw <- data.frame()
for(db in footprint.db.uris){
   fpf <- FootprintFinder(genome.db.uri, db, quiet=FALSE)
   tbl.fp.new <- getFootprintsInRegion(fpf, chrom, min.loc, max.loc)
   tbl.fp.raw <- rbind(tbl.fp.raw, tbl.fp.new)
   }


tbl.fp <- tbl.fp.raw

tbl.fp$tf <- ""
tbl.fp$tf <- as.character(tbl.fp$tf)
dim(tbl.fp)



spi1.indices  <- grep("SPI1", tbl.fp$name)
length(spi1.indices)
fli1.indices <-  grep("FLI1", tbl.fp$name)
length(fli1.indices)
lyl1.indices <-  grep("MESP1", tbl.fp$name)
length(lyl1.indices)
tal1.indices <-  grep("TAL1", tbl.fp$name)
length(tal1.indices)

tbl.fp$tf[spi1.indices] <- "SPI1"
tbl.fp$tf[fli1.indices] <- "FLI1"
tbl.fp$tf[tal1.indices] <- "TAL1"
tbl.fp$tf[lyl1.indices] <- "LYL1"


indices <- sort(c(spi1.indices, fli1.indices, tal1.indices, lyl1.indices))
tbl.fp <- tbl.fp[indices, ]
dim(tbl.fp)
table(tbl.fp$tf)

tbl.fp <- tbl.fp[, c("chrom", "fp_start", "fp_end", "tf")]
colnames(tbl.fp) <- c("chr", "start", "end", "tf")
dim(tbl.fp)

save(tbl.fp, file="../shinyApp/data/tbl.fp.RData")

