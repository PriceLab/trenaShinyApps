#f <- "~/s/data/public/human/chip-seq/DNS.Neu.50.AllAg.AllCell-chr2.bed"
#tbl <- read.table(f, skip=0, nrow=-1, sep="\t", as.is=TRUE)
#colnames(tbl) <- c("chr", "start", "end", "desc", "value", "strand", "start2", "end2", "rgb")
library(rtracklayer)
chrom <- "chr2"
tss <- 233060416
tbl.sub <- subset(tbl, start > tss - 20000 & end < tss + 20000)

f <- "~/s/data/public/human/chip-seq/SRX1023790.bw"
file.exists(f)
gr <- import(f, which=seqinfo(BigWigFile(f))["chr2"])
tbl.chr2 <- as.data.frame(gr)
dim(tbl.chr2)   # 3354818 x 6
colnames(tbl.chr2) <- c("chr", "start", "end", "width", "strand", "score")

#----------------------------------------------------------------------------------------------------
# direct ucsc queries
#----------------------------------------------------------------------------------------------------
library(RMariaDB)
driver <- RMariaDB::MariaDB()
host <- "genome-mysql.cse.ucsc.edu"
user <- "genome"
driver <- RMariaDB::MariaDB()
dbname <- "hg19"
db <- DBI::dbConnect(driver, user = user, host = host, dbname = dbname)
all.tableNames <- DBI::dbListTables(db);  # 11132

# [1] "wgEncodeRegTfbsClustered"         "wgEncodeRegTfbsClusteredInputs"
# [3] "wgEncodeRegTfbsClusteredInputsV3" "wgEncodeRegTfbsClusteredV2"
# [5] "wgEncodeRegTfbsClusteredV3"

load("~/github/trenaShinyApps/igap/inpp5d/shinyApp/data/tbl.enhancers.RData")
enhancer.loc.min <- min(tbl.enhancers$start)
enhancer.loc.max <- max(tbl.enhancers$end)

tbl <-  dbGetQuery(db, "select chrom, chromStart, chromEnd, name, score, expCount  from wgEncodeRegTfbsClusteredV3 where chrom='chr2'")
dim(tbl)  # 349645      6
    # inpp5d's hg19 proximal promoter, 7kb
   # chr2:233,923,356-233,930,491
tbl.inpp5d.prox <- subset(tbl, chromStart >= enhancer.loc.min & chromEnd <= enhancer.loc.max)
dim(tbl.inpp5d.prox)  # 695 x 6
tbl.xtab <- as.data.frame(table(tbl.inpp5d.prox$name))
tbl.xtab <- tbl.xtab[order(tbl.xtab$Freq, decreasing=TRUE),]
dim(tbl.xtab)

load("~/github/trenaShinyApps/igap/inpp5d/shinyApp/data/tbl.model.RData")
intersect(tbl.xtab$Var1, tbl.model$tf.hgnc[1:10])   # SPI1  TAL1
intersect(tbl.xtab$Var1, tbl.model$tf.hgnc[1:20])   # SPI1  TAL1 FOXP2 TBP
intersect(tbl.xtab$Var1, tbl.model$tf.hgnc[1:30])   # SPI1  TAL1 FOXP2 TBP KDM5B
intersect(tbl.xtab$Var1, tbl.model$tf.hgnc[1:40])   # SPI1  EGR1  TBP   MEF2C FOXP2 KDM5B TAL1
intersect(tbl.xtab$Var1, tbl.model$tf.hgnc[1:100])  # SPI1  EGR1  TBP   E2F1  MEF2C E2F6  FOXP2 IKZF1 IRF1  KDM5B MTA3  TAL1

