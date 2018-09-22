library(RSQLite)
#------------------------------------------------------------------------------------------------------------------------
source("constants.R")
#------------------------------------------------------------------------------------------------------------------------
print(load("~/github/trenaSGM/inst/batch/cory/AD/tbl.models.all.RData"))
head(tbl.models)
target.gene <- "MPZL1"
tbl.model <- subset(tbl.models, targetGene.hgnc==target.gene)
tbl.model <- tbl.model[order(tbl.model$rfScore, decreasing=TRUE),]
top.tfs <- tbl.model$tf.hgnc[1:10]
top.tfs


db <- dbConnect(dbDriver("SQLite"), "~/s/data/public/human/remap-2018/remap-all.sqlite")

dbGetQuery(db, "select * from chipseq limit(3)")

query <- sprintf("select * from chipseq where start >= %d and end <= %d and chr=='%s'", snp.min.loc, snp.max.loc, chromosome)
tbl.cs <- dbGetQuery(db, query)
tbl.cs <- subset(tbl.cs, tf %in% top.tfs)
save(tbl.cs, file="tbl.cs.mpzl1.RData")
cs.tfs <- unique(tbl.cs$tf)
match(cs.tfs, tbl.model$tf.hgnc)








