library(trenaSGM)
genomeName <- "hg38"
targetGene <- "COL1A1"
#f <- "trena.model.inpp5d.164kb.buildSpec.RData"
#f <- "trena.model.test01.buildSpec.Rdata"
f <- "trena.model.encode.buildSpec.RData"
print(load(f))
build.spec$regions <- build.spec$regions[, 1:3]
colnames(build.spec$regions) <- c("chrom", "start", "end")
fpBuilder <- FootprintDatabaseModelBuilder(genomeName, targetGene, build.spec, quiet=FALSE)
x <- build(fpBuilder)
