library(trenaSGM)
genomeName <- "hg38"
targetGene <- "COL1A1"
f <- "trena.model.test01.buildSpec.Rdata"
print(load(f))
fpBuilder <- FootprintDatabaseModelBuilder(genomeName, targetGene, build.spec, quiet=FALSE)
x <- build(fpBuilder)
