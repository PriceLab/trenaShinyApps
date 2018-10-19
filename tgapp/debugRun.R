library(trenaSGM)
genomeName <- "hg38"
targetGene <- "COL1A1"
print(load("trena.model.hint20.buildSpec.Rdata"))
fpBuilder <- FootprintDatabaseModelBuilder(genomeName, targetGene, build.spec, quiet=FALSE)
x <- build(fpBuilder)
