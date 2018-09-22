library(motifbreakR)
library(SNPlocs.Hsapiens.dbSNP150.GRCh38)
library(BSgenome.Hsapiens.UCSC.hg38)

print(load("tbl.RData"))
print(load("tbl.enhancers.umo463f.RData"))
enhancer.min <- min(tbl.enhancers.umo463f$enhancer_start)
enhancer.max <- max(tbl.enhancers.umo463f$enhancer_end)
enhancer.max - enhancer.min   # 7M

tbl.snps <- tbl
head(tbl.snps)
tbl.snp <- subset(tbl.snps, CHR.x == "1" & bp.hg38 >= enhancer.min & bp.hg38 <= enhancer.max)
dim(tbl.snp)  # 314 21

tbl.mb.snp <- data.frame(chr=rep("chr1", nrow(tbl.snp)),
                         start=tbl.snp$bp.hg38-1,
                         end=tbl.snp$bp.hg38,
                         name=sprintf("chr1:%d:%s:%s", tbl.snp$bp.hg38, tbl.snp$REF, tbl.snp$ALT.y),
                         score=rep(0, nrow(tbl.snp)),
                         strand=rep("+", nrow(tbl.snp)))

dim(tbl.mb.snp)
write.table(tbl.mb.snp, file="snps.bed", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)

snps.mb.frombed <- snps.from.file(file = "snps.bed",
                                  dbSNP = SNPlocs.Hsapiens.dbSNP150.GRCh38,
                                  search.genome = BSgenome.Hsapiens.UCSC.hg38,
                                  format = "bed")
snps.mb.frombed

mdb.human.jaspar.2018 <- query(MotifDb, c("sapiens", "jaspar", "2018"))

results <- motifbreakR(snpList = snps.mb.frombed, filterp = TRUE,
                       pwmList = mdb.human.jaspar.2018,
                       threshold = 1e-4,
                       method = "ic",
                       bkg = c(A=0.25, C=0.25, G=0.25, T=0.25),
                       BPPARAM = BiocParallel::bpparam())
results[results$effect=="strong"]
plotMB(results, rsid="chr1:163818172:T:C")
plotMB(results, rsid="chr1:164143663:T:C")
