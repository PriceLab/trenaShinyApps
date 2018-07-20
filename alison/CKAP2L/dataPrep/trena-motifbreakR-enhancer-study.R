snp <- "rs7594852"
chromosome <- "chr2"
tss <- 112764686
snp <- "rs7594852"
snp.bp <- 112764177
snp.loc <- sprintf("%s:%d-%d", chromosome, snp.bp, snp.bp)
print(load("~/github/genomescale_scripts/geneHancer/v4.7/tbl.enhancer.chr2.RData"))

tbl.enhancers.rs7594852 <- subset(tbl.enhancers.chr2, enhancer_start <= snp.bp & enhancer_end >= snp.bp)
target.genes <- tbl.enhancers.rs7594852$geneSymbol
length(target.genes)


library(MotifDb)
library(SNPlocs.Hsapiens.dbSNP150.GRCh38)
library(BSgenome.Hsapiens.UCSC.hg38)
library(motifbreakR)
pfms <- query(MotifDb, "sapiens", c("jsapar2018", "hocomoco", "regulon"))
length(pfms)
export(pfms, "~/github/fimoService/pfms/current.meme", format="meme")

snps.gr <- snps.from.rsid(rsid = snp, dbSNP=SNPlocs.Hsapiens.dbSNP150.GRCh38, search.genome=BSgenome.Hsapiens.UCSC.hg38)

results <- motifbreakR(snpList = snps.gr,
                       filterp = TRUE,
                       pwmList = pfms,
                       show.neutral=TRUE,
                       method = c("ic", "log", "notrans")[2],
                       bkg = c(A=0.25, C=0.25, G=0.25, T=0.25),
                       BPPARAM = BiocParallel::bpparam(),
                       verbose=TRUE)
tbl.results <- as.data.frame(mcols(results))
dim(tbl.results)
tbl.results.filtered <- subset(tbl.results, effect %in% c("weak", "strong") & pctRef > 0.8)
dim(tbl.results.filtered)

print(load("../../trenaDB/placenta.draft.RData"))
tbl.placentalModels <- placenta.rank
candidate.tfs <- unique(subset(tbl.placentalModels, target.gene %in% target.genes & Rank <= 10)$gene)
length(candidate.tfs)  # 44

# print(load("~/github/trenaSGM/inst/batch/cory/AD/tbl.models.all.RData"))
# candidate.tfs <- unique(subset(tbl.models, targetGene.hgnc %in% target.genes & rank <= 10)$tf.hgnc)
# length(candidate.tfs)  # 15
affected.tfs <- intersect(candidate.tfs, tbl.results.filtered$geneSymbol)
length(affected.tfs)
affected.tfs

tbl.hits <- subset(tbl.placentalModels, (gene %in% affected.tfs) & (target.gene %in% target.genes))

# subset(tbl.models, tf.hgnc %in% affected.tfs & targetGene.hgnc %in% target.genes & rank < 10)
# tbl.hits <- subset(tbl.models, tf.hgnc %in% affected.tfs & targetGene.hgnc %in% target.genes & rank < 10)
#         targetGene.hgnc tf.hgnc targetGene.ensembl      tf.ensembl  betaLasso  lassoPValue pearsonCoeff  rfScore   betaRidge spearmanCoeff betaSqrtLasso concordance   pcaMax bindingSites rank
# 964921           POLR1B   INSM1    ENSG00000125630 ENSG00000173404  0.0000000 8.926726e-01    0.2820912 3.298819  0.01363760     0.2909920    0.00000000   0.3762966 2.141165           NA    6
# 2059870           RGPD8    NFIA    ENSG00000169629 ENSG00000162599 -0.1205352 2.068607e-17   -0.5169928 8.425475 -0.02162032    -0.5437562   -0.05705476   0.5428285 1.291639           NA    7


tbl.tf.gene.hits <- merge(tbl.results.filtered, tbl.hits, by.x="geneSymbol", by.y="gene")[c(1,5), ]  # hocomoco over swissregulon

co.all <-  c("geneSymbol", "REF", "ALT", "snpPos", "motifPos", "dataSource", "providerName", "providerId", "seqMatch", "pctRef",
             "pctAlt", "scoreRef", "scoreAlt", "Refpvalue", "Altpvalue",  "alleleRef", "alleleAlt", "effect", "beta.lasso",
             "beta.ridge", "rf.score", "pearson.coeff", "spearman.coeff", "beta.sqrtlasso", "lasso.p.value", "concordance",
             "pcaMax", "target.gene",  "Rank")



coi <-  c("target.gene", "geneSymbol", "REF", "ALT", "snpPos", "motifPos", "dataSource", "providerName", "providerId", "seqMatch",
          "pctRef", "pctAlt",
          "effect", "beta.lasso",  "beta.ridge", "rf.score", "pearson.coeff", "spearman.coeff",  "lasso.p.value",
          "concordance", "pcaMax",  "Rank")


tbl.final <- tbl.tf.gene.hits[,coi]
colnames(tbl.final)[1:2] <- c("targetGene", "TF")

as.data.frame(t(tbl.final))
# targetGene                                            IL1A                             POLR1B
# TF                                                     MAZ                              TGIF1
# REF                                                      C                                  C
# ALT                                                      T                                  T
# snpPos                                           112764177                          112764177
# motifPos                                                 4                                  3
# dataSource                                    SwissRegulon                        HOCOMOCOv10
# providerName                              MAZ.SwissRegulon                TGIF1_HUMAN.H10MO.S
# providerId                                MAZ.SwissRegulon                TGIF1_HUMAN.H10MO.S
# seqMatch       ccccgcctcctggCaca                           ctggCac
# pctRef                                           0.8327501                          0.8121517
# pctAlt                                           0.8022636                          0.6271615
# effect                                                weak                             strong
# beta.lasso                                     -0.07072285                         0.07728170
# beta.ridge                                     -0.09708883                         0.07708446
# rf.score                                          1.387624                           1.126456
# pearson.coeff                                   -0.3132200                          0.2921831
# spearman.coeff                                  -0.3123397                          0.2802222
# lasso.p.value                                 0.0003274400                       0.0003732657
# concordance                                      0.2711878                          0.3026960
# pcaMax                                           0.5729905                          0.6934357
# Rank                                                     7                                  8

