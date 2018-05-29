source("~/github/projects/utils/geneIdMapping/symToGeneID.R"); test_assignGeneIDs()
gad.syms <- c("ABCA7", "APOE", "BIN1", "CASS4", "CD2AP", "CELF1", "CLU", "CR1", "EPHA1", "FERMT2",
              "INPP5D", "MEF2C", "MS4A6A", "NME8", "PICALM", "PTK2B", "SLC24A4", "RIN3",
              "SORL1", "ZCWPW1")
gad.risk <- c("ABCA7", "APOE", "BIN1", "CD2AP", "CELF1", "CR1", "FERMT2", "PTK2B", "INPP5D")
gad.preventive <- c("CASS4", "CLU", "EPHA1", "MEF2C", "MS4A6A", "NME8", "PICALM",
                    "SLC24A4", "RIN3", "SORL1", "ZCWPW1")

gad <- assignGeneIDs(gad.syms)$mapped


library(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
txByGene <- transcriptsBy(txdb, by="gene")
length(txByGene) # [1] 25221
unlisted <- unlist(txByGene)
TSS <- ifelse(strand(unlisted) == "+", start(unlisted), end(unlisted))
TSS <- GRanges(seqnames(unlisted), IRanges(TSS, width=1), strand(unlisted))
TSS_by_gn <- relist(TSS, txByGene)
mcols(TSS) <- mcols(unlisted)
tssByGene <- relist(TSS, txByGene)

library(igvR)
igv <- igvR()
setGenome(igv, "hg38")
showGenomicRegion(igv, "RIN3")
source("~/github/trenaShinyApps/utils/getEnhancers.R")
tbl.enhancers <- getEnhancerRegions("RIN3")
track.enhancers <- DataFrameAnnotationTrack("enhancers", tbl.enhancers, color="black")
displayTrack(igv, track.enhancers)
# build a basic model with this 10119 bp region straddling TSS
#    chrom    start      end
# 42 chr14 92510645 92520763
