# get all the genes on chromsome 2
# build a multi-gene table with all of the genehancer data
source("getAllEnhancers.R")
library(org.Hs.eg.db)
tbl.chrom <- as.data.frame(org.Hs.egCHR)
genes.chrom1 <- subset(tbl.chrom, chromosome=="1")$gene_id   # 4245
keytypes(org.Hs.eg.db)
#  [1] "ACCNUM"       "ALIAS"        "ENSEMBL"      "ENSEMBLPROT"  "ENSEMBLTRANS" "ENTREZID"
#  [7] "ENZYME"       "EVIDENCE"     "EVIDENCEALL"  "GENENAME"     "GO"           "GOALL"
# [13] "IPI"          "MAP"          "OMIM"         "ONTOLOGY"     "ONTOLOGYALL"  "PATH"
# [19] "PFAM"         "PMID"         "PROSITE"      "REFSEQ"       "SYMBOL"       "UCSCKG"
# [25] "UNIGENE"      "UNIPROT"
geneSyms.chrom1 <- select(org.Hs.eg.db, keys=genes.chrom1, keytype="ENTREZID", columns="SYMBOL")$SYMBOL  # 4245
x <- lapply(geneSyms.chrom1, getAllEnhancers)
length(x)
tbl.enhancers.chr1 <- do.call(rbind, x)
dim(tbl.enhancers.chr1)   #71916
save(tbl.enhancers.chr1, file="tbl.enhancers.chr1.RData")

min.loc <- snp.min.loc - 10000
max.loc <- snp.max.loc + 10000
tbl.enhancers.umo463f <- subset(tbl.enhancers.chr1, enhancer_start >= min.loc & enhancer_end <= max.loc)
dim(tbl.enhancers.umo463f)   # 1545 6
save(tbl.enhancers.umo463f, file="tbl.enhancers.umo463f.RData")






