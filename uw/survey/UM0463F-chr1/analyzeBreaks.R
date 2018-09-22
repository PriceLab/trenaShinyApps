library(GenomicRanges)
print(load("results.all.snps.RData"))
tbl.results <- as.data.frame(mcols(results))
dim(tbl.results)
head(tbl.results)

print(load("tbl.RData"))
tbl.snps <- tbl
print(load("tbl.enhancers.umo463f.RData"))
enhancer.min <- min(tbl.enhancers.umo463f$enhancer_start)
enhancer.max <- max(tbl.enhancers.umo463f$enhancer_end)
enhancer.max - enhancer.min   # 7M

tbl.snp <- subset(tbl.snps, CHR.x == "1" & bp.hg38 >= enhancer.min & bp.hg38 <= enhancer.max)

results.strong <- results[results$effect=="strong"]
tbl.strong <- as.data.frame(mcols(results.strong))
tfs.altered.binding <- unique(tbl.strong$geneSymbol)

target.genes.with.enhancers.in.region <- sort(unique(tbl.enhancers.umo463f$geneSymbol))
length(target.genes.with.enhancers.in.region)  # 120

print(load("~/github/trenaSGM/inst/batch/cory/AD/tbl.models.all.RData"))
enhancer.genes.with.models <- sort(unique(intersect(target.genes.with.enhancers.in.region, tbl.models$targetGene.hgnc)))
length(enhancer.genes.with.models)  #36

x <- lapply(enhancer.genes.with.models, function(enhancer.gene) subset(tbl.models, targetGene.hgnc==enhancer.gene)[1:10,])
tbl.enhancer.gene.models <- do.call(rbind, x)
dim(tbl.enhancer.gene.models)  # 360 15
head(tbl.enhancer.gene.models)
tbl.enhancer.gene.models.tfs <- sort(unique(tbl.enhancer.gene.models$tf.hgnc))
length(tbl.enhancer.gene.models.tfs) # 214
tbl.enhancer.gene.models.tfs.altered.binding <- intersect(tfs.altered.binding, tbl.enhancer.gene.models.tfs)
length(tbl.enhancer.gene.models.tfs.altered.binding) 47
tbl.oi <- subset(tbl.enhancer.gene.models, tf.hgnc %in% tbl.enhancer.gene.models.tfs.altered.binding)
head(tbl.oi)

    #--------------------------------------------------
    # look at PBX1
    #--------------------------------------------------
subset(tbl.oi, targetGene.hgnc == "PBX1")  # NFKB1, NR2C2, top 2 regulators by pcaMax score

        targetGene.hgnc tf.hgnc targetGene.ensembl      tf.ensembl   betaLasso  lassoPValue pearsonCoeff   rfScore    betaRidge spearmanCoeff betaSqrtLasso concordance   pcaMax bindingSites rank
2464027            PBX1   NFKB1    ENSG00000185630 ENSG00000109320 -0.10961711 3.500604e-26   -0.5903090 20.650174 -0.025994330    -0.5867842   -0.05335293   0.4749515 1.979561           NA    1
2464028            PBX1   NR2C2    ENSG00000185630 ENSG00000177463  0.17338876 9.109431e-16    0.4981506  5.298252  0.042919055     0.5038515    0.19143495   0.5763795 1.741863           NA    2
2464033            PBX1    ATF4    ENSG00000185630 ENSG00000128272 -0.09575517 2.005160e-19   -0.5517552  5.965769 -0.028864270    -0.5392993   -0.06079828   0.4439064 1.136085           NA    7
2464034            PBX1   MEF2C    ENSG00000185630 ENSG00000081189  0.00000000 9.732160e-01    0.4664913  6.131910  0.005928837     0.5208519    0.00000000   0.4024187 1.132957           NA    8

unique(subset(tbl.strong, geneSymbol %in% c("NFKB1", "NR2C2", "ATF4", "MEF2C")))

chr1:167145933:T:C     T   C 167145933        8       ATF4 jaspar2018     MA0833.1   MA0833.1 ggttgTctcataa              0.7861738 0.9257291 10.905422 12.84127        NA        NA 0.0009017133 0.9954914 strong
chr1:165317326:T:G.1   T   G 165317326       10      MEF2C jaspar2018     MA0497.1   MA0497.1 gaataTttttacagg            0.8914680 0.7459791 10.507320  8.79251        NA        NA 0.9665006790 0.0000000 strong
chr1:167549187:A:G.1   A   G 167549187        1      NFKB1 jaspar2018     MA0105.1   MA0105.1      Agggctttcc            0.7979674 0.9658307  8.487376 10.27281        NA        NA 0.0000000000 1.0000000 strong
chr1:165205553:T:C.2   T   C 165205553       12      NR2C2 jaspar2018     MA0504.1   MA0504.1   tgaTctctggcatga          0.7264795 0.8080439  9.334363 10.38236        NA        NA 0.0278481013 0.8506329 strong

# to generalize this:
# obtain tss of each gene
# calculate distance to high-scoring broken motifs
# unlikely, but see if any fall in this gene's enhancer regions
# build a table with distance tss, trena rank, summarize over all genes 9, distance to nearest gene-specific enhancer

    #--------------------------------------------------
    # look at POU2F1
    #--------------------------------------------------
subset(tbl.oi, targetGene.hgnc == "POU2F1")  # NFKB1, NR2C2, top 2 regulators by pcaMax score

        targetGene.hgnc tf.hgnc targetGene.ensembl      tf.ensembl   betaLasso  lassoPValue pearsonCoeff   rfScore    betaRidge spearmanCoeff betaSqrtLasso concordance   pcaMax bindingSites rank
2464027            PBX1   NFKB1    ENSG00000185630 ENSG00000109320 -0.10961711 3.500604e-26   -0.5903090 20.650174 -0.025994330    -0.5867842   -0.05335293   0.4749515 1.979561           NA    1
2464028            PBX1   NR2C2    ENSG00000185630 ENSG00000177463  0.17338876 9.109431e-16    0.4981506  5.298252  0.042919055     0.5038515    0.19143495   0.5763795 1.741863           NA    2
2464033            PBX1    ATF4    ENSG00000185630 ENSG00000128272 -0.09575517 2.005160e-19   -0.5517552  5.965769 -0.028864270    -0.5392993   -0.06079828   0.4439064 1.136085           NA    7
2464034            PBX1   MEF2C    ENSG00000185630 ENSG00000081189  0.00000000 9.732160e-01    0.4664913  6.131910  0.005928837     0.5208519    0.00000000   0.4024187 1.132957           NA    8

unique(subset(tbl.strong, geneSymbol %in% c("NFIL3", "NRL")))


subset(tbl.oi, targetGene.hgnc=="SFT2D2")
unique(subset(tbl.strong, geneSymbol %in% c("REST")))
