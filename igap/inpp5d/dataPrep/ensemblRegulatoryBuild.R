# ensemblRegulatoryBuild.R
#------------------------------------------------------------------------------------------------------------------------
library(biomaRt)
library(igvR)
#------------------------------------------------------------------------------------------------------------------------
targetGene <- "INPP5D"

if(!exists("igv")){
   igv <- igvR()
   setGenome(igv, "hg38")
   showGenomicRegion(igv, targetGene)
   }
#------------------------------------------------------------------------------------------------------------------------
mart.reg <- useMart("ENSEMBL_MART_FUNCGEN")
tbl.datasets <- listDatasets(mart.reg)
dim(tbl.datasets) # 11 3
tbl.datasets[grep("hsapiens", tbl.datasets$dataset),]

#                        dataset                                 description    version
#      hsapiens_external_feature Human Other Regulatory Regions (GRCh38.p12) GRCh38.p12
#  hsapiens_mirna_target_feature     Human miRNA Target Regions (GRCh38.p12) GRCh38.p12
#         hsapiens_motif_feature           Human Binding Motifs (GRCh38.p12) GRCh38.p12
#                  hsapiens_peak      Human Regulatory Evidence (GRCh38.p12) GRCh38.p12
#    hsapiens_regulatory_feature      Human Regulatory Features (GRCh38.p12) GRCh38.p12

datasets <- tbl.datasets[grep("hsapiens", tbl.datasets$dataset), "dataset"]

mart.reg.1 <- useMart("ENSEMBL_MART_FUNCGEN", dataset=datasets[1])
mart.reg.2 <- useMart("ENSEMBL_MART_FUNCGEN", dataset=datasets[2])
mart.reg.3 <- useMart("ENSEMBL_MART_FUNCGEN", dataset=datasets[3])
mart.reg.4 <- useMart("ENSEMBL_MART_FUNCGEN", dataset=datasets[4])
mart.reg.5 <- useMart("ENSEMBL_MART_FUNCGEN", dataset=datasets[5])

listAttributes(mart.reg)
#                        name              description               page
# 1                  activity                 Activity regulatory_feature
# 2      regulatory_stable_id     Regulatory stable ID regulatory_feature
# 3    bound_seq_region_start         Bound start (bp) regulatory_feature
# 4      bound_seq_region_end           Bound end (bp) regulatory_feature
# 5           chromosome_name Chromosome/scaffold name regulatory_feature
# 6          chromosome_start               Start (bp) regulatory_feature
# 7            chromosome_end                 End (bp) regulatory_feature
# 8         feature_type_name             Feature type regulatory_feature
# 9  feature_type_description Feature type description regulatory_feature
# 10           epigenome_name           Epigenome name regulatory_feature
# 11    epigenome_description    Epigenome description regulatory_feature
# 12             so_accession        SO term accession regulatory_feature
# 13                  so_name             SO term name regulatory_feature
# 14                   efo_id       EFO term accession regulatory_feature
listFilters(mart.reg)
#                            name                                           description
# 1               chromosome_name                                       Chromosome Name
# 2                         start                                            Start (bp)
# 3                           end                                              End (bp)
# 4            chromosomal_region                   e.g 1:100:10000000, 1:100000:200000
# 5                    band_start                                            Band Start
# 6                      band_end                                              Band End
# 7                  marker_start                                          Marker Start
# 8                    marker_end                                            Marker End
# 9                 encode_region                                         Encode region
# 10         regulatory_stable_id Regulatory stable ID (e.g. ENSR00000060894) [Max 500]
# 11 regulatory_feature_type_name                                     Feature type name
# 12               epigenome_name                                        Epigenome name
roi <- getGenomicRegion(igv)
coi <- listAttributes(mart.reg)$name #[c(1,5,6,8)]

tbl.1 <- getBM(attributes=listAttributes(mart.reg.1)$name,
               filters=list(chromosome_name="2", start=roi$start, end=roi$end),
               mart=mart.reg.1)
dim(tbl.1)

tbl.2 <- getBM(attributes=listAttributes(mart.reg.2)$name,
               filters=list(chromosome_name="2", start=roi$start, end=roi$end),
               mart=mart.reg.2)
dim(tbl.2)

tbl.3 <- getBM(attributes=listAttributes(mart.reg.3)$name,
               filters=list(chromosome_name="2", start=roi$start, end=roi$end),
               mart=mart.reg.3)
dim(tbl.3)

tbl.4 <- getBM(attributes=listAttributes(mart.reg.4)$name,
               filters=list(chromosome_name="2", start=roi$start, end=roi$end),
               mart=mart.reg.4)
dim(tbl.4)
head(tbl.4)

regions, mart=mart.reg)

ensembl_regulation = useMart(biomart="ENSEMBL_MART_FUNCGEN",host="www.ensembl.org",dataset="hsapiens_motif_feature")
binding_matrix = getBM(attributes=c('binding_matrix_id','chromosome_name', 'chromosome_start', 'chromosome_end','chromosome_strand', 'score'), filters='motif_binding_matrix_id',values="MA0005.2", mart=ensembl_regulation)
