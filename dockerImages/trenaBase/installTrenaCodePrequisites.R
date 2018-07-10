printf <- function(...) print(noquote(sprintf(...)))
my.user.library <- "~/library"
if(!my.user.library %in% .libPaths())
    .libPaths(my.user.library)

source("http://bioconductor.org/biocLite.R")


code.pkgs <- c("GenomicRanges",
               "Biostrings",
               "jsonlite",
               "httpuv",
               "BiocGenerics",
               "BiocParallel",
               "DBI",
               "RPostgreSQL",
               "RMySQL",
               "RSQLite",
               "glmnet",
               "lassopv",
               "randomForest",
               "flare",
               "vbsr",
               "stringr",
               "httpuv",
               "colorspace",
               "annotate",
               "shiny",
               "DT",
               "htmlwidgets",
               "later",
               "splitstackshape",
               "RUnit")

for(code.pkg in code.pkgs){
   suppressWarnings(
      needed <- !require(code.pkg, character.only=TRUE, lib.loc=my.user.library, quiet=TRUE)
      )
   printf("%s needed? %s", code.pkg, needed)
   if(needed)
      biocLite(code.pkg, quiet=FALSE, lib=my.user.library)
   } # for




