library(trenaSGM)
load("trem2-models.RData")
for(i in seq_len(length(all.models))){
   tbl <- all.models[[i]]$model
   printf("%s: %d %d", names(all.models)[i], nrow(tbl), ncol(tbl))
   tbl2 <- roundNumericColumns(tbl, 3, "lassoPValue")
   all.models[[i]]$model <- tbl2
   print(head(tbl2, n=3))
}
