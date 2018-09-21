library(curl)
library(RUnit)
txt <- readLines(curl("https://igv-data.systemsbiology.net/static/tair10/chromosomeAliases.txt"))
checkEquals(length(txt), 7)
tokens <- strsplit(txt[1], "\t")[[1]]
checkEquals(tokens[1], "1")
checkEquals(tokens[2], "Chr1")
checkEquals(tokens[3], "CHR1")
checkEquals(tokens[4], "chr1")
