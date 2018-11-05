args <- commandArgs(trailingOnly=TRUE)
projectName <- args[1]
library(projectName, character.only=TRUE)
initialization.command <- sprintf("trenaProject <- %s()", projectName)
eval(parse(text=initialization.command))
setTargetGene(trenaProject, getSupportedGenes(trenaProject)[1])
printf("targetGene: %s", getTargetGene(trenaProject))


