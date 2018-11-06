all:  motifDb  trena  igvShiny  cyjShiny trenaSGM

motifDb:
	(cd ~/github/MotifDb; git pull)
	R CMD build --no-build-vignettes ~/github/MotifDb

trena:
	(cd ~/github/trena; git pull)
	R CMD build --no-build-vignettes ~/github/trena

trenaSGM:
	(cd ~/github/trenaSGM; git pull)
	R CMD build --no-build-vignettes ~/github/trenaSGM


igvShiny:
	(cd ~/github/igvShiny; git pull)
	R CMD build --no-build-vignettes ~/github/igvShiny

cyjShiny:
	(cd ~/github/cyjShiny; git pull)
	R CMD build --no-build-vignettes ~/github/cyjShiny

#shiny
#DT
#htmlwidgets
#GenomicRanges
#later
##motifStack
#RUnit

