build:
	docker build -t pshannon/trenashinybase ./

bash:
	docker run -ti --rm pshannon/trenashinybase bash

host.install:
	R -f installTrenaPrequisites.R

rstudio.bioc:
	docker run -e PASSWORD=password -p 8787:8787 bioconductor/release_core2

rstudio.trena:
	docker run -e PASSWORD=password -p 8787:8787 pshannon/trenashinybase

