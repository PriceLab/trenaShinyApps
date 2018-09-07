runInShell:
	R -f runIreg1InBackground.R

runInBackground:
	nohup R -r runIreg1InBackground.R &


build:
	docker build -t ireg1shinyapp -f Dockerfile .

bash:
	docker run -ti --rm -p 60039:3838 ireg1shinyapp bash

dockerRun:
	docker run -p 60039:3838 ireg1shinyapp

