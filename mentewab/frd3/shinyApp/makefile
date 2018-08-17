runInShell:
	R -f runFrd3InBackground.R

runInBackground:
	nohup R -r runFrd3InBackground.R &


build:
	docker build -t frd3shinyapp -f Dockerfile .

bash:
	docker run -ti --rm -p 60039:3838 frd3shinyapp bash

dockerRun:
	docker run -p 60039:3838 frd3shinyapp

