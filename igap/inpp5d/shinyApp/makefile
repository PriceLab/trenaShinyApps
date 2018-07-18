runInShell:
	R -f runInpp5dInBackground.R

runInBackground:
	nohup R -r runInpp5dInBackground.R &


build:
	docker build -t inpp5dshinyapp -f Dockerfile .

bash:
	docker run -ti --rm -p 60040:3838 inpp5dshinyapp bash

dockerRun:
	docker run -p 60041:3838 inpp5dshinyapp

