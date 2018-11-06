# runInShell:
# 	R -f runInpp5dInBackground.R
#
# runInBackground:
# 	nohup R -r runInpp5dInBackground.R &
#
#
build:
	docker build -t trena_skin -f Dockerfile .

bash:
	docker run -ti --rm -p 60040:3838 trena_skin bash

dockerRun:
	docker run -p 60042:3838 inpp5dshinyapp

dockerView:
	open http://localhost:60042
