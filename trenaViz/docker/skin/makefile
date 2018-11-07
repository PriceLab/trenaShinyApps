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
	docker run -ti --rm -p 60042:3838 trena_skin bash

dockerRun:
	docker run -p 60042:3838 trena_skin

dockerView:
	open http://localhost:60042

save:
	time docker save trena_skin:latest >  ~/tmp/dockerSkin.save
	ls -lat ~/tmp/dockerSkin.save

scp:
	time scp ~/tmp/dockerSkin.save pshannon@trena:/local/tmp/
	ssh pshannon@trena ls -l /local/tmp/dockerSkin.save
	ssh pshannon@trena docker load < /local/tmp/dockerSkin.save

