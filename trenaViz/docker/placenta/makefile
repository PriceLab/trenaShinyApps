build:
	docker build -t trenaplacenta -f Dockerfile .

bash:
	docker run -ti --rm -p 60043:3838 trenaplacenta bash

dockerRun:
	docker run -p 60043:3838 trenaplacenta

dockerView:
	open http://localhost:60043

save:
	time docker save trenaplacenta:latest >  ~/tmp/dockerPLACENTA.save
	ls -lat ~/tmp/dockerPLACENTA.save

saveTest:
	time docker save trenaplacenta >  ~/tmp/dockerPLACENTA.test.save
	ls -lat ~/tmp/dockerPLACENTA.test.save

scp:
	time scp ~/tmp/dockerPLACENTA.save pshannon@trena:/local/tmp/
	ssh pshannon@trena ls -l /local/tmp/dockerPLACENTA.save
#	ssh pshannon@trena docker load < /local/tmp/dockerPLACENTA.save

scpTest:
	time scp ~/tmp/dockerPLACENTA.test.save pshannon@trena:/local/tmp/
	ssh pshannon@trena ls -l /local/tmp/dockerPLACENTA.test.save

