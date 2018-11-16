build:
	docker build -t trenaigap:test -f Dockerfile .

bash:
	docker run -ti --rm -p 60043:3838 trenaigap:test bash

dockerRun:
	docker run -p 60043:3838 trenaigap:test

dockerView:
	open http://localhost:60043

save:
	time docker save trenaigap:latest >  ~/tmp/dockerIGAP.save
	ls -lat ~/tmp/dockerIGAP.save

saveTest:
	time docker save trenaigap:test >  ~/tmp/dockerIGAP.test.save
	ls -lat ~/tmp/dockerIGAP.test.save

scp:
	time scp ~/tmp/dockerIGAP.save pshannon@trena:/local/tmp/
	ssh pshannon@trena ls -l /local/tmp/dockerIGAP.save
#	ssh pshannon@trena docker load < /local/tmp/dockerIGAP.save

scpTest:
	time scp ~/tmp/dockerIGAP.test.save pshannon@trena:/local/tmp/
	ssh pshannon@trena ls -l /local/tmp/dockerIGAP.test.save

