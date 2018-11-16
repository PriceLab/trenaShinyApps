build:
	docker build -t trenaigap:test -f Dockerfile .

bash:
	docker run -ti --rm -p 60042:3838 trenaigap bash

dockerRun:
	docker run -p 60042:3838 trenaigap:test

dockerView:
	open http://localhost:60042

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

