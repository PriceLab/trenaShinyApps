TIMEOUT = 12000

all: flaskServer chia cellphonesim inpp5d frd3 ireg1

flaskServer:
	R -f test.igv-data.flaskServer.R

chia:
	testcafe chrome chia.js --assertion-timeout 6000

cellphonesim:
	testcafe chrome cellphonesim.js --assertion-timeout 6000

frd3:
	testcafe chrome frd3.js --assertion-timeout $(TIMEOUT)

ireg1:
	testcafe chrome ireg1.js --assertion-timeout $(TIMEOUT)

inpp5d:
	testcafe chrome inpp5d.js --assertion-timeout $(TIMEOUT)


