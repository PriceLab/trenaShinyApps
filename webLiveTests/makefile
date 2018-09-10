TIMEOUT = 6000

all: frd3 ireg1

frd3:
	testcafe chrome frd3.js --assertion-timeout $(TIMEOUT)

ireg1:
	testcafe chrome ireg1.js --assertion-timeout $(TIMEOUT)
