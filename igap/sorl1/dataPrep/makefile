PORT=5558
FIMO="/Users/paul/meme/bin/fimo"
MEME_FILE=sorl1.tf.motfs.meme
default:  kill.fimo run.fimo

kill.fimo:
	- kill -9 `ps x | grep runFimoServer.py | egrep -v grep | awk '{print $$1}'`

run.fimo:
	python -i ~/github/fimoService/server/runFimoServer.py $(PORT) $(FIMO) $(MEME_FILE)
