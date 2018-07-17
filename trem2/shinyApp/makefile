runInShell:
	R -f runTrem2InBackground.R

runInBackground:
	nohup R -r runTrem2InBackground.R &


build:
	docker build -t trem2shinyapp -f Dockerfile .

bash:
	docker run -ti --rm -p 60040:3838 trem2shinyapp bash


