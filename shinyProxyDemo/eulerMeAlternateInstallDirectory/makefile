build:
	docker build -t eulershiny3 ./

bash:
	docker run -ti --rm -p 60021:3838 eulershiny3 bash
            # start R
            # shiny::runApp('/root/euler', launch.browser=FALSE)

runContainer:
	docker run -ti --rm -p 60021:3838 eulershiny3

runContainerDetached:
	docker run -ti --rm -d -p 60021:3838 eulershiny3
           # docker ps to find the containerID
           # docker kill containerID


