FROM pshannon/trenashinybase
USER root
COPY Rprofile.site /usr/local/lib/R/etc/
USER trena
RUN mkdir -p /home/trena/app/data
#RUN mkdir -p /home/trena/app/tmp
#RUN mkdir -p /home/trena/app/png
WORKDIR /home/trena/app
COPY data/* /home/trena/app/data/
#COPY png/* /home/trena/app/png/
EXPOSE 3838
#COPY message-handler.js .
#COPY readme.html .
COPY app.R .
CMD ["R", "-e", "shiny::runApp('/home/trena/app', port=3838, host='0.0.0.0')"]





