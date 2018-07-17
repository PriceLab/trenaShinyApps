FROM pshannon/trenashinybase
USER trena
RUN mkdir -p /home/trena/app/data
RUN mkdir -p /home/trena/app/tmp
RUN mkdir -p /home/trena/app/png
WORKDIR /home/trena/app
COPY data/* /home/trena/app/data/
COPY png/* /home/trena/app/png/
COPY message-handler.js .
COPY readme.html .
COPY app.R .





