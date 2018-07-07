FROM bioconductor/devel_core2:latest
USER root
RUN useradd -ms /bin/bash trena
user trena
WORKDIR /home/trena
RUN mkdir /home/trena/work
RUN mkdir /home/trena/library
COPY myBashCustomizations /home/trena/
COPY myRprofile /home/trena/
COPY installTrenaPrequisites.R /home/trena
#COPY SNPlocs.Hsapiens.dbSNP150.GRCh38_0.99.20.tar.gz .
RUN cat /home/trena/myBashCustomizations >> /home/trena/.bashrc
USER root
RUN R -f installTrenaPrequisites.R
#RUN R CMD INSTALL -l ~/library SNPlocs.Hsapiens.dbSNP150.GRCh38_0.99.20.tar.gz
#RUN rm SNPlocs.Hsapiens.dbSNP150.GRCh38_0.99.20.tar.gz
