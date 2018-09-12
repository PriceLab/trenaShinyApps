FROM pshannon/trenaprereqsdata:latest
USER root
RUN apt-get update
RUN apt-get install libxml2-dev -y
RUN apt-get -y install libmariadb-client-lgpl-dev -y
RUN apt-get install libpq-dev -y
RUN sudo apt-get install libpng-dev -y
RUN sudo apt-get install zlib1g-dev -y
COPY installTrenaCodePrequisites.R /home/trena
USER trena
COPY myRprofile /home/trena/.Rprofile
WORKDIR /home/trena
RUN R -f ./installTrenaCodePrequisites.R
COPY packages/*gz ./
RUN R CMD INSTALL -l ~/library MotifDb_1.23.11.tar.gz
RUN R CMD INSTALL -l ~/library igvShiny_0.99.48.tar.gz
RUN R CMD INSTALL -l ~/library cyjShiny_0.99.8.tar.gz
RUN R CMD INSTALL -l ~/library trena_1.3.7.tar.gz
