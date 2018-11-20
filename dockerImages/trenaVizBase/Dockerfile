FROM pshannon/trenashinybase
USER root
COPY installCodePrerequisites.R /home/trena
USER trena
WORKDIR /home/trena
RUN R -f ./installCodePrerequisites.R
COPY packages/*gz ./
RUN R CMD INSTALL -l ~/library TrenaProject_0.99.17.tar.gz
# USER root

