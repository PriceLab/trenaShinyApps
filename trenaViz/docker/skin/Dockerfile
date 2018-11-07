FROM pshannon/trenavizbase
USER root
COPY Rprofile.site /usr/local/lib/R/etc/
USER trena
WORKDIR /home/trena/app
COPY packages/*gz ./
RUN R CMD INSTALL -l ~/library TrenaProjectSkin_0.99.1.tar.gz
EXPOSE 3838
COPY trenaViz.R .
CMD ["R", "-e", "shiny::runApp('/home/trena/app/trenaViz.R', port=3838, host='0.0.0.0')"]





