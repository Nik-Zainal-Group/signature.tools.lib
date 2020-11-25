FROM rocker/r-apt:bionic

RUN apt-get update && \
    apt-get install -y -qq \
    	r-cran-devtools \
      r-cran-nmf \
      r-cran-foreach \
      r-cran-doparallel \
      r-cran-domc \
      r-cran-lpsolve \
      r-cran-ggplot2 \
      r-cran-cluster \
      r-cran-nnls \
      r-cran-gmp \
      r-cran-plyr \
      r-cran-scales \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/

WORKDIR /lib
COPY DESCRIPTION /lib/DESCRIPTION
RUN Rscript -e 'devtools::install()'
COPY . /lib

ENTRYPOINT ["Rscript"]