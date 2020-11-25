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
      r-cran-farver \
      r-cran-rcpp \
      r-cran-xtable \
      r-cran-bitops \
      r-cran-formatr \
      r-cran-futile.options \
      r-cran-snow \
      r-cran-futile.logger \
      r-cran-matrixstats \
      r-cran-generics \
      r-cran-tidyselect \
      r-cran-hms \
      r-cran-bit \
      r-cran-curl \
      r-cran-dplyr \
      r-cran-dbplyr \
      r-cran-plogr \
      r-cran-blob \
      r-cran-bit64 \
      r-cran-rappdirs \
      r-cran-progress \
      r-cran-xml \
      r-cran-rcurl \
      r-cran-rsqlite \
      r-cran-dbi \
      r-cran-readr \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/

WORKDIR /lib
COPY DESCRIPTION /lib/DESCRIPTION
RUN Rscript -e 'devtools::install()'
COPY . /lib

ENTRYPOINT ["Rscript"]