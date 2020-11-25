FROM rocker/tidyverse:4.0.3

RUN apt-get -y install libcurl4-gnutls-dev libssl-dev libgmp3-dev libxml2-dev

WORKDIR /lib
COPY DESCRIPTION /lib/DESCRIPTION
RUN Rscript -e 'devtools::install()'
COPY . /lib

ENTRYPOINT ["Rscript"]