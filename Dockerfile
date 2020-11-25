FROM rocker/tidyverse:4.0.3

WORKDIR /lib
COPY DESCRIPTION /lib/DESCRIPTION
RUN Rscript -e 'devtools::install()'
COPY . /lib

ENTRYPOINT ["Rscript"]