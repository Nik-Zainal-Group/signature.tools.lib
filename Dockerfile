FROM r-base:4.0.3

RUN apt-get update -qq \ 
  && apt-get install -t unstable -y --no-install-recommends \
    libcurl4-openssl-dev \
    libgmp3-dev \
    libssl-dev \
    libsqlite3-dev \
    libxml2-dev \
  && apt-get clean \
  && rm -rf /var/lib/apt/lists/ \
  && rm -rf /tmp/downloaded_packages/ /tmp/*.rds
RUN Rscript -e 'install.packages("devtools")'

WORKDIR /lib
COPY DESCRIPTION /lib/DESCRIPTION
RUN Rscript -e 'devtools::install()'
COPY . /lib

ENTRYPOINT ["Rscript"]