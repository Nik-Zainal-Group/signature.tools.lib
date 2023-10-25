FROM rocker/tidyverse:4.3

WORKDIR /root/signature.tools.lib
COPY . .
RUN apt-get -y update && \
  apt -y install git tar libgit2-dev libgmp-dev libbz2-dev libcurl4-openssl-dev libssh-dev libssl-dev libxml2-dev libharfbuzz-dev libfribidi-dev libfreetype-dev libpng-dev libtiff-dev libjpeg-dev && \
  rm -rf /var/lib/apt/lists/*
RUN R -e 'devtools::install()'