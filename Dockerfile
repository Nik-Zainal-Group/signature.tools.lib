FROM r-base:4.0.3

RUN Rscript -e 'install.packages("devtools")'

WORKDIR /lib
COPY DESCRIPTION /lib/DESCRIPTION
RUN Rscript -e 'devtools::install()'
COPY . /lib

ENTRYPOINT ["Rscript"]