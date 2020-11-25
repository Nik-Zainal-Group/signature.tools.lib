FROM r-base:4.0.3

COPY DESCRIPTION /lib/DESCRIPTION
WORKDIR /lib
RUN Rscript -e 'install.packages("devtools")'
RUN Rscript -e 'devtools::install()'
COPY . /lib

ENTRYPOINT ["Rscript"]