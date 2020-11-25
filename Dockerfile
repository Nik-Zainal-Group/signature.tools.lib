FROM r-base:4.0.3

WORKDIR /lib
RUN Rscript -e 'install.packages("devtools")'

COPY DESCRIPTION /lib/DESCRIPTION
RUN Rscript -e 'devtools::install()'
COPY . /lib

ENTRYPOINT ["Rscript"]