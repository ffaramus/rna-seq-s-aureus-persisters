FROM rocker/r-base:4.3.2

RUN apt-get update && apt-get install -y \
    libcurl4-openssl-dev \
    libxml2-dev \
    libssl-dev \
    && rm -rf /var/lib/apt/lists/*

RUN R -e "install.packages('BiocManager', repos='https://cloud.r-project.org')"
RUN R -e "BiocManager::install('KEGGREST')"
RUN R -e "install.packages(c('ggplot2','ggrepel','VennDiagram'), repos='https://cloud.r-project.org')"

WORKDIR /data

RUN mkdir -p /scripts
COPY *.R /scripts/

