FROM mambaorg/micromamba:1.5.8

LABEL software="deseq2"
LABEL version="1.16"
LABEL description="DESeq2 1.16 environment (R 3.4.1)"

ARG MAMBA_DOCKERFILE_ACTIVATE=1

USER root

RUN micromamba install -y -n base -c conda-forge -c bioconda \
        r-base=3.4.1 \
        bioconductor-deseq2=1.16.* \
        r-ggplot2=3.0.* \
        libgfortran=3.0.0 \
        libgcc-ng=7.* \
        libstdcxx-ng=7.* \
    && micromamba clean -a -y


RUN apt-get update && apt-get install -y procps libxml2-dev libcurl4-openssl-dev libssl-dev && rm -rf /var/lib/apt/lists/*

# Installation des packages R figés depuis les archives CRAN
RUN R -e "install.packages('https://cran.r-project.org/src/contrib/Archive/ggplot2/ggplot2_3.0.0.tar.gz', repos=NULL, type='source')"
RUN R -e "install.packages('https://cran.r-project.org/src/contrib/Archive/dplyr/dplyr_0.7.8.tar.gz', repos=NULL, type='source')"
RUN R -e "install.packages('https://cran.r-project.org/src/contrib/Archive/ggrepel/ggrepel_0.8.2.tar.gz', repos=NULL, type='source')"

ENV PATH="/opt/conda/envs/base/bin:/opt/conda/bin:$PATH"

WORKDIR /data

CMD ["bash"]
