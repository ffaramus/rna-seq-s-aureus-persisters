FROM mambaorg/micromamba:1.5.8

# Installer R 3.4.1 + DESeq2 1.16 + ggplot2
RUN micromamba install -y -n base -c conda-forge -c bioconda \
    r-base=3.4.1 \
    bioconductor-deseq2=1.16.* \
    r-ggplot2=3.0.* \
    libgfortran=3.0.0 \
    libgcc-ng=7.* \
    libstdcxx-ng=7.* \
 && micromamba clean -a -y

# Tous les RUN suivants utiliseront bash
SHELL ["/bin/bash", "-lc"]

# Installer packages R via CRAN (nécessite activation micromamba)
RUN eval "$(micromamba shell hook -s bash)" \
 && micromamba activate base \
 && echo 'install.packages(c("KEGGREST","ggrepel","VennDiagram"), repos="https://cloud.r-project.org")' > /tmp/install.R \
 && Rscript /tmp/install.R \
 && rm /tmp/install.R

WORKDIR /data

