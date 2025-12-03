#!/usr/bin/env bash

# ============================================================
# Script d'exécution du pipeline RNA-seq (Nextflow + Docker)
# Hackathon Reproductibilité – M2 AMI2B 2025
# ============================================================

set -e

echo "----------------------------------------"
echo " Build des conteneurs Docker"
echo "----------------------------------------"

docker build -t alantrbt/bowtie:latest      containers/bowtie
docker build -t alantrbt/cutadapt:1.11     containers/cutadapt
docker build -t alantrbt/sratoolkit:latest containers/sratoolkit
docker build -t alantrbt/subread:latest    containers/subread
docker build -t alantrbt/deseq2:latest     containers/deseq2
docker build -t alantrbt/tidyverse:latest  containers/tydiverse

echo ""
echo "----------------------------------------"
echo " Lancement du pipeline Nextflow"
echo "----------------------------------------"

nextflow run main.nf -with-docker -resume

echo ""
echo "Pipeline terminé avec succès."
