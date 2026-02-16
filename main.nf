#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.outdir  = params.outdir  ?: "results"

/*
 * 1) Téléchargement des lectures
 */
process DOWNLOAD_FASTQ {
    tag "$id"
    publishDir "${params.outdir}/raw_fastq", mode: 'copy'

    input:
        val id

    output:
        tuple val(id), path("${id}.fastq"), emit: fastq
        path "versions.yml", emit: versions

    script:
    """
    fasterq-dump ${id} --threads ${task.cpus} -O .
    echo "fasterq-dump: `fasterq-dump --version 2>&1 | head -1`" > versions.yml
    """

    stub:
    """
    touch ${id}.fastq
    echo "fasterq-dump: stub" > versions.yml
    """
}

/*
 * 2) Trimming (Cutadapt)
 */
process TRIM {
    tag "$id"
    publishDir "${params.outdir}/trimmed", mode: 'copy'
    container 'alantrbt/cutadapt:1.11'

    input:
        tuple val(id), path(read)

    output:
        tuple val(id), path("${id}.trimmed.fastq"), emit: fastq
        path "versions.yml", emit: versions

    script:
    """
    cutadapt -q 20 -m 4 --length 25 \
        -o ${id}.trimmed.fastq \
        ${read}
    echo "cutadapt: `cutadapt --version`" > versions.yml
    """

    stub:
    """
    touch ${id}.trimmed.fastq
    echo "cutadapt: stub" > versions.yml
    """
}

/*
 * 3) Téléchargement génome
 */
process DOWNLOAD_REFERENCE {
    publishDir "data/reference", mode: 'copy'

    output:
        path "reference.fasta", emit: fasta
        path "reference.gff3",  emit: gff
        path "versions.yml", emit: versions

    script:
    """
    wget -q -O reference.gff3 \
        "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?db=nuccore&report=gff3&id=CP000253.1"
    wget -q -O reference.fasta "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=CP000253.1&rettype=fasta"
    echo "wget: `wget --version | head -1`" > versions.yml
    """

    stub:
    """
    touch reference.fasta reference.gff3
    echo "wget: stub" > versions.yml
    """
}

/*
 * 4) Index Bowtie1
 */
process INDEX {
    tag "bowtie-index"
    publishDir "${params.outdir}/index", mode: 'copy'

    input:
        path fasta

    output:
        path "index.*", emit: index
        path "versions.yml", emit: versions

    script:
    """
    bowtie-build ${fasta} index
    echo "bowtie-build: `bowtie-build --version | head -1`" > versions.yml
    """

    stub:
    """
    touch index.1 index.2 index.3 index.4
    echo "bowtie-build: stub" > versions.yml
    """
}

/*
 * 5) Alignement
 */
process ALIGN {
    tag "$id"
    publishDir "${params.outdir}/aligned", mode: 'copy'

    input:
        tuple val(id), path(read)
        path index_files

    output:
        tuple val(id), path("${id}.bam"), emit: bam
        path "versions.yml", emit: versions

    script:
    """
    bowtie -p ${task.cpus} --sam index ${read} > ${id}.sam
    samtools view -bS ${id}.sam > ${id}.bam
    samtools sort ${id}.bam ${id}.sorted
    mv ${id}.sorted.bam ${id}.bam
    samtools index ${id}.bam
    echo "bowtie: `bowtie --version | head -1`" > versions.yml
    echo "samtools: `samtools --version | head -1`" >> versions.yml
    """

    stub:
    """
    touch ${id}.bam
    echo "bowtie: stub" > versions.yml
    echo "samtools: stub" >> versions.yml
    """
}

/*
 * 6) Comptage (featureCounts, GFF3)
 */
process COUNT {
    tag "$id"
    publishDir "${params.outdir}/counts", mode: 'copy'

    input:
        tuple val(id), path(bam)
        path gff

    output:
        path "counts_${id}.txt", emit: counts
        path "versions.yml", emit: versions

    script:
    """
    featureCounts \
        -F GFF \
        -t gene \
        -g ID \
        -T ${task.cpus} \
        -a ${gff} \
        -o counts_${id}.txt \
        ${bam}
    echo "featureCounts: `featureCounts -v 2>&1 | head -1`" > versions.yml
    """

    stub:
    """
    touch counts_${id}.txt
    echo "featureCounts: stub" > versions.yml
    """
}

/*
 * 7) DESeq2
 */
process DESEQ2 {
    publishDir "${params.outdir}/deseq2", mode: 'copy'

    input:
        path counts_files
        val samples_metadata

    output:
        path "deseq2_results.csv", emit: results
        path "versions.yml", emit: versions

    script:
    """
    Rscript /scripts/run_deseq2.R $baseDir/samples.tsv ${counts_files.join(" ")} deseq2_results.csv
    echo "Rscript: `Rscript --version | head -1`" > versions.yml
    """

    stub:
    """
    touch deseq2_results.csv
    echo "Rscript: stub" > versions.yml
    """
}

/*
 * 8) Annotation GFF → DESeq2 annoté
 */
process ANNOTATE_GENES {

    tag "annotate"
    publishDir "${params.outdir}/deseq2", mode: 'copy'
    container "bioconductor/bioconductor_docker:RELEASE_3_17"

    input:
        path deseq
        path gff

    output:
        path "deseq2_results_annotated.csv", emit: annotated
        path "versions.yml", emit: versions

    script:
    """
    Rscript /scripts/annotate_genes.R \
        $deseq \
        $gff \
        deseq2_results_annotated.csv
    echo "Rscript: `Rscript --version | head -1`" > versions.yml
    """
}

/*
 * 9) Analyse des pathways KEGG
 */
process PATHWAYS {

    tag "pathways"
    publishDir "${params.outdir}/pathways", mode: 'copy'
    container "bioconductor/bioconductor_docker:RELEASE_3_17"

    input:
        path annotated

    output:
        path "kegg_pathways_up.csv"
        path "kegg_pathways_down.csv"
        path "kegg_up_barplot.png"
        path "kegg_down_barplot.png"
        path "versions.yml", emit: versions

    script:
    """
    Rscript /scripts/run_pathways.R \
        $annotated \
        pathways_results

    cp pathways_results/* .
    echo "Rscript KEGG: `Rscript --version | head -1`" > versions.yml
    """
}

/*
 * 10) Visualisation des résultats DESeq2
 */
process PLOT_DESEQ2 {
    publishDir "${params.outdir}/deseq2/plots", mode: 'copy'
    container "alantrbt/deseq2:latest"

    input:
        path deseq_results

    output:
        path "deseq2_plot.png", emit: plot
        path "versions.yml", emit: versions

    script:
    """
    # Correction : si 'gene_name' n'existe pas, on remplace par 'gene_id' dans le script
    sed 's/gene_name/gene_id/g' /scripts/plot_deseq2.R > plot_deseq2_tmp.R
    Rscript plot_deseq2_tmp.R $deseq_results deseq2_plot.png
    echo "Rscript: `Rscript --version | head -1`" > versions.yml
    """
}

/*
 * Workflow principal
 */
workflow {

    def sample_list = file("$baseDir/samples.tsv").splitCsv(header:true, sep:'\t')

    samples_ch = Channel.from(sample_list)
    sra_ids_ch = samples_ch.map { row -> row.sra }

    reads_ch = DOWNLOAD_FASTQ(sra_ids_ch).fastq
    trimmed_ch = TRIM(reads_ch).fastq

    ref        = DOWNLOAD_REFERENCE()
    fasta_ch   = ref.fasta
    gff_ch     = ref.gff

    index_ch   = INDEX(fasta_ch).index
    aligned_ch = ALIGN(trimmed_ch, index_ch).bam
    counts_ch  = COUNT(aligned_ch, gff_ch).counts

    all_counts = counts_ch.collect()

    samples_ch.collect().set { samples_metadata }
    deseq_results = DESEQ2(all_counts, samples_metadata).results
    annotated_ch  = ANNOTATE_GENES(deseq_results, gff_ch).annotated
    PATHWAYS(annotated_ch)
    PLOT_DESEQ2(deseq_results)
}




