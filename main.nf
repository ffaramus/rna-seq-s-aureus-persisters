#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.outdir      = params.outdir ?: "results"
params.paper_table = params.paper_table ?: "data/paper/paper_results.tsv"
params.mapping     = params.mapping ?: "data/mapping/mapping.tsv"

/*
===========================================
 1) DOWNLOAD FASTQ
===========================================
*/
process DOWNLOAD_FASTQ {

    cpus { Math.max(1, (params.core_nb / 6).toInteger()) }

    input:
        val id

    output:
        tuple val(id), path("${id}.fastq"), emit: fastq
        path "versions.yml", emit: versions

    publishDir "${params.outdir}/raw_fastq", mode:'copy'

    script:
    """
    fasterq-dump ${id} --threads ${task.cpus} -O .
    echo "fasterq-dump: \$(fasterq-dump --version 2>&1 | head -1)" > versions.yml
    """
}


/*
===========================================
 2) TRIMMING
===========================================
*/
process TRIM {

    cpus { Math.max(1, (params.core_nb / 6).toInteger()) }
    container 'alantrbt/cutadapt:1.11'

    input:
        tuple val(id), path(read)

    output:
        tuple val(id), path("${id}.trimmed.fastq"), emit: fastq
        path "versions.yml", emit: versions

    publishDir "${params.outdir}/trimmed", mode: 'copy'

    script:
    """
    cutadapt -q 20 -m 4 --length 25 \
        -o ${id}.trimmed.fastq \
        ${read}

    echo "cutadapt: \$(cutadapt --version)" > versions.yml
    """
}


/*
===========================================
 3) REFERENCE GENOME
===========================================
*/
process DOWNLOAD_REFERENCE {

    publishDir "data/reference", mode: 'copy'

    output:
        path "reference.fasta", emit: fasta
        path "reference.gff3", emit: gff

    script:
    """
    wget -q -O reference.gff3 \
      "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?db=nuccore&report=gff3&id=CP000253.1"

    wget -q -O reference.fasta \
      "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=CP000253.1&rettype=fasta"
    """
}


/*
===========================================
 4) INDEX BOWTIE
===========================================
*/
process INDEX {

    publishDir "${params.outdir}/index", mode:'copy'

    input:
        path fasta

    output:
        path "index.*", emit: index

    script:
    """
    bowtie-build ${fasta} index
    """
}


/*
===========================================
 5) ALIGNMENT
===========================================
*/
process ALIGN {

    cpus { Math.max(1, (params.core_nb / 4).toInteger()) }
    publishDir "${params.outdir}/aligned", mode:'copy'

    input:
        tuple val(id), path(read)
        path idx

    output:
        tuple val(id), path("${id}.bam"), emit: bam

    script:
    """
    bowtie -p ${task.cpus} --sam index ${read} > ${id}.sam
    samtools view -bS ${id}.sam > ${id}.bam
    samtools sort ${id}.bam -o ${id}.sorted.bam
    mv ${id}.sorted.bam ${id}.bam
    samtools index ${id}.bam
    """
}


/*
===========================================
 6) FEATURECOUNTS
===========================================
*/
process FEATURECOUNTS {

    publishDir "${params.outdir}/counts", mode:'copy'
    container 'alantrbt/subread:latest'

    input:
        path bam_files
        path gff

    output:
        path "counts_matrix.txt", emit: matrix

    script:
    """
    featureCounts \
        -F GFF \
        -t gene \
        -g ID \
        -T ${task.cpus} \
        -a ${gff} \
        -o counts_matrix.txt \
        ${bam_files.join(' ')}
    """
}


/*
===========================================
 7) DESEQ2  
===========================================
*/
process DESEQ2 {

    publishDir "${params.outdir}/deseq2", mode:'copy'
    container "alantrbt/deseq2:latest"

    input:
        path counts_matrix
        path samples_file

    output:
        path "deseq2_results.csv", emit: results

    script:
    """
    Rscript /scripts/run_deseq2.R \
        ${samples_file} \
        ${counts_matrix} \
        deseq2_results.csv
    """
}


/*
===========================================
 8) KEGG TABLE  (DOWNSTREAM CONTAINER)
===========================================
*/
process CREATE_GENE_PATHWAY_TABLE {
    publishDir "${params.outdir}/downstream", mode:'copy'
    container "alantrbt/downstream:latest"

    input:
        val dummy

    output:
        path "gene_pathway_table.tsv"

    script:
    """
    Rscript /scripts/1_Create_gene_pathway.R .
    """
}


/*
===========================================
 9) DOWNSTREAM PLOTS  (DOWNSTREAM CONTAINER)
===========================================
*/
process CREATE_PLOTS {
    publishDir "${params.outdir}/downstream", mode:'copy'
    container "alantrbt/downstream:latest"

    input:
        path deseq_results
        path gene_pathway
        path mapping_file

    output:
        path "*.pdf"
        path "*.png"

    script:
    """
    Rscript /scripts/2_Create_plots.R \
        ${deseq_results} \
        ${gene_pathway} \
        ${mapping_file}
    """
}


/*
===========================================
 10) PAPER COMPARISON  (DOWNSTREAM CONTAINER)
===========================================
*/
process COMPARE_WITH_PAPER {
    publishDir "${params.outdir}/downstream/paper_comparison", mode:'copy'
    container "alantrbt/downstream:latest"

    input:
        path paper_table
        path deseq_results
        path counts_matrix
        path mapping_file

    output:
        path "*.pdf"
        path "*.tsv"

    script:
    """
    Rscript /scripts/3_Compare.R \
        ${paper_table} \
        ${deseq_results} \
        ${counts_matrix} \
        ${mapping_file}
    """
}


/*
===========================================
 11) BASIC DESEQ2 PLOT
===========================================
*/
process PLOT_DESEQ2 {
    publishDir "${params.outdir}/deseq2/plots", mode:'copy'
    container "alantrbt/deseq2:latest"

    input:
        path deseq_results

    output:
        path "deseq2_plot.png"

    script:
    """
    Rscript /scripts/plot_deseq2.R \
        ${deseq_results} \
        deseq2_plot.png
    """
}


/*
===========================================
 WORKFLOW
===========================================
*/
workflow {

    def sample_list = file("$baseDir/samples.tsv").splitCsv(header:true, sep:'\t')

    samples_ch  = Channel.from(sample_list)
    sra_ids_ch  = samples_ch.map { it.sra }

    fastq       = DOWNLOAD_FASTQ(sra_ids_ch).fastq
    trimmed     = TRIM(fastq).fastq
    ref         = DOWNLOAD_REFERENCE()
    idx         = INDEX(ref.fasta).index
    aligned     = ALIGN(trimmed, idx).bam
    matrix      = FEATURECOUNTS(aligned, ref.gff).matrix

    samples_x   = Channel.of(file("$baseDir/samples.tsv"))
    deseq_res   = DESEQ2(matrix, samples_x).results

    PLOT_DESEQ2(deseq_res)

    pathways = CREATE_GENE_PATHWAY_TABLE(Channel.value(true))

    CREATE_PLOTS(
        deseq_res,
        pathways,
        file(params.mapping)
    )

    COMPARE_WITH_PAPER(
        file(params.paper_table),
        deseq_res,
        matrix,
        file(params.mapping)
    )
}

