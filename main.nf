#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.indir ="./input/"
params.output ="./results/"

// Defaults
params.publishmode = 'symlink'

process FASTQC_v0_11_9 {
    tag "${ID} quality check"
    container 'quay.io/biocontainers/fastqc:0.11.9--hdfd78af_1'
    publishDir { "${params.output}/${ID}/QC_${ID}/fastqc_${ID}" }, pattern: "*.{zip,html}", mode: params.publishmode
    input: tuple val(ID), path(r1), path(r2)
    output:
        path("*.zip"), emit: qc
        path("*.html"), emit: html
        path("1_raw_${ID}.count"), emit: countfile
    shell:
        '''
        echo "Raw_readcount,$(expr $(zcat !{r1} | wc -l) / 4),!{ID}" >> "1_raw_!{ID}.count"
        fastqc -t !{task.cpus} -q !{r1} !{r2} --noextract
        ''' }

process TRIMMOMATIC_v0_32 {
    tag "${ID} quality trimming"
    container 'quay.io/biocontainers/trimmomatic:0.32--hdfd78af_4'
    publishDir { "${params.output}/${ID}/reads_${ID}" }, pattern: "*fastq.gz", mode: params.publishmode
    input: tuple val(ID), path(r1), path(r2)
    output:
        tuple val(ID), path("trimmed_${ID}_r1.fastq.gz"), path("trimmed_${ID}_r2.fastq.gz"), emit: reads
        path("2_trimmed_${ID}.count"), emit: countfile
    shell:
        '''
        trimmomatic PE -phred33 -threads !{task.cpus} !{r1} !{r2} \
        trimmed_!{ID}_r1.fastq.gz orphan_!{ID}_r1.fastq.gz trimmed_!{ID}_r2.fastq.gz orphan_!{ID}_r2.fastq.gz LEADING:0 TRAILING:5 SLIDINGWINDOW:4:15 MINLEN:36
        echo "Trimmed_readcount,$(expr $(zcat trimmed_!{ID}_r1.fastq.gz | wc -l) / 4),!{ID}" > "2_trimmed_!{ID}.count"
        ''' }

workflow {

    // Check inputs
    def input_dir = file(params.indir)
    if (!input_dir.exists()) {
        error("Input directory not found: ${params.indir}.\nMake sure it exists and contains paired Illumina fastq files, named like.\n  BLAH_R1.fastq.gz\n  BLAH_R2.fastq.gz")
}
    CH_fastq = channel.fromPath("${params.indir}/*.{fastq,fastq.gz,fq,fq.gz}", checkIfExists: true)
        .flatten() // emit each fastq path as its own item
        .map { file -> tuple(file.getSimpleName().replaceFirst('_R1','').replaceFirst('_R2',''), file) }    // Derive library by removing '_R1' or '_R2' suffixes. E.g. 4_12345678_R1.fastq.gz -> [4, 4_12345678_R1.fastq.gz] 
        .groupTuple(size:2) // E.g. [4, [4_12345678_R1.fastq.gz, 4_12345678_R2.fastq.gz]]
        .map { it -> tuple(it[0], it[1][0], it[1][1]) }       // Lastly, simplify.                              E.g. [X,[Y,Z]] ->  [X, Y, Z]
    
    FASTQC_v0_11_9(CH_fastq)
    TRIMMOMATIC_v0_32(CH_fastq)
}

