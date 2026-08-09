#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.indir ="./input/"
params.output ="./results/"

// MODES
params.dev = false

// Defaults
params.publishmode = 'symlink'

//# READ PROCESSING
params.phred = "33"
params.kmernorm_opts = "-k 21 -t 30 -c 3"
params.complexity_threshold = "0.05"
params.reference_threshold = "0.05"

//#     DECONTAMINATION
params.contamination_reference = "/mnt/scgc_nfs/ref/assembly_protocol/GRCh38_AG665_mm10.fa"
params.amb = "/mnt/scgc_nfs/ref/assembly_protocol/GRCh38_AG665_mm10.fa.amb"
params.ann = "/mnt/scgc_nfs/ref/assembly_protocol/GRCh38_AG665_mm10.fa.ann"
params.bwt = "/mnt/scgc_nfs/ref/assembly_protocol/GRCh38_AG665_mm10.fa.bwt"
params.pac = "/mnt/scgc_nfs/ref/assembly_protocol/GRCh38_AG665_mm10.fa.pac"
params.sa = "/mnt/scgc_nfs/ref/assembly_protocol/GRCh38_AG665_mm10.fa.sa"

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
    
    CH_fastq = params.dev ? CH_fastq.take(1) : CH_fastq

    FASTQC_v0_11_9(CH_fastq)
    TRIMMOMATIC_v0_32(CH_fastq)
    COMPLEXITY_FILTER(TRIMMOMATIC_v0_32.out.reads)
    KMERNORM_v1_1_0(COMPLEXITY_FILTER.out.reads)
    
    // LOG_COMPLEX_READS(KMERNORM_v1_0_0.out.plex_countfile.collect())
	// LOG_NORMALIZED_READS(KMERNORM_v1_0_0.out.norm_countfile.collect())
	DEINTERLEAVE(KMERNORM_v1_0_0.out.reads_gz)
	// CONTAM_READ_FINDER(DEINTERLEAVE.out.reads_gz)
	// CONTAM_READ_REMOVER(CONTAM_READ_FINDER.out.join(KMERNORM_v1_0_0.out.reads_gz))
	// LOG_CLEAN_READS(CONTAM_READ_REMOVER.out.countfile.collect())
}

process FASTQC_v0_11_9 {
    tag "${ID} quality check"
    container 'quay.io/biocontainers/fastqc:0.11.9--hdfd78af_1'
    publishDir { "${params.output}/${ID}/QC_${ID}/fastqc_${ID}" }, pattern: "*.{zip,html}", mode: params.publishmode
    publishDir { "${params.output}/sample_tracking" }, pattern: "*count", mode: params.publishmode
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
    publishDir { "${params.output}/sample_tracking" }, pattern: "*count", mode: params.publishmode
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

process COMPLEXITY_FILTER {
    tag "${ID}"
    container 'complexity-filter-env:latest'
    //container='brwnj/kmernorm:v1.0.0'
    publishDir { "${params.output}/${ID}/reads_${ID}" }, pattern: "*fastq.gz", mode: params.publishmode
    input: tuple val(ID), path(r1), path(r2)
    output: tuple val(ID), path("pe_${ID}.fastq.gz"), emit: reads
    script: template 'complexity_filter.py' }

process KMERNORM_v1_1_0 {
    tag "${ID}"
    //memory = { 16.GB * task.attempt }
    //errorStrategy = {task.attempt <= 3 ? 'retry' : 'ignore'}
    //maxRetries = 3
    container 'brwnj/kmernorm:v1.1.0'
    publishDir { "${params.output}/${ID}/reads_${ID}"}, pattern: "normalized_pe_*.fastq.gz", mode: params.publishmode
    publishDir { "${params.output}/sample_tracking" }, pattern: "*count", mode: params.publishmode
    input: tuple val(ID), path(paired)
    output:
        tuple val(ID), path("normalized_pe_${ID}.fastq.gz"), emit: reads
        path("3_pe_${ID}.count"), emit: plex_countfile
        path("4_normalized_pe_${ID}.count"), emit: norm_countfile
    shell:
        '''
        gunzip -c !{paired} > temp_paired.fastq
        echo "Complexity_filtered_readcount,$(expr $(cat temp_paired.fastq | wc -l) / 8),!{ID}" >> "3_pe_!{ID}.count"
        kmernorm !{params.kmernorm_opts} temp_paired.fastq > normalized_pe_!{ID}.fastq
        echo "Normalized_readcount,$(expr $(cat normalized_pe_!{ID}.fastq | wc -l) / 8),!{ID}" >> "4_normalized_pe_!{ID}.count"
        gzip normalized_pe_!{ID}.fastq

        # Cleanup
        rm temp_paired.fastq
        ''' }

