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
params.contam_ref_fasta = "/Users/greggavelis/Desktop/SCGC_Refcontam/GRCh38_AG665_mm10.fa"
// BWA Index is auto-detected next to contam_ref_fasta (<contam_ref_fasta>.amb/.ann/.bwt/.pac/.sa);
// BWA_INDEX only runs when one or more of those files is missing.

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
    KMERNORM_v1_0_0(COMPLEXITY_FILTER.out.reads)
    
    // LOG_COMPLEX_READS(KMERNORM_v1_0_0.out.plex_countfile.collect())
	// LOG_NORMALIZED_READS(KMERNORM_v1_0_0.out.norm_countfile.collect())
	DEINTERLEAVE(KMERNORM_v1_0_0.out.reads)

	// Build the BWA index once, reusing an existing one next to contam_ref_fasta if present.
	def bwa_fasta_file = file(params.contam_ref_fasta, checkIfExists: true)
	def existing_index = ['amb','ann','bwt','pac','sa'].collect { file("${params.contam_ref_fasta}.${it}") }
	if (existing_index.every { it.exists() }) {
		CH_bwa_fasta = channel.value(bwa_fasta_file)
		CH_bwa_index = channel.value(existing_index)
	} else {
		BWA_INDEX(channel.value(bwa_fasta_file))
		CH_bwa_fasta = BWA_INDEX.out.fasta
		CH_bwa_index = BWA_INDEX.out.index
	}

	CONTAM_READ_FINDER(DEINTERLEAVE.out.reads_gz, CH_bwa_fasta, CH_bwa_index)
	CONTAM_READ_REPORTER(CONTAM_READ_FINDER.out.aligned)
	CONTAM_READ_REMOVER(CONTAM_READ_REPORTER.out.join(KMERNORM_v1_0_0.out.reads))
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
    //container 'complexity-filter-env:latest'
    container 'brwnj/kmernorm:v1.0.0'
    publishDir { "${params.output}/${ID}/reads_${ID}" }, pattern: "*fastq.gz", mode: params.publishmode
    input: tuple val(ID), path(r1), path(r2)
    output: tuple val(ID), path("pe_${ID}.fastq.gz"), emit: reads
    script: template 'complexity_filter.py' }

process KMERNORM_v1_0_0 {
    tag "${ID}"
    //memory = { 16.GB * task.attempt }
    //errorStrategy = {task.attempt <= 3 ? 'retry' : 'ignore'}
    //maxRetries = 3
    container 'brwnj/kmernorm:v1.0.0' // 'brwnj/kmernorm:v1.1.0'
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

process DEINTERLEAVE {
	tag "${ID}"
	container 'quay.io/biocontainers/bbmap:38.90--he522d1c_3'
    publishDir { "${params.output}/${ID}/reads_${ID}"}, pattern: "r*_norm*.fastq.gz", mode: params.publishmode
    input: tuple val(ID), path(normed)
	output: tuple val(ID), path("r1_norm_${ID}.fastq.gz"), path("r2_norm_${ID}.fastq.gz"), emit: reads_gz
	script:
	"""
	reformat.sh in=${normed} out1=r1_norm_${ID}.fastq out2=r2_norm_${ID}.fastq
	gzip r1_norm_${ID}.fastq
	gzip r2_norm_${ID}.fastq
	""" }

process BWA_INDEX {
    tag "BWA index"
    container 'quay.io/biocontainers/bwa:0.7.17--h5bf99c6_8'
    // Persist the index next to the reference itself so future runs (and other
    // projects pointed at the same contam_ref_fasta) find it via the auto-detect check
    // in the workflow block instead of rebuilding it.
    publishDir { file(params.contam_ref_fasta).getParent() }, pattern: "*.{amb,ann,bwt,pac,sa}", mode: 'copy'
    input:
    path fasta

    output:
    path("${fasta}.*"), emit: index
    path(fasta),        emit: fasta

    script:
    """
    bwa index $fasta
    """
}

process CONTAM_READ_FINDER {
    tag "${ID}"
    memory '9.GB'
    // Docker Desktop's VM is capped at ~11.9GB total; this process is the only
    // thing that comes close to that ceiling, so cap concurrency at 1 to make
    // sure two samples can't stack their memory demand and OOM the VM even
    // though each individually fits under it.
    maxForks 1
    errorStrategy 'terminate'
    container 'quay.io/biocontainers/bwa:0.7.17--h5bf99c6_8'
    input:
        tuple val(ID), path(norm1), path(norm2)
        path(fasta)
        path(index)
    output: tuple val(ID), path("norm_${ID}_r1.fastq.sai"), path("norm_${ID}_r2.fastq.sai"), path("contam_${ID}.sam"), emit: aligned
    shell:
    '''
    set -e
    echo 'aligning forward reads to reference contaminant'
    bwa aln -n !{params.reference_threshold} -t !{task.cpus} !{fasta} !{norm1} > norm_!{ID}_r1.fastq.sai
    echo 'aligning reverse reads to reference contaminant'
    bwa aln -n !{params.reference_threshold} -t !{task.cpus} !{fasta} !{norm2} > norm_!{ID}_r2.fastq.sai
    echo 'retrieving hits'
    bwa sampe !{fasta} norm_!{ID}_r1.fastq.sai norm_!{ID}_r2.fastq.sai !{norm1} !{norm2} > contam_!{ID}.sam
    '''  }

process CONTAM_READ_REPORTER {
    tag "${ID}"
    container 'quay.io/biocontainers/samtools:1.24--h9dcdb79_1'
    publishDir { "${params.output}/${ID}/reads_${ID}" }, pattern: "contam_align_*.tsv", mode: params.publishmode
    input: tuple val(ID), path(sai1), path(sai2), path(sam)
    output: tuple val(ID), path("contam_align_${ID}.tsv")
    shell:
    '''
    set -e
    grep -i "^@" !{sam} > headers_only.tmp
    if [ "$(wc -l < headers_only.tmp)" -eq "$(wc -l < !{sam})" ];
    then
        echo 'no contaminant hits found'
        touch contam_align_!{ID}.tsv
    else
        echo 'contaminants were found'
        samtools view -SF0x0004 !{sam} > contam_align_!{ID}.tsv
    fi
    '''  }

process CONTAM_READ_REMOVER {
    tag "${ID}"
    errorStrategy { task.exitStatus in [0] ? 'ignore' : 'retry' }
    container 'brwnj/kmernorm:v1.0.0'
    publishDir { "${params.output}/${ID}/reads_${ID}" }, pattern: "contamfiltered_pe_*.fastq.gz", mode: 'copy'
    input: tuple val(ID), path(sam_contam), path(norm)
    output:
        tuple val(ID), path("contamfiltered_pe_${ID}.fastq.gz"), emit: reads
        path("5_contamfiltered_pe_${ID}.count"), emit: countfile
    script: template "contam_read_remover.py" }