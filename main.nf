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

//# CONTIG PROCESSING
params.assembly_minlength = '1000'
params.assembly_righttrim = '200'
params.assembly_lefttrim = '200'

//#     DECONTAMINATION
// Relative + gitignored: a fresh clone has nothing here, and BWA_INDEX/BLAST_INDEX
// download the fasta + prebuilt indexes from Zenodo (DOI 10.5281/zenodo.21682938) into
// this path automatically on first run — no manual setup needed. Override with
// --contam_ref_fasta to point at an already-populated location instead (e.g. a shared
// path on a cluster) and avoid re-downloading the ~8GB reference per clone.
params.contam_ref_fasta = "./reference/GRCh38_AG665_mm10.fa"
params.contam_min_length=100 // for BLASTn on contigs
params.contam_min_percid=95.0 // for BLASTn on contigs
// BWA/BLAST indexes are auto-detected next to contam_ref_fasta; the download only runs
// when they (or the fasta itself) are missing.

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

    def CH_stepwise_counts = "${params.output}/sample_tracking/stepwise_counts"
    def COUNT_HEADER = "Metric,Count,Sample_ID\n"

    FASTQC_v0_11_9(CH_fastq)
    CH_count_raw_reads = FASTQC_v0_11_9.out.countfile.collectFile(name: '1_raw_readcounts.csv', storeDir: CH_stepwise_counts, seed: COUNT_HEADER, cache: false)

    TRIMMOMATIC_v0_32(CH_fastq)
    CH_count_trimmed_reads = TRIMMOMATIC_v0_32.out.countfile.collectFile(name: '2_trimmed_readcounts.csv', storeDir: CH_stepwise_counts, seed: COUNT_HEADER, cache: false)
    
    COMPLEXITY_FILTER(TRIMMOMATIC_v0_32.out.reads)

    KMERNORM_v1_0_0(COMPLEXITY_FILTER.out.reads)
    CH_count_complex_reads = KMERNORM_v1_0_0.out.plex_countfile.collectFile(name: '3_complex_readcounts.csv', storeDir: CH_stepwise_counts, seed: COUNT_HEADER, cache: false)
    CH_count_normalized_reads = KMERNORM_v1_0_0.out.norm_countfile.collectFile(name: '4_normalized_readcounts.csv', storeDir: CH_stepwise_counts, seed: COUNT_HEADER, cache: false)

	DEINTERLEAVE(KMERNORM_v1_0_0.out.reads)

	// Fetch the reference fasta + its BWA index from Zenodo once, reusing them if both
	// already exist next to contam_ref_fasta. The fasta itself is no longer assumed to
	// pre-exist (BWA_INDEX downloads it too), so this can't use checkIfExists: true.
	def bwa_fasta_file = file(params.contam_ref_fasta)
	def existing_index = ['amb','ann','bwt','pac','sa'].collect { file("${params.contam_ref_fasta}.${it}") }
	if (bwa_fasta_file.exists() && existing_index.every { it.exists() }) {
		CH_bwa_fasta = channel.value(bwa_fasta_file)
		CH_bwa_index = channel.value(existing_index)
	} else {
		log.info "Reference fasta + BWA index not found at ${params.contam_ref_fasta} — downloading from Zenodo (~8GB, DOI 10.5281/zenodo.21682938). This is a one-time cost; future runs will reuse the downloaded files automatically."
		BWA_INDEX()
		CH_bwa_fasta = BWA_INDEX.out.fasta
		CH_bwa_index = BWA_INDEX.out.index
	}

	// Build the BLAST db once, reusing an existing one next to contam_ref_fasta if present.
	// Triggered here (early) rather than next to CONTAM_CONTIG_FINDER so it has the whole
	// read-processing + assembly pipeline to finish in before it's actually needed.
	// A large reference gets split into volumes by makeblastdb, so only a top-level .nal
	// alias file is written next to contam_ref_fasta — check that its listed volumes
	// actually exist too; an alias file by itself doesn't mean the db is usable.
	def nal_file = file("${params.contam_ref_fasta}.nal")
	def blast_db_ready
	if (nal_file.exists()) {
		def dblist = (nal_file.text =~ /(?m)^DBLIST\s+(.+)$/)
		def volumes = dblist ? dblist[0][1].replaceAll('"', '').trim().split(/\s+/) : []
		// DBLIST entries are bare filenames (no directory) — resolve them against
		// contam_ref_fasta's own directory, not wherever `nextflow run` was launched
		// from, or this always evaluates false and BLAST_INDEX reruns every time.
		def ref_dir = nal_file.getParent()
		blast_db_ready = volumes.size() > 0 && volumes.every { file("${ref_dir}/${it}.nhr").exists() }
	} else {
		blast_db_ready = file("${params.contam_ref_fasta}.nhr").exists()
	}
	if (blast_db_ready) {
		CH_blast_fasta = CH_bwa_fasta
		// Glob covers both single-volume (<fasta>.nhr) and multi-volume (<fasta>.nal,
		// <fasta>.00.nhr, ...) layouts without hand-maintaining BLAST's version-dependent
		// extension list, the way the BWA branch above can with its fixed 5 extensions.
		CH_blast_index = channel.value(file("${params.contam_ref_fasta}*.n??"))
	} else {
		// CH_bwa_fasta (not the raw bwa_fasta_file variable) so this properly waits on
		// BWA_INDEX's download finishing when the fasta doesn't exist yet, instead of
		// just hoping it's ready in time from being placed earlier in the script.
		log.info "BLAST db not found next to ${params.contam_ref_fasta} — building it now. This is a one-time cost; future runs will reuse it automatically."
		BLAST_INDEX(CH_bwa_fasta)
		CH_blast_fasta = BLAST_INDEX.out.fasta
		CH_blast_index = BLAST_INDEX.out.index
	}

	CONTAM_READ_FINDER(DEINTERLEAVE.out.reads_gz, CH_bwa_fasta, CH_bwa_index)
	CONTAM_READ_REPORTER(CONTAM_READ_FINDER.out.aligned)
	CONTAM_READ_REMOVER(CONTAM_READ_REPORTER.out.join(KMERNORM_v1_0_0.out.reads))
	CH_count_clean_reads = CONTAM_READ_REMOVER.out.countfile.collectFile(name: '5_clean_readcounts.csv', storeDir: CH_stepwise_counts, seed: COUNT_HEADER, cache: false)

    SPADES_v3_15_2(CONTAM_READ_REMOVER.out.reads)
    CH_count_raw_contigs =SPADES_v3_15_2.out.countfile.collectFile(name: '6_raw_contigcounts.csv', storeDir: CH_stepwise_counts, seed: COUNT_HEADER, cache: false)

    TRIM_CONTIGS(SPADES_v3_15_2.out.passed)
    CH_count_trimmed_contigs = TRIM_CONTIGS.out.countfile.collectFile(name: '7_trimmed_contigcounts.csv', storeDir: CH_stepwise_counts, seed: COUNT_HEADER, cache: false)

    CONTAM_CONTIG_FINDER(TRIM_CONTIGS.out.trimmed_contigs, CH_blast_fasta, CH_blast_index)
    CONTAM_CONTIG_REMOVER(CONTAM_CONTIG_FINDER.out.join(TRIM_CONTIGS.out.length_passing_contigs))
    CH_count_final_contigs = CONTAM_CONTIG_REMOVER.out.countfile.collectFile(name: '8_final_contigcounts.csv', storeDir: CH_stepwise_counts, seed: COUNT_HEADER, cache: false)

    MEASURE_SAG(TRIM_CONTIGS.out.trimmed_contigs)
    CH_count_sag_stats = MEASURE_SAG.out.countfile.collectFile(name: '8_sag_stats.csv', storeDir: CH_stepwise_counts, seed: COUNT_HEADER, cache: false)

    CHECKM_v1_1_9(TRIM_CONTIGS.out.trimmed_contigs)
    // .map() pulls the 'Completeness' column out of CheckM's --tab_table TSV directly in
    // Groovy (no separate PARSE_CHECKM process/container needed for three lines of pandas).
    CH_count_checkm = CHECKM_v1_1_9.out
        .map { ID, checkm_dir ->
            def lines = file("${checkm_dir}/completeness_${ID}.tsv").readLines()
            def header = lines[0].split('\t')
            def completeness = lines[1].split('\t')[header.findIndexOf { it == 'Completeness' }]
            "CheckM1_est_genome_completeness,${completeness},${ID}\n"
        }
        .collectFile(name: '9_checkm_completeness.csv', storeDir: CH_stepwise_counts, seed: COUNT_HEADER, cache: false, sort: false)

    // Combine every stepwise count file into one, earliest stage first. Each one already
    // carries its own COUNT_HEADER line (from its own collectFile seed above) — strip
    // that off per file before stacking, then let this collectFile's own seed add it back at the end.
    CH_all_counts = CH_count_raw_reads
        .concat(CH_count_trimmed_reads, CH_count_complex_reads, CH_count_normalized_reads,
                CH_count_clean_reads, CH_count_raw_contigs, CH_count_trimmed_contigs,
                CH_count_final_contigs, CH_count_sag_stats, CH_count_checkm)
        .map { it.text.readLines().drop(1).join('\n') + '\n' }
        .collectFile(name: 'all_stepwise_counts.csv', storeDir: "${params.output}/sample_tracking", seed: COUNT_HEADER, cache: false, sort: false)

    ASSEMBLY_STATS_TABULATOR(CH_all_counts)
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
        echo "Raw_readcount,$(expr $(zcat !{r1} | wc -l) / 4 + $(zcat !{r2} | wc -l) / 4),!{ID}" >> "1_raw_!{ID}.count"
        fastqc -t !{task.cpus} -q !{r1} !{r2} --noextract
        '''
    stub:
    // Neither .qc nor .html is consumed downstream (only .countfile is), so these just
    // need to exist to satisfy the declared outputs.
    """
    touch stub_fastqc.zip stub_fastqc.html
    echo "Raw_readcount,20,${ID}" >> "1_raw_${ID}.count"
    """ }

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
        echo "Trimmed_readcount,$(expr $(zcat trimmed_!{ID}_r1.fastq.gz | wc -l) / 4 + $(zcat trimmed_!{ID}_r2.fastq.gz | wc -l) / 4),!{ID}" > "2_trimmed_!{ID}.count"
        '''
    stub:
    // First 10 real reads from each of r1/r2 — becomes COMPLEXITY_FILTER's stub input,
    // which genuinely zcats/subsets it, so this needs to be valid gzipped fastq, not empty.
    """
    zcat ${r1} | head -40 | gzip > trimmed_${ID}_r1.fastq.gz
    zcat ${r2} | head -40 | gzip > trimmed_${ID}_r2.fastq.gz
    echo "Trimmed_readcount,20,${ID}" > 2_trimmed_${ID}.count
    """ }

process COMPLEXITY_FILTER {
    tag "${ID}"
    //container 'complexity-filter-env:latest'
    container 'brwnj/kmernorm:v1.0.0'
    publishDir { "${params.output}/${ID}/reads_${ID}" }, pattern: "*fastq.gz", mode: params.publishmode
    input: tuple val(ID), path(r1), path(r2)
    output: tuple val(ID), path("pe_${ID}.fastq.gz"), emit: reads
    script: template 'complexity_filter.py'
    stub:
    // First 10 real read pairs from r1/r2, properly interleaved: `paste - - - -`
    // collapses each 4-line fastq record to one line, pasting the two collapsed
    // streams side by side then expanding tabs back to newlines interleaves them.
    """
    paste <(zcat ${r1} | head -40 | paste - - - -) <(zcat ${r2} | head -40 | paste - - - -) | tr '\\t' '\\n' | gzip > pe_${ID}.fastq.gz
    """ }

process KMERNORM_v1_0_0 {
    tag "${ID}"
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
        echo "Complexity_filtered_readcount,$(expr $(cat temp_paired.fastq | wc -l) / 4),!{ID}" >> "3_pe_!{ID}.count"
        kmernorm !{params.kmernorm_opts} temp_paired.fastq > normalized_pe_!{ID}.fastq
        echo "Normalized_readcount,$(expr $(cat normalized_pe_!{ID}.fastq | wc -l) / 4),!{ID}" >> "4_normalized_pe_!{ID}.count"
        gzip normalized_pe_!{ID}.fastq

        # Cleanup
        rm temp_paired.fastq
        '''
    stub:
    // `paired` (from COMPLEXITY_FILTER's stub) is already real subsetted read data —
    // just carry it forward under the expected output name rather than re-deriving it.
    """
    cp ${paired} normalized_pe_${ID}.fastq.gz
    echo "Complexity_filtered_readcount,20,${ID}" >> 3_pe_${ID}.count
    echo "Normalized_readcount,20,${ID}" >> 4_normalized_pe_${ID}.count
    """ }

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
    tag "downloading ref contaminant + BWA index from Zenodo"
    container 'curlimages/curl:8.21.0'
    // curlimages/curl deliberately runs as a non-root user by default. Docker Desktop's
    // bind-mount layer on macOS is lenient about that mismatch against the host-owned
    // work dir, but real Linux Docker (GitHub Actions runners, HPC nodes) enforces it for
    // real and this container cannot write its own output without -u root.
    containerOptions '--entrypoint "" -u root'
    shell '/bin/sh', '-ue'
    // Persist both next to contam_ref_fasta so future runs (and other projects pointed
    // at the same path) find them via the auto-detect check instead of re-downloading.
    // enabled: !workflow.stubRun keeps stub output from landing in this shared directory.
    publishDir { file(params.contam_ref_fasta).getParent() }, mode: 'copy', enabled: !workflow.stubRun
    cache false
    output:
    path("${file(params.contam_ref_fasta).getName()}"), emit: fasta
    path("${file(params.contam_ref_fasta).getName()}.{amb,ann,bwt,pac,sa}"), emit: index

    script:
    def base = file(params.contam_ref_fasta).getName()
    // DOI 10.5281/zenodo.21682938 — "GORG Dark - Reference Contaminant Dataset": the
    // exact GRCh38_AG665_mm10.fa + prebuilt BWA index this pipeline already expects.
    // All 6 files are zip archives, each wrapping the real, already-decompressed file
    // under its real name (e.g. fetching "<base>.ann.zip" returns a zip whose sole
    // entry is literally "<base>.ann"). MD5s below are Zenodo's published checksums of
    // the .zip files themselves, checked before extracting so a corrupt/partial
    // download fails loudly instead of silently producing a broken reference.
    """
    set -e
    fetch() {
        curl -sL -o "\$1" "https://zenodo.org/api/records/21682938/files/\$1/content"
        echo "\$2  \$1" | md5sum -c -
        unzip -oq "\$1"
        rm "\$1"
    }
    fetch "${base}.zip"     262152088a82f34416c9e30e142cab9d
    fetch "${base}.amb.zip" d22f375cd7f16768f99b22c005b1b452
    fetch "${base}.ann.zip" 66da9fa8ab74c58e87f57cfa5dc087a2
    fetch "${base}.bwt.zip" cb0be188a3dfb8bac70b7d8f16dd61e5
    fetch "${base}.pac.zip" 66722c0f8f463d790e6cdb99f960e99f
    fetch "${base}.sa.zip"  c16b36eea4b1f75761bda8ed3e62bea5
    """
    stub:
    // `base` above is local to the script: closure, not visible here — recompute inline
    // (same reason output: does above rather than referencing a shared variable).
    // Explicit filenames, not brace expansion ({,.amb,...}) — that's a bash-ism, and
    // this process's shell is /bin/sh (BusyBox ash, no /bin/bash in this container).
    // Under sh, the unexpanded brace expression became one literal filename instead of
    // six, so `touch` exited 0 while never actually creating the expected output.
    """
    touch ${file(params.contam_ref_fasta).getName()} \
          ${file(params.contam_ref_fasta).getName()}.amb \
          ${file(params.contam_ref_fasta).getName()}.ann \
          ${file(params.contam_ref_fasta).getName()}.bwt \
          ${file(params.contam_ref_fasta).getName()}.pac \
          ${file(params.contam_ref_fasta).getName()}.sa
    """
}

process BLAST_INDEX {
    tag "makeblastdb on ref contaminants"
    container 'quay.io/biocontainers/blast:2.11.0--pl5262h3289130_1'
    // Persist the db next to the reference itself, same convention as BWA_INDEX, so
    // future runs find it via the auto-detect check instead of rebuilding it. Same
    // enabled: !workflow.stubRun guard as BWA_INDEX — this is exactly what leaked
    // 0-byte stub .nhr/.nin/.nsq files into the real reference directory before.
    publishDir { file(params.contam_ref_fasta).getParent() }, pattern: "${file(params.contam_ref_fasta).getName()}.*", mode: 'copy', enabled: !workflow.stubRun
    // Same reasoning as BWA_INDEX: the outer workflow-level check (not Nextflow's task
    // cache) is what avoids redundant rebuilds, so disabling caching here costs nothing
    // and closes off stub-run/real-run cache cross-contamination.
    cache false
    input:
    path fasta

    output:
    path("${fasta}.*"), emit: index
    path(fasta),        emit: fasta

    script:
    """
    makeblastdb -in $fasta -dbtype nucl -out $fasta -title $fasta
    """
    stub:
    """
    touch ${fasta}.nhr ${fasta}.nin ${fasta}.nsq
    """
}

process CONTAM_READ_FINDER {
    tag "${ID}"
    // Retries with more memory instead of a fixed ceiling, so this adapts to whatever's
    // actually available (a laptop, CI, an HPC node) rather than encoding one machine's
    // Docker Desktop allocation. 5.GB is where this genuinely OOM'd once for real against
    // this reference; attempt 2 (10GB) lands right around where it's since run reliably.
    memory { 5.GB * task.attempt }
    maxForks 1 // keeps memory-heavy retries from stacking across samples regardless of environment
    errorStrategy { task.exitStatus in [137, 140] ? 'retry' : 'terminate' }
    maxRetries 3
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
    '''
    stub:
    """
    touch norm_${ID}_r1.fastq.sai norm_${ID}_r2.fastq.sai
    printf '@HD\\tVN:1.6\\tSO:unsorted\\n' > contam_${ID}.sam
    """ }

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
    container 'brwnj/kmernorm:v1.0.0'
    publishDir { "${params.output}/${ID}/reads_${ID}" }, pattern: "contamfiltered_pe_*.fastq.gz", mode: 'copy'
    input: tuple val(ID), path(sam_contam), path(norm)
    output:
        tuple val(ID), path("contamfiltered_pe_${ID}.fastq.gz"), emit: reads
        path("5_contamfiltered_pe_${ID}.count"), emit: countfile
    script: template "contam_read_remover.py" }

process SPADES_v3_15_2 {
    tag "${ID} assembling"
    errorStrategy 'ignore'
    container 'quay.io/biocontainers/spades:3.15.2--h95f258a_1'
    publishDir { "${params.output}/${ID}/intermediate_assemblies_${ID}" }, pattern: "0_all_contigs_*.fasta", mode: params.publishmode
    input: tuple val(ID), path(clean_reads)
    output:
        tuple val(ID), path("0_all_contigs_${ID}.fasta"), emit: passed
        path("6_all_contigs_${ID}.count"), emit: countfile
    shell:
        '''
        [ ! -d spades_!{ID} ] && mkdir spades_!{ID}
        spades.py -o spades_!{ID} --careful --sc --phred-offset 33 -t !{task.cpus} --12 !{clean_reads}
        cp spades_!{ID}/contigs.fasta 0_all_contigs_!{ID}.fasta
        echo "Raw_contig_count,$(grep -c '>' spades_!{ID}/contigs.fasta),!{ID}" > 6_all_contigs_!{ID}.count
        '''
    stub:
    """
    printf '>stub_contig_1\\n' > 0_all_contigs_${ID}.fasta
    yes ACGT | head -400 | tr -d '\\n' >> 0_all_contigs_${ID}.fasta
    printf '\\n' >> 0_all_contigs_${ID}.fasta
    echo "Raw_contig_count,1,${ID}" > 6_all_contigs_${ID}.count
    """ }

process TRIM_CONTIGS {
    errorStrategy 'finish'
    tag "${ID}"
    container 'brwnj/kmernorm:v1.0.0'
    publishDir { "${params.output}/${ID}/intermediate_assemblies_${ID}" }, pattern: "*.fasta", mode: params.publishmode
    input: tuple val(ID), path(contigs)
    output:
        tuple val(ID), path("2a_long_passing_contigs_${ID}.fasta"), emit: length_passing_contigs
        tuple val(ID), path("2b_short_discarded_contigs_${ID}.fasta"), emit: short_contigs
        tuple val(ID), path("3_trimmed_contigs_${ID}.fasta"), emit: trimmed_contigs
        path("2a_long_passing_contigs_${ID}.fasta"), emit: fasta
        path("7_length_passing_contigs_${ID}.count"), emit: countfile
    script: template 'trim_and_deduplicate_contigs.py' }

process MEASURE_SAG{
    tag "${ID}"
    container 'brwnj/kmernorm:v1.0.0'
    publishDir { "${params.output}/${ID}/logs_${ID}" }
    input: tuple val(ID), path(contigs)
    output: path("sag_stats_${ID}.csv"), emit: countfile
    script:
    """
    #!/usr/bin/env python
    from Bio import SeqIO; from Bio.SeqUtils import GC
    max_contig_length = 0; STR_all_seq = ''
    for record in SeqIO.parse("${contigs}", "fasta"):
        STR_all_seq = STR_all_seq + record.seq          # read in the DNA, 1 contig at a time
        if len(record.seq) > max_contig_length:
            max_contig_length = len(record.seq)
    out=open("sag_stats_${ID}.csv", "a")
    print("Max_contig_length", str(max_contig_length), "${ID}", sep=",", file=out)
    print("Final_assembly_length", str(len(STR_all_seq)), "${ID}", sep=",", file=out)
    print("GC_content", str(round(GC(STR_all_seq),2)), "${ID}", sep=",", file=out)
    out.close()
    """  }

process CONTAM_CONTIG_FINDER {
    tag "${ID}"
    container 'quay.io/biocontainers/blast:2.11.0--pl5262h3289130_1'
    input:
        tuple val(ID), path(contigs)
        path(refcontam_fasta)
        path(index)
    output: tuple val(ID), path("contig_blast_${ID}.tsv")
    shell:
        '''
        blastn -db !{refcontam_fasta} -query !{contigs} -out tmp_1_blast.tsv -num_threads !{task.cpus} -max_target_seqs 10 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sallseqid score nident positive gaps ppos qframe sframe qseq sseq qlen slen salltitles"
        echo -e 'Query Seq-id\tSubject Seq-id\tPercentage of identical matches\tAlignment length\tNumber of mismatches\tNumber of gap openings\tStart of alignment in query\tEnd of alignment in query\tStart of alignment in subject\tEnd of alignment in subject\tExpect value\tBit score\tAll subject Seq-id(s)\tRaw score\tNumber of identical matches\tNumber of positive-scoring matches\tTotal number of gaps\tPercentage of positive-scoring matches\tQuery frame\tSubject frame\tAligned part of query sequence\tAligned part of subject sequence\tQuery sequence length\tSubject sequence length\tAll Subject Title(s)' > tmp_header.tsv
        cat tmp_header.tsv tmp_1_blast.tsv > contig_blast_!{ID}.tsv
        rm tmp_header.tsv tmp_1_blast.tsv
        '''
    stub:
    // Header-only, matching what a real run with zero blast hits produces — keeps
    // CONTAM_CONTIG_REMOVER's real (unstubbed) parser downstream working normally.
    """
    echo -e 'Query Seq-id\tSubject Seq-id\tPercentage of identical matches\tAlignment length\tNumber of mismatches\tNumber of gap openings\tStart of alignment in query\tEnd of alignment in query\tStart of alignment in subject\tEnd of alignment in subject\tExpect value\tBit score\tAll subject Seq-id(s)\tRaw score\tNumber of identical matches\tNumber of positive-scoring matches\tTotal number of gaps\tPercentage of positive-scoring matches\tQuery frame\tSubject frame\tAligned part of query sequence\tAligned part of subject sequence\tQuery sequence length\tSubject sequence length\tAll Subject Title(s)' > contig_blast_${ID}.tsv
    """ }

process CONTAM_CONTIG_REMOVER {
    tag "${ID}"
    container 'brwnj/kmernorm:v1.0.0'
    publishDir { "${params.output}/${ID}/intermediate_assemblies_${ID}" }, pattern: "4b_contam_contigs_*.fasta"
    publishDir { "${params.output}/${ID}" }, pattern: "SCGC_*_contigs.fasta" 
    input:
        tuple val(ID), path(blast_tsv), path(contigs)
    output:
        tuple val(ID), path("SCGC_${ID}_contigs.fasta"), emit: contigs
        path("SCGC_${ID}_contigs.fasta"), emit: fasta
        path("4b_contam_contigs_${ID}.fasta"), emit: contam_contigs
        path("8_final_contigs_${ID}.count"), emit: countfile
    script: template 'contam_contig_remover.py' }

process CHECKM_v1_1_9 {
    tag "${ID} estimate completeness"
    // --reduced_tree is a memory-saver for single-genome input, but it can be inaccurate for some lineages.
    // If memory is NOT limiting, run CHECKM without --reduced_tree.
    // Retries with more memory instead of a fixed ceiling — see CONTAM_READ_FINDER for
    // the same pattern and why. --reduced_tree's own guarantee is only "<16GB", so this
    // covers that range across two attempts rather than assuming a single number.
    memory { 6.GB * task.attempt }
    maxForks 1 // keeps memory-heavy retries from stacking across samples regardless of environment
    errorStrategy { task.exitStatus in [137, 140] ? 'retry' : 'terminate' }
    maxRetries 3
    container 'quay.io/biocontainers/checkm-genome:1.1.9--pyhdfd78af_0'
    publishDir { "${params.output}/${ID}/QC_${ID}" }, mode: params.publishmode
    input: tuple val(ID), path(contigs)
    output: tuple val(ID), path("checkm_${ID}")
    script:
    """
    mkdir tmp_dir; cp ${contigs} ./tmp_dir/final_contigs_${ID}.fasta
    checkm lineage_wf --reduced_tree -f checkm_${ID}/completeness_${ID}.tsv --tab_table -q -x fasta -t ${task.cpus} tmp_dir checkm_${ID}
    rm -r tmp_dir
    """
    stub:
    """
    mkdir checkm_${ID}
    printf "Bin Id\\tCompleteness\\tContamination\\n" > checkm_${ID}/completeness_${ID}.tsv
    printf "final_contigs_${ID}\\t99.9\\t0.1\\n" >> checkm_${ID}/completeness_${ID}.tsv
    """ }

process ASSEMBLY_STATS_TABULATOR {
    container 'brwnj/kmernorm:v1.0.0'
    publishDir "${params.output}"
    input: path(COUNTS_TXT)
    output: path("assembly_stats.csv")
    script:
    // Applied to every cell below so stub-run output can never be mistaken for a real
    // assembly's stats — empty string on a real run, no visible effect.
    def STUB_PREFIX = workflow.stubRun ? "STUB-RUN " : ""
    """
    #!/usr/bin/env python
    import pandas as pd
    import numpy as np
    PATH_out = "assembly_stats.csv"
    DF_log = pd.read_csv("${COUNTS_TXT}")

    LIST_col_order = ["Sample_ID", "Raw_readcount", "Trimmed_readcount", "Complexity_filtered_readcount", "Normalized_readcount", "Contam_filtered_readcount", "Raw_contig_count", "Final_clean_contig_count", "Max_contig_length", "Final_assembly_length", "GC_content", "CheckM1_est_genome_completeness"]

    print("Converting list of read/contig counts to table...")
    try:
        DF_log = DF_log.pivot(index="Sample_ID", columns="Metric", values="Count") # Pivot to table indexed by Sample_ID
    except:
        ### Deal with the rare edge case where some samples have redundant analyses (e.g. Nextflow spawned the same job twice)
        DF_log_orig = DF_log
        # Use "aggfunc=first" to discard redundant ID+metrics
        DF_log = DF_log.pivot_table(index="Sample_ID", columns="Metric", values="Count", aggfunc='first')
    
    DF_log.reset_index(inplace=True) # Make index 'Sample_ID' -> column

    # If no dirty reads were found (i.e. Contam_filtered_readcount says "NO_CHANGE"), make numeric by copying readcount from upstream.
    DF_log['Contam_filtered_readcount'] = np.where(DF_log['Contam_filtered_readcount']=='NO_CHANGE',DF_log['Normalized_readcount'],DF_log['Contam_filtered_readcount']) 

    DF_log = DF_log[LIST_col_order] # Reorder columns to match LIST_col_order
    ## Warn the user when the metrics are from a stub run, so they don't mistake them for real data.
    DF_log = DF_log.applymap(lambda x: "${STUB_PREFIX}" + str(x))
    DF_log.to_csv(PATH_out, index=False)
    """

}