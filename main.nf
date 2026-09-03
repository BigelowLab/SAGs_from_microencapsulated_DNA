#!/usr/bin/env nextflow
//
// SAGs_from_microcapsules: a Nextflow pipeline for single-cell genome assembly and annotation from Atrandi combinatorial barcoded Illumina reads
// Copyright (C) 2026  Greg Gavelis
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <https://www.gnu.org/licenses/>.
//
nextflow.enable.dsl=2

params.indir ="./input/"
params.output ="./results/"

// MODES
params.dev = false
params.dev_num_capsules = 10 // how many capsules --dev carries into assembly
params.viral = true
params.target_IDs = "" // comma-separated capsule IDs to restrict downstream processing to (e.g. "4_AACCGGTT,7_TTGGCCAA"); empty runs every capsule

// Defaults
params.publishmode = 'symlink'

//# ATRANDI DEMULTIPLEXING
params.barcode_dir = "${projectDir}/barcodes"
params.sample_random_seed = 42
params.sample_num_reads = 1000000
params.sample_hamming_dist = 1
params.split_hamming_dist = 1
params.read_threshold = 3
params.cell_threshold = 100000000000 // effectively uncapped
params.barcode_trim_length = 45 // TRIM_BARCODE's trim_galore --clip_R2
// Pheniqs' default is 2048 buffered records PER feed (input or output), so with one output
// feed pair per capsule this scales directly with capsule count — with ~1900 real capsules in
// one pool that's ~2000x the default's memory footprint. Turned down here since PHENIQS_DEMULTIPLEX
// OOM'd against real data at the default; raise it again if throughput becomes the bottleneck instead.
params.pheniqs_buffer_capacity = 64

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
params.contam_ref_fasta = "./reference/GRCh38_AG665_mm10.fa"
params.contam_min_length=100 // for BLASTn on contigs
params.contam_min_percid=95.0 // for BLASTn on contigs
// BWA/BLAST indexes are auto-detected next to contam_ref_fasta; the download only runs when they (or the fasta itself) are missing.

//#    VIRAL
params.DB_genomad_v1_11_1 = "/mnt/scgc_nfs/ref/genomad/genomad_1.11.1/genomad_db/"
params.prokka = "/mnt/scgc_nfs/ref/uniprot_swissprot_prokka.fasta"
params.PATH_hmm = "/mnt/databases/scgc/EggNOGdb/nog.hmm"
params.PATH_annot = "/mnt/databases/scgc/EggNOGdb/nog_annotation_virupdated.tsv"

//#    SSU (16S) RECOVERY + TAXONOMY
// Unlike the contaminant reference (auto-downloaded from Zenodo by BWA_INDEX), these
// databases are NOT fetched by the pipeline — supply your own and point these params at
// them (see README "Reference data"). Defaults are the paths used on the Bigelow cluster
// that produced the published results.
// TODO: publish the SILVA rRNA DB and the Prokka SwissProt DB to Zenodo and switch these
// to the same self-bootstrapping download pattern BWA_INDEX uses.
params.silva_blastdb = "/mnt/scgc_nfs/ref/silva_rrna/v128/silvamod128.fasta"
params.silva_map     = "/mnt/scgc_nfs/ref/silva_rrna/v128/silvamod128.map"
params.silva_tree    = "/mnt/scgc_nfs/ref/silva_rrna/v128/silvamod128.tre"
params.gtdb          = "/mnt/scgc_nfs/ref/gtdb/release207" // GTDB r207 — the release GTDB-Tk 2.0.0 expects; kept at 2.0.0/r207 for fidelity to the published results
params.gtdbtk_min_bp = 2500 // skip GTDB-Tk on assemblies smaller than this (total bases)

def maxContigLength(Path fasta) {
  fasta.splitFasta( record:[ seqString: true ] ) // Nextflow has a similarly named method for files which follows the same input as the channel operator
    *.seqString // Returns the values from the Map, i.e. the sequences
    *.size()    // Returns the size of each string
    .max()
}

def countBases(Path fasta) {
  fasta.splitFasta( record:[ seqString: true ] )
    *.seqString
    *.size()
    .sum()
}

workflow {

    // Check inputs
    def input_dir = file(params.indir)
    if (!input_dir.exists()) {
        error("Input directory not found: ${params.indir}.\nMake sure it exists and contains paired Illumina fastq files carrying Atrandi combinatorial barcodes, named like.\n  BLAH_R1.fastq.gz\n  BLAH_R2.fastq.gz")
    }
    // Each pair here is an Atrandi pool (many single-cell capsules multiplexed together via
    // combinatorial D/C/B/A barcodes in R2) — not yet SAG-ready. The ATRANDI DEMULTIPLEXING
    // block below splits each pool into one read pair per capsule before assembly begins.
    CH_library_fastq = channel.fromPath("${params.indir}/*.{fastq,fastq.gz,fq,fq.gz}", checkIfExists: true)
        .flatten() // emit each fastq path as its own item
        .map { file -> tuple(file.getSimpleName().replaceFirst('_R1','').replaceFirst('_R2',''), file) }    // Derive library by removing '_R1' or '_R2' suffixes. E.g. 4_12345678_R1.fastq.gz -> [4, 4_12345678_R1.fastq.gz]
        .groupTuple(size:2) // E.g. [4, [4_12345678_R1.fastq.gz, 4_12345678_R2.fastq.gz]]
        .map { library, files ->      // Pick r1/r2 out by filename rather than list position — groupTuple's
                                       // element order follows channel emission order, which isn't guaranteed
                                       // to be R1-then-R2 (unlike FASTQC/Trimmomatic, the barcode positions
                                       // this feeds into are order-sensitive enough that a silent swap here
                                       // would corrupt every downstream capsule ID).
            def r1 = files.find { it.getSimpleName().contains('_R1') }
            def r2 = files.find { it.getSimpleName().contains('_R2') }
            tuple(library, r1, r2)
        }

    def BC_D = file("${params.barcode_dir}/bcD_24.txt")
    def BC_C = file("${params.barcode_dir}/bcC_24.txt")
    def BC_B = file("${params.barcode_dir}/bcB_24.txt")
    def BC_A = file("${params.barcode_dir}/bcA_24.txt")

    SAMPLE_READS(CH_library_fastq)
    PHENIQS_MAKE_SAMPLE_CONFIG(SAMPLE_READS.out, BC_D, BC_C, BC_B, BC_A)
    PHENIQS_SAMPLE_DEMULTIPLEX(SAMPLE_READS.out.join(PHENIQS_MAKE_SAMPLE_CONFIG.out))
    PHENIQS_COUNT_SORT_SAMPLE(PHENIQS_SAMPLE_DEMULTIPLEX.out)
    PHENIQS_PLOT_HIST(PHENIQS_COUNT_SORT_SAMPLE.out)
    PHENIQS_FILTER_BC_LIST(PHENIQS_COUNT_SORT_SAMPLE.out)
    PHENIQS_NAME_CAPSULES(PHENIQS_FILTER_BC_LIST.out.join(CH_library_fastq), BC_D, BC_C, BC_B, BC_A)
    PHENIQS_DEMULTIPLEX(PHENIQS_NAME_CAPSULES.out.join(CH_library_fastq))
    PHENIQS_PARSE_REPORT(PHENIQS_DEMULTIPLEX.out.txt)

    // Each capsule fastq is named <library>_<capsuleID>_r{1,2}.fastq.gz by PHENIQS_DEMULTIPLEX.
    // Regroup into (ID, r1, r2) tuples — same shape the rest of the pipeline already expects —
    // dropping the undetermined bucket (reads whose barcode combo didn't survive filtering).
    CH_fastq = PHENIQS_DEMULTIPLEX.out.list_fastq.flatten()
        .map { file -> tuple(file.getSimpleName().replaceFirst(/_r1$/,'').replaceFirst(/_r2$/,''), file) }
        .groupTuple(size:2)
        .map { ID, files ->      // Pick r1/r2 out by filename, not list position — see the matching note
                                  // on CH_library_fastq above for why position isn't reliable here.
            def r1 = files.find { it.getSimpleName().endsWith('_r1') }
            def r2 = files.find { it.getSimpleName().endsWith('_r2') }
            tuple(ID, r1, r2)
        }
        .filter { !it[0].contains('undetermined') }

    // NOTE: unlike the pre-demux CH_fastq this replaces, --dev now caps the number of
    // *capsules* carried into assembly (params.dev_num_capsules), not the number of pools
    // demultiplexed — every pool still gets demuxed even in dev mode, since demux itself is
    // cheap relative to assembly.
    // Sorted by ID before taking, so --dev picks the same capsules on every run rather than
    // whatever order PHENIQS_DEMULTIPLEX's output glob happens to list them in (unspecified,
    // and not sorted by anything meaningful like read count either way). Only sorted when
    // --dev is actually on: toSortedList() has to collect the whole channel before re-emitting,
    // which would otherwise force every capsule through demux before any of them could start
    // assembly — fine for a small --dev subset, not something a full production run should pay.
    // --target_IDs restricts processing to an explicit, comma-separated capsule ID list
    // instead of --dev's numeric cap -- for rerunning/debugging specific known capsules
    // without waiting on the rest of the pool. Takes precedence over --dev when both are
    // set. Like --dev's take() below, this is a pure filter on an already-fully-demultiplexed
    // channel: it changes which capsules enter TRIM_BARCODE onward, not any per-capsule
    // task's own inputs, so it can't invalidate -resume caches for capsules it lets through.
    if (params.target_IDs) {
        def SET_target_IDs = params.target_IDs.tokenize(',')*.trim() as Set
        CH_fastq = CH_fastq.filter { ID, r1, r2 -> ID in SET_target_IDs }
    } else {
        CH_fastq = params.dev
            ? CH_fastq.toSortedList { a, b -> a[0] <=> b[0] }.flatMap { it }.take(params.dev_num_capsules)
            : CH_fastq
    }

    TRIM_BARCODE(CH_fastq)
    CH_fastq = TRIM_BARCODE.out

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
	def nal_file = file("${params.contam_ref_fasta}.nal")
	def blast_db_ready
	if (nal_file.exists()) {
		def dblist = (nal_file.text =~ /(?m)^DBLIST\s+(.+)$/)
		def volumes = dblist ? dblist[0][1].replaceAll('"', '').trim().split(/\s+/) : []
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
		CH_blast_index = channel.value(files("${params.contam_ref_fasta}*.n??"))
	} else {
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

    CONTAM_CONTIG_FINDER(TRIM_CONTIGS.out.trimmed_contigs.filter({ it[1].size()>0 }), CH_blast_fasta, CH_blast_index) // don't move forward with empty contigs
    // CONTAM_CONTIG_REMOVER excises contaminant regions from the *trimmed* contigs (not
    // length_passing_contigs), so its coordinates line up with the BLAST hits CONTAM_CONTIG_FINDER
    // produced from those same trimmed contigs, and SCGC_<ID>_contigs.fasta is the genuinely
    // trimmed + decontaminated assembly every downstream annotation step runs on.
    CONTAM_CONTIG_REMOVER(CONTAM_CONTIG_FINDER.out.join(TRIM_CONTIGS.out.trimmed_contigs.filter({ it[1].size()>0 })))
    CH_count_final_contigs = CONTAM_CONTIG_REMOVER.out.countfile.collectFile(name: '8_final_contigcounts.csv', storeDir: CH_stepwise_counts, seed: COUNT_HEADER, cache: false)

    // The single fully decontaminated + trimmed assembly per capsule. Everything below
    // (stats, annotation, SSU, GTDB-Tk, viral) consumes this rather than the pre-contig-decon
    // TRIM_CONTIGS output.
    CH_final_contigs = CONTAM_CONTIG_REMOVER.out.contigs.filter({ it[1].size()>0 })

    MEASURE_SAG(CH_final_contigs) //don't run on empty contigs file
    CH_count_sag_stats = MEASURE_SAG.out.countfile.collectFile(name: '8_sag_stats.csv', storeDir: CH_stepwise_counts, seed: COUNT_HEADER, cache: false)

    CHECKM_v1_1_9(CH_final_contigs) //don't run on empty contigs file
    // .map() pulls the 'Completeness' column out of CheckM's --tab_table TSV directly in Groovy
    CH_count_checkm = CHECKM_v1_1_9.out
        .map { ID, checkm_dir ->
            def lines = file("${checkm_dir}/completeness_${ID}.tsv").readLines()
            def header = lines[0].split('\t')
            def completeness = lines[1].split('\t')[header.findIndexOf { it == 'Completeness' }]
            "CheckM1_est_genome_completeness,${completeness},${ID}\n"
        }
        .collectFile(name: '9_checkm_completeness.csv', storeDir: CH_stepwise_counts, seed: COUNT_HEADER, cache: false, sort: false)

    // CheckM's lineage_wf tries both genetic codes (4 and 11) during marker gene prediction
    // and records the one it picked in storage/bin_stats.analyze.tsv, not the --tab_table completeness TSV.
    CH_translation_table = CHECKM_v1_1_9.out
        .map { ID, checkm_dir ->
            def bin_stats = file("${checkm_dir}/storage/bin_stats.analyze.tsv").text
            def translation_table = (bin_stats =~ /'Translation table':\s*(\d+)/)[0][1]
            tuple(ID, translation_table)
        }
    
    PROKKA_v1_14_6(CH_final_contigs.join(CH_translation_table))

    PROKKA_GFF_2_TSV(PROKKA_v1_14_6.out.gff.join(CH_final_contigs))
    CH_count_prokka = PROKKA_GFF_2_TSV.out.countfile.collectFile(name: '10_prokka_stats.csv', storeDir: CH_stepwise_counts, seed: COUNT_HEADER, cache: false, sort: false)

    // SSU (16S) recovery + CREST-style LCA classification against SILVA
    SSU_BLAST(CH_final_contigs)
    SSU_GET_GENE(SSU_BLAST.out.join(CH_final_contigs))
    SSU_CLASSIFIER(SSU_BLAST.out.join(SSU_GET_GENE.out))
    PARSE_CLASSIFIER(SSU_CLASSIFIER.out.classif)
    CH_count_ssu = PARSE_CLASSIFIER.out.countfile.collectFile(name: '11_ssu_classification.csv', storeDir: CH_stepwise_counts, seed: COUNT_HEADER, cache: false, sort: false)

    // GTDB-Tk taxonomy — skip assemblies below params.gtdbtk_min_bp total bases
    GTDBTK_v2_0_0(CH_final_contigs.filter({ countBases(it[1]) > (params.gtdbtk_min_bp as int) }))
    PARSE_GTDBTK(GTDBTK_v2_0_0.out)
    CH_count_gtdbtk = PARSE_GTDBTK.out.collectFile(name: '12_gtdbtk_stats.csv', storeDir: CH_stepwise_counts, seed: COUNT_HEADER, cache: false, sort: false)

    // Alaina's viral vs. cellular predictor
    //PROTEINS_VS_EGGNOG_5(PROKKA_v1_14_6.out.faa.filter({ it[1].size()>0 })) // ignore .faa containing no proteins
    //EGGNOG_HITS_TO_CELL_OR_VIRUS(PROTEINS_VS_EGGNOG_5.out.filter({ it[1].size()>2350 })) // only keep outputs where hitsTXT file is over 13 lines long (<= 13 means no hits)
    //LOG_CELL_OR_VIRUS(EGGNOG_HITS_TO_CELL_OR_VIRUS.out.countfile.collect())

    // Combine every stepwise count file into one, earliest stage first. Each one already
    // carries its own COUNT_HEADER line (from its own collectFile seed above) — strip
    // that off per file before stacking, then let this collectFile's own seed add it back at the end.
    CH_all_counts = CH_count_raw_reads
        .concat(CH_count_trimmed_reads, CH_count_complex_reads, CH_count_normalized_reads,
                CH_count_clean_reads, CH_count_raw_contigs, CH_count_trimmed_contigs,
                CH_count_final_contigs, CH_count_sag_stats, CH_count_checkm,
                CH_count_prokka, CH_count_ssu, CH_count_gtdbtk)
        .map { it.text.readLines().drop(1).join('\n') + '\n' }
        .collectFile(name: 'all_stepwise_counts.csv', storeDir: "${params.output}/sample_tracking", seed: COUNT_HEADER, cache: false, sort: false)

    // Viral classifiers run on the same fully decontaminated + trimmed assembly as the
    // annotation steps above; geNomad is further restricted to capsules whose longest contig
    // exceeds 1500 bp.
    if ( params.viral == true ) {

        CH_FINAL_CONTIGS = CH_final_contigs
        GENOMAD_v1_11_1(CH_FINAL_CONTIGS.filter({ maxContigLength(it[1]) > 1500 }))
        //PARSE_GENOMAD(GENOMAD_v1_11_1.out)
        //LOG_GENOMAD(PARSE_GENOMAD.out.collect())

        // Other viral tools
        VIRSORTER_v2_2_3(CH_FINAL_CONTIGS)
        CHECKV_v1_0_1(CH_FINAL_CONTIGS)
        DEEPVIRFINDER(CH_FINAL_CONTIGS)
    }
    ASSEMBLY_STATS_TABULATOR(CH_all_counts)
}



process SAMPLE_READS {
    tag "${library} fastq subsampled to detect barcodes"
    container 'quay.io/biocontainers/seqtk:1.2--1'
    errorStrategy 'finish'
    input: tuple val(library), path(r1), path(r2)
    output: tuple val(library), path("0_sampled_r2.fastq")
    script: "seqtk sample -s${params.sample_random_seed} ${r2} ${params.sample_num_reads} > 0_sampled_r2.fastq"
    stub:
    // FOR TESTING: seqtk's reservoir sample must stream the entire input regardless of the
    // target sample size, so it's slow against real multi-GB fastqs even under -stub-run.
    // Grabs a real prefix of R2 instead — downstream Pheniqs steps aren't stubbed and need
    // genuinely barcode-bearing sequence, not placeholder data.
    """
    zcat ${r2} | head -40000 > 0_sampled_r2.fastq
    """ }

process PHENIQS_MAKE_SAMPLE_CONFIG {
    tag "${library} barcodes recorded to json file"
    container 'quay.io/biocontainers/pandas:2.2.1'
    errorStrategy 'finish'
    publishDir { "${params.output}/sample_tracking/atrandi_demux/${library}" }, mode: params.publishmode
    input:
        tuple val(library), path(sample_fastq)
        path BC_D
        path BC_C
        path BC_B
        path BC_A
    output: tuple val(library), path('1_sample_pheniqs_config.json')
    script: "pheniqs_make_sample_config.py ${sample_fastq} ${BC_D} ${BC_C} ${BC_B} ${BC_A} ${params.sample_hamming_dist}" }

process PHENIQS_SAMPLE_DEMULTIPLEX {
    tag "${library} demultiplexing read subset to evaluate barcode frequencies"
    container 'quay.io/biocontainers/pheniqs:2.1.0--py39ha79081e_6'
    input: tuple val(library), path(sample_fastq), path(sample_config)
    output: tuple val(library), path("2_sample_demux.bam")
    shell:
    '''
    pheniqs mux --config !{sample_config} -R sample_pheniqs_report.txt -t !{task.cpus}
    mv sample_demux.bam ./2_sample_demux.bam
    ''' }

process PHENIQS_COUNT_SORT_SAMPLE {
    tag "${library}"
    container 'quay.io/biocontainers/samtools:1.24--h9dcdb79_1'
    publishDir { "${params.output}/sample_tracking/atrandi_demux/${library}" }, mode: params.publishmode
    input: tuple val(library), path(sample_demux_bam)
    output: tuple val(library), path("3_observed_bc_count.txt")
    shell:
    '''
    samtools view -h !{sample_demux_bam} | awk -F '\t' '{ for (i=1; i<=NF; i++) { if ($i ~ /^CB:Z:/) { split($i, tag, ":"); print tag[3]; } } }' | sort | uniq -c > 3_observed_bc_count.txt
    ''' }

process PHENIQS_PLOT_HIST {
    tag "${library}"
    container 'quay.io/biocontainers/mulled-v2-283013c53e9be2db71ac5442e35da355c979ca0e:db712fd85ab78376a65a1bfaa62f945d34b00413-0'
    publishDir { "${params.output}/sample_tracking/atrandi_demux/${library}" }, mode: params.publishmode
    input: tuple val(library), path(observed_bc_count)
    output: path "4_observed_bc_dist.png"
    script: "pheniqs_plot_observed.py ${observed_bc_count} ${params.read_threshold}" }

process PHENIQS_FILTER_BC_LIST {
    tag "${library}"
    container 'quay.io/biocontainers/pandas:2.2.1'
    publishDir { "${params.output}/sample_tracking/atrandi_demux/${library}" }, mode: params.publishmode
    input: tuple val(library), path(observed_bc_count)
    output: tuple val(library), path("5_filt_bc.txt")
    script: "pheniqs_filter_barcodes.py ${observed_bc_count} ${params.read_threshold} ${params.cell_threshold}" }

process PHENIQS_NAME_CAPSULES {
    tag "${library}"
    container 'quay.io/biocontainers/pandas:2.2.1'
    publishDir { "${params.output}/sample_tracking/atrandi_demux/${library}" }, mode: params.publishmode
    input:
        tuple val(library), path(filt_bc), path(r1), path(r2)
        path BC_D
        path BC_C
        path BC_B
        path BC_A
    output: tuple val(library), path('6_split_pheniqs_config.json')
    script: "pheniqs_name_capsules.py ${filt_bc} ${r1} ${r2} ${params.split_hamming_dist} ${library} ${BC_D} ${BC_C} ${BC_B} ${BC_A}" }

process PHENIQS_DEMULTIPLEX {
    tag "${library} demultiplexing by combinatorial barcode"
    // Retries with more memory instead of a fixed ceiling, same pattern as
    // CONTAM_READ_FINDER/CHECKM_v1_1_9 below — but unlike those two, 4.GB isn't an
    // empirically-observed OOM point, just a starting guess: this is the one process that
    // reads a whole (un-subsampled) pool while writing one fastq.gz pair per capsule
    // concurrently, so its footprint scales with both pool size and capsule count. Capped at
    // 3 attempts (4/8/12GB) to stop at this machine's real 12GB Docker Desktop ceiling —
    // retrying past that would just repeat the same OOM kill.
    memory 12.GB
    //memory { 4.GB * task.attempt }
    maxForks 1 // keeps memory-heavy retries from stacking across libraries regardless of environment
    errorStrategy 'finish' //{ task.exitStatus in [137, 140] ? 'retry' : 'terminate' }
    maxRetries 2
    container 'quay.io/biocontainers/pheniqs:2.1.0--py39ha79081e_6'
    publishDir { "${params.output}/sample_tracking/atrandi_demux/${library}" }, pattern: "*.json", mode: params.publishmode
    input: tuple val(library), path(split_config), path(r1), path(r2)
    output:
        tuple val(library), path('7_pheniqs_report.json'), emit: txt
        path('*.fastq.gz'), emit: list_fastq
    shell:
    '''
    pheniqs mux --config !{split_config} -R 7_pheniqs_report.json -t !{task.cpus} -B !{params.pheniqs_buffer_capacity}
    '''
    stub:
    // FOR TESTING: like SAMPLE_READS, this is slow under -stub-run regardless of --dev, since
    // --dev's take() only trims capsules *after* this process has already demuxed the whole
    // (un-subsampled) pool. Fakes params.dev_num_capsules capsules from a real slice of
    // r1/r2 instead of running pheniqs mux for real, so TRIM_BARCODE right after this (not
    // stubbed) still has real-shaped data to run trim_galore on. The report JSON is a minimal
    // but structurally valid stand-in, not real content — just enough for PHENIQS_PARSE_REPORT
    // (also not stubbed) to parse a barcode/count/index row without choking on an empty one.
    """
    for i in \$(seq 1 ${params.dev_num_capsules}); do
        zcat ${r1} | head -40 | gzip > ${library}_stub\${i}_r1.fastq.gz
        zcat ${r2} | head -40 | gzip > ${library}_stub\${i}_r2.fastq.gz
    done
    echo '{"cellular": [{"classified": [{"barcode": ["ACGTACGT"], "count": 10, "index": 1}]}]}' > 7_pheniqs_report.json
    """ }

process PHENIQS_PARSE_REPORT {
    tag "${library}"
    container 'quay.io/biocontainers/pandas:2.2.1'
    publishDir { "${params.output}/sample_tracking/atrandi_demux/${library}" }, mode: params.publishmode
    input: tuple val(library), path(split_report)
    output: tuple val(library), path('8_pheniqs_report.csv')
    script: "pheniqs_parse_report.py ${split_report}" }

process TRIM_BARCODE {
    tag "${ID} removing Atrandi barcode"
    container 'quay.io/microbiome-informatics/trim_galore:0.6.7'
    publishDir { "${params.output}/${ID}/reads_${ID}" }, pattern: "debarcoded_*fastq.gz", mode: params.publishmode
    input: tuple val(ID), path(r1), path(r2)
    output: tuple val(ID), path("debarcoded_${ID}_r1.fastq.gz"), path("debarcoded_${ID}_r2.fastq.gz")
    script:
    """
    trim_galore --clip_R2 ${params.barcode_trim_length} --paired -o output -j ${task.cpus} ${r1} ${r2}
    mv ./output/${ID}_r1_val_1.fq.gz ./debarcoded_${ID}_r1.fastq.gz
    mv ./output/${ID}_r2_val_2.fq.gz ./debarcoded_${ID}_r2.fastq.gz
    """ }

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
    // FOR TESTING: makes empty *zip and *html files so downstream processes pick up immediately without waiting for real FastQC to run (which is slow and not needed for stub testing)
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
    // FOR TESTING: makes dummy trimmed_*.fastq.gz files with 10 read pairs each (40 lines) so downstream processes can run on stub data
    """
    zcat ${r1} | head -40 | gzip > trimmed_${ID}_r1.fastq.gz
    zcat ${r2} | head -40 | gzip > trimmed_${ID}_r2.fastq.gz
    echo "Trimmed_readcount,20,${ID}" > 2_trimmed_${ID}.count
    """ }

process COMPLEXITY_FILTER {
    tag "${ID}"
    // Only real dependency is pysam (parmap was replaced with plain multiprocessing.Pool
    // below, and six was replaced with stdlib equivalents) - a stock biocontainers image
    // covers it with no extra install needed.
    container 'quay.io/biocontainers/pysam:0.24.0--py312hf5ad864_1'
    publishDir { "${params.output}/${ID}/reads_${ID}" }, pattern: "*fastq.gz", mode: params.publishmode
    input: tuple val(ID), path(r1), path(r2)
    output: tuple val(ID), path("pe_${ID}.fastq.gz"), emit: reads
    script: template 'complexity_filter.py'
    stub:
    // FOR TESTING: Makes a dummy pe_*.fastq.gz file with 10 read pairs (40 lines) so downstream processes can run on stub data
    """
    paste <(zcat ${r1} | head -40 | paste - - - -) <(zcat ${r2} | head -40 | paste - - - -) | tr '\\t' '\\n' | gzip > pe_${ID}.fastq.gz
    """ }

process KMERNORM_v1_0_0 {
    tag "${ID}"
    // Installing and providing `kmernorm` (or other normalization software) on PATH is the user's own responsibility
    container 'brwnj/kmernorm:v1.0.0'
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
    // FOR TESTING: Makes a dummy normalized_pe_*.fastq.gz file with 10 read pairs (40 lines) so downstream processes can run on stub data
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
    // curlimages/curl:8.21.0 (the previous image here) is Alpine-based with no /bin/bash at
    // all. That's fatal regardless of this process's own `shell` directive: Nextflow's Docker
    // executor always launches the container's entrypoint as `/bin/bash -c "..."` for its own
    // setup wrapper (verified directly — every task's .command.run does this, whatever image
    // or `shell` directive it uses; `shell` only controls the *inner* invocation of
    // .command.sh, nested inside that outer wrapper). No `shell` directive can work around an
    // outer wrapper the image can't even start. debian:bookworm-slim has bash (and md5sum) by
    // default, runs as root by default (no `-u root`/`--entrypoint` override needed, unlike
    // curlimages/curl), and its curl successfully reaches zenodo.org over HTTPS (verified) —
    // curl/unzip aren't preinstalled here, so the script installs them itself first.
    container 'debian:bookworm-slim'
    // enabled: !workflow.stubRun keeps stub output from interfering with actual DB files.
    publishDir { file(params.contam_ref_fasta).getParent() }, mode: 'copy', enabled: !workflow.stubRun
    cache false
    output:
    path("${file(params.contam_ref_fasta).getName()}"), emit: fasta
    path("${file(params.contam_ref_fasta).getName()}.{amb,ann,bwt,pac,sa}"), emit: index

    script:
    def base = file(params.contam_ref_fasta).getName()
    """
    set -e
    apt-get update -qq
    apt-get install -y -qq curl unzip
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
    // FOR TESTING: Makes dummy DB files so downstream processes pick up immediately without waiting for real download to run (which is slow and not needed for stub testing)
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
    // future runs find it via the auto-detect check instead of rebuilding it.
    publishDir { file(params.contam_ref_fasta).getParent() }, pattern: "${file(params.contam_ref_fasta).getName()}.*", mode: 'copy', enabled: !workflow.stubRun
    // Fancy 'enabled: !workflow.stubRun' guard above is enough to keep stub output from interfering with real DB files
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
    // FOR TESTING: Makes dummy BLAST db files so downstream processes pick up immediately without waiting for real makeblastdb to run (which is slow and not needed for stub testing)
    """
    touch ${fasta}.nhr ${fasta}.nin ${fasta}.nsq
    """
}

process CONTAM_READ_FINDER {
    tag "${ID}"
    // Retries with more memory instead of a fixed ceiling, so this adapts to whatever's
    // actually available (a laptop, CI, an HPC node) rather than encoding one machine's
    // Docker Desktop allocation. Starting attempt at 10.GB (not 5.GB): a dynamic memory
    // directive is part of Nextflow's task hash, so a task that only succeeds on a retry
    // gets cached under that retry's hash, not attempt 1's — every subsequent -resume
    // still checks attempt 1 first, misses, and re-OOMs before reaching the cached retry
    // again. Confirmed against real Atrandi capsule data: most real samples need >5GB and
    // were stuck re-running this loop on every single -resume; 10GB clears them on attempt 1.
    memory { 10.GB * task.attempt }
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
    // FOR TESTING: Makes dummy norm_*.sai and contam_*.sam files so downstream processes pick up immediately without waiting for real BWA to run (which is slow and not needed for stub testing)
    """
    touch norm_${ID}_r1.fastq.sai norm_${ID}_r2.fastq.sai
    printf '@HD\\tVN:1.6\\tSO:unsorted\\n' > contam_${ID}.sam
    """ }

process CONTAM_READ_REPORTER {
    tag "${ID}"
    container 'quay.io/biocontainers/samtools:1.24--h9dcdb79_1'
    publishDir { "${params.output}/${ID}/QC_${ID}" }, pattern: "contam_align_*.tsv", mode: params.publishmode
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
    // Only real dependency is pysam (six was replaced with stdlib equivalents).
    container 'quay.io/biocontainers/pysam:0.24.0--py312hf5ad864_1'
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
    // FOR TESTING: Makes a dummy 0_all_contigs_*.fasta file with one contig so downstream processes can run on stub data
    """
    printf '>stub_contig_1\\n' > 0_all_contigs_${ID}.fasta
    yes ACGT | head -400 | tr -d '\\n' >> 0_all_contigs_${ID}.fasta
    printf '\\n' >> 0_all_contigs_${ID}.fasta
    echo "Raw_contig_count,1,${ID}" > 6_all_contigs_${ID}.count
    """ }

process TRIM_CONTIGS {
    errorStrategy 'finish'
    tag "${ID}"
    container 'quay.io/biocontainers/biopython:1.84'
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
    container 'quay.io/biocontainers/biopython:1.84'
    publishDir { "${params.output}/${ID}/logs_${ID}" }
    input: tuple val(ID), path(contigs)
    output: path("sag_stats_${ID}.csv"), emit: countfile
    script:
    """
    #!/usr/bin/env python
    from Bio import SeqIO; from Bio.SeqUtils import gc_fraction
    max_contig_length = 0; STR_all_seq = ''
    for record in SeqIO.parse("${contigs}", "fasta"):
        STR_all_seq = STR_all_seq + record.seq          # read in the DNA, 1 contig at a time
        if len(record.seq) > max_contig_length:
            max_contig_length = len(record.seq)
    out=open("sag_stats_${ID}.csv", "a")
    print("Max_contig_length", str(max_contig_length), "${ID}", sep=",", file=out)
    print("Final_assembly_length", str(len(STR_all_seq)), "${ID}", sep=",", file=out)
    # gc_fraction() replaces the older, since-removed Bio.SeqUtils.GC(); it returns a
    # 0-1 fraction rather than a 0-100 percentage, so scale to match the old output.
    print("GC_content", str(round(gc_fraction(STR_all_seq)*100,2)), "${ID}", sep=",", file=out)
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
    // FOR TESTING: Makes a dummy contig_blast_*.tsv file with the correct header so downstream processes can run on stub data 
    """
    echo -e 'Query Seq-id\tSubject Seq-id\tPercentage of identical matches\tAlignment length\tNumber of mismatches\tNumber of gap openings\tStart of alignment in query\tEnd of alignment in query\tStart of alignment in subject\tEnd of alignment in subject\tExpect value\tBit score\tAll subject Seq-id(s)\tRaw score\tNumber of identical matches\tNumber of positive-scoring matches\tTotal number of gaps\tPercentage of positive-scoring matches\tQuery frame\tSubject frame\tAligned part of query sequence\tAligned part of subject sequence\tQuery sequence length\tSubject sequence length\tAll Subject Title(s)' > contig_blast_${ID}.tsv
    """ }

process CONTAM_CONTIG_REMOVER {
    tag "${ID}"
    // Needs Biopython + toolshed + interlap. No public image bundles all three, so this
    // takes the closest existing public image (Biopython) and installs the two small,
    // pure-Python, no-compiled-deps packages on top at task start - still an existing
    // public image, not a custom one we build/publish ourselves. Installed from inside
    // the template script itself (see templates/contam_contig_remover.py), not via
    // beforeScript - beforeScript runs on the HOST, not inside the container, for the
    // local executor this pipeline uses, so pip installed there would need to exist on
    // the host rather than in the image.
    container 'quay.io/biocontainers/biopython:1.84'
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
    // No maxForks: on a real multi-capsule pool (esp. under a grid executor) it serialized
    // CheckM to one task at a time, which also stalled PROKKA downstream (gated on
    // CH_translation_table from CHECKM_v1_1_9.out). At 6 GB/task the OOM-retry stacking it
    // was guarding against isn't a concern; cap concurrency via the executor's queueSize
    // (or a withName maxForks in a machine-specific config) if a laptop needs it.
    // Non-OOM failures ignore (was 'terminate' — a single junk-genome CheckM failure then
    // killed the whole pool run); that capsule just gets no completeness / no translation
    // table, so PROKKA skips it too.
    errorStrategy { task.exitStatus in [137, 140] ? 'retry' : 'ignore' }
    maxRetries 3
    container 'quay.io/biocontainers/checkm-genome:1.1.9--pyhdfd78af_0'
    publishDir { "${params.output}/${ID}/QC_${ID}" }, mode: params.publishmode
    input: tuple val(ID), path(contigs)
    output: tuple val(ID), path("checkm_${ID}")
    script:
    """
    # Concurrent CheckM tasks on one node hit "OSError: [Errno 98] Address already in use"
    # in mp.Manager(). The real fix is `singularity.newPidNamespace = false` in the run's
    # nextflow.config (see that file) -- without a private PID namespace CheckM's manager
    # process gets a unique host PID, so its abstract socket (\\0listener-<pid>-0) can't
    # clash. This private TMPDIR is just hygiene on top (clean matplotlib cache per task).
    mkdir -p /var/tmp/checkm_mp_${ID}_${task.attempt}
    export TMPDIR=/var/tmp/checkm_mp_${ID}_${task.attempt}
    mkdir tmp_dir; cp ${contigs} ./tmp_dir/final_contigs_${ID}.fasta
    checkm lineage_wf --reduced_tree -f checkm_${ID}/completeness_${ID}.tsv --tab_table -q -x fasta -t ${task.cpus} tmp_dir checkm_${ID}
    rm -r tmp_dir
    rm -rf /var/tmp/checkm_mp_${ID}_${task.attempt}
    """
    stub:
    """
    mkdir -p checkm_${ID}/storage
    printf "Bin Id\\tCompleteness\\tContamination\\n" > checkm_${ID}/completeness_${ID}.tsv
    printf "final_contigs_${ID}\\t99.9\\t0.1\\n" >> checkm_${ID}/completeness_${ID}.tsv
    printf "final_contigs_${ID}\\t{'Translation table': 11}\\n" > checkm_${ID}/storage/bin_stats.analyze.tsv
    """ }

process PROKKA_v1_14_6 {
    container 'quay.io/biocontainers/prokka:1.14.6--pl5262hdfd78af_1'
    errorStrategy 'finish'
    tag "${ID}"
    publishDir { "${params.output}/${ID}/annotation_${ID}" }, mode: params.publishmode
    input: tuple val(ID), path(contigs), val(translation_table)
    output:
        tuple val(ID), path("prokka_${ID}"), emit: dir
        tuple val(ID), path("prokka_${ID}/${ID}.gff"), emit: gff
        tuple val(ID), path("prokka_${ID}/${ID}.faa"), emit: faa
    shell:
    '''
    prokka --gcode !{translation_table} --outdir prokka_!{ID} --prefix !{ID} --locustag !{ID} --quiet --compliant --force --proteins !{params.prokka} --cpus !{task.cpus} !{contigs}
    ''' }

process PROKKA_GFF_2_TSV{
    tag "${ID}"
    errorStrategy 'finish'
    // pandas alone (the old container here) is missing Bio.SeqIO, which this script's
    // coding-density calc needs — mulled combo pins pandas=1.5.2, biopython=1.79, numpy=1.23.5.
    container 'quay.io/biocontainers/mulled-v2-1e9d4f78feac0eb2c8d8246367973b3f6358defc:ebca4356a18677aaa2c50f396a408343200e514b-0'
    publishDir { "${params.output}/${ID}/annotation_${ID}" }, mode: params.publishmode, pattern: "*tsv"
    input: tuple val(ID), path(gff), path(contigs)
    output:
        tuple val(ID), path("comprehensive_prokka_${ID}.tsv"), emit: tsv
        path("prokka_stats_${ID}.csv"), emit: countfile
    script: template "prokka_gff_2_tsv.py" }

process SSU_BLAST {
    tag "${ID}"
    container 'quay.io/biocontainers/blast:2.11.0--pl5262h3289130_1'
    cpus 4
    publishDir { "${params.output}/${ID}/annotation_${ID}/ssu_recovery_${ID}" }, mode: params.publishmode
    input: tuple val(ID), path(contigs)
    output: tuple val(ID), path("1_silva_blast_${ID}.tsv")
    script:
    "blastn -task megablast -query ${contigs} -db ${params.silva_blastdb} -num_alignments 10 -outfmt 6 -num_threads ${task.cpus} -out 1_silva_blast_${ID}.tsv"
    stub:
    // FOR TESTING: SILVA isn't present under -stub-run; hand SSU_GET_GENE an empty hit table
    // (a real no-hit blastn produces exactly that) so it and PARSE_CLASSIFIER still run for real.
    "touch 1_silva_blast_${ID}.tsv" }

process SSU_GET_GENE {
    tag "${ID}"
    // pysam ships the faidx/fetch API used below; same image CONTAM_READ_REMOVER already pins.
    container 'quay.io/biocontainers/pysam:0.24.0--py312hf5ad864_1'
    publishDir { "${params.output}/${ID}/annotation_${ID}/ssu_recovery_${ID}" }, mode: params.publishmode
    input: tuple val(ID), path(ssu_hits), path(contigs)
    output: tuple val(ID), path("1b_candidate_ssu_${ID}.fasta")
    script:
    """
    #!/usr/bin/env python
    import pysam
    from itertools import groupby
    BLAST6 = ["qseqid","sseqid","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore"]
    pysam.faidx("${contigs}")
    fa = pysam.FastaFile("${contigs}")
    with open("${ssu_hits}") as fh, open("1b_candidate_ssu_${ID}.fasta", "w") as fo:
        # input assumed sorted by query; take the best (first) hsp per query
        for query, qgroup in groupby(fh, key=lambda x: x.partition("\\t")[0]):
            for hsp in qgroup:
                toks = dict(zip(BLAST6, hsp.strip().split("\\t")))
                break
            if not toks.get("qseqid"):
                continue
            start = min(int(toks["qstart"]), int(toks["qend"]))
            end   = max(int(toks["qstart"]), int(toks["qend"]))
            fo.write(">{}\\n{}\\n".format(toks["qseqid"], fa.fetch(toks["qseqid"], start - 1, end)))
    """
    stub:
    "touch 1b_candidate_ssu_${ID}.fasta" }

process SSU_CLASSIFIER {
    tag "${ID}"
    container 'quay.io/biocontainers/biopython:1.84'
    // 'ignore' not 'finish' (pilot used 'finish'): on a full ~1900-capsule pool a single
    // malformed SSU shouldn't halt everything — that capsule just gets no SSU classification.
    errorStrategy 'ignore'
    publishDir { "${params.output}/${ID}/annotation_${ID}/ssu_recovery_${ID}" }, mode: params.publishmode
    input: tuple val(ID), path(ssu_hits), path(ssu)
    output:
        tuple val(ID), path("2_ssu_${ID}.fasta"), emit: fasta
        tuple val(ID), path("3_classification-16s_${ID}.tsv"), emit: classif
    script: template "ssu_classifier.py"
    stub:
    // FOR TESTING: needs the SILVA .map/.tree, absent under -stub-run. Empty outputs let
    // PARSE_CLASSIFIER (not stubbed) still run and emit its three "0" rows.
    "touch 2_ssu_${ID}.fasta 3_classification-16s_${ID}.tsv" }

process PARSE_CLASSIFIER {
    tag "${ID}"
    container 'quay.io/biocontainers/pandas:2.2.1'
    errorStrategy 'ignore'
    input: tuple val(ID), path(ssu_tsv)
    output: path("${ID}_ssu.count"), emit: countfile
    script:
    """
    #!/usr/bin/env python
    import os
    import pandas as pd
    SSU_1 = SSU_2 = SSU_3 = "0"
    if os.stat("${ssu_tsv}").st_size != 0:
        DF = pd.read_csv("${ssu_tsv}", header=None, sep="\\t")
        SSU_1 = DF.loc[0, 1]
        if len(DF) > 1: SSU_2 = DF.loc[1, 1]
        if len(DF) > 2: SSU_3 = DF.loc[2, 1]
    with open("${ID}_ssu.count", "w") as h:
        h.write("1_SSU_classification,%s,${ID}\\n2_SSU_classification,%s,${ID}\\n3_SSU_classification,%s,${ID}\\n" % (SSU_1, SSU_2, SSU_3))
    """ }

process GTDBTK_v2_0_0 {
    tag "${ID}"
    container 'quay.io/biocontainers/gtdbtk:2.0.0--pyhdfd78af_1'
    // 'ignore': a tiny/junk assembly that clears params.gtdbtk_min_bp but has no ORFs makes
    // gtdbtk classify_wf hard-exit 1 ("no genomes to process" — Prodigal called 0 genes).
    // That's deterministic, so retrying is pointless; that capsule just gets no GTDB row.
    errorStrategy 'ignore'
    // Sized for the Bigelow cluster (charlie). classify_wf's pplacer step against GTDB r207 is
    // the memory driver; tune to your own node sizes if running elsewhere.
    memory '128 GB'
    cpus 8
    publishDir { "${params.output}/${ID}/annotation_${ID}" }, mode: params.publishmode
    input: tuple val(ID), path(contigs)
    output: tuple val(ID), path("gtdbtk_classification_${ID}")
    script:
    """
    mkdir tmp_genome_dir
    cp ${contigs} tmp_genome_dir/final_contigs_${ID}.fasta
    export GTDBTK_DATA_PATH=${params.gtdb}
    gtdbtk classify_wf --genome_dir tmp_genome_dir --out_dir gtdbtk_classification_${ID} --cpus ${task.cpus} -x fasta
    rm -r tmp_genome_dir
    """
    stub:
    // FOR TESTING: the ~66GB GTDB r207 DB isn't present under -stub-run. Fake the two files
    // PARSE_GTDBTK reads so it (not stubbed) still parses a bacterial classification.
    """
    mkdir -p gtdbtk_classification_${ID}/identify
    printf 'user_genome\\tclassification\\n' > gtdbtk_classification_${ID}/gtdbtk.bac120.summary.tsv
    printf 'final_contigs_${ID}\\td__Bacteria;p__STUB;c__;o__;f__;g__;s__\\n' >> gtdbtk_classification_${ID}/gtdbtk.bac120.summary.tsv
    printf 'name\\tnumber_multiple_unique_genes\\n' > gtdbtk_classification_${ID}/identify/gtdbtk.bac120.markers_summary.tsv
    printf 'final_contigs_${ID}\\t0\\n' >> gtdbtk_classification_${ID}/identify/gtdbtk.bac120.markers_summary.tsv
    """ }

process PARSE_GTDBTK {
    tag "${ID}"
    container 'quay.io/biocontainers/pandas:2.2.1'
    input: tuple val(ID), path(gtdbtk_dir)
    output: path("gtdbtk_stats_${ID}.csv")
    script:
    """
    #!/usr/bin/env python
    from os.path import exists
    import pandas as pd
    classification_via_GTDBTk = 0
    multicopy_marker_genes = 0
    TSV_bac = "${gtdbtk_dir}/gtdbtk.bac120.summary.tsv"
    TSV_ar  = "${gtdbtk_dir}/gtdbtk.ar53.summary.tsv"
    M_bac   = "${gtdbtk_dir}/identify/gtdbtk.bac120.markers_summary.tsv"
    M_ar    = "${gtdbtk_dir}/identify/gtdbtk.ar53.markers_summary.tsv"
    if exists(TSV_bac) and not exists(TSV_ar):
        classification_via_GTDBTk = pd.read_csv(TSV_bac, sep="\\t").loc[0, "classification"]
        if exists(M_bac):
            multicopy_marker_genes = pd.read_csv(M_bac, sep="\\t").loc[0, "number_multiple_unique_genes"]
    elif exists(TSV_ar) and not exists(TSV_bac):
        classification_via_GTDBTk = pd.read_csv(TSV_ar, sep="\\t").loc[0, "classification"]
        if exists(M_ar):
            multicopy_marker_genes = pd.read_csv(M_ar, sep="\\t").loc[0, "number_multiple_unique_genes"]
    elif exists(TSV_ar) and exists(TSV_bac):
        print("warning: ambiguous whether this is archaeal or bacterial")
    else:
        print("warning: no GTDB-Tk classification files found")
    with open("gtdbtk_stats_${ID}.csv", "w") as out:
        print("classification_via_GTDBTk", classification_via_GTDBTk, "${ID}", sep=",", file=out)
        print("multicopy_marker_genes", multicopy_marker_genes, "${ID}", sep=",", file=out)
    """ }

process ASSEMBLY_STATS_TABULATOR {
    // Only real dependencies are pandas + numpy, and numpy comes bundled with this image.
    container 'quay.io/biocontainers/pandas:2.2.1'
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

    LIST_col_order = ["Sample_ID", "Raw_readcount", "Trimmed_readcount", "Complexity_filtered_readcount", "Normalized_readcount", "Contam_filtered_readcount", "Raw_contig_count", "Final_clean_contig_count", "Max_contig_length", "Final_assembly_length", "GC_content", "CheckM1_est_genome_completeness", "CDS", "tRNA", "percent_CDS_annotated", "average_CDS_length", "coding_density", "1_SSU_classification", "2_SSU_classification", "3_SSU_classification", "classification_via_GTDBTk", "multicopy_marker_genes"]

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

    # reindex (not [LIST_col_order]) so a metric that no sample produced this run — e.g. every
    # assembly fell below params.gtdbtk_min_bp, or no SSU gene was recovered — comes through as
    # an empty column instead of raising KeyError.
    DF_log = DF_log.reindex(columns=LIST_col_order)
    ## Warn the user when the metrics are from a stub run, so they don't mistake them for real data.
    DF_log = DF_log.map(lambda x: "${STUB_PREFIX}" + str(x))
    DF_log.to_csv(PATH_out, index=False)
    """ }

process GENOMAD_v1_11_1 {
  tag "${ID}"
  errorStrategy "terminate"
  // Camargo, A. P., Roux, S., Schulz, F., Babinski, M., Xu, Y., Hu, B., Chain, P. S. G., Nayfach, S., & Kyrpides, N. C. — Nature Biotechnology (2023), DOI: 10.1038/s41587-023-01953-y.
  container 'quay.io/biocontainers/genomad:1.11.1--pyhdfd78af_0'
  publishDir { "${params.output}/${ID}/annotation_${ID}" }, mode: "copy"
  cpus 4
  input: tuple val(ID), path(contigs)
  output: tuple val(ID), path("geNomad_${ID}")
  script: "genomad end-to-end --cleanup --threads ${task.cpus} --full-ictv-lineage --splits 8 ${contigs} geNomad_${ID} ${params.DB_genomad_v1_11_1}" }

process PROTEINS_VS_EGGNOG_5 {
  // container='docker://quay.io/biocontainers/hmmer:3.3.2--h87f3376_2'
  // installation notes: conda create --prefix /mnt/scgc/scgc_nfs/opt/common/anaconda3a/envs/hmmer_3.4 -c conda-forge -c bioconda hmmer=3.4 pandas numpy gzip
  beforeScript 'module load anaconda; source activate /mnt/scgc/scgc_nfs/opt/common/anaconda3a/envs/hmmer_3.4'
  conda '/mnt/scgc/scgc_nfs/opt/common/anaconda3a/envs/hmmer_3.4'
publishDir { "${params.output}/${ID}/annotation_${ID}/eggNOG_${ID}" }, mode: params.publishmode
  errorStrategy 'ignore'
  cpus 2 // 6
  memory "50.GB"
  tag "${ID}"
 
  input: tuple val(ID), path(faa)
 
  output: tuple val(ID), path("proteins_hmmsearch_v_EggNOGdb_${ID}.txt.gz")
 
  script:
  """
  hmmsearch -E 0.00001 --cpu ${task.cpus} -o stdout.log --tblout proteins_hmmsearch_v_EggNOGdb_${ID}.txt ${params.PATH_hmm} ${faa}
  rm stdout.log
  gzip proteins_hmmsearch_v_EggNOGdb_${ID}.txt
  """ }
  
process EGGNOG_HITS_TO_CELL_OR_VIRUS {
    container 'quay.io/biocontainers/pandas:2.2.1'
    publishDir { "${params.output}/${ID}/annotation_${ID}/eggNOG_${ID}" }, mode: params.publishmode
    memory '50.GB'
    errorStrategy 'terminate'
    input: tuple val(ID), path(hitsTXT)
    output:
        tuple val(ID), path("${ID}_proteins_hmmsearched_against_eggNOG.csv"), path("${ID}_eggnog.count"), emit: tsv                            
        path("${ID}_eggnog.count"), emit: countfile
    script:
    """
        #!/usr/bin/env python
    import pandas as pd
    import numpy as np
    import gzip

    TAB = "\\t"
    NEWLINE = "\\n"
    Sample_ID = "${ID}"
    
    PATH_hits = "${hitsTXT}"
    
    PATH_annot = "${params.PATH_annot}"
    DF_annot = pd.read_csv(PATH_annot, sep=TAB)
    
    PATH_out_csv = Sample_ID + "_proteins_hmmsearched_against_eggNOG.csv"
    
    PATH_out_countfile = Sample_ID + "_eggnog.count"
    
    def hmmer_to_DF(path, program="hmmsearch", format="tblout", verbose=False):   
        if format in {"tblout","domtblout"}:
            cut_index = {"tblout":18, "domtblout":22}[format]
            data = list()
            header = list()
            with gzip.open(path,'rt') as FILE:
                for line in FILE.readlines():
                    if line.startswith("#"):
                        header.append(line)
                    else:
                        row = list(filter(bool, line.strip().split(" ")))
                        row = row[:cut_index] + [" ".join(row[cut_index:])]
                        data.append(row)
            DF = pd.DataFrame(data)
            if not DF.empty:
                columns = ["target_name","target_accession","query_name","query_accession","e-value","score","bias","best_domain_e-value","best_domain_score","best_domain_bias","exp","reg","clu","ov","env","dom","rep","inc","query_description"]
                DF.columns = columns
        return DF
    
    DF_hits = hmmer_to_DF(PATH_hits,"hmmsearch","tblout")

    # Keep only the top hit for each query

    # extract query_accession from query filename, e.g.  2N57K.faa.final_tree.fa -> 2N57K
    DF_hits['query_accession'] = DF_hits['query_name'].str.split('.').str[0]
    
    # Drop unneeded cols
    DF_hits = DF_hits.drop(columns=['target_accession','query_description','bias','best_domain_e-value','best_domain_score','best_domain_bias','exp','reg','clu','ov','env','dom','rep','inc'])
    
    ## Merge along the accession
    DF = DF_hits.merge(DF_annot, left_on='query_accession', right_on='accession',how='left')
    
    DF = DF.rename(columns={'target_name':'protein',
            'query_name':'hit',
            'e-value':'evalue',
            'score':'bit',
            'description':'name',
            'category':'cat',
            'updated_viral':'virdom'})
    
    DF = DF[['protein','hit','bit','evalue','name','cat','domain','virdom']]
    
    DF.to_csv(PATH_out_csv)

    # Now, report info about the best hits specifically

    DF['bit'] = pd.to_numeric(DF['bit']) # make bitscore numeric
    # For each protein (a.k.a. query) get the hit with the highest bitsore
    idx = DF.groupby('protein')['bit'].idxmax()
    DF_best = DF.loc[idx].reset_index(drop=True)
    
    print("Summarizing hits to a logfile")
    
    LINE1 = "Viral_besthits_(eggNOG)," + str(len(DF_best.loc[DF_best['virdom'] == 'Virus' ])) + "," + Sample_ID + NEWLINE
    LINE2 = "Bacterial_besthits_(eggNOG)," + str(len(DF_best.loc[DF_best['virdom'] == 'Bacteria' ])) + "," + Sample_ID + NEWLINE
    LINE3 = "Eukaryotic_besthits_(eggNOG)," + str(len(DF_best.loc[DF_best['virdom'] == 'Eukarya' ])) + "," + Sample_ID + NEWLINE
    LINE4 = "Archaeal_besthits_(eggNOG)," + str(len(DF_best.loc[DF_best['virdom'] == 'Archaea' ]))  + "," + Sample_ID + NEWLINE
    
    with open(PATH_out_countfile, 'w') as f:
        f.writelines([LINE1, LINE2, LINE3, LINE4])
    """ }

process LOG_CELL_OR_VIRUS {
    publishDir {"${params.output}/sample_tracking/3_assemblies"}, mode: "copy"; errorStrategy 'terminate'; queue "normal"
    input: path(countfiles)
    output: path("10_cell_or_virus_stats.csv")
    shell: ''' echo "Metric,Count,Sample_ID" > 10_cell_or_virus_stats.csv; for LINE in !{countfiles}; do cat ${LINE} >> 10_cell_or_virus_stats.csv; done ''' }

process VIRSORTER_v2_2_3 {
  errorStrategy 'ignore'
  container 'docker://jiarong/virsorter:latest'
  publishDir { "${params.output}/${ID}/annotation_${ID}" }, mode: "copy"
//memory='50.GB'
  cpus 6
  tag "${ID}"

  input: tuple val(ID), path(fasta)
  output: tuple val(ID), path("virsorter_${ID}"), emit: DIR_input
  output: tuple val(ID), path("virsorter_${ID}/final-viral-score.tsv"), emit: tsv

  script:
  """
  virsorter run -w virsorter_${ID} -i ${fasta} --min-length 1500 -j ${task.cpus} all --db-dir "/mnt/scgc_nfs/ref/virsorter2/"
  """ }

process CHECKV_v1_0_1 {
  errorStrategy 'ignore'
  container='docker://quay.io/biocontainers/checkv:1.0.1--pyhdfd78af_0'
  publishDir "${params.output}/${ID}/annotation_${ID}", pattern: "checkv_${ID}", mode: "copy"
  cpus 1
  tag "${ID}"
  input: tuple val(ID), path(fasta)
  output: tuple val(ID), path("checkv_${ID}"), emit: dir
  output: tuple val(ID), path("viruses_and_proviruses_${ID}.fasta"), emit: fasta
  output: tuple val(ID), path("viruses_and_proviruses_${ID}.fasta"), path("checkv_${ID}/quality_summary.tsv"), emit: fasta_AND_tsv

  shell:
  '''
  cp !{fasta} ./!{ID}.fasta
  checkv end_to_end !{ID}.fasta checkv_!{ID} -t !{task.cpus} -d /mnt/scgc_nfs/ref/checkv/checkv-db-v1.0/
  # delete empty output fastas
  find . -type f -empty -print -delete
  # Combine the output fasta files
  for f in checkv_!{ID}/*iruses.fna; do (cat "${f}"; echo) >> viruses_and_proviruses_!{ID}.fasta; done
  rm -r checkv_!{ID}/tmp
  ''' }

process DEEPVIRFINDER {
  errorStrategy 'ignore'
  // No bioconda/biocontainers package exists for DeepVirFinder (upstream has no official
  // container either); this community image traces back to a public Dockerfile in
  // replikation/What_the_Phage (the peer-reviewed "What the Phage" pipeline) — verified by
  // pulling it and running `dvf.py --help` before adopting it. Its bundled models live at
  // /DeepVirFinder/models, not the tool's own default ./models, hence -m below.
  container 'replikation/deepvirfinder:latest'
  publishDir { "${params.output}/${ID}/annotation_${ID}/deepvirfinder_${ID}" }, mode: "copy"
  cpus 2
  tag "${ID}"

  input: tuple val(ID), path(fasta)
  output: tuple val(ID), path("deepvirfinder.tsv")

  script:
  """
  set +e
  set +o pipefail

  cp ${fasta} ./${ID}.fasta
  if dvf.py -i ${ID}.fasta -m /DeepVirFinder/models -o output -l 1500 -c ${task.cpus};
  then
    echo "DeepFirFinder worked, copying output files"
    cp output/${ID}.fasta_gt1500bp_dvfpred.txt ./deepvirfinder.tsv
  else
    echo "WARNING: DeepFirFinder failed, but making dummy files to trigger next process in pipeline."
    # Rationale: So that results from the parallel process--VirSorter--can still get analyzed by next process in pipeline.
    echo "name\tlen\tscore\tpvalue" > ./deepvirfinder.tsv  # Make empty table
  fi
  """ }

  /*
process PARSE_GENOMAD {
    tag "${ID}"
    errorStrategy: "terminate"
    container: 'brwnj/kmernorm:v1.0.0'
    input: tuple val(ID), path(DIR_geNomad)
    output: path("${ID}_geNomad_counts.csv")
    script:
    """
    #!/usr/bin/env python
    newline='\\n'
    ID="${ID}"
    PATH_out=ID+"_geNomad_counts.csv"

    import pandas as pd
    from glob import glob
    from os.path import exists

    PATH_plasmid_tsv = glob("${DIR_geNomad}/*_summary/*_plasmid_summary.tsv")[0]
    PATH_virus_tsv = glob("${DIR_geNomad}/*_summary/*_virus_summary.tsv")[0]

    NUM_virus = ''; LEN_virus = ''; NUM_plasmid = ''; LEN_plasmid = ''

    if exists(PATH_virus_tsv): 
        DF_virus = pd.read_csv(PATH_virus_tsv, sep="\t")
        NUM_virus = len(DF_virus)
        LEN_virus = DF_virus['length'].sum()
    else: print("No file matching this pattern: ${DIR_geNomad}/*_summary/*_virus_summary.tsv" )
    
    if exists(PATH_plasmid_tsv):
        DF_plasmid = pd.read_csv(PATH_plasmid_tsv, sep="\t")
        NUM_plasmid = len(DF_plasmid)
        LEN_plasmid = DF_plasmid['length'].sum()
    else: print("No file matching this pattern: ${DIR_geNomad}/*_summary/*_plasmid_summary.tsv" )
    
    #  log results
    with open(PATH_out, "w") as handle:
        handle.write("Viral_bp_geNomad,"+str(LEN_virus)+","+ID+newline+"Viruses_geNomad,"+ str(NUM_virus) +","+ID+newline+"Plasmid_bp_geNomad,"+ str(LEN_plasmid) +","+ID+newline+'Plasmids_geNomad,' +str(NUM_plasmid)+','+ID+newline)
    """ }

/*
process PARSE_DEEPVIRFINDER_AND_VIRSORTER {
  errorStrategy 'ignore'
  container 'brwnj/kmernorm:v1.0.0'
  cpus 1
  tag "${ID}"

  input: tuple val(ID), path(virsorter_TSV), path(deepvirfinder_TSV), path(fasta)
  output: tuple val(ID), path("bait_${ID}.fasta"), emit: fasta
  output: path("1_${ID}.log"), emit: log

  """
  #!/usr/bin/env python
  import pandas as pd; import shutil

  tab = "\\t"; newline = "\\n"
  PASS = False

  # Check virsorter results
  MAXvirsorter = pd.read_csv("${virsorter_TSV}",sep=tab)['max_score'].max()
  if MAXvirsorter > ${Virus_MINscore}: PASS = True

  # Check deepvirfinder results
  DF = pd.read_csv("${deepvirfinder_TSV}",sep=tab)
  DF_sig = DF.loc[DF['pvalue'] < ${Virus_MAXpvalue}]
  if len(DF_sig) > 0:
    MAXdeepvirfinder = DF_sig['score'].max()
    if len(DF_sig) > 0: PASS = True
  else: MAXdeepvirfinder = 'NA'

  if PASS == True:
    print("passed")
    shutil.copy("${fasta}", "bait_${ID}.fasta") # copy SAG to move forward with analysis
  else:
    print("Fail. No contigs had scores over ","${Virus_MINscore}"," and/or pvalues under ","${Virus_MAXpvalue}")

  # log results
  with open("1_${ID}.log", "w") as handle:
    handle.write("Bait,"+str(PASS)+",${ID}"+newline+"Maxdeepvirfinder,"+ str(MAXdeepvirfinder) +",${ID}"+newline+"MAXvirsorter,"+ str(MAXvirsorter) +",${ID}"+newline)
  """ }
*/