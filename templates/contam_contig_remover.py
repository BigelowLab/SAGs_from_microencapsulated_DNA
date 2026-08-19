#!/usr/bin/env python

STR_sample_id = "${ID}"
FLOAT_contam_min_percid = float("${params.contam_min_percid}")
INT_contam_min_length = int("${params.contam_min_length}")
PATH_blast_tsv = "${blast_tsv}"
PATH_contigs = "${contigs}"

import sys
import csv
import os.path as op
import os.path
from Bio import SeqIO
from itertools import groupby

# toolshed/interlap aren't in the biopython image this runs in - install them here rather
# than via a `beforeScript` directive, which runs on the host (not inside the container)
# for the local executor this pipeline uses.
import os
import subprocess
# Global site-packages isn't writable in this container, and $HOME (used for both pip's
# cache and its --user fallback) is read-only on some compute nodes - `--target` installs
# into an explicit, always-writable directory instead of relying on either.
PATH_pip_target = os.path.join(os.getcwd(), ".pylibs")
subprocess.run(
    [sys.executable, "-m", "pip", "install", "--quiet", "--no-cache-dir",
     "--target", PATH_pip_target, "toolshed==0.4.8", "interlap==0.2.7"],
    check=True,
)
sys.path.insert(0, PATH_pip_target)
from toolshed import nopen, reader
from interlap import InterLap

maxInt = sys.maxsize
while True:
	try:
		csv.field_size_limit(maxInt)
		break
	except OverflowError:
		maxInt =  int(maxInt/10)
print("max csv size = " + str(maxInt))
############ SUPPORT FUNCTIONS ###################
def overlap_distance(coords):
	overlap = InterLap()
	joined_intervals = InterLap()
	coords.sort(key=lambda x: x[0])
	overlap.add(coords)
	seen = set()
	for iset in coords:
		if iset in seen: continue
		seen.add(iset)
		start = iset[0]
		end = iset[1]
		for jset in overlap.find(iset):
			if jset in seen: continue
			seen.add(jset)
			start = min(start, jset[0])
			end = max([end, jset[1]])
		joined_intervals.add((start, end))
	return sum([j - i + 1 for (i, j) in joined_intervals])
def blast_results_filter(blast_tsv, fastx_file, out_file, contam_outfile, count_file):
	terms = ['homo sapiens', 'Homo sapiens', 'Homo_sapiens', 'human dna sequence', 'mus musculus', 'Delftia acidovorans', 'Achromobacter xylosoxidans', 'Stenotrophomonas maltophilia', 'Bradyrhizobium', 'AG-665_NODE_', 'AG-SIM-NEG_', 'AH-847', 'Cutibacterium', 'Malassezia']
	if isinstance(blast_tsv, str):
		blast_tsv = [blast_tsv]
	print("Filtering %s based on BLAST results." % fastx_file)
	if isinstance(terms, str):
		terms = [terms.lower()]
	else:
		terms = [t.lower() for t in terms]
	####### I. ####### pull out IDs of contaminated entries from blast results __________________________________
	contaminated = set()
	for btsv in blast_tsv:
		for qseq_id, alignments in groupby(reader(btsv, header=True, sep="\t"), key=lambda x: x['Query Seq-id']):
			###### get the coordinates of passing hits
			coords = set()
			for alignment in alignments:
				if (float(alignment['Percentage of identical matches']) > float(FLOAT_contam_min_percid) and
						any(t in alignment['All Subject Title(s)'].lower() for t in terms)):
					coords.add((int(alignment['Start of alignment in query']), int(alignment['End of alignment in query'])))
			###### sum potentially overlapping coordinates
			if coords and overlap_distance(list(coords)) > int(INT_contam_min_length):
				contaminated.add(qseq_id)
	####### II. ###### which seqIDs are not in contamlist? _________________________________________
	clean_contigs_list = []; contam_contigs_list = []
	for record in SeqIO.parse(fastx_file, "fasta"):
		if str(record.id) not in contaminated:
			clean_contigs_list.append(record)
		else:
			contam_contigs_list.append(record)
	####### III. How many were clean?
	with open(count_file, 'w') as count_handle:
		print("Final_clean_contig_count,"+str(len(clean_contigs_list))+","+STR_sample_id, file=count_handle)
	###### IV. ####### Filter out clean ___________________________________________________________________
	SeqIO.write(contam_contigs_list, contam_outfile, "fasta")
	if len(clean_contigs_list) > 0:
		SeqIO.write(clean_contigs_list, out_file, "fasta")
	else:
		print('No contigs survived decontamination')
############################################################
blast_results_filter(PATH_blast_tsv, PATH_contigs, "SCGC_"+STR_sample_id+"_contigs.fasta", "4b_contam_contigs_"+STR_sample_id+".fasta", "8_final_contigs_"+STR_sample_id+".count")


