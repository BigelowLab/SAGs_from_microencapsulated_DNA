	#!/usr/bin/env python
	import Bio; from Bio import SeqIO
	import os
	from collections import Counter
	MODE_verbose=False

	#### I. Trim and size-select sequences # -> 2a_long_passing_contigs_''', 2b_short_discarded_contigs_''', temp_trimmed.fasta

	trimmed_seqs=[]; short_seqs=[]; long_seqs=[]; num_long_seqs=0
	for record in SeqIO.parse("${contigs}", 'fasta'):
		if len(record) > int("${params.assembly_minlength}"):
			long_seqs.append(record); num_long_seqs += 1
			trimmed_seqs.append(record[int("${params.assembly_lefttrim}"):-int("${params.assembly_righttrim}")])
		else:
			short_seqs.append(record)

	SeqIO.write(long_seqs, "2a_long_passing_contigs_${ID}.fasta", 'fasta')
	SeqIO.write(short_seqs, "2b_short_discarded_contigs_${ID}.fasta", 'fasta')
	SeqIO.write(trimmed_seqs, "temp_trimmed.fasta", 'fasta')

	#### II. From surviving seqs (temp_trimmed.fasta), remove any contig that is an EXACT reverse complement of another contig. -> 3_trimmed_contigs_${ID}.fasta
	# These contigs are rare artifacts produced by SPAdes(v3) and will cause the genome to be rejected by NCBI/Genbank

	LIST_contigs_to_remove = []
	DICT_contigs_to_remove_2_length = {}

	def are_these_reverse_complements(TUP_seqIDs, DICT_contigs):
		'''
		Check two sequences (of exactly the same length) for whether they're EXACT reverse complements
		Returns True or False
		'''
		ID1 = TUP_seqIDs[0] # e.g. 'AM-124-K07_NODE_75'
		ID2 = TUP_seqIDs[1] # e.g. 'AM-124-K07_NODE_76'
		# Retrieve sequences from DICT_contigs
		seq1 = DICT_contigs[ID1].seq 
		seq2 = DICT_contigs[ID2].seq
		return(seq1 == seq2.reverse_complement())

	### Read trimmed fasta into dictionary (e.g. 'SeqID to SeqRecord')
	DICT_contigs = SeqIO.to_dict(SeqIO.parse("temp_trimmed.fasta", "fasta"))

	### For each SeqRecord, get Seq length
	DICT_seqlengths = {id:len(DICT_contigs[id].seq) for id in DICT_contigs} # e.g. {'AM-124-K07_NODE_1': 107265, 'AM-124-K07_NODE_2': 50472, ... }

	### Of these lengths, which are redundant? (I.e. >1 contigs have this exact length)
	LIST_lengths = [ DICT_seqlengths[key] for key in DICT_seqlengths.keys() ]
	LIST_redund_lengths = [length for length, count in Counter(LIST_lengths).items() if count > 1]

	#### IIa. Scan for sequences to remove. (I.e. the non-first seq of each pair that exhibits reverse complementarity)

	# Only comparing sequences that are the exact same length
	for length in LIST_redund_lengths:
		if MODE_verbose == True: print("Looking up all contigs with the length: "+str(length))
		LIST_seqIDs_with_this_length = [k for k,v in DICT_seqlengths.items() if float(v) == length]
		if MODE_verbose == True: print("Found "+str(len(LIST_seqIDs_with_this_length)))
	
		# Shorten listname for ease
		my_LIST = LIST_seqIDs_with_this_length 
	
		# Check each pair of each SeqIDs for reverse complementarity # But first pair of seqIDs like ['A','B','C'] -> [('A', 'B'), ('A', 'C'), ('B', 'C')]
		LIST_pairs_of_seqIDs=[(my_LIST[i], my_LIST[j]) for i in range(len(my_LIST)-1) for j in range(i+1,len(my_LIST))]

		for TUP_seqIDs in LIST_pairs_of_seqIDs:
			if are_these_reverse_complements(TUP_seqIDs, DICT_contigs) == True:
				if MODE_verbose == True: print(TUP_seqIDs[0]  + " and " + TUP_seqIDs[1] + " are exact reverse complements. Flagging second sequence for removal.")
				LIST_contigs_to_remove.append(TUP_seqIDs[1]) # Mark the second member of the pair for removal
				DICT_contigs_to_remove_2_length[TUP_seqIDs[1]] = length

	##### IIb. Remove the seqs and write out as new file:

	if MODE_verbose == False: print("Removing " + str(len(LIST_contigs_to_remove)) + " contigs that are exact reverse complements and keeping the rest.")

	LIST_records_to_write = []

	for record in SeqIO.parse("temp_trimmed.fasta", "fasta"):
		if record.id not in LIST_contigs_to_remove:
			LIST_records_to_write.append(record)
	### If seqID not in LIST_contigs_to_remove, write to new file

	SeqIO.write(LIST_records_to_write, "3_trimmed_contigs_${ID}.fasta", "fasta")

	if MODE_verbose == False: print("Done.")

	# How many contigs survived size-selection, trimming, and reverse-complement removal?
	with open("7_length_passing_contigs_${ID}.count", 'w') as count_handle:
		print("Contigs_over_min_length,"+str(len(LIST_records_to_write))+",${ID}", file=count_handle)

	# Cleanup
	if MODE_verbose == False: os.remove("temp_trimmed.fasta")

