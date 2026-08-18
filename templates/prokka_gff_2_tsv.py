#!/usr/bin/env python

STR_sample_id = "${ID}"
PATH_gff="${gff}"
PATH_contigs = "${contigs}" # Need this to get total assembly length
TAB = "\\t"

import pandas as pd
from Bio import SeqIO

with open(PATH_gff) as gff_fh, open("comprehensive_prokka_"+STR_sample_id+".tsv","w") as output:
	## print header to output tsv
	print("contig_id","locus_tag","ftype","start","end","strand","prokka_gene","prokka_EC_number","prokka_product","swissprot_gene","swissprot_EC_number","swissprot_product","swissprot_eggNOG","swissprot_KO","swissprot_Pfam","swissprot_CAZy","swissprot_TIGRFAMs",sep=TAB,file=output)

	for line in gff_fh:

		### initialize column values
		locus_tag=''; Name=''; prokka_gene=''; prokka_EC_number=''; prokka_product=''; swissprot_gene=''; swissprot_EC_number=''; swissprot_product=''; swissprot_eggNOG=''; swissprot_KO=''; swissprot_Pfam=''; swissprot_CAZy='';  swissprot_TIGRFAMs=''

		if line.startswith("##FASTA"):
			break
		if line.startswith("#"):
			continue
		LIST_of_fields = line.strip().split("\t")

		### Fill required columns
		#contig_id = LIST_of_fields[0]
		# change contig_id from "gnl|Prokka|AM-276-A01_1" -> "SCGC_AM-276-A01_contig1"
		contig_id = LIST_of_fields[0].replace('gnl|Prokka|','SCGC!').replace('_','_contig').replace('!','_')
		ftype = LIST_of_fields[2]
		start = LIST_of_fields[3]
		end = LIST_of_fields[4]
		strand = LIST_of_fields[6]

		## Skip features that are just genes
		if ftype == "gene": continue

		## Read annots into a dictionary
		DICT_of_annots = {}
		LIST_of_annots = LIST_of_fields[-1].split(';')   #annots are always in last field output by prokka
		for i in LIST_of_annots:
			split_keyvalue = i.split('=') # key is on right of '=' value is onleft
			if len(split_keyvalue) != 2 and MODE_verbose == True:
				print(split_keyvalue)
			if len(split_keyvalue) == 2:
				key = split_keyvalue[0]
				value = split_keyvalue[1]
				DICT_of_annots[key] = value # load into dictionary

		## The 'product' annotation may contain additional annots (eggNOG, KO, Pfam, CAZy, & TIGRFAMs) all sep by '^^'
		if 'product' in DICT_of_annots.keys():
			STR_product = DICT_of_annots.get('product')
			if "^^" not in STR_product:
				prokka_product = STR_product

			else:
				## Load in those additional annotations
				STR_product = STR_product.replace("%2",",")
				for LIST_additional_annots in STR_product.split("^^"):
					items = LIST_additional_annots.split("::")
					DICT_of_annots[items[0]] = items[1]

				prokka_product = DICT_of_annots.get('product')

		## Fill optional columns
		if 'locus_tag' in DICT_of_annots.keys(): locus_tag = DICT_of_annots.get('locus_tag')
		if 'gene' in DICT_of_annots.keys(): prokka_gene = DICT_of_annots.get('gene')
		if 'eC_number' in DICT_of_annots.keys(): prokka_EC_number = DICT_of_annots.get('eC_number')
		if 'eggNOG' in DICT_of_annots.keys(): swissprot_eggNOG = DICT_of_annots.get('eggNOG')
		if 'KO' in DICT_of_annots.keys(): swissprot_KO = DICT_of_annots.get('KO')
		if 'Pfam' in DICT_of_annots.keys(): swissprot_Pfam = DICT_of_annots.get('Pfam')
		if 'CAZy' in DICT_of_annots.keys(): swissprot_CAZy = DICT_of_annots.get('CAZy')
		if 'TIGRFAMs' in DICT_of_annots.keys(): swissprot_TIGRFAMs = DICT_of_annots.get('TIGRFAMs')

		print(contig_id,locus_tag,ftype,start,end,strand,prokka_gene,prokka_EC_number,prokka_product,swissprot_gene,swissprot_EC_number,swissprot_product,swissprot_eggNOG,swissprot_KO,swissprot_Pfam,swissprot_CAZy,swissprot_TIGRFAMs,sep="\t",file=output)

###### Collect core stats for logging / sample summary
DF=pd.read_table("comprehensive_prokka_"+STR_sample_id+".tsv")
count_handle = open("prokka_stats_"+STR_sample_id+".csv","w")

## basic counts
CDS = len(DF.loc[DF['ftype']=='CDS'])
print("CDS",CDS,STR_sample_id,sep=",",file=count_handle)
## print("rRNA",str(len(DF.loc[DF['ftype']=='rRNA'])),STR_sample_id,sep=",",file=count_handle) # prokka's rRNA counts are unreliable
print("tRNA",str(len(DF.loc[DF['ftype']=='tRNA'])),STR_sample_id,sep=",",file=count_handle)

## calc percent_CDS_annotated and average_CDS_length
DF_CDS=DF.loc[DF['ftype']=='CDS']
if len(DF_CDS) > 0:
	if 'hypothetical protein' not in DF_CDS['prokka_product'].tolist():
		percent_CDS_annotated = 100
	else:
		NUM_annotated = CDS - DF_CDS['prokka_product'].value_counts()['hypothetical protein']
		percent_CDS_annotated = round(NUM_annotated / CDS * 100, 2)
	print("percent_CDS_annotated",str(percent_CDS_annotated),STR_sample_id,sep=",",file=count_handle)

	DF_CDS['length']=DF_CDS['end']-DF_CDS['start']
	average_CDS_length=round(DF_CDS['length'].mean(), 2)
	print("average_CDS_length",str(average_CDS_length),STR_sample_id,sep=",",file=count_handle)

## calc coding density
STR_all_seq = ''
for record in SeqIO.parse(PATH_contigs, "fasta"): STR_all_seq = STR_all_seq + record.seq          # read in the DNA, 1 contig at a time
final_assembly_length = len(STR_all_seq)

if CDS != 0:
	coding_density = round(CDS * average_CDS_length / final_assembly_length*100,2)
else:
	coding_density = 0

print("coding_density",coding_density,STR_sample_id,sep=",",file=count_handle)


print(DF_CDS['prokka_product']) ###gremlin
count_handle.close()
