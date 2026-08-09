#!/usr/bin/env python

STR_sample_id = "${ID}"
PATH_sam_contam = "${sam_contam}"
PATH_norm = "${norm}"

import gzip
import six
from six.moves import zip_longest
from pysam import FastxFile
import os
import os.path as op
import shutil
from os.path import exists

################## ACCESSORY FUNCTIONS ####################
def file_exists(fnames):
	if isinstance(fnames, six.string_types):
		fnames = [fnames]
	for f in fnames:
		if not os.path.exists(f) or os.path.getsize(f) == 0:
			return False
	return True
def grouper(n, iterable, fillvalue=None):
	args = [iter(iterable)] * n
	return zip_longest(fillvalue=fillvalue, *args)
def readfx(fastx):
	if not file_exists(fastx):
		print("File Not Found: %s" % fastx)
		#print "File Not Found: %s" % fastx
		raise IOError(2, "No such file:", fastx)
	fx = ""
	try:
		fx = FastxFile(fastx)
		for f in fx:
			yield f.name, f.sequence, f.quality
	finally:
		if fx:
			fx.close()
def check_sync(name1, name2):
	try:
		assert name1 == name2
	except AssertionError:
		assert name1[:-2] == name2[:-2], \
			"Run failed due to read sync: %s %s" % (name1, name2)
################### CORE SCRIPT ###########################
#### First check if there are any contaminants at all
observed=0; survived = 0
if os.stat(PATH_sam_contam).st_size == 0:
	print("Contamlist is empty. There are no contaminants to remove")
	src=PATH_norm; dst="contamfiltered_pe_"+STR_sample_id+".fastq.gz"
	#src=os.readlink(PATH_norm); dst="contamfiltered_pe_"+STR_sample_id+".fastq.gz"
	shutil.copy(src,dst)
	with open("5_contamfiltered_pe_"+STR_sample_id+".count", 'w') as count_file:
		print("Contam_filtered_readcount,NO_CHANGE,"+STR_sample_id, file=count_file)

else:
	### Read in list of contaminants
	with open(PATH_sam_contam) as f:
		hits = set([row.split()[0] for row in f])
	### Print out non-contaminant reads
	with open("contamfiltered_pe_"+STR_sample_id+".fastq", 'w') as fq_fh:
		for (n1, s1, q1), (n2, s2, q2) in grouper(2, readfx(PATH_norm)):
			check_sync(n1, n2)
			observed += 1
			# bwa strips read number from name so this works on R1 and R2
			read_name = n1.partition("/")[0]
			if read_name in hits: continue
			survived += 1
			print('@' + n1, s1, '+', q1, '@' + n2, s2, '+', q2, sep='\\n', file=fq_fh)
			#print >>fq_fh, '@' + n1 + '\\n' + s1 + '\\n' + '+\\n' + q1 + '\\n' + '@' + n2 + '\\n' + s2 + '\\n' + '+\\n' + q2
		#print("Pairs observed: %d" % observed); print("Pairs survived: %d" % survived)
		print(str(observed) + " paired observed")
		print(str(survived) + " pairs survived")
		print(str(observed - survived) + " contaminant pairs removed")
		with open("5_contamfiltered_pe_"+STR_sample_id+".count", 'w') as count_file:
			print("Contam_filtered_readcount,"+str(survived * 2)+","+STR_sample_id, file=count_file)

	############## Zip and cleanup #######

	if exists("contamfiltered_pe_"+STR_sample_id+".fastq"):
		with open("contamfiltered_pe_"+STR_sample_id+".fastq", 'rb') as f_in:
			with gzip.open("contamfiltered_pe_"+STR_sample_id+".fastq.gz", 'wb') as f_out:
				shutil.copyfileobj(f_in, f_out)

		### delete the non-zipped version to save space:
		print("Removing the non-zipped file to save space.")
		os.remove("contamfiltered_pe_"+STR_sample_id+".fastq")

	if survived == 0:
		os.remove("contamfiltered_pe_"+STR_sample_id+".fastq.gz"); print('No reads survived decontamination')

print("Done")


