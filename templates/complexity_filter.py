#!/usr/bin/env python

PATH_r1="${r1}"
PATH_r2="${r2}"
STR_sample_id="${ID}"
FLOAT_complexity_threshold=float("${params.complexity_threshold}")
NUM_threads=int("${task.cpus}")


import os
import os.path as op
import shutil
import six
import tempfile
import gzip
from contextlib import contextmanager
from itertools import groupby
import contextlib
import itertools
import parmap
import time
from pysam import FastxFile
import multiprocessing

##################### SUPPORT FUNCTIONS ##################################
def multiprocess(f, iterable, *args, **kwargs):
	chunksize = kwargs.pop('chunksize', 1000)
	key = kwargs.pop('key', lambda k, l=itertools.count(): next(l)//chunksize)
	for k, g in itertools.groupby(iterable, key=key):
		yield parmap.map(f, g, *args, **kwargs)
@contextlib.contextmanager
def file_transaction(*rollback_files):
	exts = {".vcf": ".idx", ".bam": ".bai", "vcf.gz": ".tbi", ".fastq.gz": ".count"}
	safe_names, orig_names = _flatten_plus_safe(rollback_files)
	# remove any half-finished transactions
	remove_files(safe_names)
	try:
		if len(safe_names) == 1:
			yield safe_names[0]
		else:
			yield tuple(safe_names)
	# failure -- delete any temporary files
	except:
		remove_files(safe_names)
		remove_tmpdirs(safe_names)
		raise
	# worked -- move the temporary files to permanent location
	else:
		for safe, orig in zip(safe_names, orig_names):
			if os.path.exists(safe):
				shutil.move(safe, orig)
				for check_ext, check_idx in six.iteritems(exts):
					if safe.endswith(check_ext):
						safe_idx = safe + check_idx
						if os.path.exists(safe_idx):
							shutil.move(safe_idx, orig + check_idx)
		remove_tmpdirs(safe_names)
def remove_tmpdirs(fnames):
	for x in fnames:
		xdir = os.path.dirname(os.path.abspath(x))
		if xdir and os.path.exists(xdir):
			shutil.rmtree(xdir, ignore_errors=True)
def remove_files(fnames):
	for x in fnames:
		if x and os.path.exists(x):
			if os.path.isfile(x):
				os.remove(x)
			elif os.path.isdir(x):
				shutil.rmtree(x, ignore_errors=True)
def _flatten_plus_safe(rollback_files):
	tx_files, orig_files = [], []
	for fnames in rollback_files:
		if isinstance(fnames, six.string_types):
			fnames = [fnames]
		for fname in fnames:
			basedir = safe_makedir(os.path.dirname(fname))
			tmpdir = safe_makedir(tempfile.mkdtemp(dir=basedir))
			tx_file = os.path.join(tmpdir, os.path.basename(fname))
			tx_files.append(tx_file)
			orig_files.append(fname)
	return tx_files, orig_files
def safe_makedir(dname):
	if not dname:
		return dname
	num_tries = 0
	max_tries = 5
	while not os.path.exists(dname):
		try:
			os.makedirs(dname)
		except OSError:
			if num_tries > max_tries:
				raise
			num_tries += 1
			time.sleep(2)
	return dname
def readfx(fastx):
	fx = ""
	try:
		fx = FastxFile(fastx)
		for f in fx:
			yield f.name, f.sequence, f.quality
	finally:
		if fx:
			fx.close()
def complexity_filter_sequence(seq, threshold=0.05, alphabet='ACGT'):
	count_threshold = len(seq) * threshold
	for elem in alphabet:
		if seq.count(elem) < count_threshold:
			return False
	return True
def complexity_filter_pe_record(xxx_todo_changeme, threshold, alphabet):
	((n1, s1, q1), (n2, s2, q2)) = xxx_todo_changeme
	try:
		assert n1 == n2
	except AssertionError:
		assert n1[:-2] == n2[:-2], "Pairing failed due to sync: %s %s" % (n1, n2)
	pass1 = complexity_filter_sequence(s1, threshold, alphabet)
	pass2 = complexity_filter_sequence(s2, threshold, alphabet)
	if pass1 and pass2:
		if not n1.endswith('/1'): n1 += '/1'
		if not n2.endswith('/2'): n2 += '/2'
		records = "@%s\\n%s\\n+\\n%s\\n@%s\\n%s\\n+\\n%s" % (n1, s1, q1, n2, s2, q2)
		return records
	return ""
############################### CORE FUNCTIONS #############################################
def complexity_filter(r1, r2, out_file, threshold=0.05, cores=1):
	count_file = out_file + '.count'
	total_pairs = 0
	filtered_pairs = 0
	p = multiprocessing.Pool(cores)
	print("Complexity filtering. Threshold is %0.02f" % threshold)
	with file_transaction(out_file) as tx_out_file:
		with open(tx_out_file, 'w') as out_handle:
			for results in multiprocess(complexity_filter_pe_record, zip(readfx(r1), readfx(r2)), threshold, 'ACGT', pool=p):
				for record in results:
					total_pairs += 1
					if not record:
						filtered_pairs += 1
						continue
					print(record, file=out_handle)
					#print >>out_handle, record
	surviving_pairs = total_pairs - filtered_pairs
	print("Filtering complete. Observed pairs: %d. Total surviving: %d" % (total_pairs, surviving_pairs))
	return out_file, surviving_pairs
###########################################
complexity_filter(PATH_r1, PATH_r2, "pe_"+STR_sample_id+".fastq", FLOAT_complexity_threshold, NUM_threads)
with open("pe_"+STR_sample_id+".fastq", 'rb') as f_in:
	with gzip.open("pe_"+STR_sample_id+".fastq.gz", 'wb') as f_out:
		shutil.copyfileobj(f_in, f_out)

### delete the non-zipped version to save space:
print("Removing the non-zipped file to save space.")
os.remove("pe_"+STR_sample_id+".fastq")

print("Done")

#complexity_filter( "${r1}", "${r2}", "pe_${ID}.fastq", ${params.complexity_threshold}, ${task.cpus})
#with open("pe_${ID}.fastq", 'rb') as f_in:
	# with gzip.open("pe_${ID}.fastq.gz", 'wb') as f_out:
	# 	shutil.copyfileobj(f_in, f_out)

