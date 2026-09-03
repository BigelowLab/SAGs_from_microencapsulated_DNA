#!/usr/bin/env python
#
# CREST-style lowest-common-ancestor classifier for a recovered SSU (16S/18S) gene,
# run against a SILVA blast + the SILVA .map/.tree. Ported from
# BigelowLab/pilot_boothbay_capsules; the Python 2/3 `six` shim it used to carry was
# removed (this only ever runs on Python 3) and `Counter` is now imported explicitly
# (the original relied on it in BlastHits.majority() without importing it).

STR_sample_id="${ID}"
PATH_ssu="${ssu}"
PATH_ssu_hits="${ssu_hits}"
PATH_silva_map="${params.silva_map}"
PATH_silva_tree="${params.silva_tree}"

from collections import defaultdict, deque, OrderedDict, Counter
from Bio import Phylo; from Bio import SeqIO
import re; import os;
import bisect

RANKS = {8: "s", 7: "g", 6: "f", 5: "o", 4: "c", 3: "p", 2: "k"}
BLAST6 = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore",]
class OTU(object):
	def __init__(self, name, sequence, classification=None):
		self.name = name
		self.sequence = sequence
		self.classification = classification
class BlastHits(object):
	def __init__(self, names=None, max_hits=10, top_fraction=None):
		if names is None:
			# increasing bitscore sorted
			self.names = deque()
			self.percent_ids = deque()
			self.bitscores = deque()
		else:
			self.names = names
		self.max_hits = max_hits
		self.top_fraction = top_fraction
	def __repr__(self):
		return "{cls}[{tax}]".format(cls=self.__class__.__name__, tax=self.names)
	def __len__(self):
		return len(self.names)
	def add(self, name, percent_id, bitscore):
		bitscore = float(bitscore)
		if self.top_fraction and self.bitscores:
			# the filter
			if bitscore < (self.bitscores[-1] * self.top_fraction):
				bitscore = None
			# new best
			elif bitscore > self.bitscores[-1]:
				score = self.bitscores[0]
				while score < bitscore * self.top_fraction:
					self.names.popleft()
					self.percent_ids.popleft()
					self.bitscores.popleft()
					score = self.bitscores[0]
		if bitscore:
			# insert into sorted list
			idx = bisect.bisect_left(self.bitscores, bitscore)
			self.bitscores.insert(idx, bitscore)
			self.percent_ids.insert(idx, percent_id)
			self.names.insert(idx, name)
			if len(self.names) > self.max_hits:
				# remove lowest bitscore
				self.names.popleft()
				self.percent_ids.popleft()
				self.bitscores.popleft()
	def best_hit(self):
		return self.names[-1]
	def majority(self):
		# no repeated names
		if len(self.names) == len(set(self.names)):
			return self.best_hit()
		else:
			# count each taxonomy, grab top taxonomy
			most_common = Counter(self.names).most_common(1)[0][0]
			# need to flip to grab best bitscore
			names_reversed = self.names.copy()
			names_reversed.reverse()
			# left most index match
			idx = names_reversed.index(most_common)
			return names_reversed[idx]
class Tree(object):
	ROOT = 0; META = 1; DOMAIN = 2; SUPERKINGDOM = 3; KINGDOM = 4; PHYLUM = 5; CLASS = 6; ORDER = 7; FAMILY = 8; GENUS = 9; SPECIES = 10; SUBSPECIES = 11
	depths = {ROOT: "root", META: "meta", DOMAIN: "domain", SUPERKINGDOM: "superkingdom", KINGDOM: "kingdom", PHYLUM: "phylum", CLASS: "class", ORDER: "order", FAMILY: "family", GENUS: "genus", SPECIES: "species", SUBSPECIES: "strain",}
	def __init__(self, mapfile, trefile):
		self.tree = Phylo.read(
			trefile, "newick", values_are_confidence=True, rooted=True
		)
		self.root = self.tree.root; self.names = {}; self.node_names = {}; self.node_ids = {}; self.assignment_min = {}; self.parents = {}
		self.no_hits = self.add_node("No hits", self.root)
		for child in self.get_all_children(self.root):
			self.node_ids[child.name] = child
		# why not just remove these from the map?
		accession_re = [
			re.compile(r"\D\D\d\d\d\d\d\d\Z"),
			re.compile(r"\D\d\d\d\d\d\Z"),
			re.compile(r"\D\D\D\D\d\d\d\d\d\d\d\d\d\Z"),
			re.compile(r"\D\D\D\D\d\d\d\d\d\d\d\d\Z"),
		]
		# Read nodes from .map file (id\t name\t cutoff)
		with open(mapfile) as fh:
			for line in fh:
				toks = line.strip().split("\t")
				node_id = toks[0]
				name = toks[1]
				similarity_cutoff = float(toks[3])
				# Find node and map name or accession to it
				n = self.node_ids.get(node_id)
				if n:
					self.node_names[name] = n
					# Unless this is just an accession, update node name and assignment min.
					if similarity_cutoff >= 0 and not (
						accession_re[0].match(name)
						or accession_re[1].match(name)
						or accession_re[2].match(name)
						or accession_re[3].match(name)
					):
						self.assignment_min[name] = similarity_cutoff
						n.name = name
	def verify_node(self, node):
		if isinstance(node, Phylo.BaseTree.Clade):
			return node
		else:
			nn = node
			node = self.node_names.get(nn)
			return node
	def add_node(self, nodename, parent, assignment_min=0):
		# insert instead of append?
		node = Phylo.Newick.Clade(name=nodename)
		parent = self.verify_node(parent)
		if not parent:
			return None
		parent.clades.append(node)
		self.node_names[nodename] = node
		self.assignment_min[node.name] = assignment_min
		self.parents[node] = parent
		return node
	def get_immediate_children(self, node):
		node = self.verify_node(node)
		children = []
		for c in node.clades:
			children.append(c)
		return children
	def get_parent(self, node):
		node = self.verify_node(node)
		if not node:
			return None
		if node in self.parents:
			return self.parents[node]
		else:
			p = self.tree.get_path(node)
			if p and len(p) > 1:
				parent = p[-2]
				self.parents[node] = parent
				return parent
			else:
				return self.tree.root
	def get_rank(self, node):
		depth = self.get_depth(node)
		if depth > Tree.SUBSPECIES:
			depth = Tree.SUBSPECIES
		return Tree.depths[depth]
	def get_depth(self, node):
		node = self.verify_node(node)
		pth = self.get_path(node)
		if len(pth) == 1 and pth[0] is self.tree.root:
			return 0
		else:
			return len(pth)
	def get_path(self, node):
		plist = [node]
		if node is self.tree.root:
			return plist
		parent = self.get_parent(node)
		while parent and parent is not self.tree.root:
			plist = [parent] + plist
			parent = self.get_parent(parent)
		if not parent:
			print("Cannot find parent beyond %s" % plist)
		else:
			return plist
	def get_common_ancestor(self, node_names):
		nodes = list(set([self.node_names.get(n) for n in node_names]))
		if len(nodes) == 1:
			return nodes[0]
		paths = sorted([self.get_path(n) for n in nodes], key=len)
		lca_path = paths[0]
		for path in paths[1:]:
			while not lca_path[-1] in path:
				if len(lca_path) == 1:
					return self.tree.root
				else:
					lca_path.pop()
		return lca_path[-1]
	def get_all_children(self, node, children=None, parent=None):
		if children is None:
			children = []
		for c in node.clades:
			self.get_all_children(c, children, parent=node)
		children.append(node)
		self.parents[node] = parent
		return children
	def get_taxonomy(self, node):
		taxonomy = OrderedDict()
		for i in "kpcofgs":
			taxonomy[i] = "?"
		if node is not None:
			for clade in self.get_path(node):
				depth = self.get_depth(clade)
				if (
					depth > Tree.META
					and depth < Tree.SUBSPECIES
					and not depth == Tree.SUPERKINGDOM
					and not depth == Tree.KINGDOM
				):
					if Tree.depths[depth][0] == "d":
						abb = "k"
					else:
						abb = Tree.depths[depth][0]
					taxonomy[abb] = clade.name.replace(" ", "_")
		return taxonomy
def parse_blasthits(blasthits, otus, tre, min_score=155, top_fraction=0.98):
	hsps = defaultdict(lambda: BlastHits(top_fraction=top_fraction))
	with open(blasthits) as blast_hits_fh:
		# TODO iterator as some may not be able to process very large files
		for hsp in blast_hits_fh:
			toks = dict(zip(BLAST6, hsp.strip().split("\t")))
			if float(toks["bitscore"]) < min_score:
				continue
			hsps[toks["qseqid"]].add(
				toks["sseqid"], float(toks["pident"]) / 100, toks["bitscore"]
			)
	for otu_name, hits in hsps.items():
		otu = otus[otu_name]
		lca_node = tre.get_common_ancestor(hits.names)
		while (
			lca_node.name in tre.assignment_min
			and hits.percent_ids[-1] < tre.assignment_min[lca_node.name]
			and lca_node is not tre.root
		):
			lca_node = tre.get_parent(lca_node)
		otu.classification = lca_node
	return otus
def run_crest_classifier(fasta, blasthits, mapfile, trefile, outfasta, outtab, min_score=155, top_fraction=0.98, min_length=500,):
	otus = OrderedDict()
	for record in SeqIO.parse(fasta, "fasta"):
		otus[record.id] = OTU(record.id, record.seq)
	# No recovered SSU gene (empty blast / empty candidate fasta) — write empty outputs and stop.
	if not otus or os.stat(blasthits).st_size == 0:
		open(outfasta, "w").close(); open(outtab, "w").close()
		return
	tre = Tree(mapfile, trefile)
	otus = parse_blasthits(blasthits, otus, tre, min_score, top_fraction)
	with open(outfasta, "w") as fasta_out, open(outtab, "w") as tsv_out:
		for otu_id, otu in otus.items():
			if len(otu.sequence) < min_length:
				continue
			taxonomy = tre.get_taxonomy(otu.classification)
			full_name = "{name} {taxonomy}".format(
				name=otu.name,
				taxonomy=";".join(
					["%s__%s" % (abb, tax) for abb, tax in taxonomy.items()]
				),
			)
			print(full_name, file=fasta_out)
			print(otu.sequence, file=fasta_out)
			print(
				otu_id,
				";".join(["%s__%s" % (abb, tax) for abb, tax in taxonomy.items()]),
				sep="\t",
				file=tsv_out,
				)
###################################
run_crest_classifier(PATH_ssu, PATH_ssu_hits, PATH_silva_map, PATH_silva_tree, "2_ssu_"+STR_sample_id+".fasta", "3_classification-16s_"+STR_sample_id+".tsv")
