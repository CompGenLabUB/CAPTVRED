#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
File   		: rbbh.py
Author 		: Dominik R. Laetsch, dominik.laetsch at gmail dot com

Modifier by     : Josep F Abril
"""

from __future__ import division
import sys
import os
import re

class ProteomeObj():
	def __init__(self):
		self.seq_by_id = {}
		self.count_by_species = {}

	def parse_fasta(self, fastas):
		for fasta in fastas:
			aa_id = ''
			aa_species = ''
			aa_seq = ''
			with open(fasta) as fh:
				for line in fh:
					if line.startswith(">"):
						if (aa_seq):
							self.seq_by_id[aa_id] = seq
							self.count_by_species[species] = self.count_by_species.get(species, 0) + 1
							aa_id, aa_seq, aa_species = '', '', ''

						aa_id = line.rstrip("\n").lstrip(">")
						aa_species = line.rstrip("\n").lstrip(">").split('|')[0]
					else:
						seq += line.rstrip("\n")
				self.seq_by_id[aa_id] = seq

class RbbhCollection():
	def __init__(self):
		self.rbbh = []
		self.seen = set()

	def parse_hits(self, rbhCollection):
		for rbhObj in sorted(rbhCollection.rbhs.values(), key=lambda x: x.mean_of_bitscores, reverse=True):
			a, b = rbhObj.pair
			if rbhObj.is_valid():
				if a in self.seen:
					pass
				elif b in self.seen:
					pass
				else:
					self.rbbh.append(rbhObj)
			self.seen.update(rbhObj.proteins)

	def print_hits(self):
		print "#%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % ("A", "B", "MeanBitscore", "A_bitscore", "B_bitscore", "A_eval", "B_eval", "A_qcov", "B_qcov", "A_pid", "B_pid")
		for rbbhObj in self.rbbh:
			a, b = rbbhObj.pair
			if b.startswith('GROS'):
				a, b = b, a
			print "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (a, b, rbbhObj.mean_of_bitscores, rbbhObj.bitscores[0], rbbhObj.bitscores[1], rbbhObj.evalues[0], rbbhObj.evalues[1], rbbhObj.qcovs[0], rbbhObj.qcovs[1], rbbhObj.pidents[0], rbbhObj.pidents[1])



class RbhObj():
	def __init__(self, pair):
		self.pair = tuple(pair)
		self.hits = []
		self.hsps = []
		self.bitscores = []
		self.evalues = []
		self.sum_of_bitscores = 0
		self.qcovs = []
		self.qcovhsp = []
		self.sum_of_qcovs = 0
		self.sum_of_qcovhsp = 0
		self.proteins = set()
		self.oneway = 0
		self.mean_of_bitscores = 0
		self.pidents = []

	def bless(self):
		for hsp in self.hsps:
			if not hsp.qseqid in [x.qseqid for x in self.hits]:
				self.hits.append(hsp)
		self.proteins = {hitObj.qseqid for hitObj in self.hits}
		self.oneway = True if len(self.proteins) == 1 else False
		self.bitscores = [x.bitscore for x in self.hits]
		self.pidents = [x.pident for x in self.hits]
		self.evalues = [x.evalue for x in self.hits]
		self.sum_of_bitscores = sum(self.bitscores)
		self.mean_of_bitscores = sum(self.bitscores)/len(self.bitscores)
		self.qcovs = [x.qcovs for x in self.hits]
		self.qcovhsp = [x.qcovhsp for x in self.hits]
		self.sum_of_qcovs = sum(self.qcovs)
		self.sum_of_qcovhsp = sum(self.qcovhsp)

	def is_valid(self):
		if self.oneway == True:
			return False
		if not any(x <= evalue_theshold for x in self.evalues):
			return False
		if not all(x >= cov_theshold for x in self.qcovs):
			return False
		return True

class HitObj():
	def __init__(self, blasthit):


		# qseqid				sseqid					pident	mismatch	gapopen	qstart	qend	sstart	send	evalue	bitscore	qcovs	qcovhsp	qlen	slen	length
		# GPAL|GPLIN_000000100  GROS|GROS_g13706.t1     55.77   41      	3       81      184     3       101     2e-18   77.4    	55      55      190     108     104
		# GPAL|GPLIN_000001400  GROS|GROS_g06997.t1     73.58   20      	5       1       304     1       386     0.0     546   		100     100     304     386     386
		# qcovs 	= percent query length covered by subject
		# qcovhsp 	= percent query length covered by this HSP
                #
		# self.pair = frozenset([blasthit[0], blasthit[1]])
		# self.qseqid = blasthit[0]
		# self.sseqid = blasthit[1]
		# self.pident = float(blasthit[2])
		# self.mismatch = int(blasthit[3])
		# self.gapopen = int(blasthit[4])
		# self.qstart = int(blasthit[5])
		# self.qend = int(blasthit[6])
		# self.sstart = int(blasthit[7])
		# self.send = int(blasthit[8])
		# self.evalue = float(blasthit[9])
		# self.bitscore = float(blasthit[10])
		# self.qcovs = float(blasthit[11])
		# self.qcovhsp = float(blasthit[12])
		# self.qlen = int(blasthit[13])
		# self.slen = int(blasthit[14])
		# self.length = int(blasthit[15])

                #    dd_Smed_v6_10001_0_1    930     SMEST055666001.1        930
                #    10      925     14      929     9161792     0.0     1617    99.127  908
                #    8       908     0       0       99.13   1       1   TTTTTTTCA    TTTTTTTCA 
                # $BLOUTFMT
                # 6 qseqid.0 qlen.1 sseqid.2 slen.3 qstart.4 qend.5 sstart.6 send.7 length.8 score.9
                #   evalue.10 bitscore.11 pident.12 nident.13 mismatch.14 positive.15 gapopen.16
                #   gaps.17 ppos.18 qframe.19 sframe.20 qcovs.21 qcovhsp.22 qseq.23 sseq.24
		self.pair = frozenset([blasthit[0], blasthit[2]])
		self.qseqid = blasthit[0]
		self.sseqid = blasthit[2]
		self.pident = float(blasthit[12])
		self.mismatch = int(blasthit[14])
		self.gapopen = int(blasthit[16])
		self.qstart = int(blasthit[4])
		self.qend = int(blasthit[5])
		self.sstart = int(blasthit[6])
		self.send = int(blasthit[7])
		self.evalue = float(blasthit[10])
		self.bitscore = float(blasthit[11])
		self.qcovs = float(blasthit[21])
		self.qcovhsp = float(blasthit[22])
		self.qlen = int(blasthit[1])
		self.slen = int(blasthit[3])
		self.length = int(blasthit[8])

class RbhCollection():
	def __init__(self):
		self.rbhs = {} # dict of lists of HitObjs

	def add_hsp(self, hitObj):
		if hitObj.evalue < evalue_theshold:
			if not hitObj.pair in self.rbhs:
				self.rbhs[hitObj.pair] = RbhObj(hitObj.pair)
			self.rbhs[hitObj.pair].hsps.append(hitObj)

	def parse_blast(self, blasts):
		for blast in blasts:
			with open(blast) as fh:
				for line in fh:
					blasthit = [x.strip() for x in line.rstrip("\n").split("\t")] # cleaning leading and trailing whitespaces, and turning it into a list ...
					try:
						self.add_hsp(HitObj(blasthit))
					except TypeError:
						sys.exit("Format should be : qseqid, sseqid, pident, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore, qcovs, qcovhsp, qlen, slen, length (16)\n%s (%s)" % (line.rstrip("\n"), len(line.rstrip("\n").split("\t"))))

	def bless_rbhs(self):
		for rbhObj in self.rbhs.values():
			rbhObj.bless()


if __name__ == "__main__":


	fastas = []
	blasts = []

	try:
		blasts.append(sys.argv[1])
		fastas.append(sys.argv[2])
		blasts.append(sys.argv[3])
		fastas.append(sys.argv[4])
		evalue_theshold = float(sys.argv[5])
		cov_theshold = float(sys.argv[6])
	except:
		sys.exit("./rbbh.py [BLAST1] [FASTA1] [BLAST2] [FASTA2] [EVAL] [COV]")

	rbhCollection = RbhCollection()
	rbhCollection.parse_blast(blasts)
	rbhCollection.bless_rbhs()

	rbbhCollection = RbbhCollection()
	rbbhCollection.parse_hits(rbhCollection)
	rbbhCollection.print_hits()
