
from sys import stdout, stdin
from typing import *					#type hints
from contextlib import contextmanager
from os.path import isfile				#file exists checks
from itertools import batched


def attr_to_dict(attrs: str) -> dict[str:str]:
	"""
	Convert a GFF record's attributes field to a dictionary
	"""
	attr_dict = {}
	if not attrs or attrs == ".":
		return attr_dict
		
	#Split into individual attributes
	for attr in attrs.split(';'):
		#Split key-value pairs
		key, sep, value = attr.partition('=')
		attr_dict[key] = value
	
	return attr_dict


def read_tsv(path: str, comment_char: str = None, hlen: int = 0) -> list[str]:
	"""
	Take a TSV file and yield each line, broken into fields
	The first `hlen' lines, and lines starting with `comment_char' are skipped
	"""
	with ropen(path) as fd:
		# Skip initial lines if requested
		for _ in range(hlen):
			next(fd)
		
		for line in fd:
			line = line.rstrip()
			
			#Skip comments
			if comment_char and line.startswith(comment_char):
				continue
			
			#Split to fields
			yield line.split('\t')


def require_files(*args):
	"""
	Given one or more paths, check if all of them exist as files, and raise an exception if not
	"""
	for file in args:
		#Skip None and "-"
		if file and file != "-" and not isfile(file):
			raise FileNotFoundError(f"File {file} does not exist")


def refuse_files(*args):
	"""
	Given one or more paths, check if none of them exist as files, and raise an exception if yes
	"""
	for file in args:
		#Skip None and "-"
		if file and file != "-" and not isfile(file):
			raise FileExistsError(f"File {file} already exists")


@contextmanager
def ropen(path: str) -> TextIO:
	"""
	Open a file for reading, or return a handle to stdin if `path' is "-"
	"""
	if path == "-":
		yield stdin
	else:
		with open(path, 'r') as fd:
			yield fd


@contextmanager
def wopen(path: str) -> TextIO:
	"""
	Open a file for writing, or return a handle to stdout if `path' is "-"
	"""
	if path == "-":
		yield stdout
	else:
		with open(path, 'w') as fd:
			yield fd


def are_altseq(s1: str, s2: str, n: int = 1) -> bool:
	"""
	Given two transcripts' sequences, detect if they both originate from the same transcript
	"""
	#What follows is a description of how this test works
	#
	#Let S be the complete ground truth sequence of the transcript, and
	#let S1 & S2 be sequences obtained from two separate sequencings of that transcript.
	#This test assumes that sequencing the transcript can yield any continuous substring of S
	#(this implicitly assumes insertions, deletions & subsitutions), and so S1 & S2 are such
	#substrings
	#
	#
	#S, S1 & S2 can take one of the following configurations:
	#
	#S:		-------------------------------------------
	#S1:	 ---------------
	#S2:	                  ----------------
	#The two sequencings sequenced separate areas of the transcript
	#In this case, without access to S, it's impossible to determine whether S1 & S2 are truly
	#from separate sequencings of S, or if they come from two different transcripts.
	#
	#S:		-------------------------------------------
	#S1:	 ---------------------
	#S2:	                  ----------------
	#S1 & S2 overlap partially, with mutual overhangs.
	#In this case, there exists a suffix of S1 which is also a prefix of S2.
	#
	#S:		-------------------------------------------
	#S1:	 ------------------------------
	#S2:	               -------------
	#S1 overlaps S2 completely, with unilateral overhang(s).
	#In this case, S2 is a substring of S2.
	#
	#
	#This test can yield false negatives in the following case:
	#1. in the above pictured case #1, where S1 & S2 are separate areas of the transcript
	#2. where S1 & S2 come from overlapping aread of the same transcript, but were sequenced with
	#   error(s)
	#
	#This test can yield false positives in the following case:
	#1. if S1 & S2 come from two different but similar transcripts (e.g. two alternative
	#   transcripts of the same gene), and S1 & S2 both are subsequences of the regions where the
	#   two transcripts are identical
	#2. if S1 & S2 come from disparate transcripts but are very short, the sequences may meet the
	#   test conditions coincidentally - to prevent these cases, the `n' parameter might be used to
	#   pass a minimal length that either a subsequence or a shared suffix-prefix must be for the
	#   test to yield a positive result
	
	#Get the shorter and longer of s1 & s2
	shorter = min(s1, s2, key=len)
	longer =  max(s2, s1, key=len)
	
	
	#If the shorter sequence is shorter than n, there's no way to pass
	if len(shorter) < n:
		return False
	
	
	#Check subsequence
	if shorter in longer:
		return True
	
	#Check shared prefix-suffix (and vice versa)
	if shorter[:n] == longer[-n:] or longer[:n] == shorter[-n:]:
		return True
	
	#If neither conditions passed, return False
	return False


def truncate_seqids(fasta: dict[str:str]) -> dict[str:str]:
	"""
	Given a deserialized FASTA file, truncate the sequence IDs to the first word
	"""
	return { seqid.partition(' ')[0]:seq for seqid,seq in fasta.items() }


def max_orf_len(s: str) -> int:
	"""
	Examine the ORFs of a nucleotide sequence and return the length of the longest ORF
	Assumes standard translation table and only examines forward ORFs
	Also considers partial ORFs, i.e. assumes there's a start codon somewhere prior to the start of
	the sequence
	"""
	#ORF opening & closing codons
	START: set[str] = { "ATG" }
	STOP:  set[str] = { "TAA", "TAG", "TGA" }
	
	
	#Will keep track of the longest open frame
	max_len = 0
	

	#Process all three forward ORFS
	for orf in range(3):
		frame_open = True
		cur_len = 0
		
		for codon in batched('N'*orf + s, 3):
			codon = ''.join(codon).ljust(3, 'N')
			
			
			if frame_open:
				if codon in STOP:
					frame_open = False
					max_len = cur_len if cur_len > max_len else max_len
					cur_len = 0
				else:
					cur_len += 1
			else:
				if codon in START:
					frame_open = True
					cur_len = 1
		
	
		if frame_open:
			max_len = cur_len if cur_len > max_len else max_len
	
	
	return max_len