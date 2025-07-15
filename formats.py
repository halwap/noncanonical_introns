###################################################################################################
#   IMPORTS & SETUP
###################################################################################################

#Standard library
#Type hints
from __future__ import annotations
from typing import *
#I/O
from sys import stdout, stdin, stderr
from contextlib import contextmanager
from os.path import isfile
from pickle import load
from Bio.SeqIO.FastaIO import SimpleFastaParser
#Phase timing reporting
from time import time
#Convenient iteration
from itertools import chain


###################################################################################################
#	BASE I/O
###################################################################################################

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


def require_files(*paths):
	"""
	Given one or more paths, check if all of them exist as files, and raise an exception if not
	"""
	for file in paths:
		#Check if file exists; Nones, empty strings and "-" are skipped
		if file not in { '-', "", None } and not isfile(file):
			raise FileNotFoundError(f"File {file} does not exist")


def refuse_files(*paths):
	"""
	Given one or more paths, check if none of them exist as files, and raise an exception if yes
	"""
	for file in paths:
		#Check if file exists; Nones, empty strings and "-" are skipped
		if file not in { '-', "", None } and isfile(file):
			raise FileExistsError(f"File {file} already exists")



###################################################################################################
#	TSV PROCESSING
###################################################################################################

def to_tsv(*fields: Any, newline: bool = False) -> str:
	"""
	Given any number of values, return those values formatted as a line of a TSV file
	Does not include a newline at the end
	"""
	return '\t'.join( map( str, fields ) )


def parse_tsv(path: str, *, comment_char: str|None = '#', header: int = 0) -> Iterator[list[str]]:
	"""
	Take a TSV file and yield each line, broken into fields
	The first `header' lines, and lines starting with `comment_char' are skipped
	"""
	with ropen(path) as fd:
		# Skip initial lines if requested
		for _ in range(header):
			next(fd)
		
		for line in fd:
			#Skip comments
			if comment_char and line.startswith(comment_char):
				continue
			
			#Split to fields
			line = line.rstrip()
			yield line.split('\t')


def flatten_tsv(path: str, **kwargs) -> Iterator[str]:
	"""
	Yields all fields of every line of a TSV file as part of a single iterator
	"""
	return chain.from_iterable(parse_tsv(path, **kwargs))


def parse_list(path: str, **kwargs) -> Iterator[str]:
	"""
	Given a single-column TSV file, yield each line's contents as part of a single iterator
	"""
	for line in parse_tsv(path, **kwargs):
		assert len(line) == 1
		yield line[0]


###################################################################################################
#	GFF PARSING
###################################################################################################

#Type alias for GFF strand field
Strand: TypeAlias = Literal['+', '-', '.']


class GFF:
	"""
	Class for representing a singular GFF file entry
	"""
	seqid:  str
	source: Optional[str]
	type_:  str
	start:  int
	end:    int
	score:  Optional[float]
	strand: Strand
	phase:  Optional[int]
	attrs:  dict[str,str]
	
	
	def __init__(
		self,
		seqid:  str,
		source: str,
		type_:  str,
		start:  int,
		end:    int,
		score:  Optional[float] = None,
		strand: Strand          = '.',
		phase:  Optional[int]   = None,
		attrs:  dict[str,str]   = {}
	):
		"""
		Initialize a GFF entry
		"""
		#Validate input
		assert end >= start, \
			f"GFF entry has invalid coords: {start} start, {end} end"
		assert strand in "+-.", \
			f"GFF entry has invalid strand: {strand}"
		assert phase is None or phase in {0,1,2}, \
			f"GFF entry has invalid phase: {phase}"
		
		#Submit values to instance
		self.seqid  = seqid
		self.source = source
		self.type_  = type_
		self.start  = start
		self.end    = end
		self.score  = score
		self.strand = strand
		self.phase  = phase
		self.attrs  = attrs
	
	
	@classmethod
	def from_fields(cls, fields: list[str]) -> GFF:
		"""
		Initialize an instance based on a tab-split line taken from a GFF file
		"""
		assert len(fields) == 9, \
			f"GFF entry has wrong number of fields: {line}"
		
		return cls(
					seqid  = fields[0],
					source = fields[1],
					type_  = fields[2],
					start  = int(fields[3]),
					end    = int(fields[4]),
					score  = float(fields[5]) if fields[5] != '.' else None,
					strand = fields[6],
					phase  = int(fields[7]) if fields[7] != '.' else None,
					attrs  = GFF.attr_to_dict( fields[8] )
				)
	
	
	def __repr__(self) -> str:
		"""
		Convert entry to GFF format
		Does not include a trailing newline
		"""
		return to_tsv(
					self.seqid,
					self.source,
					self.type_,
					self.start,
					self.end,
					self.score if self.score != None else '.',
					self.strand,
					self.phase if self.phase != None else '.',
					GFF.dict_to_attr( self.attrs )
				)
	
	
	def __str__(self) -> str:
		return self.__repr__()
	
	
	@staticmethod
	def attr_to_dict(attrs: str) -> dict[str,str]:
		"""
		Convert an attributes field to a dictionary
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
	
	
	@staticmethod
	def dict_to_attr(attrs: dict[str,str]) -> str:
		"""
		Convert a dictionary to an attributes field as a single string
		"""
		return ';'.join( key+'='+value for key,value in attrs.items() ) if attrs else '.'
	
	
	def __len__(self) -> int:
		return self.end - self.start + 1
	
	
	#def __getattr__(self, key: str) -> str:
	#	"""
	#	Returns an attribute from self.attrs
	#	"""
	#	return self.attrs[key]
	#
	#
	#def __setattr__(self, key: str, value: Any):
	#	if hasattr(self, key):
	#		setattr(self, key, value)
	#	else:
	#		self.attrs[key] = value
	
	
	def orphan(self, *, backup: bool = True):
		"""
		Remove the Parent attribute of an entry, if it has one
		If `backup' is True, the parent's ID will be backed up as the "former_parent" attribute
		"""
		if "Parent" in self.attrs:
			if backup:
				self.attrs["former_parent"] = self.attrs["Parent"]
			del self.attrs["Parent"]
	
	
	def adopt(self, id_: str, *, backup: bool = True):
		"""
		Set the Parent attribute of an entry
		If `backup' is True, the previous parent's ID will be backed up as the "former_parent"
		attribute
		"""
		if "Parent" in self.attrs and backup:
			self.attrs["former_parent"] = self.attrs["Parent"]
		self.attrs["Parent"] = id_
	
	
	def overlaps(self, other: Self, *, unstranded: bool = False) -> bool:
		"""
		Determine if a pair of entries have any overlap between each other
		If `unstranded' is False, the entries must have the same strand to count as overlapping
		"""
		#Check that the entries have the same strand, if `unstranded' is False
		if not unstranded and self.strand != other.strand:
			return False
		
		#Check for overlap
		return self.seqid == other.seqid and self.start <= other.end and other.start <= self.end


def parse_gff(path: str) -> Iterator[GFF]:
	"""
	Read GFF file and yield each each contained entry
	`path' can be a path to a file, or '-' to signify stdin
	"""
	for fields in parse_tsv(path):
		#Instantiate a given entry
		yield GFF.from_fields(fields)



###################################################################################################
#	OTHER FORMATS
###################################################################################################

def deserialize_fasta(path: str, trunc: bool = False) -> dict[str,str]:
	"""
	Given a path to a FASTA file, deserialize it to a dictionary
	If `trunc' is True, the sequence identifiers will be truncated to only the first word
	"""
	#Load FASTA
	with ropen(path) as fd:
		fasta = dict( SimpleFastaParser(fd) )
	
	#Truncate sequence identifiers if requested
	if trunc:
		return { seqid.partition(' ')[0]:seq for seqid,seq in fasta.items() }
	else:
		return fasta


def load_model(path: str):
	"""
	Load a pickled sklearn model from filename
	Note that unlike other functions which take file paths, this one will NOT take '-' as stdin
	"""
	with open(path, 'rb') as fd:
		return load(fd)



###################################################################################################
#	STDERR OUTPUT
###################################################################################################

def eprint(*args, **kwargs):
	"""
	Print a message to stderr
	Sans the "file" parameter, semantics are identical to builtin print()
	"""
	print(*args, file=stderr, **kwargs)


def phase(new_phase: str|None = None):
	"""
	Function used to keep track of phases of the analyses. Calling this function marks the end of
	the previous phase (if there is one), and begins a new one (if a name for it was passed)
	"""
	#Report the timing of the last phase, if there was one
	if hasattr(phase, "timer"):
		#Move the printer cursor one line up, to the 40th character of the line
		eprint(f"\033[0;A\033[40;C{time() - phase.timer:.3f} sec")
	
	
	#If a new phase was requested
	if new_phase != None:
		#Print name of new phase
		eprint(new_phase.ljust(40))	
		
		#Reset the phase timer
		phase.timer = time()
	
	
	#If no new phase was requested
	else:
		#Remove the phase timer to effectively disable it
		if hasattr(phase, "timer"):
			del phase.timer
