
from sys import stdout, stdin
from typing import *					#type hints
from contextlib import contextmanager
from os.path import isfile				#file exists checks


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


def require_files(files: list[str]):
	"""
	Given an iterable of paths, check if all of them exist as files, and raise an exception if not
	"""
	for file in files:
		#Skip None and "-"
		if file and file != "-" and not isfile(file):
			raise FileNotFoundError(f"File {file} does not exist")


def refuse_files(files: list[str]):
	"""
	Given an iterable of paths, check if none of them exist as files, and raise an exception if yes
	"""
	for file in files:
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



