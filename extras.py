#Type hints
from typing import *

from collections import defaultdict
from itertools import chain, batched



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
	
	
	#Check shared prefix-suffix (both ways)
	if shorter[:n] == longer[-n:] or longer[:n] == shorter[-n:]:
		return True
	
	#Check subsequence
	if shorter in longer:
		return True
	
	#If neither conditions passed, return False
	return False



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




def merge_sets(sets: Collection[Collection[Hashable]]) -> list[set[Hashable]]:
	"""
	Given a set of sets, merge ones which have elements in common, and return a set of completely
	disjoint sets
	"""
	#This dictionary's keys are all the unique items across all sets in `sets'
	#Each value is the given key's parent
	parents: dict[Hashable,Hashable] = {}
	

	#Populate `parents' by iterating through all items across all sets
	for item in chain.from_iterable(sets):
		#Initially, each item is its own parent
		parents[item] = item
	
	
	def find_root(item: Hashable) -> Hashable:
		"""
		Find the root of an item, i.e. its great-great-great-...-grandparent
		Also automatically shortens the path through `parents' from `item' to its root
		"""
		while parents[item] != item:
			parents[item] = parents[parents[item]]
			item = parents[item]
		return item
	
	
	#Modify `parents' such that if two items come from the same set, they have the same root
	for set_ in sets:
		#Skip single-element & empty sets
		if len(set_) > 1:
			
			set_: Iterator[Hashable] = iter(set_)
			
			
			#Get the first (in practice, arbitrary) item of `set_', along with its root
			item1: Hashable = next(set_)
			root1: Hashable = find_root(item1)
			
			
			#Iterate through the remaining items in `set_'
			for item2 in set_:
				#Get the root of the other item
				root2: Hashable = find_root(item2)
				
				#Daisy chain the tree that `item2' is in, to the tree that `item1' is in
				#This effectively merges these two trees, albeit the resultant tree is potentially
				#unbalanced
				if root1 != root2:
					parents[root2] = root1
	
	
	#Group items into new sets, such that all items with the same root are in one set
	merged_sets: dict[Hashable,set[Hashable]] = defaultdict(set)
	for item in parents:
		merged_sets[ find_root(item) ].add(item)
	
	return list( merged_sets.values() )



def max_of_each(*iters: Iterable[Any], key: Callable = lambda x: x) -> Iterator[Any]:
	"""
	Given multiple iterables, find the max-valued item in each iterable, and yield that item
	Like with max(), a custom function to apply on each may be supplied with the `key' argument
	To use with an iterable of iterables, unpack it: max_of_each(*iter_of_iters)
	"""
	for iter_ in iters:
		yield max(iter_, key=key)


def min_of_each(*iters: Iterable[Any], key: Callable = lambda x: x) -> Iterator[Any]:
	"""
	Given multiple iterables, find the min-valued item in each iterable, and yield that item
	Like with min(), a custom function to apply on each may be supplied with the `key' argument
	To use with an iterable of iterables, unpack it: min_of_each(*iter_of_iters)
	"""
	for iter_ in iters:
		yield min(iter_, key=key)


def bifurcate_gff(entries: Iterable['GFF'], mask: Iterable[str]) -> tuple[set[str], set[str]]:
	"""
	Given a list of GFF entries and a list of entry IDs, split the IDs of all fetures in `entries'
	into two sets, such that one set includes the entries in `mask', along with all their lineages,
	and the other set excludes the entries in `mask', and as much of their lineage as possible
	"""
	#The two sets to be returned
	inc: set[str] = set()
	exc: set[str] = set()
	
	#Keys: IDs of features which have children; values: IDs of a given feature's children
	children: dict[str:set[str]] = defaultdict(set)
	#Keys: IDs of child features; values: ID of a given feature's parent
	parents: dict[str:str] = dict()
	
	
	#Populate `exc' (with all IDs in `entries', for now) and `children'
	for entry in entries:
		assert "ID" in entry.attrs
		exc.add(entry.attrs["ID"])
		if "Parent" in entry.attrs:
			children[ entry.attrs["Parent"] ].add( entry.attrs["ID"] )
			parents[ entry.attrs["ID"] ] = entry.attrs["Parent"]
	
	
	def descendants(id_) -> set[str]:
		"""
		Get all descendants of a given feature
		"""
		out = set()
		
		for child in children[id_]:
			#Add immediate children
			out.add(child)
			#Recursively add further descendants
			out.update(descendants(child))
		
		return out
	
	
	def singular_ancestors(id_) -> set[str]:
		"""
		Get the ancestors of a given feature, up to the oldest ancestor which does not have any
		additional children
		"""
		out = set()
		
		while id_ in parents and len(children[ parents[id_] ]) == 1:
			id_ = parents[id_]
			out.add(id_)
		
		return out
	
	
	def ancestors(id_) -> set[str]:
		"""
		Get all ancestors of a given feature
		"""
		out = set()
		
		while id_ in parents:
			id_ = parents[id_]
			out.add(id_)
		
		return out
	
	
	def lineage(id_) -> set[str]:
		"""
		Get the entire lineage on which a given feature lies, i.e. a set containing the (IDs of)
		the feature, all its descendants, and all its ancestors
		"""
		out: set[str] = descendants(id_) | ancestors(id_)
		out.add(id_)
		return out
	
	
	def delete(id_):
		"""
		Delete a feature from `exc' and unlink it from its parent & children
		"""
		exc.remove(id_)
		#Unlink this feature from its children
		if id_ in children:
			for child in children[id_]:
				del parents[child]
			del children[id_]
		#Unlink this feature from its parent
		if id_ in parents:
			children[ parents[id_] ].remove(id_)
			del parents[id_]
	
	
	#Move features from `exc' to `inc'
	for id_ in mask:
		#If the feature was moved already
		if id_ not in exc:
			continue
		
		
		#Move children
		for child in descendants(id_):
			#Add feature to `inc'
			inc.add(child)
			#Remove from `exc' and unlink
			delete(child)
		
		
		#Get singular and all ancestors
		sanc = singular_ancestors(id_)
		anc = ancestors(id_)
		
		#Move feature
		inc.add(id_)
		delete(id_)
		
		#Move singular ancestors, copy non-singular ancestors
		for ancestor in anc:
			inc.add(ancestor)
			if ancestor in sanc:
				delete(ancestor)
	
	
	return inc, exc

