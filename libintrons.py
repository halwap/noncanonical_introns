###################################################################################################
#   IMPORTS & SETUP
###################################################################################################

#Standard library 
#Type hints
from __future__ import annotations
from typing import *
#Efficient lambdas
from operator import attrgetter
#Convenient iteration
from itertools import pairwise, batched
#Diagnostics
from warnings import filterwarnings
#Statistics
from statistics import mean

#Third-party imports
#Retrieving seqs from the '-' strand
from Bio.Seq import reverse_complement
#Array processing
from numpy import array, append, empty, argmax

#Local imports
from formats import *

#Hide warnings
filterwarnings("ignore", category=UserWarning)



###################################################################################################
#   CLASSES
###################################################################################################

class GenomicSequence:
	scaffold: str
	start:    int
	end:      int
	strand:   Strand
	seq:      str|None
	
	
	def __init__(
		self,
		scaffold: str,
		start:    int,
		end:      int,
		*,
		strand:   Strand,
		seq:      str|None = None
	):
		self.scaffold = scaffold
		
		
		#Validate coordinates
		assert start < end
		
		self.start = start
		self.end   = end
		
		
		#Validate strand
		assert strand in "+-."
		
		self.strand = strand
		
		
		#Validate sequence length vis-a-vis passed coords
		if seq:
			assert len(seq) == end - start
		
		self.seq = seq
	
	
	def valid_len(self):
		"""
		Check if the `start' and `end' attributes "agree" with the `seq' attribute about the length
		of the sequence
		"""
		assert self.seq
		return len(self.seq) == self.end - self.start
	
	
	def __bool__(self) -> bool:
		return True

	
	
	#Report type and position of sequence

	def __repr__(self):
		#Example representation: '[GenomicSequence @ scaffold_0(+): 1153-2093]'
		return f"[{type(self).__name__} @ {self.scaffold}({self.strand}): {self.start}-{self.end}]"

	def __str__(self):
		#Use that same formatting for string conversion
		return self.__repr__()
	
	
	#Shortcuts for string-like ops on sequences
	
	def __len__(self):
		"""
		Get length of sequence
		Assumes the sequence is present
		"""
		return len(self.seq)
	

	def __getitem__(self, index: int|slice):
		"""
		Index and slice nucleotides of sequence
		"""
		return self.seq[index]
	

	def __add__(self, other: Self|str) -> str:
		"""
		Concatenate a GenomicSequence's sequence with another one's, or with a str
		"""
		if isinstance(other, str):
			return self.seq + other
		elif isinstance(other, GenomicSequence):
			return self.seq + other.seq
		else:
			raise TypeError(f"Can't concatenate {type(self).__name__} & {type(other).__name__}")
	

	def __radd__(self, other: Self|str) -> str:
		"""
		Concatenate a GenomicSequence's sequence with another one's, or with a str
		"""
		if isinstance(other, str):
			return other + self.seq
		elif isinstance(other, GenomicSequence):
			return other.seq + self.seq 
		else:
			raise TypeError(f"Can't concatenate {type(other).__name__} & {type(self).__name__}")
	
	
	#Methods for checking the relative position of two GenomicSequence objects
	def comparable(self, other: Self) -> bool:
		"""
		Check whether two GenomicSequence objects are from the same strand & scaffold
		This is a necessary conditions for further comparisons to even be considered
		"""
		return self.strand == other.strand and self.scaffold == other.scaffold
	
	def __lt__(self, other: Self) -> bool:
		"""
		Check whether a GenomicSequence is entirely to the right of `self', without flushed borders
		Comparison is only relevant if both objects are from the same scaffold & strand
		"""
		return self.comparable(other) and self.end < other.start
	
	def __le__(self, other: Self) -> bool:
		"""
		Check whether a GenomicSequence is entirely to the right of `self', with flushed borders
		Comparison is only relevant if both objects are from the same scaffold & strand
		"""
		return self.comparable(other) and self.end == other.start
	
	def __contains__(self, other: Self) -> bool:
		"""
		Check whether `self' completely overlaps another GenomicSequence
		Comparison is only relevant if both objects are from the same scaffold & strand
		"""
		return self.comparable(other) and self.start <= other.start and other.end <= self.end
	
	def __eq__(self, other: Self) -> bool:
		"""
		Check whether two GenomicSequence objects have the same coordinates, scaffold & strand
		"""
		return self.comparable(other) and self.start == other.start and other.end == self.end
	

	
	
	def to_gff(self, type_: str = "sequence_feature", **attrs: str) -> GFF:
		"""
		Convert a GenomicSequence object to a GFF entry
		
		Params:
			type_: str
				Type of feature
				The default type is "sequence_feature"
			**attrs: str
				Attributes of the sequence to report, passed as keywords
				All values should be strings
		"""
		return GFF(
					seqid  = self.scaffold,
					source = "libintrons",
					type_  = type_,
					start  = self.start+1,
					end    = self.end,
					score  = None,
					strand = self.strand,
					phase  = None,
					attrs  = attrs
			)



class Gene(GenomicSequence):
	scaffold:    str
	start:       int
	end:         int
	strand:      Strand
	seq:         str|None
	transcript:  Transcript
	exons:       list[Exon]
	introns:     list[Intron]
	name:        str|None
	
	
	def __init__(self,
		scaffold:   str,
		start:      int,
		end:        int,
		*,
		strand:     Strand       = '.',
		seq:        str|None     = None,
		transcript: Transcript   = None,
		exons:      list[Exon]   = [],
		introns:    list[Intron] = [],
		name:       str|None     = None
	):
		super().__init__(scaffold, start, end, strand=strand, seq=seq)
		self.exons           = exons
		self.introns         = introns
		self.name            = name
		#Autogenerate a transcript if it was not passed
		self.transcript      = transcript or Transcript(scaffold, start, end, strand=strand)
	
	
	def fixup_exons(self):
		"""
		Appropriately sort and link the exons in self.exons
		"""
		#If there's no exons, or only 1 exon, there's nothing to do
		if len(self.exons) < 2:
			return
		
		#Sort exons by start position
		#If the gene is on the positive strand, this is also sorting them by their within-gene
		#order; if negative, however, they are sorted in reversed within-gene order
		self.exons.sort(key=attrgetter("start"))
		
		#Validate coordinates
		assert self.valid_coords()
		
		
		#Link exons
		for exon1, exon2 in pairwise(self.exons):
			exon1.next_exon = exon2
			exon2.prev_exon = exon1
		
		
		#Validate linkage
		assert self.valid_links()
	
	
	def add_seqs(self, genome: dict[str,str], which: str = "gtei"):
		"""
		Given a deserialized fasta file, build the sequences of the gene/exons/introns
		`which' determines which specific features should obtain new sequences
		By default, all 3 features (gene, exons & introns) will get new seqs
		"""
		#Define lambda for obtaining the proper orientation of a sequence
		if self.strand == '-':
			orient: Callable[[str],str] = lambda seq: reverse_complement(seq)
		else:
			orient: Callable[[str],str] = lambda seq: seq
		
		
		#Validate coordinates
		assert self.valid_coords()
		
		
		#Get sequence of appropriate scaffold
		#This requires the gene's scaffold to be present in the genome
		assert self.scaffold in genome
		scaffold_seq: str = genome[self.scaffold]
		
		
		#Push sequence to gene
		if 'g' in which:
			self.seq = orient(scaffold_seq[self.start:self.end])
			assert self.valid_len()
		
		
		#Push sequence to transcript
		if 't' in which:
			self.transcript.seq = ""
			for exon in self.order_exons():
				self.transcript.seq += orient(scaffold_seq[exon.start:exon.end])
		
		
		#Push sequences to exons
		if 'e' in which:
			for exon in self.exons:
				exon.seq = orient(scaffold_seq[exon.start:exon.end])
				assert exon.valid_len()
		
		
		#Push sequences to introns
		if 'i' in which:
			for intron in self.introns:
				intron.seq = orient(scaffold_seq[intron.start:intron.end])
				assert intron.valid_len()

		
		#Validate created sequences
		assert self.valid_seqs()
	
	
	def add_introns(self):
		"""
		Populate self.introns by inserting an intron between every pair of exons
		"""
		#Check that there are exons at all
		assert len(self.exons) > 0
		
		#Validate coordinates & linkage
		assert self.valid_coords()
		assert self.valid_links()
		

		#If there's only one exon, there are no introns to add; exit early
		if len(self.exons) == 1:
			return
		
		
		self.introns = []
		
		
		for exon1, exon2 in pairwise(self.exons):

			#Get coordinates of the intron
			start = exon1.end
			end   = exon2.start
			
			#Instantiate intron & add it to the list
			self.introns.append(
				Intron(
					self.scaffold, start, end, strand = self.strand,
					gene = self, prev_exon = exon1, next_exon = exon2
				)
			)
			
			#Validate linkage & coords between the new intron & exons
			#There's no linkage check between exon1 & exon2, but that's ok, since the
			#self.valid_links() asserts above already did that
			assert exon1 <= self.introns[-1] <= exon2
			assert exon1 @ self.introns[-1] @ exon2
		
		
		#Check exon & intron count
		assert len(self.exons) == len(self.introns) + 1
		
		#Validate coords & linkage, this time with introns
		assert self.valid_coords()
		assert self.valid_links()
	
	
	def add_intron_variants(self, min_exon_len: int):
		"""
		Build a list of variants for each intron in the gene
		"""
		for intron in self.introns:
			intron.add_variants(min_exon_len)
	
	
	def rectify_introns(self):
		"""
		Rectify all introns in the gene by selecting the highest-scored variant for each one.
		Assumes introns have been scored already.
		"""
		#This method should do nothing if the gene has no introns
		if len(self.introns) == 0:
			return
		
		#Check that all introns have been scored
		assert all( intron.scored() for intron in self.introns )
		
		#Validate coords & linkage
		assert self.valid_coords()
		assert self.valid_links()
		
		
		for intron_idx, intron in enumerate(self.introns):
			#Skip introns without variants
			if intron.variants:
				
				#Check that all variants of this intron have been scored
				assert all( variant.scored() for variant in intron.variants )
				
				
				#Get index of variant with highest score
				best_variant_idx = argmax([ variant.unif_score for variant in intron.variants ])
				
				
				#If the highest-scored variant supersedes the base intron, supplant it
				if intron.variants[best_variant_idx].unif_score > intron.unif_score:
					self.swap_intron_and_variant(intron_idx, best_variant_idx)
		
		
		#With the best variant selected for each intron, the exons need also be switched, so that
		#all exon-intron boundaries match
		
		#Take care of terminal exons
		self.exons[0] = self.introns[0].prev_exon
		self.exons[-1] = self.introns[-1].next_exon

		
		#The left and right boundaries of each medial exon were nudged independently, and so may
		#need to be reconciled
		for exon_idx, (intron1, intron2) in enumerate(pairwise(self.introns), 1):
			if not intron1 < intron2:
				raise ValueError("Vanishing exon found, try increasing minimum exon length")
			
			
			#If the introns disagree on the exon between them, reconciliation is necessary
			if intron1.next_exon != intron2.prev_exon:
				
				#Reconciliation is done by having the two introns dictate the new start & end
				#of the exon, respectively
				intron1.next_exon.end = intron2.prev_exon.end
				
				
				#Link new exon & secons intron (first exon & intron are already linked)
				intron2.prev_exon = intron1.next_exon
				intron1.next_exon.next_intron = intron2

			
			#At this point, the introns should agree on the exon between them, either naturally, or
			#thanks to reconciliation
			#Note that `intron1.next_exon' and `intron2.prev_exon' might not necessarily be the
			#same object, they just need to have identical coordinates
			assert intron1.next_exon == intron2.prev_exon
			
			#Replace exon in list
			self.exons[exon_idx] = intron1.next_exon
			

			#Check linkage & coordinates between the exon and the introns
			assert intron1 @ self.exons[exon_idx]  @ intron2
			assert intron1 <= self.exons[exon_idx] <= intron2
	
	
	def swap_intron_and_variant(self, intron_idx: int, variant_idx: int):
		"""
		Take self.introns[intron_idx] and self.introns[intron_idx].variants[variant_idx], and swap
		them around.
		"""
		#Validate input
		assert intron_idx < len(self.introns)
		
		#Get the base intron & the variant
		intron: Intron = self.introns[intron_idx]
		variant: Intron = intron.variants[variant_idx]
		
		#Validate input
		assert variant_idx < len(intron.variants)
		

		#Replace `variant' with `intron' in both lists
		intron.variants[variant_idx] = intron
		self.introns[intron_idx] = variant
		
		#Repossess the list of variants
		variant.variants = intron.variants
		del intron.variants
	
	
	def to_gff(self) -> list[GFF]:
		"""
		Convert a Gene and all its children features to a list of GFF instances
		"""
		#Validate coordinates
		assert self.valid_coords()
		
		#Check that both the gene and the transcript have names to use as their IDs in the GFF
		assert self.name
		assert self.transcript.name
		
		
		#Will hold converted entries
		entries: list[GFF] = []
		
		#Detemine whether all introns of the gene have been scored
		introns_scored = self.introns and all(intron.scored() for intron in self.introns)
		
		
		#Gene
		attrs = { "ID": self.name }
		#If all introns have been scored, report the average score in attributes
		if introns_scored:
			attrs["intron_score"] = str( mean( map( attrgetter("unif_score" ), self.introns ) ))
		
		entries.append( super().to_gff("gene", **attrs) )
		
		
		#Transcript
		mrna_id = self.transcript.name
		entries.append( self.transcript.to_gff("mRNA", Parent=self.name, ID=mrna_id) )
		
		
		#Exons
		for n, exon in enumerate(self.exons, 1):
			entries.append( exon.to_gff("exon", Parent=mrna_id, ID=f"{mrna_id}.exon{n}") )
		
		
		#Introns
		for n, intron in enumerate(self.introns, 1):
			attrs = { "Parent": mrna_id, "ID": f"{mrna_id}.intron{n}" }
			
			#Report splice site if available
			if intron.splice_site:
				attrs["splice_site"] = intron.splice_site
			
			#Report number of variants, if adding variants was requested
			if intron.variants is not None:
				attrs["variant_cnt"] = str(len(intron.variants))
			
			#Report scores if available
			if introns_scored:
				attrs["conv_score"] = str(intron.c_score)
				attrs["nonconv_score"] = str(intron.nc_score)
			
			entries.append( intron.to_gff("intron", **attrs) )
		
		
		return entries
	
	
	def frags(self) -> list[Exon|Intron]:
		"""
		Returns a list of interpersed exons & introns
		"""
		out = []
		for exon, intron in zip(self.exons, self.introns):
			out.append(exon)
			out.append(intron)
		out.append(self.exons[-1])
		return out
	
	
	def order_exons(self) -> Iterable[Exon]:
		"""
		Returns an iterable of exons in within-gene order
		"""
		return reversed(self.exons) if self.strand == '-' else self.exons
	
	
	def order_introns(self) -> Iterable[Intron]:
		"""
		Returns an iterable of introns in within-gene order
		"""
		return reversed(self.introns) if self.strand == '-' else self.introns
	
	
	def order_frags(self) -> Iterable[Exon|Intron]:
		"""
		Returns an iterable of interspersed exons & introns in within-gene order
		"""
		return reversed(self.frags()) if self.strand == '-' else self.frags()
	
	
	def valid_coords(self) -> bool:
		"""
		Check that the coordinates of the gene and all its sub-features are in order:
		1. The transcript, all exons & all introns are within the bounds of the gene
		2. Exons are properly ordered and non-contiguous on the scaffold
		3. Introns are properly ordered and non-contiguous on the scaffold
		4. Exons and introns are properly ordered and contiguous on the scaffold
		Checks which do not apply are skipped
		"""
		#Check that the transcript is within the bounds of the gene
		if self.transcript not in self:
			return False

		
		#If there are exons to check
		if self.exons:
			#Check that all exons are within bounds of the gene
			if any( exon not in self for exon in self.exons ):
				return False

			#Check that all exons are properly ordered, and there are gaps between each one
			if any( not exon1 < exon2 for exon1, exon2 in pairwise(self.exons) ):
				return False
		
		
		#If there are introns to check
		if self.introns:
			#Check that introns exons are within bounds of the gene
			if any( intron not in self for intron in self.introns ):
				return False
			
			#Check that all introns are properly ordered, and there are gaps between each one
			if any( not intron1 < intron2 for intron1, intron2 in pairwise(self.introns) ):
				return False
		
		
		#If there are both exons & introns
		if self.exons and self.introns:
			#Check that consecutive exon-intron pairs are properly ordered, and flush
			if any( not frag1 <= frag2 for frag1, frag2 in pairwise(self.frags()) ):
				return False
		
		
		#If all relevant checks passed, the coordinates are valid
		return True
	
	
	def valid_links(self) -> bool:
		"""
		Check that all exons & introns are properly linked
		"""
		#If there are exons to check
		if self.exons:
			#Check that each pair of exons is mutually linked
			if any( not exon1 @ exon2 for exon1, exon2 in pairwise(self.exons) ):
				return False
		
		
		#If there are both exons & introns
		if self.exons and self.introns:
			#Check that consecutive exon-intron pairs are mutually linked
			if any( not frag1 @ frag2 for frag1, frag2 in pairwise(self.frags()) ):
				return False
		
		
		#If all relevant checks passed, the links are valid
		return True
	
	
	def valid_seqs(self) -> bool:
		"""
		Check that the sequences of the gene and all its sub-features are in order:
		1. Either all or no exons have sequences
		2. Either all or no introns have sequences
		3. The transcript's sequence matches the combined sequence of the exons
		4. The gene's sequence matches the combined sequence of the exons and introns
		Checks which do not apply are skipped
		"""
		#Check which specific features actually have sequences
		gene_seq:       bool = bool(self.seq)
		transcript_seq: bool = bool(self.transcript.seq)
		exons_seq:      bool = any(exon.seq for exon in self.exons)
		introns_seq:    bool = any(intron.seq for intron in self.introns)
		
		
		#If at least one exon has a sequence, check that all exons do
		if exons_seq:
			if not all(exon.seq for exon in self.exons):
				return False

		#Likewise for introns
		if introns_seq:
			if not all(intron.seq for intron in self.introns):
				return False
		

		#If both the transcript & exons have sequences, check that the exons' sequnences match the
		#transcript's
		if transcript_seq and exons_seq \
		and concat_genseq(*self.order_exons()) != self.transcript.seq:
			return False
		
		#Likewise, if the gene, the introns & exons all have sequences, check that they match
		if gene_seq and exons_seq and introns_seq \
		and concat_genseq(*self.order_frags()) not in self.seq:
			return False
		
		
		#If all checks passed, the sequences are valid
		return True



class Transcript(GenomicSequence):
	scaffold: str
	start:    int
	end:      int
	strand:   Strand
	seq:      str|None
	name:     str|None
	
	
	def __init__(
		self,
		scaffold: str,
		start:    int,
		end:      int,
		*,
		strand:   Strand   = '.',
		seq:      str|None = None,
		name:     str|None = None
	):
		super().__init__(scaffold, start, end, strand=strand, seq=seq)
		self.name = name


class Exon(GenomicSequence):
	scaffold:    str
	start:       int
	end:         int
	strand:      Strand
	seq:         str|None
	gene:        Gene|None
	prev_exon:   Exon|None
	next_exon:   Exon|None
	prev_intron: Intron|None
	next_intron: Intron|None
	
	
	def __init__(
		self,
		scaffold:    str,
		start:       int,
		end:         int,
		*,
		strand:      Strand      = '.',
		seq:         str|None    = None,
		gene:        Gene|None   = None,
		prev_exon:   Exon|None   = None,
		next_exon:   Exon|None   = None,
		prev_intron: Intron|None = None,
		next_intron: Intron|None = None
	):
		super().__init__(scaffold, start, end, strand=strand, seq=seq)
		
		#If preceeding intron was passed
		if prev_intron:
			#Validate input
			assert prev_intron <= self
			
			#Link exon to intron
			prev_intron.next_exon = self
		
		#If following intron was passed
		if next_intron:
			#Validate input
			assert self <= next_intron
		
			#Link exon to intron
			next_intron.prev_exon = self
		
		#Prospectively link introns to exon
		self.prev_intron = prev_intron
		self.next_intron = next_intron
		
		
		#If gene was passed
		if gene:
			#Validate input
			assert self in gene
		
		#Prospectiely link gene to exon
		self.gene = gene
		
		
		#If preceeding exon was passed
		if prev_exon:
			#Validate input
			assert prev_exon < self
			
			#Link this exon to that one
			prev_exon.next_exon = self
		
		#If following exon was passed
		if next_exon:
			#Validate input
			assert self < next_intron

			#Link this exon to that one
			next_exon.prev_exon = self
		
		#Prospectively links exons to this one
		self.prev_exon = prev_exon
		self.next_exon = next_exon
	
	
	def __rmatmul__(self, other: Exon|Intron|None) -> Self|None:
		"""
		Check if this exon and a preceding exon or intron are linked together correctly, to be
		used via the '@' operator, e.g. "ft1 @ ft2"
		"ft1 @ ft2" will return ft2 if the linkage was found to be correct, and None if incorrect
		The result can then be converted to a bool, since GenomicSequence (and subclass) objects
		are truthy by default, while None is falsy
		Also, the '@' operator can be chained: "ft1 @ ft2 @ ft3" is equivalent to
		"ft1 @ ft2 and ft2 @ ft3", and will check the linkage between ft1-ft2 and ft2-ft3 as
		intended (with ft3 or None as the result)
		Note that ft1-ft3 isn't checked in this case; crucially, this means that something like
		"exon1 @ intron @ exon2" might also need a separate "exon1 @ exon2"
		"""
		#If None was passed, this comparison is further down a chain of '@' operators, and an
		#earlier comparison failed; pass the result further to short circuit
		if other is None:
			return None
		
		if not isinstance(other, (Exon, Intron)):
			raise TypeError(f"Can't test links for {type(other).__name__} & {type(self).__name__}")
		
		
		#Get the features pointed to
		prev = self.prev_exon if isinstance(other, Exon) else self.prev_intron
		next_ = other.next_exon
		
		#Both features must link to each other for the linkage to be correct
		if prev and next_ and prev == other and next_ == self:
			return self
		else:
			return None
	
	
	#Thin wrapper around __rmatmul__, needed for it to work correctly
	def __matmul__(self, other):
		return other.__rmatmul__(self)
	
	
	def former_exon(self) -> Exon|None:
		"""
		Returns the exon preceding this exon in within-gene order
		"""
		return self.prev_exon if self.strand == '+' else self.next_exon
	
	
	def latter_exon(self) -> Exon|None:
		"""
		Returns the exon following this exon in within-gene order
		"""
		return self.next_exon if self.strand == '+' else self.prev_exon
	
	
	def former_intron(self) -> Intron|None:
		"""
		Returns the intron preceding this intron in within-gene order
		"""
		return self.prev_intron if self.strand == '+' else self.next_intron
	
	
	def latter_intron(self) -> Intron|None:
		"""
		Returns the intron following this exon in within-gene order
		"""
		return self.next_intron if self.strand == '+' else self.prev_intron


class Intron(GenomicSequence):
	scaffold:    str
	start:       int
	end:         int
	strand:      Strand
	seq:         str|None
	gene:        Gene|None
	prev_exon:   Exon|None
	next_exon:   Exon|None
	variants:    list[Intron]|None
	c_score:     float|None
	nc_score:    float|None
	unif_score:  float|None
	splice_site: str|None
	traits:      dict[str,bool|float]


	def __init__(
		self,
		scaffold:  str,
		start:     int,
		end:       int,
		*,
		strand:    Strand    = '.',
		seq:       str|None  = None,
		gene:      Gene|None = None,
		prev_exon: Exon|None = None,
		next_exon: Exon|None = None,
	):
		super().__init__(scaffold, start, end, strand=strand, seq=seq)
		
		
		#Validate `gene'
		if gene:
			assert self in gene
		self.gene = gene
		
		
		#If preceeding exon was passed
		if prev_exon:
			#Validate input
			assert prev_exon <= self
			
			#Link this exon to that one
			prev_exon.next_intron = self
		
		#If following exon was passed
		if next_exon:
			#Validate input
			assert self <= next_exon
			#Link this exon to that one
			next_exon.prev_intron = self
		
		#Prospectively links exons to this one
		self.prev_exon = prev_exon
		self.next_exon = next_exon
		
		
		#Attributes whose specific values are to be computed later
		self.variants    = None
		self.c_score     = None
		self.nc_score    = None
		self.unif_score  = None
		self.splice_site = None
		self.traits      = {}
	
	
	def __rmatmul__(self, other: Exon|None) -> Self|None:
		"""
		Check if this intron and a preceding exon are linked together correctly, to be
		used via the '@' operator, e.g. "ft1 @ ft2"
		"ft1 @ ft2" will return ft2 if the linkage was found to be correct, and None if incorrect
		The result can then be converted to a bool, since GenomicSequence (and subclass) objects
		are truthy by default, while None is falsy
		Also, the '@' operator can be chained: "ft1 @ ft2 @ ft3" is equivalent to
		"ft1 @ ft2 and ft2 @ ft3", and will check the linkage between ft1-ft2 and ft2-ft3 as
		intended (with ft3 or None as the result)
		Note that ft1-ft3 isn't checked in this case; crucially, this means that something like
		"exon1 @ intron @ exon2" might also need a separate "exon1 @ exon2"
		"""
		#If None was passed, this comparison is further down a chain of '@' operators, and an
		#earlier comparison failed; pass the result further to short circuit
		if other is None:
			return None
		
		#Refuse comparisons with incompatible types
		if not isinstance(other, Exon):
			raise TypeError(f"Can't test links for {type(other).__name__} & {type(self).__name__}")
		
		
		#Get the features pointed to
		prev = self.prev_exon
		next_ = other.next_intron
		
		#Both features must link to each other for the linkage to be correct
		if prev and next_ and prev == other and next_ == self:
			return self
		else:
			return None
	
	
	def scored(self) -> bool:
		"""
		Check if this intron has been scored
		"""
		return None not in (self.c_score, self.nc_score, self.unif_score)
	
	
	def has_exons(self) -> bool:
		"""
		Check if this intron has surrounding exons
		"""
		return bool(self.prev_exon and self.next_exon)
	
	
	def has_seqs(self) -> bool:
		"""
		Check if this intron and its surrounding exons all have sequences
		Implicitly also tests whether the surrounding exons are at all present
		"""
		return self.has_exons() and bool(self.seq and self.prev_exon.seq and self.next_exon.seq)
	
	
	def former_exon(self) -> Exon|None:
		"""
		Returns the exon preceding this intron in within-gene order
		"""
		return self.prev_exon if self.strand == '+' else self.next_exon
	
	
	def latter_exon(self) -> Exon|None:
		"""
		Returns the exon following this intron in within-gene order
		"""
		return self.next_exon if self.strand == '+' else self.prev_exon
	
	
	def add_variants(self, min_exon_len: int):
		"""
		Find possible variants of an intron by examining the nucleotides at exon-intron boundaries
		and add all of them to self.variants
		"""
		#Check that the intron and surrounding exons have sequences
		assert self.has_seqs()
		
		
		#Get surrounding exons as local variables, to save on typing
		prev_exon = self.prev_exon
		next_exon = self.next_exon
		
		
		#Check coordinates & linkage
		assert prev_exon <= self <= next_exon
		assert prev_exon @ self @ next_exon and prev_exon @ next_exon
		
		
		#Project what the intron & exons would look like before & after splicing, and validate
		#these sequences vis-a-vis gene & transcript sequences
		assert ( pre_mrna_seq := self.former_exon() + self + self.latter_exon() ) in self.gene.seq
		assert ( mrna_seq := self.former_exon() + self.latter_exon() ) in self.gene.transcript.seq
		#An additional purpose of these asserts is to initialize `pre_mrna_seq' & `mrna_seq'
		#The actual purpose of these asserts is to get the values of `pre_mrna_seq' & `mrna_seq'
		#These two variables are used exclusively in asserts, so they don't need to be initialized
		#if asserts are disabled - putting the assignments in asserts achieves that
		

		#prev_exon & next_exon are previous & next in the coordinate sense, which, as long as
		#the gene in on the positive strand, is also their within-gene order
		#If the gene is on the negative strand, this is opposite to within-gene order, however,
		#which would cause variant detection to yield incorrect results
		#To remedy this, if the gene is on the negative strand, the sequences of the intron & exons
		#are temporarily reversed
		if self.strand == '-':
			#Backup correctly-oriented sequences
			self.seq_      = self.seq
			prev_exon.seq_ = prev_exon.seq
			next_exon.seq_ = next_exon.seq
			#Reverse the sequences
			self.seq      = self.seq[::-1]
			prev_exon.seq = prev_exon.seq[::-1]
			next_exon.seq = next_exon.seq[::-1]
		
		
		#What shifts are possible depends on the nucleotides at the ends of the
		#intron and the previous exon (for shifts towards the 5' end), or at
		#the starts of the intron and the next exon (for shifts towards 3')
		#
		#The size of the largest possible shift towards 5' is equal to the
		#length of the longest suffix shared between the intron and the
		#previous exon
		#
		#The size of the largest possible shift towards 3' is equal to the
		#length of the longest prefix shared between the intron and the
		#next exon
		#
		#Consider the input intron looks like this: TTCAG|TGGTC...TGCAG|TGACT
		#
		#The longest shared suffix between the intron and the previous exon is
		#CAG (length 3), so the intron can be shifted towards 5' by up to 3
		#nucleotides
		#
		#The longest shared prefix between the intron and the next exon is
		#TG (length 2), so the intron can be shifted towards 3' by up to 2
		#nucleotides
		#
		#Thus the possible variants are as follows:
		#Shift of -3: TT|CAGTGGTC...TG|CAGTGACT
		#Shift of -2: TTC|AGTGGTC...TGC|AGTGACT
		#Shift of -1: TTCA|GTGGTC...TGCA|GTGACT
		#No shift:    TTCAG|TGGTC...TGCAG|TGACT
		#Shift of +1: TTCAGT|GGTC...TGCAGT|GACT
		#Shift of +2: TTCAGTG|GTC...TGCAGTG|ACT
		#(Negative shift values mean a shift towards the 5' end,
		#positive shift values mean a shift towards the 3' end)
		#The unshifted variant is already available as `self' and is not to be
		#included in self.variants
		
		
		#Get the shared prefix and suffix's length to determine the range of shift values
		
		#Intron-prev exon shared suffix
		suff_len = 0
		#`min_exon_len' is used to put a lower bound on how short exons can become upon a shift
		left_bound = min(len(self), len(prev_exon) - min_exon_len)
		while suff_len < left_bound and self[-1 - suff_len] == prev_exon[-1 - suff_len]:
			suff_len += 1
		
		#Intron-next exon shared prefix
		pref_len = 0
		right_bound = min(len(self), len(next_exon) - min_exon_len)
		#This check also performs a length check for both sequences
		while pref_len < right_bound and self[pref_len] == next_exon[pref_len]:
			pref_len += 1
		
		
		#Instantiate/reset list of variants
		self.variants = []
		
		
		#Iterate over all possible shift values
		for shift in range(-suff_len, pref_len+1):
			#Skip a shift value of 0, since it means no shift
			if shift == 0:
				continue
			
			
			#Get coordinates of intron & exon after shifts
			new_prev_exon_start = prev_exon.start
			new_prev_exon_end   = prev_exon.end + shift

			new_intron_start    = self.start + shift
			new_intron_end      = self.end   + shift

			new_next_exon_start = next_exon.start + shift
			new_next_exon_end   = next_exon.end
			
			
			#Instantiate the variant's flanking exons
			new_prev_exon = Exon(
				self.scaffold, new_prev_exon_start, new_prev_exon_end, strand = self.strand,
				gene = self.gene
			)
			new_next_exon = Exon(
				self.scaffold, new_next_exon_start, new_next_exon_end, strand = self.strand,
				gene = self.gene, prev_exon = new_prev_exon
			)
			
			#Instantiate the variant
			new_intron = Intron(
				self.scaffold, new_intron_start, new_intron_end, strand = self.strand,
				gene = self.gene, prev_exon = new_prev_exon, next_exon = new_next_exon
			)
			
			#Passing the prev/next exons to the constructors mean that the introns & exons are
			#automatically linked
			#Note that the variant exons are linked to each other, but not to any further preceding
			#or following exons
			
			
			#Build the new exons' & intron's sequences
			#Shift towards 5'
			if shift < 0:
				new_prev_exon.seq = prev_exon[:shift]
				new_intron.seq    = prev_exon[shift:] + self[:shift]
				new_next_exon.seq =                     self[shift:] + next_exon
			
			#Shift towards 3'
			else:
				new_prev_exon.seq = prev_exon + self[:shift]
				new_intron.seq    =             self[shift:] + next_exon[:shift]
				new_next_exon.seq =                            next_exon[shift:]
			
			
			#Check sequence lengths
			assert new_prev_exon.valid_len()
			assert new_intron.valid_len()
			assert new_next_exon.valid_len()
			
			
			#Undo the new exons' & intron's sequences being reversed if needed
			if self.strand == '-':
				new_intron.seq    = new_intron.seq[::-1]
				new_prev_exon.seq = new_prev_exon.seq[::-1]
				new_next_exon.seq = new_next_exon.seq[::-1]
			
			
			#Check whether the shift leaves the unspliced sequence untouched
			assert pre_mrna_seq == new_intron.former_exon() + new_intron + new_intron.latter_exon()
			assert mrna_seq == new_intron.former_exon() + new_intron.latter_exon()
			
			
			#Validate coords & linkage
			assert new_prev_exon <= new_intron <= new_next_exon
			assert new_prev_exon @ new_intron @ new_next_exon and new_prev_exon @ new_next_exon
			
			
			#Add variant to list
			self.variants.append(new_intron)
		
		
		#Unreverse the sequences if needed
		if self.strand == '-':
			#Restore backups
			self.seq      = self.seq_
			prev_exon.seq = prev_exon.seq_
			next_exon.seq = next_exon.seq_
			#Remove backups
			del self.seq_
			del prev_exon.seq_
			del next_exon.seq_
	
	
	def add_splice_site(self):
		"""
		Compute the splice site an intron, and save it under the `splice_site' attribute
		Will be None if the intron is less than 4 nt long
		"""
		assert self.seq
		
		self.splice_site: str|None = self[:2] + self[-2:] if len(self) >= 4 else None
	
	
	
	def add_traits(self, weighted: bool = True, pairing_len: int = 21):
		"""
		Add attributes describing the intron's conventional and nonconventional traits
		If `weighted' is True, the pairing scores will be weighted
		`pairing_len' is the max number of nucleotides to examine during pairing score calculation
		"""
		#Check that the intron and surrounding exons have sequences
		assert self.has_seqs()
		
		#Check coordinates & linkage
		assert self.prev_exon <= self <= self.next_exon
		assert self.prev_exon @ self @ self.next_exon and self.prev_exon @ self.next_exon
		
		
		#Get surrounding exons as local variables, to save on typing
		#Note that the exons are taken in within-gene order
		prev_exon = self.former_exon()
		next_exon = self.latter_exon()
		
		
		#Re-instantiate trait dictionary
		#Trait values are bools (informing if a given trait is present), with the exception of
		#pairing scores, which are floats
		self.traits: dict[str,bool|float] = {}
		
		
		#Last nucleotide of previous exon is a pyrimidine
		self.traits["prev_exon_y"] = prev_exon[-1] in PYRIMIDINES
		
		#First nucleotide of intron is a purine
		self.traits["intron_r"]    = self[0]       in PURINES
		
		#Last nucleotide of intron is a pyrimidine
		self.traits["intron_y"]    = self[-1]      in PYRIMIDINES
		
		#First nucleotide of next exon is a purine
		self.traits["next_exon_r"] = next_exon[0]  in PURINES
		
		
		#CAG & CTG at appropriate positions
		self.traits["intron_cag"] = len(self) >= 6 and self[3:6]   == "CAG"
		self.traits["intron_ctg"] = len(self) >= 8 and self[-8:-5] == "CTG"
		self.traits["intron_cagctg"] = self.traits["intron_cag"] and self.traits["intron_ctg"]
		
		
		#Add splice site
		#The splice site is its own attribute, not a trait
		self.add_splice_site()
		
		#Intron has a conventional splice site
		self.traits["ss_is_conv"] = self.splice_site in CONV_SS
		
		
		
		#Determine pairing scores
		#The intron's start and end (plus a few nucleotides from both exons) are examined to find
		#the best contiguous run of pairing nucleotides
		#This is done in 3 different configurations:
		#
		#
		#1. With the intron's 3' end jutting out 1 nucleotide past the 5' end
		#
		#→ ┈┈──────────────────────────────────────────────┬───***─────────────────────────────┈┈ →
		#    ← 5' prev exon 3' →             ... AGTCAAGG T│AGTCAGGCTAGTC ...              ← 5'
		#  ┈┈────────────────────────────────────────────┬─┘                              intron
		#    ← 3' next exon 5' →             ... AGGTCCAA│C AAGCGTCTCGCGTC ...             ← 3'
		#← ┈┈────────────────────────────────────────────┴──────***────────────────────────────┈┈ ←
		#(The intron is folded in on itself, to project its pre-mRNA secondary structure;
		# the characteristic CAG & CTG are marked with asterisks)
		#
		#
		#2. With the 3' end jutting out 2 nucleotides past the 5' end - this way, if the
		#   characteristic CAG & CTG are at the appropriate positions, they're aligned and pair up:
		#
		#→ ┈┈──────────────────────────────────────────────┬───***─────────────────────────────┈┈ →
		#    ← 5' prev exon 3' →             ... AGTCAAG GT│AGTCAGGCTAGTC ...              ← 5'
		#  ┈┈───────────────────────────────────────────┬──┘                              intron
		#    ← 3' next exon 5' →             ... GGTCCAA│CA AGCGTCTCGCGTC ...              ← 3'
		#← ┈┈───────────────────────────────────────────┴──────***─────────────────────────────┈┈ ←
		#
		#
		#3. With the 3' end jutting out 3 nucleotides past the 5' end:
		#
		#→ ┈┈──────────────────────────────────────────────┬───***─────────────────────────────┈┈ →
		#    ← 5' prev exon 3'→              ... AGTCAA GGT│AGTCAGGCTAGTC ...              ← 5'
		#  ┈┈──────────────────────────────────────────┬───┘                              intron
		#    ← 3' next exon 5'→              ... GTCCAA│CAA GCGTCTCGCGTC ...               ← 3'
		#← ┈┈──────────────────────────────────────────┴──────***──────────────────────────────┈┈ ←
		#
		#
		#Separately for each configuration, runs of nucleotides which pair up (A-T, C-G, G-T) are
		#found. Each run is scored, and the pairing score for a given configuration is the score of
		#its highest-scoring run.
		#For unweighted pairing scores, the score of a run is simply its length.
		#For weighted scores, each pair of nucleotides is given a weight (C-G > A-T > G-T), and the
		#score of a run is the sum of the weights of its constituent nucleotides.
		

		#Get sequences for the above described alignments
		#
		#Sequence of intron forward, plus the last 5 nucleotides of previous exon
		#The previous exon is padded with Ns, if needed, to ensure it is in fact 5 nucleotides
		#The intron is limited to be `pairing_len'-5 nucleotides, so that this whole sequence is
		#`pairing_len' long
		#Limiting this sequence's length serves to apply `pairing_len''s length limit
		forward: str = prev_exon[-5:].rjust(5, 'N') + self[:pairing_len-5]
		#In the above illustrated example, this sequence would be
		#"AAGGTAGTCAGGCTAGTC..." (truncated to be `pairing_len' long)
		
		
		#Get sequence of intron, plus a bit of the following exon, reversed
		#The number of nucleotides from the following exon are such that, upon alignment, the 3'
		#end of the intron juts out past the 5' end 
		#3' end juts out by 1 nucleotide past 5' end
		reverse_1: str = self + next_exon[:4].ljust(4, 'N')
		reverse_1 = reverse_1[::-1]
		#3' end juts out by 2 nucleotides past 5' end
		reverse_2 = reverse_1[1:]
		#3' end juts out by 3 nucleotides past 5' end
		reverse_3 = reverse_1[2:]
		#In the above illustrated example, these sequences would be
		#reverse_1: "CCAACAAGCGTCTCGCGTC..."
		#reverse_2: "CAACAAGCGTCTCGCGTC..."
		#reverse_3: "AACAAGCGTCTCGCGTC..."
		
		
		#For each configuration, check if each pair of corresponding nucleotides is pairing
		#Non-pairing pairs yield 0.0, and pairing pairs yield their weight if `weighted' is True,
		#or a flat 1.0 if False
		weights_1: list[float] = test_pairs(forward, reverse_1, weighted)
		weights_2: list[float] = test_pairs(forward, reverse_2, weighted)
		weights_3: list[float] = test_pairs(forward, reverse_3, weighted)
		#At this point, each `weights_*' contains the weight (or a flat 1.0/0.0) for each pair of
		#corresponding nucleotides from `forward' and the respective `reverse_*'
		#In the above illustrated example, the values in these lists would be
		#weights_1: 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.375, 0.0, 0.0, 0.375, 1.0,   0.0, 0.0, ...
		#weights_2: 0.0, 0.0, 0.0, 1.0, 0.5, 0.0, 0.0, 0.0,   1.0, 0.5, 1.0,   0.375, 0.0, 0.375...
		#weights_3: 0.0, 0.0, 1.0, 0.0, 0.5, 0.0, 1.0, 0.375, 0.0, 0.0, 0.375, 1.0,   1.0, 0.0, ...
		#if `weighted' is True, or
		#weights_1: 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0,   0.0, 0.0, 1.0,   1.0,   0.0, 0.0, ...
		#weights_2: 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0,   1.0, 1.0, 1.0,   1.0,   0.0, 1.0, ...
		#weights_3: 0.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 1.0,   0.0, 0.0, 1.0,   1.0,   1.0, 0.0, ...
		#if `weighted' is False
		
		
		#For each configuration, identify runs of pairing pairs, find the run whose weights sum up
		#to the highest value, and use that sum as the pairing score for that configuration
		self.traits["pair_score_1"] = max_sum_of_run(weights_1)
		self.traits["pair_score_2"] = max_sum_of_run(weights_2)
		self.traits["pair_score_3"] = max_sum_of_run(weights_3)
		
		
		
		#Check if pairs of nucleotides at specific positions of the intron are pairing
		#These positions are where the characteristic CAG & CTG are typically located in
		self.traits["pair_3_6"] = len(self) >= 6 and (self[3], self[-6]) in PAIR_WEIGHTS
		self.traits["pair_4_7"] = len(self) >= 7 and (self[4], self[-7]) in PAIR_WEIGHTS
		self.traits["pair_5_8"] = len(self) >= 8 and (self[5], self[-8]) in PAIR_WEIGHTS
	
	
	def get_ml_traits(self) -> array:
		"""
		Retrieve intron traits relevant to ML predictions. These are, in order:
		1. Does the preceding exon end with a pyrimidine?
		2. Does the intron start with a purine?
		3. Does the intron have the characteristic CAG & CTG at the right position?
		4. Does the intron end with a pyrimidine?
		5. Does the following exon start with a purine?
		6. Pairing score for configuration #2 (3' end juts out by 2 nucleotides)
		7. Pairing score for configuration #1 (3' end juts out by 1 nucleotide)
		8. Pairing score for configuration #3 (3' end juts out by 3 nucleotides)
		The first 5 traits are float-encoded booleans (1.0: True, 0.0: False),
		the latter 3 traits are ordinary floats.
		"""
		#Check that traits are actually present
		assert self.traits
		
		
		#Return dummy traits for very short introns
		if len(self) < 6:
			return array([0.0] * 8)
		
		
		#Extract the traits in the appropriate order, converting the bool traits to floats
		return array([
			self.traits["prev_exon_y"],
			self.traits["intron_r"],
			self.traits["intron_cagctg"],
			self.traits["intron_y"],
			self.traits["next_exon_r"],
			self.traits["pair_score_2"],
			self.traits["pair_score_1"],
			self.traits["pair_score_3"]
		])
	
	
	def get_stats(self, idx: int = 0) -> list[list[str|int|float|bool]]:
		"""
		Get statistics for an intron & all of its variants. The statistics are, in order:
		
		Per-intron statistics:
			the name of the gene of the intron
			the name of the transcript of the intron
			the index of the intron within the gene (taken from `idx')
			the number of variants the intron has
			the scaffold
			the strand
			the difference between the unified scores of the base intron and the best variant,
				as a measure of certainty of variant selection for this intron
		
		Per-variant statistics:
			the start coordinate
			the end coordinate
			the rank of the variant (the base intron is rank #0, variants are ranked by score)
			unified score
			conv score
			nonconv score
			the difference between the conv & nonconv score,
				as a measure of certainty of this variant's classification
			splice site
			whether the splice site is conventional
			last 5 nucleotides of the preceding exon
			first 10 nucleotides of the intron/variant
			last 10 nucleotides of the intron/variant
			first 5 nucleotides of the following exon
			whether the last nucleotide of the preceding exon is a pyrimidine
			whether the first nucleotide is a purine
			whether the last nucleotide is a pyrimidine
			whether the first nucleotide of the preceding exon is a purine
			whether the characteristic CAG & CTG are present at the approptiate positions
			whether the fourth and sixth-to-last nucleotide pair up
			whether the fifth and seventh-to-last nucleotide pair up
			whether the sixth and eighth-to-last nucleotide pair up
		"""
		#Sort variants by descending score
		self.variants.sort(key = lambda v: (v.unif_score, v.c_score, v.nc_score), reverse=True)
		
		
		#Build list containing the base intron and all variants for the purpose of other asserts
		assert ( allvars := [self] + self.variants )
		
		#Validate that the intron and all the variants all have the required information
		#Name of gene & transcript
		assert all( var.gene.name            for var in allvars )
		assert all( var.gene.transcript.name for var in allvars )
		#Scores
		assert all( var.scored()             for var in allvars )
		#Traits
		assert all( var.traits               for var in allvars )
		#Sequences
		assert all( var.has_seqs()           for var in allvars )
		
		
		#Measure of certainty of variant selection for this intron
		if self.variants:
			certainty: float = self.unif_score - self.variants[0].unif_score
		#Empty string if there are no variants
		else:
			certainty: str = ""
		
		
		#Will hold the returned stats
		stats: list[list[str|int|float|bool]] = []
		

		#Get statistics for the base intron
		stats.append([
			#Per-intron stats
			self.gene.name,							#Name of gene
			self.gene.transcript.name,				#Name of transcript
			idx,									#Index of intron within gene
			len(self.variants),						#Number of variants of intron
			self.scaffold,							#Scaffold
			self.strand,							#Strand
			certainty,								#Certainty of variant selection for intron
			#Per-variant stats
			self.start+1,							#Start position
			self.end,								#End position
			0,										#Variant's rank
			self.unif_score,						#Unified score
			self.c_score,							#Conventionality score
			self.nc_score,							#Nonconventionality score
			abs(self.c_score - self.nc_score),		#Certainty of variant classification
			self.splice_site or "",					#Splice site, if available
			self.traits["ss_is_conv"],				#Is the splice site conventional
			self.former_exon()[-5:],				#Last  5  nt of prev exon
			self[:10],								#First 10 nt of intron
			self[-10:],								#Last  10 nt of intron
			self.latter_exon()[:5],					#First 5  nt of next exon
			self.traits["prev_exon_y"],				#Is the prev exon's last  nt a pyrimidine
			self.traits["intron_r"],				#Is the intron's    first nt a purine
			self.traits["intron_y"],				#Is the intron's    last  nt a pyrimidine
			self.traits["next_exon_r"],				#Is the next exon's first nt A purine
			self.traits["intron_cagctg"],			#Are the characteristic CAG & CTG present
			self.traits["pair_score_1"],			#Pairing score
			self.traits["pair_score_2"],			#Pairing score
			self.traits["pair_score_3"],			#Pairing score
			self.traits["pair_3_6"],				#Do the nts at positions 3 and -6 pair up
			self.traits["pair_4_7"],				#Do the nts at positions 4 and -7 pair up
			self.traits["pair_5_8"]					#Do the nts at positions 5 and -8 pair up
		])
		
		
		#Get statistics for each variant
		for n, var in enumerate(self.variants, 1):
			stats.append([
				#Per-intron stats are only reported for optimal variants
				"",									#Name of gene
				"",									#Name of transcript
				"",									#Index of intron within gene
				"",									#Number of variants of intron
				"",									#Scaffold
				"",									#Strand
				"",									#Certainty of variant selection for intron
				#Per-variant stats
				var.start+1,						#Start position
				var.end,							#End position
				n,									#Variant's rank
				var.unif_score,						#Unified score
				var.c_score,						#Conventionality score
				var.nc_score,						#Nonconventionality score
				abs(var.c_score - var.nc_score),	#Certainty of variant classification
				var.splice_site or "",				#Splice site, if available
				var.traits["ss_is_conv"],			#Is the splice site conventional
				var.former_exon()[-5:],				#Last  5  nt of prev exon
				var[:10],							#First 10 nt of intron
				var[-10:],							#Last  10 nt of intron
				var.latter_exon()[:5],				#First 5  nt of next exon
				var.traits["prev_exon_y"],			#Is the prev exon's last  nt a pyrimidine
				var.traits["intron_r"],				#Is the intron's    first nt a purine
				var.traits["intron_y"],				#Is the intron's    last  nt a pyrimidine
				var.traits["next_exon_r"],			#Is the next exon's first nt A purine
				var.traits["intron_cagctg"],		#Are the characteristic CAG & CTG present
				var.traits["pair_score_1"],			#Pairing score
				var.traits["pair_score_2"],			#Pairing score
				var.traits["pair_score_3"],			#Pairing score
				var.traits["pair_3_6"],				#Do the nts at positions 3 and -6 pair up
				var.traits["pair_4_7"],				#Do the nts at positions 4 and -7 pair up
				var.traits["pair_5_8"]				#Do the nts at positions 5 and -8 pair up
			])
		
		
		return stats
	
	
	#List of descriptors for every stat returned by Intron.get_stats()
	STATS: list[str] = [
		#Per-intron stats
		"gene",
		"transcript",
		"intron_idx",
		"variant_cnt",
		"scaffold",
		"strand",
		"score_outpace",
		#Per-variant stats
		"start",
		"end",
		"variant_rank",
		"unif_score",
		"c_score",
		"nc_score",
		"score_delta",
		"splice_site",
		"splice_site_is_conv",
		"e-5",
		"i10",
		"i-10",
		"e5",
		"prev_exon_y",
		"start_r",
		"end_y",
		"next_exon_r",
		"cagctg",
		"pairing_score_1",
		"pairing_score_2",
		"pairing_score_3",
		"pair_3_-6",
		"pair_4_-7",
		"pair_5_-8"
	]


def concat_genseq(*seqs: GenomicSequence, sep: str = '') -> str:
	"""
	Given GenomicSequence objects (subclasses also ok), extract the sequence of each one and
	concatenate them, with an optional separator between each sequence
	"""
	return sep.join( map( attrgetter("seq"), seqs ) )


def hash_genseq(seq: GenomicSequence) -> int:
	"""
	Hash function for GenomicSequence (& subclass) objects
	The intention of this hashing function is that if a pair of objects would evaluate as equal
	using the '==' operator, they have identical hashes
	This function can be hooked up to the GenomicSequence class to make it hashable:
	"GenomicSequence.__hash__ = hash_genseq"
	"""
	return hash(f"{seq.scaffold} {seq.start} {seq.end} {seq.strand}")


def hash_genseq_unstranded(seq: GenomicSequence) -> int:
	"""
	Hash function for GenomicSequence (& subclass) objects
	Unlike hash_genseq(), two objects will have identical hashes if they have the same scaffold
	and coordinates, but not necessarily the same strand
	This function can be hooked up to the GenomicSequence class to make it hashable:
	"GenomicSequence.__hash__ = hash_genseq"
	"""
	return hash(f"{seq.scaffold} {seq.start} {seq.end}")


###################################################################################################
#   DESERIALIZATION
###################################################################################################

def deserialize_gff(path: str) -> dict[str,Gene]:
	"""
	Given a path to a GFF file, deserialize it to a dictionary containing all "gene" features
	listed therein, with each gene's "mRNA" & "exon" feature(s) linked to it
	"""	
	#Dict of created genes
	genes: dict[str,Gene] = {}
	#List containing info about exons to be created after parsing the whole file
	deferred_exons: list[tuple] = []
	#Table for converting from transcript IDs to gene IDs
	mrna_to_gene: dict[str,str] = {}
	
	
	#for scaffold, source, type_, start, end, strand, 
	for entry in parse_gff(path):
		#Ensure each entry has an ID attribute
		assert "ID" in entry.attrs
		
		#Reject features with a '.' strand
		if entry.strand not in '-+':
			continue
		
		
		scaffold: str           = entry.seqid
		type_:    str           = entry.type_
		start:    int           = entry.start-1
		end:      int           = entry.end
		strand:   Strand        = entry.strand
		attrs:    dict[str,str] = entry.attrs
		id_:      str           = attrs["ID"]
		
		
		#Process "gene", "mRNA" & "exon" features appropriately
		if type_ == "gene":
			#Instantiate new gene
			genes[id_] = Gene(scaffold, start, end, name=id_, strand=strand, exons=[])
		
		
		elif type_ == "mRNA":
			#Ensure that mRNA feature has a Parent
			assert "Parent" in attrs
			
			#Submit to conversion table
			mrna_to_gene[id_] = attrs["Parent"]
		

		elif type_ == "exon":
			#Ensure that exon feature has a Parent
			assert "Parent" in attrs
			
			#Keep track of this exon to add to a gene later
			deferred_exons.append((attrs["Parent"], scaffold, start, end, strand))
	
	
	#Instantiate all exons and link them to their parent gene
	for mrna_id, scaffold, start, end, strand in deferred_exons:
		#Get the name of the exon's grandparent gene
		assert mrna_id in mrna_to_gene
		gene_id: str = mrna_to_gene[mrna_id]
		assert gene_id in genes
		
		
		#Get the grandparent gene
		gene: Gene = genes[gene_id]
		
		#Instantiate the exon and add it to the gene's exon list
		gene.exons.append(
			Exon(
				scaffold, start, end, strand=strand, gene=gene
			)
		)
	
	
	#Submit the name of each transcript to its gene
	for mrna_id, gene_id in mrna_to_gene.items():
		assert gene_id in genes
		
		genes[gene_id].transcript.name = mrna_id
	
	
	return genes



###################################################################################################
#   INTRON SCORING
###################################################################################################

#Nucleotide pair weights.
#Each pair of nucleotides which actually pair up in pre-mRNA have the following weights:
#A-T: 0.5, C-G: 1.0, G-T: 0.375.
#Non-pairing nucleotide pairs have an implicit weight of 0.0.
#For commutativity, each pairing pair is present in this dictionary in both orders.
PAIR_WEIGHTS: dict[tuple[str,str],float] = {
	('A','T'): 0.5,   ('T','A'): 0.5,
	('G','C'): 1.0,   ('C','G'): 1.0,
	('G','T'): 0.375, ('T','G'): 0.375
}


#Splice sites considered conventional
CONV_SS: set[str] = { "GTAG", "GCAG", "CTAC", "CTGC" }


#Purines & pyrimidines
PURINES:     str = "AGR"
PYRIMIDINES: str = "CTY"


def test_pairs(seq1: str, seq2: str, weighted: bool = True) -> list[float]:
	"""
	Given two pairs of sequences, check if each pair of corresponding nucleotides actually pair up
	and return a numeric value for that pair.
	If `weighted' is True, the numeric value for each pair is that pair's weight.
	If `weighted' is False, the value is 1.0 for pairing pairs, and 0.0 for non-pairing pairs.
	"""
	#Will hold the result for each pair
	out: list[float] = []
	

	#If the results should be weighted
	if weighted:
		#Iterate over pairs of corresponding nucleotides
		for pair in zip(seq1,seq2):
			#Get weight for this pair
			out.append( PAIR_WEIGHTS.get( pair, 0.0 ) )
	
	#If the results should not be weighted
	else:
		#Iterate over pairs of corresponding nucleotides
		for pair in zip(seq1,seq2):
			#1.0 if the this is a pairing pair, 0.0 otherwise
			out.append( 1.0 if pair in PAIR_WEIGHTS else 0.0 )
	

	return out


def max_sum_of_run(iterable: Iterable[float|int]) -> float:
	"""
	Given an iterable of floats or ints, identify runs of positive values, find the run which sums
	up to the largest value, and return that sum.
	"""
	#Sum of the highest-sum run found so far
	max_sum = 0.0
	#Sum of the current run
	cur_sum = 0.0
	
	for i in iterable:
		#If this is the start or continuation of a run
		if i > 0.0:
			cur_sum += i
			max_sum = cur_sum if cur_sum >= max_sum else max_sum
		
		#If this is the end of a run, or not part of a run
		else:
			cur_sum = 0.0
	
	return max_sum


def score_introns(
	genes:              Iterable[Gene],
	nc_model,
	c_model,
	*,
	unif_score_from_ss: bool           = False,
	batch_size:         int            = 0,
	weighted:           bool           = True
):
	"""
	Calculate conventional & nonconventional structure compatibility scores for all introns
	(& variants) of the passed genes, using the supplied models.

	Each scored intron is given a conventionality score and a nonconventionality score, plus a
	unified score, which is the higher score of the two.

	Alternatively, if `unif_score_from_ss' is True, introns with conventional splice sites are
	instead	forced to have the maximum possible unified score, and introns without take the
	nonconventional score as the unified score

	By default, all introns are scored at once, but if number of introns is very large, an OOM
	error might occur. To remedy this, the introns can be processed in batches by passing a
	positive integer to `batch_size'. This value will be the size of a single batch of introns.
	Smaller values reduce the risk of an OOM error, but slow down scoring.
	
	`weighted' controls whether the introns' pairing scores are calculated as weighted/unweighted.
	"""
	#Build list of introns to assess
	phase("Score introns: gather")
	introns: list[Introns] = []
	for gene in genes:
		for intron in gene.introns:
			#Include the intron and all its variants in assessment
			introns.append(intron)
			introns += intron.variants
	
	assert len(introns) > 0
	
	
	#Compute traits for each intron
	phase("Score introns: compute traits")
	for intron in introns:
		intron.add_traits(weighted)
	
	
	#Retrieve ML traits and splice site
	phase("Score introns: retrieve traits")
	ml_traits: list[array] =  [ intron.get_ml_traits() for intron in introns ]
	
	#Check number of trait sets
	assert len(introns) == len(ml_traits)
	
	
	#Nonconv scores
	phase("Score introns: nonconv scores")
	#Instantiate nonconv scores as a vector of floats
	nc_scores: array = empty([0], dtype="float64")
	
	#Split introns' traits into batches of size `batch'
	#If `batch' is <= 0, the traits are all put into one big batch and effectively unbatched
	for batch in batched(ml_traits, batch_size if batch_size > 0 else len(ml_traits)):
		#nonconv_model.predict_proba(batch) returns a matrix with 2 columns, and `batch_size' rows
		#Each row contains the conv & nonconv scores (first & second item, respectively), according
		#to the nonconv model
		#The nonconv scores are extracted into a vector (using [:,1]), and that vector is appended
		#to the overall nonconv score vector (with append())
		nc_scores = append(nc_scores, nc_model.predict_proba(batch)[:,1])
	
	#Check number of nonconv scores
	assert len(introns) == len(nc_scores)
	
	
	#Conv scores
	phase("Score introns: conv scores")
	#Instantiate conv scores as a vector of floats
	c_scores: array = empty([0], dtype="float64")
	
	#Split traits into batches
	for batch in batched(ml_traits, batch_size if batch_size > 0 else len(ml_traits)):
		#Just like in the previous for loop, except the overall conv score vector is appended to,
		#the conv model is used for prediction, and conv scores are extracted (with [:,0])
		c_scores = append(c_scores, c_model.predict_proba(batch)[:,0])
	
	#Check number of nonconv scores
	assert len(introns) == len(c_scores)
	
	
	#Convert NumPy arrays to Python list, which also converts NumPy floats to regular Python floats
	c_scores = c_scores.tolist()
	nc_scores = nc_scores.tolist()

	
	#Submit each intron's scores to the actual object
	phase("Score introns: save scores")
	for intron, c_score, nc_score in zip(introns, c_scores, nc_scores):
		#Per-class score
		intron.c_score  = c_score
		intron.nc_score = nc_score
		#Unfied score
		if unif_score_from_ss:
			intron.unif_score = 1.0 if intron.traits["ss_is_conv"] else nc_score
		else:
			intron.unif_score = max( c_score, nc_score )




