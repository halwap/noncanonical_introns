###################################################################################################
#   IMPORTS & SETUP
###################################################################################################

#Standard library 
#Efficient lambdas
from operator import attrgetter
#Convenient iteration
from itertools import pairwise, batched
#Type hints
from typing import *
#Model deserialization
from pickle import load
#Diagnostics
from warnings import filterwarnings
#Phase duration timing
from time import time
#Stderr printing
from sys import stderr
#Statistics
from statistics import mean
#Enum support
from enum import Enum


#Third-party imports
#FASTA deserialization
from Bio.SeqIO.FastaIO import SimpleFastaParser
#Retrieving seqs from the '-' strand
from Bio.Seq import reverse_complement
#Convenient iteration
from more_itertools import roundrobin
#Array processing
from numpy import array, append, empty, argmax


#Hide warnings
filterwarnings("ignore", category=UserWarning)




###################################################################################################
#   CLASSES
###################################################################################################

class GenomicSequence:
    scaffold: str
    start: int
    end: int
    seq: str|None
    strand: Literal['+', '-', '.']|None

    def __init__(self, scaffold, start, end, seq=None, strand=None):
        self.scaffold = scaffold
        

        #Validate coordinates
        assert start < end, \
            f"Invalid coords: start: {start}, end: {end}"

        self.start = start
        self.end   = end
        
        
        #Validate strand
        assert strand, \
            f"Bad strand: no value passed"
        assert len(strand) == 1 and strand in "+-.", \
            f"Bad strand: {strand} is not one of '+', '-', or '.'"
        
        self.strand = strand
        
        
        if seq:
            #Validate sequence length vis-a-vis passed coords
            assert len(seq) == end - start, \
                f"Length mismatch: {len(seq)} by seq, {end - start} by coords"
        self.seq = seq
    
    
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
        """
        return len(self.seq)
    
    def __getitem__(self, index: int|slice):
        """
        Index and slice nucleotides of sequence
        """
        return self.seq[index]
    

    def __add__(self, other: Type['GenomicSequence']|str) -> str:
        """
        Concatenate a GenomicSequence (or subclass)'s sequence with another one's, or with a str
        """
        if not ( issubclass( type(other), GenomicSequence ) or isinstance(other, str) ):
            raise TypeError(f"Can't concatenate {type(self).__name__} & {type(other).__name__}")
        
        if isinstance(other, str):
            return self.seq + other
        else:
            return self.seq + other.seq
    

    def __radd__(self, other: Type['GenomicSequence']|str) -> str:
        """
        Concatenate a GenomicSequence (or subclass)'s sequence with another one's, or with a str
        """
        if not ( issubclass( type(other), GenomicSequence ) or isinstance(other, str) ):
            raise TypeError(f"Can't concatenate {type(self).__name__} & {type(other).__name__}")
        
        if isinstance(other, str):
            return other + self.seq
        else:
            return other.seq + self.seq 
    
    
    #Comparisons for checking strand, scaffold and relative position of two GenomicSequence objects
    def same_scaff_and_strand(self, other: Type['GenomicSequence']) -> bool:
        """
        Check whether two GenomicSequence objects are from the same strand & scaffold
        """
        return self.strand == other.strand and self.scaffold == other.scaffold
    
    def __lt__(self, other: Type['GenomicSequence']) -> bool:
        """
        Check whether a GenomicSequence is entirely to the right of `self', without flushed borders
        Comparison is only relevant if both objects are from the same scaffold & strand
        """
        return self.same_scaff_and_strand(other) \
           and self.end < other.start
    
    def __le__(self, other: Type['GenomicSequence']) -> bool:
        """
        Check whether a GenomicSequence is entirely to the right of `self', with flushed borders
        Comparison is only relevant if both objects are from the same scaffold & strand
        """
        return self.same_scaff_and_strand(other) \
           and self.end == other.start
    
    def __contains__(self, other: Type['GenomicSequence']) -> bool:
        """
        Check whether `self' completely overlaps another GenomicSequence
        Comparison is only relevant if both objects are from the same scaffold & strand
        """
        return self.same_scaff_and_strand(other) \
           and self.start <= other.start \
           and other.end <= self.end
    
    def __eq__(self, other: Type['GenomicSequence']) -> bool:
        """
        Check whether two GenomicSequence objects have the same coordinates, scaffold & strand
        """
        return self.same_scaff_and_strand(other) \
           and self.start == other.start \
           and other.end == self.end
    
    
    #Serialization
    def emit_gff(self, fd: TextIO, ft: str, attr: dict[str:str]|None = None):
        """Serialize a GenomicSequence object to file `fd'
        
        Params:
            fd: TextIO
                File descriptor to write to
            ft: str
                Feature type to save sequence as, e.g. "gene", "mRNA" etc.
            attr: dict[str:str] | None
                A dictionary containing the attributes of the sequence, or None
                if it should not have any attributes
        """
        assert fd and fd.mode == 'w', \
            f"Cannot serialize {self} to given file"
        
        assert self.strand in "+-.", \
            f"{self} has invalid strand"
        
        
        #Convert attributes from a dict to a str
        
        
        fd.write(
            '\t'.join((
                self.scaffold,      #seqid
                ".",                #source
                ft,                 #type
                #+1 because GenomicSequence objects store their coords differently than in GFFs
                str(self.start+1),  #start
                str(self.end),      #end
                '.',                #score
                self.strand,        #strand
                '.',                #phase
                #Convert attributes dictionary to a GFF attribute strin
                ';'.join( map ( '='.join, attr.items() ) ) if attr else '.'
            )) + '\n'
        )


class Gene(GenomicSequence):
    scaffold: str
    start: int
    end: int
    seq: str|None
    strand: Literal['+', '-', '.']|None

    transcript: 'Transcript'
    exons: list['Exon']
    introns: list['Intron']
    name: str

    def __init__(self, scaffold, start, end, seq='', strand='',
                 transcript=None, exons=None, introns=None, name=''):
        GenomicSequence.__init__(self, scaffold, start, end, seq=seq, strand=strand)
        self.transcript      = transcript or Transcript(scaffold, start, end, strand=strand)
        self.exons           = exons or []
        self.introns         = introns or []
        self.name            = name
    
    
    def fixup_exons(self):
        '''
        When exons are pushed as-is to Gene.exons, they might not necessarily be sorted
        ordered and linked to one another -  this methods applies the necessary fixes
        '''
        #If there's no exons, or only 1 exon, there's nothing to do
        if len(self.exons) < 2:
            return
        
        #Sort exons by position
        self.exons.sort(key=attrgetter("start"))
        
        #Validate coords of exons
        assert self.valid_coords(), \
            f"Exons of {self} have bad coords"
        
        
        #Link exons
        for left_exon, right_exon in pairwise(self.exons):
            left_exon.next_exon = right_exon
            right_exon.prev_exon = left_exon
        

        #Validate linkage of exons
        assert self.valid_links(), \
            f"Exons of {self} unlinked"
    
    
    def add_seqs(self, genome: dict[str:str], which: str|set[str] = "gtei"):
        """
        Given a deserialized fasta file, build the sequences of the gene/exons/introns
        `which' determines which specific features should obtain new sequences
        By default, all 3 features (gene, exons & introns) will get new seqs
        """
        def maybe_rc(seq: str) -> str:
            """
            Get the reverse complement of a seq, if the strand of the gene necessitates it
            """
            return reverse_complement(seq) if self.strand == '-' else seq
        
        
        #Validate coordinates of exons, and introns, if requested
        assert self.valid_coords('i' in which), \
            f"Exons/introns of {self} have bad coords"
        
        
        #Get sequence of appropriate scaffold
        #This requires the gene's scaffold to be present in the genome
        assert self.scaffold in genome, \
            f"Scaffold of {self} not in genome"
        scaffold_seq = genome[self.scaffold]
        
        
        #Push sequence to gene
        if 'g' in which:
            self.seq = maybe_rc(scaffold_seq[self.start:self.end])
            assert len(self) == self.end - self.start, \
                f"{self}: bad length, {len(self)} vs {self.end - self.start}"
        
        
        #Push sequence to transcript
        if 't' in which:
            self.transcript.seq = ""
            for exon in self.ordered_exons():
                self.transcript.seq += maybe_rc(scaffold_seq[exon.start:exon.end])
        
        
        #Push sequences to exons
        if 'e' in which:
            for exon in self.exons:
                exon.seq = maybe_rc(scaffold_seq[exon.start:exon.end])
                assert len(exon) == exon.end - exon.start, \
                    f"{exon} in {self}: bad length, {len(exon)} vs {exon.end - exon.start}"
        
        
        #Push sequences to introns
        if 'i' in which:
            for intron in self.introns:
                intron.seq = maybe_rc(scaffold_seq[intron.start:intron.end])
                assert len(intron) == intron.end - intron.start, \
                    f"{intron} in {self}: bad length, {len(intron)} vs {intron.end - intron.start}"
        
        
        #Check that all created sequences match
        assert self.valid_seqs('i' in which), \
            f"Seq of {self} differs from joint seq of constituent features"

    
    def add_introns(self):
        """
        Populate self.introns by inserting an intron between every pair of exons
        """
        #Check that there are exons at all
        assert len(self.exons) > 0, \
            f"{self} has no exons"
        
        #Validate coordinates & linkage of exons
        assert self.valid_coords(), \
            f"Exons of {self} have bad coords"
        assert self.valid_links(), \
            f"Exons of {self} unlinked"
        

        #If there's only one exon, there are no introns to add; exit early
        if len(self.exons) == 1:
            return
        
        
        self.introns = []
        
        
        for left_exon, right_exon in pairwise(self.exons):

            #Get coordinates of the intron
            start = left_exon.end
            end = right_exon.start
            
            #Instantiate intron & add it to the list
            self.introns.append(
                Intron(
                    self.scaffold, start, end, strand = self.strand,
                    gene = self, prev_exon = left_exon, next_exon = right_exon
                )
            )
        
        
        #Check exon & intron count
        assert len(self.exons) == len(self.introns) + 1, \
            f"{self} has {len(self.exons)} exons, but {len(self.introns)} introns"
        
        #Check coords & linkage, this time with introns
        assert self.valid_coords(True), \
            f"Exons/introns of {self} have bad coords"
        assert self.valid_links(True), \
            f"Exons/introns of {self} unlinked"
    

    def add_intron_variants(self, min_exon_len: int):
        """
        Build a list of variants for each intron in the gene
        """
        for intron in self.introns:
            intron.get_variants(min_exon_len)
    
    
    def rectify_introns(self):
        """
        Rectify all introns in the gene by selecting the highest-scored variant for each one.
        Assumes introns have been scored already.
        """
        
        #This method should do nothing if the gene has no introns
        if len(self.introns) == 0:
            return
        
        #Check that all introns have been scored
        assert all( intron.unif_score != None for intron in self.introns ), \
            f"{self} has unscored intron(s)"
        

        for intron_idx, intron in enumerate(self.introns):
            #Skip introns without variants
            if intron.variants:
                
                #Check that all variants of this intron have been scored
                assert all( variant.unif_score != None for variant in intron.variants ), \
                    f"{intron} of {self} has unscored variant(s)"
                
                
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
        for exon_idx, (left_intron, right_intron) in enumerate(pairwise(self.introns), 1):
            if not left_intron < right_intron:
                raise ValueError("Vanishing exon found, try increasing minimum exon length")
            
            
            #If the introns disagree on the exon between them, reconciliation is necessary
            if left_intron.next_exon != right_intron.prev_exon:
                
                #Reconciliation is done by having the left & right introns dictate the start & end
                #of the exon, respectively
                left_intron.next_exon.end = right_intron.prev_exon.end
                
                
                #Link new exon & right intron (left exon & intron are already linked)
                right_intron.prev_exon = left_intron.next_exon
                left_intron.next_exon.next_intron = right_intron

            
            #At this point, the introns should agree on the exon between them, either naturally, or
            #thanks to reconciliation
            assert left_intron.next_exon == right_intron.prev_exon, \
                f"{left_intron} & {right_intron} disagree on their in-between exon, respectively"
            
            #Replace exon in list
            self.exons[exon_idx] = left_intron.next_exon
            

            #Check linkage
            assert left_intron @ self.exons[exon_idx], \
                f"{left_intron} & {self.exons[exon_idx]} of {self} unlinked"
            assert self.exons[exon_idx] @ right_intron, \
                f"{self.exons[exon_idx]} & {right_intron} of {self} unlinked"
            
            #Check coordinates (intron-exon only, since introns were checked already)
            assert left_intron <= self.exons[exon_idx], \
                f"{left_intron} & {self.exons[exon_idx]} of {self} not flush"
            assert self.exons[exon_idx] <= right_intron , \
                f"{self.exons[exon_idx]} & {right_intron} of {self} not flush"
    
    
    def swap_intron_and_variant(self, intron_idx: int, variant_idx: int):
        """
        Take self.introns[intron_idx] and self.introns[intron_idx].variants[variant_idx], and swap
        them around.
        """
        #Validate input
        assert intron_idx < len(self.introns), \
            f"Attempting to swap nonexistent intron of {self}"
                
        #Get the base intron
        intron: Intron = self.introns[intron_idx]
        
        
        #Validate input
        assert variant_idx < len(intron.variants), \
            f"Attempting to swap nonexistent variant of {intron} of {self}"
        
        #Get the variant
        variant: Intron = intron.variants[variant_idx]
        
        
        #Make `intron' a variant
        intron.variants[variant_idx] = intron
        
        #Repossess the list of variants
        variant.variants = intron.variants
        del intron.variants
        
        #Make `variant' a base intron
        self.introns[intron_idx] = variant
    
    
    def serialize(self, fd: TextIO):
        """
        Serialize a gene and all its children features (transcript, exons, introns) to a GFF file
        """
        #Validate coordinates
        assert self.valid_coords(True), \
            f"Exons/introns of {self} have bad coords"
        
        
        #Check that both the gene and the transcript have names to use as their IDs in the GFF
        assert self.name, \
            f"{self} lacks a name"
        assert self.transcript.name, \
            f"{self.transcript} lacks a name"
        
        
        #GFF attributes dictionary
        attr: dict[str:str] = {}
        
        
        #Use gene name as ID
        attr["ID"] = self.name 
        
        #Get mean intron score, if there are scored introns
        if self.introns and all( intron.unif_score != None for intron in self.introns ):
            attr["avg_intron_score"] = str( mean( map( attrgetter("unif_score" ), self.introns ) ))
        
        #Serialize gene
        self.emit_gff(fd, "gene", attr)
        
        
        #Modify attributes for transcript - make the gene its parent, and get its ID
        attr["Parent"] = attr["ID"]
        attr["ID"] = self.transcript.name
        
        #Remove mean intron score for sub-gene features
        if "avg_intron_score" in attr:
            del attr["avg_intron_score"]
        
        #Serialize transcript
        self.transcript.emit_gff(fd, "mRNA", attr)
        
        
        #Modify attributes for exons & introns - make the transcript their parent
        attr["Parent"] = attr["ID"]
        
        #Serialize exons
        for n, exon in enumerate(self.exons, 1):
            #Generate ID attribute
            attr["ID"] = f"{attr["Parent"]}.exon{n}"

            exon.emit_gff(fd, "exon", attr)
        
        #Serialize introns
        for n, intron in enumerate(self.introns, 1):
            #Generate ID attribute
            attr["ID"] = f"{attr["Parent"]}.intron{n}"
            
            #Get splice site attribute
            if intron.splice_site:
                attr["splice_site"] = intron.splice_site + '/' + \
                                      reverse_complement(intron.splice_site)
            
            #Get conventional score attribute
            if intron.c_score != None:
                attr["conv_score"] = str(intron.c_score)
            
            #Get nonconventional score attribute
            if intron.nc_score != None:
                attr["nonconv_score"] = str(intron.nc_score)
            
            
            intron.emit_gff(fd, "intron", attr)
    
    
    def ordered_exons(self) -> Iterable['Exon']:
        """
        Returns an iterable of exons, sorted by order within gene, not position
        """
        #Check all exons are ordered
        assert self.valid_coords(), \
            f"Exons of {self} have bad coords"
        
        #Return the exons, ordered as they are within the gene
        return reversed( self.exons ) if self.strand == '-' else self.exons
    
    
    def ordered_exons_and_introns(self) -> Iterable[Union['Exon','Intron']]:
        """
        Returns an iterable of interspersed exons & introns, sorted by order within gene, not
        position
        """
        assert self.valid_coords(True), \
            f"Exons/introns of {self} have bad coords"
        
        #Intersperse exons & introns, reversing their order if needed
        if self.strand == '-':
            return roundrobin( reversed(self.exons), reversed(self.introns) )
        else:
            return roundrobin( self.exons, self.introns)
    
    
    def valid_coords(self, introns_too: bool = False) -> bool:
        """
        Check if the coordinates of all the exons and introns of the gene are valid: this includes
        that exons/introns are within the bounds of the gene, are properly ordered, and are
        appropriately flush/unflush.

        By default only exons are checked; `introns_too' should be set to True to also check
        introns. Additionally, exons & introns are checked iff the gene has any.
        """
        #Check exons are within bound of gene
        if any( exon not in self for exon in self.exons ):
            return False
        #Check order and non-flushness of exons
        if any( not l_exon < r_exon for l_exon, r_exon in pairwise(self.exons) ):
            return False
        
        #If intron should be checked
        if introns_too:
            #Check introns are within bound of gene
            if any( intron not in self for intron in self.introns ):
                return False
            #Check order and non-flushness of introns
            if any( not l_intron < r_intron for l_intron, r_intron in pairwise(self.introns) ):
                return False
            #Check each exon is before and flush with its following intron
            if any( not exon <= intron for exon, intron in zip(self.exons, self.introns) ):
                return False
            #Check each intron is before and flush with its following exon
            if any( not intron <= exon for intron, exon in zip(self.introns, self.exons[1:]) ):
                return False
        
        #If all relevant checks passed, the coordinates are valid
        return True
    
    
    def valid_links(self, introns_too: bool = False) -> bool:
        """
        Check if the linkage of all the exons (and optionally, introns) of the gene are valid
        """
        #Check that exons are mutually linked
        if any( not l_exon @ r_exon for l_exon, r_exon in pairwise(self.exons) ):
            return False
        
        #If intron should be checked
        if introns_too:
            #Check linkage between each exon and its following intron
            if any( not exon @ intron for exon, intron in zip(self.exons, self.introns) ):
                return False
            #Check linkage between each intron and its following exon
            if any( not intron @ exon for intron, exon in zip(self.introns, self.exons[1:]) ):
                return False
        
        #If all relevant checks passed, the links are valid
        return True
    
    
    def valid_seqs(self, introns_too: bool = False) -> bool:
        """
        Check if the exons' combined sequence matches the transcript's, and optionally, the exons
        and introns' combined sequence matches the gene's
        """
        #Check that the transcript & exons all have sequences
        if not self.transcript.seq:
            return False
        if any( not exon.seq for exon in self.exons ):
            return False
        
        #Check if the exons' combined seq matches the transcript's seq
        if concat_genseq(self.ordered_exons()) not in self.transcript.seq:
            return False
        
        
        #If intron should be checked
        if introns_too:
            #Check that the gene & introns all have sequences
            if not self.seq:
                return False
            if any( not intron.seq for intron in self.introns ):
                return False
            

            #Check if the exons' & introns' combined seq matches the gene's seq
            if concat_genseq(self.ordered_exons_and_introns()) not in self.seq:
                return False
        

        #If all relevant checks passed, the sequences are valid
        return True



class Transcript(GenomicSequence):
    scaffold: str
    start: int
    end: int
    seq: str
    strand: Literal['+', '-', '.']|None
    
    name: str

    def __init__(self, scaffold, start, end, seq='', strand='', name=''):
        GenomicSequence.__init__(self, scaffold, start, end, seq=seq, strand=strand)
        self.name = name


class Exon(GenomicSequence):
    scaffold: str
    start: int
    end: int
    seq: str|None
    strand: Literal['+', '-', '.']|None

    gene: Gene|None
    prev_exon: 'Exon'|None
    next_exon: 'Exon'|None
    prev_intron: 'Intron'|None
    next_intron: 'Intron'|None


    def __init__(self, scaffold, start, end, seq='', strand='',
                 gene=None, prev_exon=None, next_exon=None, prev_intron=None, next_intron=None):
        GenomicSequence.__init__(self, scaffold, start, end, seq=seq, strand=strand)
        
        #If preceeding intron was passed
        if prev_intron:
            #Validate input
            assert prev_intron <= self, \
                f"Intron-exon misalignment: {prev_intron} & {self}"
            
            #Link exon to intron
            prev_intron.next_exon = self
        
        #If following intron was passed
        if next_intron:
            #Validate input
            assert self <= next_intron, \
                f"Exon-intron misalignment: {self} & {next_intron}"
        
            #Link exon to intron
            next_intron.prev_exon = self
        
        #Prospectively link introns to exon
        self.prev_intron = prev_intron
        self.next_intron = next_intron
        
        
        #If gene was passed
        if gene:
            #Validate input
            assert self in gene, \
                f"{self} outside of {gene}"
        
        #Prospectiely link gene to exon
        self.gene = gene
        
        
        #If preceeding exon was passed
        if prev_exon:
            #Validate input
            assert prev_exon < self,  \
                f"Exon-exon misalignment: {prev_exon} & {self}"
            
            #Link this exon to that one
            prev_exon.next_exon = self
        
        #If following exon was passed
        if next_exon:
            #Validate input
            assert self < next_exon, \
                f"Exon-exon misalignment: {self} & {next_exon}"
            #Link this exon to that one
            next_exon.prev_exon = self
        
        #Prospectively links exons to this one
        self.prev_exon = prev_exon
        self.next_exon = next_exon
    
    
    def __matmul__(self, other: Type[GenomicSequence]) -> bool:
        """
        Check if this exon and a following exon or intron are linked together
        Call with "self @ other_exon" etc.
        This test is non-commutative: the left operand is expected to be "before" the right operand
        """
        if not isinstance(other, (Exon, Intron)):
            raise TypeError(f"Can't test links for {type(self).__name__} & {type(other).__name__}")
                
        #Get the linked-to exons/intron
        left = other.prev_exon
        right = self.next_exon if isinstance(other, Exon) else self.next_intron
        
        return left and right and (left == self) and (right == other)
    
    
    def __rmatmul__(self, other: Type[GenomicSequence]) -> bool:
        """
        Check if this exon and a preceedgin exon or intron are linked together
        Call with "intron @ self" etc.
        This test is non-commutative: the left operand is expected to be "before" the right operand
        """
        if not isinstance(other, (Exon, Intron)):
            raise TypeError(f"Can't test links for {type(other).__name__} & {type(self).__name__}")
                
        left = self.prev_exon if isinstance(other, Exon) else self.prev_intron
        right = other.next_exon
                
        return left and right and (left == other) and (right == self)


class Intron(GenomicSequence):
    """
    This is a class for representing introns.
    
    Attributes:
        scaffold: str
            Name of sequence on which the intron is located
        start: int
            Coordinate of first nucleotide of intron within scaffold
        end: int
            Coordinate of last nucleotide of intron within scaffold
        gene: Gene|None
            The gene this intron is a part of, if one was provided
        prev_exon: Exon|None
            The exon preceding the intron (by position in scaffold, not in-gene order), if provided
        next_exon: Exon|None
            The exon following the intron (by position in scaffold, not in-gene order), if provided
        strand: str|None
            Which strand (+ or -) the intron is located on, if it was provided
        seq: str|None
            Nucleotide sequence of the intron, if one was provided
        variants: list[Intron]
            Variants (alt positions) of the intron; empty if it has none, or if not yet computed
        c_score: float|None
            Conventionality score of the intron, or None if not yet computed
        nc_score: float|None
            Nonconventionality score of the intron, or None if not yet computed
        unif_score: float|None
            Unified score of the intron, or None if not yet computed
        splice_site: str|None
            The splice site of the intron, or None if not yet computed
        traits: dict[str:bool|float]
            Structural traits of the intron; empty if not yet computed
    """
    scaffold: str
    start: int
    end: int
    seq: str|None
    strand: Literal['+', '-', '.']|None

    gene: Gene|None
    prev_exon: 'Exon'|None
    next_exon: 'Exon'|None
    variants: list['Intron']
    c_score: float|None
    nc_score: float|None
    unif_score: float|None
    splice_site: str|None
    trais: dict[str:bool|float]


    def __init__(
                self,
                scaffold:  str,
                start:     int,
                end:       int,
                seq:       str|None  = None,
                strand:    str|None  = None,
                gene:      Gene|None = None,
                prev_exon: Exon|None = None,
                next_exon: Exon|None = None,
            ):
        GenomicSequence.__init__(self, scaffold, start, end, seq=seq, strand=strand)
        

        #Validate `gene'
        if gene:
            assert self in gene, \
                f"{self} outside of {gene}"
        self.gene: Gene|None = gene
        

        #If preceeding exon was passed
        if prev_exon:
            #Validate input
            assert prev_exon <= self, \
                f"Exon-intron misalignment: {prev_exon} & {self}"
            
            #Link this exon to that one
            prev_exon.next_intron = self

        #If following exon was passed
        if next_exon:
            #Validate input
            assert self <= next_exon, \
                f"Intron-exon misalignment: {self} & {next_exon}"
            #Link this exon to that one
            next_exon.prev_intron = self
        
        #Prospectively links exons to this one
        self.prev_exon: Exon|None = prev_exon
        self.next_exon: Exon|None = next_exon
        
        
        #Attributes whose specific values are to be computed later
        self.variants:    list[Intron]         = []
        self.c_score:     float|None           = None
        self.nc_score:    float|None           = None
        self.unif_score:  float|None           = None
        self.splice_site: str|None             = None
        self.traits:      dict[str:bool|float] = {}
    
    
    def get_variants(self, min_exon_len: int):
        """
        Find possible variants of an intron by examining the nucleotides at exon-intron boundaries
        and add all of them to self.variants
        """
        #Get surrounding exons as local variables, to save on typing
        prev_exon = self.prev_exon
        next_exon = self.next_exon
        
        
        #Check presence of flanking exons, linkage & coords
        assert prev_exon, \
            f"{self} does not have a preceeding exon"
        assert next_exon, \
            f"{self} does not have a following exon"
        assert prev_exon @ next_exon, \
            f"{prev_exon} & {next_exon} are unlinked"
        assert prev_exon @ self, \
            f"{prev_exon} & {self} are unlinked"
        assert self @ next_exon, \
            f"{self} & {next_exon} are unlinked"
        assert prev_exon < next_exon, \
            f"{prev_exon} & {next_exon} are misordered"
        assert prev_exon <= self, \
            f"{prev_exon} & {self} are misordered/unflush"
        assert self <= next_exon, \
            f"{self} & {next_exon} are misordered/unflush"
        assert self.seq, \
            f"{self.seq} has no sequence"
        assert prev_exon.seq, \
            f"{prev_exon.seq} has no sequence"
        assert next_exon.seq, \
            f"{next_exon.seq} has no sequence"
        
        
        #prev_exon & next_exon are previous & next in the coordinate sense, which, as long as
        #the gene in on the positive strand, is also their within-gene order.
        #If the gene is on the negative strand, this is opposite to within-gene order, however,
        #which would cause variant detection to yield incorrect results.
        #To remedy this, if the gene is on the negative strand, the sequences of the intron & exons
        #are temporarily reversed.
        if self.strand == '-':
            #Backup correctly-oriented sequences
            self.seq_      = self.seq
            prev_exon.seq_ = prev_exon.seq
            next_exon.seq_ = next_exon.seq
            #Reverse the sequences
            self.seq      = self.seq[::-1]
            prev_exon.seq = prev_exon.seq[::-1]
            next_exon.seq = next_exon.seq[::-1]
        
        
        #Project what the intron & exons would look like before & after splicing 
        assert ( pre_mrna_seq := prev_exon + self + next_exon )
        assert ( mrna_seq     := prev_exon + next_exon )
        #The previous asserts already checked that the exons & intron have sequences, so these
        #asserts should always pass
        #The actual purpose of this stanza is to get the values of `pre_mrna_seq' & `mrna_seq'
        #These two variables are used exclusively in asserts, so they don't need to be initialized
        #if asserts are disabled - wrapping the assignments in asserts achieves that
        

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
        

        #Iterate over all possible shift values
        for shift in range(-suff_len, pref_len+1):
            #Skip a shift value of 0, since it means no shift
            if shift == 0:
                continue
            
            
            #Get shifted intron coordinates
            #These are also the coordinates for the end of the previous exon, and the start of the
            #next exon, respectively
            new_intron_start = self.start + shift
            new_intron_end   = self.end   + shift
            
            
            #Instantiate the variant's flanking exons
            new_prev_exon = Exon(
                self.scaffold, self.prev_exon.start, new_intron_start, strand = self.strand,
                gene = self.gene
            )
            new_next_exon = Exon(
                self.scaffold, new_intron_end, self.next_exon.end, strand = self.strand,
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
            assert len(new_prev_exon) == new_prev_exon.end - new_prev_exon.start, \
                f"{new_prev_exon} length mismatch: {len(new_prev_exon)} by seq"
            assert len(new_intron)    == new_intron.end    - new_intron.start, \
                f"{new_intron} length mismatch: {len(new_intron)} by seq"
            assert len(new_next_exon) == new_next_exon.end - new_next_exon.start, \
                f"{new_next_exon} length mismatch: {len(new_next_exon)} by seq"
            
            
            #Check whether the shift leaves the unspliced sequence untouched
            assert new_prev_exon.seq + new_intron.seq + new_next_exon.seq == pre_mrna_seq, \
                f"Variant {shift} of {self} yields different pre-mRNA"
            
            #Likewise for the spliced sequence
            assert new_prev_exon.seq + new_next_exon.seq == mrna_seq, \
                f"Variant {shift} of {self} yields different mRNA"
            

            #Check alignment & linkage of exons & intron
            assert new_prev_exon < new_next_exon, \
                f"{new_prev_exon} & {new_next_exon} out of order"
            assert new_prev_exon <= new_intron, \
                f"{new_prev_exon} & {new_intron} out of order"
            assert new_intron <= new_next_exon, \
                f"{new_intron} & {new_next_exon} out of order"
            assert new_prev_exon @ new_next_exon, \
                f"{new_prev_exon} & {new_next_exon} unlinked"
            assert new_prev_exon @ new_intron, \
                f"{new_prev_exon} & {new_intron} unlinked"
            assert new_intron @ new_next_exon, \
                f"{new_intron} & {new_next_exon} unlinked"
            
            
            #Undo the new exons' & intron's sequences being reversed if needed
            if self.strand == '-':
                new_intron.seq    = new_intron.seq[::-1]
                new_prev_exon.seq = new_prev_exon.seq[::-1]
                new_next_exon.seq = new_next_exon.seq[::-1]
            
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
    
    
    def add_traits(self, weighted: bool = True, pairing_len: int = 21):
        """
        Add attributes describing the intron's conventional and nonconventional traits
        If `weighted' is True, the pairing scores will be weighted
        `pairing_len' is the max number of nucleotides to examine during pairing score calculation
        """
        #Get surrounding exons as local variables, to save on typing
        prev_exon = self.prev_exon
        next_exon = self.next_exon
        

        #Check presence of flanking exons, linkage, coords & presence of sequences
        assert prev_exon, \
            f"{self} does not have a preceeding exon"
        assert next_exon, \
            f"{self} does not have a following exon"
        assert prev_exon @ next_exon, \
            f"{prev_exon} & {next_exon} are unlinked"
        assert prev_exon @ self, \
            f"{prev_exon} & {self} are unlinked"
        assert self @ next_exon, \
            f"{self} & {next_exon} are unlinked"
        assert prev_exon < next_exon, \
            f"{prev_exon} & {next_exon} are misordered"
        assert prev_exon <= self, \
            f"{prev_exon} & {self} are misordered/unflush"
        assert self <= next_exon, \
            f"{self} & {next_exon} are misordered/unflush"
        assert self.seq, \
            f"{self.seq} has no sequence"
        assert prev_exon.seq, \
            f"{prev_exon.seq} has no sequence"
        assert next_exon.seq, \
            f"{next_exon.seq} has no sequence"
        
        
        #prev_exon & next_exon are previous & next in the coordinate sense, which, as long as
        #the gene in on the positive strand, is also their within-gene order.
        #If the gene is on the negative strand, this is opposite to within-gene order, and so
        #the exons need to be switched in that case.
        if self.strand == '-':
            prev_exon, next_exon = next_exon, prev_exon
        
        
        #Instantiate trait dictionary
        #Traits are bools (informing if a given trait is present), with the exception of pairing
        #scores, which are floats
        self.traits: dict[str:bool|float] = {}
        
        
        #Last nucleotide of previous exon is a pyrimidine
        self.traits["prev_exon_y"] = prev_exon[-1] in "CT"
        
        #First nucleotide of intron is a purine
        self.traits["intron_r"]    = self[0]       in "AG"
        
        #Last nucleotide of intron is a pyrimidine
        self.traits["intron_y"]    = self[-1]      in "CT"
        
        #First nucleotide of next exon is a purine
        self.traits["next_exon_r"] = next_exon[0]  in "AG"
        
        
        #CAG & CTG at appropriate positions
        self.traits["intron_cag"] = len(self) >= 6 and self[3:6]   == "CAG"
        self.traits["intron_ctg"] = len(self) >= 8 and self[-8:-5] == "CTG"
        
        
        #Add splice site; will be None, if the intron is < 4 nucleotides
        #The splice site is its own attribute, not a trait
        self.splice_site: str = self[:2] + self[-2:] if len(self) >= 4 else None
        
        #Intron has a conventional splice site
        self.traits["ss_is_conv"] = self.splice_site in { "GTAG", "GCAG", "CTAC", "CTGC" }
        
        
        
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
    
    
    class PairScoresReport(Enum):
        """
        Enum used to control how pairing scores should be reported by the following method, to be
        used in intron scoring. Does not affect actual calculation of pairing scores, only how they
        are used in intron scoring.
        """
        #All pairing scores are reported as-is; default behavior
        KEEP_ALL = 0
        #Pairing scores are all reported as 0.0
        DROP_ALL = 1
        #All three pairing scores are reported to be the value of the single highest pairing score
        PROPAGATE_BEST = 2
        #The single highest pairing score is reported as-is; the other two are reported as 0.0
        #In case of ties, the precedence is pair_score_2 > pair_score_1 > pair_score_3
        KEEP_BEST = 3

        #Let's say the pairing scores were calculated to be:
        #pair_score_1: 1.0, pair_score_2: 3.0, pair_score_3: 2.0
        #
        #For KEEP_ALL, they will be reported as follows:
        #pair_score_1: 1.0, pair_score_2: 3.0, pair_score_3: 2.0
        #
        #For DROP_ALL, they will be reported as follows:
        #pair_score_1: 0.0, pair_score_2: 0.0, pair_score_3: 0.0
        #
        #For PROPAGATE_BEST, they will be reported as follows:
        #pair_score_1: 3.0, pair_score_2: 3.0, pair_score_3: 3.0
        #
        #For KEEP_BEST, they will be reported as follows:
        #pair_score_1: 0.0, pair_score_2: 3.0, pair_score_3: 0.0
    
    
    
    def get_ml_traits(self, pair_scores: PairScoresReport = PairScoresReport.KEEP_ALL) -> array:
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
        The first 5 traits are booleans encoded as float (1.0: True, 0.0: False),
        the latter 3 traits are floats.

        If `pair_scores' is False, the last three traits are forced to be 0.0, regardless of
        the actual pairing scores' values.
        """
        #Check that traits are actually present
        assert self.traits, \
            f"{self} does not have traits computed"
        
        
        #Return dummy traits for very short introns
        if len(self) < 6:
            return array([0.0] * 8)
        
        
        #Get the appropriate pairing score values to report
        match pair_scores:
            #Report all pairing scores as 0.0
            case Intron.PairScoresReport.DROP_ALL:
                ps1 = 0.0
                ps2 = 0.0
                ps3 = 0.0
            
            #Use best pairing score as all pairing scores
            case Intron.PairScoresReport.PROPAGATE_BEST:
                max_ps = max(self.traits["pair_score_1"],
                             self.traits["pair_score_2"],
                             self.traits["pair_score_3"])
                ps1 = max_ps
                ps2 = max_ps
                ps3 = max_ps
            
            #Report best pairing score as-is, and the other ones as 0.0
            case Intron.PairScoresReport.KEEP_BEST:
                #Get best pairing score
                max_ps = max(self.traits["pair_score_1"],
                             self.traits["pair_score_2"],
                             self.traits["pair_score_3"])
                
                #Report best pairing score, suppress other ones
                if max_ps == self.traits["pair_score_2"]:
                    ps1 = 0.0
                    ps2 = max_ps
                    ps3 = 0.0
                elif max_ps == self.traits["pair_score_1"]:
                    ps1 = max_ps
                    ps2 = 0.0
                    ps3 = 0.0
                else:
                    ps1 = 0.0
                    ps2 = 0.0
                    ps3 = max_ps

            #By default, report all pairing scores as-is (KEEP_ALL)
            case _:
                ps1 = self.traits["pair_score_1"]
                ps2 = self.traits["pair_score_2"]
                ps3 = self.traits["pair_score_3"]
        
        
        #Extract the traits in the appropriate order, converting the bool traits to floats
        return array([
            self.traits["prev_exon_y"],
            self.traits["intron_r"],
            self.traits["intron_cag"] and self.traits["intron_ctg"],
            self.traits["intron_y"],
            self.traits["next_exon_r"],
            ps2,
            ps1,
            ps3
        ])
    
    
    def write_stats(self, fd: TextIO, idx: int = 0):
        """
        Write statistics for an intron & for all of its variants to a file. The statistics are:
         1. the name of the gene of the intron
         2. the name of the transcript of the intron
         3. the index of the intron within the gene (taken from `idx')
         4. the number of variants the intron has
         5. the rank of the variant (the base intron is rank #0, variants are ranked by score)
         6. the scaffold
         7. the strand
         8. the start coordinate
         9. the end coordinate
        10. unified score
        11. conv score
        12. nonconv score
        13. the difference between the conv & nonconv score
        14. the difference between the unified scores of the base intron and the best variant
            (base intron only)
        15. splice site
        16. whether the splice site is conventional
        17. last 5/first 10 nucleotides of the preceding exon and the intron/variant, respectively
        18. last 10/first 10 nucleotides of the intron/variant and the following exon, respectively
        19. whether the last nucleotide of the preceding exon is a pyrimidine
        20. whether the first nucleotide is a purine
        21. whether the last nucleotide is a pyrimidine
        22. whether the first nucleotide of the preceding exon is a purine
        23. whether the characteristic CAG & CTG are present at the approptiate positions
        24. whether the fourth and sixth-to-last nucleotide pair up
        25. whether the fifth and seventh-to-last nucleotide pair up
        26. whether the sixth and eighth-to-last nucleotide pair up
        """
        #Sort variants by descending score
        self.variants.sort(key = lambda v: (v.unif_score, v.c_score, v.nc_score), reverse=True)
        
        
        #Get index of intron as a string
        idx: str = str(idx)
        
        
        #Build list containing the base intron and all variants for the purpose of other asserts
        assert ( eachvar := [self] + self.variants )
        
        #Validate that the intron and all the variants all have the required information
        assert all( var.gene                 for var in eachvar), \
            f"{self} (or variant(s) thereof) lacks parent gene"
        assert all( var.gene.name            for var in eachvar), \
            f"{self} (or variant(s) thereof) lacks parent gene's name"
        assert all( var.gene.transcript      for var in eachvar), \
            f"{self} (or variant(s) thereof) lacks parent transcript"
        assert all( var.gene.transcript.name for var in eachvar), \
            f"{self} (or variant(s) thereof) lacks parent transcript's name"
        assert all( var.unif_score != None   for var in eachvar), \
            f"{self} (or variant(s) thereof) lacks unified score"
        assert all( var.c_score    != None   for var in eachvar), \
            f"{self} (or variant(s) thereof) lacks conventional score"
        assert all( var.nc_score   != None   for var in eachvar), \
            f"{self} (or variant(s) thereof) lacks nonconventional score"
        assert all( var.traits               for var in eachvar), \
            f"{self} (or variant(s) thereof) lacks computed traits"
        assert all( var.prev_exon            for var in eachvar), \
            f"{self} (or variant(s) thereof) lacks preceding exon"
        assert all( var.next_exon            for var in eachvar), \
            f"{self} (or variant(s) thereof) lacks following exon"
        assert all( var.seq                  for var in eachvar), \
            f"{self} (or variant(s) thereof) lacks sequence"
        assert all( var.prev_exon.seq        for var in eachvar), \
            f"{self} (or variant(s) thereof) lacks sequence of preceding exon"
        assert all( var.next_exon.seq        for var in eachvar), \
            f"{self} (or variant(s) thereof) lacks sequence of following exon"
        
        
        #Check whether the inton is on the positive strand
        pos: bool = self.strand == '+'
        
        
        #Write statistics for the base intron
        fd.write(
            '\t'.join([
                self.gene.name,
                self.gene.transcript.name,
                idx,
                str(len(self.variants)),
                '0',
                self.scaffold,
                self.strand,
                str(self.start+1),
                str(self.end),
                str(self.unif_score),
                str(self.c_score),
                str(self.nc_score),
                str( abs(self.c_score - self.nc_score) ),
                str( abs(self.unif_score - self.variants[0].unif_score) ) if self.variants else "",
                self.splice_site or "too short",
                str(self.traits["ss_is_conv"]),
                ( self.prev_exon if pos else self.next_exon )[-5:],
                self[:10],
                self[-10:],
                ( self.next_exon if pos else self.prev_exon )[:5],
                str(self.traits["prev_exon_y"]),
                str(self.traits["intron_r"]),
                str(self.traits["intron_y"]),
                str(self.traits["next_exon_r"]),
                str(self.traits["intron_cag"] and self.traits["intron_ctg"]),
                str(self.traits["pair_score_1"]),
                str(self.traits["pair_score_2"]),
                str(self.traits["pair_score_3"]),
                str(self.traits["pair_3_6"]),
                str(self.traits["pair_4_7"]),
                str(self.traits["pair_5_8"])
            ]) + '\n'
        )
        
        
        #Write statistics for each variant
        for n, var in enumerate(self.variants, 1):
            fd.write(
                '\t'.join([
                    self.gene.name,
                    self.gene.transcript.name,
                    idx,
                    str(len(self.variants)),
                    str(n),
                    self.scaffold,
                    self.strand,
                    str(var.start+1),
                    str(var.end),
                    str(var.unif_score),
                    str(var.c_score),
                    str(var.nc_score),
                    str( abs(var.c_score - var.nc_score) ),
                    "",
                    var.splice_site or "too short",
                    str(var.traits["ss_is_conv"]),
                    ( var.prev_exon if pos else var.next_exon )[-5:],
                    var[:10],
                    var[-10:],
                    ( var.next_exon if pos else var.prev_exon )[:5],
                    str(var.traits["prev_exon_y"]),
                    str(var.traits["intron_r"]),
                    str(var.traits["intron_y"]),
                    str(var.traits["next_exon_r"]),
                    str(var.traits["intron_cag"] and var.traits["intron_ctg"]),
                    str(var.traits["pair_score_1"]),
                    str(var.traits["pair_score_2"]),
                    str(var.traits["pair_score_3"]),
                    str(var.traits["pair_3_6"]),
                    str(var.traits["pair_4_7"]),
                    str(var.traits["pair_5_8"])
                ]) + '\n'
            )




###################################################################################################
#   DESERIALIZATION
###################################################################################################

def deserialize_gff(gff: str, invert: bool = False) -> dict[str:Gene]:
    def get_id_and_parent(attr: str) -> tuple[str, str]:
        """
        Given an attribute field, extract the "ID" and "Parent" fields, in that order.
        Returns an empty string if a given field is absent from the attribute field.
        If `invert' is true, the strand field's value will be inverted for all processed records.
        """
        return attr.partition("ID=")[2].partition(';')[0], \
               attr.partition("Parent=")[2].partition(';')[0]
    
    
    fd = open(gff)
    
    
    #Dict of created genes
    genes: dict[str:Gene] = {}
    #List containing info about exons to be created after parsing the whole file
    deferred_exons: list[tuple] = []
    #Table for converting from transcript IDs to gene IDs
    mrna_to_gene: dict[str:str] = {}
    
    
    #Iterate over GFF records
    for lineno, line in enumerate(fd,1):
        #Skip comments & pragmas
        if line.startswith('#'):
            continue
        
        line = line.strip()
        fields: list[str] = line.split('\t')
        
        assert len(fields) == 9, \
            f"Line {lineno} of GFF file does not have 9 fields:\n{line}"
        
        
        scaffold: str = fields[0]
        ft:       str = fields[2]
        start:    int = int(fields[3])-1
        end:      int = int(fields[4])
        strand:   str = fields[6]
        attr:     str = fields[8]
        
        
        if strand not in '-+':
            continue
        
        
        #Invert strand if requested
        if invert:
            strand = '-' if strand == '+' else '+'
        
        
        if ft == "gene":
            #Get name of gene (from its ID) and instantiate it
            gene_name: str = get_id_and_parent(attr)[0]

            assert gene_name, \
                f"Gene @ line {lineno} of GFF has no ID:\n{line}"
            
            #Instantiate new gene
            genes[gene_name] = Gene(scaffold, start, end, name=gene_name, strand=strand, exons=[])
        

        elif ft == "exon":
            #Get name of parent transcript
            mrna_name: str = get_id_and_parent(attr)[1]

            assert mrna_name, \
                f"Exon @ line {lineno} of GFF has no Parent:\n{line}"
            
            #Keep track of this exon to add to a gene later
            deferred_exons.append((mrna_name, scaffold, start, end, strand))
        

        elif ft == "mRNA":
            #Get name of transcript & its parent gene
            mrna_name, gene_name = get_id_and_parent(attr)

            assert mrna_name, \
                f"Transcript @ line {lineno} of GFF has no ID:\n{line}"
            assert gene_name, \
                f"Transcript @ line {lineno} of GFF has no Parent:\n{line}"

            #Submit to conversion table
            mrna_to_gene[mrna_name] = gene_name
    
    
    #Instantiate all exons and link them to their parent gene
    for mrna_name, scaffold, start, end, strand in deferred_exons:
        #Get the name of the exon's grandparent gene
        assert mrna_name in mrna_to_gene, \
            f"GFF file contains an exon referencing a missing transcript '{gene_name}'"
        gene_name: str = mrna_to_gene[mrna_name]
        assert gene_name in genes, \
            f"GFF file contains transcript {mrna_name}, referencing a missing gene '{gene_name}'"
        
        
        #Get the grandparent gene
        gene: Gene = genes[gene_name]
        
        #Instantiate the exon and add it to the gene's exon list
        gene.exons.append(
            Exon(
                scaffold, start, end, strand=strand, gene=gene
            )
        )
    
    
    #Submit the name of each transcript to its gene
    for mrna_name, gene_name in mrna_to_gene.items():
        assert gene_name in genes, \
            f"GFF file contains transcript {mrna_name}, referencing a missing gene '{gene_name}'"
        
        genes[gene_name].transcript.name = mrna_name
    
    
    fd.close()
    return genes


def deserialize_fasta(fasta: str) -> dict[str:str]:
    """
    Given a path to a FASTA file, deserialize it to a dictionary.
    Keys are sequence names, values are sequences.
    """
    with open(fasta) as fd:
        return dict( SimpleFastaParser(fd) )


def load_model(model: str):
    """
    Load a pickled sklearn model from filename
    """
    with open(model, 'rb') as fd:
        return load(fd)




###################################################################################################
#   INTRON SCORING
###################################################################################################

#Nucleotide pair weights.
#Each pair of nucleotides which actually pair up in pre-mRNA have the following weights:
#A-T: 0.5, C-G: 1.0, G-T: 0.375.
#Non-pairing nucleotide pairs have an implicit weight of 0.0.
#For commutativity, each pairing pair is present in this dictionary forward & backwards.
PAIR_WEIGHTS: dict[tuple[str,str]:float] = {
    ('A','T'): 0.5,   ('T','A'): 0.5,
    ('G','C'): 1.0,   ('C','G'): 1.0,
    ('G','T'): 0.375, ('T','G'): 0.375
}


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


def score_introns(genes: Iterable[Gene], nc_model, c_model, *,
                  unif_score_from_ss: bool = False,
                  batch_size: int = 0,
                  weighted: bool = True,
                  pair_scores: Intron.PairScoresReport = Intron.PairScoresReport.KEEP_ALL):
    """
    Calculate conventional & nonconventional structure compatibility scores for all introns
    (& variants) of the passed genes, using the supplied models.

    Each scored intron is given a conventionality score and a nonconventionality score, plus a
    unified score, which is the higher score of the two.

    Alternatively, if `unif_by_ss' is True, introns with conventional splice sites are instead
    forced to have the maximum possible unified score, and introns without copy the unified score
    from their nonconventional score.

    By default, all introns are scored at once, but if number of introns is very large, an OOM
    error might occur. To remedy this, the introns can be processed in batches by passing a
    positive integer to `batch_size'. This value will be the size of the intron batches.
    
    `weighted' controls whether the introns' pairing scores are weighted or unweighted.
    `pair_scores' controls how the pairing scores are used in intron scoring.
    The pairing scores are calculated regardless of the value `pair_scores' - the parameter
    controls only how the pairing scores are used in intron scoring.
    """
    #Build list of introns to assess
    phase("Score introns: gather")
    introns: list[Introns] = []
    for gene in genes:
        for intron in gene.introns:
            #Include the intron and all its variants in assessment
            introns.append(intron)
            introns += intron.variants
    
    assert len(introns) > 0, \
        "Intron gathering went wrong"
    
    
    #Compute traits for each intron
    phase("Score introns: compute traits")
    for intron in introns:
        intron.add_traits(weighted)
    
    
    #Retrieve ML traits and splice site
    phase("Score introns: retrieve traits")
    ml_traits: list[array] =  [ intron.get_ml_traits(pair_scores) for intron in introns ]
    
    #Check number of trait sets
    assert len(introns) == len(ml_traits), \
        "Intron trait gathering went wrong"
    
    
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
    assert len(introns) == len(nc_scores), \
        "Nonconv scoring went wrong"
    
    
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
    assert len(introns) == len(c_scores), \
        "Conv scoring went wrong"
    
    
    #Submit each intron's scores to the actual object
    phase("Score introns: save scores")
    for intron, c_score, nc_score in zip(introns, c_scores, nc_scores):
        #Per-class score
        intron.c_score  = float(c_score)
        intron.nc_score = float(nc_score)
        #Unfied score
        if unif_score_from_ss:
            intron.unif_score = 1.0 if intron.traits["ss_is_conv"] else float(nc_score)
        else:
            intron.unif_score = float(max( c_score, nc_score ))




###################################################################################################
#   MISCELLANEOUS
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
    #If this is the end of a phase, end it by reporting its duration
    if hasattr(phase, "timer"):
        #Move the printer cursor one line up, to the 40th character of the line
        eprint(f"\033[0;A\033[40;C{time() - phase.timer:.3f} sec")
    
    #If a new_phase was requested, start it by printing its name
    if new_phase != None:
        eprint(new_phase.ljust(40))
    
    #Reset the timer
    phase.timer = time()


def concat_genseq(iterable: Iterable[Type[GenomicSequence]], sep: str = '') -> str:
    """
    Given an iterable of GenomicSequence (and/or subclass) objects, extract the sequences of each
    one and concatenate them, with an optional separator between each sequence
    """
    return sep.join( map( attrgetter("seq"), iterable ) )

