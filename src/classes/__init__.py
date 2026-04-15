"""
Core data model classes for genomic sequences and their components.
"""

from .genomic_sequence import GenomicSequence
from .gene import Gene
from .intron import Intron
from .exon import Exon
from .transcript import Transcript

__all__ = [
    "GenomicSequence",
    "Gene", 
    "Intron",
    "Exon",
    "Transcript"
]