"""
Transcript class representing a gene transcript.
"""

class Transcript():
    """
    Represents a transcript sequence derived from a gene.
    
    Attributes:
        scaffold_name (str): Name of the scaffold
        start (int): Start position on scaffold
        end (int): End position on scaffold  
        sequence (str): The transcript sequence
        strand (str): Strand orientation ('+', '-', or '.')
    """
    def __init__(self, scaffold_name, start, end, sequence='', strand=''):
        self.scaffold_name = scaffold_name
        self.start = start
        self.end = end
        self.sequence = sequence
        self.strand = strand