class GenomicSequence:
    '''
    Base class representing a genomic sequence with coordinates and optional sequence data.
    
    Attributes:
        scaffold_name (str): Name of the scaffold/chromosome
        scaffold_start (int): Start position on scaffold (0-based)
        scaffold_end (int): End position on scaffold (0-based, exclusive)
        sequence (str, optional): The actual sequence string
        strand (str): Strand orientation ('+', '-', or '.')
    '''
    def __init__(self, scaffold_name, scaffold_start, scaffold_end, sequence=None, strand=None):
        self.scaffold_name = scaffold_name
        self.scaffold_start = scaffold_start
        self.scaffold_end = scaffold_end
        if not self.scaffold_start<self.scaffold_end:
            raise ValueError(f"Incorrect start & end coordinates: start={self.scaffold_start}, end={self.scaffold_end} in scaffold {scaffold_name}")
        if sequence and len(sequence) != self.length():
            raise ValueError(f'Incorrect sequence length in {self}: len(seq) = {len(sequence)}, self.length={self.length()}')
        self.sequence = sequence
        if not strand or (strand and strand not in ['+', '-', '.']):
            raise ValueError(f'Strand can only be +, - or ., is {strand}')
        self.strand = strand

    def length(self):
        '''Returns length of the sequence according to the start and end positions.'''
        return self.scaffold_end - self.scaffold_start
    
    def __repr__(self):
        return ' '.join([self.scaffold_name, str(self.scaffold_start), str(self.scaffold_end)])
        
    def __str__(self):
        return self.__repr__()