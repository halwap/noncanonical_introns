from .genomic_sequence import GenomicSequence

class Exon(GenomicSequence):
    """
    Represents an exon within a gene.
    
    Attributes:
        gene: Reference to the parent Gene object
        prev_exon: Reference to the previous exon
        next_exon: Reference to the next exon
        prev_intron: Reference to the previous intron
        next_intron: Reference to the next intron
    """
    def __init__(self, scaffold_name, scaffold_start, scaffold_end, sequence='', strand='',
                 gene=None, prev_exon=None, next_exon=None, prev_intron=None, next_intron=None):
        GenomicSequence.__init__(self, scaffold_name, scaffold_start, scaffold_end,
                                 sequence=sequence, strand=strand)
        if (prev_intron and self.scaffold_name!=prev_intron.scaffold_name) or \
            (next_intron and self.scaffold_name!=next_intron.scaffold_name):
            raise ValueError(f"""[EXON] Neighbouring introns are in different scaffolds:
                             prev i: {prev_intron}, current: {self}, next i: {next_intron}""")
        self.prev_intron = prev_intron
        self.next_intron = next_intron
        self.gene = gene
        if gene:
            self.check_gene()

        if (prev_exon and prev_exon.scaffold_name!=self.scaffold_name) \
            or (next_exon and next_exon.scaffold_name!=self.scaffold_name):
            raise ValueError(f"""[EXON] Neighbouring exons are in different scaffolds:
                             prev e: {prev_exon}, current: {self}, next e: {next_exon}""")
        if prev_exon and (prev_exon.scaffold_end!=scaffold_start+1):
            raise ValueError(f"[EXON] self.prev_exon.scaffold_end={prev_exon.scaffold_end} != self.scaffold_start+1={scaffold_start}+1")
        self.prev_exon = prev_exon
        self.next_exon = next_exon

        if self.next_exon and (self.scaffold_end+1!=self.next_exon.scaffold_start):
            raise ValueError(f"[EXON] self.scaffold_end+1={self.scaffold_end}+1 != \
                self.next_exon.scaffold_start={self.next_exon.scaffold_start}")

    def __str__(self):
        return f"""
    Exon in gene: {self.scaffold_name}
    scaff loc: {self.scaffold_start}-{self.scaffold_end}"""
    
    def check_gene(self):
        if not self.scaffold_name == self.gene.scaffold_name:
            print(f"Intron: {self}, its gene: {self.gene}")
            raise ValueError(f"[EXON] Wrong scaffolds, gene is {self.gene.scaffold_name}, exon is {self.scaffold_name}")
        if not (self.gene.scaffold_start<=self.scaffold_start<self.scaffold_end<=self.gene.scaffold_end):
            print(f"gene {self.gene}, intron {self}")
            raise ValueError(f"""[EXON] Exon not in gene! gene={self.gene}:({self.gene.scaffold_start} {self.gene.scaffold_end}),
                            exon={self}:({self.scaffold_start} {self.scaffold_end})""")
                
    def assign_gene(self, gene):
        self.gene = gene
        self.check_gene()
    
    def merge_exons(self, next_exon_to_merge):
        print("Merging exons", self, next_exon_to_merge)
        self.sequence += next_exon_to_merge.sequence
        self.scaffold_end = next_exon_to_merge.scaffold_end
        self.next_exon = next_exon_to_merge.next_exon
