
from collections import OrderedDict as OD


__all__ = ['gene', 'aORF']


class Gene:

    def __init__(self, ID, chromosome, start, end, multi):

        self.ID = ID
        self.chromosome = chromosome
        self.start = start
        self.end = end
        self.multi = multi
        self.aORFs = OD()
        self.exons = OD()

    def add_ORF(self, key, value):

        if type(value) != OD:
            raise TypeError('aORF must be an OrderedDict')
        
        elif key in self.aORFs:
            raise KeyError(f'aORF {key} already exists')
        
        else:
            self.aORFs[key] = value

    def add_exon(self, key, value):

        if type(value) != OD:
            raise TypeError('exon must be an OrderedDict')
        
        elif key in self.exons:
            raise KeyError(f'exon {key} already exists')
        
        else:
            self.exons[key] = value

class Orf: 

    def __init__(self, ID, start, end, gene ):

        self.ID = ID
        self.start = start
        self.end = end
        self.gene = gene
        self._frame = None
        self._ribospike = None

    @property
    def ribospike(self):
        return self._ribospike
    
    @ribospike.setter
    def ribospike(self, value):
        if type(value) != int:
            raise TypeError('ribospike must be int')
        else:
            self._ribospike = value
            self._frame = self.gene.start + self._ribospike

    @property
    def frame(self):
        return self._frame