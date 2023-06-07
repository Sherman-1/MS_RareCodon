from Bio.Seq import Seq 
from collections import OrderedDict as OD



class GenomicFeature:

    def __init__(self, ID : str, start : int, end : int):
        self.ID = ID
        self.start = start - 1
        self.end = end - 1
        self.counter = 0


class Gene(GenomicFeature):

    counter = 0

    def __init__(self, ID, start, end, chromosome : str,  multi : bool, sense : str, frame : int):
        super().__init__(ID, start, end)
        self.chromosome = chromosome
        self.frame = frame
        self.sense = sense
        self.multi = multi
        self.aORFs = OD()
        self._exons = OD()
        Gene.counter += 1
        

    @property
    def exons(self):
        return self._exons
    
    @property
    def exons_list(self):
        return list(self._exons.values())

    @property
    def orfs_list(self):
        return list(self.aORFs.values())

    @exons.setter
    def exons(self, value):
    
        self._exons = value

    def add_exon(self, key, value):
        
        if key in self._exons:
            raise KeyError(f'exon {key} already exists')

        elif type(value) != Exon:

            raise TypeError('exon must be an Exon')
        
        else:
            self._exons[key] = value

            
    def add_orf(self, key, value):

        if type(value) != Orf:
            raise TypeError('value must be Orf')

        elif key in self.aORFs:
            raise KeyError(f'orf {key} already exists')
        else:
            self.aORFs[key] = value

    
    def sort(self):

        if self.sense == "+":

            if int(self.start) > int(self.end) : self.start, self.end = self.end, self.start
            self._exons = OD(sorted(self._exons.items(), key=lambda x: x[1].start))
            

        elif self.sense == "-":

            if int(self.end) > int(self.start) : self.start, self.end = self.end, self.start
            self._exons = OD(sorted(self._exons.items(), key=lambda x: x[1].start, reverse=True))

    def sort_aorfs(self):

        if self.sense == "+":

            self.aORFs = OD(sorted(self.aORFs.items(), key=lambda x: x[1].start))

        elif self.sense == "-":

            self.aORFs = OD(sorted(self.aORFs.items(), key=lambda x: x[1].start, reverse=True))

    def get_adjacent_nucleotides(self, sequence, nb_codons, nb_nt):

        sequence = Seq(sequence)

        for orf in self.orfs_list:

            if orf.ribostartLocalisation == "exon":

                
                position = orf.ribostart
                strand = self.sense
                codons_upstream = []
                codons_downstream = []
                
                if strand == "+":

                    # Get upstream codons
                    for i in range(nb_codons):
                        start = position -1 - 3 * (i + 1)
                        end = position - 1 - 3 * i
                        if start < 0:  # check if start of sequence is reached
                            break
                        codon = sequence[start:end]
                        codons_upstream.append(codon)
                    
                    # Get downstream codons
                    for j in range(nb_codons):
                        start = position + 2 + 3 * j 
                        end = position + 2 + 3 * (j + 1)
                        if end > len(sequence):  # check if end of sequence is reached
                            break
                        codon = sequence[start:end]
                        codons_downstream.append(codon) 

                elif strand == "-":

                    # Get downstream codons
                    for i in range(nb_codons):
                        start = position - 2 - 3 * ( i + 1 )
                        end = position - 2 - 3 * i
                        codon = str(Seq(sequence[start:end]).reverse_complement())
                        codons_downstream.append(codon) 

                    # Get upstream codons
                    for j in range(nb_codons):
                        start = position + 1 + 3 * j
                        end = position + 1 + 3 * ( j + 1 )
                        codon = str(Seq(sequence[start:end]).reverse_complement())
                        codons_upstream.append(codon)


                orf.upstream = codons_upstream
                orf.downstream = codons_downstream
                    
            else: 
                
                orf.upstream = "NA"
                orf.downstream = "NA"
                



class Exon(GenomicFeature):

    def __init__(self, ID, start, end, gene : Gene, abs_frame : int):
        super().__init__(ID, start, end )
        self.abs_frame = abs_frame
        self.gene = gene

        if self.gene.sense == "-":
            self.start, self.end = self.end, self.start
 
class Orf(GenomicFeature):

    def __init__(self, ID, start, end, gene, ribospike : int):
        super().__init__(ID, start, end)
        self.gene = gene
        self.ribospike = ribospike

        if self.gene.sense == "-":
            self.start, self.end = self.end, self.start
            self.ribostart = int(self.start - 3*self.ribospike)
        else:
            self.ribostart = int(self.start + 3*self.ribospike)
            
        self.upstream = None
        self.downstream = None
        self.frame = None

        self.ribostartLocalisation = None
        self.exon = None
        
    def locRiboStart(self):

        if self.gene.sense == "+":
            for exon in self.gene.exons_list:
                if exon.start <= self.ribostart <= exon.end:
                    self.ribostartLocalisation = "exon"
                    self.exon = exon
                    self.frame = abs(self.exon.start - self.ribostart) % 3
                    return
                
            if self.gene.start <= self.ribostart < self.gene.exons_list[0].start:
                self.ribostartLocalisation = "5UTR"
                return
            
            elif self.gene.exons_list[-1].end < self.ribostart <= self.gene.end:
                self.ribostartLocalisation = "3UTR"
                return

            elif self.ribostart < self.gene.start:
                self.ribostartLocalisation = "upstream"
                return
            
            elif self.ribostart > self.gene.end:
                self.ribostartLocalisation = "downstream"
                return
            
            else:
                self.ribostartLocalisation = "intron"
                return
            
        elif self.gene.sense == "-":
            for exon in self.gene.exons_list:
                if exon.end <= self.ribostart <= exon.start:
                    self.exon = exon
                    self.ribostartLocalisation = "exon"
                    self.frame = abs(self.exon.start - self.ribostart) % 3

                    return
                
            if self.gene.exons_list[0].start < self.ribostart <= self.gene.start:
                self.ribostartLocalisation = "5UTR"
                return
            
            elif self.gene.end <= self.ribostart < self.gene.exons_list[-1].end:
                self.ribostartLocalisation = "3UTR"
                return

            elif self.ribostart > self.gene.start:
                self.ribostartLocalisation = "upstream"
                return
            
            elif self.ribostart < self.gene.end:
                self.ribostartLocalisation = "downstream"
                return
            
            else:
                self.ribostartLocalisation = "intron"
                return
            

class GeneStructureError(Exception):
    def __init__(self, message):
        self.message = message
        super().__init__(self.message)
