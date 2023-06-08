from Bio.Seq import Seq 
from collections import OrderedDict as OD



class GenomicFeature:

    def __init__(self, ID : str, start : int, end : int):
        self.ID = ID
        self.start = start
        self.end = end
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

    def get_adjacent_codons(self, sequence, nb_codons): # Deprecated for now

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
                        start = position - 3 - 3 * ( i + 1 )
                        end = position - 3 - 3 * i
                        codon = str(Seq(sequence[start:end]).reverse_complement())
                        codons_downstream.append(codon) 

                    # Get upstream codons
                    for j in range(nb_codons):
                        start = position  + 3 * j
                        end = position + 3 * ( j + 1 )
                        codon = str(Seq(sequence[start:end]).reverse_complement())
                        codons_upstream.append(codon)


                orf.upstream = codons_upstream
                orf.downstream = codons_downstream

    def get_adjacent_nucleotides(self, sequence, nb_nt : int):

        for orf in self.orfs_list:

            position = int(orf.ribostart - 1) # -1 because python starts at 0
            strand = self.sense

            if strand == "+":

                # Get upstream codons
                start = position - nb_nt 
                end = position 
                if start < 0:  # check if start of sequence is reached
                    codons_upstream = sequence["seq"][:end]

                else:
                    codons_upstream = sequence["seq"][start:end]

                # Get downstream codons
                start = position + 3  
                end = position + 3  + nb_nt
                if end > sequence["len"]: 
                    codons_downstream = sequence["seq"][start:]
                else:
                    codons_downstream = sequence["seq"][start:end]

            elif strand == "-":

                # Get upstream codons
                start = position + 1  
                end = position + 1  + nb_nt
                if end > sequence["len"]:
                    codons_upstream = Seq(sequence["seq"][start:]).reverse_complement()
                else:
                    codons_upstream = Seq(sequence["seq"][start:end]).reverse_complement()

                # Get downstream codons
                start = position - 2 - nb_nt
                end = position - 2 
                if start < 0:  # check if start of sequence is reached
                    codons_downstream = Seq(sequence["seq"][:end]).reverse_complement() 
                else:
                    codons_downstream = Seq(sequence["seq"][start:end]).reverse_complement()

            orf.upstream = str(codons_upstream)
            orf.downstream = str(codons_downstream)

class Exon(GenomicFeature):

    def __init__(self, ID, start, end, gene : Gene, abs_frame : int):
        super().__init__(ID, start, end )
        self.abs_frame = abs_frame
        self.gene = gene

        if self.gene.sense == "-":
            self.start, self.end = self.end, self.start
 
class Orf(GenomicFeature):

    counter = 0
    # Allows to retrieve instance attributes from class 
    items = ["gene","ribospike", "MSMS", "MS_sds", "MS_cond", "age_rel", 
             "upstream", "downstream", "rel_frame", "ribostartLocalisation", "exon",
             "start", "end", "ribostart", "ID"]
    def __init__(self, ID, start, end, gene, MSMS, MS_sds, MS_cond, age_rel, start_seq, ribospike : int):
        super().__init__(ID, start, end)
        self.gene = gene
        self.ribospike = ribospike
        self.MSMS = MSMS
        self.MS_sds = MS_sds
        self.MS_cond = MS_cond
        self.age_rel = age_rel
        self.start_seq = start_seq
        Orf.counter += 1
        self.upstream = None
        self.downstream = None
        self.rel_frame = None
        self.ribostartLocalisation = None
        self.exon = None

        if self.gene.sense == "-":
            self.start, self.end = self.end, self.start
            self.ribostart = int(self.start - 3*self.ribospike)
        else:
            self.ribostart = int(self.start + 3*self.ribospike)
            
        
        
    def locRiboStart(self): 

        if self.gene.sense == "+":
            for exon in self.gene.exons_list:
                if exon.start <= self.ribostart <= exon.end:
                    self.ribostartLocalisation = "exon"
                    self.exon = exon
                    self.rel_frame = abs(self.exon.start - self.ribostart) % 3
                    if self.rel_frame == 0:
                        Warning(f"RiboStart {self.ID} is in frame with exon {self.exon.ID}")
                        return 1
                    else:
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
                    self.rel_frame = abs(self.exon.start - self.ribostart) % 3
                    if self.rel_frame == 0:
                        Warning(f"RiboStart {self.ID} is in frame with exon {self.exon.ID}")
                        return 1
                    else:
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
