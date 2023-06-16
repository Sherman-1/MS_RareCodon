from Bio.Seq import Seq 
from collections import OrderedDict as OD
from data import MAX_CODONS
from data import NB_NT
from warnings import warn

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

        for i,exon in enumerate(self.exons_list):

                exon.number = i 


    def returnSplicedExons(self, sequence):


            tmpSeq = str()
            for exon in self.exons_list:

                START = exon.start - 1
                END = exon.end - 1

                if self.sense == "+":
                    tmpSeq += sequence[START-exon.abs_frame:END+1]  
                elif self.sense == "-":
                    tmpSeq += str(Seq(sequence[END:START+exon.abs_frame+1]).reverse_complement())
                    
            if len(tmpSeq)%3 != 0:
                warn(f"{self.ID} Multi : {self.multi} Sense : {self.sense}")
            return tmpSeq

    def getAdjTest(self, SeqRecord):
        
        CODONS = MAX_CODONS * 3 + 2 # Enough nucleotides to get MAX_CODONS codons regardless of relative frame 
        splicedExons = self.returnSplicedExons(SeqRecord["seq"])

        for orf in self.orfs_list:

            if orf.ribostartLocalisation == "exon":
                POSITION = orf.distExon
                SPLICED_START = 0
                SPLICED_END = len(splicedExons)

                # upstream sequence
                if POSITION < CODONS:
                    upstream = splicedExons[SPLICED_START:POSITION]
                else:
                    upstream = splicedExons[POSITION-CODONS:POSITION]

                # downstream sequence
                if POSITION+3+CODONS > SPLICED_END:
                    downstream = splicedExons[POSITION+3:]
                else:
                    downstream = splicedExons[POSITION+3:POSITION+3+CODONS]

                if len(upstream) > CODONS or len(downstream) > CODONS:
                    warn(f"upstream : {len(upstream)} downstream : {len(downstream)} for {orf.ID}")

                orf.upstream = upstream
                orf.downstream = downstream
                



    def sort_aorfs(self):

        if self.sense == "+":

            self.aORFs = OD(sorted(self.aORFs.items(), key=lambda x: x[1].start))

        elif self.sense == "-":

            self.aORFs = OD(sorted(self.aORFs.items(), key=lambda x: x[1].start, reverse=True))

    def getAdjacentNucleotides(self, sequence):

        for orf in self.orfs_list:

            # Here we declare the constants that are indexed from 0
            # to avoid confusion with +1/-1 corrections related to 
            # the biological indexing of the sequence
            POSITION = int(orf.ribostart - 1) # -1 for python indexing
             # -1 for python indexing
            STRAND = self.sense
            CODONS = MAX_CODONS * 3 + 2 # This way, we get enough nucleotides to get MAX_CODONS codons regardless of the relative frame of the ribostart

            
            if orf.ribostartLocalisation == "exon":
                EXON_END = orf.exon.end - 1 # -1 for python indexing
                EXON_START = orf.exon.start - 1 # -1 for python indexing

                if STRAND == "+":

                    ## Get upstream nucleotides
                    if abs( EXON_START - POSITION + 1 ) > CODONS:
                        start = POSITION - CODONS
                    else:
                        start = EXON_START
                
                    end = POSITION # Not position - 1 because in seq[a:b] b is not included
                    if end < start:
                        print(f"start : {start}")
                        print(f"end : {end}")
                        print(f"orf.exon.end : {orf.exon.end}")
                        print(f"position : {POSITION}")
                        print(f"orf.ID : {orf.ID}")
                        print(f"codons : {CODONS}")
                        print(f"orf.ribostart : {orf.ribostart}")
                        raise ValueError
                    codons_upstream = sequence["seq"][start:end]
                    if len(codons_upstream) == 0:
                        pass
                    if len(codons_upstream) > CODONS:
                        raise ValueError
                    
                    
                
                    ## Get downstream nucleotides
                    start = POSITION + 3 

                    if abs(EXON_END - start) > CODONS:
                        end = start + CODONS 
                    else:
                        end = EXON_END + 1 # +1 because in seq[a:b] b is not included, -1 for python indexing
                    if end < start:
                        print(f"start : {start}")
                        print(f"end : {end}")
                        print(f"orf.exon.end : {orf.exon.end}")
                        print(f"position : {POSITION}")
                        print(f"orf.ID : {orf.ID}")
                        print(f"codons : {CODONS}")
                        print(f"orf.ribostart : {orf.ribostart}")
                        raise ValueError
                    codons_downstream = sequence["seq"][start:end]
                    if len(codons_downstream) == 0:
                       pass
                    if len(codons_downstream) > CODONS:
                        raise ValueError("Too many downstream nucleotides for orf {orf.ID}")

                elif STRAND == "-":

                    ## Get upstream nucleotides
                    # +1 to get upstreams nucleotides
                    start = POSITION + 1
                    if abs(EXON_START - POSITION + 1 ) > CODONS:
                        end = start + CODONS 
                    else:
                        end = EXON_START + 1 # +1 because in seq[a:b] b is not included


                    # If end < start, print the details of the operations that led to this error
                    if end < start:
                        print(f"start : {start}")
                        print(f"end : {end}")
                        print(f"orf.exon.end : {orf.exon.end}")
                        print(f"position : {POSITION}")
                        print(f"orf.ID : {orf.ID}")
                        print(f"codons : {CODONS}")
                        print(f"orf.ribostart : {orf.ribostart}")
                        raise ValueError


                    codons_upstream = str(Seq(sequence["seq"][start:end]).reverse_complement())
                   
                    if len(codons_upstream) > CODONS:
                        print(f"{orf.ID} length up : {len(codons_upstream)}")

                    ## Get downstream nucleotides
                    end = POSITION - 2 # -2 because in seq[a:b] b is not included
                    if abs(POSITION - 2  - EXON_END) > CODONS:
                        start = POSITION - 2 - CODONS
                    else:
                        start = EXON_END

                    if end < start:
                        print(f"start : {start}")
                        print(f"end : {end}")
                        print(f"orf.exon.end : {orf.exon.end}")
                        print(f"position : {POSITION}")
                        print(f"orf.ID : {orf.ID}")
                        print(f"codons : {CODONS}")
                        print(f"orf.ribostart : {orf.ribostart}")
                        raise ValueError
                    
                    codons_downstream = str(Seq(sequence["seq"][start:end]).reverse_complement())
                    
                    if len(codons_downstream) > CODONS:
                        print(f"{orf.ID} length down : {len(codons_downstream)}")
                    

            else:

                if STRAND == "+":

                    # Get upstream nucleotides
                    # Example : position = 20, NB_NT = 4 => seq[16:20] got us the nt 16, 17, 18, 19 hence 4 nt => good
                    start = POSITION - NB_NT 
                    if start < 0:
                        start = 0
                    end = POSITION # Not position - 1 because in seq[a:b] b is not included
                    codons_upstream = sequence["seq"][start:end]

                    # Get downstream nucleotides
                    start = POSITION + 3
                    end = start + NB_NT 
                    if end > sequence["len"]:
                        end = sequence["len"]+1 # +1 because in seq[a:b] b is not included and in python even if b is out of range, it doesn't raise an error
                    codons_downstream = sequence["seq"][start:end]

                elif STRAND == "-":

                    # Get upstreams nucleotides
                    start = POSITION + 1
                    end = start + NB_NT + 1 # +1 because in seq[a:b] b is not included
                    if end > sequence["len"]:
                        end = sequence["len"]+1
                    codons_upstream = str(Seq(sequence["seq"][start:end]).reverse_complement())

                    # Get downstream nucleotides
                    start = POSITION - 2 - NB_NT
                    if start < 0:
                        start = 0
                    end = POSITION - 2 # -2 because in seq[a:b] b is not included
                    codons_downstream = str(Seq(sequence["seq"][start:end]).reverse_complement())

            orf.upstream = str(codons_upstream)
            orf.downstream = str(codons_downstream)
    
    

     

class Exon(GenomicFeature):

    def __init__(self, ID, start, end, gene : Gene, abs_frame : int):
        super().__init__(ID, start, end )
        self.abs_frame = abs_frame
        self.gene = gene
        self.number = None
        
        if self.gene.sense == "-":
            self.start, self.end = self.end, self.start
        self.length = abs(self.end - self.start) + 1
        if self.gene.sense == "+":
            self.start += self.abs_frame
        elif self.gene.sense == "-":
            self.start -= self.abs_frame

 
class Orf(GenomicFeature):

    counter = 0
    # Allows to retrieve instance attributes from class 
    items = ["gene","ribospike", "MSMS", "MS_sds", "MS_cond", "age_rel", 
             "upstream", "downstream", "rel_frame", "ribostartLocalisation", "exon",
             "start", "end", "ribostart", "ID"]
    def __init__(self, ID, start, end, gene, MSMS, reads, p0, p1, p2,
                 MS_sds, MS_cond, age_rel, start_seq, ribospike : int):
        super().__init__(ID, start, end)
        self.gene = gene
        self.ribospike = ribospike
        self.MSMS = MSMS
        self.MS_sds = MS_sds
        self.MS_cond = MS_cond
        self.age_rel = age_rel
        self.reads = reads
        self.p0 = p0
        self.p1 = p1
        self.p2 = p2
        self.start_seq = start_seq
        Orf.counter += 1
        self.upstream = None
        self.downstream = None
        self.rel_frame = None
        self.ribostartLocalisation = None
        self.exon = None
        self.distExon = None

        if self.gene.sense == "-":
            self.start, self.end = self.end, self.start
            self.ribostart = int(self.start - 3*self.ribospike)
        else:
            self.ribostart = int(self.start + 3*self.ribospike)
        
        self.relRiboStart = self.ribostart 
            
        
        
    def locateRiboStart(self): 

        if self.gene.sense == "+":
            for exon in self.gene.exons_list:
                if exon.start <= self.ribostart <= exon.end:
                    self.exon = exon
                    self.ribostartLocalisation = "exon"
                    self.rel_frame = abs(self.exon.start - self.ribostart) % 3
                    if self.rel_frame == 0:
                        raise ValueError(f"rel_frame = 0 for {self.ID}")
                    return 
                
            if self.gene.start <= self.ribostart < self.gene.exons_list[0].start:
                self.ribostartLocalisation = "5UTR"
                self.rel_frame = 0
                self.exon = "NA"
                return
            
            elif self.gene.exons_list[-1].end < self.ribostart <= self.gene.end:
                self.ribostartLocalisation = "3UTR"
                self.rel_frame = 0
                self.exon = "NA"
                return

            elif self.ribostart < self.gene.start:
                self.ribostartLocalisation = "upstream"
                self.rel_frame = 0
                self.exon = "NA"
                return
            
            elif self.ribostart > self.gene.end:
                self.ribostartLocalisation = "downstream"
                self.rel_frame = 0
                self.exon = "NA"
                return
            
            else:
                self.ribostartLocalisation = "intron"
                self.rel_frame = 0
                self.exon = "NA"
                return
            
        elif self.gene.sense == "-":
            for exon in self.gene.exons_list:
                if exon.end <= self.ribostart <= exon.start:
                    self.exon = exon
                    self.ribostartLocalisation = "exon"
                    self.rel_frame = abs(self.exon.start - self.ribostart) % 3
                    if self.rel_frame == 0:
                        raise ValueError(f"rel_frame = 0 for {self.ID}")
                    return 
                
            if self.gene.exons_list[0].start < self.ribostart <= self.gene.start:
                self.ribostartLocalisation = "5UTR"
                self.rel_frame = 0
                self.exon = "NA"
                return
            
            elif self.gene.end <= self.ribostart < self.gene.exons_list[-1].end:
                self.ribostartLocalisation = "3UTR"
                self.rel_frame = 0
                self.exon = "NA"
                return

            elif self.ribostart > self.gene.start:
                self.ribostartLocalisation = "upstream"
                self.rel_frame = 0
                self.exon = "NA"
                return
            
            elif self.ribostart < self.gene.end:
                self.ribostartLocalisation = "downstream"
                self.rel_frame = 0
                self.exon = "NA"
                return
            
            else:
                self.ribostartLocalisation = "intron"
                self.rel_frame = 0
                self.exon = "NA"
                return
            
    def getExonDist(self):

        if self.exon != "NA":

            if self.gene.multi == False or self.exon.number == 0:
                
                self.distExon = abs(self.exon.start - self.ribostart)
            
            else:
                    
                dist = 0
                for i in range(self.exon.number):

                    dist += self.gene.exons_list[i].length

                dist += abs(self.exon.start - self.ribostart) + self.exon.abs_frame 
                self.distExon = dist
                    


class GeneStructureError(Exception):
    def __init__(self, message):
        self.message = message
        super().__init__(self.message)
