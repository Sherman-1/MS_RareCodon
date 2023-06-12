import polars as pl
import gff3_parser
from files import multifasta_to_dict
from halo import Halo

spinner = Halo(text='Loading data . . .', spinner='dots')
spinner.start()


ORFS = pl.from_pandas(gff3_parser.parse_gff3("input/mapping_orf_Scer_SGD_noMT.gff", parse_attributes=True, verbose = False))

GFF = gff3_parser.parse_gff3("input/Scer.gff", parse_attributes=True, verbose = False)

GFF_POLARS = pl.from_pandas(GFF)

RIBO = pl.read_csv("input/aORFs_complete_info.csv", infer_schema_length = 10000).filter(pl.col("first.codon.idx") != "NA")

PATTERN = r'\b([\w-]+)_CDS\b'

FASTA_DICT = multifasta_to_dict("input/Scer.fna", genome = True)

NB_NT = 10

MAX_CODONS = 3

spinner.stop()

ORF_DF_COLUMNS = ['Seqid',
 'Source',
 'Type',
 'Start',
 'End',
 'Score',
 'Strand',
 'Phase',
 'Status',
 'color',
 'Parent',
 'ID',
 'Ovp_with',
 'Ovp_gene']

RIBO_DF_COLUMNS = ['Seq_ID',
 'Num_reads',
 'Num_p0',
 'Num_p1',
 'Num_p2',
 'Perc_p0',
 'Perc_p1',
 'Perc_p2',
 'HCA',
 'Disord',
 'Aggreg',
 'Relative.species',
 'Relative.age',
 'prot.size',
 'GGC',
 'GCC',
 'GTC',
 'GAC',
 'GGG',
 'CCC',
 'GGA',
 'TCC',
 'GAG',
 'CTC',
 'GGT',
 'ACC',
 'GTT',
 'AAC',
 'GCG',
 'CGC',
 'CCG',
 'CGG',
 'TCG',
 'CGA',
 'ACG',
 'CGT',
 'GAT',
 'ATC',
 'GCT',
 'AGC',
 'ACT',
 'AGT',
 'CCT',
 'AGG',
 'TCT',
 'AGA',
 'CTT',
 'AAG',
 'CTG',
 'CAG',
 'CAA',
 'TTG',
 'GCA',
 'TGC',
 'ACA',
 'TGT',
 'GAA',
 'TTC',
 'AAA',
 'TTT',
 'GTG',
 'CAC',
 'CAT',
 'ATG',
 'GTA',
 'TAC',
 'CTA',
 'TAG',
 'ATA',
 'TAT',
 'AAT',
 'ATT',
 'TTA',
 'TAA',
 'CCA',
 'TGG',
 'TCA',
 'TGA',
 'msms.ynb.sds',
 'msms.ynb.sds.mg132',
 'coverage',
 'TPM',
 'MSMS',
 'RPKM',
 'firstATG',
 'first_near_cognate',
 'first.codon',
 'first.codon.idx',
 'chromosome_ID']
