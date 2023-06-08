import polars as pl
import gff3_parser
from files import multifasta_to_dict

ORFS = pl.from_pandas(gff3_parser.parse_gff3("input/mapping_orf_Scer_SGD_noMT.gff", parse_attributes=True))

GFF = gff3_parser.parse_gff3("input/Scer.gff", parse_attributes=True)

GFF_POLARS = pl.from_pandas(GFF)

RIBO = pl.read_csv("input/aORFs_complete_info.csv").filter(pl.col("first.codon.idx") != "NA")

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

PATTERN = r'\b([\w-]+)_CDS\b'

FASTA_DICT = multifasta_to_dict("input/Scer.fna", genome = True)