import polars as pl
import gff3_parser

orfs = pl.from_pandas(gff3_parser.parse_gff3("input/mapping_orf_Scer_SGD_noMT.gff", parse_attributes=True))

GFF = gff3_parser.parse_gff3("input/Scer.gff", parse_attributes=True)

COLUMNS = [["seq_id", "start", "end", "strand", "phase", "attributes"]]

GFF_POLARS = pl.from_pandas(GFF)

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
