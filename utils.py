import polars as pl
import re
from Bio.Seq import Seq
from Bio.SeqIO import parse
from data import ORF_DF_COLUMNS, GFF_POLARS
from classes import Gene, Exon


def multifasta_to_dict(path, genome = False):

    """
    Reads a FASTA file and returns a dictionary where each record in the file
    is a key-value pair with the record identifier as the key and the record sequence (in uppercase) as the value.

    Args:
    path (str): The path to the FASTA file to be read.

    Returns:
    dict: A dictionary where the keys are the record identifiers and the values
          are the corresponding record sequences as Seq objects in uppercase letters.
    """
    
    records = parse(path, "fasta")

    if genome:
        dico = {}
        for record in records:
            dico[record.id] = {}
            dico[record.id]["seq"] = Seq(str(record.seq).upper())
            dico[record.id]["len"] = len(record.seq)
            
        return dico
    
    else:

        return {record.id: Seq(str(record.seq).upper()) for record in records}
    
    
def get_good_column_names(columns):


    parent_col = [col for col in columns if re.match(r"[Pp]arent", col)]
    name_col = [col for col in columns if re.match(r"[Nn]ame", col)]

    if len(parent_col) == 1 and len(name_col) == 1:

        return str(parent_col[0]), str(name_col[0])
    
    else:

        raise AttributeError("Problem with GFF columns : Parent or Name not found. See get_good_column_names()")
        

def return_gene_infos(gene_infos) -> dict:

    chromosome = gene_infos["Seqid"].unique().to_list()[0]
    multi = gene_infos["Type"].to_list().count("CDS") > 1

    gene_infos = gene_infos.filter(pl.col("Type") == "gene")

    start = gene_infos["Start"].to_list()[0]
    end = gene_infos["End"].to_list()[0]
    sense = gene_infos["Strand"].to_list()[0]
    
    return {

        "chromosome" : chromosome,
        "start" : int(start),
        "end" : int(end),
        "sense" : sense,
        "multi" : multi
    } 

def init_gene_object(gene_id, gff_dataframe):


    """
    
    On veut pour la gene_id donnée initiliaser un objet Gene avec ses attributs
    tirés du gff_dataframe.
    
    """

    pattern = fr".*{gene_id}.*"

    parent, name = get_good_column_names(gff_dataframe.columns)

    # Polars does not support regex filtering : pandas is used instead
    gene_rows = pl.from_pandas(gff_dataframe[
        gff_dataframe[parent].str.contains(pattern, regex=True, na=False) |
        gff_dataframe[name].str.contains(pattern, regex=True, na=False)
    ]) # Get all rows related to the gene being initialized


    
    gene_infos = return_gene_infos(gene_rows)


    gene = Gene(

        ID = gene_id,
        chromosome = gene_infos["chromosome"],
        start = gene_infos["start"],
        end = gene_infos["end"],
        sense = gene_infos["sense"],
        multi = gene_infos["multi"],
        frame = gene_infos["start"] % 3
        
    )
    

    cds_counter = 1
    for exon in gene_rows.filter(pl.col("Type") == "CDS").iter_rows(named = True):

        key = f'{exon["ID"]}-{cds_counter}' if exon["ID"] else f'{exon["Name"]}-{cds_counter}'
        gene.add_exon(

            key = key,
            value = Exon(
                        ID = key,
                        gene = gene,
                        start = int(exon["Start"]),
                        end = int(exon["End"]),
                        abs_frame = int(exon["Start"]) % 3)
        )

        cds_counter = cds_counter + 1

    gene.sort() 

    return gene

def check_double_overlap(row : tuple):

    orf = dict(zip(ORF_DF_COLUMNS, row))
    
    overlaps = [match for item in orf["Ovp_with"].split("|") for match in re.findall(r"\b([\w-]+)_mRNA\b", item)]

    if len(overlaps) == 0:

        overlaps = [match for item in orf["Ovp_with"].split("|") for match in re.findall(r"\b([\w-]+)_CDS\b", item)]
    
    buffer = []

    if len(overlaps) != 1:

        for overlap in overlaps:

            
            if GFF_POLARS.filter(
                
                (pl.col("ID") == overlap)
                )["Strand"].unique().to_list()[0] == orf["Strand"]:

                buffer.append(overlap)
                

        if len(buffer) == 1: # Several genes are overlapped by the ORF, but only one is on the same strand

            orf["Ovp_gene"] = buffer[0]
            return tuple(orf.values())
        
    
        elif len(buffer) == 0: # No gene found on the same strand as the ORF

            orf["Ovp_gene"] = "NA"
            return tuple(orf.values())
        
        else: # Several genes found on the same strand as the ORF

            orf["Ovp_gene"] = "Two_or_more_genes"
            return tuple(orf.values())


    else: # If there is only one gene found in the overlapping information given by ORFMine ID

        orf["Ovp_gene"] = overlaps[0]
        return tuple(orf.values())
