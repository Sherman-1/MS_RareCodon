from Bio.SeqIO import parse
from Bio.Seq import Seq

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
    