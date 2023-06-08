from classes import Orf
import csv


def write_data(gene_list, filename = "output/output.csv"):

    attribute_names = ["ID","ribospike", "MSMS", "MS_sds", "MS_cond", "age_rel", 
             "upstream", "downstream", "rel_frame", "ribostartLocalisation",
             "start", "end", "ribostart","start_seq"] 
    
    colnames = attribute_names.copy()
    colnames.extend(["gene","sense","chromosome"])
    
    
    with open(filename, mode="w", newline="") as file:
        writer = csv.writer(file)
        writer.writerow(colnames)

        for gene in gene_list:

            for orf in gene.orfs_list:

                row = [orf.__dict__.get(attr_name, "") for attr_name in attribute_names]
                row.extend([gene.ID,gene.sense,gene.chromosome])

                writer.writerow(row)

            