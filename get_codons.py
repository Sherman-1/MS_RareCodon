
import polars as pl   
import csv



def get_upstream(data_path = "output/output.csv", output_path = "output/upstream_codons.csv", nb_codons : int = 3):

    data = pl.read_csv(data_path, infer_schema_length = 10000).filter(pl.col("ribostartLocalisation") == "exon")

    threshold_f1 = ( nb_codons * 3 ) + 1
    threshold_f2 = ( nb_codons * 3 ) + 2

    handle = open(output_path, 'w')
    writer = csv.writer(handle)
    writer.writerow(['-3','-2','-1','ID'])

    sorted = data.filter(
        
        ((pl.col('upstream').apply(lambda x: len(x)) >= threshold_f1) & (pl.col('rel_frame') == 1)) |
        ((pl.col('upstream').apply(lambda x: len(x)) >= threshold_f2) & (pl.col('rel_frame') == 2))

    )

    sorted = sorted.with_columns((pl.col('upstream') + pl.col('start_seq')).alias('sequence'))




    # Case F1 validated
    for orf in sorted.filter(pl.col('rel_frame') == 1).iter_rows(named = True):
        codons = []
        for i in range(nb_codons):

            ribostart = len(orf['sequence']) - 3
            start = ribostart - 1 - 3 * i
            end = ribostart + 2 - 3 * i
            codons.append(orf['sequence'][start:end])
        
        codons.reverse()
        codons.append(orf['ID'])
        writer.writerow(codons)
        

    # Case F2 validated
    for orf in sorted.filter(pl.col('rel_frame') == 2).iter_rows(named = True):
        codons = []
        for i in range(nb_codons):

            ribostart = len(orf['sequence']) - 3
            start = ribostart - 2 - 3 * i
            end = ribostart + 1 - 3 * i
            codons.append(orf['sequence'][start:end])

        codons.reverse()
        codons.append(orf['ID'])
        writer.writerow(codons)


    handle.close()

    return 0

        
def get_downstream(data_path = "output/output.csv",  output_path = "output/downstream_codons.csv", nb_codons : int = 3):

    data = pl.read_csv(data_path, infer_schema_length = 10000).filter(pl.col("ribostartLocalisation") == "exon")

    threshold_f1 = ( nb_codons * 3 )  + 2
    threshold_f2 = ( nb_codons * 3 ) + 1

    handle = open(output_path, 'w')
    writer = csv.writer(handle)
    writer.writerow(['+1','+2','+3','ID'])
    

    sorted = data.filter(
        
        ((pl.col('downstream').apply(lambda x: len(x)) >= threshold_f1) & (pl.col('rel_frame') == 1)) |
        ((pl.col('downstream').apply(lambda x: len(x)) >= threshold_f2) & (pl.col('rel_frame') == 2))

    )

    sorted = sorted.with_columns((pl.col('start_seq') + pl.col('downstream')).alias('sequence'))


    # Case F1 validated
    for orf in sorted.filter(pl.col('rel_frame') == 1).iter_rows(named = True):

        codons = []
        for i in range(nb_codons):

            ribostart = 0
            start = ribostart + 2 + 3 * i
            end = ribostart + 5 + 3 * i
            codons.append(orf['sequence'][start:end])

        codons.append(orf['ID'])
        writer.writerow(codons)


    # Case F2 validated
    for orf in sorted.filter(pl.col('rel_frame') == 2).iter_rows(named = True):

        codons = []
        for i in range(nb_codons):

            ribostart = 0
            start = ribostart + 1 + 3 * i
            end = ribostart + 4 + 3 * i
            codons.append(orf['sequence'][start:end])
            
        codons.append(orf['ID'])
        writer.writerow(codons)

    return 0    

if __name__ == "__main__":
    
    print("test 1")
    get_downstream()
    get_upstream()    
    print("test 2")