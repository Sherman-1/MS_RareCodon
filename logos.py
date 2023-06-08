import polars as pl 
import logomaker as lm 



data = pl.read_csv("output/output.csv", infer_schema_length = 10000)


data = data.filter(

    pl.col("ribostartLocalisation") == "exon"

)


downstream = data["downstream"].to_list()
upstream = data["upstream"].to_list()

down_counts = lm.alignment_to_matrix(downstream)
up_counts = lm.alignment_to_matrix(upstream)

down_logo = lm.Logo(down_counts)
up_logo = lm.Logo(up_counts)

down_logo.ax.set_title("Exon downstream")
up_logo.ax.set_title("Exon upstream")

# Draw logos into output folder

down_logo.fig.savefig("output/downstream.png")
up_logo.fig.savefig("output/upstream.png")
