# Imports 
from pathlib import Path
import polars as pl
import matplotlib.pyplot as plt
from load_data import load_data

# Use final_comp.csv for Data 
pairs, confidences = load_data("data/")
sizeCoreceted = (pairs.collect())["chain_pair_iptm_mean_corrected"]

# filtering data
df_filtered = pairs.filter(pl.col("chain_pair_iptm_mean_corrected") >= 0.061) 

## Min and Max Values 
max = (df_filtered.collect())["chain_pair_iptm_mean_corrected"].max()
min = (df_filtered.collect())["chain_pair_iptm_mean_corrected"].min()
x = (df_filtered.collect())["chain_pair_iptm_mean_corrected"]

# Spras ready format 
normalization = (x - min)/ (max - min)
top_edges3_normal = df_filtered.with_columns(
    (pl.col("chain_pair_iptm_mean_corrected")/max)
    )
# Sanity Check
print(len(top_edges3_normal.collect()))
print((top_edges3_normal.collect())["chain_pair_iptm_mean_corrected"].describe())

# plot 
fig, ax = plt.subplots()
plt.hist((top_edges3_normal.collect())["chain_pair_iptm_mean_corrected"], bins=50000, color='teal', edgecolor='black', alpha=0.5)
ax.ticklabel_format(style='plain', axis='y')
ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, loc: "{:,}".format(int(x))))
plt.show()

# Saving file in spras format using final_comp,csv
ready_file = top_edges3_normal.select([
    pl.col("chain_pair_iptm_mean_corrected").alias("weight"),
    pl.col("af3_id1").alias("protein1"),
    pl.col("af3_id2").alias("protein2"),
])

ready_file = ready_file.with_columns(
        pl.col("protein1").str.to_uppercase(),
        pl.col("protein2").str.to_uppercase()
    )

ready_file = ready_file.with_columns(direction=pl.lit("U"))
ready_file = ready_file.select(["protein1","protein2","weight","direction"])

print(ready_file.collect().head())

Path("../outputs/spras").mkdir(parents=True, exist_ok=True)
ready_file.collect().write_csv("../outputs/spras/yeast_spras.txt", include_header=False)