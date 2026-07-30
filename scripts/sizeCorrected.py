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
ready_file_df = ready_file.collect()

# Save file
repo_root = Path(__file__).resolve().parents[1]
output_dir = repo_root / "outputs" / "spras"
output_dir.mkdir(parents=True, exist_ok=True)
output_path = output_dir / "yeast_spras.txt"
ready_file_df.write_csv(output_path, include_header=False)
print(f"CSV saved successfully at: {output_path}")