# May crash when running...dont run
from pathlib import Path
import polars as pl
import matplotlib.pyplot as plt
from scipy import stats
import numpy as np
import seaborn as sns
from load_data import load_data
from sklearn.metrics import  roc_curve, auc, roc_auc_score
from sklearn.metrics import precision_recall_curve, auc, average_precision_score

pairs, confidences = load_data("data/")
sizeCoreceted = (pairs.collect())["chain_pair_iptm_mean_corrected"]

# fig, ax = plt.subplots()
# plt.hist(sizeCoreceted, bins=100000, color='teal', edgecolor='black', alpha=0.5)
# ax.ticklabel_format(style='plain', axis='y')
# ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, loc: "{:,}".format(int(x))))
# plt.xlim(-0.015, 0.8)
# # Display the plot
# plt.show()

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

# More figures below for visualization 
sns.set_theme(style="whitegrid")
fig,ax = plt.subplots(1,2, figsize=(14,5))
sns.histplot(
    sizeCoreceted,
    bins = 1000,
    kde=True,
    color="#5157a1",
    edgecolor = None, 
    stat="density",)
ax.set_title(
    "iptm scoresize corrected (unfiltered)"
)
ax.legend()