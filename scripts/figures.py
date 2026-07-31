from data_analyzer import DataAnalyzer as da
from load_data import load_data
import seaborn as sns
import matplotlib.pyplot as plt
from pathlib import Path
import os

pairs, confidences = load_data("data/")

# distributionYeastScoreValue.png figure code 
df_metrics = pairs.unique(subset=["af3_id1", "af3_id2"]).select(["chain_pair_iptm_best"]).collect()
data = df_metrics["chain_pair_iptm_best"].to_numpy()
fig, ax = plt.subplots(figsize=(8, 5))
plt.hist(data, bins=30, edgecolor="#3B4A62", color="#8399C8")

ax.hist(data, bins=30, edgecolor="#3B4A62", color="#8399C8")

ax.ticklabel_format(style='plain', axis='y')
ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, loc: f"{int(x):,}"))

ax.set_xlabel('AlphaFold ipTM Score')
ax.set_ylabel('Number of Pairs (Frequency)')
ax.set_title('Distribution of Yeast Interactome Score Values')

output_dir = '../output/figures'
output_path = os.path.join(output_dir, 'distributionYeastScoreValue.png') #already in folder
plt.savefig(output_path, dpi=600, bbox_inches='tight')

plt.show()


# interactionCount.png figure 


# confidentInteraction.png figure
da = da()
bias_analysis_df, _, _ = da.run_pipeline() 
sns.boxplot(
    bias_analysis_df,
    y="protein_avg_iptm", 
    x="confident_interaction_count", 
    palette= "Set2", 
    boxprops=dict(alpha = 0.8), 
    fliersize=2)

output_dir = '../output/figures'
output_path = os.path.join(output_dir, 'confidentInteraction.png') #already in folder
plt.savefig(output_path, dpi=600, bbox_inches='tight')

plt.show()

# iptmBest_stringScore.png figure 
sns.set_theme(style="whitegrid")

g2 = sns.jointplot(
    data=df_final_comparison_sig2, 
    x='combined_score', 
    y='chain_pair_iptm_best', 
    kind='reg',
    scatter=False)

scatter = g2.ax_joint.scatter(
    data=df_final_comparison_sig2,
    x='combined_score',
    y='chain_pair_iptm_best',
    c=df_final_comparison_sig2['combined_score'], 
    cmap='viridis',
    edgecolor='w',                            
    linewidth=0.5
)

plt.show()

# iptmMean_stringScore.png figure 
sns.set_theme(style="whitegrid")

g = sns.jointplot(
    data=df_final_comparison_sig, 
    x='combined_score', 
    y='chain_pair_iptm_mean', 
    kind='reg',
    scatter=False)

scatter = g.ax_joint.scatter(
    data=df_final_comparison_sig,
    x='combined_score',
    y='chain_pair_iptm_mean',
    c=df_final_comparison_sig['combined_score'], 
    cmap='viridis',
    edgecolor='w',                            
    linewidth=0.5
)

plt.show()