### Cytoscape
---
All files here (after running relevant scripts) will aid in Cytoscape visualization. 

`yeast_network.cys`: Contains only high-scoring PPIs (iptm $\geq$ 0.6 __AND__ STRING $\geq$ 400)

Files used for cytoscape: 

**df_final_comparison**:

* `df_final_comparison_mean.csv`: Looks at the "chain_pair_iptm_mean" score where it's threshold to be $\geq$ 0.6
 
* `df_final_comparison_best.csv`: Looks at the "chain_pair_iptm_best" score where it's threshold to be $\geq$ 0.6

**df_final_comparison_string**: 

* `df_final_comparison_mean.csv`: Looks at the "combined_score" score where it's threshold to be $\geq$ 400, and uses "chain_pair_iptm_mean" as a base line to decide how many edges there should be. Such that this file looks at where both the STRING scores the iptm_mean are considered good. 

* `df_final_comparison_best.csv`: Looks at the "combined_score" score where it's threshold to be $\geq$ 400, and uses "chain_pair_iptm_best" as a base line to decide how many edges there should be. Such that this file looks at where both the STRING scores the iptm_best are considered good. 