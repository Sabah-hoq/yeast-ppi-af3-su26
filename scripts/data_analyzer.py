from load_data import load_data
import polars as pl

class DataAnalyzer:
    def __init__(self, data_dir="data/"):
        self.data_dir = data_dir

    def clean_data(self):
        pairs, confidences = load_data(self.data_dir)

        unique_proteins = pl.concat(
            [
                pairs.select(pl.col("af3_id1").alias("protein_id")).unique(),
                pairs.select(pl.col("af3_id2").alias("protein_id")).unique(),
            ]
        ).unique()
        print(
            f"Total Unique Proteins: {(unique_proteins.select(pl.len()).collect(engine='streaming').item()):,}"
        )

        lf_merged = pairs.join(confidences, on="name", how="inner")
        lf_confident = lf_merged.filter(pl.col("chain_pair_iptm_best") >= 0.6)

        print("Removing duplicate structural runs...")
        lf_clean = lf_confident.unique(subset=["af3_id1", "af3_id2"])

        lf_final = lf_clean.select(
            [
                "af3_id1",
                "af3_id2",
                "chain_pair_iptm_best",
                "chain_pair_iptm_mean",
                "ptm",
                "iptm",
            ]
        )

        print("Processing optimized data stream...")
        df_result = lf_final.collect(streaming=True)

        print(f"Total Unique Confident Interactions Found: {len(df_result):,}")
        print(df_result.head())

        unique_p1 = set(df_result["af3_id1"])
        unique_p2 = set(df_result["af3_id2"])
        total_unique_proteins = len(unique_p1.union(unique_p2))

        print(f"\nTotal Unique Proteins in the Confident Network: {total_unique_proteins:,}")
        return df_result, total_unique_proteins

    def stats(self, df_result):
        df_proteins_long = df_result.melt(
            id_vars=["chain_pair_iptm_best"],
            value_vars=["af3_id1", "af3_id2"],
            value_name="protein_id",
        )

        protein_stats = (
            df_proteins_long.group_by("protein_id")
            .agg(
                pl.len().alias("confident_interaction_count"),
                pl.col("chain_pair_iptm_best").mean().alias("protein_avg_iptm"),
            )
        )
        mean_of_protein_means = protein_stats["protein_avg_iptm"].mean()
        std_of_protein_means = protein_stats["protein_avg_iptm"].std()

        return mean_of_protein_means, std_of_protein_means, protein_stats

    def analysis(self, mean_of_protein_means, std_of_protein_means, protein_stats):
        bias_analysis_df = (
            protein_stats.with_columns(
                (pl.col("protein_avg_iptm") - mean_of_protein_means).alias("iptm_diff_from_avg"),
                ((pl.col("protein_avg_iptm") - mean_of_protein_means) / std_of_protein_means).alias(
                    "standard_deviation"
                ),
            )
            .sort("protein_avg_iptm", descending=True)
        )
        return bias_analysis_df

    def run_pipeline(self):
        df_result, total_unique_proteins = self.clean_data()
        mean_of_protein_means, std_of_protein_means, protein_stats = self.stats(df_result)
        bias_analysis_df = self.analysis(
            mean_of_protein_means,
            std_of_protein_means,
            protein_stats,
        )
        return bias_analysis_df, protein_stats, total_unique_proteins