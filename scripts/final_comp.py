import os
import argparse
from pathlib import Path
import polars as pl
from string__downloader import download_string_data  
from load_data import load_data
from data_analyzer import DataAnalyzer

def normalize_string_id_column(lf: pl.LazyFrame) -> pl.LazyFrame:
    string_col = "string_protein_id" if "string_protein_id" in lf.columns else "#string_protein_id"
    return (
        lf.rename({string_col: "string_id"})
        .with_columns(
            pl.col("string_id").str.replace(r"^4932\.", "", literal=True).alias("string_id")
        )
    )


def load_and_map_data(data_dir_path):
    data_dir = Path(data_dir_path)
    
    df2 = download_string_data(data_id="protein.physical.links.detailed", organism_id=4932,cols_to_clean=["protein1", "protein2"])
    string_aliases = download_string_data(data_id="protein.aliases")
    string_info = download_string_data(data_id="protein.info")

    pairs, confideences = load_data(str(data_dir))
    unique_proteins =(
        pl.concat([
            pairs.select(pl.col("af3_id1").alias("protein_id")),
            pairs.select(pl.col("af3_id2").alias("protein_id"))
        ]).drop_nulls().unique().sort("protein_id"))

    map = pl.concat([
        pairs.select([
            pl.col("af3_id1").alias("uniprot_id"),
            pl.col("michaelis2023:Source Gene names (SGD/UniProt-primary or ordered locus)").alias("gene_name"),
        ]),
        pairs.select([
            pl.col("af3_id2").alias("uniprot_id"),
            pl.col("michaelis2023:Target Gene names  (SGD/UniProt-primary or ordered locus)").alias("gene_name"),
        ])
    ]).drop_nulls().unique(subset=["uniprot_id"])

    already_mapped_ids = map.select("uniprot_id").collect().get_column("uniprot_id").to_list()
    need_api = unique_proteins.filter(pl.col("protein_id").is_in(already_mapped_ids).not_())
    
    clean_aliases = (
        string_aliases
        .rename({"alias": "protein_id"})
        .pipe(normalize_string_id_column)
        .with_columns(pl.col("protein_id").str.to_lowercase())
        .filter(pl.col("source").str.contains("UniProt"))
        )
    
    clean_info = normalize_string_id_column(string_info)
    
    matched_missing = (
        need_api.join(clean_aliases, on="protein_id", how="inner")
        .select([pl.col("protein_id").alias("uniprot_id"), pl.col("string_id")])
        .unique()
    )

    df_matches = (
        matched_missing.join(clean_info, on="string_id", how="inner")
        .select([
            pl.col("uniprot_id"),
            pl.col("string_id"),
            pl.col("preferred_name").alias("string_gene_name"),
            pl.col("annotation"),
        ])
        .collect()
    )

    unique_proteins_mapped = (
        unique_proteins.join(clean_aliases, on="protein_id", how="inner")
        .select([pl.col("protein_id").alias("uniprot_id"), pl.col("string_id")])
        .join(clean_info, on="string_id", how="inner")
        .select([
            pl.col("uniprot_id"),
            pl.col("string_id"),
            pl.col("preferred_name").alias("string_gene_name"),
            pl.col("annotation"),
        ])
        .unique(subset=["uniprot_id"])
        .collect()
    )

    df_matches_unique = df_matches.unique(subset=["uniprot_id"], keep="first")
    df_final_mapping = (
        unique_proteins_mapped.unique(subset=["uniprot_id"], keep="first")
        .join(df_matches_unique, on="uniprot_id", how="left")
        .select(pl.exclude(r"^.*_right$"))
    )
    id_map_dict = dict(zip(df_final_mapping["uniprot_id"], df_final_mapping["string_id"]))
    
    pairs_2 = pairs.select([
            pl.col("af3_id1").alias("protein1"), pl.col("af3_id2").alias("protein2"),
            pl.col("chain_pair_iptm_best"), pl.col("chain_pair_iptm_mean")
        ]).drop_nulls().unique()

    
    
        # using pairs_2 protein1 and protein2 to make a new map
    df_alphafold_mapped = pairs_2.with_columns([
            pl.col("protein1").replace(id_map_dict, default=None).alias("string_id1"),
            pl.col("protein2").replace(id_map_dict, default=None).alias("string_id2")
        ])
    
        # FIX 3: drop rows where the mapping didn't resolve, so they don't
        # inflate the final row count with unmatchable garbage rows.
    df_alphafold_mapped = df_alphafold_mapped.drop_nulls(subset=["string_id1", "string_id2"])
    
    df_all_mapped = (
        df_alphafold_mapped.with_columns(
            pl.col("string_id1").str.replace(r"^4932\.", "", literal=True).alias("string_id1_norm"),
            pl.col("string_id2").str.replace(r"^4932\.", "", literal=True).alias("string_id2_norm"),
        )
        .with_columns(
            pl.min_horizontal("string_id1_norm", "string_id2_norm").alias("pair_key1"),
            pl.max_horizontal("string_id1_norm", "string_id2_norm").alias("pair_key2"),
        )
        .drop(["string_id1_norm", "string_id2_norm"])
    )

    matched = (
    map.join(
        string_aliases,
        on="protein_id",
        how="inner"
    )
    .select([
        pl.col("protein_id").alias("uniprot_id"),
        pl.col("string_id")
    ])
    .unique()
    )
    string_info2 = (
        pl.scan_csv(
            str(data_dir / "4932.protein.info.v12.0.txt"),
            separator="\t"
        )
        .rename({"#string_protein_id": "string_id"})
    ) 
    df2_2 = df2.select(["protein1", "protein2", "combined_score"]).unique()
    print(df2_2.collect().head())

    df_string_ordered = (
        df2_2.with_columns(
            pl.col("protein1").str.replace(r"^4932\.", "", literal=True).alias("protein1_norm"),
            pl.col("protein2").str.replace(r"^4932\.", "", literal=True).alias("protein2_norm"),
        )
        .with_columns(
            pl.min_horizontal("protein1_norm", "protein2_norm").alias("pair_key1"),
            pl.max_horizontal("protein1_norm", "protein2_norm").alias("pair_key2"),
        )
        .drop(["protein1_norm", "protein2_norm"])
    )
    df_string_unique_pairs = df_string_ordered.select(["pair_key1", "pair_key2", "combined_score"]).unique(subset=["pair_key1", "pair_key2"])
    print(df_string_unique_pairs.collect().head())
    print(df_all_mapped.collect().head())
    df_final_comparison = df_all_mapped.join(
            df_string_unique_pairs, on=["pair_key1", "pair_key2"], how="left"
        ).drop("pair_key1", "pair_key2")
    
    print(df_final_comparison.collect()["combined_score"])

    return df_final_comparison.collect()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Global Score Pipeline Wrapper")
    parser.add_argument("--data-dir", type=str, required=True, help="Path to local data folder (with summary_pairs.parquet)")
    parser.add_argument("--output", type=str, default="final_comp.csv", help="output filename")
    args = parser.parse_args()

    final_df = load_and_map_data(args.data_dir)
    
    data_folder = Path(args.data_dir)
    csv_path = data_folder / args.output
    final_df.write_csv(csv_path)
    print(f"CSV file saved at: {csv_path} (Total rows: {final_df.height})") #this is should be 7.5mil 

    repo_root = data_folder.parent
    folder_path = repo_root / "notebooks/PR_ROC"
    folder_path.mkdir(parents=True, exist_ok=True)
    csv_file_path = folder_path / "final_comp.csv"

    final_df.write_csv(csv_file_path)
    print(f"CSV saved successfully at: {csv_file_path}")

