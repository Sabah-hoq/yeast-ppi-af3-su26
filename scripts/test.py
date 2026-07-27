import os
import argparse
from pathlib import Path
import polars as pl
from string__downloader import download_string_data  
from load_data import load_data

def load_and_map_data(data_dir_path):
    data_dir = Path(data_dir_path)
    
    df2 = download_string_data(data_id="protein.physical.links.detailed", cols_to_clean=["protein1", "protein2"])
    string_aliases = download_string_data(data_id="protein.aliases")
    string_info = download_string_data(data_id="protein.info")

    pairs, confidences = load_data("data/")

    unique_protiens = (
    pl.concat([
        pairs.select(pl.col("af3_id1").alias("protein_id")),
        pairs.select(pl.col("af3_id2").alias("protein_id"))
    ])
    .drop_nulls()
    .unique()
    .sort("protein_id")
    )

    map = (
        pl.concat([
            pairs.select([
                pl.col("af3_id1").alias("uniprot_id"),
                pl.col("michaelis2023:Source Gene names (SGD/UniProt-primary or ordered locus)") #some of these are empty 
                    .alias("gene_name"),
        ]),
        pairs.select([
        pl.col("af3_id2").alias("uniprot_id"),
        pl.col("michaelis2023:Target Gene names  (SGD/UniProt-primary or ordered locus)") #some of these are empty
            .alias("gene_name"),
            ])
        ])
        .drop_nulls()
        .unique(subset =["uniprot_id"])
    )

    coverage = map.collect().height
    total = unique_protiens.collect().height
    print(f"Already mapped: {coverage}/{total} ({100 *coverage/total:.1f}%)")
    #What we are able to map ^

    already_mapped_ids = map.collect()["uniprot_id"].to_list()
    need_api = unique_protiens.filter(
        pl.col("protein_id").is_in(already_mapped_ids).not_()
    )
    print(f"Need API lookup: {need_api.collect().height}")

    string_aliases = (
    pl.scan_csv(str(data_dir / "4932.protein.aliases.v12.0.txt"), separator="\t")
    .rename({"#string_protein_id": "string_id", "alias": "protein_id"})
    .with_columns(pl.col("protein_id").str.to_lowercase())
    .filter(pl.col("source").str.contains("UniProt"))
    )

    matched_missing = (
        need_api.join(
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

    unique_protiens = (
        unique_protiens.join(
            string_aliases,
            on="protein_id",
            how="inner"
        )
        .select(
            pl.col("protein_id").alias("uniprot_id"),
            pl.col("string_id")
        )
    )

    string_info = (
        pl.scan_csv(str(data_dir / "4932.protein.info.v12.0.txt"), separator="\t")
        .rename({"#string_protein_id": "string_id"})
    )

    full_matches = (
        matched_missing.join(
            string_info,
            on="string_id",
            how="inner"
        )
        .select([
            pl.col("uniprot_id"),
            pl.col("string_id"),
            pl.col("preferred_name").alias("string_gene_name"),
            pl.col("annotation")
        ])
    )

    unique_protiens = (
        unique_protiens.join(
            string_info,
            on="string_id",
            how="inner"
        )
        .select([
            pl.col("uniprot_id"),
            pl.col("string_id"),
            pl.col("preferred_name").alias("string_gene_name"),
            pl.col("annotation")
        ])
    )

    # see your results
    df_matches = full_matches.collect()

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

    full_matches = (
        matched_missing.join(
            string_info2,
            on="string_id",
            how="inner"
        )
        .select([
            pl.col("uniprot_id"),
            pl.col("string_id"),
            pl.col("preferred_name").alias("string_gene_name"),
            pl.col("annotation")
        ])
    )

    df_matches = full_matches.collect()

    df_matches_unique_lazy = df_matches.unique(subset=['uniprot_id'], keep='first').lazy()

    df_final_lazy = unique_protiens.unique(subset=['uniprot_id'], keep='first').join(
        other=df_matches_unique_lazy,
        on='uniprot_id',
        how='left'
        ).select(pl.exclude("^.*_right$"))

    df_final = df_final_lazy.collect()

    df2_2 = (
    df2.select([
        pl.col("protein1"),
        pl.col("protein2"),
        pl.col("combined_score")
    ])
    .collect()
    .unique()
    )

    # Now to do this with pairs 
    pairs_2 = (
        pairs.select([
            pl.col("af3_id1").alias("protein1"),
            pl.col("af3_id2").alias("protein2"),
            pl.col("chain_pair_iptm_best"),
            pl.col("chain_pair_iptm_mean"),
            pl.col("chain_pair_iptm_best_corrected"),
            pl.col("chain_pair_iptm_mean_corrected")
        ])
        .collect()
        .drop_nulls()
        .unique()
    )

    id_map_dict = dict(zip(df_final['uniprot_id'], df_final['string_id']))

    # using pairs_2 protein1 and protein2 to make a new map
    df_alphafold_mapped = pairs_2.with_columns([
        pl.col("protein1").replace_strict(id_map_dict, default=None).alias("string_id1"),
        pl.col("protein2").replace_strict(id_map_dict, default=None).alias("string_id2")
    ])
    # Drop any rows where the mapping returned None 
    df_alphafold_mapped = df_alphafold_mapped.drop_nulls(subset=["string_id1", "string_id2"])

    df_all_mapped = df_alphafold_mapped.with_columns(
    pl.min_horizontal("string_id1", "string_id2").alias("pair_key1"),
    pl.max_horizontal("string_id1", "string_id2").alias("pair_key2")
    )

    df_string_ordered = df2_2.with_columns(
        pl.min_horizontal("protein1", "protein2").alias("pair_key1"),
        pl.max_horizontal("protein1", "protein2").alias("pair_key2")
    )

    df_string_unique_pairs = df_string_ordered.select(
        ["pair_key1", "pair_key2", "combined_score"]
    ).unique(subset=["pair_key1", "pair_key2"])

    df_final_comparison = df_all_mapped.join(
        df_string_unique_pairs,
        on=["pair_key1", "pair_key2"],
        how="left"
    ).drop("pair_key1", "pair_key2")

    print(df_final_comparison["combined_score"].count()) # Should be 102658...not 0

    return df_final_comparison


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Global Score Pipeline Wrapper")
    parser.add_argument("--data-dir", type=str, required=True, help="Path to local data folder (with summary_pairs.parquet)")
    parser.add_argument("--output", type=str, default="final_comp.csv", help="output filename")
    args = parser.parse_args()
    
    print("Running mapping pipeline standalone...")

    final_df = load_and_map_data(args.data_dir)
    print(final_df)