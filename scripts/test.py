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

    pairs, confidences = load_data(str(data_dir))



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Global Score Pipeline Wrapper")
    parser.add_argument("--data-dir", type=str, required=True, help="Path to local data folder (with summary_pairs.parquet)")
    parser.add_argument("--output", type=str, default="final_comp.csv", help="output filename")
    args = parser.parse_args()
    
    print("Running mapping pipeline standalone...")

    final_df = load_and_map_data(args.data_dir)
    print(final_df)