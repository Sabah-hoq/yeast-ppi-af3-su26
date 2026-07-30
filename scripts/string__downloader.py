## Copying from here: https://github.com/jurgjn/pooled-ppi/blob/main/src/pooled_ppi/string_db.py
#!/usr/bin/env python3
from pathlib import Path
import polars as pl
from cached_path import cached_path

def download_string_data(
    data_id: str = "protein.physical.links.detailed", 
    organism_id: int = 4932, 
    cols_to_clean: list = None, 
    string_db_version: str = "v12.0"
) -> pl.LazyFrame:
    """

    data_id : str
        The type of dataset to download (e.g., 'protein.physical.links.detailed', 'protein.aliases', 'protein.info')
    organism_id : int
        NCBI tax ID (default: 4932 for Saccharomyces cerevisiae)
    cols_to_clean : list
        Columns where you want to strip out the organism prefix (e.g., ['protein1', 'protein2'])
    string_db_version : str
        The STRING version string (default: 'v12.0')
    """
    separator = " "
    if data_id in {"protein.aliases", "protein.info"}:
        separator = "\t"
    url = f"https://stringdb-downloads.org/download/{data_id}.{string_db_version}/{organism_id}.{data_id}.{string_db_version}.txt.gz"
    print(f"Fetching/Loading from cache: {url}")
    local_path = cached_path(url)

    lf = pl.scan_csv(str(local_path), separator=separator)

    schema_cols = lf.collect_schema().names()
    rename_map = {col: col.lstrip("#") for col in schema_cols if col.startswith("#")}
    if rename_map:
        lf = lf.rename(rename_map)

    if cols_to_clean:
        prefix_pattern = fr"^{organism_id}\."
        lf = lf.with_columns([
            pl.col(col).str.replace(prefix_pattern, "") for col in cols_to_clean
        ])

    return lf


def save_lazyframe_to_data(lf: pl.LazyFrame, filename: str, data_dir: Path | None = None) -> Path:
    if data_dir is None:
        data_dir = Path(__file__).resolve().parent.parent / "data"

    data_dir.mkdir(parents=True, exist_ok=True)
    output_path = data_dir / filename
    lf.collect().write_csv(output_path)
    print(f"Saved {output_path}")
    return output_path


if __name__ == "__main__":
    print("Testing yeast data download with Polars...")

    data_dir = Path(__file__).resolve().parent.parent / "data"
    data_dir.mkdir(parents=True, exist_ok=True)

    physical_links_lf = download_string_data(
        data_id="protein.physical.links.detailed",
        organism_id=4932,
        cols_to_clean=["protein1", "protein2"]
    )
    protein_aliases_lf = download_string_data(
        data_id="protein.aliases",
        organism_id=4932,
    )

    protein_info_lf = download_string_data(
            data_id="protein.info",
            organism_id=4932,
        )

    save_lazyframe_to_data(physical_links_lf, "4932.protein.physical.links.detailed.v12.0.txt", data_dir)
    save_lazyframe_to_data(protein_aliases_lf, "4932.protein.aliases.v12.0.txt", data_dir)
    save_lazyframe_to_data(protein_info_lf, "4932.protein.info.v12.0.txt", data_dir)

    print(protein_aliases_lf.head().collect())
    print(len(protein_aliases_lf.collect()))