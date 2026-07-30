from pathlib import Path

from scripts.load_data import resolve_data_dir


def test_resolve_data_dir_uses_repo_data_folder_when_running_from_repo_root(tmp_path):
    repo_root = tmp_path / "repo"
    repo_root.mkdir()
    (repo_root / "data").mkdir()

    resolved = resolve_data_dir(repo_root)

    assert resolved == repo_root / "data"


def test_resolve_data_dir_uses_supplied_path_when_explicit(tmp_path):
    explicit_dir = tmp_path / "custom-data"
    explicit_dir.mkdir()

    resolved = resolve_data_dir(explicit_dir)

    assert resolved == explicit_dir
