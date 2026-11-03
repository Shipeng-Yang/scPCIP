"""End-to-end check: the bundled example reproduces the shipped calls exactly."""
from pathlib import Path

import pandas as pd

from pcip.cli import main

HERE = Path(__file__).resolve().parents[1] / "examples"


def test_example_reproduces_expected_calls(tmp_path):
    rc = main(["predict", "--input", str(HERE / "maize_root_example.h5ad"), "--species", "Zea mays",
               "--out", str(tmp_path)])
    assert rc == 0
    got = pd.read_csv(tmp_path / "maize_root_example.predictions.tsv", sep="\t")
    exp = pd.read_csv(HERE / "expected_predictions.tsv", sep="\t")
    pd.testing.assert_frame_equal(got[exp.columns].reset_index(drop=True), exp)


def test_unsupported_organ_fails(tmp_path):
    import pytest

    with pytest.raises(SystemExit, match="not supported"):
        main(["predict", "--input", str(HERE / "maize_root_example.h5ad"), "--species", "Zea mays",
              "--organ", "stem", "--out", str(tmp_path)])


def test_annotate_matches_command_line(tmp_path):
    import anndata as ad
    import numpy as np

    import pcip

    rc = main(["predict", "--input", str(HERE / "maize_root_example.h5ad"), "--species", "Zea mays",
               "--out", str(tmp_path)])
    assert rc == 0
    cli = pd.read_csv(tmp_path / "maize_root_example.predictions.tsv", sep="\t")
    cli_scores = pd.read_csv(tmp_path / "maize_root_example.type_scores.tsv", sep="\t", index_col=0)

    adata = ad.read_h5ad(HERE / "maize_root_example.h5ad")
    obs_before = list(adata.obs.columns)
    pcip.annotate(adata, species="Zea mays")
    for col in cli.columns.drop(["cell_id", "dataset_id"]):
        np.testing.assert_array_equal(adata.obs[f"scpcip_{col}"].to_numpy(), cli[col].to_numpy(), err_msg=col)
    np.testing.assert_array_equal(adata.obsm["scpcip_type_scores"].to_numpy(), cli_scores.to_numpy())
    assert list(adata.obs.columns[: len(obs_before)]) == obs_before
    assert adata.uns["scpcip"]["n_cells"] == adata.n_obs


def test_annotate_rejects_bad_organ():
    import anndata as ad
    import pytest

    import pcip

    adata = ad.read_h5ad(HERE / "maize_root_example.h5ad")
    with pytest.raises(ValueError, match="not supported"):
        pcip.annotate(adata, species="Zea mays", organ="stem")


def test_sample_consistency_flags_disagreeing_replicates():
    import numpy as np

    from pcip.report import sample_consistency

    rng = np.random.default_rng(0)
    clusters = np.repeat(["a", "b"], 400)
    samples = np.tile(np.repeat(["s1", "s2"], 200), 2)
    same = np.where(clusters == "a", "cortex", "xylem")
    _, ok = sample_consistency(same, samples, clusters)
    assert ok["status"] == "consistent" and ok["weighted_mean_tv"] == 0.0
    flipped = np.where((clusters == "a") == (samples == "s1"), "cortex", "xylem")
    _, bad = sample_consistency(flipped, samples, clusters)
    assert bad["status"] == "not_reproducible"
    noisy = np.where(rng.random(800) < 0.9, same, "phloem")
    _, mid = sample_consistency(noisy, samples, clusters)
    assert mid["status"] == "consistent"
