import pytest
from pathlib import Path
from conllu import parse
import polars as pl


@pytest.fixture
def current_conllu_data():
    r = []
    for splt in ["train", "dev", "test"]:
        data = parse(Path(f"sl_sst-ud-{splt}.conllu").read_text())
        for i in data:
            r.append(
                {
                    "wordlen": len(i),
                    "sent_id": i.metadata["sent_id"],
                    "speaker": i.metadata["speaker_id"],
                    "split": splt,
                    "doc": i.metadata["sent_id"].split(".")[0],
                }
            )
    return pl.DataFrame(r)


@pytest.fixture
def old_conllu_data():
    r = []
    for splt in ["train", "dev", "test"]:
        data = parse(
            Path(f"/home/peter/UD_Slovenian-SST/sl_sst-ud-{splt}.conllu").read_text()
        )
        for i in data:
            r.append(
                {
                    "wordlen": len(i),
                    "sent_id": i.metadata["sent_id"],
                    "speaker": i.metadata["speaker_id"],
                    "split": splt,
                    "doc": i.metadata["sent_id"].split(".")[0],
                }
            )
    return pl.DataFrame(r)


@pytest.fixture
def final_split_data():
    return pl.read_csv("final_split.csv")


def test_shape_data_vs_metadata(final_split_data, current_conllu_data):
    assert final_split_data.shape[0] == current_conllu_data.shape[0]


def test_shape_data_vs_old_data(current_conllu_data, old_conllu_data):
    assert current_conllu_data.shape[0] == old_conllu_data.shape[0]


def test_sent_id_consistency(current_conllu_data, old_conllu_data):
    assert set(current_conllu_data["sent_id"].to_list()) == set(
        old_conllu_data["sent_id"].to_list()
    )


def test_document_unique_split(current_conllu_data):
    gb = (
        current_conllu_data.group_by("doc")
        .agg(pl.col("split").n_unique().alias("in_splits"), pl.col("split").unique())
        .filter(pl.col("in_splits") > 1)
    )
    assert gb.shape[0] == 0


def test_token_count_consistency(current_conllu_data, old_conllu_data):
    left = old_conllu_data.select(pl.col("sent_id"), pl.col("wordlen")).sort("sent_id")
    right = current_conllu_data.select(pl.col("sent_id"), pl.col("wordlen")).sort(
        "sent_id"
    )
    assert left["wordlen"].to_list() == right["wordlen"].to_list()


def test_missing_sent_ids(current_conllu_data, old_conllu_data):
    missing = current_conllu_data.filter(
        ~pl.col("sent_id").is_in(old_conllu_data["sent_id"])
    )
    assert missing.shape[0] == 0


def test_nonleakage(current_conllu_data):
    df = current_conllu_data
    docs_occur = (
        df.group_by("doc")
        .agg(pl.col("split").n_unique())
        .select(pl.col("split").max())["split"]
        .to_list()
    )
    assert docs_occur == [1]


def test_token_counts(current_conllu_data):
    counts = (
        current_conllu_data.with_columns(
            pl.when(pl.col("doc").str.contains("Artur"))
            .then(pl.lit("Artur"))
            .otherwise(
                pl.when(pl.col("doc").str.contains("Gos"))
                .then(pl.lit("Gos"))
                .otherwise(None)
            )
            .alias("subcorpus")
        )
        .group_by(["subcorpus", "split"])
        .agg(pl.col("wordlen").sum())
        .pivot(on="subcorpus", index="split")
        .sort("split")
    )
    print("Counts:", counts)
