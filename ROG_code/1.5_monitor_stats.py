from pathlib import Path
import polars as pl
from conllu import parse
from os import environ

environ["POLARS_FMT_MAX_ROWS"] = "15"
artursplit = (
    pl.read_csv("/home/peter/mezzanine_resources/ROG-Artur-train-dev-test-split.csv")
    .with_columns(pl.col("Split").str.to_lowercase().alias("split"))
    .select(pl.exclude("Split"))
)
artur_metadata = pl.read_csv(
    "../Gos_2.1_ Speech-to-Discourse-Type - Speech-to-type.tsv", separator="\t"
)
artur_metadata = pl.concat(
    [
        artur_metadata,
        artur_metadata.with_columns(
            pl.col("TEXT-ID").str.replace(r"-G(\d+-)", r"-Gvecg-")
        ),
    ]
)


def get_letters_from_metadata(sent_id: str) -> str:
    recording = sent_id.split(".")[0]
    if "Artur" in sent_id:
        return artur_metadata.filter(pl.col("TEXT-ID") == recording).unique()["TYPE"][0]
    files = list(Path("../../SST").glob(f"*{recording}*.xml")) + list(
        Path("../../SPOG").glob(f"*{recording}*.xml")
    )
    letters = files[0].name.split("-")[2]
    return letters

from joblib import Memory
memory = Memory("brisi")
@memory.cache
def in_corpus(s: str, corpus: str) -> bool:
    from subprocess import run

    r = run(
        f"""cat ../../{corpus}/*.xml | grep '{s.strip()}"'""",
        shell=True,
        capture_output=True,
    )
    return r.returncode == 0


r = list()
for splt in "train dev test".split(" "):
    data = parse(Path(f"../../UD_Slovenian-SST/sl_sst-ud-{splt}.conllu").read_text())
    for i in data:
        r.append(
            {
                "wordlen": len(i),
                "sent_id": i.metadata["sent_id"],
                "speaker": i.metadata["speaker_id"],
                "old_split": splt,
                "in_sst": in_corpus(i.metadata["sent_id"], "SST"),
                "in_spog": in_corpus(i.metadata["sent_id"], "SPOG"),
            }
        )

data = (
    pl.DataFrame(r)
    .with_columns(
        pl.col("sent_id")
        .map_elements(get_letters_from_metadata, return_dtype=pl.String)
        .alias("old_name")
    )
    .with_columns(
        pl.col("old_name").str.head(2).alias("letters"),
        pl.col("sent_id").str.split(".").list[0].alias("doc"),
    )
)
data = data.join(artursplit, how="left", left_on="doc", right_on="Recording ID")

assigned = pl.concat(
    [
        pl.read_csv("test_candidates.csv").drop_nulls(subset="split"),
        pl.read_csv("dev_candidates.csv").drop_nulls(subset="split"),
    ]
)
print("Currently assigned:")
k = "split letters".split()
print(assigned.group_by(k).agg(pl.col("wordlen").sum()).sort(k))
data = (
    data.join(assigned.select(pl.col("doc"), pl.col("split")), how="left", on="doc")
    .with_columns(
        split=pl.when(pl.col("split").is_not_null())
        .then(pl.col("split"))
        .otherwise(pl.col("split_right"))
    )
    .drop(pl.col("split_right"))
)


data = data.with_columns(
    split=pl.when(pl.col("split").is_null())
    .then(pl.lit("train"))
    .otherwise(pl.col("split"))
)
k = "split letters".split()
gb = data.group_by(k).agg(pl.col("wordlen").sum()).sort(k)
print("Current split distribution:")
print(gb)




print("Overall result:\n", data.group_by("split").agg(pl.col("wordlen").sum()))

k = "letters"
gb = data.group_by(k).agg(pl.col("wordlen").sum()).sort(k)
print("Current split distribution:")
print(gb)


print("Final stats:")
print("Overall:")
gb = data.group_by("letters").agg(pl.col("wordlen").sum(),
                                (pl.col("wordlen").sum()/data["wordlen"].sum()).alias("ratio")).sort("letters")
print(gb)

for split in "dev test train".split():
    print(f"For split {split}:")
    gb = data.filter(pl.col("split")==split).group_by("letters").agg(pl.col("wordlen").sum(),
                                (pl.col("wordlen").sum()/data.filter(pl.col("split")==split)["wordlen"].sum()).alias("ratio")).sort("letters")
    print(gb)
data.write_csv("final_split.csv")

2 + 2
