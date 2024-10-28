from pathlib import Path
import polars as pl
from conllu import parse
from os import environ

environ["POLARS_FMT_MAX_ROWS"] = "20"
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

memory = Memory("../brisi")


r = []


@memory.cache
def in_corpus(s: str, corpus: str) -> bool:
    from subprocess import run

    r = run(
        f"""cat ../{corpus}/*.xml | grep '{s.strip()}"'""",
        shell=True,
        capture_output=True,
    )
    return r.returncode == 0


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
    .filter(~pl.col("sent_id").str.contains("Artur-"))
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
darinka = (
    data.filter(pl.col("doc").str.contains("Gos"))
    .with_columns(pl.col("doc").str.replace("Gos", "").cast(pl.Int16).alias("gos_id"))
    .filter(
        pl.col("gos_id").is_between(lower_bound=160, upper_bound=278)
        & (~pl.col("gos_id").is_in([217, 218]))
    )
)

gb = (
    darinka.group_by("doc")
    .agg(pl.col("wordlen").sum(), pl.col("letters").first())
    .with_columns(split=None)
)


test = pl.read_csv("test_candidates.csv").drop_nulls("split")
darinka = data.filter(
    pl.col("doc").str.contains("Gos") & ~pl.col("doc").is_in(test["doc"])
).with_columns(pl.col("doc").str.replace("Gos", "").cast(pl.Int16).alias("gos_id"))

gb = (
    darinka.group_by("doc")
    .agg(pl.col("wordlen").sum(), pl.col("letters").first())
    .with_columns(split=None)
)
gb.write_csv("dev_candidates.csv")
