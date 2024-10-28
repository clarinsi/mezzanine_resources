from pathlib import Path
from conllu import parse
import polars as pl

df = pl.read_csv("conllu_splits.csv")


def key(item):
    from string import digits

    sid = item.metadata["sent_id"]
    l = sid.split(".")
    l = ["".join([j for j in i if j in digits]) for i in l]
    t = tuple([float(i) if bool(i) else 0 for i in l])
    return t


# data = (
#     parse(Path("../UD_Slovenian-SST/sl_sst-ud-train.conllu").read_text())
#     + parse(Path("../UD_Slovenian-SST/sl_sst-ud-train.conllu").read_text())
#     + parse(Path("../UD_Slovenian-SST/sl_sst-ud-train.conllu").read_text())
# )
r = []
for splt in "train dev test".split():
    data = parse(Path(f"../UD_Slovenian-SST/sl_sst-ud-{splt}.conllu").read_text())
    for i in data:
        r.append(
            {
                "split": splt,
                "wordlen": len(i),
                "subcorpus": i.metadata["sent_id"].split(".")[0],
            }
        )

df = (
    pl.DataFrame(r)
    .with_columns(
        pl.when(pl.col("subcorpus").str.starts_with("Gos"))
        .then(pl.lit("Gos"))
        .when(pl.col("subcorpus").str.starts_with("Artur-J"))
        .then(pl.lit("Artur-J"))
        .when(pl.col("subcorpus").str.starts_with("Artur-P"))
        .then(pl.lit("Artur-P"))
        .when(pl.col("subcorpus").str.starts_with("Artur-N"))
        .then(pl.lit("Artur-N"))
        .otherwise(None)
        .alias("subcorpus"),
    )
    .with_columns(
        pl.when(pl.col("subcorpus").str.contains("Artur"))
        .then(pl.lit("Artur"))
        .otherwise(pl.lit("Gos"))
        .alias("supersubcorpus")
    )
)
k = "split subcorpus".split()
print(df.group_by(k).agg(pl.col("wordlen").sum()).sort(k))
k = "split supersubcorpus".split()
print(df.group_by(k).agg(pl.col("wordlen").sum()).sort(k))
2 + 2
