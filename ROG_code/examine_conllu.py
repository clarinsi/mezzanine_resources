import polars as pl
from os import environ
environ["POLARS_FMT_MAX_ROWS"] = "100"
df = pl.read_csv("conllu_splits_prop_2.csv")
gb = (
    df.group_by(["split"])
    .agg(pl.col("wordlen").sum())
    .with_columns((pl.col("wordlen") / df["wordlen"].sum()).alias("ratio"))
).sort(["split"])
print(gb)
gb = df.group_by(["split", "letters"]).agg(pl.col("wordlen").sum()).sort(["split", "letters"])
print(gb.pivot(on="letters", index="split"))
gb = df.group_by(["letters"]).agg(pl.col("wordlen").sum()).sort("letters")
print(gb)
2 + 2
