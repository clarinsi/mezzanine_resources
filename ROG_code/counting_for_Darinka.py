from pathlib import Path
import polars as pl

df = pl.read_csv(
    "/home/peter/mezzanine_resources/ROG/ROG-Art/ROG-Art-speakers.tsv", separator="\t"
).filter(pl.col("SUBCORPUS").is_in(["Artur-J", "Artur-P"]))

for c in "SEX AGE".split():
    print("Grouping by", c)
    print(
        df.group_by(c).agg(
            pl.col("PRS-ID").n_unique().alias("num_speakers"),
            pl.col("WORDS").sum().alias("num_words"),
        )
    )
print("Total speakers", df["PRS-ID"].n_unique(), "total num words", df["WORDS"].sum())
2 + 2
