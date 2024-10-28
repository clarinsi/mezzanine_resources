from pathlib import Path
from conllu import parse
import polars as pl


r = list()
for splt in "train dev test".split(" "):
    data = parse(Path(f"sl_sst-ud-{splt}.conllu").read_text())
    for i in data:
        r.append(
            {
                "wordlen": len(i),
                "sent_id": i.metadata["sent_id"],
                "speaker": i.metadata["speaker_id"],
                "split": splt,
                "doc": i.metadata["sent_id"].split(".")[0]
            }
        )
# Data from current conllu
fc = pl.DataFrame(r)

r = list()
for splt in "train dev test".split(" "):
    data = parse(Path(f"../UD_Slovenian-SST/sl_sst-ud-{splt}.conllu").read_text())
    for i in data:
        r.append(
            {
                "wordlen": len(i),
                "sent_id": i.metadata["sent_id"],
                "speaker": i.metadata["speaker_id"],
                "split": splt,
                "doc": i.metadata["sent_id"].split(".")[0]
            }
        )
# Data from OLD splits
fo = pl.DataFrame(r)

# Data from file
ff = pl.read_csv("final_split.csv")




# Test that we are not losing data
assert ff.shape[0] == fc.shape[0]
assert ff.shape[0] == fo.shape[0]
assert set(fc["sent_id"].to_list()) == set(fo["sent_id"].to_list())
# Test that all documents belong to only one split
gb = fc.group_by("doc").agg(pl.col("split").n_unique().alias("in_splits"), pl.col("split").unique()).filter(pl.col("in_splits")>1)
assert gb.shape[0] == 0


# Test that we are not losing tokes:
left = fo.select(pl.col("sent_id"), pl.col("wordlen")).sort("sent_id")
right = fc.select(pl.col("sent_id"), pl.col("wordlen")).sort("sent_id")
assert left["wordlen"].to_list() == right["wordlen"].to_list()



# Test that we have all the sent_ids that were in the old data
missing = fc.filter(~pl.col("sent_id").is_in(fo["sent_id"]))
assert missing.shape[0] == 0
2+2