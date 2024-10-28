from pathlib import Path
from conllu import parse
import polars as pl

df = pl.read_csv("final_split.csv")


def key(item):
    from string import digits

    sid = item.metadata["sent_id"]
    l = sid.split(".")
    l = ["".join([j for j in i if j in digits]) for i in l]
    t = tuple([float(i) if bool(i) else 0 for i in l])
    return t


results = dict(train=list(), dev=list(), test=list())

data = (
    parse(Path("../UD_Slovenian-SST/sl_sst-ud-train.conllu").read_text())
    + parse(Path("../UD_Slovenian-SST/sl_sst-ud-dev.conllu").read_text())
    + parse(Path("../UD_Slovenian-SST/sl_sst-ud-test.conllu").read_text())
)
for i in data:
    sent_id = i.metadata["sent_id"]
    destination = df.filter(pl.col("sent_id") == sent_id)["split"][0]
    results[destination].append(i)
for k, l in results.items():
    l = sorted(l, key=key)
    s = ""
    for i in l:
        s = s + i.serialize()
    Path(f"sl_sst-ud-{k}.conllu").write_text(s)

2 + 2
