try:
    infiles = snakemake.input
    outfile = snakemake.output[0]
except NameError:
    infiles = ["Iriss-N-G5036-P600034.exb.json"]
    outfile = "stats.json"


import pandas as pd

df = pd.concat([pd.read_json(i, orient="records") for i in infiles])
df.reset_index(drop=True).to_json(outfile, orient="records", indent=4)
from pathlib import Path

for f in infiles:
    Path(f).unlink()
