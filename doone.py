try:
    infile = snakemake.input[0]
    outfile = snakemake.output[0]
except:
    infile = "SST/SPOG-Gos164-JRtvslpovp-rd0904201545.xml"
    outfile = "snakedir/SST/SPOG-Gos164-JRtvslpovp-rd0904201545.csv"

import pandas as pd
import lxml.etree as ET
from pathlib import Path

doc = ET.fromstring(Path(infile).read_bytes())
# desc_type = doc.find(".//{*}desc[@type='topic']").text
textdesc = doc.find(".//{*}textDesc").find(".//{*}domain").text


us = doc.findall(".//{*}u")
whos = [u.get("who", "").replace("#", "") for u in us]
num_ws = [len(list(u.findall(".//{*}w"))) for u in us]

df = pd.DataFrame(dict(who=whos, numw=num_ws)).groupby("who").numw.sum().reset_index()
df["file"] = str(Path(infile).name)
# df["desc_type"] = desc_type
df["textdesc"] = textdesc

df[df.who != ""].to_csv(outfile, index=False)
