try:
    infile = snakemake.input[0]
    outfile = snakemake.output[0]
except:
    infile = "SPOG/SPOG-Gos136-JRrajzcajp-rd0907171301.xml"
    outfile = "snakedir/SPOG/SPOG-Gos136-JRrajzcajp-rd0907171301.csv"

import pandas as pd
import lxml.etree as ET
from pathlib import Path

doc = ET.fromstring(Path(infile).read_bytes())
# desc_type = doc.find(".//{*}desc[@type='topic']").text
textdesc = doc.find(".//{*}textDesc").find(".//{*}domain").text

excluded_w_ids = []
if "SPOG/" in infile:
    sstpath = Path(infile.replace("SPOG/", "SST/"))
    if sstpath.exists():
        print("SST will be filtered out")
        sst = ET.fromstring(sstpath.read_bytes())
        for w in sst.findall(".//{*}w"):
            excluded_w_ids.append(w.get("{http://www.w3.org/XML/1998/namespace}id"))
        print(f"Found {len(excluded_w_ids)} to filter out.")

us = list(doc.findall(".//{*}u"))
whos = [u.get("who", "").replace("#", "") for u in us]
num_ws = [
    len(
        [
            w
            for w in u.findall(".//{*}w")
            if (w.get("{http://www.w3.org/XML/1998/namespace}id") not in excluded_w_ids)
            and (len(list(w.iter())) == 1)
        ]
        # + [w for w in u.findall(".//{*}pc")]
    )
    for u in us
]

df = pd.DataFrame(dict(who=whos, numw=num_ws)).groupby("who").numw.sum().reset_index()
df["file"] = str(Path(infile).name)
# df["desc_type"] = desc_type
df["textdesc"] = textdesc

df[df.who != ""].to_csv(outfile, index=False)
