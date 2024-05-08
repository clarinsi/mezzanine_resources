try:
    infile = snakemake.input[0]
    outfile = snakemake.output.csv
    out_html = snakemake.output.html
    templatepath = snakemake.params.template
except:
    infile = "snakedir/SPOG"
    outfile = "snakedir/SPOG.csv"
    out_html = "snakedir/SPOG.html"
    templatepath = "template.html"
import pandas as pd
import lxml.etree as ET
from pathlib import Path

files = list(Path(infile).glob("*.csv"))
df = pd.concat([pd.read_csv(i) for i in files])
df["speech"] = (
    df.file.str.replace(".xml", "").str.replace("SPOG-", "").str.replace("Iriss-", "")
)
speakers = pd.read_csv("Gos.TEI/Gos-speakers.tsv", sep="\t").dropna()
speeches = (
    pd.read_csv("Gos.TEI/Gos-speeches.tsv", sep="\t")
    .dropna()
    .reset_index()
    .rename(columns={"index": "speech"})
)
speeches["speech"] = speeches.speech.str.replace("Artur-", "")
df = df.merge(
    speakers["PRS-ID PERM-RESD CHILD-RESD SEX AGE SUBCORPUS".split()],
    left_on="who",
    right_on="PRS-ID",
    how="left",
)
df = df.merge(
    speeches["speech".split()],
    left_on="speech",
    right_on="speech",
    how="left",
)
c1 = df["PERM-RESD"].str.strip() == "-"
df.loc[c1, "PERM-RESD"] = df["CHILD-RESD"][c1]
c2 = df["CHILD-RESD"].str.strip() == "-"
df.loc[c2, "CHILD-RESD"] = df["PERM-RESD"][c2]
df = df.drop_duplicates().reset_index(drop=True)
df.to_csv(outfile, index=False)
2 + 2


from jinja2 import Environment, FileSystemLoader

env = Environment(loader=FileSystemLoader("."))
template = env.get_template(name=templatepath)
tip_diskurza = (
    df.groupby(["SUBCORPUS", "textdesc"])
    .agg({"file": "nunique", "numw": "sum"})
    .rename(columns={"numw": "Število besed", "file": "Število posnetkov"})
)
spol = (
    df.groupby("SEX")
    .agg({"who": "nunique", "numw": "sum"})
    .rename(columns={"who": "Število govorcev", "numw": "Število besed"})
)
starost = (
    df.groupby("AGE")
    .agg({"who": "nunique", "numw": "sum"})
    .rename(columns={"who": "Število govorcev", "numw": "Število besed"})
)
statisticnaregija = (
    df.groupby("PERM-RESD")
    .agg({"who": "nunique", "numw": "sum"})
    .rename(columns={"who": "Število govorcev", "numw": "Število besed"})
)

html = template.render(
    subcorpus=str(Path(infile).name),
    tip_diskurza=tip_diskurza,
    spol=spol,
    starost=starost,
    statisticnaregija=statisticnaregija,
)

Path(out_html).write_text(html)
