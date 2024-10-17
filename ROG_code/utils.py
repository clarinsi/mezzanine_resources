def rogify(inpath, outpath, public=False):
    fixDA(inpath, outpath, public=public)
    prettify(outpath, outpath, public=public)


def get_tier_name(s: str) -> str:
    import re

    match = re.match(r"^(.*)\s\[.*\]$", s)
    return match.group(1).strip()


def fixDA(inpath: str, outpath: str, public=False):
    """Moves DA comments to separate tier. If public, delete comment tier and
    secondary DA tier.

    :param str inpath
    :param str outpath
    :param bool public
    """
    from pathlib import Path
    from lxml import etree as et

    doc = et.fromstring(Path(inpath).read_bytes())
    tiers = doc.findall(".//tier")
    tier_names = [i.attrib.get("display-name") for i in tiers]
    speakers = set([get_tier_name(i) for i in tier_names])
    for speaker in speakers:
        # Move comments to comment tier:
        t = doc.find(f".//tier[@display-name='{speaker} [dialogActsPrimary]']")
        attribs = {k: v.replace("Primary", "Comment") for k, v in t.attrib.items()}
        newtier = et.Element("tier", attrib=attribs)
        for event in t.findall(".//event"):
            text = event.text
            if not " " in text:
                continue
            primary, *comments = text.split(" ")
            comments = " ".join(comments)
            event.text = primary
            newevent = et.Element("event", attrib=event.attrib)
            newevent.text = comments
            newtier.append(newevent)
        t.getparent().insert(t.getparent().index(t) + 1, newtier)
        if public:
            # Delete dialogActsSecondary:
            t = doc.find(f".//tier[@display-name='{speaker} [dialogActsSecondary]']")
            t.getparent().remove(t)
            # Delete dialogActsComment
            t = doc.find(f".//tier[@display-name='{speaker} [dialogActsComment]']")
            t.getparent().remove(t)
    et.indent(doc)
    Path(outpath).write_bytes(et.tostring(doc, encoding="utf-8"))


def prettify(inpath: str, outpath: str, public=False):
    from pathlib import Path
    from lxml import etree as et
    import datetime

    doc = et.fromstring(Path(inpath).read_bytes())
    t = doc.find(".//transcription-name")
    t.text = t.text.replace("Iriss", "Rog-Art")
    t = doc.find(".//referenced-file")
    t.set(
        "url",
        str(Path("../AVD/", Path(inpath).with_suffix("").with_suffix(".wav"))).replace(
            "Iriss", "RogArt"
        ),
    )
    t = doc.find(".//comment")
    t.text = (
        t.text
        + f" \n\n {datetime.datetime.now(datetime.UTC).isoformat()}, PR: Import annotations from other sources (Dialog Acts, Prosodic Units, CONLLU...) and prepare for release. Release: {'public' if public else 'internal'}."
    )

    # Find hidden tiers and unhide:
    ts = doc.findall(".//ud-information[@attribute-name='exmaralda:hidden']")
    for t in ts:
        tt = t.getparent()
        tt.getparent().remove(tt)
    ts = doc.findall(".//ud-information[@attribute-name='AutoSave']")
    for t in ts:
        tt = t.getparent()
        tt.getparent().remove(tt)

    et.indent(doc)
    Path(outpath).write_bytes(et.tostring(doc, encoding="utf-8"))


# rogify(
#     "/home/peter/mezzanine_resources/Iriss-DA-disfl-conll-pros/Iriss-J-Gvecg-P500001.exb.xml",
#     "brisi.xml",
#     public=False,
# )
def get_conll_ids():
    from pathlib import Path
    from conllu import parse

    ids = []
    for splt in "train dev test".split(" "):
        data = parse(Path(f"../UD_Slovenian-SST/sl_sst-ud-{splt}.conllu").read_text())
        ids.extend(list(set([i.metadata["sent_id"].split(".")[0] for i in data])))
    return ids


# get_conll_ids()


def make_conll_splits():
    from pathlib import Path
    from conllu import parse
    from subprocess import run

    r = []
    for splt in "train dev test".split(" "):
        data = parse(Path(f"../UD_Slovenian-SST/sl_sst-ud-{splt}.conllu").read_text())
        for i in data:
            r.append({"Sent_id": i.metadata["sent_id"], "Split": splt})
    return r


def in_corpus(s: str, corpus: str) -> bool:
    from subprocess import run

    r = run(
        f"""cat ../{corpus}/*.xml | grep '{s.strip()}"'""",
        shell=True,
        capture_output=True,
    )
    return r.returncode == 0


def do_conllus():
    import polars as pl
    from pathlib import Path
    from conllu import parse
    from string import digits

    r = make_conll_splits()
    df = pl.DataFrame(r).with_columns(
        pl.col("Sent_id")
        .map_elements(lambda i: in_corpus(i, "SST"), return_dtype=pl.Boolean)
        .alias("in_sst"),
        pl.col("Sent_id")
        .map_elements(lambda i: in_corpus(i, "SPOG"), return_dtype=pl.Boolean)
        .alias("in_spog"),
        pl.col("Sent_id").str.contains("Artur").alias("in_artur"),
        pl.col("Sent_id").str.split(".").list[0].alias("file"),
    )
    print("SST and SPOG segments assigned.")

    def key(item):
        sid = item.metadata["sent_id"]
        l = sid.split(".")
        l = ["".join([j for j in i if j in digits]) for i in l]
        t = tuple([float(i) if bool(i) else 0 for i in l])
        return t

    data = (
        parse(Path("../UD_Slovenian-SST/sl_sst-ud-train.conllu").read_text())
        + parse(Path("../UD_Slovenian-SST/sl_sst-ud-train.conllu").read_text())
        + parse(Path("../UD_Slovenian-SST/sl_sst-ud-train.conllu").read_text())
    )
    print("Conllu parsed.")
    from tqdm import tqdm

    for file in tqdm(df["file"].unique(), total=df["file"].unique().shape[0]):
        # SST -> GO1
        subset = df.filter((pl.col("file") == file) & (pl.col("in_sst") == True))
        if subset.shape[0] != 0:
            subdata = [i for i in data if i.metadata["sent_id"] in subset["Sent_id"]]
            subdata = sorted(subdata, key=key)

            path = Path(f"../ROG/CONLLU/Rog-Go1-{file}.conllu")
            path.parent.mkdir(exist_ok=True)
            path.write_text("\n".join([i.serialize() for i in data]))
            # print("wrote")
        # SPOG -> GO2
        subset = df.filter((pl.col("file") == file) & (pl.col("in_spog") == True))
        if subset.shape[0] != 0:
            subdata = [i for i in data if i.metadata["sent_id"] in subset["Sent_id"]]
            subdata = sorted(subdata, key=key)

            path = Path(f"../ROG/CONLLU/Rog-Go2-{file}.conllu")
            path.parent.mkdir(exist_ok=True)
            path.write_text("\n".join([i.serialize() for i in data]))
            # print("wrote")
        # Artur -> Rog-Art
        subset = df.filter((pl.col("file") == file) & (pl.col("in_artur") == True))
        if subset.shape[0] != 0:
            subdata = [i for i in data if i.metadata["sent_id"] in subset["Sent_id"]]
            subdata = sorted(subdata, key=key)

            path = Path(f"../ROG/CONLLU/Rog-Art{file.replace('Artur-', '-')}.conllu")
            path.parent.mkdir(exist_ok=True)
            path.write_text("\n".join([i.serialize() for i in data]))
            # print("wrote")
    2 + 2


def fix_trs(inpath: str, outpath: str):
    from pathlib import Path
    from lxml import etree as et
    import datetime

    doc = et.fromstring(Path(inpath).read_bytes())

    doc.set(
        "audio_filename",
        "../AVD/"
        + str(Path(inpath).with_suffix("").with_suffix(".wav").name)
        .replace("Iriss", "Rog-Art")
        .replace("-pog", "")
        .replace("-std", ""),
    )
    import datetime

    doc.set("version_date", f"{datetime.datetime.now(datetime.UTC).isoformat()}")
    doc.set("version", str(int(doc.get("version")) + 1))
    et.indent(doc)
    Path(outpath).write_bytes(
        et.tostring(
            doc,
            encoding="UTF-8",
            xml_declaration=True,
            doctype='<!DOCTYPE Trans SYSTEM "trans-14.dtd">',
        )
    )


# fix_trs(
#     "/home/peter/mezzanine_resources/ROG/ROG-Art/TRS/Rog-Art-J-Gvecg-P500001-pog.trs",
#     "brisi.xml",
# )
# do_conllus()
