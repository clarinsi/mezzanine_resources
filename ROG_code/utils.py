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


rogify(
    "/home/peter/mezzanine_resources/Iriss-DA-disfl-conll-pros/Iriss-J-Gvecg-P500001.exb.xml",
    "brisi.xml",
    public=False,
)
