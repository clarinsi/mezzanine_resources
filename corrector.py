import pandas as pd
from pathlib import Path
from conllu import parse
from lxml import etree
from numpy import median
from loguru import logger

files = (
    pd.read_html("/home/peter/mezzanine_resources/rog_salvage/LONG_XPOS-1.html")[0]
    .Comm.unique()
    .tolist()
)
# files = [
#     i.with_suffix("").with_suffix("").name
#     for i in Path(f"/home/peter/mezzanine_resources/ROG/ROG-Art/EXB/").glob("*.exb.xml")
# ]
for file in files:
    conllufile = Path(f"/home/peter/mezzanine_resources/ROG/CONLLU/{file}.conllu")
    exbfile = Path(f"/home/peter/mezzanine_resources/ROG/ROG-Art/EXB/{file}.exb.xml")
    outpath = exbfile
    conlludata = parse(Path(conllufile).read_text())
    doc = etree.fromstring(Path(exbfile).read_bytes())
    logger.debug(f"Opening file {Path(exbfile).name}")
    timeline = {i.get("id"): float(i.get("time")) for i in doc.findall(".//tli")}
    tiers = doc.findall(".//tier")
    traceability_tiers = [i for i in tiers if "[traceability]" in i.get("display-name")]
    lens_of_events = []
    for t in traceability_tiers:
        for e in t.findall(".//event"):
            lens_of_events.append(len(e.text))

    medianlentgh = median(lens_of_events)

    suspect_events = []
    for t in traceability_tiers:
        for e in t.findall(".//event"):
            if len(e.text) > 1.5 * medianlentgh or " " in e.text.strip():
                suspect_events.append(e)
                logger.info(f"Found traceability: {e.text}")
    for e in suspect_events:
        try:
            speaker = e.getparent().get("speaker")
            wordtier = doc.find(f".//tier[@display-name='{speaker} [word]']")
            assert wordtier is not None
            start_idx = e.get("start")
            end_idx = e.get("end")
            top_tier_events = [
                i
                for i in wordtier.findall(".//event")
                if (timeline[i.get("start")] >= timeline[start_idx])
                and (timeline[i.get("end")] <= timeline[end_idx])
            ]
            traceability_tokens = e.text.strip().split()
            assert len(traceability_tokens) == len(top_tier_events)
            insert_index = e.getparent().index(e)
            for i, (top_tier_event, traceability_token) in enumerate(
                zip(top_tier_events, traceability_tokens)
            ):
                newevent = etree.Element("event", attrib=top_tier_event.attrib)
                newevent.text = traceability_token
                e.getparent().insert(insert_index + 1 + i, newevent)
            fixalso = "lemma upos xpos feats head deprel conllu"
            for tier_designator in fixalso.split():
                tier = doc.find(
                    f".//tier[@display-name='{speaker} [{tier_designator}]']"
                )
                assert tier is not None
                event_to_drop = tier.find(
                    f".//event[@start='{e.get('start')}'][@end='{e.get("end")}']"
                )
                if event_to_drop is not None:
                    event_to_drop.getparent().remove(event_to_drop)
                for token in traceability_tokens:
                    i = [
                        d
                        for sent in conlludata
                        for d in sent
                        if d["misc"]["Gos2.1_token_id"] == token
                    ]
                    assert len(i) == 1
                    i = i[0]
                    if tier_designator == "conllu":
                        payload = i.__repr__()
                    else:
                        payload = i[tier_designator]
                    event_with_the_same_timestamps = [
                        i
                        for i in doc.findall(".//event")
                        if i.text.strip() == token.strip()
                    ][0]
                    newevent = etree.Element(
                        "event", attrib=event_with_the_same_timestamps.attrib
                    )
                    newevent.text = str(payload)
                    tier.append(newevent)
            e.getparent().remove(e)
        except RuntimeWarning:
            logger.warning(f"Unfixable token id: {e.text}")
            continue
    for tier in doc.findall(".//tier"):
        tier[:] = sorted(
            tier,
            key=lambda elem: timeline.get(elem.attrib.get("start", "notfound"), -1),
        )
    etree.indent(doc)
    Path(outpath).write_bytes(etree.tostring(doc, encoding="utf-8"))
