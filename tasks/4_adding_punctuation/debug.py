from pathlib import Path

original_gos_tei = Path(
    "/home/peter/mezzanine_resources/Gos.TEI/Artur-J/Artur-J-Gvecg-P500002.xml"
)
conllu_path = Path(
    "/home/peter/mezzanine_resources/UD-SST-split/Artur-J-Gvecg-P500002.conllu"
)
exb_path = Path(
    "/home/peter/mezzanine_resources/Iriss-disfl-anno-phase5-fin-corr/Iriss-J-Gvecg-P500002.exb.xml"
)
outpath = Path("brisi.xml")

import EXBUtils
import conllu

exb = EXBUtils.EXB(exb_path)
cnl = conllu.parse(conllu_path.read_text())


# Add punctuation:
punctuation_to_be_added = [
    i for sentence in cnl for i in sentence if i["upos"] == "PUNCT"
]

for sentence in cnl:
    speaker = sentence.metadata["speaker_id"]
    # Is the first token punctuation?
    if sentence[0]["upos"] == "PUNCT":
        raise NotImplementedError(
            "This was not anticipated... Who starts a sentence with a punctuation? ¿What are you, Spanish?"
        )
    else:
        for previous, current in zip(sentence, sentence[1:]):
            if current["upos"] == "PUNCT":
                punctuation = current["form"]
                previous_token_id = previous["misc"]["Gos2.1_token_id"]
                traceability_event = [
                    e
                    for e in exb.doc.findall(".//event")
                    if e.text.strip() == previous_token_id.strip()
                ][0]
                end_id = traceability_event.get("end")
                time_s = exb.timeline.get(end_id)
                time_str = exb.timeline_str.get(end_id)

                # Let's insert a new tli:
                previous_tli = exb.doc.find(f".//tli[@id='{end_id}']")
                ind = previous_tli.getparent().index(previous_tli)
                previous_tli.getparent().insert(
                    ind + 1,
                    EXBUtils.ET.Element("tli", id=end_id + "_punct", time=time_str),
                )
                # Let's change all annotation tiers so that they end at the new timestamp:
                for n in [2, 3, 4]:
                    tier_to_check = exb.get_all_nth_tiers(n=n)[speaker]
                    for event in tier_to_check.findall(f".//event[@end='{end_id}']"):
                        event.set("end", end_id + "_punct")

                # Let's insert the new punctuation in its proper place:
                for n in [0, 1]:
                    tier_to_check = exb.get_all_nth_tiers(n=n)[speaker]
                    for event in tier_to_check.findall(f".//event[@end='{end_id}']"):
                        ind = event.getparent().index(event)
                        newevent = EXBUtils.ET.Element(
                            "event", start=end_id, end=end_id + "_punct"
                        )
                        newevent.text = punctuation
                        event.getparent().insert(ind + 1, newevent)

                    # Fix also the next item:
                    for event in tier_to_check.findall(f".//event[@start='{end_id}']"):
                        if event.get("end") == end_id + "_punct":
                            continue
                        event.set("start", end_id + "_punct")
                # Fix traceability tiers too:
                traceability_tier = exb.get_traceability_tiers()[speaker]
                for event in traceability_tier.findall(f".//event[@start='{end_id}']"):
                    event.set("start", end_id + "_punct")
                    # Add a new traceability event:
                    ind = event.getparent().index(event)
                    newevent = EXBUtils.ET.Element(
                        "event", start=end_id, end=end_id + "_punct"
                    )
                    newevent.text = current["misc"]["Gos2.1_token_id"] + " "
                    event.getparent().insert(ind, newevent)

                2 + 2
exb.save(exb.doc, outpath)

# for sentence in cnl:
#     sent_id = sentence.metadata["sent_id"]
#     token_ids = [t["misc"]["Gos2.1_token_id"] for t in sentence]

2 + 2
