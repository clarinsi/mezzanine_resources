from pathlib import Path

original_gos_tei = Path(
    "/home/peter/mezzanine_resources/Gos.TEI/Artur-J/Artur-J-Gvecg-P580002.xml"
)
conllu_path = Path(
    "/home/peter/mezzanine_resources/UD-SST-split/Artur-J-Gvecg-P580002.conllu"
)
exb_path = Path(
    "/home/peter/mezzanine_resources/Iriss-disfl-anno-phase5-fin-corr/Iriss-J-Gvecg-P580002.exb.xml"
)
outpath = Path("brisi.xml")

import EXBUtils
import conllu

exb = EXBUtils.EXB(exb_path)
cnl = conllu.parse(conllu_path.read_text())


# Add punctuation:
for sentence in cnl:
    speaker = sentence.metadata["speaker_id"]
    # Is the first token punctuation?
    if sentence[0]["upos"] == "PUNCT":
        raise NotImplementedError(
            "This was not anticipated... Who starts a sentence with a punctuation? ¿What are you, Spanish?"
        )

    for previous, current in zip(sentence, sentence[1:]):
        if current["upos"] == "PUNCT":
            punctuation = current["form"]
            previous_token_id = EXBUtils.trimns(previous["misc"]["Gos2.1_token_id"])
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
                for event in tier_to_check.findall(f".//event[@start='{end_id}']"):
                    event.set("start", end_id + "_punct")

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

            # If no traceability corrections were performed, then we also didn't introduce a new event.
            # Lettuce fix that:
            if (
                len(
                    traceability_tier.findall(
                        f".//event[@start='{end_id}'][@end='{end_id}_punct']"
                    )
                )
                == 0
            ):
                event = traceability_tier.findall(f".//event[@end='{end_id}']")[0]
                ind = event.getparent().index(event)
                newevent = EXBUtils.ET.Element(
                    "event", start=end_id, end=end_id + "_punct"
                )
                newevent.text = current["misc"]["Gos2.1_token_id"] + " "
                event.getparent().insert(ind, newevent)

# Let's add sentence boundaries for every speaker:
for speaker in exb.speakers:
    tier = EXBUtils.ET.Element(
        "tier",
        attrib={
            "id": f"{speaker}SENT",
            "category": "sentenceid",
            "type": "a",
            "display-name": f"{speaker} [sentenceid]",
        },
    )
    for sentence in cnl:
        sent_id = sentence.metadata["sent_id"]
        if sentence.metadata["speaker_id"] != speaker:
            continue
        token_ids = [t["misc"]["Gos2.1_token_id"] for t in sentence]
        # Let's find min-max token_ids in the exb:
        # start:
        token_id = token_ids[0]
        found_events = [
            e for e in exb.doc.findall(".//event") if e.text.strip() == token_id.strip()
        ]
        assert (
            len(found_events) == 1
        ), f"Token id {token_id} was not found. Fix it or die trying."
        start = found_events[0].get("start")

        # end:

        token_id = token_ids[-1]
        found_events = [
            e for e in exb.doc.findall(".//event") if e.text.strip() == token_id.strip()
        ]
        assert (
            len(found_events) == 1
        ), f"Token id {token_id} was not found. Fix it or die trying."
        end = found_events[0].get("end")
        newevent = EXBUtils.ET.Element("event", start=start, end=end)
        newevent.text = sent_id
        tier.append(newevent)
    exb.doc.findall(".//tier")[-1].getparent().append(tier)


# Fixing words and norms:
for sentence in cnl:
    speaker = sentence.metadata["speaker_id"]
    # Let's test if we have split words present:
    tokids = [EXBUtils.trimns(t["misc"]["Gos2.1_token_id"]) for t in sentence]
    from collections import Counter

    c = Counter(tokids)

    # First we fix "monowords", so tokens that map 1:1 to EXB tokens.
    monoword_tokens = [
        t for t in sentence if c[EXBUtils.trimns(t["misc"]["Gos2.1_token_id"])] == 1
    ]
    for t in monoword_tokens:
        tracevent = [
            e
            for e in exb.doc.findall(".//event")
            if e.text.strip() == t["misc"]["Gos2.1_token_id"].strip()
        ][0]
        word_event = exb.get_all_nth_tiers(0)[speaker].find(
            f".//event[@start='{tracevent.get('start')}'][@end='{tracevent.get('end')}']"
        )
        assert word_event is not None
        word_event.text = t["misc"]["pronunciation"]

        norm_event = exb.get_all_nth_tiers(1)[speaker].find(
            f".//event[@start='{tracevent.get('start')}'][@end='{tracevent.get('end')}']"
        )
        assert norm_event is not None
        norm_event.text = t["form"]

    # And now we move to multiword tokens, where we map N words from conllu to 1 cell in EXB.
    multiword_tokens = [
        t for t in sentence if c[EXBUtils.trimns(t["misc"]["Gos2.1_token_id"])] > 1
    ]
    supertokenids = set(
        [EXBUtils.trimns(t["misc"]["Gos2.1_token_id"]) for t in multiword_tokens]
    )
    for sti in supertokenids:
        tokens = [
            t for t in sentence if EXBUtils.trimns(t["misc"]["Gos2.1_token_id"]) == sti
        ]
        tracevent = [
            e
            for e in exb.doc.findall(".//event")
            if e.text.strip()
            == EXBUtils.trimns(tokens[0]["misc"]["Gos2.1_token_id"]).strip()
        ][0]
        word_event = exb.get_all_nth_tiers(0)[speaker].find(
            f".//event[@start='{tracevent.get('start')}'][@end='{tracevent.get('end')}']"
        )
        assert word_event is not None
        word_event.text = " ".join([t["misc"]["pronunciation"] for t in tokens[0:1]])

        norm_event = exb.get_all_nth_tiers(1)[speaker].find(
            f".//event[@start='{tracevent.get('start')}'][@end='{tracevent.get('end')}']"
        )
        assert norm_event is not None
        norm_event.text = " ".join([t["form"] for t in tokens])


# Now we automate it for upos, lemma, feats, head, deprel and deps
for _speaker in exb.speakers:
    for feature in ["upos", "lemma", "feats", "head", "deprel", "misc"]:
        featuretier = EXBUtils.ET.Element(
            "tier",
            attrib={
                "id": f"{_speaker} [{feature.upper()}]",
                "category": feature.upper(),
                "type": "a",
                "display-name": f"{_speaker} [{feature.upper()}]",
            },
        )

        for sentence in cnl:
            speaker = sentence.metadata["speaker_id"]
            if speaker != _speaker:
                continue
            # Let's test if we have split words present:
            tokids = [EXBUtils.trimns(t["misc"]["Gos2.1_token_id"]) for t in sentence]
            from collections import Counter

            c = Counter(tokids)

            # First we fix "monowords", so tokens that map 1:1 to EXB tokens.
            monoword_tokens = [
                t
                for t in sentence
                if c[EXBUtils.trimns(t["misc"]["Gos2.1_token_id"])] == 1
            ]
            for t in monoword_tokens:
                tracevent = [
                    e
                    for e in exb.doc.findall(".//event")
                    if e.text.strip() == t["misc"]["Gos2.1_token_id"].strip()
                ][0]
                newevent = EXBUtils.ET.Element(
                    "event", start=tracevent.get("start"), end=tracevent.get("end")
                )
                newevent.text = str(t[feature])
                featuretier.append(newevent)

            # And now we move to multiword tokens, where we map N words from conllu to 1 cell in EXB.
            multiword_tokens = [
                t
                for t in sentence
                if c[EXBUtils.trimns(t["misc"]["Gos2.1_token_id"])] > 1
            ]
            supertokenids = set(
                [
                    EXBUtils.trimns(t["misc"]["Gos2.1_token_id"])
                    for t in multiword_tokens
                ]
            )
            for sti in supertokenids:
                tokens = [
                    t
                    for t in sentence
                    if EXBUtils.trimns(t["misc"]["Gos2.1_token_id"]) == sti
                ]
                tracevent = [
                    e
                    for e in exb.doc.findall(".//event")
                    if e.text.strip()
                    == EXBUtils.trimns(tokens[0]["misc"]["Gos2.1_token_id"]).strip()
                ][0]
                newevent = EXBUtils.ET.Element(
                    "event", start=tracevent.get("start"), end=tracevent.get("end")
                )
                newevent.text = " ".join([str(t[feature]) for t in tokens])
                featuretier.append(newevent)
        exb.doc.findall(".//tier")[-1].getparent().append(featuretier)

# 2 + 2


exb.save(exb.doc, outpath)

2 + 2


# for sentence in cnl:
#     sent_id = sentence.metadata["sent_id"]
#     token_ids = [t["misc"]["Gos2.1_token_id"] for t in sentence]

2 + 2
