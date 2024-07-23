from pathlib import Path

error_path = Path("errors.txt")
try:
    conllu_path = Path(snakemake.input.conllu)
    exb_path = Path(snakemake.input.exb)
    textgrid_path = Path(snakemake.input.pu)
    outpath = Path(snakemake.output[0])
except NameError:
    searchstring = "-J-Gvecg-P500014"
    conllu_path = list(
        Path(
            "/home/peter/mezzanine_resources/UD-SST-split/Artur-J-Gvecg-P500016.conllu"
        ).parent.glob(f"*{searchstring}*")
    )[0]
    exb_path = list(
        Path(
            "/home/peter/mezzanine_resources/Iriss-disfl-anno-phase5-fin-corr/Iriss-J-Gvecg-P500016.exb.xml"
        ).parent.glob(f"*{searchstring}*")
    )[0]
    textgrid_path = list(
        Path(
            "/home/peter/mezzanine_resources/iriss-prosodic-units/Iriss-J-Gvecg-P500016-avd_lr_sm-no-arg.TextGrid"
        ).parent.glob(f"*{searchstring}*")
    )[0]
    outpath = Path("brisi.xml")

import EXBUtils
import conllu
import textgrids

exb = EXBUtils.EXB(exb_path)
cnl = conllu.parse(conllu_path.read_text())
grid_pu = textgrids.TextGrid(textgrid_path)["PU"]

for s, tier in exb.get_traceability_tiers().items():
    for e in tier.findall(".//ud-tier-information"):
        e.getparent().remove(e)


def f(token):
    if (
        token["misc"]["Gos2.1_token_id"] == "Artur-N-G5044-P600044.s5.element1"
        or token["misc"]["Gos2.1_token_id"] == "Artur-N-G5044-P600044.s5.element2"
    ):
        return False
    return (token["upos"] == "PUNCT") or (token["form"].startswith("[name:"))


def get_starting_tokens(sentence, f=f) -> list[int]:
    indices = list()
    for i, t in enumerate(sentence):
        if f(t):
            indices.append(i)
        else:
            break
    return indices


def get_ending_tokens(sentence, f=f) -> list[int]:
    indices = list()
    elements_that_fit = []
    for i, t in enumerate(sentence[::-1]):
        if f(t):
            elements_that_fit.append(t)
        else:
            break
    for i, t in enumerate(sentence):
        if t in elements_that_fit:
            indices.append(i)
    return sorted(indices)


def get_middle_tokens(sentence, f=f) -> list[list[int]]:
    lists_of_indices = list()
    current_list = list()
    for i, t in enumerate(sentence[:-1]):
        if f(t):
            current_list.append(i)
        if not f(sentence[i + 1]):
            lists_of_indices.append(current_list)
            current_list = []
    lists_of_indices = [i for i in lists_of_indices if len(i) > 0]
    lists_of_indices = [
        i for i in lists_of_indices if i != get_starting_tokens(sentence, f=f)
    ]
    lists_of_indices = [
        i for i in lists_of_indices if i != get_ending_tokens(sentence, f=f)
    ]
    return lists_of_indices


def insert_before(token, end_id):
    time_str = exb.timeline_str[end_id]
    start_id = end_id + "_"
    # Do tli stuff
    new_tli = EXBUtils.ET.Element("tli", id=start_id, time=time_str)
    prev_tli = exb.doc.find(f".//tli[@id='{end_id}']")
    prev_tli.getparent().insert(
        prev_tli.getparent().index(prev_tli),  # Not +1, for we are prepending.
        new_tli,
    )
    exb.update_timeline()
    # Move preceeding neighbouring events around
    for event in exb.doc.findall(f".//event[@end='{end_id}']"):
        event.set("end", start_id)
    # # Let's change all annotation tiers so that they end at the new timestamp:
    # # Inhibited for consistency with prosodic units.
    # for n in [2, 3, 4]:
    #     tier_to_check = exb.get_all_nth_tiers(n=n)[speaker]
    #     for event in tier_to_check.findall(f".//event[@start='{end_id}']"):
    #         event.set("start", start_id)

    # Let's insert the new punctuation in its proper place:
    for n in [0, 1]:
        tier_to_check = exb.get_all_nth_tiers(n=n)[speaker]
        # Insert punctuation in the proper place:
        for event in tier_to_check.findall(f".//event[@start='{end_id}']"):
            ind = event.getparent().index(event)
            newevent = EXBUtils.ET.Element("event", start=start_id, end=end_id)
            newevent.text = token["form"] if n == 2 else token["misc"]["pronunciation"]
            event.getparent().insert(ind, newevent)
    # Fix traceability tiers too:
    traceability_tier = exb.get_traceability_tiers()[speaker]
    for event in traceability_tier.findall(f".//event[@start='{end_id}']"):
        # Add a new traceability event:
        ind = event.getparent().index(event)
        newevent = EXBUtils.ET.Element("event", start=start_id, end=end_id)
        newevent.text = token["misc"]["Gos2.1_token_id"] + " "
        event.getparent().insert(ind, newevent)

    return start_id


def insert_after(token, start_id):
    time_str = exb.timeline_str[start_id]
    end_id = start_id + "_"
    # Do tli stuff
    new_tli = EXBUtils.ET.Element("tli", id=end_id, time=time_str)
    prev_tli = exb.doc.find(f".//tli[@id='{start_id}']")
    prev_tli.getparent().insert(prev_tli.getparent().index(prev_tli) + 1, new_tli)
    exb.update_timeline()
    # Move neighbouring events around
    for event in exb.doc.findall(f".//event[@start='{start_id}']"):
        event.set("start", end_id)
    # # Let's change annotation tiers so that they end at the new timestamp:
    # # This was inhibited so that it's consistent with prosodic units.
    # for n in [2, 3, 4]:
    #     tier_to_check = exb.get_all_nth_tiers(n=n)[speaker]
    #     for event in tier_to_check.findall(f".//event[@end='{start_id}']"):
    #         event.set("end", end_id)
    # Let's insert the new element in its proper place:
    for n in [0, 1]:
        tier_to_check = exb.get_all_nth_tiers(n=n)[speaker]
        # Insert punctuation in the proper place:
        for event in tier_to_check.findall(f".//event[@end='{start_id}']"):
            ind = event.getparent().index(event)
            newevent = EXBUtils.ET.Element("event", start=start_id, end=end_id)
            newevent.text = token["form"] if n == 2 else token["misc"]["pronunciation"]
            event.getparent().insert(ind + 1, newevent)
    # Fix traceability tiers too:
    traceability_tier = exb.get_traceability_tiers()[speaker]
    for event in traceability_tier.findall(f".//event[@end='{start_id}']"):
        # Add a new traceability event:
        ind = event.getparent().index(event) + 1
        newevent = EXBUtils.ET.Element("event", start=start_id, end=end_id)
        newevent.text = token["misc"]["Gos2.1_token_id"] + " "
        event.getparent().insert(ind, newevent)
    return end_id


# First pass: assign pronunciation to names:
names_token_ids = []  # These tokens will be overlooked in supertoken_joining
for ii, sentence in enumerate(cnl):
    for i in range(len(sentence)):
        if sentence[i]["form"].startswith("[name:"):
            cnl[ii][i]["misc"]["pronunciation"] = cnl[ii][i]["form"]
            names_token_ids.append(sentence[i]["misc"]["Gos2.1_token_id"])
# Add punctuation and names:
for i, sentence in enumerate(cnl):
    speaker = sentence.metadata["speaker_id"]
    starting = get_starting_tokens(
        sentence,
        f=f,
    )
    middle = get_middle_tokens(
        sentence,
        f=f,
    )
    ending = get_ending_tokens(
        sentence,
        f=f,
    )
    if starting:
        next_word_token_id = EXBUtils.trimns(
            sentence[starting[-1] + 1]["misc"]["Gos2.1_token_id"]
        )
        # Find start:
        tracevent = [
            e
            for e in exb.doc.findall(".//event")
            if str(e.text).strip() == next_word_token_id
        ][0]
        end_id = tracevent.get("start")
        time_str = exb.timeline_str[end_id]
        for t in starting[::-1]:
            token = sentence[t]
            end_id = insert_before(token, end_id)
            exb.update_timeline()
    if bool(middle != [] or ending != []):
        for mlist in middle + [ending]:
            if not mlist:
                continue
            previous_word_token_id = EXBUtils.trimns(
                sentence[min(mlist) - 1]["misc"]["Gos2.1_token_id"]
            )
            tracevent = [
                e
                for e in exb.doc.findall(".//event")
                if e.text.strip() == previous_word_token_id.strip()
            ][0]
            start_id = tracevent.get("end")
            for m in mlist:
                token = sentence[m]
                start_id = insert_after(token, start_id)
                exb.update_timeline()


# Let's add sentence boundaries for every speaker:
for speaker in exb.speakers:
    tier = EXBUtils.ET.Element(
        "tier",
        attrib={
            "id": f"{speaker}sentenceId",
            "category": "sentenceId",
            "type": "a",
            "display-name": f"{speaker} [sentenceId]",
        },
    )
    for sentence in cnl:
        sent_id = sentence.metadata["sent_id"]
        if sentence.metadata["speaker_id"] != speaker:
            continue
        token_ids = [
            t["misc"]["Gos2.1_token_id"]
            for t in sentence
            if t["misc"]["Gos2.1_token_id"] not in names_token_ids
        ]
        # Let's find min-max token_ids in the exb:
        # start:
        token_id = EXBUtils.trimns(token_ids[0])
        found_events = [
            e for e in exb.doc.findall(".//event") if e.text.strip() == token_id.strip()
        ]
        assert (
            len(found_events) > 0
        ), f"Token id {token_id} was not found. Fix it or die trying."
        start = found_events[0].get("start")

        # end:

        token_id = EXBUtils.trimns(token_ids[-1])
        found_events = [
            e for e in exb.doc.findall(".//event") if e.text.strip() == token_id.strip()
        ]
        assert (
            len(found_events) > 0
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
        tiers = [
            t
            for t in exb.doc.findall(".//tier")
            if "traceability" in t.get("display-name")
        ]
        tracevent = [
            e
            for tier in tiers
            for e in tier.findall(".//event")
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


# Now we add conllu stuff
for _speaker in exb.speakers:
    for feature in "conllu lemma upos xpos feats head deprel".split():
        featuretier = EXBUtils.ET.Element(
            "tier",
            attrib={
                "id": f"{_speaker} [{feature}]",
                "category": f"{feature}",
                "type": "a",
                "display-name": f"{_speaker} [{feature}]",
                "speaker": f"{_speaker}",
            },
        )
        # featuretier.append(
        #     EXBUtils.ET.fromstring(
        #         """<ud-tier-information><ud-information attribute-name="exmaralda:hidden">true</ud-information></ud-tier-information>"""
        #     )
        # )
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
                tiers = [
                    tier
                    for tier in exb.doc.findall(".//tier")
                    if "traceability" in tier.get("display-name")
                ]
                tracevent = [
                    e
                    for tier in tiers
                    for e in tier.findall(".//event")
                    if e.text.strip() == t["misc"]["Gos2.1_token_id"]
                ][0]
                newevent = EXBUtils.ET.Element(
                    "event", start=tracevent.get("start"), end=tracevent.get("end")
                )
                newevent.text = (
                    str((dict(t),)) if feature == "conllu" else str(t[feature])
                )
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
                newevent.text = (
                    str([dict(t) for t in tokens])
                    if feature == "conllu"
                    else " ".join(
                        [
                            str(
                                t.get(
                                    feature,
                                )
                            )
                            for t in tokens
                        ]
                    )
                )
                featuretier.append(newevent)
        exb.doc.findall(".//tier")[-1].getparent().append(featuretier)


# Lettuce add Simona's prosodic units:
for _speaker in exb.speakers:
    pu_tier = EXBUtils.ET.Element(
        "tier",
        attrib={
            "id": f"{_speaker} [prosodicUnits]",
            "category": "prosodicUnits",
            "type": "a",
            "display-name": f"{_speaker} [prosodicUnits]",
        },
    )
    for interval in grid_pu:
        interval.text = interval.text.replace("#", "").strip()
        if interval.text.strip() == "":
            continue
        if not interval.text.strip().startswith("Artur"):
            continue
        start_tokenid = "Artur" + interval.text.split("Artur")[1]
        end_tokenid = "Artur" + interval.text.split("Artur")[-1]
        startracevent = [
            e
            for e in exb.doc.findall(".//event")
            if e.text.strip() == EXBUtils.trimns(start_tokenid.strip())
        ][0]
        if not _speaker in startracevent.getparent().get("display-name"):
            continue
        endtracevent = [
            e
            for e in exb.doc.findall(".//event")
            if e.text.strip() == EXBUtils.trimns(end_tokenid.strip())
        ][0]
        e = EXBUtils.ET.Element(
            "event",
            attrib={
                "start": startracevent.get("start"),
                "end": endtracevent.get("end"),
            },
        )
        e.text = interval.text
        pu_tier.append(e)
    exb.doc.find(".//tier").getparent().append(pu_tier)

# A fix for the missing sentence ID...
sentence_id_event = [
    e
    for e in exb.doc.findall(".//event")
    if str(e.text).strip() == "Artur-N-G5044-P600044.s5_reseg.1196"
]
if sentence_id_event:
    should_start_from = [
        e
        for e in exb.doc.findall(".//event")
        if str(e.text).strip() == "Artur-N-G5044-P600044.s5.element1"
    ][0].get("start")
    sentence_id_event[0].set("start", should_start_from)


speakers = exb.speakers
# Rename top two tiers with explicit suffices:
for s in speakers:
    tiers = exb.doc.findall(f".//tier[@speaker='{s}']")
    tiers[0].set("display-name", f"{s} [word]")
    tiers[1].set("display-name", f"{s} [norm]")

# Add dialog acts tiers
# Lettuce add Simona's prosodic units:
for _speaker in exb.speakers:
    for tiersuffix in "Primary Secondary".split():
        tier = EXBUtils.ET.Element(
            "tier",
            attrib={
                "id": f"{_speaker} [dialogActs{tiersuffix}]",
                "category": f"dialogActs{tiersuffix}",
                "type": "a",
                "display-name": f"{_speaker} [dialogActs{tiersuffix}]",
                "speaker": f"{_speaker}",
            },
        )  # type: ignore
        exb.doc.find(".//tier").getparent().append(tier)

# Explicitly encode speakers in tiers
for s in speakers:
    tiers = [
        tier
        for tier in exb.doc.findall(f".//tier")
        if s.strip() in tier.get("display-name")
    ]
    for tier in tiers:
        tier.set("speaker", s.strip())
# Reorder tiers so that they are as we expect them.
tier_suffices = EXBUtils.tier_suffices
order = [s + t for s in speakers for t in tier_suffices]
mapping = {o: i for i, o in enumerate(order)}
parent = exb.doc.find(".//tier").getparent()
parent[:] = sorted(
    parent, key=lambda child: mapping.get(child.get("display-name"), -100)
)

exb.save(exb.doc, outpath)
