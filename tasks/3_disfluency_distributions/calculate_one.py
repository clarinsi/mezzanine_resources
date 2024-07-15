try:
    infile = snakemake.input[0]
    outfile = snakemake.output[0]
except NameError:
    infile = "/home/peter/mezzanine_resources/Iriss-disfl-anno-phase5-fin-corr/Iriss-N-G6100-P610002.exb.xml"
    outfile = "brisi.json"


from lxml import etree as ET
import pandas as pd
from pathlib import Path

target_tiers = [
    # "[additional]",
    "[nonverbalDisfluency]",
    "[verbalDisfluency]",
    "[disfluencyStructure]",
]


def test_target(name: str) -> bool:
    return any([i in name for i in target_tiers])


doc = ET.fromstring(Path(infile).read_bytes())
tiers = [i for i in doc.findall(".//{*}tier") if test_target(i.get("display-name"))]
speakers = {i.get("speaker") for i in tiers}
timeline = {tli.get("id"): float(tli.get("time")) for tli in doc.findall(".//tli")}
results = []
for speaker in speakers:
    filtered_tiers = [i for i in tiers if i.get("speaker") == speaker]
    for tier in filtered_tiers:
        if "disfluencyStructure" in tier.get("display-name"):
            continue
            # For disfluencyStructure tier will be counted elsewhere.
        for event in tier.findall(".//event"):
            results.append(
                {
                    "type": event.text.strip(),
                    "tier": tier.get("category"),
                    "speaker": speaker,
                    "file": Path(infile).name,
                }
            )
            # Is there also a disfluency structure tier sub_event pertaining to this superevent?
            event_start = event.get("start")
            event_end = event.get("end")
            potential_timestamps = [
                k
                for k, v in timeline.items()
                if ((v >= timeline[event_start]) and (v <= timeline[event_end]))
            ]
            ds_tier = [
                i
                for i in filtered_tiers
                if "disfluencyStructure" in i.get("display-name")
            ][0]
            for s in potential_timestamps:
                for e in potential_timestamps:
                    found_event = ds_tier.find(f"event[@start='{s}'][@end='{e}']")
                    if found_event is not None:
                        results.append(
                            {
                                "type": found_event.text.strip()
                                + "_"
                                + event.text.strip(),
                                "tier": ds_tier.get("category"),
                                "speaker": speaker,
                                "file": Path(infile).name,
                            }
                        )
pd.DataFrame(data=results).to_json(
    outfile, orient="records", indent=4, force_ascii=False
)
2 + 2
