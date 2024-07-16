from pathlib import Path
from lxml import etree as ET
import pandas as pd

tier_suffices = [
    "",
    "",
    " [additional]",
    " [nonverbalDisfluency]",
    " [verbalDisfluency]",
    " [disfluencyStructure]",
    " [traceability]",
]
tier_roles = [
    "word",
    "norm",
    "additional",
    "nonverbalDisfluency",
    "verbalDisfluency",
    "disfluencyStructure",
    "traceability",
]


class EXB:
    def __init__(self, input: bytes | str | Path):
        if isinstance(input, bytes):
            self.doc = ET.fromstring(input)
        elif isinstance(input, str) or isinstance(input, Path):
            self.doc = ET.fromstring(Path(input).read_bytes())
        else:
            raise NotImplementedError(
                f"Unknown input of type {type(input)}. Supported: bytes for file content or strings or Paths for file location."
            )
        self.tiers = self.doc.findall(".//{*}tier")
        self.speakers = self.get_speakers()
        self.timeline = {
            tli.get("id"): float(tli.get("time"))
            for tli in self.doc.findall(".//{*}tli")
        }

    def get_speakers(self) -> set[str]:
        speakers = [
            speaker
            for t in self.tiers
            if (speaker := t.get("speaker", None)) is not None
        ]
        return set(speakers)

    def get_all_nth_tiers(self, n: int = 1) -> dict[str, ET.Element]:
        return {
            s: self.doc.findall(f".//{{*}}tier[@speaker='{s}']")[n]
            for s in self.speakers
        }

    def pandalize(self, filter_pauses: bool = True) -> pd.DataFrame:
        results = []
        for s, tier in self.get_all_nth_tiers().items():
            for event in tier.findall(".//event"):
                trac_tier = self.doc.find(
                    f".//tier[@display-name='{s} [traceability]']"
                )

                trac_token_candidates = set(
                    list(trac_tier.findall(f".//event[@start='{event.get('start')}']"))
                    + list(trac_tier.findall(f".//event[@end='{event.get('end')}']"))
                )
                if len(trac_token_candidates) == 1:
                    trac_token = list(trac_token_candidates)[0]
                elif event.text.strip().startswith("["):
                    trac_token = None
                else:
                    EXB.print(event)
                    2 + 2
                results.append(
                    {
                        "speaker": s,
                        **event.attrib,
                        "norm": event.text.strip(),
                        "start_s": self.timeline.get(event.get("start")),
                        "end_s": self.timeline.get(event.get("end")),
                        "traceability": trac_token,
                    }
                )
        df = pd.DataFrame(data=results)
        if filter_pauses:
            c = df.norm.str.strip().str.startswith(
                "["
            ) & df.norm.str.strip().str.endswith("]")
            df = df[~c]
        return df.sort_values("start_s").reset_index(drop=True)

    @staticmethod
    def to_str(e: ET.Element) -> str:
        return ET.tostring(e).decode()

    @staticmethod
    def print(e: ET.Element) -> None:
        print(EXB.to_str(e))
