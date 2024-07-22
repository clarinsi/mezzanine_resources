from pathlib import Path
from lxml import etree as ET
import pandas as pd

tier_suffices = [
    " [word]",
    " [norm]",
    " [additional]",
    " [nonverbalDisfluency]",
    " [verbalDisfluency]",
    " [disfluencyStructure]",
    " [traceability]",
    " [sentenceId]",
    " [lemma]",
    " [upos]",
    " [xpos]",
    " [feats]",
    " [head]",
    " [deprel]",
    " [conllu]",
    " [prosodicUnits]",
]
tier_roles = [
    "word",
    "norm",
    "additional",
    "nonverbalDisfluency",
    "verbalDisfluency",
    "disfluencyStructure",
    "traceability",
    "sentence_id",
    "lemma",
    "upos",
    "xpos",
    "feats",
    "head",
    "deprel",
    "conllu",
    "prosodicUnits",
]


class EXB:
    """A class for handling EXB files."""

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
        self.timeline = self.get_timeline()
        self.timeline_str = {
            tli.get("id"): tli.get("time") for tli in self.doc.findall(".//{*}tli")
        }
        self.df = self.pandalize()

    def update_timeline(self):
        self.timeline = self.get_timeline()
        self.timeline_str = {
            tli.get("id"): tli.get("time") for tli in self.doc.findall(".//{*}tli")
        }

    def get_timeline(self):
        return {
            tli.get("id"): float(tli.get("time"))
            for tli in self.doc.findall(".//{*}tli")
        }

    def get_speakers(self) -> set[str]:
        """Get a set of all speakers present in the doc.

        :return set[str]: speakers
        """
        speakers = [
            speaker
            for t in self.tiers
            if (speaker := t.get("speaker", None)) is not None
        ]
        return set(speakers)

    def get_display_names(self) -> list[str]:
        """Get display names for all tiers

        :return list[str]: list of display names
        """
        display_names = [
            display_name
            for t in self.tiers
            if (display_name := t.get("display-name", None)) is not None
        ]
        return display_names

    def get_all_nth_tiers(self, n: int = 1) -> dict[str, ET.Element]:
        """Get a dictionary of {speaker: nth tier} for all speakers present.

        :param int n: Which tier to take, defaults to 1
        :return dict[str, ET.Element]: dictionary of nth tiers per speaker.
        """
        return {
            s: self.doc.findall(f".//{{*}}tier[@speaker='{s}']")[n]
            for s in self.speakers
        }

    def get_traceability_tiers(self) -> dict[str, ET.Element]:
        """Get a dictionary of {speaker: traceability tier}

        :return dict[str, ET.Element]: dictionary.
        """
        return {
            s: self.doc.find(f".//{{*}}tier[@display-name='{s} [traceability]']")
            for s in self.speakers
        }

    def pandalize(self) -> pd.DataFrame:
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
                2 / 2
                if len(trac_token_candidates) == 1:
                    trac_token = list(trac_token_candidates)[0].text.strip()
                elif event.text.strip().startswith("["):
                    trac_token = None
                else:
                    print("For event:")
                    EXB.print(event)
                    print(f"there are {len(trac_token_candidates)} choices")
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
        return df.sort_values("start_s").reset_index(drop=True)

    @staticmethod
    def to_str(e: ET.Element) -> str:
        return ET.tostring(e, encoding="unicode", pretty_print=True, with_tail=False)

    @staticmethod
    def print(e: ET.Element) -> None:
        print(EXB.to_str(e))

    @staticmethod
    def save(e: ET.Element, path: str | Path) -> None:
        Path(path).parents[0].mkdir(parents=True, exist_ok=True)
        Path(path).write_text(EXB.to_str(e))


def trimns(s: str) -> str:
    return (
        s.replace(".n1", "")
        .replace(".n2", "")
        .replace(".n3", "")
        .replace(".element1", "")
        .replace(".element2", "")
        .replace(".element3", "")
        .replace(".element4", "")
        .replace(".element5", "")
    )
