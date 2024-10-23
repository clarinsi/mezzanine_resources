from pathlib import Path
import polars as pl
from conllu import parse

r = []


def get_letters_from_metadata(sent_id: str) -> str:
    recording = sent_id.split(".")[0]
    files = list(Path("../SST").glob(f"*{recording}*.xml")) + list(
        Path("../SPOG").glob(f"*{recording}*.xml")
    )
    letters = files[0].name.split("-")[2][:2]
    return letters


for splt in "train dev test".split(" "):
    data = parse(Path(f"../UD_Slovenian-SST/sl_sst-ud-{splt}.conllu").read_text())
    for i in data:
        r.append(
            {
                "wordlen": len(i),
                "sent_id": i.metadata["sent_id"],
                "speaker": i.metadata["speaker_id"],
                "old_split": splt,
                "letter": get_letters_from_metadata(i.metadata["sent_id"]),
            }
        )
df = pl.DataFrame(r).filter(~pl.col("sent_id").str.contains("Artur-"))


2 + 2
