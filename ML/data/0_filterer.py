import pandas as pd
import json
from pathlib import Path

data = [json.loads(i) for i in Path("prosodic_units.jsonl").read_text().splitlines()]
TARGET = "prosodicUnit"

# Ballancing stuff
REPEAT_TARGETS = False  # Should the segments with positive examples be repeated?
REPEAT_FACTOR = 5  # How many times should the positive segments be repeated?


new_data = []
for line in data:
    file = line["file"]
    for content in line["strips"]:
        new_data.append({"file": file, **content})
df = pd.DataFrame(data=new_data).round(3)
df = df[df.duration <= 30]
df["fileroot"] = df.file.apply(lambda s: str(Path(s).with_suffix("").with_suffix("")))
df["wavpath"] = [f"iriss_wavs/{i}-avd.wav" for i in df.fileroot]
assert df.wavpath.apply(lambda s: Path(s).exists()).all(), "Files are missing!"
df["segment_path"] = df.apply(
    lambda row: f"segments30s/{Path(row['wavpath']).with_suffix('').name}_{row['start_ms']}_{row['end_ms']}.wav",
    axis="columns",
)


def assign_labels(row):
    labels = []
    L0, L1 = "0", TARGET
    i = pd.interval_range(start=row["start_ms"], end=row["end_ms"], freq=0.02)
    df = pd.DataFrame(
        data={
            "left": [ie.left for ie in i],
            "right": [ie.right for ie in i],
            "label": [L0 for _ in i],
        }
    )
    for a in row["annotations"]:
        if TARGET in a["label"]:
            start_ms = a["start_ms"]
            end_ms = a["end_ms"]
            c = (df.right > start_ms) & (df.right < end_ms)
            df.loc[c, "label"] = L1
    # mapper = {L0: [1,0], L1: [0,1]}
    # return df.label.apply(lambda s: mapper.get(s)).tolist()
    return df.label.tolist()
    2 + 2


df["labels"] = df.apply(assign_labels, axis="columns")
df.to_csv(f"filtered_{TARGET}.csv", index=False)
df.tail(530)[["segment_path", "labels"]].to_csv(f"test_{TARGET}.csv", index=False)


if not REPEAT_TARGETS:
    df.head(530)[["segment_path", "labels"]].to_csv(f"train_{TARGET}.csv", index=False)
else:
    df = df.head(530)
    df["positive"] = df.labels.map(lambda l: TARGET in l)

    ps = df.index[df.positive].tolist()
    ns = df.index[~df.positive].tolist()
    indices = ns + ps * REPEAT_FACTOR
    df = (
        df.loc[indices, :]
        .sample(frac=1.0)
        .drop(columns=["positive"])
        .reset_index(drop=True)
    )
    df[["segment_path", "labels"]].to_csv(f"train_{TARGET}.csv", index=False)
2 + 2
