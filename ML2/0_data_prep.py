from pathlib import Path
import pandas as pd
import json

segmentpath = Path("data/filled_pauses.jsonl")
splitspath = Path("../IRISS-train-dev-test-splits.csv")
data = [json.loads(i) for i in segmentpath.read_text().splitlines()]
TARGET = "filledPause"
new_data = []
for line in data:
    file = line["file"]
    for content in line["strips"]:
        new_data.append({"file": file, **content})
segdf = pd.DataFrame(new_data)
segdf["matcher"] = segdf.file.str.split(".").str[0].str.replace("Iriss", "Artur")
splitdf = pd.read_csv(splitspath, sep=";", nrows=57, encoding="cp1250").rename(
    columns={"Recording ID": "ID", "Train/Dev/Eval": "split"}
)
splitdf["split"] = splitdf.split.fillna("Train")
splitdf["subcorpus"] = splitdf.ID.str.split("-").str[1]
segdf["fileroot"] = segdf.file.apply(
    lambda s: str(Path(s).with_suffix("").with_suffix(""))
)
segdf["wavpath"] = [f"iriss_wavs/{i}-avd.wav" for i in segdf.fileroot]
segdf["segment_path"] = segdf.apply(
    lambda row: f"segments30s/{Path(row['wavpath']).with_suffix('').name}_{row['start_ms']:.3f}_{row['end_ms']:.3f}.wav",
    axis="columns",
)
segdf = segdf.merge(
    splitdf[["ID", "split", "subcorpus"]], left_on="matcher", right_on="ID", how="left"
).drop(columns=["ID", "matcher"])


def assign_labels(row):
    labels = []
    L0, L1 = 0, 1
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


segdf["y_true"] = segdf.apply(assign_labels, axis=1)


segdf.loc[segdf.split.isin(["Eval"]), ["segment_path", "y_true"]].to_json(
    f"test_{TARGET}.jsonl", lines=True, orient="records"
)
segdf.loc[segdf.split.isin(["Dev", "Train"]), ["segment_path", "y_true"]].to_json(
    f"traindev_{TARGET}.jsonl", lines=True, orient="records"
)
2 + 2
