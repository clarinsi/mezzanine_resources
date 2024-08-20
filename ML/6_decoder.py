import pandas as pd
from pydub import AudioSegment
from itertools import pairwise

TARGET = "prosodicUnit"
f = f"model_{TARGET}_1e-5_30_1.csv"
df = pd.read_csv(f).head(5)
segments = []
for i, row in df.iterrows():
    i = 0
    AS = AudioSegment.from_file("data/" + row["segment_path"])
    row["y_pred"] = [str(i) for i in eval(row["y_pred"])]
    ndf = pd.DataFrame(
        data={
            "millisecond": [20 * i for i in range(len(row["y_pred"]))],
            "y_pred": row["y_pred"],
            "y_true": eval(row["labels"]),
        }
    )
    ndf["y_pred"] = ndf.y_pred.str.replace(TARGET, "1").astype(int)
    ndf["y_true"] = ndf.y_true.str.replace(TARGET, "1").astype(int)
    ndf["millisecond"] = ndf.millisecond.astype(int)
    ndf = ndf.dropna()
    indices_of_change = ndf.y_pred.diff()[ndf.y_pred.diff() != 0].index.values
    print("Predicted segments:")
    for si, ei in pairwise(indices_of_change):
        if ndf.loc[si:ei, "y_pred"].mode()[0] != 1:
            continue
        else:
            print(
                f"""{ndf.loc[si, "millisecond"]} : {ndf.loc[ei, "millisecond"]} (PU_pred_{i}.mp3)"""
            )
            current = AS[ndf.loc[si, "millisecond"] : ndf.loc[ei, "millisecond"]]
            segments.append(current)
            current.export(f"PU_pred_{i}.mp3", format="mp3")
            i += 1
    indices_of_change = ndf.y_true.diff()[ndf.y_true.diff() != 0].index.values
    print("True segments:")
    for si, ei in pairwise(indices_of_change):
        if ndf.loc[si:ei, "y_true"].mode()[0] != 1:
            continue
        else:
            print(
                f"""{ndf.loc[si, "millisecond"]} : {ndf.loc[ei, "millisecond"]} (PU_true_{i}.mp3)"""
            )
            current = AS[ndf.loc[si, "millisecond"] : ndf.loc[ei, "millisecond"]]
            segments.append(current)
            current.export(f"PU_true_{i}.mp3", format="mp3")
            i += 1

    AS.export("PU_whole.mp3")
if not segments:
    raise LookupError("No disfluencies for this type.")
out = sum(segments)
out.export(f.replace(".csv", ".mp3").replace("model", "output"), format="mp3")
