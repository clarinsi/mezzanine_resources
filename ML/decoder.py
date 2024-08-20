import pandas as pd
from pydub import AudioSegment
from itertools import pairwise

df = pd.read_csv("model_5e-5_20.csv")
c = 0
for i, row in df.iterrows():
    # if not "N-G6060-P606001-avd_53.264_83.061.wav" in row["segment_path"]: continue
    AS = AudioSegment.from_file("data/" + row["segment_path"])
    ndf = pd.DataFrame(
        data={
            "millisecond": [20 * i for i in range(len(eval(row["y_pred"])))],
            "label": eval(row["y_pred"]),
        }
    )
    # ndf["label"] = ndf.label.str.replace("filledPause","1").astype(int)
    ndf["millisecond"] = ndf.millisecond.astype(int)
    ndf = ndf.dropna()
    indices_of_change = ndf.label.diff()[ndf.label.diff() != 0].index.values
    for si, ei in pairwise(indices_of_change):
        if ndf.loc[si:ei, "label"].mode()[0] != 1:
            continue
        else:
            AS[ndf.loc[si, "millisecond"] : ndf.loc[ei, "millisecond"]].export(
                f"{c}.wav",
                format="wav",
                parameters=[
                    "-ac",
                    "1",
                    "-ar",
                    "16000",
                ],
            )
            c += 1
            2 + 2
