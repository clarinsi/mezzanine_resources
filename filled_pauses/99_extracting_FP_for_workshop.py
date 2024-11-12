import pandas as pd

df = pd.read_json(
    "/cache/peterr/mezzanine_resources/filled_pauses/model_filledPause_3e-5_20_4/checkpoint-900_predictions.jsonl",
    lines=True,
)
df = df[df.segment_path.str.contains("Iriss-J-Gvecg-P500014")]


def frames_to_intervals(frames: list) -> list[pd.Interval]:
    from itertools import pairwise

    return_list = []
    ndf = pd.DataFrame(
        data={
            "millisecond": [20 * i for i in range(len(frames))],
            "frames": frames,
        }
    )

    ndf["millisecond"] = ndf.millisecond.astype(int)
    ndf = ndf.dropna()
    indices_of_change = ndf.frames.diff()[ndf.frames.diff() != 0].index.values
    for si, ei in pairwise(indices_of_change):
        if ndf.loc[si : ei - 1, "frames"].mode()[0] == 0:
            pass
        else:
            return_list.append(
                pd.Interval(ndf.loc[si, "millisecond"], ndf.loc[ei - 1, "millisecond"])
            )
    return return_list


df["y_pred"] = df.y_pred.apply(frames_to_intervals)
df["y_true"] = df.y_true.apply(frames_to_intervals)
df["start"] = df.segment_path.str.split("_").str[1].astype(float)

intervals = []
for i, row in df.iterrows():
    start = row["start"]
    print(start)
    for i in row["y_pred"]:
        intervals.append(
            (round(start + i.left / 1000, 2), round(start + i.right / 1000, 2))
        )
# for i, inter in enumerate(intervals):
#     print(i, inter)
timestamps = sorted([i for j in intervals for i in j])
indices = [i for i in range(len(timestamps))]
for i, t in enumerate(timestamps):
    print(f"""<tli id="{i}" time="{t}"/>""")

for inter in intervals:
    left, right = inter
    left_index, right_index = None, None
    for i, t in enumerate(timestamps):
        if left == t:
            left_index = i
        if right == t:
            right_index = i
    print(f"""<event start="{left_index}" end="{right_index}">predicted FP</event>""")
2 + 2
