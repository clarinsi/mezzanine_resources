import pandas as pd
from itertools import pairwise
from sklearn.metrics import classification_report, confusion_matrix
from itertools import chain
from pathlib import Path
from pydub import AudioSegment

f = "model_prosodicUnit_1e-5_30_1.csv"
o = f.replace("model", "instances").replace(".csv", "")
# csvs = [str(i) for i in Path(".").glob("model_*.csv")]
# for f in csvs:

# o = Path(f).replace("model", "stats").replace(".csv", ".txt")
df = pd.read_csv(f).head(1)
# df["labels"] = df.labels.apply(eval).apply(lambda l: [0 if i == "0" else 1 for i in l])
# df["y_pred"] = df.y_pred.apply(eval)


# Reverse labelling:
df["labels"] = df.labels.apply(eval).apply(lambda l: [0 if i == "0" else 1 for i in l])
df["y_pred"] = df.y_pred.apply(eval).apply(lambda l: [0 if i == 0 else 1 for i in l])


def frames_to_intervals(frames: list) -> list[pd.Interval]:
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


df["gold_intervals"] = df.labels.apply(frames_to_intervals)
df["pred_intervals"] = df.y_pred.apply(frames_to_intervals)


def extract(gold_intervals: list[pd.Interval], pred_intervals: list[pd.Interval]):
    TP, FN, FP = [], [], []
    pred_intervals = [i + 1 for i in pred_intervals]
    events = set(gold_intervals + pred_intervals)
    inhibited_events = []
    for event in events:
        if event in inhibited_events:
            continue
        others = [e for e in events if ((e not in inhibited_events) and (e != event))]
        overlapping = [other for other in others if events_overlap(event, other)]
        if not overlapping:
            if event in gold_intervals:
                FN.append(event)
                inhibited_events.append(event)
            else:
                FP.append(event)
                inhibited_events.append(event)
        else:
            other = overlapping[0]
            TP.append(
                pd.Interval(
                    min(event.left, other.left),
                    max(event.right, other.right),
                )
            )
            inhibited_events.append(event)
            inhibited_events.append(other)

    assert len(TP) * 2 + len(FN) + len(FP) == len(gold_intervals) + len(pred_intervals)
    return TP, FN, FP


def events_overlap(this, other):
    return (this.right >= other.left) and (this.left <= other.right)


df["TP"] = df.apply(
    lambda row: extract(row["gold_intervals"], row["pred_intervals"])[0], axis=1
)
df["FN"] = df.apply(
    lambda row: extract(row["gold_intervals"], row["pred_intervals"])[1], axis=1
)
df["FP"] = df.apply(
    lambda row: extract(row["gold_intervals"], row["pred_intervals"])[2], axis=1
)

TPs, FNs, FPs = [], [], []
for i, row in df.iterrows():
    AS = AudioSegment.from_file("data/" + row["segment_path"])
    for interval in row["TP"]:
        TPs.append(AS[interval.left : interval.right])
    for interval in row["FN"]:
        FNs.append(AS[interval.left : interval.right])
    for interval in row["FP"]:
        FPs.append(AS[interval.left : interval.right])
print(f"""{len(TPs)=}
{len(FPs)=}
{len(FNs)=}""")

if TPs:
    sum(TPs).export(f"{o}_TP.mp3", format="mp3")
if FNs:
    sum(FNs).export(f"{o}_FN.mp3", format="mp3")
if FPs:
    sum(FPs).export(f"{o}_FP.mp3", format="mp3")
2 + 2
