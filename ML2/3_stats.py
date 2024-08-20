import pandas as pd
from pathlib import Path
from sklearn.metrics import classification_report
import numpy as np
from itertools import pairwise


def elementwise_metrics(y_true, y_pred):
    M = min(len(y_true), len(y_pred))
    y_true = np.array(y_true)[:M]
    y_pred = np.array(y_pred)[:M]
    return classification_report(
        y_true, y_pred, output_dict=True, zero_division=np.nan
    ).get(
        "1",
        {"precision": np.nan, "recall": np.nan, "f1-score": np.nan, "support": 0},
    )


def overlap_metrics(y_true, y_pred):
    M = min(len(y_true), len(y_pred))
    y_true = np.array(y_true)[:M]
    y_pred = np.array(y_pred)[:M]
    true_intervals = frames_to_intervals(y_true)
    pred_intervals = frames_to_intervals(y_pred)
    TP, FN, FP = extract(true_intervals, pred_intervals)
    TP, FN, FP = len(TP), len(FN), len(FP)
    support = len(true_intervals)
    if (support == 0) or (TP == 0):
        return {
            "precision": np.nan,
            "recall": np.nan,
            "f1-score": np.nan,
            "support": support,
        }
    precision = TP / (TP + FP)
    recall = TP / (TP + FN)
    F_1 = 2 * precision * recall / (precision + recall)
    return {
        "precision": precision,
        "recall": recall,
        "f1-score": F_1,
        "support": support,
    }


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


def events_overlap(this: pd.Interval, other: pd.Interval):
    return (this.right >= other.left) and (this.left <= other.right)


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


f = Path("model_filledPause_3e-5_20_4/checkpoint-211_predictions.jsonl")
for f in Path("model_filledPause_3e-5_20_4/").glob("checkpoint-*_predictions.jsonl"):
    df = pd.read_json(f, lines=True)


    df["elementwise"] = df.apply(
        lambda row: elementwise_metrics(row["y_true"], row["y_pred"]), axis="columns"
    )
    df["overlap"] = df.apply(
        lambda row: overlap_metrics(row["y_true"], row["y_pred"]), axis="columns"
    )
    dfe = pd.json_normalize(df.elementwise)
    dfo = pd.json_normalize(df.overlap)

    dfboth = pd.merge(
        dfe,
        dfo,
        left_index=True,
        right_index=True,
        suffixes = ("_elementwise", "_overlap")
    )
    dfboth["segment_path"] = df.segment_path.values
    dfboth.to_json(str(f).replace("_predictions", "_stats"), lines=True, orient="records")
2 + 2
