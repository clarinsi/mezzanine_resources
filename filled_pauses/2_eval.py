from transformers import AutoFeatureExtractor, Wav2Vec2BertForAudioFrameClassification
from datasets import load_dataset, Dataset, Audio
import torch
import numpy as np
import soundfile as sf
import tqdm
import os

device = torch.device("cuda")
import numpy as np
import os
import pandas as pd
import datasets
from datasets import load_dataset, load_metric, Audio
from itertools import zip_longest

os.environ["WANDB_DISABLED"] = "true"
os.environ["CUDA_VISIBLE_DEVICES"] = "2"
from pathlib import Path

TARGET = "filledPause"
checkpoints = Path(f"./model_{TARGET}_3e-5_20_4/").glob("checkpoint-*")
TESTFILE = f"test_{TARGET}.jsonl"

for checkpoint in checkpoints:
    if checkpoint.name.endswith(".jsonl"):
        continue

    def label_processor(l: list):
        mapper = {"0": [1, 0], TARGET: [0, 1]}
        return [mapper[i] for i in l]

    def inverse_label_processor(l: list):
        return [1 if i == [0, 1] else 0 for i in l]

    df = pd.read_csv("data/filled_pauses.csv").drop(
        columns="wavpath duration annotations fileroot".split()
    )
    df["wavpath"] = df.apply(
        lambda row: "data/test/"
        + row["file"].replace(".exb.xml", "")
        + f"-avd_{row['start_ms']:0.3f}_{row['end_ms']:0.3f}.wav",
        axis=1,
    )
    wavs_to_test_on = [str(i) for i in Path("data/test").glob("*.wav")]

    df["label"] = df.labels.apply(eval).apply(label_processor)
    print("Before selecting test split:", df.shape)
    df = df[df.wavpath.isin(wavs_to_test_on)].reset_index(drop=True)
    print("After selecting test split:", df.shape)

    df["audio"] = df.wavpath.values
    df = df.drop(
        columns=["file", "start_ms", "end_ms", "segment_path", "labels", "wavpath"]
    )
    y_true = [inverse_label_processor(i) for i in df.label.values]
    ds = datasets.Dataset.from_pandas(df).cast_column(
        "audio",
        Audio(sampling_rate=16_000, mono=True),
    )
    feature_extractor = AutoFeatureExtractor.from_pretrained(str(checkpoint))
    model = Wav2Vec2BertForAudioFrameClassification.from_pretrained(str(checkpoint)).to(
        device
    )
    from sklearn.metrics import classification_report, confusion_matrix

    def evaluator(chunks):
        sampling_rate = chunks["audio"][0]["sampling_rate"]
        with torch.no_grad():
            inputs = feature_extractor(
                [i["array"] for i in chunks["audio"]],
                return_tensors="pt",
                sampling_rate=sampling_rate,
            ).to(device)
            logits = model(**inputs).logits
        y_pred = np.array(logits.cpu()).argmax(axis=-1)
        return {
            "y_pred": [i.tolist()[: len(l)] for i, l in zip(y_pred, chunks["label"])]
        }

    n_test = ds.map(evaluator, batch_size=15, batched=True, desc="Running inference")
    y_preds = [[i for i in row] for row in n_test["y_pred"]]
    df["y_true"] = y_true
    df["y_pred"] = y_preds
    checkpoint.mkdir(parents=True, exist_ok=True)
    df = df.rename(columns={"audio": "segment_path"}).drop(columns="label")
    df.to_json(
        checkpoint.with_name(checkpoint.name + "_predictions.jsonl"),
        orient="records",
        lines=True,
    )
