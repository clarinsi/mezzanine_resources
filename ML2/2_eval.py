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
os.environ["CUDA_VISIBLE_DEVICES"] = "1"
from pathlib import Path

TARGET = "filledPause"
checkpoints = Path(f"./model_{TARGET}_3e-5_20_4/").glob("checkpoint-*")
TESTFILE = f"test_{TARGET}.jsonl"

for checkpoint in checkpoints:
    if checkpoint.name.endswith(".jsonl"):
        continue

    def label_processor(l: list):
        mapper = {0: [1, 0], 1: [0, 1]}
        return [mapper[i] for i in l]

    test = pd.read_json(TESTFILE, lines=True).rename(columns={"segment_path": "audio"})
    test["label"] = test.y_true.apply(label_processor)
    test["audio"] = test.audio.apply(lambda s: "data/" + s)
    test = datasets.Dataset.from_pandas(test).cast_column(
        "audio", Audio(sampling_rate=16_000, mono=True)
    )

    feature_extractor = AutoFeatureExtractor.from_pretrained(str(checkpoint))
    model = Wav2Vec2BertForAudioFrameClassification.from_pretrained(str(checkpoint)).to(
        device
    )
    from sklearn.metrics import classification_report, confusion_matrix

    y_trues = []
    y_preds = []

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

    n_test = test.map(evaluator, batch_size=15, batched=True, desc="Running inference")
    y_preds = [[i for i in row] for row in n_test["y_pred"]]
    n_test = pd.read_json(TESTFILE, lines=True)
    n_test["y_pred"] = y_preds
    checkpoint.mkdir(parents=True, exist_ok=True)
    n_test.to_json(
        checkpoint.with_name(checkpoint.name + "_predictions.jsonl"),
        orient="records",
        lines=True,
    )
