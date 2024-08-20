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
os.environ["CUDA_VISIBLE_DEVICES"] = "0"
from pathlib import Path

TARGET = "lengthening"
checkpoint = Path(f"./model_{TARGET}_2e-5_20_1/checkpoint-2660")
TESTFILE = f"data/test_{TARGET}.csv"


def label_processor(l: list):
    mapper = {"0": [1, 0], TARGET: [0, 1]}
    return [mapper[i] for i in eval(l)]


test = pd.read_csv(TESTFILE).rename(
    columns={"segment_path": "audio", "labels": "label"}
)
test["label"] = test.label.apply(label_processor)
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
    return {"y_pred": [i.tolist()[: len(l)] for i, l in zip(y_pred, chunks["label"])]}


n_test = test.map(evaluator, batch_size=53, batched=True)
y_preds = [[i for i in row] for row in n_test["y_pred"]]
n_test = pd.read_csv(TESTFILE)
n_test["y_pred"] = y_preds
n_test.to_csv(f"{checkpoint.parent}.csv")
