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

MAX_LENGTH_s = 30
FRAME_RATE_hz = 50

def label_processor(l: list):
    mapper = {
        "0": [1,0],
        "filledPause": [0,1]
    }
    return [mapper[i] for i in eval(l)]

# train = pd.read_csv("data/train.csv").rename(columns={"segment_path": "audio", "labels": "label"})
# train["label"] = train.label.apply(label_processor)
# train["audio"] = train.audio.apply(lambda s: "data/"+s)
# train = datasets.Dataset.from_pandas(train).cast_column("audio", Audio(sampling_rate = 16_000, mono=True))

test = pd.read_csv("data/test.csv").rename(columns={"segment_path": "audio", "labels": "label"})
test["label"] = test.label.apply(label_processor)
test["audio"] = test.audio.apply(lambda s: "data/"+s)
test = datasets.Dataset.from_pandas(test).cast_column("audio", Audio(sampling_rate = 16_000, mono=True))

feature_extractor = AutoFeatureExtractor.from_pretrained("./brisi/checkpoint-266")
model = Wav2Vec2BertForAudioFrameClassification.from_pretrained(
    "./brisi/checkpoint-266" ).to(device)
from sklearn.metrics import classification_report, confusion_matrix
y_trues = []
y_preds = []

def evaluate(chunks):
    sampling_rate = chunks["audio"][0]["sampling_rate"]
    with torch.no_grad():
        inputs = feature_extractor([i["array"] for i in chunks["audio"]],
                return_tensors="pt", sampling_rate = sampling_rate).to(device)
        logits = model(**inputs).logits
    y_pred = np.array(logits.cpu()).argmax(axis=-1)
    return {"y_pred": [i.tolist() for i in y_pred]}


test = test.map(evaluate, batched=True, batch_size=50)
y_preds = [[i for i in row] for row in test["y_pred"]]
test = pd.read_csv("data/test.csv")
test["y_pred"] = y_preds
test.to_csv("pred30epoch.csv")