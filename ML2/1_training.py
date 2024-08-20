import numpy as np
import os
import pandas as pd
import datasets
from datasets import load_dataset, load_metric, Audio
from itertools import zip_longest

os.environ["WANDB_DISABLED"] = "true"
# os.environ["CUDA_VISIBLE_DEVICES"] = "1"

MAX_LENGTH_s = 30
FRAME_RATE_hz = 50
TARGET = "filledPause"


def label_processor(l: list):
    mapper = {0: [1, 0], 1: [0, 1]}
    return [mapper[i] for i in l]


train = pd.read_json(f"traindev_{TARGET}.jsonl", lines=True).rename(
    columns={"segment_path": "audio"}
)
train["label"] = train.y_true.apply(label_processor).drop(columns=["y_true"])
train["audio"] = train.audio.apply(lambda s: "data/" + s)
train = datasets.Dataset.from_pandas(train).cast_column(
    "audio", Audio(sampling_rate=16_000, mono=True), 
)

import transformers
from transformers import AutoFeatureExtractor, Wav2Vec2BertForAudioFrameClassification

import torch

device = torch.device("cuda")
from tqdm.auto import tqdm
import gc

torch.cuda.empty_cache()
gc.collect()

feature_extractor = AutoFeatureExtractor.from_pretrained(
    "facebook/w2v-bert-2.0",
)


def preprocess_function(examples):
    sampling_rate = examples["audio"][0]["sampling_rate"]

    inputs = feature_extractor(
        [x["array"] for x in examples["audio"]],
        sampling_rate=sampling_rate,
        return_tensors="pt",
    )
    inputs = feature_extractor.pad(
        inputs,
        padding="max_length",  # pad to max_length, not just to the longest sequence
        max_length=int(FRAME_RATE_hz * MAX_LENGTH_s),
        truncation=False,
        return_tensors="pt",
    )
    # M = inputs["input_features"].shape[1]
    # inputs["label"] = [i[0:M] for i in examples["label"]]
    inputs["label"] = torch.Tensor(
        [
            [x for x, y in zip_longest(e, i, fillvalue=[1, 0])]
            for e, i in zip(examples["label"], inputs["input_features"])
        ]
    )
    inputs = inputs.convert_to_tensors(tensor_type="pt").to(device)
    return inputs


train = train.map(
    preprocess_function,
    batched=True,
    batch_size=70,
    remove_columns="audio",
    desc="Extracting features for train",
)
from transformers import Trainer as Trainer, TrainingArguments as TrainingArguments
from itertools import product

LRs = [
    "3e-5",
    # "1e-5",
    # "8e-6",
]
EPs = [
    "20",
    # "10",
]
gases = [
    # 1,
    4
]
for LR, EP, gas in product(LRs, EPs, gases):
    model = Wav2Vec2BertForAudioFrameClassification.from_pretrained(
        "facebook/w2v-bert-2.0", num_labels=2
    ).cuda()
    print(
        f"Model will be saved to model_{TARGET}_{LR}_{EP}_{gas}",
    )
    D = 4
    training_args = TrainingArguments(
        output_dir=f"model_{TARGET}_{LR}_{EP}_{gas}",
        overwrite_output_dir=True,
        learning_rate=float(LR),
        per_device_train_batch_size=4 // D,
        gradient_accumulation_steps=gas * D,
        num_train_epochs=int(EP),
        weight_decay=0.01,
        save_strategy="epoch",
        logging_steps=10,
        # save_total_limit=2,
    )
    trainer = Trainer(
        model=model,
        args=training_args,
        train_dataset=train,
        tokenizer=feature_extractor,
    )
    trainer.train()
