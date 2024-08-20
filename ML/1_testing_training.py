import numpy as np
import os

os.environ["WANDB_DISABLED"] = "true"
os.environ["CUDA_VISIBLE_DEVICES"] = "4"
import pandas as pd
import datasets
from datasets import load_dataset, load_metric, Audio

# Dummy wav
wavpath = "MP_01.wav"
frame_rate_hz = 50  # Hz
train_duration_s = 2
device = "cuda"
num_samples = 5
dataset = datasets.Dataset.from_pandas(
    pd.DataFrame(
        data={
            "audio": [wavpath for i in range(num_samples)],
        }
    )
)
dataset = dataset.cast_column("audio", Audio(mono=True, sampling_rate=16000))


def trimmer(example):
    example["audio"]["array"] = example["audio"]["array"][
        0 : int((train_duration_s-1)* example["audio"]["sampling_rate"])
    ]
    return example


dataset = dataset.map(trimmer)
dataset = dataset.add_column(
    "label",
    np.random.choice(
        [0, 1], (num_samples, int(train_duration_s * frame_rate_hz), 3)
    ).tolist(),
)


import transformers
from transformers import AutoFeatureExtractor, Wav2Vec2BertForAudioFrameClassification

import torch

torch.device(device)
from tqdm.auto import tqdm


feature_extractor = AutoFeatureExtractor.from_pretrained("facebook/w2v-bert-2.0")
model = Wav2Vec2BertForAudioFrameClassification.from_pretrained(
    "facebook/w2v-bert-2.0", num_labels=3
)

model.freeze_base_model()


def preprocess_function(examples):
    sampling_rate = examples["audio"][0]["sampling_rate"]

    inputs = feature_extractor(
        [x["array"] for x in examples["audio"]],
        sampling_rate=sampling_rate,
        # padding="max_length",  # pad to max_length, not just to the longest sequence
        # max_length=int(frame_rate_hz * train_duration_s),
        # truncation=False,
        return_tensors= "pt"

    )
    inputs = feature_extractor.pad(
        inputs,
        padding="max_length",  # pad to max_length, not just to the longest sequence
        max_length=int(frame_rate_hz * train_duration_s),
        truncation=False,
        return_tensors= "pt"
    )
    M = inputs["input_features"].shape[1]
    inputs["label"] = [i[0:M] for i in examples["label"]]
    return inputs


dataset = dataset.map(preprocess_function, batched=True, remove_columns="audio")

from transformers import Trainer as Trainer, TrainingArguments as TrainingArguments

training_args = TrainingArguments(
    output_dir="brisi",
    learning_rate=5e-5,
    per_device_train_batch_size=16,
    num_train_epochs=10,
    weight_decay=0.01,
    save_strategy="epoch",
)
trainer = Trainer(model=model, args=training_args, train_dataset=dataset,
    tokenizer=feature_extractor,
                  )
trainer.train()

model.save_pretrained("./brisi")