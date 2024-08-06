import numpy as np
import os
import pandas as pd
import datasets
from datasets import load_dataset, load_metric, Audio
from itertools import zip_longest

os.environ["WANDB_DISABLED"] = "true"
os.environ["CUDA_VISIBLE_DEVICES"] = "4"

MAX_LENGTH_s = 30
FRAME_RATE_hz = 50

def label_processor(l: list):
    mapper = {
        "0": [1,0],
        "filledPause": [0,1]
    }
    return [mapper[i] for i in eval(l)]

train = pd.read_csv("data/train.csv").rename(columns={"segment_path": "audio", "labels": "label"})
train["label"] = train.label.apply(label_processor)
train["audio"] = train.audio.apply(lambda s: "data/"+s)
train = datasets.Dataset.from_pandas(train).cast_column("audio", Audio(sampling_rate = 16_000, mono=True))

# test = pd.read_csv("data/test.csv").rename(columns={"segment_path": "audio", "labels": "label"})
# test["label"] = test.label.apply(label_processor)
# test["audio"] = test.audio.apply(lambda s: "data/"+s)
# test = datasets.Dataset.from_pandas(test).cast_column("audio", Audio(sampling_rate = 16_000, mono=True))



import transformers
from transformers import AutoFeatureExtractor, Wav2Vec2BertForAudioFrameClassification

import torch

device = torch.device("cuda")
from tqdm.auto import tqdm


feature_extractor = AutoFeatureExtractor.from_pretrained("facebook/w2v-bert-2.0",)
model = Wav2Vec2BertForAudioFrameClassification.from_pretrained(
    "facebook/w2v-bert-2.0", num_labels=2
).cuda()

model.freeze_base_model()

def preprocess_function(examples):
    sampling_rate = examples["audio"][0]["sampling_rate"]

    inputs = feature_extractor(
        [x["array"] for x in examples["audio"]],
        sampling_rate=sampling_rate,
        return_tensors= "pt"

    )
    inputs = feature_extractor.pad(
        inputs,
        padding="max_length",  # pad to max_length, not just to the longest sequence
        max_length=int(FRAME_RATE_hz * MAX_LENGTH_s),
        truncation=False,
        return_tensors= "pt"
    )
    # M = inputs["input_features"].shape[1]
    # inputs["label"] = [i[0:M] for i in examples["label"]]
    inputs["label"] = torch.Tensor([[x for x, y in zip_longest(e, i, fillvalue = [1,0])] for e, i in zip(examples["label"], inputs["input_features"])])
    inputs=inputs.convert_to_tensors(tensor_type="pt").to(device)
    return inputs


train = train.map(preprocess_function, batched=True, batch_size=50, remove_columns="audio", desc="Extracting features for train")
# test = test.map(preprocess_function, batched=True, batch_size=50, remove_columns="audio", desc="Extracting features for test")

from transformers import Trainer as Trainer, TrainingArguments as TrainingArguments

training_args = TrainingArguments(
    output_dir="brisi",
    overwrite_output_dir=True,
    learning_rate=2e-5,
    per_device_train_batch_size=16,
    num_train_epochs=10,
    weight_decay=0.01,
    save_strategy="epoch",
    logging_steps = 10,

)
trainer = Trainer(model=model, args=training_args, train_dataset=train,
tokenizer=feature_extractor
                  )
trainer.train()