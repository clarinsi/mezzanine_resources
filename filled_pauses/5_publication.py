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
model_path = "./model_filledPause_3e-5_20_4/checkpoint-900"
feature_extractor = AutoFeatureExtractor.from_pretrained(model_path)
model = Wav2Vec2BertForAudioFrameClassification.from_pretrained(model_path)
model.push_to_hub("5roop/wav2vecbert2-filledPause")
