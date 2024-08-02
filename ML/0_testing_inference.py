from transformers import AutoFeatureExtractor, Wav2Vec2BertForAudioFrameClassification
from datasets import load_dataset, Dataset, Audio
import torch
import numpy as np
import soundfile as sf

torch.device("cpu")

# from subprocess import run
# run("wget https://huggingface.co/classla/wav2vec2-large-slavic-parlaspeech-hr-lm/raw/main/00020578b.flac.wav", shell=True)

ds = load_dataset("classla/ParlaSpeech-RS", split="train", streaming=True)


for row in ds:
    sampling_rate = row["audio"]["sampling_rate"]

    feature_extractor = AutoFeatureExtractor.from_pretrained("facebook/w2v-bert-2.0")
    model = Wav2Vec2BertForAudioFrameClassification.from_pretrained(
        "facebook/w2v-bert-2.0"
    )
    for i in range(1, 300, 15):
        audio = row["audio"]["array"]
        for _ in range(i - 1):
            audio = np.concatenate([audio, row["audio"]["array"]])
        print(f"Duration: {audio.shape[0]/sampling_rate:0.1f}")
        # audio file is decoded on the fly
        inputs = feature_extractor(
            audio, return_tensors="pt", sampling_rate=sampling_rate
        )

        with torch.no_grad():
            logits = model(**inputs).logits

        probabilities = torch.sigmoid(logits[0])

        # labels is a one-hot array of shape (num_frames, num_speakers)

        labels = (probabilities > 0.5).long()
        print("\t", "Len labels: ", len(labels))
    # print(
    #     sampling_rate / (row["audio"]["array"].shape[0] / len(labels)),
    #     sampling_rate / (row["audio"]["array"].shape[0] / (len(labels) + 1)),
    # )
    2 + 2
