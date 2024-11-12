from transformers import AutoFeatureExtractor, Wav2Vec2BertForAudioFrameClassification
from datasets import Dataset, Audio
import torch
import numpy as np
from pathlib import Path

device = torch.device("cuda")
model_name = "5roop/wav2vecbert2-filledPause"
feature_extractor = AutoFeatureExtractor.from_pretrained(model_name)
model = Wav2Vec2BertForAudioFrameClassification.from_pretrained(model_name).to(device)

ds = Dataset.from_dict(
    {
        "audio": [
            "/cache/peterr/mezzanine_resources/filled_pauses/data/dev/Iriss-J-Gvecg-P500001-avd_2082.293_2112.194.wav"
        ],
    }
).cast_column("audio", Audio(sampling_rate=16_000, mono=True))


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
    return {"y_pred": y_pred.tolist()}


ds = ds.map(evaluator, batched=True)
print(ds["y_pred"][0])
