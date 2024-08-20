from transformers import AutoFeatureExtractor, Wav2Vec2BertForAudioFrameClassification
from datasets import load_dataset, Dataset, Audio
import torch
import numpy as np
import pandas as pd
import soundfile as sf
from scipy.io import wavfile
import tqdm
import os

device = torch.device("cuda")
import numpy as np
import os
import pandas as pd
import datasets
from datasets import load_dataset, load_metric, Audio
from itertools import zip_longest, pairwise

os.environ["WANDB_DISABLED"] = "true"
os.environ["CUDA_VISIBLE_DEVICES"] = "0"
checkpoint = "model_filledPause_2e-5_20_4/checkpoint-660"


feature_extractor = AutoFeatureExtractor.from_pretrained(str(checkpoint))
model = Wav2Vec2BertForAudioFrameClassification.from_pretrained(str(checkpoint)).to(
    device
)


def evaluator(chunk):
    sampling_rate = chunk["audio"]["sampling_rate"]
    with torch.no_grad():
        inputs = feature_extractor(
            chunk["audio"]["array"],
            return_tensors="pt",
            sampling_rate=sampling_rate,
        ).to(device)
        logits = model(**inputs).logits
    y_pred = np.array(logits.cpu()).argmax(axis=-1).reshape(-1)
    return {"y_pred": y_pred}


# # Fake dataset:
# ds = Dataset.from_dict(
#     {
#         "audio": [
#             "data/segments30s/Iriss-J-Gvecg-P500001-avd_2082.293_2112.194.wav",
#             "data/segments30s/Iriss-J-Gvecg-P500002-avd_341.769_371.477.wav",
#         ]
#     }
# ).cast_column("audio", Audio(mono=True, sampling_rate=16000))
# ds = ds.to_iterable_dataset()

# Real dataset
print("Start loading in streaming mode")
ds = load_dataset("classla/ParlaSpeech-HR", split="train", streaming=True)
print("Loading done.")


c = 0
for datum in ds:
    y_pred = evaluator(datum)["y_pred"]
    ndf = pd.DataFrame(
        data={
            "millisecond": [20 * i for i in range(len(y_pred))],
            "y_pred": y_pred,
        }
    )
    ndf["y_pred"] = ndf.y_pred.astype(int)
    ndf["millisecond"] = ndf.millisecond.astype(int)
    ndf = ndf.dropna()
    indices_of_change = ndf.y_pred.diff()[ndf.y_pred.diff() != 0].index.values
    for si, ei in pairwise(indices_of_change):
        if ndf.loc[si:ei, "y_pred"].mode()[0] != 1:
            continue
        else:
            # print("will save from ", si, "to", ei)
            tosave = datum["audio"]["array"][
                int(si / 50 * datum["audio"]["sampling_rate"]) : int(
                    ei / 50 * datum["audio"]["sampling_rate"]
                )
            ]
            wavfile.write(
                f"ParlaSpeech-HR_{c}.wav",
                datum["audio"]["sampling_rate"],
                np.int16(tosave / np.max(np.abs(tosave)) * 32767),
            )
            c += 1
