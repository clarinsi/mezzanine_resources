from pathlib import Path
import lxml.etree as ET
import pandas as pd

dev = list(Path("../../UD_Slovenian-SST").glob("*dev*.conllu"))[0]
train = list(Path("../../UD_Slovenian-SST").glob("*train*.conllu"))[0]
test = list(Path("../../UD_Slovenian-SST").glob("*test*.conllu"))[0]

outdir = "../../SST-split"
Path(outdir).mkdir(exist_ok=True, parents=True)
for file in Path(outdir).glob("*conllu"):
    file.unlink()

import conllu

dev_sentences = conllu.parse(dev.read_text())
for i in range(len(dev_sentences)):
    dev_sentences[i].metadata["split"] = "dev"

test_sentences = conllu.parse(test.read_text())
for i in range(len(test_sentences)):
    test_sentences[i].metadata["split"] = "test"

train_sentences = conllu.parse(train.read_text())
for i in range(len(train_sentences)):
    train_sentences[i].metadata["split"] = "train"

# This shall hold all the TokenLists with the added metadata on split.
combined = list(dev_sentences) + list(test_sentences) + list(train_sentences)
df = pd.DataFrame(
    {
        "tokenlist": combined,
        "first_token_id": [i[0]["misc"]["Gos2.1_token_id"] for i in combined],
        "sent_id": [i.metadata["sent_id"] for i in combined],
        "split": [i.metadata["split"] for i in combined],
        "sent_id": [i.metadata["sent_id"] for i in combined],
    }
)
### # Test for leakage between splits:
# print(
#     "Train intersection dev:",
#     set(df[df.split == "train"].first_token_id.tolist())
#     & set(df[df.split == "dev"].first_token_id.tolist()),
# )
# print(
#     "Train intersection test:",
#     set(df[df.split == "train"].first_token_id.tolist())
#     & set(df[df.split == "test"].first_token_id.tolist()),
# )
# print(
#     "dev intersection test:",
#     set(df[df.split == "dev"].first_token_id.tolist())
#     & set(df[df.split == "test"].first_token_id.tolist()),
# )


df["file"] = df.first_token_id.str.split(".").str[0].apply(lambda s: s + ".conllu")
pattern = r"(\d+)"
# Extract all matches and convert to a list of tuples
df["extracted_digits"] = (
    df["sent_id"]
    .str.split(".")
    .str[-1]
    .str.extractall(pattern)
    .groupby(level=0)[0]
    .apply(lambda x: tuple(x.astype(int)))
)

df = df.sort_values(by="extracted_digits").reset_index(drop=True)

from tqdm import tqdm

for i, row in tqdm(df.iterrows(), total=df.shape[0]):
    with open(str(Path(outdir) / row["file"]), "a") as f:
        f.write(row["tokenlist"].serialize() + "\n")
2 + 2
