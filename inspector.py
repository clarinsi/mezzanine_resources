from pathlib import Path
from conllu import parse
from lxml import etree
import pandas as pd
from numpy import median
from loguru import logger

connllufile = (
    "/home/peter/mezzanine_resources/ROG/CONLLU/Rog-Art-J-Gvecg-P580047.conllu"
)
target = "Artur-N-G6100-P610002.s25-s28_reseg.1661"

data = parse(Path(conllufile).read_text())
data = [i for i in data if i.metadata["sent_id"] == target][0]
for i in data:
    print(
        i,
        "::::",
        # i["deprel"],
        # i["misc"]["Gos2.1_token_id"],
        i.__repr__(),
        end="\n\n",
    )
