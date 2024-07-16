from pathlib import Path
import EXBUtils
import conllu

ddir = Path("../../Iriss-disfl-anno-phase5-fin-corr")
darinkas = list(ddir.glob("*.exb.xml"))
kdir = Path("../../UD-SST-split")
for d in darinkas:
    print("\n\nChecking file ", d.name)
    doc = EXBUtils.EXB(d)
    df = doc.pandalize()
    pattern = d.name.replace("Iriss", "Artur").replace(".exb.xml", "")
    k = list(kdir.glob(f"*{pattern}*"))[0]

    cnl = conllu.parse(Path(k).read_text())
    conllu_words = []
    for tokenlist in cnl:
        for token in tokenlist:
            if (token["upos"] == "PUNCT") or token["form"].startswith("[name:"):
                continue
            conllu_words.append(token["form"])

    for c, (i, row) in zip(conllu_words, df.iterrows()):
        if c != row["norm"].replace("…", "-"):
            print(f"Conllu: {c:^15}, exb: {row["norm"]:^15}")
            # print(row)
        else:
            print("\rMatched\r", end="")
