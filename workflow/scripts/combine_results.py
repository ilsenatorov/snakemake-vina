import os
import re
import numpy as np
import pandas as pd


class VinaParser:
    def __init__(self, filename):
        self.filename = filename
        basename = filename.split("/")[-1].split(".")[0]
        self.prot, self.lig = basename.split("_")
        self.read_file()

    def read_file(self):
        with open(self.filename, "r") as file:
            self.text = file.read()

    def get_runtime(self):
        return float(re.search("in\s(\d+\.\d+)\sseconds\n", self.text).groups()[0])

    def get_top_score(self):
        return float(re.search("1\s+(\-*\d+\.\d+)\s+", self.text).groups()[0])

    def run(self):
        return dict(
            prot=self.prot,
            lig=self.lig,
            score=self.get_top_score(),
            runtime=self.get_runtime(),
        )


if __name__ == "__main__":
    df = []
    for i in snakemake.input:
        df.append(VinaParser(i).run())
    df = pd.DataFrame(df, columns=["prot", "lig", "score", "runtime"])
    df = df[df["prot"] != 0]
    df["pair"] = df["prot"] + "_" + df["lig"]
    df.to_csv(snakemake.output[0], index=False)
