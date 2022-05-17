import re

import pandas as pd


class VinaParser:
    """Parses the logs of vina"""

    def __init__(self, filename: str):
        self.filename = filename
        basename = filename.split("/")[-1].split(".")[0]
        self.prot, self.lig = basename.split("_")
        self.read_file()

    def read_file(self):
        """Reads the file and stores the text"""
        with open(self.filename, "r") as file:
            self.text = file.read()

    def get_runtime(self) -> float:
        """Returns the runtime of the vina run"""
        return float(re.search(r"in\s(\d+\.\d+)\sseconds\n", self.text).groups()[0])

    def get_top_score(self) -> float:
        """Returns the top score of the vina run"""
        return float(re.search(r"1\s+(\-*\d+\.\d+)\s+", self.text).groups()[0])

    def run(self) -> dict:
        """Returns a dictionary with the results of the vina run"""
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
