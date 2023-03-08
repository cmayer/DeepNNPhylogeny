#!/usr/bin/env python3

import os
from configparser import ConfigParser

config = ConfigParser()

home_dir = os.path.expanduser("~")

nn_search_path = {
    0: os.path.join(home_dir, "PhylNNsaved", "DNA-ModelPred"),
    1: os.path.join(home_dir, "PhylNNsaved", "AA-ModelPred"),
    2: os.path.join(home_dir, "PhylNNsaved", "TopPredJC-normalisation-layer"),
    3: os.path.join(home_dir, "PhylNNsaved", "TopPredK2P-normalisation-layer"),
    4: os.path.join(home_dir, "PhylNNsaved", "TopPredF81-normalisation-layer"),
    5: os.path.join(home_dir, "PhylNNsaved", "TopPredF84-normalisation-layer"),
    6: os.path.join(home_dir, "PhylNNsaved", "TopPredHKY-normalisation-layer"),
    7: os.path.join(home_dir, "PhylNNsaved", "TopPredGTR-normalisation-layer"),
    8: os.path.join(home_dir, "PhylNNsaved", "TopPredJTT"),
    9: os.path.join(home_dir, "PhylNNsaved", "TopPredLG"),
    10: os.path.join(home_dir, "PhylNNsaved", "TopPredWAG"),
    11: os.path.join(home_dir, "PhylNNsaved", "TopPredDayhoff"),
    12: os.path.join(home_dir, "PhylNNsaved")
}

config["NN-Search-Path"] = nn_search_path

with open("DeepNNPhylogeny.config", "w") as f:
    config.write(f)
