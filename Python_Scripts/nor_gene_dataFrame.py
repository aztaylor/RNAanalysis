import pandas as pd
import numpy as np

genes = {
    "SOS responsive, LexA-dependent genes":
    np.array(["polB", "recA", "recN", "sbmC", "ssb", "sulA", "uvrA", "uvrB"]),

    "LexA-independent genes":
    np.array(["dnaA", "nrdA", "nrdB", "cydA", "tdcB", "tdcC", "cstA", "sspA", "sspB",
    "deoA", "deoB", "deoC", "dnaG", "gmk", "gyrA", "gyrB", "intZ", "mcrB", "oraA", "pyrG"]),

    "Anaerobic induction":
    np.array(["cydA", "tdcB", "tdcC"]),

    "Starvation induction":
    np.array(["cstA", "sspA", "sspB"]),

    "Other regulation":
    np.array(["deoA", "deoB", "deoC", "dnaG", "gmk", "gyrA", "gyrB", "intZ", "mcrB", "oraA", "pyrG"]),

    "Transcription-related genes":
    np.array(["fliA", "fliZ", "greA", "rho", "rpoA", "rpoD"]),

    "Translation-related genes":
    np.array(["fusA", "infA", "rnpA", "rpmB", "rpmG", "rpmH", "rpsG", "rpsU", "trmD"]),

    "Hypothetical or unknown function genes":
    np.array(["yabO", "yagP", "ychB", "ydiY", "yebE", "yebF", "yeeN", "yfgB", "yfhN",
    "yfhO", "yghB", "yhjG", "yi52–10", "yibB", "yjeS", "ylaC", "yleA"])
}

gene_df = pd.DataFrame(dict([( k,pd.Series(v)) for k, v in genes.items() ]))
gene_df.to_csv("nor_genes.csv")
