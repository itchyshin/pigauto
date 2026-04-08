"""
Standalone TabPFN imputation for use from R via system2().

Usage:
    python run_tabpfn.py <input.csv> <output.csv>

Reads a numeric matrix (missing values as empty cells) from <input.csv>,
runs tabimpute.ImputePFN, and writes the imputed matrix to <output.csv>.
CSV files have no row names; the first row is a header (trait names).
"""

import sys
import warnings
import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

from tabimpute.interface import ImputePFN


def main():
    if len(sys.argv) != 3:
        print("Usage: run_tabpfn.py <input.csv> <output.csv>", file=sys.stderr)
        sys.exit(2)

    infile, outfile = sys.argv[1], sys.argv[2]

    X = pd.read_csv(infile)
    cols = list(X.columns)
    X_np = X.values.astype(float)
    # Ensure NaN representation (pandas already uses NaN for empty cells)
    imp = ImputePFN(device="cpu")
    out_np = imp.impute(X_np)
    out_df = pd.DataFrame(out_np, columns=cols)
    out_df.to_csv(outfile, index=False)
    print(f"Wrote {out_df.shape[0]} x {out_df.shape[1]} matrix to {outfile}")


if __name__ == "__main__":
    main()
