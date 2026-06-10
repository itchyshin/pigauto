#!/usr/bin/env python3
"""Small TabPFN regression runner for pigauto benchmark scripts.

This intentionally lives in script/ rather than the package API. TabPFN is an
optional external benchmark dependency, not a pigauto runtime dependency.
"""

from __future__ import annotations

import argparse
import inspect
import json
import sys
import time

import numpy as np
import pandas as pd


def _make_regressor(device: str, seed: int):
    try:
        from tabpfn import TabPFNRegressor
    except Exception as exc:  # pragma: no cover - exercised on user machines
        raise RuntimeError(
            "Could not import tabpfn. Install it in the Python environment used "
            "by PIGAUTO_TABPFN_PYTHON, for example `python -m pip install "
            "tabpfn`, then rerun the benchmark."
        ) from exc

    params = inspect.signature(TabPFNRegressor).parameters
    kwargs = {}
    if "device" in params:
        kwargs["device"] = device
    if "random_state" in params:
        kwargs["random_state"] = seed
    if "ignore_pretraining_limits" in params:
        kwargs["ignore_pretraining_limits"] = True

    return TabPFNRegressor(**kwargs)


def _read_train(path: str, target: str):
    frame = pd.read_csv(path)
    if target not in frame.columns:
        raise ValueError(f"Target column {target!r} is missing from {path}")
    y = frame.pop(target).to_numpy(dtype=np.float32)
    x = frame.to_numpy(dtype=np.float32)
    return x, y, list(frame.columns)


def _read_predict(path: str):
    frame = pd.read_csv(path)
    return frame.to_numpy(dtype=np.float32), list(frame.columns)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--train", required=True)
    parser.add_argument("--predict", required=True)
    parser.add_argument("--target", default=".target")
    parser.add_argument("--out", required=True)
    parser.add_argument("--metadata", default=None)
    parser.add_argument("--device", default="auto")
    parser.add_argument("--seed", type=int, default=1)
    args = parser.parse_args()

    t0 = time.perf_counter()
    x_train, y_train, train_cols = _read_train(args.train, args.target)
    x_pred, pred_cols = _read_predict(args.predict)
    if train_cols != pred_cols:
        raise ValueError("Training and prediction feature columns differ.")

    model = _make_regressor(args.device, args.seed)
    model.fit(x_train, y_train)
    pred = np.asarray(model.predict(x_pred), dtype=np.float32).reshape(-1)

    pd.DataFrame({"prediction": pred}).to_csv(args.out, index=False)

    if args.metadata:
        with open(args.metadata, "w", encoding="utf-8") as handle:
            json.dump(
                {
                    "n_train": int(x_train.shape[0]),
                    "n_predict": int(x_pred.shape[0]),
                    "n_features": int(x_train.shape[1]),
                    "device": args.device,
                    "seed": int(args.seed),
                    "wall_sec": time.perf_counter() - t0,
                },
                handle,
                indent=2,
            )

    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"run_tabpfn_phylo.py: {exc}", file=sys.stderr)
        raise SystemExit(1)
