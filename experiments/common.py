# -*- coding: utf-8 -*-
"""实验公共工具。"""

import os
import sys

import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

DATA_DIR = os.path.join(ROOT, "data")
OUTPUT_DIR = os.path.join(ROOT, "outputs")


def load_sites(name):
    return pd.read_csv(os.path.join(DATA_DIR, name))


def ensure_output():
    os.makedirs(OUTPUT_DIR, exist_ok=True)


def save(df, name):
    ensure_output()
    df.to_csv(os.path.join(OUTPUT_DIR, name), index=False, encoding="utf-8-sig")


def show_ranking(title, res, score_col="isvm"):
    tmp = res[["site", score_col]].sort_values(score_col, ascending=False).reset_index(drop=True)
    tmp["rank"] = tmp[score_col].rank(ascending=False, method="min").astype(int)
    print(f"\n=== {title} ===")
    for _, r in tmp.iterrows():
        print(f"  {int(r['rank'])}. {r['site']}: {r[score_col]:.4f}")
    return tmp
