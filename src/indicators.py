# -*- coding: utf-8 -*-
"""四张科学价值图的指标聚合。"""

import numpy as np
import pandas as pd

from .config import INDICATOR_PARAMETERS


def indicator_scores(score_df):
    rows = []
    for _, r in score_df.iterrows():
        scores = {"site": r["site"]}
        for ind in ["C", "M", "G", "T"]:
            scores[ind] = float(np.mean([r[p] for p in INDICATOR_PARAMETERS[ind]]))
        rows.append(scores)
    return pd.DataFrame(rows).set_index("site")
