from __future__ import annotations

import numpy as np
from scipy import sparse


def mean_axis0(x):
    if sparse.issparse(x):
        return np.asarray(x.mean(axis=0)).ravel()
    return np.asarray(x).mean(axis=0)


def mean_axis1(x):
    if sparse.issparse(x):
        return np.asarray(x.mean(axis=1)).ravel()
    return np.asarray(x).mean(axis=1)
