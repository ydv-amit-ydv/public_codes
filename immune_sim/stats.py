#!/usr/bin/env python3
"""
Statistical core for validating the pipeline.

Everything that produces a number a reader might act on should come through
here, so that the number arrives with an interval and a test attached.

Contents
--------
  bootstrap_ci     percentile bootstrap for any metric
  delong_auc_var   DeLong variance / covariance for ROC AUC
  delong_test      paired DeLong test: is AUC_a > AUC_b on the same cases?
  auc_ci_delong    analytic AUC confidence interval
  permutation_auc  exact-ish permutation p-value for AUC > 0.5
  brier            Brier score and its calibration/refinement decomposition
  calibration      reliability curve and expected calibration error
  bh_fdr           Benjamini-Hochberg correction across many alleles
  prob_direction   posterior probability that an effect has a given sign
"""
import numpy as np
from scipy import stats as sps

SEED = 0


# --------------------------------------------------------------- bootstrap
def bootstrap_ci(y, p, metric, n=2000, alpha=0.05, seed=SEED, paired_with=None):
    """
    Percentile bootstrap CI for metric(y, p).

    paired_with: a second score vector resampled on the SAME indices, in
    which case the metric is called metric(y, p, p2) and the interval is for
    the difference. Resampling cases (not scores) keeps the pairing.
    """
    rng = np.random.default_rng(seed)
    y = np.asarray(y); p = np.asarray(p)
    idx = np.arange(len(y))
    vals = []
    for _ in range(n):
        s = rng.choice(idx, size=len(idx), replace=True)
        if len(np.unique(y[s])) < 2:
            continue
        try:
            vals.append(metric(y[s], p[s]) if paired_with is None
                        else metric(y[s], p[s], np.asarray(paired_with)[s]))
        except Exception:
            continue
    if not vals:
        return float("nan"), float("nan"), []
    lo = float(np.percentile(vals, 100 * alpha / 2))
    hi = float(np.percentile(vals, 100 * (1 - alpha / 2)))
    return lo, hi, vals


# ------------------------------------------------------------------ DeLong
def _midrank(x):
    J = np.argsort(x)
    Z = x[J]
    N = len(x)
    T = np.zeros(N, dtype=float)
    i = 0
    while i < N:
        j = i
        while j < N and Z[j] == Z[i]:
            j += 1
        T[i:j] = 0.5 * (i + j - 1) + 1
        i = j
    out = np.empty(N, dtype=float)
    out[J] = T
    return out


def _structural(scores, m):
    """DeLong V10/V01 components. scores: (k, n) with positives first."""
    k, n = scores.shape
    nn = n - m
    tx = np.array([_midrank(scores[r, :m]) for r in range(k)])
    ty = np.array([_midrank(scores[r, m:]) for r in range(k)])
    tz = np.array([_midrank(scores[r, :]) for r in range(k)])
    aucs = (tz[:, :m].sum(axis=1) / m / nn - (m + 1.0) / (2.0 * nn))
    v01 = (tz[:, :m] - tx) / nn
    v10 = 1.0 - (tz[:, m:] - ty) / m
    s01 = np.cov(v01) if k > 1 else np.array([[np.var(v01, ddof=1)]])
    s10 = np.cov(v10) if k > 1 else np.array([[np.var(v10, ddof=1)]])
    s = np.atleast_2d(s01) / m + np.atleast_2d(s10) / nn
    return aucs, s


def _order(y, *score_sets):
    y = np.asarray(y).astype(int)
    order = np.argsort(-y, kind="mergesort")
    m = int(y.sum())
    return m, np.vstack([np.asarray(s)[order] for s in score_sets])


def auc_ci_delong(y, p, alpha=0.05):
    """Analytic AUC with a DeLong confidence interval."""
    m, S = _order(y, p)
    if m == 0 or m == len(y):
        return float("nan"), (float("nan"), float("nan"))
    auc, var = _structural(S, m)
    se = float(np.sqrt(np.clip(var[0, 0], 0, None)))
    z = sps.norm.ppf(1 - alpha / 2)
    a = float(auc[0])
    return a, (float(np.clip(a - z * se, 0, 1)), float(np.clip(a + z * se, 0, 1)))


def delong_test(y, p_a, p_b):
    """
    Paired DeLong test on the same cases.
    Returns (auc_a, auc_b, difference, two-sided p).
    """
    m, S = _order(y, p_a, p_b)
    if m == 0 or m == len(y):
        return (float("nan"),) * 4
    aucs, s = _structural(S, m)
    L = np.array([1.0, -1.0])
    var = float(L @ s @ L)
    diff = float(aucs[0] - aucs[1])
    if var <= 0:
        return float(aucs[0]), float(aucs[1]), diff, 1.0
    z = diff / np.sqrt(var)
    return float(aucs[0]), float(aucs[1]), diff, float(2 * sps.norm.sf(abs(z)))


def permutation_auc(y, p, n=5000, seed=SEED):
    """Permutation p-value for H0: AUC = 0.5 (label shuffling)."""
    from sklearn.metrics import roc_auc_score
    rng = np.random.default_rng(seed)
    y = np.asarray(y); obs = roc_auc_score(y, p)
    cnt = sum(roc_auc_score(rng.permutation(y), p) >= obs for _ in range(n))
    return float(obs), float((cnt + 1) / (n + 1))


# ------------------------------------------------- calibration of probabilities
def brier(y, p):
    """
    Brier score with Murphy's decomposition.
    BS = reliability - resolution + uncertainty  (lower reliability is better)
    """
    y = np.asarray(y, float); p = np.asarray(p, float)
    bs = float(np.mean((p - y) ** 2))
    bins = np.clip((p * 10).astype(int), 0, 9)
    base = y.mean()
    rel = res = 0.0
    for b in range(10):
        m = bins == b
        if not m.any():
            continue
        w = m.mean()
        rel += w * (p[m].mean() - y[m].mean()) ** 2
        res += w * (y[m].mean() - base) ** 2
    return {"brier": bs, "reliability": float(rel), "resolution": float(res),
            "uncertainty": float(base * (1 - base))}


def calibration(y, p, n_bins=10):
    """Reliability curve and expected / maximum calibration error."""
    y = np.asarray(y, float); p = np.asarray(p, float)
    edges = np.linspace(0, 1, n_bins + 1)
    rows, ece, mce = [], 0.0, 0.0
    for i in range(n_bins):
        m = (p >= edges[i]) & (p < edges[i + 1] if i < n_bins - 1 else p <= 1.0)
        if not m.any():
            continue
        conf, obs, w = p[m].mean(), y[m].mean(), m.mean()
        rows.append({"bin_lo": edges[i], "bin_hi": edges[i + 1], "n": int(m.sum()),
                     "mean_predicted": float(conf), "observed_frequency": float(obs)})
        ece += w * abs(conf - obs)
        mce = max(mce, abs(conf - obs))
    return {"bins": rows, "ece": float(ece), "mce": float(mce)}


# --------------------------------------------------------- multiple testing
def bh_fdr(pvals, alpha=0.05):
    """Benjamini-Hochberg. Returns (adjusted p, significant mask)."""
    p = np.asarray(pvals, float)
    n = len(p)
    order = np.argsort(p)
    adj = np.empty(n)
    prev = 1.0
    for rank, i in enumerate(order[::-1]):
        r = n - rank
        prev = min(prev, p[i] * n / r)
        adj[i] = prev
    return adj, adj <= alpha


# -------------------------------------------------------- effect direction
def prob_direction(samples, greater_than=0.0):
    """
    Share of a bootstrap / Monte-Carlo sample on one side of a threshold.

    Read this as: given the uncertainty actually sampled, how often does the
    effect come out in the stated direction. It is NOT a probability that the
    underlying biology is true - only the sampled sources of uncertainty are
    represented.
    """
    s = np.asarray(samples, float)
    s = s[np.isfinite(s)]
    if s.size == 0:
        return float("nan")
    return float((s > greater_than).mean())


def fmt_p(p):
    if not np.isfinite(p):
        return "n/a"
    if p < 1e-4:
        return "<1e-4"
    return f"{p:.3g}"
