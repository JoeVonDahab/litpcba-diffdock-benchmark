"""EF / BEDROC / ROC-AUC — shared metrics, matching the paper's definitions.
ponytail: one flat module; BEDROC delegated to rdkit's vetted CalcBEDROC.
"""
import numpy as np


def ef_at(y_true, scores, frac=0.01, higher_better=True):
    """Enrichment factor at top `frac`. y_true in {0,1}."""
    y = np.asarray(y_true, float)
    s = np.asarray(scores, float)
    n = len(y)
    order = np.argsort(-s if higher_better else s)
    k = max(1, int(frac * n))          # floor, matching the paper's int(total*0.01)
    n_act = y.sum()
    if n_act == 0 or k == 0:
        return np.nan
    # standard EF: (actives in top-k / k) / (total actives / total)
    return float((y[order][:k].sum() / k) / (n_act / n))


def roc_auc(y_true, scores, higher_better=True):
    from sklearn.metrics import roc_auc_score
    s = np.asarray(scores, float)
    y = np.asarray(y_true, int)
    if y.sum() == 0 or y.sum() == len(y):
        return np.nan
    return float(roc_auc_score(y, s if higher_better else -s))


def bedroc(y_true, scores, alpha=20.0, higher_better=True):
    """BEDROC via rdkit (Truchon & Bayly 2007), 0-1 scale."""
    from rdkit.ML.Scoring.Scoring import CalcBEDROC
    s = np.asarray(scores, float)
    y = np.asarray(y_true, int)
    order = np.argsort(-s if higher_better else s)
    scored = [(int(v),) for v in y[order]]  # sorted best-first, col 0 = activity
    return float(CalcBEDROC(scored, 0, alpha))


def all_metrics(y_true, scores, higher_better=True):
    return {
        "EF1%": ef_at(y_true, scores, 0.01, higher_better),
        "EF10%": ef_at(y_true, scores, 0.10, higher_better),
        "ROC_AUC": roc_auc(y_true, scores, higher_better),
        "BEDROC20": bedroc(y_true, scores, 20.0, higher_better),
        "n": int(len(y_true)),
        "n_actives": int(np.asarray(y_true).sum()),
    }


if __name__ == "__main__":
    y = np.array([1, 1, 0, 0, 0, 0, 0, 0, 0, 0])
    perfect = np.arange(10, 0, -1)
    assert abs(ef_at(y, perfect, 0.1) - 5.0) < 1e-9, ef_at(y, perfect, 0.1)  # top-1 of 10 has 1 active -> (1/2)/0.1=5
    assert abs(roc_auc(y, perfect) - 1.0) < 1e-9
    print("metrics self-check OK")
