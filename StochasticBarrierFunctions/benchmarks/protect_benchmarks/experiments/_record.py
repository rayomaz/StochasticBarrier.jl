"""Helper to record PRoTECT experiment results to a unified CSV file."""

import csv
import os
import time
import traceback

_HEADER = [
    "tool", "table", "class", "system",
    "barrier_type", "pwc_method", "degree", "num_partitions",
    "time_s", "eta", "beta", "Ps", "status",
]


def _results_path():
    out_dir = os.environ.get("SB_RESULTS_DIR", "/results")
    os.makedirs(out_dir, exist_ok=True)
    return os.path.join(out_dir, "protect_results.csv")


def _ensure_header(path):
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        with open(path, "w", newline="") as f:
            csv.writer(f).writerow(_HEADER)


def _extract(result):
    """Map PRoTECT dt_SS result dict → (eta, beta, Ps, status).

    PRoTECT keys: gamma (≈ η), lambda (≈ β), confidence (≈ Ps).
    Missing keys or an error field → FAILED.
    """
    if not isinstance(result, dict) or not result or "error" in result:
        return "", "", "", "FAILED"
    eta = result.get("gamma", "")
    beta = result.get("lambda", "")
    Ps = result.get("confidence", "")
    if eta == "" or beta == "":
        return "", "", "", "FAILED"
    return eta, beta, Ps, "OK"


def record(dt_ss_callable, *, system, cls, table, degree,
           barrier_type="SOS", pwc_method="", num_partitions="NA",
           **dt_ss_kwargs):
    """Run PRoTECT's dt_SS (or equivalent) and append a row to the CSV.

    `dt_ss_callable` is invoked with **dt_ss_kwargs; the first positional
    argument is the degree. Timing is measured here to stay consistent
    across tools.
    """
    path = _results_path()
    _ensure_header(path)

    start = time.time()
    try:
        result = dt_ss_callable(degree, **dt_ss_kwargs)
        elapsed = time.time() - start
        eta, beta, Ps, status = _extract(result)
    except Exception:
        elapsed = time.time() - start
        traceback.print_exc()
        eta, beta, Ps, status = "", "", "", "FAILED"
        result = None

    with open(path, "a", newline="") as f:
        csv.writer(f).writerow([
            "protect", str(table), cls, system,
            barrier_type, pwc_method, str(degree), str(num_partitions),
            f"{elapsed:.4f}", str(eta), str(beta), str(Ps), status,
        ])

    return result
