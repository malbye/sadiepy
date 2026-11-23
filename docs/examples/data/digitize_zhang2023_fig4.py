"""Digitize the gauge markers of Zhang et al. (2023), Fig. 4 (p. 7).

Provenance script for ``zhang2023_fig4_digitized.csv``; it is not run by the docs build.
Requires poppler's ``pdftoppm`` and Pillow, and the paper at ``refs/Zhang_et_al_2023.pdf``.

Method: render the page at 400 dpi; the axes frames were located from long dark line
segments (pixel constants below). For each marker, take the longest vertical run of the
marker's fill colour in a 17-px column at the gauge position, and convert its centre to
H/H0 with the axis calibration (0 to 1.2 over ``HH`` px). Markers partly hidden by a
fitted line fall back to the extent of the fill-plus-outline colour family. Two markers
are hidden behind the legend or another marker and were read by eye from a zoomed crop
(``MANUAL``). Each row records which method was used.

Usage::

    python digitize_zhang2023_fig4.py refs/Zhang_et_al_2023.pdf out.csv
"""

import csv
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
from PIL import Image

X0 = [374, 1351, 2327]  # left edge of the axes in each panel column, px
W = 743  # axes width (x/B from 0 to 1), px
YB = [2148, 3004, 3859]  # bottom edge (H/H0 = 0) of each panel row, px
HH = 615  # axes height (H/H0 from 0 to 1.2), px
CASES = [
    [(0.06, 1.4, 0.3), (0.08, 1.4, 0.3), (0.08, 1.6, 0.3)],
    [(0.08, 1.4, 0.4), (0.08, 1.6, 0.4), (0.10, 1.6, 0.4)],
    [(0.08, 1.4, 0.5), (0.08, 1.6, 0.5), (0.16, 1.8, 0.5)],
]
CONFIGS = [("M1", "M3"), ("M1", "M2", "M3"), ("M1", "M2", "M3")]
FILL = {"M1": (144, 192, 208), "M2": (192, 32, 32), "M3": (240, 162, 128)}
TOL = {"M1": 38, "M2": 45, "M3": 30}
X_OVER_B = [0.25, 0.5, 0.75, 1.0]
MANUAL = {((0.16, 1.8, 0.5), "M1", 0.25): 1.01, ((0.16, 1.8, 0.5), "M2", 0.25): 0.891}


def _longest_run(mask: np.ndarray, lo: int, hi: int, min_len: int) -> tuple[int, int] | None:
    best, start = None, None
    for y in range(lo, hi + 1):
        if mask[y] and start is None:
            start = y
        if (not mask[y] or y == hi) and start is not None:
            if y - start >= min_len and (best is None or y - start > best[1] - best[0]):
                best = (start, y)
            start = None
    return best


def digitize(page: np.ndarray) -> list[tuple[float, float, float, str, float, float, str]]:
    r, g, b = page[..., 0], page[..., 1], page[..., 2]
    family = {
        "M1": (b - r > 50) & (b > 150),
        "M2": (r - b > 100) & (g < 110),
        "M3": (r > 200) & (g > 110) & (g < 200) & (b < 170),
    }
    rows = []
    for i in range(3):
        for j in range(3):
            case = CASES[i][j]
            for cfg in CONFIGS[i]:
                fill = np.linalg.norm(page - np.array(FILL[cfg]), axis=2) < TOL[cfg]
                rows.append((*case, cfg, 0.0, 1.0, "reference"))
                for xb in X_OVER_B:
                    if (case, cfg, xb) in MANUAL:
                        rows.append((*case, cfg, xb, MANUAL[(case, cfg, xb)], "manual-occluded"))
                        continue
                    x = round(X0[j] + xb * W)
                    lo = YB[i] - int(1.06 / 1.2 * HH)  # below the legend
                    hi = YB[i] - int((0.47 if xb <= 0.6 else 0.05) / 1.2 * HH)  # above labels
                    run = _longest_run(fill[:, x - 8 : x + 9].sum(1) >= 5, lo, hi, 12)
                    method = "colour"
                    if run is None:
                        run = _longest_run(family[cfg][:, x - 3 : x + 4].sum(1) >= 3, lo, hi, 15)
                        method = "colour-outline"
                    if run is None:
                        raise RuntimeError(f"marker not found: {case} {cfg} x/B={xb}")
                    centre = (run[0] + run[1] - 1) / 2
                    rows.append((*case, cfg, xb, (YB[i] - centre) / HH * 1.2, method))
    return rows


def main(pdf: str, out: str) -> None:
    with tempfile.TemporaryDirectory() as tmp:
        stem = Path(tmp) / "p"
        subprocess.run(
            ["pdftoppm", "-f", "9", "-l", "9", "-r", "400", "-png", pdf, stem], check=True
        )
        png = next(Path(tmp).glob("p*.png"))
        page = np.asarray(Image.open(png).convert("RGB")).astype(float)
    with open(out, "w", newline="") as f:
        w = csv.writer(f, lineterminator="\n")
        w.writerow(["H_m", "T_s", "h_m", "config", "x_over_B", "H_over_H0", "method"])
        for row in digitize(page):
            w.writerow([*row[:5], f"{row[5]:.3f}", row[6]])


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
