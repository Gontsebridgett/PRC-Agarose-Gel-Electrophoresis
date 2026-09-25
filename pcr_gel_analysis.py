"""
pcr_gel_analysis.py

Fits a DNA ladder standard curve (log10(size in bp) vs. migration distance)
and uses it to estimate PCR product sizes from measured migration distance,
flagging whether each product matches its expected amplicon size.
"""

import csv
import math
from pathlib import Path

DATA_PATH = Path(__file__).parent / "sample_data" / "gel_migration_data.csv"
MATCH_TOLERANCE_PCT = 15  # within +/-15% of expected size counts as a match


def load_data(path: Path):
    ladder, samples = [], []
    with open(path, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            if row["type"] == "ladder":
                ladder.append((float(row["migration_mm"]), float(row["known_size_bp"])))
            else:
                expected = float(row["expected_size_bp"]) if row["expected_size_bp"] else None
                samples.append((row["band_id"], float(row["migration_mm"]), expected))
    return ladder, samples


def linear_regression(x_vals, y_vals):
    n = len(x_vals)
    mean_x = sum(x_vals) / n
    mean_y = sum(y_vals) / n
    ss_xy = sum((x - mean_x) * (y - mean_y) for x, y in zip(x_vals, y_vals))
    ss_xx = sum((x - mean_x) ** 2 for x in x_vals)
    m = ss_xy / ss_xx
    b = mean_y - m * mean_x
    ss_tot = sum((y - mean_y) ** 2 for y in y_vals)
    ss_res = sum((y - (m * x + b)) ** 2 for x, y in zip(x_vals, y_vals))
    r_squared = 1 - (ss_res / ss_tot) if ss_tot else 1.0
    return m, b, r_squared


def main():
    ladder, samples = load_data(DATA_PATH)

    distances = [d for d, _ in ladder]
    log_sizes = [math.log10(size) for _, size in ladder]
    m, b, r_squared = linear_regression(distances, log_sizes)

    print(f"Standard curve fit: log10(bp) = {m:.4f} * distance + {b:.3f}   (R-squared = {r_squared:.3f})\n")

    for band_id, distance, expected in samples:
        estimated_size = 10 ** (m * distance + b)
        line = f"{band_id}: migration {distance} mm  ->  Estimated size: ~{estimated_size:.0f} bp"
        if expected:
            diff_pct = abs(estimated_size - expected) / expected * 100
            if diff_pct <= MATCH_TOLERANCE_PCT:
                line += f"  (expected: {expected:.0f} bp, match)"
            else:
                line += f"  (expected: {expected:.0f} bp, mismatch \u2014 check primers)"
        print(line)


if __name__ == "__main__":
    main()
