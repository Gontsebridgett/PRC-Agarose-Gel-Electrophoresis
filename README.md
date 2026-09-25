# PCR Product Analysis by Agarose Gel Electrophoresis

A protocol for confirming PCR product size using agarose gel electrophoresis,
paired with a Python script that builds a DNA ladder standard curve and
estimates PCR product size(s) from measured migration distance — the DNA
equivalent of the SDS-PAGE protein sizing workflow.

## Overview

After PCR amplification, agarose gel electrophoresis is used to confirm that
the expected product was generated, at the expected size. DNA fragments are
negatively charged and migrate toward the anode through an agarose gel matrix
under an electric field; smaller fragments migrate faster (travel further)
than larger fragments. Running a DNA ladder of known fragment sizes alongside
PCR products allows the size of each product band to be estimated by
comparison.

## Principle

- Agarose gel acts as a molecular sieve, similarly to polyacrylamide in
  SDS-PAGE, but typically used for larger, longer nucleic acid fragments
- Migration distance is inversely related to **log10(fragment size in bp)**
  for the ladder standards, forming a standard curve
- A single sharp band at the expected size confirms successful, specific
  amplification; multiple bands or smearing indicate non-specific
  amplification, primer dimers, or degraded template

## Materials & Reagents

- PCR product(s) and a DNA ladder (e.g. 100 bp or 1 kb ladder)
- Agarose powder and TAE or TBE running buffer
- DNA loading dye
- DNA-safe stain (e.g. SYBR Safe, GelRed, or ethidium bromide)
- Horizontal gel electrophoresis apparatus and power supply
- Gel imaging system (UV or blue-light transilluminator)

## Method (Summary)

| Step | Action |
|---|---|
| 1 | Prepare agarose gel (typically 1&ndash;2%, higher % for smaller fragments) with stain |
| 2 | Pour gel with comb inserted; allow to set |
| 3 | Load PCR products (mixed with loading dye) and DNA ladder into wells |
| 4 | Run at constant voltage (~80&ndash;120V) until dye front has migrated appropriately |
| 5 | Image gel under UV/blue light |
| 6 | Measure migration distances and compare to the ladder standard curve |

## Result Interpretation

| Observation | Interpretation |
|---|---|
| Single sharp band at expected size | Successful, specific amplification |
| Band at unexpected size | Non-specific amplification or primer mis-annealing |
| Multiple bands | Non-specific amplification or multiple targets amplified |
| Smearing | Degraded template or non-optimal PCR conditions |
| No band | Failed amplification (primer, template, or reaction issue) |
| Very low band (~50&ndash;100 bp), present without template | Primer dimers |

## Analysis Script

`pcr_gel_analysis.py` fits a standard curve from DNA ladder migration
distances (log10(size) vs. distance) and uses it to estimate the size of PCR
product bands from their measured migration distance — confirming whether the
product matches the expected amplicon size.

### Usage

```bash
pip install -r requirements.txt
python pcr_gel_analysis.py
```

### Sample output

```
Standard curve fit: log10(bp) = -0.0170 * distance + 3.397   (R-squared = 0.998)

PCR_Product_1: migration 58.0 mm  ->  Estimated size: ~258 bp  (expected: 230 bp, match)
PCR_Product_2: migration 75.0 mm  ->  Estimated size: ~133 bp  (expected: 230 bp, mismatch — check primers)
```

## Repository Structure

```
pcr-agarose-gel-electrophoresis/
├── README.md
├── pcr_gel_analysis.py
├── requirements.txt
└── sample_data/
    └── gel_migration_data.csv
```

## Relationship to Other Projects in This Portfolio

This project builds on [`pcr-simulation`](../pcr-simulation), which predicts
the expected amplicon size from primer positions — the two together mirror
the full real-world workflow of running PCR, then confirming the product by
gel electrophoresis.
