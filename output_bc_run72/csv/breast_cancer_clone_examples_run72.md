# Breast Cancer Clone/VAF Examples from `run72`

This file highlights a small set of representative patients from [breast_cancer_clone_proportions.csv](/home/rkisleva/SDG_tools/synthea_developer3.3.0/synthea/output_runs/output_bc_run72/csv/breast_cancer_clone_proportions.csv) to show how the clone-building and VAF logic works across different scenarios.

## At a Glance

| Patient ID | Why included | Timepoints | Clones | Total genes |
| --- | --- | --- | --- | --- |
| `61eec2ea-b8e1-fbce-2342-62cb62811619` | Simplest case: 2 genes, 2 clones, clear decreasing clone | 2 | 2 | 2 |
| `5c15c240-389e-b191-cd61-3b38b9c132b4` | 2 timepoints with unknown, increasing, and stable clone behavior | 2 | 3 | 5 |
| `1fbc5bb3-8af5-989c-5dc0-6759cc90215c` | Higher mutation burden with 4 clones | 2 | 4 | 13 |
| `60f8ea69-2e25-5ea3-44f7-150dff7e38a2` | 3 timepoints with decreasing and stable-then-decreasing trajectories | 3 | 3 | 3 |
| `9125f158-2fb2-7b53-4dee-0fa41b47e685` | 3 timepoints with both increasing and decreasing late clones | 3 | 3 | 5 |

## Patient 1: Minimal two-clone case

Patient: `61eec2ea-b8e1-fbce-2342-62cb62811619`  
Profile: `ER+ / PR- / HER2-`, Stage `1A`
Sequencing dates:
- `t0 = 2023-04-11T08:54:17Z`
- `t1 = 2023-11-01T12:34:14Z`

| Clone | Type | Genes | Signature | t0 VAF % | t1 VAF % |
| --- | --- | --- | --- | ---: | ---: |
| `clone_1` | founding | `TP53` | `t0:present|t1:unknown` | 46.46 | 97.30 |
| `clone_2` | late | `PGR` | `t0:present|t1:decreasing` | 53.54 | 2.70 |

Why it is useful:
- simplest example of a founding clone plus one late decreasing clone
- easy to see the pruning threshold effect at the second sequencing

## Patient 2: Two timepoints with increasing and stable clones

Patient: `5c15c240-389e-b191-cd61-3b38b9c132b4`  
Profile: `ER+ / PR+ / HER2-`, Stage `4`
Sequencing dates:
- `t0 = 2008-03-17T13:09:18Z`
- `t1 = 2008-09-19T22:13:29Z`

| Clone | Type | Genes | Signature | t0 VAF % | t1 VAF % |
| --- | --- | --- | --- | ---: | ---: |
| `clone_1` | founding | `FOXQ1;NF1;NOTCH2` | `t0:present|t1:unknown` | 80.36 | 40.91 |
| `clone_2` | late | `NCOR2` | `t0:present|t1:increasing` | 2.13 | 44.09 |
| `clone_3` | late | `RB1` | `t0:present|t1:stable` | 17.51 | 15.00 |

Why it is useful:
- shows an increasing clone starting very small at baseline
- shows a stable clone staying near its initial VAF
- shows the founding clone losing relative share as another clone expands

## Patient 3: Higher mutation burden with four clones

Patient: `1fbc5bb3-8af5-989c-5dc0-6759cc90215c`  
Profile: `ER- / PR- / HER2+`, Stage `1A`
Sequencing dates:
- `t0 = 2017-01-26T22:19:52Z`
- `t1 = 2017-11-16T07:09:56Z`

| Clone | Type | Genes | Signature | t0 VAF % | t1 VAF % |
| --- | --- | --- | --- | ---: | ---: |
| `clone_1` | founding | `CDH1;PTEN;TP53` | `t0:present|t1:unknown` | 27.53 | 36.10 |
| `clone_2` | branch | `CCND1;KMT2C;MAP3K1` | `t0:present|t1:unknown` | 19.44 | 30.33 |
| `clone_3` | branch | `GATA3;MAP2K4;ZNF703` | `t0:present|t1:unknown` | 19.44 | 30.33 |
| `clone_4` | late | `ERBB3;ESR2;HER2;YWHAZ` | `t0:present|t1:decreasing` | 33.59 | 3.24 |

Why it is useful:
- example with many total drivers split across founding, branch, and late clones
- late treatment-related clone starts substantial and then collapses
- useful to explain how the script avoids putting all mutations into one clone

## Patient 4: Three timepoints with stable-then-decreasing behavior

Patient: `60f8ea69-2e25-5ea3-44f7-150dff7e38a2`  
Profile: `ER+ / PR- / HER2-`, Stage `3C`
Sequencing dates:
- `t0 = 2017-07-29T22:37:09Z`
- `t1 = 2018-01-24T22:06:39Z`
- `t2 = 2018-06-26T14:30:47Z`

| Clone | Type | Genes | Signature | t0 VAF % | t1 VAF % | t2 VAF % |
| --- | --- | --- | --- | ---: | ---: | ---: |
| `clone_1` | founding | `LINC00290` | `t0:present|t1:unknown|t2:unknown` | 32.60 | 51.34 | 95.04 |
| `clone_2` | late | `ABCB1` | `t0:present|t1:decreasing|t2:decreasing` | 29.83 | 11.14 | 3.33 |
| `clone_3` | late | `PGR` | `t0:present|t1:stable|t2:decreasing` | 37.57 | 37.52 | 1.63 |

Why it is useful:
- shows a stable midpoint followed by later collapse
- shows a 3-timepoint trajectory rather than only baseline/follow-up
- illustrates why clone signatures are tracked across all available timepoints

## Patient 5: Three timepoints with both increasing and decreasing late clones

Patient: `9125f158-2fb2-7b53-4dee-0fa41b47e685`  
Profile: `ER+ / PR+ / HER2-`, Stage `4`
Sequencing dates:
- `t0 = 2003-12-23T01:56:37Z`
- `t1 = 2004-06-15T20:16:26Z`
- `t2 = 2004-10-10T01:16:26Z`

| Clone | Type | Genes | Signature | t0 VAF % | t1 VAF % | t2 VAF % |
| --- | --- | --- | --- | ---: | ---: | ---: |
| `clone_1` | founding | `ARID1A;TP53` | `t0:present|t1:unknown|t2:unknown` | 46.63 | 22.05 | 46.37 |
| `clone_2` | late | `NCOR2` | `t0:present|t1:increasing|t2:increasing` | 3.18 | 31.26 | 51.79 |
| `clone_3` | late | `NF1;RB1` | `t0:present|t1:stable|t2:decreasing` | 50.19 | 46.69 | 1.84 |

Why it is useful:
- clean example with both expansion and collapse in the same patient
- the increasing clone starts under 5% at baseline and becomes dominant later
- the stable-then-decreasing clone stays high initially, then falls below the threshold

## Short Interpretation

- Founding clones usually carry the early/trunk events and often remain present across all timepoints.
- Late increasing clones are seeded small and then expand.
- Late decreasing clones begin larger and are pushed low at later sequencing.
- Stable clones are kept close to their previous value until the signature says they should fall.
- The total clone VAF within each patient/timepoint is normalized to sum to 100%.
