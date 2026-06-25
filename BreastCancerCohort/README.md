# Representative Breast Cancer Cohort

This directory contains the representative synthetic cohort generated and distributed with this repository. The clinical CSV files are under
`Clinical Data/`, and the corresponding patient- and clone-level mutation files
are under `MAF files/`.

## Omitted intermediate file

`breast_cancer_assigned_passenger_mutations.tsv` is a large intermediate table
and is not included in this packaged cohort. The passenger assignments from that
table have already been split across clones and are present in the clone-level
MAF files under `MAF files/`. The omission therefore does not affect use of the
released final MAF files, but the intermediate table is required to regenerate
those MAF files from scratch.

Regeneration requires the PCAWG/ICGC resources described in the main
repository `README.md`.
