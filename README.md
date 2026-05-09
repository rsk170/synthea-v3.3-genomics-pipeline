# Synthea<sup>TM</sup> Patient Generator ![Build Status](https://github.com/synthetichealth/synthea/workflows/.github/workflows/ci-build-test.yml/badge.svg?branch=master) [![codecov](https://codecov.io/gh/synthetichealth/synthea/branch/master/graph/badge.svg)](https://codecov.io/gh/synthetichealth/synthea)

Synthea<sup>TM</sup> is a Synthetic Patient Population Simulator. The goal is to output synthetic, realistic (but not real), patient data and associated health records in a variety of formats.

Read our [wiki](https://github.com/synthetichealth/synthea/wiki) and [Frequently Asked Questions](https://github.com/synthetichealth/synthea/wiki/Frequently-Asked-Questions) for more information.

Currently, Synthea<sup>TM</sup> features include:
- Birth to Death Lifecycle
- Configuration-based statistics and demographics (defaults with Massachusetts Census data)
- Modular Rule System
  - Drop in [Generic Modules](https://github.com/synthetichealth/synthea/wiki/Generic-Module-Framework)
  - Custom Java rules modules for additional capabilities
- Primary Care Encounters, Emergency Room Encounters, and Symptom-Driven Encounters
- Conditions, Allergies, Medications, Vaccinations, Observations/Vitals, Labs, Procedures, CarePlans
- Formats
  - HL7 FHIR (R4, STU3 v3.0.1, and DSTU2 v1.0.2)
  - Bulk FHIR in ndjson format (set `exporter.fhir.bulk_data = true` to activate)
  - C-CDA (set `exporter.ccda.export = true` to activate)
  - CSV (set `exporter.csv.export = true` to activate)
  - CPCDS (set `exporter.cpcds.export = true` to activate)
- Rendering Rules and Disease Modules with Graphviz

## Developer Quick Start

These instructions are intended for those wishing to examine the Synthea source code, extend it or build the code locally. Those just wishing to run Synthea should follow the [Basic Setup and Running](https://github.com/synthetichealth/synthea/wiki/Basic-Setup-and-Running) instructions instead.

### Installation

**System Requirements:**
Synthea<sup>TM</sup> requires Java JDK 11 or newer. We strongly recommend using a Long-Term Support (LTS) release of Java, 11 or 17, as issues may occur with more recent non-LTS versions.

To clone the Synthea<sup>TM</sup> repo, then build and run the test suite:
```
git clone https://github.com/synthetichealth/synthea.git
cd synthea
./gradlew build check test
```

### Changing the default properties


The default properties file values can be found at `src/main/resources/synthea.properties`.
By default, synthea does not generate CCDA, CPCDA, CSV, or Bulk FHIR (ndjson). You'll need to
adjust this file to activate these features.  See the [wiki](https://github.com/synthetichealth/synthea/wiki)
for more details, or use our [guided customizer tool](https://synthetichealth.github.io/spt/#/customizer).



### Generate Synthetic Patients
Generating the population one at a time...
```
./run_synthea
```

Command-line arguments may be provided to specify a state, city, population size, or seed for randomization.
```
run_synthea [-s seed] [-p populationSize] [state [city]]
```

Full usage info can be printed by passing the `-h` option.
```
$ ./run_synthea -h     

> Task :run
Usage: run_synthea [options] [state [city]]
Options: [-s seed]
         [-cs clinicianSeed]
         [-p populationSize]
         [-r referenceDate as YYYYMMDD]
         [-g gender]
         [-a minAge-maxAge]
         [-o overflowPopulation]
         [-c localConfigFilePath]
         [-d localModulesDirPath]
         [-i initialPopulationSnapshotPath]
         [-u updatedPopulationSnapshotPath]
         [-t updateTimePeriodInDays]
         [-f fixedRecordPath]
         [-k keepMatchingPatientsPath]
         [--config*=value]
          * any setting from src/main/resources/synthea.properties

Examples:
run_synthea Massachusetts
run_synthea Alaska Juneau
run_synthea -s 12345
run_synthea -p 1000
run_synthea -s 987 Washington Seattle
run_synthea -s 21 -p 100 Utah "Salt Lake City"
run_synthea -g M -a 60-65
run_synthea -p 10 --exporter.fhir.export=true
run_synthea --exporter.baseDirectory="./output_tx/" Texas
```

Some settings can be changed in `./src/main/resources/synthea.properties`.

Synthea<sup>TM</sup> will output patient records in C-CDA and FHIR formats in `./output`.

### Synthea<sup>TM</sup> GraphViz
Generate graphical visualizations of Synthea<sup>TM</sup> rules and modules.
```
./gradlew graphviz
```

### Concepts and Attributes
Generate a list of concepts (used in the records) or attributes (variables on each patient).
```
./gradlew concepts
./gradlew attributes
```

## Breast Cancer Genomics Workflow

This repository also contains a custom breast cancer genomics pipeline built on
top of the standard Synthea CSV export. The final outputs of this workflow are:

- `observations_pruned_by_clone_vaf.csv`
- `breast_cancer_assigned_passenger_mutations.tsv`
- per-patient, per-clone MAF files under `maf_files/`

The workflow assumes that CSV export is enabled and that a run folder already
exists under `output_runs/<run_name>/csv/`.

### 1. Generate the Synthea run

First generate the breast cancer cohort with Synthea so that at minimum the run
contains:

- `observations.csv`
- `medications.csv`
- `patients.csv`

Example:

```bash
./run_synthea -p 30 --exporter.csv.export=true --exporter.baseDirectory="./output_runs/output_bc_run1/"
```

In the examples below, replace `output_bc_run1` with your run folder name.

### 2. Reconstruct clone groups

This step reads baseline genomic observations and later post-treatment clonal
trend observations from `observations.csv` and builds patient-specific
`founding`, `branch`, and `late` clone groups.

```bash
python3 scripts/build_breast_cancer_clones.py \
  --observations output_runs/output_bc_run1/csv/observations.csv \
  --output output_runs/output_bc_run1/csv/breast_cancer_clone_groups.csv
```

Output:

- `output_runs/output_bc_run1/csv/breast_cancer_clone_groups.csv`

### 3. Assign clone proportions over time

This step assigns relative abundance trajectories to the reconstructed clones
across sequencing timepoints and normalizes them so that clone proportions sum
to 100% at each timepoint.

```bash
python3 scripts/build_breast_cancer_clone_proportions.py \
  --clone-groups output_runs/output_bc_run1/csv/breast_cancer_clone_groups.csv \
  --output output_runs/output_bc_run1/csv/breast_cancer_clone_proportions.csv
```

Output:

- `output_runs/output_bc_run1/csv/breast_cancer_clone_proportions.csv`

### 4. Create the pruned variant-oriented observations file

This step rewrites genomic rows in `observations.csv` into standardized
variant-oriented observation groups. It keeps only driver genes whose estimated
clonal prevalence remains above the configured threshold at each sequencing
event, and assigns exact variants from the CIViC, generic driver, or
non-disruptive variant libraries.

```bash
python3 scripts/build_breast_cancer_pruned_observations.py \
  --clone-proportions output_runs/output_bc_run1/csv/breast_cancer_clone_proportions.csv \
  --observations output_runs/output_bc_run1/csv/observations.csv \
  --medications output_runs/output_bc_run1/csv/medications.csv \
  --civic-driver-variants scripts/civic_breast_cancer_driver_variants_from_maf.csv \
  --driver-variants scripts/breast_cancer_driver_variants_from_maf.csv \
  --non-disruptive-variants scripts/breast_cancer_non_disruptive_variants_from_maf.csv \
  --output output_runs/output_bc_run1/csv/observations_pruned_by_clone_vaf.csv
```

Output:

- `output_runs/output_bc_run1/csv/observations_pruned_by_clone_vaf.csv`

### 5. Assign passenger mutations

This step creates a breast-cancer passenger-only mutation pool from the ICGC
MAF by removing Table S3 driver genes, then deterministically assigns each
synthetic patient with genomic sequencing to one source breast tumor sample and
 transfers that sample’s passenger mutation set to the patient.

```bash
python3 scripts/build_breast_cancer_passenger_mutations.py \
  --patients output_runs/output_bc_run1/csv/patients.csv \
  --observations output_runs/output_bc_run1/csv/observations.csv \
  --passenger-maf scripts/breast_cancer_passenger_only_from_maf.maf \
  --output output_runs/output_bc_run1/csv/breast_cancer_assigned_passenger_mutations.tsv
```

Outputs:

- `scripts/breast_cancer_passenger_only_from_maf.maf`
- `output_runs/output_bc_run1/csv/breast_cancer_assigned_passenger_mutations.tsv`

### 6. Create complete per-patient clone MAF files

This step combines:

- the driver variants present in `observations_pruned_by_clone_vaf.csv`, grouped
  by reconstructed clone membership
- a non-overlapping split of the assigned passenger mutation set for the same
  patient

It writes one MAF file per patient per clone under a new `maf_files/`
directory in the run root.

```bash
python3 scripts/build_breast_cancer_complete_maf_files.py \
  --pruned-observations output_runs/output_bc_run1/csv/observations_pruned_by_clone_vaf.csv \
  --assigned-passengers output_runs/output_bc_run1/csv/breast_cancer_assigned_passenger_mutations.tsv \
  --clone-groups output_runs/output_bc_run1/csv/breast_cancer_clone_groups.csv \
  --clone-proportions output_runs/output_bc_run1/csv/breast_cancer_clone_proportions.csv
```

Outputs:

- `output_runs/output_bc_run1/maf_files/<patient_uuid>/clone_*.maf`

Each clone-specific MAF contains:

- the driver variants for that clone that are recoverable from
  `observations_pruned_by_clone_vaf.csv`
- a randomly split, non-overlapping subset of the patient's passenger mutation
  rows, so passenger mutations do not repeat across clone files

### Minimal command sequence

After the Synthea run has been generated, the full genomics post-processing
sequence is:

```bash
python3 scripts/build_breast_cancer_clones.py \
  --observations output_runs/output_bc_run1/csv/observations.csv \
  --output output_runs/output_bc_run1/csv/breast_cancer_clone_groups.csv

python3 scripts/build_breast_cancer_clone_proportions.py \
  --clone-groups output_runs/output_bc_run1/csv/breast_cancer_clone_groups.csv \
  --output output_runs/output_bc_run1/csv/breast_cancer_clone_proportions.csv

python3 scripts/build_breast_cancer_pruned_observations.py \
  --clone-proportions output_runs/output_bc_run1/csv/breast_cancer_clone_proportions.csv \
  --observations output_runs/output_bc_run1/csv/observations.csv \
  --medications output_runs/output_bc_run1/csv/medications.csv \
  --civic-driver-variants scripts/civic_breast_cancer_driver_variants_from_maf.csv \
  --driver-variants scripts/breast_cancer_driver_variants_from_maf.csv \
  --non-disruptive-variants scripts/breast_cancer_non_disruptive_variants_from_maf.csv \
  --output output_runs/output_bc_run1/csv/observations_pruned_by_clone_vaf.csv

python3 scripts/build_breast_cancer_passenger_mutations.py \
  --patients output_runs/output_bc_run1/csv/patients.csv \
  --observations output_runs/output_bc_run1/csv/observations.csv \
  --passenger-maf scripts/breast_cancer_passenger_only_from_maf.maf \
  --output output_runs/output_bc_run1/csv/breast_cancer_assigned_passenger_mutations.tsv

python3 scripts/build_breast_cancer_complete_maf_files.py \
  --pruned-observations output_runs/output_bc_run1/csv/observations_pruned_by_clone_vaf.csv \
  --assigned-passengers output_runs/output_bc_run1/csv/breast_cancer_assigned_passenger_mutations.tsv \
  --clone-groups output_runs/output_bc_run1/csv/breast_cancer_clone_groups.csv \
  --clone-proportions output_runs/output_bc_run1/csv/breast_cancer_clone_proportions.csv
```

# License

Copyright 2017-2023 The MITRE Corporation

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
