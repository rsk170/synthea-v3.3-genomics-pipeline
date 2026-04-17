## Treatment Biomarker Realization Layer

### Summary
Keep Synthea’s current treatment-group assignment and treatment-budget logic as the internal scaffold, and add a separate post-generation realization layer that turns each assigned gene into a clinically realistic genomic event.

This plan avoids rewriting the core module now, preserves the existing clone/trend workflow, and fixes the main mismatch in your current pipeline:
- Synthea internal logic can still decide *which treatment-related slot exists*
- the realization layer decides *what exact event that slot becomes* in clinical output and genome input

The key design choice is:
- **small-variant MAF output is only for genome-embeddable SNV/MNV/indel events**
- **CNA/expression-style treatment biomarkers stay in clinical output as separate biomarker events**
- **exact administered drug is used post hoc**, not just the broad treatment group

### Key Changes
- Add a **two-tier evidence model** instead of overloading `breast_cancer_driver_variants_from_maf.csv` with everything.
  - **Raw evidence sheet**: literature/CIViC curation, including exact drug, group, gene, sensitive/resistant role, named variant/event, and source citation.
  - **Resolved runtime catalog**: machine-consumable rows that the pipeline can actually use, with one row per allowed treatment biomarker realization.

- The **resolved runtime catalog** should be keyed by:
  - exact drug
  - treatment group
  - gene
  - role (`sensitive` / `resistant`)
  - event type (`small_variant`, `cna`, `expression`)
  - realization mode (`maf_exact`, `group_level_nonvariant`)
  - exact MAF match when available (`helper_variant_key` or full row reference)
  - fallback metadata only when explicitly curated

- Use **actual administered drug after generation** to choose treatment biomarker evidence.
  - Read the patient’s real therapy from generated treatment data / `medications.csv`
  - For each treatment-related gene, choose the **earliest qualifying drug in that treatment group** as the canonical drug for realization
  - Match against the resolved catalog by exact drug first
  - If no exact-drug entry exists, only use a **curated group-level entry** for that same gene/role; do not automatically borrow evidence from unrelated drugs

- Keep the current post-processing split:
  - **random / non-treatment genes** keep using the current MAF-backed variant selection workflow
  - **treatment-related genes** use the new treatment realization catalog
  - one assigned treatment-related gene still consumes **one treatment-budget slot**, regardless of whether it becomes a small variant, CNA, or expression biomarker

- Build two downstream exports from the same realized patient events:
  - **Clinical genomics output**
    - small variants use the current mutation-block format
    - non-variant biomarkers use a separate biomarker block, not the gHGVS-based variant block
  - **Tumor genome input**
    - includes only small-variant MAF-like rows
    - merges realized treatment small variants with passenger mutations
    - excludes expression events
    - excludes CNAs for now

### Implementation Changes
- Extend the current observation-rewrite workflow so it first classifies each assigned gene as:
  - treatment-related with resolved event from the treatment catalog
  - or non-treatment/random with existing MAF-backed variant selection

- Preserve one underlying event per patient-gene across all retained timepoints.
  - The event identity stays fixed
  - only VAF / clonal role changes by timepoint

- Add a dedicated patient-level **realized events table** as the intermediate truth.
  - One row per patient-gene-event
  - Includes gene, event type, drug, role, chosen evidence source, exact variant reference if present, and whether it is genome-exportable
  - This becomes the source for both clinical observations and tumor-genome MAF assembly

- Keep `breast_cancer_driver_variants_from_maf.csv` as the small-variant source for:
  - random driver genes
  - treatment genes only when the chosen realization is an exact MAF-backed small variant
  - Do not use this file to represent CNAs or expression biomarkers

### Test Plan
- **CDK4/6 resistance case**: patient receives `palbociclib`, has `RB1`; realization selects a CDK4/6-resistant event and exports it as a small variant if catalog-backed by MAF.
- **Endocrine resistance case**: patient receives `tamoxifen` / AI / `fulvestrant`, has `NF1`; realization uses endocrine-resistance evidence and does not fall back to generic unknown behavior.
- **HER2 resistance case**: patient receives trastuzumab-class therapy, has `MET`; realization uses HER2-resistance evidence and exports it consistently.
- **Expression/CNA-only evidence case**: gene has no acceptable small variant but does have curated amplification/overexpression evidence; it still fills the treatment slot in clinical output, but does not appear in the genome MAF.
- **Unresolved named CIViC variant**: exact literature variant has no local MAF row; pipeline does not invent a synthetic row under the current policy, and instead uses only an explicitly curated group-level fallback or marks the event as non-variant clinical-only.
- **Random gene control**: non-treatment genes continue to use the current MAF-backed selection unchanged.
- **Multi-timepoint control**: the same patient-gene keeps the same realized event at baseline and follow-up timepoints; only VAF changes.

### Assumptions and Defaults
- **No rewrite of core Synthea treatment assignment now**; the plan assumes an overlay model.
- **MAF-first policy**: no synthetic small-variant rows are created in this phase.
- **Actual drug is preferred** over broad group when selecting treatment biomarker evidence.
- **Any genomic event counts toward the treatment budget**.
- **Genome simulator input is small-variant only** for now.
- **CNA/expression biomarkers remain clinical-only** until you explicitly add a separate CNA simulation path.

### Core Suggestion
Do not try to force every treatment-related gene into a disruptive small variant. That is what is making the workflow incoherent.

Instead:
- let Synthea keep deciding **which biologic slot exists**
- let the realization layer decide **what exact genomic event best represents that slot**
- let the genome export include only **events the simulator can actually consume**
