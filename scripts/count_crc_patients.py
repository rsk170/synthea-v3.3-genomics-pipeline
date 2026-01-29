import pandas as pd
from pathlib import Path

OUTDIR = Path("./output_crc_run1/csv/")

patients = pd.read_csv(OUTDIR / "patients.csv")
conditions = pd.read_csv(OUTDIR / "conditions.csv")

# Your CRC codes from keep_crc.json
crc_snomed = {
    "93761005",   # Primary malignant neoplasm of colon
    "109838007",  # Overlapping malignant neoplasm of colon
    "363406005",  # Malignant neoplasm of colon
    "94260004",   # Metastatic malignant neoplasm to colon
}

# Make sure CODE is treated as string (sometimes it gets read as int)
conditions["CODE"] = conditions["CODE"].astype(str)

exported_patient_ids = set(patients["Id"].astype(str))

crc_patients = set(
    conditions.loc[
        conditions["PATIENT"].astype(str).isin(exported_patient_ids)
        & conditions["CODE"].isin(crc_snomed),
        "PATIENT"
    ].astype(str)
)

missing_crc = sorted(exported_patient_ids - crc_patients)

print(f"Exported patients.csv unique patients: {len(exported_patient_ids)}")
print(f"Patients with >=1 CRC SNOMED code in conditions.csv: {len(crc_patients)}")
print(f"Patients missing those CRC codes: {len(missing_crc)}")

# Show a few missing, if any
print("Example missing patient IDs (up to 20):", missing_crc[:20])

# Optional: save the missing list
pd.Series(missing_crc, name="PATIENT").to_csv(OUTDIR / "patients_missing_crc.csv", index=False)
print(f"Wrote: {OUTDIR / 'patients_missing_crc.csv'}")
