# CIViC Breast Cancer Variant Summary

Source:
- `scripts/civic_breast_cancer_driver_variants_from_maf.csv`

Variants are listed using the `helper_variant_key` format, for example:
- `TP53|17|7578406|7578406|C|T|Missense_Mutation`

## Summary Table

| Medication | Gene | Sensitivity Variants | Resistance Variants |
|---|---|---|---|
| `Doxorubicin` | `TP53` | `TP53\|17\|7578406\|7578406\|C\|T\|Missense_Mutation` |  |
| `Fulvestrant` | `AKT1` | `AKT1\|14\|105246551\|105246551\|C\|T\|Missense_Mutation` |  |
| `Fulvestrant` | `PIK3CA` | `PIK3CA\|3\|178916876\|178916876\|G\|A\|Missense_Mutation`; `PIK3CA\|3\|178921553\|178921553\|T\|A\|Missense_Mutation`; `PIK3CA\|3\|178927980\|178927980\|T\|C\|Missense_Mutation`; `PIK3CA\|3\|178928079\|178928079\|G\|A\|Missense_Mutation`; `PIK3CA\|3\|178936082\|178936082\|G\|A\|Missense_Mutation`; `PIK3CA\|3\|178936091\|178936091\|G\|A\|Missense_Mutation`; `PIK3CA\|3\|178936094\|178936094\|C\|A\|Missense_Mutation`; `PIK3CA\|3\|178936095\|178936095\|A\|C\|Missense_Mutation`; `PIK3CA\|3\|178936095\|178936095\|A\|G\|Missense_Mutation`; `PIK3CA\|3\|178952084\|178952084\|C\|T\|Missense_Mutation`; `PIK3CA\|3\|178952085\|178952085\|A\|G\|Missense_Mutation`; `PIK3CA\|3\|178952085\|178952085\|A\|T\|Missense_Mutation`; `PIK3CA\|3\|178952090\|178952090\|G\|C\|Missense_Mutation` |  |
| `Lapatinib` | `ERBB2` | `ERBB2\|17\|37879658\|37879658\|G\|A\|Missense_Mutation` |  |
| `Lapatinib` | `PIK3CA` | `PIK3CA\|3\|178952085\|178952085\|A\|G\|Missense_Mutation` | `PIK3CA\|3\|178916946\|178916946\|G\|C\|Missense_Mutation`; `PIK3CA\|3\|178936091\|178936091\|G\|A\|Missense_Mutation`; `PIK3CA\|3\|178952085\|178952085\|A\|G\|Missense_Mutation` |
| `Neratinib` | `ERBB2` | `ERBB2\|17\|37879658\|37879658\|G\|A\|Missense_Mutation`; `ERBB2\|17\|37880261\|37880261\|G\|T\|Missense_Mutation`; `ERBB2\|17\|37881332\|37881332\|G\|A\|Missense_Mutation`; `ERBB2\|17\|37881414\|37881414\|T\|G\|Missense_Mutation` |  |
| `Palbociclib` | `PIK3CA` | `PIK3CA\|3\|178916876\|178916876\|G\|A\|Missense_Mutation`; `PIK3CA\|3\|178921553\|178921553\|T\|A\|Missense_Mutation`; `PIK3CA\|3\|178928079\|178928079\|G\|A\|Missense_Mutation`; `PIK3CA\|3\|178936082\|178936082\|G\|A\|Missense_Mutation`; `PIK3CA\|3\|178936091\|178936091\|G\|A\|Missense_Mutation`; `PIK3CA\|3\|178936094\|178936094\|C\|A\|Missense_Mutation`; `PIK3CA\|3\|178952085\|178952085\|A\|G\|Missense_Mutation`; `PIK3CA\|3\|178952090\|178952090\|G\|C\|Missense_Mutation` |  |
| `Tamoxifen` | `TP53` |  | `TP53\|17\|7577120\|7577120\|C\|T\|Missense_Mutation` |
| `Trastuzumab` | `PIK3CA` | `PIK3CA\|3\|178952085\|178952085\|A\|G\|Missense_Mutation` | `PIK3CA\|3\|178916946\|178916946\|G\|C\|Missense_Mutation`; `PIK3CA\|3\|178927980\|178927980\|T\|C\|Missense_Mutation`; `PIK3CA\|3\|178936082\|178936082\|G\|A\|Missense_Mutation`; `PIK3CA\|3\|178936091\|178936091\|G\|A\|Missense_Mutation`; `PIK3CA\|3\|178952085\|178952085\|A\|G\|Missense_Mutation` |

## Notes

- Some variants appear in both sensitivity and resistance contexts for the same gene and medication.
- The CIViC file currently contains rows for these medications:
  - `Doxorubicin`
  - `Fulvestrant`
  - `Lapatinib`
  - `Neratinib`
  - `Palbociclib`
  - `Tamoxifen`
  - `Trastuzumab`
