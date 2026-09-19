# Data access and governance

Obtain MIMIC-IV, eICU-CRD and SICdb through their respective access, training and data-use procedures. This repository cannot grant access or redistribute governed data.

- MIMIC-IV: https://physionet.org/content/mimiciv/
- eICU-CRD: https://physionet.org/content/eicu-crd/
- SICdb: https://www.sicdb.com/

The source study used local PostgreSQL databases named `mimiciv31` and `eicu`, and SICdb CSV archives including `cases.csv.gz`, `data_float_h.csv.gz`, `data_range.csv.gz`, and `d_references.csv.gz`. Local database names do not independently authenticate releases. Verify source versions, schemas, code dictionaries and archive provenance before re-extraction. Source-hour aggregation and outcome ascertainment differ between databases; preserve the source-specific definitions in the code.

Extraction queries use read-only database sessions. Generated CSV, pickle and NumPy files may contain patient or stay identifiers, individual clinical measurements and dates. They must remain in a governed local environment and must not be committed. `.gitignore` is a safeguard, not a substitute for access control or manual review. Never use force-add to publish excluded outputs. Keep database passwords in the process environment or an approved local secret store.
