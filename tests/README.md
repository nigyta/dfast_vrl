# DFAST_VRL / vadr2mss.py operational tests

End-to-end (subprocess) tests that drive the two CLI scripts at the repo root
and verify the produced files. They mirror the manual smoke tests we used to
run for each release.

## What is covered

- `test_dfast_vrl.py`
  - `--help`, `--version`, missing-input handling, `--force` semantics.
  - Full pipeline on a SARS-CoV-2 sample
    (`examples/nig_vrl_sample.fasta` + `examples/metadata.txt`).
- `test_vadr2mss.py`
  - `--help`, invalid `-M`, missing-input handling.
  - Full pipeline for **mpox**, **sarscov2**, **RSV**, **Flu (A)**, **Flu (B)**.
- `test_cox1.py`
  - `--help`, missing-input handling, `--force` semantics.
  - Full **COX1** pipeline on the two-record multi-FASTA
    (`examples/cox1/multi.fa` + `examples/cox1/metadata_example.tsv`), asserting
    the MSS `COMMON` block, the per-entry `source`/`CDS`/`DBLINK` qualifiers and
    the report JSON.
  - One `xfail` (`test_no_empty_qualifier_values`) pins an open defect: blank
    `##SPECIFIC` columns are emitted as empty-valued qualifiers. Drop the marker
    when that is fixed.

Tests that need a particular VADR model are skipped automatically when the
corresponding `*.minfo` file is missing under `$VADRMODELDIR` (default
`/vadr_models`). Tests are also skipped when the example fasta / metadata
files are missing — useful outside the container, where `/data/flu/...`
might not be mounted.

## Running

From the repo root:

```bash
pip install pytest        # if not already installed
pytest tests
```

Run a subset:

```bash
pytest tests/test_vadr2mss.py -k mpox
pytest tests/test_vadr2mss.py -k "flu_A or flu_B"
pytest tests/test_dfast_vrl.py::test_full_pipeline_scov2
```

Inside the Docker image:

```bash
docker run -it --rm \
    -v $PWD/vadr_models:/vadr_models \
    -v $PWD:/dfast_vrl \
    nigyta/dfast_vrl:latest \
    bash -lc "pip install pytest && pytest /dfast_vrl/tests"
```

## Runtime notes

- Per-test timeout is 30 minutes; mpox (~200 kb) and Flu can take several
  minutes each.
- Each test writes into pytest's `tmp_path`; nothing is left in the repo.
- The tests treat the original `examples/...` and `/data/flu/...` files as
  read-only (they are passed by path; the scripts copy metadata into the
  output dir before mutating it).
