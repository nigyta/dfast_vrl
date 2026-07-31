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
    (`examples/cox1/multi.fa` + `common_example.json` + `specific_example.tsv`),
    asserting the MSS `COMMON` block, the per-entry `source`/`CDS`/`DBLINK`
    qualifiers and the report JSON.
  - Metadata is in the ddbj_mss_tools `batch_wgs_builder` format: a common JSON
    plus a two-header-row TSV keyed by `_entry` or `_file_path`. Both row keys
    are covered, as are the refusals (mandatory columns, entry mismatch,
    `COMMENT` in the TSV).

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
pytest tests/test_cox1.py
```

### Inside the dev container (docker compose)

`docker-compose.yml` mounts the repo at `/dfast_vrl` and the models at
`/vadr_models`, so a running `app` service can execute the suite as-is. This is
the only environment where every model is present and nothing is skipped.

**pytest is not installed in the image.** Install it once per container — the
image ships a PEP 668 "externally managed" Python, so plain `pip install` is
refused and `--break-system-packages` is required:

```bash
docker compose exec app pip install --break-system-packages pytest
docker compose exec app bash -lc 'cd /dfast_vrl && pytest tests'
```

The install does not survive `docker compose down` / a rebuild, so repeat it
after recreating the container.

### Inside the published image

```bash
docker run -it --rm \
    -v $PWD/vadr_models:/vadr_models \
    -v $PWD:/dfast_vrl \
    nigyta/dfast_vrl:latest \
    bash -lc "pip install --break-system-packages pytest && pytest /dfast_vrl/tests"
```

## Runtime notes

- Per-test timeout is 30 minutes. The whole suite takes ~2 min in the dev
  container (23 tests, all models present); mpox (~200 kb) and Flu dominate.
  Allow considerably longer on a slower host or under emulation.
- Each test writes into pytest's `tmp_path`; nothing is left in the repo.
- The tests treat the original `examples/...` and `/data/flu/...` files as
  read-only (they are passed by path; the scripts copy metadata into the
  output dir before mutating it).
