# Reproducing the current milestone

Run commands from the repository root. Python 3.10.13 was used. All bulk input, indexes, and temporary files belong under ignored `data/raw/`, `data/cache/`, or `work/`.

```bash
EXP=experiments/1_viral_antibody_tcr
python -m venv "$EXP/work/venv"
"$EXP/work/venv/bin/pip" install --require-hashes -r "$EXP/requirements.lock"
```

The original environment was populated with `uv pip`, using the same lock. The local BLAST+ 2.12.0 binaries are at `/home/bizon/anaconda3/envs/imprint/bin`; MUSCLE 3.8.31 is at `/home/bizon/git/ANARCI/bin/muscle`. Their paths and hashes are in the manifests. On another host, install those versions and adjust `BLAST_DIR` in `scripts/search.py` and `MUSCLE` in `scripts/build_references.py`. Do not silently substitute a newer binary when reproducing a result.

## Inputs

`config/acquisition_pinned.json` contains successful download URLs and expected SHA-256 hashes. Original request groups preserve failed attempts and the order of discovery. Public current-release endpoints can change: a future mismatch must be resolved by locating the recorded release, not by accepting different data under the same result label.

```bash
python "$EXP/scripts/acquire.py" "$EXP/config/acquisition_pinned.json"
```

The INPHARED host was slow for single-stream downloads. `acquire_ranges.py` used eight workers, 4 MiB ranges, `If-Match`, exact `Content-Range` checks, and a final checksum. Initial single-stream work was resumed only against the observed ETag. No full MetaVR download was started.

## Preparation and searches

Preparation dependencies are sequential. Do not search a FASTA until its preparation process has exited successfully.

1. Run `prepare_corpus.py` to produce the pilot and the full downloaded RefSeq protein FASTA.
2. Use `build_references.py` and `prepare_benchmark.py`, or restore the exact small reference inputs from `references/` into `data/cache/references/`. The frozen HMM avoids timestamp-dependent changes to regenerated HMM files. Complete the benchmark searches before `summarize_calibration.py`; the final calibration also reads the INPHARED protein counts.
3. Run the pilot HMM, BLASTP, and independent TBLASTN searches, followed by the full RefSeq protein searches.
4. Run `prepare_expansion.py`, then full RefSeq TBLASTN. Run `prepare_translated_candidates.py` to recover translated intervals and deposited features.
5. Run `prepare_inphared.py protein` and `prepare_inphared.py genome` after the corresponding downloads finish, then HMM, BLASTP, and independent TBLASTN.
6. `prepare_auxiliary.py swissprot` produces the reviewed comparison database; `prepare_auxiliary.py annotation` prepares the annotation-selected records and, after their HMM search, the candidate FASTA. `prepare_auxiliary.py phage_candidates` extracts the phage rescue candidate after BLASTP. Compare candidate proteins/translated intervals with the reviewed database.
7. Run `summarize_calibration.py`, `summarize_candidates.py`, and `validate_results.py`.

For example:

```bash
"$EXP/work/venv/bin/python" "$EXP/scripts/search.py" hmm \
  "$EXP/data/cache/references/receptor_profiles.hmm" \
  "$EXP/data/cache/refseq/all_proteins.faa" \
  "$EXP/results/refseq_all_hmm.tsv" --cpus 8
```

Every executed search has its **exact command**, input hashes, binary/version information, parameters, output hash, start/end timestamps, elapsed time, and maximum process RSS in the matching `results/*.run.json`. These commands cover benchmark, pilot, full-corpus, and candidate-comparison searches. Keep total search threads at or below 16. Runtime measurements are local observations, not a capacity estimate for MetaVR.

The script SHA in each search manifest is the source used for that search; `scripts/archive/search_initial.py` preserves the exact initial version before the optional HMM/BLAST effective database-size parameters were added. Its hash matches the three early benchmark manifests; its Git commit is the repository HEAD at invocation, before this milestone was committed. Preparation commands are identified above, with input/output hashes in the preparation summaries. `prepare_auxiliary.py` captures extraction steps initially executed interactively. It is not a scheduler or a one-command implementation of all future stages in the original plan.

## Validation and known limits

`validate_results.py` rehashes completed downloads and all search inputs/outputs, verifies frozen reference copies, and checks the known-positive and five annotation-only CDS validations. The positive-validation script checks the actual FeLV CDS translation. Full search target hashes were checked against completed preparation outputs.

The planned full calibration gate remains incomplete. Alternative genetic codes, BLASTN germline screening, dedicated splice/frameshift rescue, raw-read analysis, and phylogeny are future work. Unresolved candidates and incomplete corpora are recorded in `RESULTS.md`; no script automatically calls a novel natural acquisition.

## MetaVR scale-up

The [scale-up plan](SCALE_UP.md) extends the same issue. Install the optional index-reader dependency with the additional `requirements-scale.lock`. The original core environment lock remains unchanged.

```bash
EXP=experiments/1_viral_antibody_tcr
"$EXP/work/venv/bin/python" "$EXP/scripts/acquire_large.py" "$EXP/config/metavr_scale.json"
"$EXP/work/venv/bin/python" "$EXP/scripts/scan_bgzf.py" \
  "$EXP/data/raw/metavr/METAVR.faa.bgz" "$EXP/results/metavr" \
  --workers 4 --cpus 4 --shard-mib 256
```

Run acquisition to completion before scanning. Repeating either command validates its existing inputs/checkpoints and resumes incomplete work. The downloader writes HTTP ranges into one partial file and checks every completed range on restart; its manifest is `manifests/metavr_download.json`. It reserves 50 GiB of free disk. Search checkpoints contain contiguous BGZF virtual offsets, uncompressed shard hashes, counts, and output hashes. Four workers each use four search threads. The input archive is hashed before a search or restart. Candidate sequences and alignments are small retained outputs; bulk archive and genome-ID coverage files stay ignored.

Shard E-values use reference `Z=1,000,000`; final sequence and independent-domain E-values are multiplied by the actual total protein count divided by that reference. Main candidates require independent-domain E <= 0.001 and >=50 aligned residues; rescue candidates extend to E <= 0.01. These are per-profile thresholds, not an FDR estimate across the 14 models. The relationship between `Z` and E-values follows the [PyHMMER significance documentation](https://pyhmmer.readthedocs.io/en/v0.12.1/guide/significance.html). Only a `complete` scan manifest establishes full protein coverage.

`test_scale_pipeline.py` compares a multi-shard resumed scan against an unsharded scan and the original scanner at `Z=1,000,000,000`, tests a record crossing BGZF blocks, rejects a damaged checkpoint, and checks the natural FeLV positive. `test_large_download.py` injects a download interruption, corrupts a previously completed range, and checks exact recovery and source-change rejection. Their validation reports are under `results/`.

`metavr_context.py <candidate.faa> <output-directory>` retrieves candidate assemblies through the pinned nucleotide `.ffx` index and HTTP byte ranges, without acquiring the full nucleotide archive. It validates index layout, compressed-block checksums, sequence length/alphabet, and agreement with MetaVR metadata. Candidate-to-genome mapping searches for the exact translated protein using the deposited genetic code, allowing the usual initiating-methionine rule. Failure to map is retained as an unresolved provenance problem. This candidate-specific check is not a whole-corpus translated-DNA screen. The small `metavr_retrieval_control` artifact tests the route on an ordinary MetaVR protein, not on an antibody/TCR positive.
