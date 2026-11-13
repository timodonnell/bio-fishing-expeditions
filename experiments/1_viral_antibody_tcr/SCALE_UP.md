# MetaVR scale-up plan — 2026-11-13

This extends [issue #1](https://github.com/timodonnell/bio-fishing-expeditions/issues/1). The previous milestone searched 3,207,211 protein records and 48,286 genome records, with overlap. This stage expands environmental-virus coverage; it does not change the evidence standard for natural antibody/TCR capture.

## Planned work

1. Freeze the public [MetaVR protein archive](https://portal.nersc.gov/cfs/m342/METAVR/METAVR.faa.bgz) by URL, response headers, byte count, retrieval time, and final SHA-256. The observed archive is 49,996,313,962 bytes, last modified 2026-08-26, ETag `"ba403356a-659fa4fc0dec0"`. Retrieve all protein records, not only species representatives or selected annotations.
2. Use resumable, checksum-verified HTTP ranges written directly into one partial file. Revalidate completed ranges on restart and reject changed source identity. Keep large data ignored by Git. Increase the stage's local cache allowance to 75 GiB while preserving at least 50 GiB free disk; use at most eight download connections and 16 search threads. Current free disk is about 159 GiB. Do not download the 72 GB nucleotide archive or the 115 GB prebuilt DIAMOND index in parallel with this stage.
3. Stream the BGZF FASTA into bounded shards and search all proteins with the existing 14 host-trained IG/TR HMMs. Save restart offsets, shard sequence/residue counts, hashes, parameters, software versions, candidate sequences, and alignments. Avoid a whole-file decompressed FASTA or permanent search index. Record observed protein and distinct source-genome identifier counts; advertised catalog counts are not executed coverage.
4. Make E-values comparable across shards: retain permissive hits using a fixed reference search-space size, then rescale independent-domain and sequence E-values to the final observed number of targets. Apply the existing main threshold (independent-domain E <= 0.001, aligned length >= 50 aa); retain weaker hits through E <= 0.01 separately. A checkpoint alone is not a completed search.
5. Verify shard/resume equivalence against a small unsharded control and rerun the held-out/known-FeLV calibration at the final search-space size. Preserve the previously documented calibration gaps. Confirm the natural TCR control remains detectable without adding it to training.
6. Compare retained candidates with the reviewed UniProt database and the receptor references, recover MetaVR genome/source metadata, and assess competing ancestry and contamination. Profile similarity alone cannot establish an antibody or TCR acquisition. Report unresolved candidates honestly.
7. Publish exact coverage, run manifests, candidate evidence, and limitations to `main`; update the existing issue. A second, independent whole-corpus protein-homology route and translated nucleotide screening require their own capacity/recall checks and remain subsequent work. Candidate-only BLAST is not a second independent corpus screen.

## Interpretation and stopping rules

This is an exploratory protein-profile expansion. It cannot detect an uncalled ORF, an intron-disrupted or frameshifted capture, or sequences excluded by catalog assembly/virus classification. Generic viral Ig folds, host contamination, and engineered inserts remain competing explanations. New natural captures require genomic context and independent support under the original plan. A negative result only bounds this method on the actual searched release.

Stop on source identity changes, checksum/parser failures, resource limits, or failed positive-control/equivalence checks. Save the last valid checkpoint and diagnose before continuing. Do not silently skip corrupt records. Interrupted downloads or searches must be reported as incomplete.

## Executed so far

- Rechecked the archive's HTTP headers and retrieved an 8 MiB range; the server returned the exact requested range with the same ETag. BGZF content contains protein FASTA records with source UViG identifiers in their headers.
- No full MetaVR search has completed at plan creation. Results and actual coverage will be recorded separately.
