# Viral-context IgV-like candidate follow-up — 2026-11-16

This is follow-up within [issue #1](https://github.com/timodonnell/bio-fishing-expeditions/issues/1), prompted by a request for the gene/allele match, an illustration, and novelty assessment.

## Planned work

1. Recover the exact 360-aa candidate, its existing profile and reviewed-protein matches, and the 27 translated assembly features. Separate the receptor-family HMM from the best individual sequence match.
2. Inspect the primary database record for that antibody match. Determine whether it identifies a germline gene/allele, a historical antibody sequence, or something else; do not infer a donor allele from a short distant match.
3. Recover the pending public-database comparison if still available. Check the exact viral-context candidate separately if the earlier query was a related sequence. Record failed/expired requests and database/version limitations.
4. Generate a reproducible exact pairwise alignment against the recorded reviewed-protein hit, preserving parameters and scores. Illustrate the assembly neighborhood, matched portion of the protein, and residue-level alignment using scientific plotting software.
5. Search primary literature and exact identifiers for prior reports. Distinguish public availability of the sequence, novelty of an Ig-like viral protein, and evidence for a newly recognized antibody/TCR acquisition. A failure to find a paper is not proof of novelty.
6. Publish the figure, evidence tables, conclusions, and limitations in this experiment directory and update the existing issue if the interpretation changes.

## Evidence criteria and limits

Viral flanking homology can support viral context but cannot establish antibody ancestry. A local match to one antibody does not identify its germline allele, rearrangement, CDR specificity, or function. Catalog classification and assembled junctions require read-level validation before a strong natural-capture claim. Domain architecture, competing Ig-superfamily ancestry, and phylogeny remain important unresolved checks if this follow-up cannot settle them.

## Executed work

The exact reviewed-protein alignment and reference annotations were verified, the genome/protein/alignment figure was generated, and primary literature and exact identifiers were checked. See [IGV_CANDIDATE.md](IGV_CANDIDATE.md) and its linked manifests. No donor gene/allele or novel natural capture was established. The older related-sequence nr request returned UNKNOWN; a new exact-candidate request remains pending and is excluded from completed evidence.
