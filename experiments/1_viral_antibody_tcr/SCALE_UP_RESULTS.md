# MetaVR protein expansion — completed 2026-11-13

We searched **369,319,989 protein records** from the complete acquired MetaVR protein archive: **115.15 times** the preceding RefSeq/INPHARED protein-record total. The search retained 26 candidate records. Two match host immunoglobulin loci, 14 have competing MHC-family homology, one matches CD300, and **nine environmental IgV-like records remain unresolved**. One unresolved candidate has substantial evidence of viral genomic context, but antibody/TCR ancestry has not been established.

**No new supported natural antibody or TCR acquisition was established.** The previously documented natural FeLV TCR beta-chain acquisition remains the positive answer for TCRs; see the [evidence packet](known_natural_tcr.md). Natural antibody acquisition remains an open question. This is a completed protein-profile expansion under [issue #1](https://github.com/timodonnell/bio-fishing-expeditions/issues/1), not completion of every search in the [original plan](PLAN.md).

## Actual coverage

| Measurement | Executed value |
| --- | ---: |
| Protein records | 369,319,989 |
| Amino-acid residues | 82,096,374,322 |
| Distinct source-genome identifiers in protein headers | 24,429,973 |
| Compressed archive bytes | 49,996,313,962 |
| Uncompressed FASTA bytes | 110,257,106,776 |
| Completed contiguous shards | 411 |
| Final retained domain rows | 32 |
| Retained protein records | 26: 15 main, 11 rescue |
| Candidate assemblies retrieved and checked | 21 |
| Candidate proteins exactly mapped to translated assembly sequence | 26/26 |

The source is the [NERSC MetaVR protein BGZF](https://portal.nersc.gov/cfs/m342/METAVR/METAVR.faa.bgz), retrieved on 2026-11-13, with Last-Modified **2026-08-26 21:57:55 GMT**, ETag `"ba403356a-659fa4fc0dec0"`, and SHA-256 **`30ac6f74eceae3490ddd26e169d63efe66c1ae410775890ab029653afad49a8d`**. The [download manifest](manifests/metavr_download.json) records every range and its checksum. The [scan manifest](results/metavr/scan.json) and [411 checkpoints](results/metavr/shards/) establish contiguous coverage through EOF.

Counts refer to records and header identifiers, not independent viral populations. No full cross-catalog deduplication was performed. Eight retained protein records exactly match a sequence in the earlier RefSeq/INPHARED inputs; this is **candidate-only overlap accounting**. The protein-header count is not a count of genomes searched by translated DNA, and is not inferred from the catalog website's advertised total.

The archive contains repeated identifiers. An initial run stopped after eight shards; [that attempt](results/metavr_attempt1/) is preserved and is not additional coverage. The completed scanner gives every input record a stable ordinal suffix, preserving all occurrences. It observed 38,048 duplicate identifiers **within shards**, which is not a global duplicate count.

## Methods and checks

The 14 frozen human/mouse-trained IG/TR profiles were unchanged (HMM SHA-256 `e64d4c4729d58529b3cd993d6a2b79d8fac539e16e63ec03221288b9f0339893`). PyHMMER 0.12.3 searched 256 MiB uncompressed shards with four workers and four threads per worker. Scores used reference `Z=1,000,000`, followed by sequence and independent-domain E-value rescaling to the observed **369,319,989** targets. Main candidates require independent-domain E ≤ 0.001 and ≥50 aligned amino acids; rescue candidates extend through E ≤ 0.01. These are retrieval thresholds, not a false-discovery rate across 14 models.

The final scan took 1,212 seconds after the initial archive hash; this is an observed local runtime, not a general benchmark. The parent process's peak RSS was 3.06 GiB, excluding worker RSS, which is recorded separately in checkpoints. Bulk data and temporary indexes remain outside Git.

Candidate BLASTP used the frozen reviewed UniProt **2026_04** database (576,073 proteins) and BLAST+ **2.12.0**. All nine environmental candidates were also compared against the downloaded RefSeq viral protein set, with no reported match at E ≤ 0.01. These candidate-specific comparisons are not independent whole-MetaVR searches. The full commands, inputs, output hashes, versions, and parameters are in the corresponding `*.run.json` files.

Candidate assemblies were obtained using the pinned nucleotide `.ffx` index and byte ranges, transferring about 24 MB of retained ranges rather than downloading the full nucleotide archive. Every candidate maps exactly to a translated assembly sequence using its deposited genetic code, allowing initiating methionine at a valid alternative start codon. This establishes assembly linkage; it does not independently validate the assembly or prove natural viral incorporation.

Integrity evidence:

- [Bulk validation](results/scale_integrity_validation.json): **1,241 hash checks**, contiguous archive coverage, global E-value rescaling, complete retained FASTA, and pipeline controls passed.
- [Post-scan review validation](results/scale_review_validation.json): **146 hash checks** across eight completed search runs; all 26 candidate mappings, all 27 proteins around the viral-context lead, host-reference inputs, and assessment counts checked.
- [Core validation](results/validation.json): **129 hash checks** and 19 earlier search runs passed, including the newly acquired small reference inputs.
- [Pipeline test](results/scale_pipeline_validation.json): sharded/unsharded equivalence, large-block-crossing records, restart preservation, damaged-checkpoint rejection, duplicate-identifier preservation, and known-positive retrieval passed. [Download interruption/source-change tests](results/scale_download_validation.json) also passed.

## Calibration at the full database size

The complete control set was rerun directly at the observed final `Z`; prior near-identity train/holdout exclusions were retained.

| Control class | Retrieved / evaluated at main threshold |
| --- | ---: |
| Held-out receptor domains | 79/79 |
| Held-out receptor proteins | 72/72 |
| 90-aa fragments | 77/77 |
| 60-aa fragments | 74/79 |
| 30-aa fragments | 0/77 |
| Hard negatives | 8/231 |
| Shuffled controls | 0/83 |
| Known natural FeLV TCR positive | 1/1 |
| Camelid VHH structural-construct spot check | 1/1 |
| Shark VNAR structural-construct spot check | 1/1 |

The 30-aa class cannot meet the 50-aa main alignment gate; its zero is a structural exclusion, not an estimate of short-fragment recognition. Hard-negative recovery confirms that a profile match is insufficient for ancestry. The added controls are deposited constructs **PDB 1MVF entity 1** and **1T6V entity 2**, frozen in [references](references/scale_antibody_controls.faa), with [source and exclusion details](results/scale_antibody_controls.json). Neither is a natural viral carrier. Notably, the shark VNAR is recovered by a **TRA_V** profile: profile names must not be interpreted as ancestry labels.

These spot checks do not complete the planned engineered-carrier, broader VHH/VNAR, divergence/family-holdout, or splice/frameshift controls. [Full calibration table](results/metavr/calibration_targets.tsv).

## Candidate assessment

| Assessment | Main | Rescue | Total |
| --- | ---: | ---: | ---: |
| Host immunoglobulin loci | 2 | 0 | 2 |
| Competing MHC-family homology | 10 | 4 | 14 |
| Competing CD300-family homology | 0 | 1 | 1 |
| Unresolved environmental IgV-like | 3 | 6 | 9 |

All identifiers, profile coordinates, scores, genome metadata, source samples, competing matches, and exact prior-corpus matches are in [candidate_assessments.tsv](results/metavr/candidate_assessments.tsv). These are evidence-based triage categories; MHC and CD300 are distinct from antibody/TCR chains, even when found in a virus-associated record.

### Host immunoglobulin matches

**Mouse: `IMGVR_UViG_3300051405_000263`.** The entire **18,461 bp** assembly is **100% identical** to the reverse complement of **GRCm39 chr16:19,008,356–19,026,816**. Direct sequence equality was checked independently of the BLAST alignment. The candidate overlaps **Iglc4 / ENSMUSG00000106039**, an immunoglobulin lambda constant gene. Source metadata identifies a mouse fecal microbiome. MetaVR labels the assembly low viral confidence with zero hallmark genes, despite a Retroviridae prediction. This supports host sequence misclassification, not viral capture. [Evidence](results/metavr/host_comparison/mouse_identity.json).

**Human: `IMGVR_UViG_3300032031_000646`.** All **10,225 query bases** are covered by two **99.507% and 99.413% identity** alignments to the hg38 immunoglobulin heavy-chain region. The alignments involve duplicated blocks at chr14:106,746,132–106,752,220 and 106,706,423–106,710,678, so this is not one exact contiguous reference copy. The source is a human scalp metagenome; viral confidence is low with zero hallmark genes. Host IGH sequence/assembly is strongly supported, with no independent viral flanks demonstrated. [Evidence](results/metavr/host_comparison/human_identity.json).

### Environmental lead with viral context

The nine unresolved records represent **seven exact distinct protein sequences** from soil metagenomes at Hubbard Brook and Harvard Forest. Similarity to reviewed antibodies is limited to a variable-like domain; this does not establish an antibody-derived coding sequence.

The strongest viral-context lead is **`IMGVR_UViG_3300079393_000044`**, a **17,620 bp** assembly from a Harvard Forest NEON soil sample. Its 360-aa candidate, source identifier **`Ga0643767_0002625_14604_15686`**, matches IGK_V at global independent-domain **E = 3.95 × 10⁻⁴**. The translated coding sequence maps to bases **14,604–15,683** on the positive strand (the deposited identifier includes the stop codon). MetaVR predicts Iridoviridae/Ranavirus, high viral confidence, one hallmark, and approximately 9.15% genome completeness.

We retrieved all **27 deposited proteins** and verified their translation against the same assembly. Nearby genes independently support viral context:

| Neighbor's deposited CDS interval | Reviewed viral homolog | BLASTP E-value |
| --- | --- | ---: |
| 14,195–14,554 | IIV3 transcription elongation factor S-II-like protein, Q197A5 | 3.73 × 10⁻²⁰ |
| 15,715–16,191 | Frog virus 3 protein 045L, Q6GZT1 | 8.72 × 10⁻⁹ |
| 17,133–17,618 | IIV6 Bro-N domain protein 289L, Q91FN5 | 2.89 × 10⁻⁷⁰ |

Additional genes elsewhere on the assembly match iridovirus Bro-N and phosphatase proteins. [All neighbor matches](results/metavr/soil_neighbors/gene_review.tsv), [translations](results/metavr/soil_neighbors/locations.json), and [summary](results/metavr/soil_neighbors/summary.json) are retained. These observations support a viral assembly containing an IgV-like protein. They **do not yet distinguish antibody/TCR ancestry from another Ig-superfamily protein**. Ig-like domains were already reported in the [Singapore grouper iridovirus genome study](https://doi.org/10.1128/JVI.78.22.12576-12590.2004), so viral context plus an Ig-like domain cannot by itself be called a new receptor capture.

The other eight environmental assemblies have low catalog viral confidence and no hallmark genes; they remain unresolved rather than being automatically classified as host contamination. Shared sequence across assemblies is not equivalent to independent read support.

## Pending work and limits

- Establish the environmental family's ancestry with broader homologs, competing Ig-superfamily references, domain architecture, and phylogeny; validate candidate junctions using source reads where available.
- Complete an independent whole-MetaVR protein-homology route and translated-nucleotide search. The present expansion searched predicted proteins only; missing ORFs, introns, frameshifts, and catalog-excluded sequences can be missed.
- Complete the remaining calibration controls and whole-GenBank/unclassified-sequence work. Previously reported annotation-only leads remain open.
- Deduplicate larger catalogs before claiming independent added sequence coverage.

Two optional NCBI remote comparisons were submitted: **FF0V9BZ1016** (mouse assembly against core_nt) and **FF1AP9JS014** (one environmental candidate against nr). They were still queued at the last recorded poll and are **not counted as completed searches or evidence**. The mouse conclusion instead rests on the completed local reference comparison and exact sequence equality. Polling commands and request parameters are retained under [host_comparison](results/metavr/host_comparison/remote_blast.json) and [soil_nr_comparison](results/metavr/soil_nr_comparison/remote_blast.json); remote RIDs may expire.

The protein expansion is complete and reproducible from the recorded inputs. Its negative findings do not establish that natural antibody acquisition has never occurred.
