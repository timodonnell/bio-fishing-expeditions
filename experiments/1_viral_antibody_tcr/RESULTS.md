# Core RefSeq/INPHARED results — 2026-11-13

**Subsequent milestone:** the [completed MetaVR expansion](SCALE_UP_RESULTS.md) searched 369,319,989 protein records and retained nine unresolved environmental IgV-like records after candidate triage. The counts and pending-work statements below describe the earlier core milestone; use the scale-up report for current MetaVR status.

**Yes: natural incorporation of a TCR beta-chain gene into feline leukemia virus is documented.** The original report is [Fulton et al., 1987](https://pubmed.ncbi.nlm.nih.gov/3029597/); a [1992 follow-up](https://pmc.ncbi.nlm.nih.gov/articles/PMC241135/) supports its natural tumor-associated provenance and packaging. We recovered the deposited locus **X05155.1 / CAA28801.1 / UniProt P11364**, verified its translation, and retrieved it using host-only references. See the [evidence packet](known_natural_tcr.md).

This is a known acquisition of one receptor chain, not a new discovery or a complete functional receptor. **Natural antibody capture remains unresolved.** No new supported natural acquisition was established by the screens below. Weak and annotation-only candidates remain open; this is an interim result, not a completed search of all viral sequence databases.

## Executed searches

| Corpus | Genome records searched by TBLASTN | Bases searched | Protein records searched by HMM and BLASTP | Residues searched |
| --- | ---: | ---: | ---: | ---: |
| RefSeq viral FTP snapshot | 19,621 | 583,706,107 | 722,740 | 175,002,185 |
| INPHARED 14Apr2025, excluding RefSeq identifiers | 28,665 | 1,713,788,886 | 2,484,471 | 517,146,495 |
| Total records, including cross-corpus overlap | **48,286** | **2,297,494,993** | **3,207,211** | **692,148,680** |

The RefSeq release marker read **237**; the acquired viral files had directory timestamps of 2026-09-03. Their hashes, rather than a mutable current-release URL alone, define the snapshot. Four RefSeq GenBank records had undefined sequence and were excluded from nucleotide searches. All protein files listed in the acquired viral directory were searched.

INPHARED's “excluding RefSeq” designation does not eliminate sequence overlap. We measured **5,691 genome records** and **1,204,734 protein records** identical to a downloaded RefSeq sequence. INPHARED contributed **21,528 distinct genome sequences** and **764,276 distinct protein sequences** not present there. Exact genome comparisons did not normalize rotation or reverse complement, so these are sequence-content counts, not independent viral populations. INPHARED accessions omit versions; the archived input hash and release identify the actual sequences searched.

A preceding pilot searched 10,000 genomes, 306,564,279 bases, and 340,419 annotated proteins. It is a subset of RefSeq and is **not additional coverage**. Selection used 614 family/molecule/length/annotation strata and seed 20261112. Host qualifiers were retained but host-taxonomic stratification was not implemented.

Separately, two Entrez Protein queries selected **960 records** under viral taxonomy with receptor/antibody terminology. These were profile-screened and inspected for provenance. This is an **annotation-selected survey**, not a whole-GenBank homology screen. Queries, record accessions, and download hashes are retained.

## Methods and calibration

The reference survey yielded 534 reviewed host IG/TR proteins from 17 species. Human/mouse references supplied 460 protein queries and 14 family/domain HMMs; other species were held out. P11364 was excluded from training. Frozen small inputs are in [references/](references/); source records are from UniProt **2026_04**. Profile construction used MUSCLE 3.8.31 and PyHMMER 0.12.3; sequence searches used BLAST+ 2.12.0. The Python environment is locked with package hashes.

Retrieval retained sequence E-values or independent-domain E-values up to 0.01. The main queue used E ≤ 0.001 and at least 50 aligned amino acids. SEG was enabled for BLAST protein queries. TBLASTN searched both strands independently of the protein results, using the standard genetic code. Raw E-values are corpus-specific retrieval statistics, not a pooled significance test across profiles/corpora. Calibration HMM E-values were rescaled to 3,207,211 target records; BLAST calibration used a conservative effective database of one billion residues. Candidates were retained permissively for review rather than requiring agreement among methods.

After excluding near-identical train/holdout matches (≥90% identity spanning ≥80% of the target), retrieval was:

| Control | HMM | BLASTP | Union |
| --- | ---: | ---: | ---: |
| Held-out intact domains | 79/79 | 79/79 | 79/79 |
| Held-out protein records | 72/72 | 72/72 | 72/72 |
| 60-aa fragments | 78/79 | 74/79 | 78/79 |
| 90-aa fragments | 77/77 | 77/77 | 77/77 |
| Known natural FeLV TCR | 1/1 | 1/1 | 1/1 |
| Hard negatives retrieved | 20/231 | 61/231 | 71/231 |
| Shuffled domains retrieved | 0/83 | 0/83 | 0/83 |

Thirty-amino-acid fragments cannot pass the main queue's 50-aa requirement; their zero main-queue retrieval is not a measurement of short-fragment sensitivity. Related domains and species are not independent evolutionary trials. The benchmark establishes useful retrieval for this panel, not universal sensitivity or specificity. The hard-negative matches demonstrate that receptor-like scores alone cannot establish receptor ancestry.

**The full calibration gate in the original plan is not complete.** Verified engineered-carrier controls, dedicated VHH/IgNAR panels, family-holdout tests, simulated divergence, and systematic frameshift/splice controls remain unfinished. The core extension was run as an explicitly exploratory screen; no broad negative-result sensitivity claim is made. The known natural FeLV case served as an independently withheld positive, including a reverse-complement nucleotide control.

## Candidate review

The [route counts](results/search_route_counts.tsv) retain both main and rescue hits. Full RefSeq searches yielded 17 main HMM protein targets, 14 main BLASTP protein targets, and 14 main TBLASTN genome targets; these overlap and must not be added as independent candidates. Including pilot and rescue results, **59 protein candidates** and **51 distinct translated intervals** were reviewed against 576,073 reviewed UniProt proteins and the deposited annotations.

Many candidates have stronger support for MHC-family proteins, CD200/OX-2 homologs, cytokine-binding proteins, viral surface proteins, polymerases, or phage nucleases. For example, the strongest RefSeq HMM hit, YP_004207833.1, is annotated as MHC-like and matches a reviewed MHC protein much more strongly than the receptor queries. The translated search also recovered MHC-related sequence at a locus annotated as truncated by a frameshift, illustrating why nucleotide screening adds information. Assessments and residual uncertainty are in [protein candidates](results/refseq_protein_candidates.tsv) and [translated candidates](results/refseq_translated_candidates.tsv); they are triage assessments, not complete evolutionary reconstructions.

INPHARED produced **no main-queue protein or translated-DNA hits** with these references/settings. One weak BLASTP rescue match, GU071097_00087, is retained at E=0.005 over 50 aa, without supporting profile or reviewed-protein matches. Its source is Synechococcus phage S-SSM5, and its genomic location is recorded in [the rescue record](results/inphared_rescue_candidate.json). This does not establish antibody ancestry.

The annotation-selected survey makes the provenance problem concrete. Its 42 reviewed candidates comprise three representations of the **same FeLV precedent**, 20 protein-structure components without viral genomic evidence, 14 older protein-only records with unresolved taxonomic provenance, and five annotation-only genomic leads. An antibody chain in a virus–antibody structure is not evidence that the viral genome encodes it.

| Annotation-only lead | Deposited locus | Protein length | Current assessment |
| --- | --- | ---: | --- |
| DAJ54870.1 | BK043519.1 | 46 aa | Phage metagenomic gene with antibody-like annotation based on a PDB motif |
| DAW92981.1 | BK025828.1 | 87 aa | Phage metagenomic gene with the same annotation source |
| DBA58143.1 | BK068258.1 | 460 aa | Pithovirus gene annotated as TCR-like |
| XOK31829.1 | PQ067554.1 | 114 aa | Chlorovirus gene annotated as TCR beta |
| XYW61038.1 | PV288767.1 | 65 aa | Chlorovirus gene annotated as TCR gamma-like |

All five source CDS translations were verified against their protein records. None produced a retained receptor-profile match or a reviewed-protein BLASTP match at E ≤ 0.01 in these searches. Short length and divergence limit what that means. They remain **unresolved annotation-only leads**, not validated acquisitions and not proven artifacts. Their exact locations, source qualifiers, read/study links where present, and hashes are in [annotation_only_loci.json](results/annotation_only_loci.json). Raw-read support and evolutionary ancestry have not been assessed.

## Coverage gaps and next work

MetaVR bulk access works. HEAD responses measured its protein BGZF at **49,996,313,962 bytes**, nucleotide BGZF at **71,552,362,108 bytes**, and DIAMOND index at **114,983,287,451 bytes**. No MetaVR sequences were searched. A resumable streaming/sharded search is the next large-corpus stage; the present cache stayed below the initial 20 GiB limit and retained the required disk headroom.

Whole GenBank homology screening, MetaVR, supplementary gut/marine catalogs, unclassified-metagenome rescue, germline nucleotide BLASTN, alternative genetic codes, read-based junction validation, and phylogenetic adjudication remain pending. The original planning search did not find the FeLV precedent; the reference survey did. This is why neither a preliminary literature search nor absence from RefSeq supports “never.”

Issue #1 remains open. The immediate priorities are to strengthen the reference/control panel, assess the five annotation-only loci, and execute the broad metagenomic stage with explicit coverage accounting. The current defensible answer is **yes for natural TCR beta-chain capture; natural antibody capture unresolved**. Engineered antibody/TCR precedents remain documented separately in the [literature preflight](literature_preflight.md).

## Reproducibility

See [RUNNING.md](RUNNING.md), [download manifests](manifests/downloads.json), [run manifests and outputs](results/), and [validation results](results/validation.json). Each search records commands, parameters, software/binary versions, input/output hashes, elapsed time, and peak process memory. Bulk databases and indexes are ignored by Git. No new experimental biological work was performed.
