# Research plan: antibody/TCR-derived sequence in viral genomes

Issue [#1](https://github.com/timodonnell/bio-fishing-expeditions/issues/1). Written 2026-11-12, before sequence screening. All numerical screening thresholds below are **pilot starting points**, not validated discovery criteria. Changes after calibration must be recorded before the broad search.

## 1. Define what would answer the question

Primary question: is there credible evidence that a naturally occurring viral genome acquired coding sequence specifically derived from a jawed-vertebrate antibody/immunoglobulin (IG) or T-cell receptor (TR) locus?

Include animal, plant, fungal, protist, bacterial, and archaeal viruses; DNA and RNA genomes; free viruses and proviruses. Do not restrict discovery to human viruses or to viruses already predicted to infect vertebrates. Record fragments and pseudogenized acquisitions as well as intact coding sequences. An acquired constant domain is a gene-fragment finding, not evidence for a complete antigen-binding antibody.

Report two independent dimensions:

| Dimension | Categories |
| --- | --- |
| Molecular identity | Complete antibody/TCR chain; variable domain; constant domain; germline V/J/C or other identifiable fragment; generic Ig-superfamily protein; unresolved |
| Genomic provenance | Supported natural acquisition; provisional natural candidate; engineered viral genome/vector; packaged phagemid or host nucleic acid; host/assembly contamination; protein association only; unresolved |

Antibody heavy and light chains and TCR alpha, beta, gamma, and delta chains are in scope. A full conventional antibody or TCR requires the appropriate chains and architecture; one acquired chain does not demonstrate a complete functional receptor. Single-domain antibodies are a separate architecture. Function and antigen specificity cannot be inferred solely from sequence incorporation.

The central competing explanations are: genuine IG/TR ancestry; another Ig-superfamily lineage with a similar fold; an engineered insert; host contamination; an incorrect viral boundary or chimeric assembly; or an annotation mistake. Proteins that bind antibodies, CD4/CD8 or other host receptor homologs, and RAG/transposon ancestry are not themselves antibody/TCR gene captures.

## 2. Establish the literature baseline

Search PubMed and Europe PMC using the planned query families in `search_queries.tsv`; supplement with publisher searches and backward/forward citation tracing. Include natural acquisition, gene capture, horizontal transfer, transduction, immunoglobulin/TCR homologs, phage display, and viral vectors. Search both IG and TR terminology and relevant viral families. A search for antibodies *against* viruses answers a different question.

For each relevant paper record DOI/PMID, year, viral taxon, claimed receptor type, natural/engineered provenance, genomic versus protein-only evidence, sequence accession if available, and exactly what was inspected (abstract, full text, supplement, sequence record). Retrieve full methods and deposited constructs before choosing engineered sequence controls. In phage-display papers distinguish a complete engineered phage genome from a packaged phagemid and helper-phage proteins.

Initial evidence is in `literature_preflight.md`. It establishes useful controls and false-positive categories, but is not a systematic literature review and does not establish that a natural example is absent. The final review must retain queries, search dates, result identifiers/counts, exclusions, and inaccessible full texts.

## 3. Freeze the sequence search space

Use the following order. Published catalog sizes describe scale, not the number searched in this experiment. Exact release, accessions, bytes, sequence counts, nucleotide/amino-acid totals, checksums, and retrieval date must be measured at acquisition.

| Priority | Corpus | Purpose and caveat |
| --- | --- | --- |
| A | NCBI viral RefSeq genomes and proteins | Curated, tractable starting corpus across viral groups. Include incomplete records in a labeled stratum. |
| A | NCBI GenBank viral nucleotide records and linked proteins | Expand beyond RefSeq and beyond annotated/complete genomes. Use Virus/Datasets plus Entrez/FTP exports as needed; reconcile endpoint coverage. Preserve taxid, accession.version, FEATURES, segment, BioProject/BioSample, publication, and read links. |
| A | INPHARED genomes, proteins, and metadata | Dedicated cultivated-phage coverage and useful provenance. The inspected README advertises a 14Apr2025 snapshot; do not invent a newer release or treat it as current without verification. Substantial overlap with GenBank is expected. |
| B | MetaVR | Broad environmental and uncultivated-virus discovery. The 2026 paper reports 24,435,662 UViGs; this is the successor to IMG/VR. Cover both viral-confidence tiers, with separate adjudication standards and counts. |
| B fallback | IMG/VR v4 | Use as an explicitly older fallback or to recover missing provenance. Its published collection exceeds 15 million genomes/fragments. Do not add overlapping MetaVR/IMG/VR counts as independent coverage. |
| B supplement | Gut Phage Database and Metagenomic Gut Virus catalog | Add gut-associated phages and original study/read links; published collections contain 142,809 and 189,680 sequences, respectively. Measure unique added content. |
| B supplement | Global Ocean Viromes 2.0 | Environmental breadth and marine phages; retain the distinction between all supplied fragments and the publication's 195,728 populations at least 10 kb long. |
| C | Unclassified metagenomic contigs/proteins, including accessible MGnify studies | Rescue sequences missed by viral taxonomic filters. Require protein-to-contig mapping and separate viral evidence. A metagenomic protein hit alone cannot establish viral residence. |
| Validation | NCBI nr/nt, host genomes/transcriptomes, original SRA/ENA studies | Test competing homology/provenance explanations and verify promising loci. These searches are not a substitute for the primary viral census. |

NCBI documents downloads of viral genomes with proteins, CDS, and metadata ([documentation](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/how-tos/virus/virus-download/)). The broader catalog choices are supported by the [MetaVR](https://pubmed.ncbi.nlm.nih.gov/41312645/), [IMG/VR](https://pmc.ncbi.nlm.nih.gov/articles/PMC9825611/), [GPD](https://pmc.ncbi.nlm.nih.gov/articles/PMC7895897/), [MGV](https://pmc.ncbi.nlm.nih.gov/articles/PMC8241571/), and [GOV2](https://pmc.ncbi.nlm.nih.gov/articles/PMC6525058/) papers.

Search all distinct sequences where practical. Exact-sequence deduplication may save computation, provided every accession, locus, and study remains recoverable. Do not search only one genome per species/vOTU: an insertion may occur in just one strain. Protein-cluster representatives can support an initial triage but cannot stand in for a complete member search without measured recall and an explicit coverage limitation. Identical proteins in different genomes still need separate context checks.

Preserve potentially engineered records for a separate analysis instead of deleting them by keyword. Explicitly track taxonomy-unassigned records, proviral segments of cellular assemblies, missing CDS, short fragments, and excluded records. Taxonomy and viral-confidence scores are evidence to inspect, not ground truth.

## 4. Build references and controls before discovery

Prepare a versioned panel from documented IMGT/NCBI IgBLAST germline references, curated public protein records, and selected full-length rearranged receptors. Record source terms and accession-level provenance; do not assume an entire repertoire collection can be redistributed.

- Cover IGH, IGK, IGL and TRA, TRB, TRG, TRD; include variable and constant domains separately and nucleotide V/D/J/C segments where available.
- Include representative mammals, birds, reptiles, amphibians, bony fish, and cartilaginous fish; include camelid VHH and shark IgNAR. Document gaps rather than assigning every hit to human/mouse by default.
- Keep rearranged receptors and germline segments distinct. V(D)J junctions can support identity but must not be required for acquisition of an unrearranged gene segment.
- Select experimentally documented engineered inserts as positive controls only after their actual nucleotide records and carrier type are verified. No accession has yet been selected as a validated positive control.
- Add hard negatives: phage Hoc and other bacterial/phage Ig-like proteins, viral IgSF homologs, host CD4/CD8/CD200, adhesion proteins, MHC-related proteins, and Fc-binding proteins. Include length/composition-matched shuffled controls as an additional statistical check, not the sole negative panel.
- Hold out whole receptor families or species from query/profile construction for calibration. Remove exact and near-identical leakage between training and holdout sets.

Build separate family/domain profiles rather than one generic Ig model. Existing broad Ig-domain annotations are a recall aid and a source of decoys. Keep competing-family models and a broad protein comparison database. [IgBLAST](https://ncbi.github.io/igblast/) can annotate receptor variable regions and germline matches, but a successful assignment is not proof of viral provenance; failure on a divergent or underrepresented species is not definitive exclusion.

## 5. Calibrate a pilot

Use up to 10,000 accession-pinned viral genomes, stratified by host group, genome type, sequence length, and annotation completeness; supplement strata if a simple random sample misses them. Use seed `20261112`. Add the positive and hard-negative panels. Benchmark known receptor sequences and their truncated/divergent versions as **software test inputs**, with unmodified biological records kept separate.

Test full domains and 30/60/90-amino-acid fragments, reverse strands, missed ORFs, polyprotein fusions, and selected frameshift/splice cases. Measure recall by receptor family, species, fragment length, and identity range; measure hard-negative misclassification and random-background hits. Simulated divergence is only a sensitivity probe, not a complete model of ancient evolution.

Pilot gate: recover all verified near-exact engineered IG/TR controls with usable sequence context; obtain at least 95% retrieval recall on held-out intact IG/TR domains; and correctly retain hard negatives as decoys or unresolved rather than natural acquisitions. Report denominators and uncertainty. If the gate fails, revise the panel/method and rerun the pilot, or explicitly narrow the claims; do not quietly proceed with an unsupported negative-result sensitivity claim.

Record wall time, peak RAM, disk amplification, and hit counts. Freeze the chosen thresholds and software/database versions in the run configuration before the broad screen. Pilot criteria concern retrieval/classification, not the existence of natural acquisitions.

## 6. Run complementary discovery routes

Take the **union** of candidates from independent routes. Requiring every method to agree at discovery would lose divergent or fragmented acquisitions.

1. **Protein sequence similarity:** search all viral proteins against the reference panel with a sensitive local aligner such as DIAMOND or MMseqs2; validate important alignments with BLASTP. Retain domain-level hits inside large proteins, not just full-length matches.
2. **Profile search:** run IG/TR family/domain HMMs over all viral proteins, independently of sequence-similarity hits. Apply broad Ig-family profiles as an auxiliary route and compare candidates with competing IgSF families. HMM and sequence scores are correlated evidence, not two independent experiments.
3. **Translated nucleotide search:** search viral nucleotide sequences independently of protein-screen outcomes, using a calibrated translated search and TBLASTN confirmation. Search both strands and retain coordinates/frame. Use annotated genetic codes and assess plausible alternative phage codes where necessary. This route addresses missing genes, partial coding regions, and annotation errors. Avoid one prokaryotic gene caller as the only annotation source for eukaryotic viruses.
4. **Nucleotide homology:** use germline and rearranged nucleotide references to find recent captures and gene remnants; inspect interrupted/spliced candidate loci with suitable alignments. Short D/J or motif-only matches cannot establish specific receptor ancestry.

Initial permissive candidate thresholds: sequence-search E-value at most `1e-3` or HMM independent-domain E-value at most `1e-3`, with at least 50 aligned amino acids for the main queue. Maintain a separate fragment/rescue queue for shorter matches and near-threshold findings (`1e-2`), calibrated against matched backgrounds. Record identity, bit score, aligned length, query/domain coverage, target coordinates, and low-complexity overlap; a universal percentage-identity cutoff is inappropriate.

These are retrieval settings, not significance tests. Account for total search space across shards and multiple profiles; use a common effective database size where supported or recompute comparable statistics. Do not interpret per-shard E-values as whole-corpus evidence. Keep sufficiently many reported alignments to avoid top-hit truncation; document limits and rerun saturated cases. Compare masked searches with controlled unmasked rescue and flag low-complexity-driven matches.

Full-corpus translated searches may be expensive. Batch and benchmark them; if only a subset is searched, report its exact scope and leave the nucleotide milestone incomplete. A protein-only search cannot support a claim about all viral genomic sequence.

## 7. Adjudicate identity, provenance, and genomic linkage

For each candidate, retrieve the original nucleotide record and its complete available contig/segment, not just the hit peptide. Build a compact evidence packet:

1. **Specific ancestry:** align to diverse IG/TR and competing IgSF families; inspect framework/CDR or constant-domain features and architecture. Conserved cysteines or an Ig fold alone are insufficient. Perform reciprocal broad-database searches, retaining alternative best explanations. Annotate receptor regions with IgBLAST and an appropriate numbering/profile method where applicable.
2. **Engineering provenance:** examine organism/source qualifiers, patents/publications, vector annotations, synthetic-construct labels, cloning adapters, tags/linkers, selection markers, and phage-display context. A missing engineering keyword is not proof of natural origin; a linker-like motif alone is not proof of engineering.
3. **Viral context:** inspect coordinates and flanks, neighboring hallmark genes, assembly graph if available, read coverage, and plausible viral/proviral boundaries. Use geNomad/CheckV where appropriate, with explicit caution outside their well-supported virus groups. Do not remove host-like loci automatically: the proposed acquisition is itself host-like. Distinguish a viral insert from a virus integrated next to a host IG/TR locus.
4. **Physical linkage:** where reads exist, map competitively to candidate, host/receptor loci, viral relatives, and relevant vector sequences. Seek independent molecules spanning both insert–virus junctions with unique anchors; remove duplicate/PCR support and inspect ambiguous mappings. Long reads or concordant independent assemblies strengthen the case. Uniform coverage, circularity, or a contig label alone does not prove integration. Packaged host DNA and expressed/chimeric RNA are separate outcomes.
5. **Replication of the observation:** search independent isolates/studies for the same locus and viral synteny. Several database copies of one assembly are one observation. If original reads are unavailable, retain the limitation; independent credible finished genomes may supply alternative support, but cannot create nonexistent junction-read evidence.
6. **Evolutionary interpretation:** compare viral relatives with and without the locus and construct domain-level alignments/trees with diverse host and viral homologs. Assess alignment quality and branch support. Short divergent domains may not resolve direction or timing of transfer; report that uncertainty. Optional structural comparison can support a fold assignment but cannot establish antibody/TCR ancestry or genomic incorporation by itself.

The roles of boundary validation are supported by [CheckV](https://www.nature.com/articles/s41587-020-00774-7) and [geNomad](https://pmc.ncbi.nlm.nih.gov/articles/PMC11324519/). Their predictions inform manual adjudication rather than replacing it.

Use the following outcome rules:

| Outcome | Required interpretation |
| --- | --- |
| Supported natural IG/TR acquisition | Specific IG/TR ancestry, credible natural-source provenance, convincing linkage within a viral genome, and independent linkage/assembly corroboration; no better supported engineering or contamination explanation. State whether only a fragment or chain was acquired. |
| Provisional candidate | Identity or linkage is suggestive but material evidence is missing or conflicting. Give the exact unresolved test. |
| Engineered insert | Documented or strongly supported laboratory construction; report complete viral genome, viral vector, and phagemid separately. |
| Other Ig-like protein | Viral locus is credible but receptor-specific ancestry is unsupported or another IgSF lineage fits better. |
| Artifact/non-genomic association | Host flank, chimeric assembly, misannotation, packaging without genomic integration, or protein-only association. |

## 8. Resources and staged execution

Local preflight found 64 logical CPUs, roughly 504 GiB total RAM, and about 169 GiB free disk on a shared filesystem. These are host observations, not a dedicated allocation. BLAST, HMMER, DIAMOND/MMseqs2, IgBLAST, NCBI Datasets, CheckV, and geNomad were not on the active PATH; this does not prove they are absent from every environment.

Start with an isolated, version-pinned environment and at most 16 CPU threads. Limit the initial pilot/cache to 20 GiB and keep at least 50 GiB filesystem headroom. Estimate downloaded, decompressed, indexed, and temporary sizes before each stage. Process one shard at a time where possible, retaining manifests and small evidence outputs. Measure performance before estimating full runtime; do not promise that a multi-million-genome translated search fits in one session or in available disk.

MetaVR's home page and API tutorial returned HTTP 403 through the browsing tool during planning. Direct API/bulk access has not been tested. Test the documented access route during acquisition; if unavailable, record the failure and continue the accessible corpora, using IMG/VR as a disclosed fallback. An inaccessible corpus remains a coverage gap. Access checks and capacity decisions are execution dependencies, not additional biological experiments.

Stages and deliverables:

1. **Plan/preflight (this milestone):** this protocol, issue, sources, inventory, and explicit pending status.
2. **Acquisition/calibration:** small frozen reference panel, real control accessions, software lock, accession manifests, pilot report, frozen search settings, measured resource estimate.
3. **Core search:** RefSeq, expanded GenBank, and INPHARED; per-corpus protein and nucleotide completion reports and preliminary candidates.
4. **Broad search:** MetaVR plus nonredundant supplemental catalogs, followed by the unclassified-sequence rescue; disclose unavailable sources or incomplete routes.
5. **Validation/report:** evidence packets, adjudication table, evolutionary interpretation where justified, and a bounded answer. Keep the issue open while required searches or candidate review are unfinished.

## 9. Reproducibility and stopping criteria

Implement a resumable workflow under this experiment directory after the plan milestone. Keep bulk inputs in ignored `data/raw/`, indexes in `data/cache/`, and temporary outputs in `work/`. Commit small manifests, scripts, parameter files, logs, and evidence tables. Future artifacts should include:

| Artifact | Minimum fields/content |
| --- | --- |
| Download manifest | Source/release, stable URL or query, retrieval UTC, accession.version, bytes, SHA-256, sequence count/total length, license/access notes, completion/error status |
| Run manifest | Git commit, software versions, exact command and parameters, input hashes, seed, start/end UTC, exit status, wall time, peak memory/disk, output hashes |
| Coverage table | Records and unique sequences supplied/searched/excluded/failed by corpus and route; bases/residues; host/genome-type/confidence strata; overlap accounting |
| Control results | Reference provenance, training/holdout assignment, expected outcome, retrieval/classification result, family/length/divergence sensitivity and uncertainty |
| Candidate table | Stable candidate ID, source accession/contig/segment, protein ID, 1-based inclusive coordinates and strand, proposed IG/TR family, alignment statistics, domain evidence, context, engineering flags, read/study links, decision and rationale |
| Evidence packet | Alignments, locus diagram, relevant metadata, competing explanations, junction-read summary, independent support, citations, reviewer/date, unresolved questions |

Completion requires accounting for every planned corpus and search route as completed, explicitly unavailable, or explicitly curtailed; reviewing all candidates above the frozen criteria; auditing a sample of near-threshold rejected hits; and reporting calibration and limitations. Do not silently label a curtailed or inaccessible search complete.

A positive answer must specify what was acquired and why natural viral residence is supported. A negative result must say **“no supported natural acquisition found in these releases with these methods and measured sensitivity”**, accompanied by exact search totals. It cannot establish “never.” Separate lack of detectable homology from lack of evolutionary occurrence, and preserve unresolved cases. Keep engineered precedents and generic Ig-like viral proteins in the final report even if the natural search is negative.
