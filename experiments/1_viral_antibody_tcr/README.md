# Have viruses acquired antibody or TCR genes?

Issue: [#1](https://github.com/timodonnell/bio-fishing-expeditions/issues/1)

Status, 2026-11-13: **natural TCR beta-chain capture established from published FeLV evidence; core RefSeq/INPHARED screens complete; broader search and unresolved candidates remain open.**

The answer is yes for a naturally captured **TCR beta chain**: GenBank **X05155.1**, protein **CAA28801.1**, UniProt **P11364**, reported in 1987. The source CDS translation was verified and the sequence recovered with host-only references. This is a previously known case, not a novel discovery or a complete alpha/beta receptor. Natural antibody capture remains unresolved.

We searched **3,207,211 protein records** and translated **48,286 genome records** across RefSeq viruses and INPHARED phages. These totals include overlap; exact-sequence overlap is reported separately. A targeted annotation survey added 960 viral-taxonomy protein records and exposed structure-component and legacy-taxonomy pitfalls. Five annotation-only genomic leads remain unresolved.

- [Results and limitations](RESULTS.md)
- [Natural FeLV TCR evidence packet](known_natural_tcr.md)
- [Original research plan](PLAN.md)
- [Reproduction instructions](RUNNING.md)
- [Data-source inventory](data_sources.tsv)
- [Download manifest](manifests/downloads.json), [run outputs](results/), and [integrity validation](results/validation.json)
- [Frozen references](references/), [scripts](scripts/), and [software lock](requirements.lock)
- [Historical literature preflight](literature_preflight.md), [sources](sources.tsv), and [queries](search_queries.tsv)

The full calibration gate, whole-GenBank search, MetaVR, supplemental metagenomic catalogs, and read/phylogenetic validation are unfinished. The issue remains open. No absence result here establishes that a natural antibody acquisition has never occurred.
