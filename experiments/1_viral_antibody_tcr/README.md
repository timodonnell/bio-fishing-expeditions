# Have viruses acquired antibody or TCR genes?

Issue: [#1](https://github.com/timodonnell/bio-fishing-expeditions/issues/1)

Status, 2026-11-13: **MetaVR expansion complete: 369,319,989 protein records screened. Natural TCR beta-chain capture is established from published FeLV evidence; natural antibody capture remains unresolved.**

The answer is yes for a naturally captured **TCR beta chain**: GenBank **X05155.1**, protein **CAA28801.1**, UniProt **P11364**, reported in 1987. The source CDS translation was verified and the sequence recovered with host-only references. This is a previously known case, not a novel discovery or a complete alpha/beta receptor. Natural antibody capture remains unresolved.

The new MetaVR profile search covers **369,319,989 protein records** associated with **24,429,973 source-genome identifiers**, about **115×** the prior protein-record total. Of 26 retained candidates, two match host immunoglobulin loci, 14 have competing MHC homology, one matches CD300, and **nine environmental IgV-like records remain unresolved**. One has neighboring genes with strong viral homology; receptor ancestry still requires validation. No new supported natural capture was established.

Earlier searches covered **3,207,211 protein records** by HMM/BLASTP and **48,286 genome records** by TBLASTN across RefSeq viruses and INPHARED phages. These totals include overlap; the MetaVR expansion used HMMs and did not repeat whole-corpus independent homology or translated-DNA routes. A targeted annotation survey added 960 records and five unresolved annotation-only genomic leads.

- [MetaVR scale-up results and limitations](SCALE_UP_RESULTS.md)
- [Earlier RefSeq/INPHARED results](RESULTS.md)
- [Natural FeLV TCR evidence packet](known_natural_tcr.md)
- [Original research plan](PLAN.md)
- [MetaVR scale-up plan](SCALE_UP.md) — acquisition, full protein-profile search, and candidate review complete
- [Reproduction instructions](RUNNING.md)
- [Data-source inventory](data_sources.tsv)
- [Download manifest](manifests/downloads.json), [run outputs](results/), and [integrity validation](results/validation.json)
- [Frozen references](references/), [scripts](scripts/), and [software lock](requirements.lock)
- [Historical literature preflight](literature_preflight.md), [sources](sources.tsv), and [queries](search_queries.tsv)

The full calibration gate, whole-GenBank search, independent whole-MetaVR homology/translated-DNA routes, supplemental catalogs, and read/phylogenetic validation are unfinished. The issue remains open. No absence result here establishes that a natural antibody acquisition has never occurred.
