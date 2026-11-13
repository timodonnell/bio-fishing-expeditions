# Natural TCR capture by feline leukemia virus

Assessment: **supported natural acquisition of a TCR beta-chain coding sequence**. This answers the TCR part of issue [#1](https://github.com/timodonnell/bio-fishing-expeditions/issues/1). It is a previously published example, recovered during the reference/annotation survey on 2026-11-13, not a novel discovery from the bulk viral screen.

## Evidence and scope

[Fulton et al., Nature (1987)](https://pubmed.ncbi.nlm.nih.gov/3029597/), DOI [10.1038/326190a0](https://doi.org/10.1038/326190a0), reported retroviral transduction of a T-cell antigen receptor beta-chain gene in a naturally occurring feline tumor. The associated deposited sequence is **X05155.1**, clone **T17T-22**, encoding **CAA28801.1**, also represented by reviewed UniProt **P11364**.

[Terry et al., Journal of Virology (1992)](https://pmc.ncbi.nlm.nih.gov/articles/PMC241135/), DOI [10.1128/jvi.66.6.3538-3549.1992](https://doi.org/10.1128/jvi.66.6.3538-3549.1992), further characterized the T17 material. It contained separate receptor-bearing, myc-bearing, and helper FeLV proviruses. The receptor and myc acquisitions had distinct junctions, and the genomes were packaged. Recombinant forms were lost preferentially during passage. The study did not establish a readily transmissible oncogenic function for the captured receptor sequence.

The acquired material is **one TCR beta chain**, not a complete alpha/beta receptor and not an antibody. Its presence does not establish an antigen specificity, a complete functional receptor, or an adaptive advantage. The literature describes a natural tumor-associated acquisition; later experiments with cloned material do not make that original acquisition an engineered event.

## Deposited locus

[GenBank X05155.1](https://www.ncbi.nlm.nih.gov/nuccore/X05155.1) is a **1,862-base partial viral record**, not a complete viral genome. Its source specifies feline leukemia virus and thymic lymphosarcoma T17. Coordinates below are 1-based, inclusive, on the forward strand; fuzzy endpoints retain the record's uncertainty.

| Feature | Coordinates |
| --- | --- |
| FeLV pol fragment | <1–86 |
| TCR beta CDS, including stop | 94–1059 |
| Leader | 94–174 |
| V segment | 175–459 |
| D segment | 460–477 |
| J segment | 478–525 |
| Constant region | 526–1056 |
| Viral long terminal repeat fragment | 1823–>1862 |

The deposited record notes a historical sequence correction at bases 567–568. The retrieved sequence, rather than a transcription from the original printed article, was used here.

The extracted CDS translates exactly to the **321-amino-acid** GenBank and UniProt proteins. The validation script checks the stop codon and equality against both records. Machine-readable features and hashes are in [known_positive_validation.json](results/known_positive_validation.json). Biopython feature strings in that JSON use 0-based, half-open coordinates; the explicit CDS fields and this table use 1-based coordinates.

## Independent computational recovery

P11364 and its feline sequence were withheld from the human/mouse training references. The host-derived profiles recover both variable and constant receptor regions. BLASTP retrieves the protein, and TBLASTN retrieves the source nucleotide locus on both the deposited orientation and an exact reverse-complement software control. The reverse-complement record is only a test input, not a second biological observation.

The source accession and the exact P11364 amino-acid sequence are absent from the downloaded RefSeq corpus, as checked in [refseq_expansion_summary.json](results/refseq_expansion_summary.json). A negative RefSeq-only search therefore cannot exclude even this known acquisition.

## Provenance and limits

Inspected: the two PubMed abstracts/metadata, the accessible PMC abstract, the deposited GenBank record, and the UniProt JSON record. Full-text PDF/XML retrieval attempts failed or returned HTML and are marked as failed in [downloads.json](manifests/downloads.json). The full historical experimental methods and raw data have not been independently reanalyzed. No new raw-read junction validation was performed.

The natural-origin conclusion rests on the primary reports together with their deposited viral locus. A database label alone would not establish it. This evidence resolves an existence question for TCR capture; it does not measure its frequency or resolve whether a natural antibody capture exists elsewhere.
