# Literature preflight — 2026-11-12

This is a targeted discovery pass to design the experiment. It is not the completed systematic review or a sequence-database search. No natural-acquisition candidate has been validated.

## What the initial sources establish

| Observation | Source and inspected evidence | Consequence for the experiment |
| --- | --- | --- |
| Antibody variable domains have been displayed on engineered fd bacteriophage. | McCafferty et al., 1990, [Nature](https://www.nature.com/articles/348552a0). Publisher-indexed abstract inspected; direct full-text opening failed. | An engineered precedent and lead for positive controls. Retrieve the construct sequence and full methods before treating a particular accession as a genomic control. |
| TCR alpha chains have been displayed using an engineered M13 system. | Onda et al., 1995, [PubMed](https://pubmed.ncbi.nlm.nih.gov/8643108/). Abstract inspected. | TCR engineering precedent. The abstract alone does not settle complete phage genome versus phagemid carriage, nor establish a complete alpha/beta receptor. |
| Antibody genes have been delivered by recombinant AAV. | Lewis et al., 2002, [primary paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC136414/), [publisher](https://journals.asm.org/doi/full/10.1128/jvi.76.17.8769-8775.2002). Indexed methods describe an antibody expression cassette between AAV terminal repeats. | Engineered viral-vector genome evidence; distinguish vectors from natural viral genomes and from their production plasmids. |
| Naturally occurring phage structural proteins can have Ig-like domains. | Fraser et al., 2006, [primary paper](https://doi.org/10.1016/j.jmb.2006.03.043), abstract; Sathaliyawala et al., 2010, [Hoc study](https://pmc.ncbi.nlm.nih.gov/articles/PMC2909354/), indexed discussion. | Essential hard negatives. An Ig-like fold or profile match does not identify an antibody/TCR-derived gene. |
| Human immunoglobulins can be carried on HCMV particles. | [Primary experimental study](https://www.frontiersin.org/journals/microbiology/articles/10.3389/fmicb.2022.1106401/full), indexed results/discussion. | Protein incorporation and genomic incorporation must be separate outcomes; this observation does not demonstrate acquisition of immunoglobulin coding sequence. |

The targeted natural-capture queries did not yield a source that this preflight could establish as a natural antibody/TCR gene acquisition. This is an observation about this limited search, not a claim that no such publication or genome exists. `search_queries.tsv` records selected executed queries and a separate planned systematic query set.

## Database and access findings

The [2026 MetaVR paper](https://pubmed.ncbi.nlm.nih.gov/41312645/) identifies MetaVR as the successor of IMG/VR and reports 24,435,662 uncultivated viral genomes with programmatic access. Use it as the preferred broad catalog, with older IMG/VR coverage documented separately. This published count is not an experiment search total.

The [INPHARED README](https://github.com/RyanCook94/inphared/blob/main/README.md), retrieved through GitHub, linked a 14Apr2025 snapshot. This is an observed downloadable-release lead, not proof that no newer release exists. Its phage data substantially overlap GenBank.

Direct browsing of MetaVR's home page and API tutorial returned HTTP 403. Several publisher/PMC direct opens also failed or returned a browser challenge; indexed abstracts and excerpts were available as indicated above. Direct sequence downloads and MetaVR API requests have not been attempted. These access limits must be resolved or recorded during acquisition, and no download/checksum/coverage claims are made here.

## Current answer

There are documented engineered antibody/TCR display and antibody viral-vector precedents, as well as natural viral Ig-like proteins and protein-only antibody association. **Whether a natural viral genome has acquired sequence specifically derived from an antibody or TCR remains unresolved by this experiment.** The planned sequence search and provenance checks address that narrower question.
