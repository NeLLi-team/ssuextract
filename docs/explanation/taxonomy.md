# Taxonomy policy

SILVA and PR2 serve different domains and use different rank systems.
SSUextract retains those native systems instead of flattening both into a shared
set of rank names.

## Native authority

- SILVA 138.2 supplies taxonomy for Bacteria and Archaea.
- PR2 5.1.1 supplies eukaryotic nuclear and nucleomorph taxonomy.
- PR2 host taxonomy describes plastid, apicoplast, and mitochondrial records;
  the compartment is stored separately.

This policy preserves PR2's fixed eukaryotic hierarchy and SILVA's prokaryotic
lineage without inventing unsupported rank equivalences.

## Exact cross-domain identity

An exact nucleotide sequence can occur in SILVA as a bacterial record and in
PR2 as a plastid or other organellar record. The curated build contains 1,650
such exact-sequence conflicts. SSUextract records:

- `domain = ambiguous`
- `compartment = mixed`
- an empty preferred taxonomy
- all native alternatives in structured evidence

The runtime annotation layer preserves this state when equal-best hits include
the conflicted sequence.

## IMG-derived assignments

Exact matches to a current SILVA or PR2 sequence retain the complete native
taxonomy. Other IMG sequences can receive similarity-derived evidence. Cluster
propagation is capped at domain because the source cluster construction and
within-cluster taxonomic coherence are not documented well enough to support
lower-rank propagation.

The classifier evaluates candidate ties within 98% of the best bit score and
requires at least 80% query coverage. It fetches one hit beyond the 500-candidate
policy limit; a tied overflow backs the result off to a defensible rank or leaves
it unclassified.

