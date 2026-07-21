# Taxonomy sources and selection

SSUextract stores the native SILVA and PR2 taxonomy strings. It does not map
their ranks onto a shared hierarchy.

## Preferred source

- SILVA 138.2 supplies taxonomy for Bacteria and Archaea.
- PR2 5.1.1 supplies eukaryotic nuclear and nucleomorph taxonomy.
- PR2 supplies the host taxonomy for plastid, apicoplast, and mitochondrial
  records. The organelle compartment is stored separately.

## Exact cross-domain identity

An exact nucleotide sequence can occur in SILVA as a bacterial record and in
PR2 as a plastid or other organellar record. The curated build contains 1,650
such exact-sequence conflicts. SSUextract records:

- `domain = ambiguous`
- `compartment = mixed`
- an empty preferred taxonomy
- all native alternatives in structured evidence

Equal-best runtime hits to such a sequence retain the ambiguous state.

## IMG-derived assignments

Exact matches to SILVA or PR2 retain the complete native taxonomy. Other IMG
sequences can receive similarity-derived assignments. Cluster-only assignments
stop at domain; lower ranks are not propagated to the member sequence. Detailed
output reports the matched sequence in `blast_sseqid`, its cluster centroid names
in `centroid_names`, and the calibrated taxonomy supported for those centroids in
`centroid_taxonomy`.

The classifier retains candidate ties within 98% of the best bit score and
requires at least 80% query coverage. It checks one hit beyond the 500-candidate
limit; tied overflow results are assigned at a higher shared rank or left
unclassified.

## Assignment and ranked evidence

`cmsearch_summary.tsv` reports the selected assignment. `blast_top_hits.tsv`
reports ranked reference evidence and does not vote across lower-scoring hits.
Equal-best subjects can reduce the summary assignment to their lowest common
taxonomy; the top-hit table keeps every equal-best subject as an individual row.
An exact, full-length IMG match can remain unclassified when neither the member
nor its centroid has calibrated taxonomy. A shorter PR2 or SILVA match may have
a deeper lineage, but that lineage remains a separate evidence row.

The top-hit table retains the best subject from each available reference source
among the fetched candidates, including subjects outside the requested overall
cutoff. This makes native PR2 and SILVA taxonomy visible without transferring it
to an IMG sequence that lacks supporting classification.

## Tree-neighbor assignment

Tree classification is disabled by default. When enabled, SSUextract searches
the extracted gene against both marker indexes and compares the best 100 unique
subjects. The majority route selects the 16S rRNA gene or 18S rRNA gene model;
best bit score and the accepted Infernal model are recorded tie-breaks. The 16S
rRNA gene route includes bacterial, archaeal, and organellar references.

The selected references and query are aligned with the chosen covariance model.
Covariance-model insert columns and match columns with more than 90% gaps are
removed before IQ-TREE 3 estimates branch lengths. References are then ordered
by patristic distance from the query. The selected taxonomy is the lowest common
taxonomy of the nearest named references. For an IMG reference, calibrated
centroid taxonomy is used when it is deeper than the member assignment.

All named references tied by patristic distance at the neighbor cutoff
contribute to the LCA. BLAST rank therefore cannot resolve an unsupported tree
tie in favor of a narrower lineage.

`taxonomy_mode` identifies the selected method. In tree mode, `taxonomy` holds
the tree-neighbor result while the `blast_taxonomy` fields preserve the normal
BLAST result. The neighbor table records every reference distance, lineage
source, and whether that reference contributed to the assignment.

At least three reference subjects are required to infer a tree. A query below
that threshold keeps its BLAST assignment, reports `taxonomy_mode=blast`, and
records `tree_skipped_insufficient_references` in `tree_assignment_method`.
