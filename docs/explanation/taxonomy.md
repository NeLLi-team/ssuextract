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
