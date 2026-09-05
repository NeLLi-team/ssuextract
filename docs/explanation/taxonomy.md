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

Runtime candidates that include such a sequence retain the ambiguous state.

## IMG-derived assignments

During database construction, exact matches to SILVA or PR2 retain the native
taxonomy. Other IMG sequences can receive similarity-derived assignments.
Cluster-only assignments stop at domain; lower ranks are not propagated to the member sequence. Detailed
output reports the matched sequence in `blast_sseqid`, its cluster centroid names
in `centroid_names`, and the calibrated taxonomy supported for those centroids in
`centroid_taxonomy`.

The classifier uses subjects with bit scores at least 98% of the best score and
requires at least 80% query coverage. It checks one hit beyond the 500-candidate
limit; tied overflow results are assigned at a higher shared rank or left
unclassified.

## Query assignment and ranked evidence

`cmsearch_summary.tsv` reports the query assignment. Reference taxonomy describes
the database sequence. A local match does not transfer that label to the query.

Runtime BLAST searches use the `blastn` task, an E-value limit of `1e-5`, and one
HSP per subject. The resolver requires at least 80% query coverage, calculated
from the query endpoints. It selects unique subjects with bit scores at least
98% of the best eligible score and finds their lowest common taxonomy. The
extra target detects candidate overflow. Unknown records cannot establish an
exact assignment. Explicit cross-domain conflicts remain ambiguous.

An exact native assignment requires every candidate to have the same native
SILVA or PR2 lineage, 100% identity over the complete query and subject, and no
mismatches or gaps. A PR2 species call also requires all nine taxonomy ranks.
The reference must have a native assignment method; a derived IMG label cannot
use this exception. A legacy 12-column BLAST file lacks subject length and
cannot prove the exact-match condition.

Nonexact query calls require schema-3 runtime calibration for the selected
profile and search policy. A missing calibration leaves the query unclassified.
The candidate LCA, best reference taxonomy, and ranked evidence remain visible.
Database v1.0.2 profiles do not contain this runtime calibration. Their schema-2
IMG centroid calibration does not validate runtime query assignment.

Schema-2 centroid calibration groups calls by the true query class. It does not
establish precision for the predicted class. Stored v1.0.2 centroid labels use
that older artifact; corrected labels require a rebuilt database profile.

Runtime calibration groups calls by predicted marker, source set, and domain.
False calls from other true classes count against the predicted group. The
calibration file records the reference content hash and policy. Accepted rank
limits or identity rules constrain the query call. Precision depends on the test class
mixture and reference coverage; the calibration report must also state call
rates and held-out sample counts.

Runtime calibration fixes the BLAST target limit. A larger limit can admit calls
that the calibrated policy rejected as truncated. Taxonomy assignment stops if
`--max_blast_targets` differs from the calibrated `max_targets`.

Nonexact PR2 calls cannot reach species. Even an exact SSU match describes
agreement with the available reference records; it does not establish
genome-wide species identity. For example, bacterial 16S similarity does not
define a general species boundary based on whole genomes.
[Primary study](https://pmc.ncbi.nlm.nih.gov/articles/PMC11264914/).

`blast_top_hits.tsv` retains assignment candidates as individual evidence rows.
A shorter match outside the coverage rule can retain a deep native reference
lineage in that table, but it cannot supply the selected query taxonomy.

The top-hit table retains the best subject from each available reference source
among the fetched candidates, including subjects outside the requested overall
cutoff. Those rows expose native PR2 and SILVA taxonomy when an IMG sequence has
no supported assignment.

## Tree-neighbor assignment

Tree classification is disabled by default. When enabled, SSUextract searches
the extracted gene against both marker indexes and compares the best 100 unique
subjects that pass the coverage rule. The majority route selects the 16S rRNA
gene or 18S rRNA gene model; best bit score and the accepted Infernal model are
recorded tie-breaks. The 16S
rRNA gene route includes bacterial, archaeal, and organellar references.

The selected references and query are aligned with the chosen covariance model.
Covariance-model insert columns and match columns with more than 90% gaps are
removed before IQ-TREE 3 estimates branch lengths. References are then ordered
by patristic distance from the query. Direct reference labels determine the
neighbor LCA. Centroid taxonomy remains evidence and cannot add ranks.

All references tied at the neighbor distance boundary contribute. An ambiguous
reference within that boundary retains its alternatives. The final tree call
is the shared prefix of the neighbor LCA and the supported BLAST query call on
the selected marker. Conflicting domains produce ambiguity. Missing named tree
evidence leaves the supported BLAST call in place. Branch distance and SH-aLRT
support are reported as evidence; they do not supply an independent taxonomic
confidence threshold.

`taxonomy_mode` identifies the selected method. In tree mode, `taxonomy` holds
the tree-neighbor result while the `blast_taxonomy` fields preserve the normal
BLAST result. The neighbor table records every reference distance, lineage
source, and whether that reference contributed to the assignment.

At least three reference subjects are required to infer a tree. A query below
that threshold keeps its BLAST assignment, reports `taxonomy_mode=blast`, and
records `tree_skipped_insufficient_references` in `tree_assignment_method`.
