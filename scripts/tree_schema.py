from __future__ import annotations


REFERENCE_FIELDS = [
    "leaf_id",
    "blast_sseqid",
    "hit_rank",
    "reference_identifiers",
    "reference_versions",
    "reference_source",
    "taxonomy",
    "taxonomy_source",
    "taxonomy_domain",
    "compartment",
    "taxonomy_assignment_method",
    "centroid_names",
    "centroid_taxonomy",
    "centroid_taxonomy_source",
    "blast_pident",
    "blast_length",
    "blast_evalue",
    "blast_bitscore",
]

TREE_ASSIGNMENT_FIELDS = [
    "name",
    "sample",
    "model",
    "tree_model",
    "tree_marker",
    "tree_route_decision",
    "tree_route_16s_votes",
    "tree_route_18s_votes",
    "tree_route_16s_best_bitscore",
    "tree_route_18s_best_bitscore",
    "tree_taxonomy",
    "tree_taxonomy_source",
    "tree_taxonomy_domain",
    "tree_compartment",
    "tree_assignment_method",
    "tree_basis_neighbors",
    "tree_nearest_sseqid",
    "tree_nearest_reference_identifiers",
    "tree_nearest_distance",
    "tree_query_edge_support",
    "tree_inference_model",
]

TREE_NEIGHBOR_FIELDS = [
    "name",
    "sample",
    "model",
    "tree_model",
    "tree_marker",
    "tree_neighbor_rank",
    "tree_distance",
    "used_for_assignment",
    "tree_lineage",
    "tree_lineage_source",
    "tree_lineage_basis",
    *REFERENCE_FIELDS,
]

SUMMARY_TREE_FIELDS = [
    "taxonomy_mode",
    "blast_taxonomy",
    "blast_taxonomy_source",
    "blast_taxonomy_domain",
    "blast_compartment",
    "blast_taxonomy_assignment_method",
    "blast_taxonomy_alternatives",
    *TREE_ASSIGNMENT_FIELDS[3:],
]
