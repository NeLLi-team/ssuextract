#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import sys
import tempfile
from collections.abc import Iterable, Sequence
from pathlib import Path

import duckdb


SCRIPTS = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(SCRIPTS))

from database_contracts import PR2_TAXONOMY_SUFFIX_TO_COMPARTMENT


TABLE_NAMES = (
    "sequences",
    "source_records",
    "taxonomy_assignments",
    "preferred_taxonomy",
    "img_location",
)
SOURCE_DISPLAY = (
    {
        "source": "SILVA",
        "reference_source": "SILVA",
        "marker": None,
        "content": "16S rRNA gene references; bacterial and archaeal taxonomy",
        "config_key": "silva_ssu_nr99_fasta",
    },
    {
        "source": "PR2",
        "reference_source": "PR2",
        "marker": None,
        "content": "18S rRNA gene and organellar 16S rRNA gene references; eukaryotic taxonomy",
        "config_key": "pr2_ssu_fasta",
    },
    {
        "source": "IMG 16S rRNA genes",
        "reference_source": "IMG",
        "marker": "16S",
        "content": "Additional IMG 16S rRNA gene sequences",
        "config_key": "eukcensus_16s_fasta",
    },
    {
        "source": "IMG 18S rRNA genes",
        "reference_source": "IMG",
        "marker": "18S",
        "content": "Additional IMG 18S rRNA gene sequences",
        "config_key": "eukcensus_18s_fasta",
    },
)
DOMAIN_RANKS = {
    "Bacteria": ("phylum", "SILVA"),
    "Archaea": ("phylum", "SILVA"),
    "Eukaryota": ("supergroup", "PR2"),
}
PR2_COMPARTMENT_SUFFIXES = tuple(PR2_TAXONOMY_SUFFIX_TO_COMPARTMENT)


class CompositionError(RuntimeError):
    pass


def _reject_duplicate_keys(pairs: list[tuple[str, object]]) -> dict:
    result: dict[str, object] = {}
    for key, value in pairs:
        if key in result:
            raise CompositionError(f"duplicate JSON key: {key!r}")
        result[key] = value
    return result


def read_json(path: Path) -> dict:
    try:
        value = json.loads(path.read_text(encoding="utf-8"), object_pairs_hook=_reject_duplicate_keys)
    except (OSError, ValueError) as error:
        raise CompositionError(f"cannot read JSON {path}: {error}") from error
    if not isinstance(value, dict):
        raise CompositionError(f"JSON root must be an object: {path}")
    return value


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def atomic_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(prefix=f".{path.name}.", dir=path.parent)
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8", newline="") as handle:
            handle.write(text)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
        path.chmod(0o644)
    except BaseException:
        try:
            os.unlink(temporary)
        except FileNotFoundError:
            pass
        raise


def atomic_tsv(path: Path, fieldnames: Sequence[str], records: Iterable[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(prefix=f".{path.name}.", dir=path.parent)
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(records)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
        path.chmod(0o644)
    except BaseException:
        try:
            os.unlink(temporary)
        except FileNotFoundError:
            pass
        raise


def table_paths(profile_root: Path) -> dict[str, Path]:
    return {name: profile_root / "tables" / f"{name}.parquet" for name in TABLE_NAMES}


def artifact_map(manifest: dict) -> dict[str, dict]:
    artifacts = manifest.get("artifacts")
    if not isinstance(artifacts, list):
        raise CompositionError("profile manifest artifacts must be a list")
    mapped = {item.get("path"): item for item in artifacts if isinstance(item, dict)}
    if len(mapped) != len(artifacts):
        raise CompositionError("profile manifest artifact paths must be unique")
    return mapped


def validate_profile(
    profile: str,
    root: Path,
    qa_path: Path,
    catalog: dict,
) -> tuple[dict, dict, dict, dict[str, Path]]:
    manifest_path = root / "manifest.json"
    provenance_path = root / "provenance.json"
    manifest = read_json(manifest_path)
    provenance = read_json(provenance_path)
    qa = read_json(qa_path)
    paths = table_paths(root)

    if manifest.get("profile") != profile or provenance.get("profile") != profile:
        raise CompositionError(f"profile identity mismatch for {profile}")
    version = catalog["profiles"][profile]["version"]
    if manifest.get("version") != version or provenance.get("version") != version:
        raise CompositionError(f"profile version mismatch for {profile}")
    if qa.get("profile_manifest_sha256") != sha256(manifest_path):
        raise CompositionError(f"QA manifest digest mismatch for {profile}")
    archive = catalog["profiles"][profile]["archive"]
    if qa.get("archive", {}).get("sha256") != archive["sha256"]:
        raise CompositionError(f"catalog archive digest mismatch for {profile}")
    if qa.get("archive", {}).get("bytes") != archive["bytes"]:
        raise CompositionError(f"catalog archive size mismatch for {profile}")

    artifacts = artifact_map(manifest)
    for name, path in paths.items():
        if not path.is_file():
            raise CompositionError(f"missing {profile} table: {name}")
        relative = f"tables/{name}.parquet"
        artifact = artifacts.get(relative)
        if artifact is None:
            raise CompositionError(f"manifest does not list {profile} {relative}")
        if artifact.get("bytes") != path.stat().st_size or artifact.get("sha256") != sha256(path):
            raise CompositionError(f"manifest does not bind {profile} {relative}")
    return manifest, provenance, qa, paths


def fetch_dicts(connection: duckdb.DuckDBPyConnection, query: str, parameters=()) -> list[dict]:
    cursor = connection.execute(query, parameters)
    names = [item[0] for item in cursor.description]
    return [dict(zip(names, row, strict=True)) for row in cursor.fetchall()]


def scalar(connection: duckdb.DuckDBPyConnection, query: str, parameters=()):
    return connection.execute(query, parameters).fetchone()[0]


def composition_rows(
    connection: duckdb.DuckDBPyConnection,
    profile: str,
    version: str,
    paths: dict[str, Path],
    qa: dict,
) -> list[dict]:
    sequence_path = str(paths["sequences"])
    source_path = str(paths["source_records"])
    preferred_path = str(paths["preferred_taxonomy"])
    rows: list[dict] = [
        {
            "profile": profile,
            "version": version,
            "dimension": "total_unique_sequences",
            "category": "All",
            "marker": "",
            "count": qa["counts"]["sequences"],
        }
    ]
    for record in fetch_dicts(
        connection,
        """SELECT markers AS category, count(*) AS count
           FROM read_parquet(?) GROUP BY markers ORDER BY markers""",
        [sequence_path],
    ):
        rows.append(
            {
                "profile": profile,
                "version": version,
                "dimension": "marker_membership",
                "category": record["category"],
                "marker": "",
                "count": record["count"],
            }
        )
    for record in fetch_dicts(
        connection,
        """SELECT reference_source AS category, marker, count(*) AS count
           FROM read_parquet(?) GROUP BY ALL ORDER BY reference_source, marker""",
        [source_path],
    ):
        rows.append(
            {
                "profile": profile,
                "version": version,
                "dimension": "source_records",
                **record,
            }
        )
    for record in fetch_dicts(
        connection,
        """SELECT CASE WHEN cross_domain_conflict THEN 'Ambiguous' ELSE domain END AS category,
                  count(*) AS count
           FROM read_parquet(?) GROUP BY category ORDER BY category""",
        [preferred_path],
    ):
        rows.append(
            {
                "profile": profile,
                "version": version,
                "dimension": "taxonomy_domain",
                "category": record["category"],
                "marker": "",
                "count": record["count"],
            }
        )
    if profile == "img":
        location = qa["img_location"]
        location_path = str(paths["img_location"])
        observed = fetch_dicts(
            connection,
            """SELECT count(*) AS rows,
                      count(latitude) AS latitude_rows,
                      count(longitude) AS longitude_rows,
                      count(*) FILTER (WHERE latitude IS NOT NULL AND longitude IS NOT NULL)
                        AS complete_rows
               FROM read_parquet(?)""",
            [location_path],
        )[0]
        for key in ("rows", "latitude_rows", "longitude_rows"):
            if observed[key] != location[key]:
                raise CompositionError(f"IMG location QA count mismatch for {key}")
        complete = observed["complete_rows"]
        for category, count in (
            ("Located", complete),
            ("Coordinates unavailable", location["rows"] - complete),
        ):
            rows.append(
                {
                    "profile": profile,
                    "version": version,
                    "dimension": "img_taxon_coordinates",
                    "category": category,
                    "marker": "",
                    "count": count,
                }
            )
    return rows


def taxonomy_rank_rows(
    connection: duckdb.DuckDBPyConnection,
    profile: str,
    version: str,
    preferred_path: Path,
) -> tuple[list[dict], list[dict]]:
    observed_suffixes = {
        record["suffix"]
        for record in fetch_dicts(
            connection,
            """WITH tokens AS (
                   SELECT trim(list_extract(string_split(taxonomy, ';'), 2)) AS token
                   FROM read_parquet(?)
                   WHERE domain = 'Eukaryota'
               )
               SELECT DISTINCT regexp_extract(token, ':([^:]+)$', 1) AS suffix
               FROM tokens WHERE token LIKE '%:%'""",
            [str(preferred_path)],
        )
    }
    if observed_suffixes:
        raise CompositionError(
            f"PR2 compartment suffixes remain in preferred taxonomy for {profile}: "
            f"{sorted(observed_suffixes)}"
        )
    connection.execute(
        """CREATE OR REPLACE TEMP TABLE rank_assignments AS
           WITH parsed AS (
               SELECT domain,
                      trim(list_extract(string_split(taxonomy, ';'), 2)) AS raw_lineage
               FROM read_parquet(?)
               WHERE domain IN ('Bacteria', 'Archaea', 'Eukaryota')
           )
           SELECT domain, raw_lineage AS lineage
           FROM parsed""",
        [str(preferred_path)],
    )
    connection.execute(
        """CREATE OR REPLACE TEMP TABLE resolved_rank_assignments AS
           SELECT domain, lineage
           FROM rank_assignments
           WHERE lineage IS NOT NULL
             AND lineage <> ''
             AND lower(lineage) NOT IN ('unclassified', 'incertae sedis')"""
    )

    domain_totals = {
        record["domain"]: record["total_sequences"]
        for record in fetch_dicts(
            connection,
            """SELECT domain, count(*) AS total_sequences
               FROM rank_assignments GROUP BY domain""",
        )
    }
    resolved = fetch_dicts(
        connection,
        """SELECT domain, count(*) AS resolved_sequences,
                  count(DISTINCT lineage) AS named_lineages
           FROM resolved_rank_assignments
           GROUP BY domain""",
    )
    resolved_by_domain = {record["domain"]: record for record in resolved}
    resolution_rows = []
    for domain, (target_rank, _taxonomy_source) in DOMAIN_RANKS.items():
        total = int(domain_totals.get(domain, 0))
        details = resolved_by_domain.get(
            domain, {"resolved_sequences": 0, "named_lineages": 0}
        )
        resolved_count = int(details["resolved_sequences"])
        resolution_rows.append(
            {
                "profile": profile,
                "version": version,
                "domain": domain,
                "target_rank": target_rank,
                "total_sequences": total,
                "resolved_sequences": resolved_count,
                "unresolved_sequences": total - resolved_count,
                "resolution_percent": f"{resolved_count / total * 100:.3f}" if total else "0.000",
                "named_lineages": int(details["named_lineages"]),
            }
        )

    lineage_rows = []
    for record in fetch_dicts(
        connection,
        """SELECT domain, lineage, count(*) AS count
           FROM resolved_rank_assignments
           GROUP BY domain, lineage
           ORDER BY CASE domain
                      WHEN 'Bacteria' THEN 1
                      WHEN 'Archaea' THEN 2
                      WHEN 'Eukaryota' THEN 3
                    END,
                    count DESC, lineage""",
    ):
        lineage_rows.append(
            {
                "profile": profile,
                "version": version,
                "domain": record["domain"],
                "target_rank": DOMAIN_RANKS[record["domain"]][0],
                "lineage": record["lineage"],
                "count": record["count"],
            }
        )
    if sum(row["count"] for row in lineage_rows) != sum(
        row["resolved_sequences"] for row in resolution_rows
    ):
        raise CompositionError(f"taxonomy rank totals do not reconcile for {profile}")
    return resolution_rows, lineage_rows


def write_location_tables(
    connection: duckdb.DuckDBPyConnection,
    source_path: Path,
    location_path: Path,
    output_directory: Path,
) -> tuple[Path, Path, dict]:
    source = str(source_path)
    locations = str(location_path)
    columns = [
        row[0]
        for row in connection.execute("DESCRIBE SELECT * FROM read_parquet(?)", [locations]).fetchall()
    ]
    if columns != ["taxon_oid", "latitude", "longitude"]:
        raise CompositionError(f"IMG location privacy allowlist mismatch: {columns}")

    connection.execute(
        """CREATE OR REPLACE TEMP TABLE img_taxa AS
           SELECT s.taxon_oid, l.latitude, l.longitude,
                  count(*) FILTER (WHERE s.marker = '16S') AS records_16s,
                  count(*) FILTER (WHERE s.marker = '18S') AS records_18s,
                  count(DISTINCT s.sequence_id) FILTER (WHERE s.marker = '16S') AS sequences_16s,
                  count(DISTINCT s.sequence_id) FILTER (WHERE s.marker = '18S') AS sequences_18s,
                  CASE
                    WHEN count(*) FILTER (WHERE s.marker = '16S') > 0
                     AND count(*) FILTER (WHERE s.marker = '18S') > 0 THEN '16S + 18S'
                    WHEN count(*) FILTER (WHERE s.marker = '16S') > 0 THEN '16S only'
                    WHEN count(*) FILTER (WHERE s.marker = '18S') > 0 THEN '18S only'
                    ELSE 'No IMG SSU record'
                  END AS marker_group
           FROM (
               SELECT * FROM read_parquet(?) WHERE reference_source = 'IMG'
           ) s
           LEFT JOIN read_parquet(?) l ON s.taxon_oid = l.taxon_oid
           GROUP BY s.taxon_oid, l.latitude, l.longitude
           ORDER BY s.taxon_oid""",
        [source, locations],
    )
    if scalar(connection, "SELECT count(*) - count(DISTINCT taxon_oid) FROM img_taxa") != 0:
        raise CompositionError("IMG taxon location output has duplicate taxon_oids")
    if scalar(connection, "SELECT count(*) FROM img_taxa WHERE marker_group = 'No IMG SSU record'") != 0:
        raise CompositionError("IMG location table contains taxa without IMG SSU records")
    img_source_taxa = scalar(
        connection,
        """SELECT count(DISTINCT taxon_oid)
           FROM read_parquet(?) WHERE reference_source = 'IMG'""",
        [source],
    )
    img_source_taxa_without_location_row = scalar(
        connection,
        """SELECT count(DISTINCT s.taxon_oid)
           FROM read_parquet(?) s
           ANTI JOIN read_parquet(?) l USING (taxon_oid)
           WHERE s.reference_source = 'IMG'""",
        [source, locations],
    )
    if scalar(connection, "SELECT count(*) FROM img_taxa") != img_source_taxa:
        raise CompositionError("IMG taxon denominator does not match the source records")

    connection.execute(
        """CREATE OR REPLACE TEMP TABLE img_sites AS
           SELECT latitude, longitude,
                  count(*) AS taxa,
                  sum(records_16s) AS records_16s,
                  sum(records_18s) AS records_18s,
                  sum(sequences_16s) AS sequences_16s,
                  sum(sequences_18s) AS sequences_18s,
                  CASE
                    WHEN sum(records_16s) > 0 AND sum(records_18s) > 0 THEN '16S + 18S'
                    WHEN sum(records_16s) > 0 THEN '16S only'
                    ELSE '18S only'
                  END AS marker_group,
                  array_to_string(list_slice(list(taxon_oid ORDER BY taxon_oid), 1, 5), ', ') AS taxon_oid_preview
           FROM img_taxa
           WHERE latitude IS NOT NULL AND longitude IS NOT NULL
           GROUP BY latitude, longitude
           ORDER BY latitude, longitude"""
    )

    taxon_output = output_directory / "img_taxon_locations.parquet"
    site_output = output_directory / "img_location_sites.parquet"
    output_directory.mkdir(parents=True, exist_ok=True)
    for table, destination in (("img_taxa", taxon_output), ("img_sites", site_output)):
        temporary = destination.with_name(f".{destination.name}.tmp")
        temporary.unlink(missing_ok=True)
        connection.execute(
            f"COPY (SELECT * FROM {table}) TO ? (FORMAT PARQUET, COMPRESSION ZSTD)",
            [str(temporary)],
        )
        os.replace(temporary, destination)
        destination.chmod(0o644)

    coverage = fetch_dicts(
        connection,
        """SELECT count(*) AS taxa,
                  count(*) FILTER (WHERE latitude IS NOT NULL AND longitude IS NOT NULL) AS located_taxa,
                  count(DISTINCT (latitude, longitude))
                    FILTER (WHERE latitude IS NOT NULL AND longitude IS NOT NULL) AS sites,
                  sum(records_16s + records_18s) AS source_records,
                  sum(records_16s + records_18s)
                    FILTER (WHERE latitude IS NOT NULL AND longitude IS NOT NULL)
                    AS located_source_records
           FROM img_taxa""",
    )[0]
    coverage["location_table_rows"] = scalar(
        connection, "SELECT count(*) FROM read_parquet(?)", [locations]
    )
    coverage["taxa_without_location_row"] = img_source_taxa_without_location_row
    coverage["marker_groups"] = fetch_dicts(
        connection,
        "SELECT marker_group, count(*) AS sites FROM img_sites GROUP BY marker_group ORDER BY marker_group",
    )
    return taxon_output, site_output, coverage


def source_rows(source_config: dict, composition: list[dict]) -> list[dict]:
    totals: dict[tuple[str, str, str], int] = {}
    for row in composition:
        if row["dimension"] == "source_records":
            key = (row["profile"], row["category"], row["marker"])
            totals[key] = totals.get(key, 0) + int(row["count"])

    output: list[dict] = []
    for display in SOURCE_DISPLAY:
        raw = source_config["sources"][display["config_key"]]
        reference_source = display["reference_source"]
        marker = display["marker"]

        def profile_total(profile: str) -> int:
            return sum(
                count
                for (row_profile, row_source, row_marker), count in totals.items()
                if row_profile == profile
                and row_source == reference_source
                and (marker is None or row_marker == marker)
            )

        output.append(
            {
                "source": display["source"],
                "version": raw["version"],
                "content": display["content"],
                "curated_source_records": profile_total("curated"),
                "img_source_records": profile_total("img"),
                "doi": raw.get("doi", raw.get("citation_doi", "")),
                "license": raw.get("license", ""),
                "url": raw["url"],
                "sha256": raw["sha256"],
            }
        )
    return output


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Prepare database-composition documentation data.")
    parser.add_argument("--curated-profile", type=Path, required=True)
    parser.add_argument("--curated-qa", type=Path, required=True)
    parser.add_argument("--img-profile", type=Path, required=True)
    parser.add_argument("--img-qa", type=Path, required=True)
    parser.add_argument("--catalog", type=Path, default=Path("config/database_catalog.json"))
    parser.add_argument("--sources", type=Path, default=Path("config/database_sources.json"))
    parser.add_argument("--output-directory", type=Path, default=Path("docs/data"))
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    catalog = read_json(args.catalog)
    release_state = catalog.get("release_state", "published")
    if release_state not in {"published", "unpublished_candidate"}:
        raise CompositionError(f"invalid database release state: {release_state!r}")
    sources = read_json(args.sources)
    validated = {
        "curated": validate_profile("curated", args.curated_profile, args.curated_qa, catalog),
        "img": validate_profile("img", args.img_profile, args.img_qa, catalog),
    }

    connection = duckdb.connect(":memory:")
    connection.execute("SET threads = 1")
    composition: list[dict] = []
    taxonomy_resolution: list[dict] = []
    taxonomy_lineages: list[dict] = []
    input_tables: list[dict] = []
    for profile in ("curated", "img"):
        manifest, provenance, qa, paths = validated[profile]
        if provenance.get("counts") != qa.get("counts"):
            raise CompositionError(f"provenance and QA counts differ for {profile}")
        composition.extend(
            composition_rows(connection, profile, manifest["version"], paths, qa)
        )
        profile_resolution, profile_lineages = taxonomy_rank_rows(
            connection,
            profile,
            manifest["version"],
            paths["preferred_taxonomy"],
        )
        taxonomy_resolution.extend(profile_resolution)
        taxonomy_lineages.extend(profile_lineages)
        input_tables.extend(
            {
                "profile": profile,
                "table": name,
                "bytes": path.stat().st_size,
                "sha256": sha256(path),
            }
            for name, path in sorted(paths.items())
        )

    output_directory = args.output_directory
    composition_path = output_directory / "database_composition.tsv"
    sources_path = output_directory / "database_sources.tsv"
    resolution_path = output_directory / "database_taxonomy_resolution.tsv"
    lineages_path = output_directory / "database_taxonomy_lineages.tsv"
    atomic_tsv(
        composition_path,
        ("profile", "version", "dimension", "category", "marker", "count"),
        composition,
    )
    atomic_tsv(
        sources_path,
        (
            "source",
            "version",
            "content",
            "curated_source_records",
            "img_source_records",
            "doi",
            "license",
            "url",
            "sha256",
        ),
        source_rows(sources, composition),
    )
    atomic_tsv(
        resolution_path,
        (
            "profile",
            "version",
            "domain",
            "target_rank",
            "total_sequences",
            "resolved_sequences",
            "unresolved_sequences",
            "resolution_percent",
            "named_lineages",
        ),
        taxonomy_resolution,
    )
    atomic_tsv(
        lineages_path,
        ("profile", "version", "domain", "target_rank", "lineage", "count"),
        taxonomy_lineages,
    )

    _, _, _, img_paths = validated["img"]
    taxon_path, site_path, coverage = write_location_tables(
        connection,
        img_paths["source_records"],
        img_paths["img_location"],
        output_directory,
    )
    profiles = {}
    for profile in ("curated", "img"):
        manifest, provenance, qa, _ = validated[profile]
        profiles[profile] = {
            "version": manifest["version"],
            "archive": catalog["profiles"][profile]["archive"],
            "manifest_sha256": qa["profile_manifest_sha256"],
            "counts": qa["counts"],
            "markers": qa["markers"],
            "sources": qa["sources"],
        }
        if profile == "img":
            profiles[profile]["coordinate_coverage"] = coverage

    taxonomy_policy = validated["curated"][1]["release_build"]["taxonomy_policy"]
    if taxonomy_policy != validated["img"][1]["release_build"]["taxonomy_policy"]:
        raise CompositionError("profiles do not share one taxonomy policy")
    build_tree_hashes = {validated[profile][2]["source_tree_sha256"] for profile in validated}
    if len(build_tree_hashes) != 1:
        raise CompositionError("profiles were not built from the same release source tree")
    taxonomy_summary = {
        domain: {
            "target_rank": target_rank,
            "taxonomy_source": taxonomy_source,
            "taxonomy_source_version": next(
                row["version"]
                for row in source_rows(sources, composition)
                if row["source"] == taxonomy_source
            ),
            "taxonomy_position": 2,
        }
        for domain, (target_rank, taxonomy_source) in DOMAIN_RANKS.items()
    }
    taxonomy_summary["Eukaryota"]["normalization"] = (
        "PR2 compartment suffixes are removed from every taxonomy rank during import"
    )
    taxonomy_summary["Eukaryota"]["normalized_suffixes"] = list(PR2_COMPARTMENT_SUFFIXES)
    taxonomy_summary["unresolved_definition"] = (
        "no non-empty second taxonomy segment, or a second segment labeled unclassified "
        "or incertae sedis"
    )
    outputs = [
        composition_path,
        sources_path,
        resolution_path,
        lineages_path,
        taxon_path,
        site_path,
    ]
    release_provenance = {
        "state": release_state,
        "version": catalog["profiles"]["curated"]["version"],
        "zenodo_concept_record_id": catalog["zenodo"]["concept_record_id"],
    }
    if release_state == "published":
        record_id = catalog["zenodo"]["record_id"]
        release_provenance.update(
            {
                "zenodo_record_id": record_id,
                "doi": f"10.5281/zenodo.{record_id}",
                "url": f"https://zenodo.org/records/{record_id}",
            }
        )
    else:
        release_provenance["parent_zenodo_record_id"] = catalog["zenodo"][
            "record_id"
        ]

    provenance_output = {
        "schema_version": 2,
        "database_release": release_provenance,
        "build_source_tree": {
            "sha256": build_tree_hashes.pop(),
            "scope": "fingerprinted database-build code, configuration, environment lock, and notices",
        },
        "profiles": profiles,
        "taxonomy_policy": taxonomy_policy,
        "taxonomy_summary": taxonomy_summary,
        "privacy": {
            "img_location_columns": ["taxon_oid", "latitude", "longitude"],
            "excluded_source_metadata": "all non-allowlisted IMG metadata fields",
        },
        "inputs": input_tables,
        "outputs": [
            {
                "path": str(path.relative_to(output_directory.parent.parent)),
                "bytes": path.stat().st_size,
                "sha256": sha256(path),
            }
            for path in outputs
        ],
    }
    atomic_text(
        output_directory / "database_composition_provenance.json",
        json.dumps(provenance_output, indent=2, sort_keys=True) + "\n",
    )


if __name__ == "__main__":
    main()
