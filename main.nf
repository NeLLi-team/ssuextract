#!/usr/bin/env nextflow

import groovy.json.JsonSlurper

nextflow.enable.dsl = 2


if (params.help) {
    helpMessage()
    exit 0
}

if (params.version) {
    log.info "SSUextract ${workflow.manifest.version}"
    exit 0
}

validateNonNegativeInteger(params.min_extract_length, 'min_extract_length')
validatePositiveInteger(params.threads_per_job, 'threads_per_job')
validatePositiveInteger(params.max_blast_targets, 'max_blast_targets')
validateIdentifier(params.database_profile.toString(), 'database profile')
database_config = loadDatabaseConfig(params.database_path, params.database_profile)
model_markers = loadModelMarkers(params.model_marker_map)


workflow {
    fna_files = Channel
        .fromPath("${params.querydir}/*.{fna,fa,fasta}", checkIfExists: true)
        .map { file ->
            sample_id = file.baseName
            validateIdentifier(sample_id, 'sample')
            tuple(sample_id, file)
        }

    cm_models = Channel
        .fromPath("${params.modeldir}/*.cm", checkIfExists: true)
        .map { file ->
            model_id = file.baseName
            validateIdentifier(model_id, 'model')
            marker = markerForModel(model_id, model_markers)
            db_prefix = databasePrefixForMarker(database_config, marker)
            tuple(
                model_id,
                file,
                marker,
                db_prefix,
                database_config.taxonomy_file ?: '',
                database_config.legacy
            )
        }

    sample_model_combinations = fna_files.combine(cm_models)

    CMSEARCH(sample_model_combinations)
    cmsearch_files_by_sample = CMSEARCH.out
        .map { sample_id, model_id, cmsearch_out, fna_file, cm_model, marker, db_prefix, taxonomy_file, legacy_database ->
            tuple(sample_id, cmsearch_out)
        }
        .groupTuple()
    RESOLVE_MODEL_HITS(cmsearch_files_by_sample)
    extraction_inputs = CMSEARCH.out
        .map { sample_id, model_id, cmsearch_out, fna_file, cm_model, marker, db_prefix, taxonomy_file, legacy_database ->
            tuple(sample_id, model_id, fna_file, cm_model, marker, db_prefix, taxonomy_file, legacy_database)
        }
        .combine(RESOLVE_MODEL_HITS.out, by: 0)
    EXTRACT_HITS(extraction_inputs)
    BLAST_ANNOTATE(EXTRACT_HITS.out)

    summary_files = BLAST_ANNOTATE.out
        .map { sample_id, model_id, m8, summary, metadata -> summary }
        .collect()
    metadata_files = BLAST_ANNOTATE.out
        .map { sample_id, model_id, m8, summary, metadata -> metadata }
        .collect()
    m8_files = BLAST_ANNOTATE.out
        .map { sample_id, model_id, m8, summary, metadata -> m8 }
        .collect()

    FINALIZE_SUMMARIES(summary_files, metadata_files, m8_files)
}


process CMSEARCH {
    tag "${sample_id}_${model_id}"
    publishDir "${params.outdir}/out", mode: 'copy', pattern: '*.out'
    cpus params.threads_per_job

    input:
    tuple \
        val(sample_id), \
        path(fna_file), \
        val(model_id), \
        path(cm_model), \
        val(marker), \
        val(db_prefix), \
        val(taxonomy_file), \
        val(legacy_database)

    output:
    tuple \
        val(sample_id), \
        val(model_id), \
        path("${sample_id}_${model_id}.out"), \
        path(fna_file), \
        path(cm_model), \
        val(marker), \
        val(db_prefix), \
        val(taxonomy_file), \
        val(legacy_database)

    script:
    """
    cmsearch \
        --anytrunc \
        --cpu "${task.cpus}" \
        -o /dev/null \
        --tblout "${sample_id}_${model_id}.out" \
        "${cm_model}" \
        "${fna_file}"
    """
}


process RESOLVE_MODEL_HITS {
    tag "${sample_id}"

    input:
    tuple val(sample_id), path(cmsearch_files)

    output:
    tuple val(sample_id), path("${sample_id}.accepted-hits.tsv")

    script:
    cmsearch_arguments = cmsearch_files
        .collect { "--cmsearch ${shellQuote(it)}" }
        .join(' ')
    """
    python3 "${projectDir}/scripts/resolve_model_hits.py" \
        ${cmsearch_arguments} \
        --output "${sample_id}.accepted-hits.tsv"
    """
}


process EXTRACT_HITS {
    tag "${sample_id}_${model_id}"
    publishDir "${params.outdir}/extracted", mode: 'copy', pattern: '*.fna'
    publishDir "${params.outdir}/stats", mode: 'copy', pattern: '*.hits.tsv'

    input:
    tuple \
        val(sample_id), \
        val(model_id), \
        path(fna_file), \
        path(cm_model), \
        val(marker), \
        val(db_prefix), \
        val(taxonomy_file), \
        val(legacy_database), \
        path(accepted_hits)

    output:
    tuple \
        val(sample_id), \
        val(model_id), \
        path("${sample_id}_${model_id}.fna"), \
        path("${sample_id}_${model_id}.hits.tsv"), \
        path("${sample_id}_${model_id}.meta.tsv"), \
        val(marker), \
        val(db_prefix), \
        val(taxonomy_file), \
        val(legacy_database)

    script:
    """
    python3 "${projectDir}/scripts/extract_hits.py" \
        --model-file "${cm_model}" \
        --fasta "${fna_file}" \
        --sample "${sample_id}" \
        --model "${model_id}" \
        --accepted-hits "${accepted_hits}" \
        --minimum-length "${params.min_extract_length}" \
        --fasta-output "${sample_id}_${model_id}.fna" \
        --hits-output "${sample_id}_${model_id}.hits.tsv" \
        --metadata-output "${sample_id}_${model_id}.meta.tsv"
    """
}


process BLAST_ANNOTATE {
    tag "${sample_id}_${model_id}"
    publishDir "${params.outdir}/m8", mode: 'copy', pattern: '*.m8'
    cpus params.threads_per_job

    input:
    tuple \
        val(sample_id), \
        val(model_id), \
        path(extracted_fna), \
        path(hits_table), \
        path(metadata), \
        val(marker), \
        val(db_prefix), \
        val(taxonomy_file), \
        val(legacy_database)

    output:
    tuple \
        val(sample_id), \
        val(model_id), \
        path("${sample_id}_${model_id}.m8"), \
        path("${sample_id}_${model_id}.summary.tsv"), \
        path(metadata)

    script:
    db_prefix_argument = shellQuote(db_prefix)
    taxonomy_argument = legacy_database \
        ? '' \
        : "--taxonomy-db ${shellQuote(taxonomy_file)}"
    blast_fetch_targets = (params.max_blast_targets as int) + 1
    """
    blastn \
        -outfmt 6 \
        -db ${db_prefix_argument} \
        -query "${extracted_fna}" \
        -max_target_seqs "${blast_fetch_targets}" \
        -max_hsps 1 \
        -num_threads "${task.cpus}" \
        -out "${sample_id}_${model_id}.m8"

    python3 "${projectDir}/scripts/annotate_hits.py" \
        --hits "${hits_table}" \
        --m8 "${sample_id}_${model_id}.m8" \
        ${taxonomy_argument} \
        --max-targets "${params.max_blast_targets}" \
        --output "${sample_id}_${model_id}.summary.tsv"
    """
}


process FINALIZE_SUMMARIES {
    publishDir "${params.outdir}", mode: 'copy', pattern: 'cmsearch_summary.*'
    publishDir "${params.outdir}/m8", mode: 'copy', pattern: 'merged.m8'

    input:
    path(summary_files)
    path(metadata_files)
    path(m8_files)

    output:
    path('cmsearch_summary.tsv')
    path('cmsearch_summary.tab')
    path('merged.m8')

    script:
    """
    python3 "${projectDir}/scripts/finalize_summaries.py" \
        --summary-output cmsearch_summary.tsv \
        --category-output cmsearch_summary.tab \
        --merged-m8-output merged.m8
    """
}


workflow.onComplete {
    println "Pipeline completed at: $workflow.complete"
    println "Execution status: ${workflow.success ? 'OK' : 'failed'}"
    println "Results directory: ${params.outdir}"
}


def resolveProjectPath(value) {
    def candidate = new File(value.toString())
    return candidate.isAbsolute() ? candidate.toString() : "${projectDir}/${value}"
}


def loadJsonFile(pathValue, label) {
    def path = new File(resolveProjectPath(pathValue))
    if (!path.isFile()) {
        throw new IllegalArgumentException("Missing ${label}: ${path}")
    }
    return new JsonSlurper().parse(path)
}


def loadModelMarkers(pathValue) {
    def config = loadJsonFile(pathValue, 'model-marker map')
    if (config.schema_version != 1 || !(config.models instanceof Map)) {
        throw new IllegalArgumentException(
            "Invalid model-marker map: expected schema_version 1 and a models object"
        )
    }
    return config.models
}


def markerForModel(modelId, modelMarkers) {
    def marker = modelMarkers[modelId]
    if (!marker) {
        throw new IllegalArgumentException(
            "No database marker is configured for model '${modelId}'. " +
            "Add it to --model_marker_map."
        )
    }
    return marker.toString()
}


def loadDatabaseConfig(databasePath, profile) {
    def root = new File(resolveProjectPath(databasePath))
    def profileDir = new File(root, profile.toString())
    def manifestFile = new File(profileDir, 'manifest.json')
    if (!manifestFile.isFile()) {
        if (profile.toString() != 'curated') {
            throw new IllegalArgumentException(
                "Database profile '${profile}' is not installed under ${root}. " +
                "Run 'pixi run setup -- --database_profile ${profile}'."
            )
        }
        def legacyPrefix = new File(root, 'silva-138-1_pr2-4-12').toString()
        ['nhr', 'nin', 'nsq'].each { suffix ->
            if (!new File("${legacyPrefix}.${suffix}").isFile()) {
                throw new IllegalArgumentException(
                    "Database profile '${profile}' is not installed under ${root}. " +
                    "Run 'pixi run setup -- --database_profile ${profile}'."
                )
            }
        }
        log.warn 'Using deprecated legacy SILVA 138.1/PR2 4.12 database layout.'
        return [
            legacy: true,
            prefixes: ['16S': legacyPrefix, '18S': legacyPrefix],
            taxonomy_file: null
        ]
    }

    def manifest = new JsonSlurper().parse(manifestFile)
    if (manifest.schema_version != 1 || manifest.profile != profile) {
        throw new IllegalArgumentException(
            "Invalid database manifest ${manifestFile}: schema/profile mismatch"
        )
    }
    if (!(manifest.blast_databases instanceof Map) || !manifest.blast_databases) {
        throw new IllegalArgumentException(
            "Invalid database manifest ${manifestFile}: no BLAST databases"
        )
    }
    def prefixes = manifest.blast_databases.collectEntries { marker, database ->
        def prefix = resolveContainedPath(
            profileDir,
            database.prefix.toString(),
            "BLAST prefix for ${marker}"
        )
        validateBlastPrefix(prefix.toString())
        [marker.toString(), prefix.toString()]
    }
    def taxonomyRelative = manifest.taxonomy_database?.preferred
    if (!taxonomyRelative) {
        throw new IllegalArgumentException(
            "Invalid database manifest ${manifestFile}: missing taxonomy_database.preferred"
        )
    }
    def taxonomyFile = resolveContainedPath(
        profileDir,
        taxonomyRelative.toString(),
        'preferred-taxonomy Parquet'
    )
    if (!taxonomyFile.isFile() || taxonomyFile.length() == 0) {
        throw new IllegalArgumentException(
            "Missing preferred-taxonomy Parquet: ${taxonomyFile}"
        )
    }
    return [legacy: false, prefixes: prefixes, taxonomy_file: taxonomyFile.toString()]
}


def resolveContainedPath(root, relativePath, label) {
    def rootPath = root.canonicalFile.toPath()
    def candidate = new File(root, relativePath).canonicalFile
    if (!candidate.toPath().startsWith(rootPath)) {
        throw new IllegalArgumentException("${label} escapes database profile: ${relativePath}")
    }
    return candidate
}


def databasePrefixForMarker(databaseConfig, marker) {
    def prefix = databaseConfig.prefixes[marker]
    if (!prefix) {
        throw new IllegalArgumentException(
            "Database profile has no BLAST database for marker '${marker}'"
        )
    }
    return prefix.toString()
}


def validateBlastPrefix(prefix) {
    def process = new ProcessBuilder('blastdbcmd', '-db', prefix, '-info')
        .redirectErrorStream(true)
        .start()
    def output = process.inputStream.text
    if (process.waitFor() != 0) {
        throw new IllegalArgumentException(
            "Invalid BLAST database prefix '${prefix}': ${output.trim()}"
        )
    }
}


def shellQuote(value) {
    return "'" + value.toString().replace("'", "'\"'\"'") + "'"
}


def validateIdentifier(identifier, kind) {
    if (!(identifier ==~ /[A-Za-z0-9][A-Za-z0-9._-]*/)) {
        throw new IllegalArgumentException(
            "Invalid ${kind} identifier '${identifier}'. " +
            "Use letters, numbers, '.', '_', or '-'."
        )
    }
}


def validatePositiveInteger(value, name) {
    try {
        if ((value as int) < 1) {
            throw new IllegalArgumentException("--${name} must be at least 1")
        }
    } catch (NumberFormatException ignored) {
        throw new IllegalArgumentException("--${name} must be an integer")
    }
}


def validateNonNegativeInteger(value, name) {
    try {
        if ((value as int) < 0) {
            throw new IllegalArgumentException("--${name} must be non-negative")
        }
    } catch (NumberFormatException ignored) {
        throw new IllegalArgumentException("--${name} must be an integer")
    }
}


def helpMessage() {
    log.info """
    Usage:

      nextflow run main.nf --querydir data/example --modeldir resources/models

    Mandatory arguments:
      --querydir [path]           Directory containing FASTA files (.fna, .fa, .fasta)
      --modeldir [path]           Directory containing covariance models (.cm)

    Optional arguments:
      --outdir [path]             Output directory (default: results/[querydir_name])
      --min_extract_length [int]  Minimum extracted sequence length (default: 500)
      --threads_per_job [int]     Threads per cmsearch or BLAST task (default: 2)
      --max_blast_targets [int]   BLAST subjects; ties at limit back off to domain (default: 500)
      --database_path [path]      BLAST database directory (default: resources/database)
      --database_profile [name]   Database profile: curated or img (default: curated)
      --model_marker_map [path]   JSON mapping covariance models to 16S/18S markers
      --version                   Print the SSUextract version
      --help                      Print this help message
    """.stripIndent()
}
