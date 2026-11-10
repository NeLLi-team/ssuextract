#!/usr/bin/env nextflow

workflow {
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
    validatePositiveInteger(params.top_hits, 'top_hits')
    validateBoolean(params.tree_classification, 'tree_classification')
    validateMinimumInteger(params.tree_reference_count, 'tree_reference_count', 3)
    validatePositiveInteger(params.tree_assignment_neighbors, 'tree_assignment_neighbors')
    validateFraction(params.tree_trim_gap_fraction, 'tree_trim_gap_fraction')
    if ((params.tree_assignment_neighbors as int) > (params.tree_reference_count as int)) {
        throw new IllegalArgumentException(
            '--tree_assignment_neighbors cannot exceed --tree_reference_count'
        )
    }
    validateIdentifier(params.database_profile.toString(), 'database profile')
    def database_config = loadDatabaseConfig(params.database_path, params.database_profile)
    def runtime_code_digest = runtimeCodeDigest()
    def model_markers = loadModelMarkers(params.model_marker_map)
    def tree_classification = params.tree_classification.toString() == 'true'
    def tree_model_ids = tree_classification \
        ? loadTreeModelIds(model_markers) \
        : [:]

    def query_input = file(params.query, checkIfExists: true)
    def query_files = query_input.isDirectory() \
        ? channel.fromPath("${query_input}/*.{fna,fa,fasta}", checkIfExists: true) \
        : channel.of(validateQueryFasta(query_input))
    def fna_files = query_files
        .map { fasta ->
            def sample_id = fasta.baseName
            validateIdentifier(sample_id, 'sample')
            tuple(sample_id, fasta)
        }

    def cm_models = channel
        .fromPath("${params.modeldir}/*.cm", checkIfExists: true)
        .map { model_file ->
            def model_id = model_file.baseName
            validateIdentifier(model_id, 'model')
            def marker = markerForModel(model_id, model_markers)
            tuple(
                model_id,
                model_file,
                marker,
                databasePrefixForMarker(database_config, marker),
                database_config.taxonomy_file,
                database_config.source_records_file
            )
        }

    def sample_model_combinations = fna_files.combine(cm_models)

    CMSEARCH(sample_model_combinations)
    def cmsearch_files_by_sample = CMSEARCH.out
        .map { sample_id, model_id, cmsearch_out, fna_file, cm_model, marker, db_prefix, taxonomy_file, source_records_file ->
            tuple(sample_id, cmsearch_out)
        }
        .groupTuple()
    RESOLVE_MODEL_HITS(cmsearch_files_by_sample, runtime_code_digest)
    def extraction_inputs = CMSEARCH.out
        .map { sample_id, model_id, cmsearch_out, fna_file, cm_model, marker, db_prefix, taxonomy_file, source_records_file ->
            tuple(sample_id, model_id, fna_file, cm_model, marker, db_prefix, taxonomy_file, source_records_file)
        }
        .combine(RESOLVE_MODEL_HITS.out, by: 0)
    EXTRACT_HITS(extraction_inputs, runtime_code_digest)
    BLAST_ANNOTATE(EXTRACT_HITS.out, database_config, runtime_code_digest)

    def tree_assignment_files
    def tree_neighbor_files
    if (tree_classification) {
        PREPARE_TREE_TASKS(EXTRACT_HITS.out, database_config, runtime_code_digest, tree_model_ids)
        def tree_task_inputs = PREPARE_TREE_TASKS.out.tasks.flatMap { sample_id, model_id, task_directories ->
            def directories = task_directories instanceof Collection \
                ? task_directories \
                : [task_directories]
            directories.collect { task_directory ->
                def task_config = new groovy.json.JsonSlurper().parse(
                    new File(task_directory.toString(), 'task.json')
                )
                def marker = task_config.tree_marker.toString()
                def model = tree_model_ids[marker]
                tuple(
                    sample_id,
                    model_id,
                    task_config.query_key.toString(),
                    task_directory,
                    file(resolveProjectPath("${params.modeldir}/${model}.cm"), checkIfExists: true),
                    databasePrefixForMarker(database_config, marker)
                )
            }
        }
        TREE_CLASSIFY(tree_task_inputs, database_config, runtime_code_digest)
        tree_assignment_files = TREE_CLASSIFY.out
            .map { sample_id, model_id, query_key, tree_directory ->
                file("${tree_directory}/${query_key}.tree_assignment.tsv")
            }
            .mix(PREPARE_TREE_TASKS.out.skipped)
            .collect()
            .map { files ->
                files ?: [file("${projectDir}/config/empty.tree_assignment.tsv")]
            }
            .ifEmpty([file("${projectDir}/config/empty.tree_assignment.tsv")])
        tree_neighbor_files = TREE_CLASSIFY.out
            .map { sample_id, model_id, query_key, tree_directory ->
                file("${tree_directory}/${query_key}.tree_neighbors.tsv")
            }
            .collect()
            .map { files ->
                files ?: [file("${projectDir}/config/empty.tree_neighbors.tsv")]
            }
            .ifEmpty([file("${projectDir}/config/empty.tree_neighbors.tsv")])
    } else {
        tree_assignment_files = channel.value(
            [file("${projectDir}/config/empty.tree_assignment.tsv")]
        )
        tree_neighbor_files = channel.value(
            [file("${projectDir}/config/empty.tree_neighbors.tsv")]
        )
    }

    def summary_files = BLAST_ANNOTATE.out
        .map { sample_id, model_id, m8, summary, top_hits, metadata -> summary }
        .collect()
    def metadata_files = BLAST_ANNOTATE.out
        .map { sample_id, model_id, m8, summary, top_hits, metadata -> metadata }
        .collect()
    def m8_files = BLAST_ANNOTATE.out
        .map { sample_id, model_id, m8, summary, top_hits, metadata -> m8 }
        .collect()
    def top_hit_files = BLAST_ANNOTATE.out
        .map { sample_id, model_id, m8, summary, top_hits, metadata -> top_hits }
        .collect()

    FINALIZE_SUMMARIES(
        summary_files,
        metadata_files,
        m8_files,
        top_hit_files,
        tree_assignment_files,
        tree_neighbor_files,
        tree_classification ? 'tree' : 'blast',
        database_config,
        runtime_code_digest
    )

    workflow.onComplete = {
        println "Pipeline completed at: ${workflow.complete}"
        println "Execution status: ${workflow.success ? 'OK' : 'failed'}"
        println "Results directory: ${params.outdir}"
    }
}


process CMSEARCH {
    tag "${sample_id}_${model_id}"
    publishDir "${params.outdir}/out", mode: 'copy', pattern: '*.out'

    input:
    tuple \
        val(sample_id), \
        path(fna_file), \
        val(model_id), \
        path(cm_model), \
        val(marker), \
        val(db_prefix), \
        val(taxonomy_file), \
        val(source_records_file)

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
        val(source_records_file)

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
    val(runtime_code_digest)

    output:
    tuple val(sample_id), path("${sample_id}.accepted-hits.tsv")

    script:
    def cmsearch_arguments = cmsearch_files
        .collect { cmsearch_file -> "--cmsearch ${shellQuote(cmsearch_file)}" }
        .join(' ')
    """
    # Runtime code SHA-256: ${runtime_code_digest}
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
        val(source_records_file), \
        path(accepted_hits)
    val(runtime_code_digest)

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
        val(source_records_file)

    script:
    """
    # Runtime code SHA-256: ${runtime_code_digest}
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
    publishDir "${params.outdir}/m8", mode: 'copy', pattern: '*.top_hits.tsv'

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
        val(source_records_file)
    val(database_config)
    val(runtime_code_digest)

    output:
    tuple \
        val(sample_id), \
        val(model_id), \
        path("${sample_id}_${model_id}.m8"), \
        path("${sample_id}_${model_id}.summary.tsv"), \
        path("${sample_id}_${model_id}.top_hits.tsv"), \
        path(metadata)

    script:
    def db_prefix_argument = shellQuote(db_prefix)
    def taxonomy_argument = shellQuote(taxonomy_file)
    def source_records_argument = shellQuote(source_records_file)
    def blast_fetch_targets = Math.max(
        (params.max_blast_targets as int) + 1,
        params.top_hits as int
    )
    def calibration_argument = database_config.runtime_calibration \
        ? "--runtime-calibration ${shellQuote(database_config.runtime_calibration)}" \
        : ''
    """
    # Database SHA-256: ${database_config.content_digest}
    # Runtime code SHA-256: ${runtime_code_digest}
    blastn \
        -task blastn \
        -evalue 1e-5 \
        -outfmt '6 std qlen slen' \
        -db ${db_prefix_argument} \
        -query "${extracted_fna}" \
        -max_target_seqs "${blast_fetch_targets}" \
        -max_hsps 1 \
        -num_threads "${task.cpus}" \
        -out "${sample_id}_${model_id}.m8"

    python3 "${projectDir}/scripts/annotate_hits.py" \
        --hits "${hits_table}" \
        --m8 "${sample_id}_${model_id}.m8" \
        --taxonomy-db ${taxonomy_argument} \
        --source-records-db ${source_records_argument} \
        ${calibration_argument} \
        --marker ${shellQuote(marker.toString())} \
        --reference-digest "${database_config.reference_digest}" \
        --query-fasta "${extracted_fna}" \
        --top-hits "${params.top_hits}" \
        --top-hits-output "${sample_id}_${model_id}.top_hits.tsv" \
        --max-targets "${params.max_blast_targets}" \
        --output "${sample_id}_${model_id}.summary.tsv"
    """
}


process PREPARE_TREE_TASKS {
    tag "${sample_id}_${model_id}"

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
        val(source_records_file)
    val(database_config)
    val(runtime_code_digest)
    val(tree_model_ids)

    output:
    tuple \
        val(sample_id), \
        val(model_id), \
        path('tree_inputs/*'), \
        optional: true, \
        emit: tasks
    path("${sample_id}_${model_id}.skipped.tree_assignment.tsv"), emit: skipped

    script:
    def marker_model_arguments = tree_model_ids
        .collect { tree_marker, tree_model ->
            "--marker-model ${shellQuote("${tree_marker}=${tree_model}")}"
        }
        .join(' ')
    def fetch_targets = Math.max(
        Math.max(params.tree_reference_count as int, params.max_blast_targets as int),
        100
    ) + 1
    def database_16s = shellQuote(databasePrefixForMarker(database_config, '16S'))
    def database_18s = shellQuote(databasePrefixForMarker(database_config, '18S'))
    def taxonomy_argument = shellQuote(taxonomy_file)
    def source_records_argument = shellQuote(source_records_file)
    def calibration_argument = database_config.runtime_calibration \
        ? "--runtime-calibration ${shellQuote(database_config.runtime_calibration)}" \
        : ''
    """
    # Database SHA-256: ${database_config.content_digest}
    # Runtime code SHA-256: ${runtime_code_digest}
    blastn \
        -task blastn \
        -evalue 1e-5 \
        -outfmt '6 std qlen slen' \
        -db ${database_16s} \
        -query "${extracted_fna}" \
        -max_target_seqs "${fetch_targets}" \
        -max_hsps 1 \
        -num_threads "${task.cpus}" \
        -out 16S.tree_route.m8

    blastn \
        -task blastn \
        -evalue 1e-5 \
        -outfmt '6 std qlen slen' \
        -db ${database_18s} \
        -query "${extracted_fna}" \
        -max_target_seqs "${fetch_targets}" \
        -max_hsps 1 \
        -num_threads "${task.cpus}" \
        -out 18S.tree_route.m8

    python3 "${projectDir}/scripts/tree_reference_selection.py" prepare \
        --query-fasta "${extracted_fna}" \
        --blast 16S=16S.tree_route.m8 \
        --blast 18S=18S.tree_route.m8 \
        --taxonomy-db ${taxonomy_argument} \
        --source-records-db ${source_records_argument} \
        ${calibration_argument} \
        --reference-digest "${database_config.reference_digest}" \
        --sample "${sample_id}" \
        --detected-model "${model_id}" \
        --detected-marker ${shellQuote(marker.toString())} \
        ${marker_model_arguments} \
        --reference-count "${params.tree_reference_count}" \
        --max-targets "${params.max_blast_targets}" \
        --route-hits 100 \
        --output-directory tree_inputs \
        --skipped-assignments-output "${sample_id}_${model_id}.skipped.tree_assignment.tsv"
    """
}


process TREE_CLASSIFY {
    tag "${sample_id}_${query_key}"
    publishDir path: { "${params.outdir}/phylogeny/${sample_id}/${model_id}" }, mode: 'copy'

    input:
    tuple \
        val(sample_id), \
        val(model_id), \
        val(query_key), \
        path(tree_task, name: 'tree_input'), \
        path(cm_model), \
        val(db_prefix)
    val(database_config)
    val(runtime_code_digest)

    output:
    tuple \
        val(sample_id), \
        val(model_id), \
        val(query_key), \
        path("${query_key}")

    script:
    def db_prefix_argument = shellQuote(db_prefix)
    """
    # Database SHA-256: ${database_config.content_digest}
    # Runtime code SHA-256: ${runtime_code_digest}
    mkdir "${query_key}"
    cp -R "${tree_task}/." "${query_key}/"

    blastdbcmd \
        -db ${db_prefix_argument} \
        -entry_batch "${query_key}/reference_ids.txt" \
        -outfmt %f \
        -out "${query_key}/references.fna"

    python3 "${projectDir}/scripts/tree_reference_selection.py" alignment-input \
        --task-directory "${query_key}" \
        --reference-fasta "${query_key}/references.fna" \
        --output "${query_key}/cmalign_input.fna"

    cmalign \
        --cpu "${task.cpus}" \
        --outformat AFA \
        -o "${query_key}/cmalign.fna" \
        "${cm_model}" \
        "${query_key}/cmalign_input.fna"

    python3 "${projectDir}/scripts/tree_phylogeny.py" trim \
        --input "${query_key}/cmalign.fna" \
        --output "${query_key}/cmalign.trimmed.fna" \
        --qc "${query_key}/alignment_qc.json" \
        --maximum-gap-fraction "${params.tree_trim_gap_fraction}"

    outgroup=\$(awk 'NR > 1 {value=\$1} END {print value}' "${query_key}/references.tsv")
    test -n "\${outgroup}"
    iqtree3 \
        -s "${query_key}/cmalign.trimmed.fna" \
        -st DNA \
        -m GTR+F+R4 \
        -fast \
        -alrt 1000 \
        -keep-ident \
        -o "\${outgroup}" \
        -T 1 \
        -seed 1 \
        -pre "${query_key}/iqtree" \
        -redo \
        -quiet

    python3 "${projectDir}/scripts/tree_phylogeny.py" classify \
        --tree "${query_key}/iqtree.treefile" \
        --references "${query_key}/references.tsv" \
        --task "${query_key}/task.json" \
        --assignment-output "${query_key}/${query_key}.tree_assignment.tsv" \
        --neighbors-output "${query_key}/${query_key}.tree_neighbors.tsv" \
        --assignment-neighbors "${params.tree_assignment_neighbors}" \
        --inference-model GTR+F+R4

    cmalign -h 2>&1 | sed -n '2p' > "${query_key}/tool_versions.txt"
    iqtree3 --version >> "${query_key}/tool_versions.txt"
    """
}


process FINALIZE_SUMMARIES {
    publishDir "${params.outdir}", mode: 'copy', pattern: 'cmsearch_summary.*'
    publishDir "${params.outdir}", mode: 'copy', pattern: 'blast_top_hits.tsv'
    publishDir "${params.outdir}", mode: 'copy', pattern: 'tree_nearest_neighbors.tsv'
    publishDir "${params.outdir}/m8", mode: 'copy', pattern: 'merged.m8'
    publishDir "${params.outdir}", mode: 'copy', pattern: 'run_provenance.json'

    input:
    path(summary_files)
    path(metadata_files)
    path(m8_files)
    path(top_hit_files)
    path(tree_assignment_files)
    path(tree_neighbor_files)
    val(taxonomy_mode)
    val(database_config)
    val(runtime_code_digest)

    output:
    path('cmsearch_summary.tsv')
    path('cmsearch_summary.tab')
    path('blast_top_hits.tsv')
    path('tree_nearest_neighbors.tsv')
    path('merged.m8')
    path('run_provenance.json')

    script:
    def provenance = groovy.json.JsonOutput.prettyPrint(groovy.json.JsonOutput.toJson([
        schema_version: 1,
        application_version: workflow.manifest.version,
        database_profile: params.database_profile,
        database_version: database_config.version,
        database_content_sha256: database_config.content_digest,
        reference_search_contract_sha256: database_config.reference_digest,
        database_manifest_sha256: database_config.manifest_digest,
        runtime_code_sha256: runtime_code_digest,
        assignment_policy: 'runtime_lca_v1',
        max_blast_targets: params.max_blast_targets as int,
        taxonomy_mode: taxonomy_mode,
        runtime_calibration_sha256: database_config.calibration_digest
    ]))
    """
    # Runtime code SHA-256: ${runtime_code_digest}
    printf '%s\\n' ${shellQuote(provenance)} > run_provenance.json
    python3 "${projectDir}/scripts/finalize_summaries.py" \
        --summary-output cmsearch_summary.tsv \
        --category-output cmsearch_summary.tab \
        --top-hits-output blast_top_hits.tsv \
        --taxonomy-mode "${taxonomy_mode}" \
        --tree-neighbor-output tree_nearest_neighbors.tsv \
        --merged-m8-output merged.m8
    """
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
    return new groovy.json.JsonSlurper().parse(path)
}


def loadModelMarkers(pathValue) {
    def config = loadJsonFile(pathValue, 'model-marker map')
    if (config.schema_version != 1 || !(config.models instanceof Map)) {
        throw new IllegalArgumentException(
            "Invalid model-marker map: expected schema_version 1 and a models object"
        )
    }
    config.models.each { model, marker ->
        if (!(marker instanceof String) || !(marker in ['16S', '18S'])) {
            throw new IllegalArgumentException(
                "Invalid model-marker map: model '${model}' must map to 16S or 18S"
            )
        }
    }
    return config.models
}


def loadTreeModelIds(modelMarkers) {
    def byMarker = modelMarkers.groupBy { model, marker -> marker.toString() }
    def required = ['16S', '18S']
    required.each { marker ->
        def models = byMarker[marker]?.keySet()?.toList()?.sort() ?: []
        if (models.size() != 1) {
            throw new IllegalArgumentException(
                "--tree_classification requires exactly one covariance model " +
                "mapped to ${marker}; found ${models.size()}"
            )
        }
    }
    return required.collectEntries { marker ->
        [(marker): byMarker[marker].keySet().first().toString()]
    }
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
        throw new IllegalArgumentException(
            "Database profile '${profile}' is not installed under ${root}. " +
            "Run 'pixi run setup --database_profile ${profile}'."
        )
    }

    def manifest = new groovy.json.JsonSlurper().parse(manifestFile)
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
        if (!(marker in ['16S', '18S'])) {
            throw new IllegalArgumentException(
                'Unsupported database marker: expected 16S or 18S'
            )
        }
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
    def sourceRecordsRelative = manifest.taxonomy_database?.source_records
    if (!sourceRecordsRelative) {
        throw new IllegalArgumentException(
            "Invalid database manifest ${manifestFile}: " +
            "missing taxonomy_database.source_records"
        )
    }
    def sourceRecordsFile = resolveContainedPath(
        profileDir,
        sourceRecordsRelative.toString(),
        'source-record Parquet'
    )
    if (!sourceRecordsFile.isFile() || sourceRecordsFile.length() == 0) {
        throw new IllegalArgumentException(
            "Missing source-record Parquet: ${sourceRecordsFile}"
        )
    }
    def calibrationFile = manifest.runtime_calibration \
        ? resolveContainedPath(profileDir, manifest.runtime_calibration.toString(), 'runtime calibration') \
        : null
    def referenceFiles = prefixes.values().collectMany { prefix -> blastFiles(prefix) }
    referenceFiles.add(taxonomyFile)
    def runtimeFiles = referenceFiles + [sourceRecordsFile]
    if (calibrationFile) {
        runtimeFiles.add(calibrationFile)
    }
    def artifacts = (manifest.artifacts ?: []).collectEntries { artifact ->
        [(artifact.path.toString()): artifact]
    }
    def hashes = runtimeFiles.unique().collectEntries { path ->
        def relative = profileDir.canonicalFile.toPath().relativize(path.canonicalFile.toPath()).toString()
        def artifact = artifacts[relative]
        if (!artifact || !path.isFile() || path.length() != artifact.bytes) {
            throw new IllegalArgumentException("Missing or changed database artifact: ${path}")
        }
        def digest = sha256File(path)
        if (digest != artifact.sha256) {
            throw new IllegalArgumentException("Database artifact checksum mismatch: ${path}")
        }
        [(relative): digest]
    }
    def referenceHashes = referenceFiles.collectEntries { path ->
        def relative = profileDir.canonicalFile.toPath().relativize(path.canonicalFile.toPath()).toString()
        [(relative): hashes[relative]]
    }
    def manifestDigest = sha256File(manifestFile)
    def calibrationRelative = calibrationFile \
        ? profileDir.canonicalFile.toPath().relativize(calibrationFile.toPath()).toString() \
        : null
    return [
        prefixes: prefixes,
        taxonomy_file: taxonomyFile.toString(),
        source_records_file: sourceRecordsFile.toString(),
        runtime_calibration: calibrationFile?.toString(),
        calibration_digest: calibrationRelative ? hashes[calibrationRelative] : null,
        manifest_digest: manifestDigest,
        content_digest: contentDigest(hashes + ['manifest.json': manifestDigest]),
        reference_digest: contentDigest(referenceHashes),
        version: manifest.version
    ]
}


def sha256File(path) {
    def digest = java.security.MessageDigest.getInstance('SHA-256')
    path.withInputStream { input ->
        new java.security.DigestInputStream(input, digest)
            .transferTo(java.io.OutputStream.nullOutputStream())
    }
    digest.digest().encodeHex().toString()
}


def contentDigest(hashes) {
    def text = hashes.keySet().sort().collect { name -> "${name}\t${hashes[name]}\n" }.join('')
    java.security.MessageDigest.getInstance('SHA-256').digest(text.getBytes('UTF-8')).encodeHex().toString()
}


def blastFiles(prefixValue) {
    def prefix = new File(prefixValue.toString()).canonicalFile
    (prefix.parentFile.listFiles() ?: []).findAll { candidate ->
        candidate.isFile() && candidate.name.startsWith(prefix.name + '.')
    }.sort { candidate -> candidate.name }
}


def runtimeCodeDigest() {
    def names = [
        'resolve_model_hits.py', 'hit_processing.py', 'extract_hits.py',
        'annotate_hits.py', 'top_hit_reporting.py', 'runtime_taxonomy.py',
        'taxonomy_utils.py', 'tree_reference_selection.py', 'tree_phylogeny.py',
        'tree_schema.py', 'finalize_summaries.py'
    ]
    contentDigest(names.collectEntries { name ->
        [(name): sha256File(new File("${projectDir}/scripts/${name}"))]
    })
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
    def blastdbcmd = new ProcessBuilder('blastdbcmd', '-db', prefix, '-info')
        .redirectErrorStream(true)
        .start()
    def output = blastdbcmd.inputStream.text
    if (blastdbcmd.waitFor() != 0) {
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
    } catch (NumberFormatException _ignored) {
        throw new IllegalArgumentException("--${name} must be an integer")
    }
}


def validateMinimumInteger(value, name, minimum) {
    try {
        if ((value as int) < minimum) {
            throw new IllegalArgumentException(
                "--${name} must be at least ${minimum}"
            )
        }
    } catch (NumberFormatException _ignored) {
        throw new IllegalArgumentException("--${name} must be an integer")
    }
}


def validateNonNegativeInteger(value, name) {
    try {
        if ((value as int) < 0) {
            throw new IllegalArgumentException("--${name} must be non-negative")
        }
    } catch (NumberFormatException _ignored) {
        throw new IllegalArgumentException("--${name} must be an integer")
    }
}


def validateQueryFasta(path) {
    if (!path.isFile()) {
        throw new IllegalArgumentException(
            "--query must be a FASTA file or directory: ${path}"
        )
    }
    if (!(path.name ==~ /.+\.(fna|fa|fasta)/)) {
        throw new IllegalArgumentException(
            "--query file must end in .fna, .fa, or .fasta: ${path}"
        )
    }
    path
}


def validateBoolean(value, name) {
    if (!(value.toString() in ['true', 'false'])) {
        throw new IllegalArgumentException(
            "--${name} must be supplied as a bare flag"
        )
    }
}


def validateFraction(value, name) {
    try {
        def parsed = value as double
        if (parsed < 0 || parsed >= 1) {
            throw new IllegalArgumentException(
                "--${name} must be at least 0 and less than 1"
            )
        }
    } catch (NumberFormatException _ignored) {
        throw new IllegalArgumentException("--${name} must be numeric")
    }
}


def helpMessage() {
    log.info """
    Usage:

      pixi run ssuextract -q data/example --modeldir resources/models
      nextflow run main.nf --query data/example --modeldir resources/models

    Input arguments:
      --query [path]              FASTA file or directory of .fna, .fa, or .fasta files
                                  (default: data/example)
      -q [path]                   Same input through the Pixi wrapper
      --modeldir [path]           Directory containing covariance models (.cm)
                                  (default: resources/models)

    Optional arguments:
      --outdir [path]             Output directory (default: results/[input_name])
      --min_extract_length [int]  Minimum union of hit spans in nucleotides (default: 500)
      --threads_per_job [int]     Threads per Infernal, BLAST, or cmalign task (default: 2)
      --max_blast_targets [int]   Assignment candidate limit; overflow abstains (default: 500)
      --top_hits [int]            Overall BLAST hits; assignment evidence is retained (default: 5)
      --tree_classification       Classify with a query-neighbor SSU tree
      --tree_reference_count [n]  BLAST references per tree (default: 100)
      --tree_assignment_neighbors [n]
                                 Nearest reference count for tree evidence (default: 5)
      --tree_trim_gap_fraction [n]
                                 Remove columns above this gap fraction (default: 0.9)
      --database_path [path]      BLAST database directory (default: resources/database)
      --database_profile [name]   Database profile: curated or img (default: curated)
      --model_marker_map [path]   JSON mapping models to 16S rRNA gene or 18S rRNA gene markers
      --max_cpus [int]            Maximum CPUs for one task (default: 16)
      --max_memory [size]         Maximum memory for one task (default: 128.GB)
      --max_time [duration]       Maximum wall time for one task (default: 240.h)
      --version                   Print the SSUextract version
      --help                      Print this help message
    """.stripIndent()
}
