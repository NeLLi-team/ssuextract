#!/usr/bin/env nextflow

/*
========================================================================================
    SSUextract: Small Subunit rRNA Extraction Pipeline
========================================================================================
    Github : https://github.com/NeLLi-team/ssuextract
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
========================================================================================
    VALIDATE INPUTS
========================================================================================
*/

// Show help message
if (params.help) {
    helpMessage()
    exit 0
}

/*
========================================================================================
    CONFIG FILES
========================================================================================
*/

/*
========================================================================================
    IMPORT LOCAL MODULES/SUBWORKFLOWS
========================================================================================
*/

/*
========================================================================================
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
========================================================================================
*/

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

workflow {

    // Get input files
    fna_files = Channel.fromPath("${params.querydir}/*.{fna,fa,fasta}")
        .map { file -> [file.baseName, file] }
    
    // Get covariance models
    cm_models = Channel.fromPath("${params.modeldir}/*.cm")
        .map { file -> [file.baseName, file] }

    // Process headers
    CHECK_FNA_HEADERS(fna_files)
    
    // Create cross product of processed files and models
    processed_files = CHECK_FNA_HEADERS.out
    sample_model_combinations = processed_files.combine(cm_models)
    
    // Run cmsearch
    CMSEARCH(sample_model_combinations)
    
    // Get stats
    GET_CMSTATS(CMSEARCH.out)
    
    // Extract sequences
    EXTRACT_SEQUENCES(
        GET_CMSTATS.out.map { meta, seqmap, seqsumt, fna -> [meta, seqmap, fna] }
    )
    
    // Annotate with BLAST
    BLAST_ANNOTATE(EXTRACT_SEQUENCES.out)
    
    // Merge BLAST results
    MERGE_BLAST_RESULTS(BLAST_ANNOTATE.out.map { it[1] }.collect())
    
    // Process results
    PROCESS_RESULTS(MERGE_BLAST_RESULTS.out)
    
    // Build final summary
    BUILD_SUMMARY(
        PROCESS_RESULTS.out,
        EXTRACT_SEQUENCES.out.map { it[1] }.collect(),
        BLAST_ANNOTATE.out.map { it[1] }.collect(),
        GET_CMSTATS.out.map { it[2] }.collect()
    )
}

/*
========================================================================================
    PROCESSES
========================================================================================
*/

process CHECK_FNA_HEADERS {
    tag "$sample_id"
    publishDir "${params.outdir}/fna", mode: 'copy'
    
    input:
    tuple val(sample_id), path(fna_file)
    
    output:
    tuple val(sample_id), path("${sample_id}_processed.fna")
    
    script:
    """
    # Simple copy for now to test pipeline
    cp ${fna_file} ${sample_id}_processed.fna
    """
}

process CMSEARCH {
    tag "${sample_id}_${model_id}"
    publishDir "${params.outdir}/out", mode: 'copy'
    cpus params.threads_per_job
    
    input:
    tuple val(sample_id), path(fna_file), val(model_id), path(cm_model)
    
    output:
    tuple val("${sample_id}_${model_id}"), path("${sample_id}_${model_id}.out"), path(fna_file), path(cm_model)
    
    script:
    """
    cmsearch --anytrunc --cpu ${task.cpus} -o /dev/null --tblout ${sample_id}_${model_id}.out ${cm_model} ${fna_file}
    """
}

process GET_CMSTATS {
    tag "$meta"
    publishDir "${params.outdir}/stats", mode: 'copy'
    
    input:
    tuple val(meta), path(cmsearch_out), path(fna_file), path(cm_model)
    
    output:
    tuple val(meta), path("${meta}.seqmap"), path("${meta}.seqsumt"), path(fna_file)
    
    script:
    """
    python3 ${projectDir}/scripts/get_cmstats.py ${cmsearch_out} ${meta} ${cm_model} > ${meta}.log
    """
}

process EXTRACT_SEQUENCES {
    tag "$meta"
    publishDir "${params.outdir}/extracted", mode: 'copy'
    
    input:
    tuple val(meta), path(seqmap), path(fna_file)
    
    output:
    tuple val(meta), path("${meta}.fna")
    
    script:
    """
    python3 ${projectDir}/scripts/get_cmsequences.py ${fna_file} ${seqmap} ${meta}.fna ${params.min_extract_length}
    """
}

process BLAST_ANNOTATE {
    tag "$meta"
    publishDir "${params.outdir}/m8", mode: 'copy'
    cpus params.threads_per_job
    
    input:
    tuple val(meta), path(extracted_fna)
    
    output:
    tuple val(meta), path("${meta}.m8")
    
    script:
    db_path = "${projectDir}/${params.database_path}/silva-138-1_pr2-4-12"
    """
    blastn -outfmt 6 -db ${db_path} -query ${extracted_fna} -max_target_seqs 5 -num_threads ${task.cpus} -out ${meta}.m8
    """
}

process MERGE_BLAST_RESULTS {
    publishDir "${params.outdir}/m8", mode: 'copy'
    
    input:
    path(m8_files)
    
    output:
    path("merged.m8")
    
    script:
    """
    cat ${m8_files} > merged.m8
    """
}

process PROCESS_RESULTS {
    publishDir "${params.outdir}", mode: 'copy'
    
    input:
    path(merged_m8)
    
    output:
    path("cmsearch_summary.tab")
    
    script:
    """
    python ${projectDir}/scripts/cmprocessing.py fna ${merged_m8} cmsearch_summary.tab
    """
}

process BUILD_SUMMARY {
    publishDir "${params.outdir}", mode: 'copy'
    
    input:
    path(summary_tab)
    path(extracted_files)
    path(m8_files)
    path(seqsumt_files)
    
    output:
    path("cmsearch_summary.tsv")
    
    script:
    """
    # Create temporary symlink structure for compatibility
    mkdir -p temp_${params.querydir_name}/cmsearch_out
    
    # Create subdirectories and organize files
    mkdir -p temp_${params.querydir_name}/cmsearch_out/extracted
    mkdir -p temp_${params.querydir_name}/cmsearch_out/m8  
    mkdir -p temp_${params.querydir_name}/cmsearch_out/stats
    
    # Copy files to expected locations
    cp *.fna temp_${params.querydir_name}/cmsearch_out/extracted/ 2>/dev/null || true
    cp *.m8 temp_${params.querydir_name}/cmsearch_out/m8/ 2>/dev/null || true
    cp *.seqsumt temp_${params.querydir_name}/cmsearch_out/stats/ 2>/dev/null || true
    cp cmsearch_summary.tab temp_${params.querydir_name}/cmsearch_out/ 2>/dev/null || true
    
    # Create dummy .fna files for each unique sample
    for f in *.fna; do
        if [ -f "\$f" ]; then
            filename=\$(basename "\$f")
            sample=\$(echo "\$filename" | sed 's/_RF[0-9]*\\.fna//' | sed 's/_processed//')
            touch "temp_${params.querydir_name}/\${sample}.fna"
        fi
    done
    
    # Debug: show what we created
    echo "Created dummy files:"
    ls -la temp_${params.querydir_name}/*.fna || echo "No .fna files created"
    echo "Files in extracted:"
    ls -la temp_${params.querydir_name}/cmsearch_out/extracted/ || echo "No extracted files"
    
    # Run the table generation script with absolute paths
    python ${projectDir}/scripts/get_table.py temp_${params.querydir_name} ${projectDir}/${params.modeldir}
    
    # Check if output was created
    if [ -f temp_${params.querydir_name}/cmsearch_out/cmsearch_summary.tsv ]; then
        mv temp_${params.querydir_name}/cmsearch_out/cmsearch_summary.tsv .
    else
        echo "Error: cmsearch_summary.tsv was not created"
        echo "Contents of temp directory:"
        find temp_${params.querydir_name} -type f
        exit 1
    fi
    
    # Cleanup
    rm -rf temp_${params.querydir_name}
    """
}

/*
========================================================================================
    COMPLETION EMAIL AND SUMMARY
========================================================================================
*/

workflow.onComplete {
    println "Pipeline completed at: $workflow.complete"
    println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
    println "Results directory: ${params.outdir}"
}

/*
========================================================================================
    THE END
========================================================================================
*/

def helpMessage() {
    log.info"""
    
    Usage:
    
    The typical command for running the pipeline is as follows:
    
    nextflow run main.nf --querydir data/example --modeldir resources/models
    
    Mandatory arguments:
      --querydir [path]           Path to directory containing input FASTA files (.fna, .fa, .fasta)
      --modeldir [path]           Path to directory containing covariance models (.cm files)
    
    Optional arguments:
      --outdir [path]             Path to output directory (default: results/[querydir_name])
      --min_extract_length [int]  Minimum sequence length for extraction (default: 30)
      --threads_per_job [int]     Number of threads per job (default: 2)
      --database_path [path]      Path to BLAST database directory (default: resources/database)
      --help                      Show this help message and exit
    
    """.stripIndent()
}