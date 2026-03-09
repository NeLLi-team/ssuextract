#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
PROJECT_DIR=$(cd -- "${SCRIPT_DIR}/.." && pwd)
CONFIG_FILE="${PROJECT_DIR}/config/local.config"
DEFAULT_DB_DIR="${PROJECT_DIR}/resources/database"
DB_BASE_URL="https://portal.nersc.gov/cfs/nelli/ssuextract_db"
DB_PREFIX="silva-138-1_pr2-4-12"
SMOKE_FASTA="${PROJECT_DIR}/data/example/LKH565_P11_Ci.fna"
DB_SUFFIXES=(fasta nhr nin nsq)

main() {
    local command="${1:-run}"
    shift || true

    case "${command}" in
        setup) setup_database "$@" ;;
        run) run_pipeline "$@" ;;
        smoke) run_smoke "$@" ;;
        *)
            echo "Usage: scripts/pipeline_cli.sh {setup|run|smoke} [nextflow args]" >&2
            exit 1
            ;;
    esac
}

has_help_flag() {
    local arg
    for arg in "$@"; do
        [[ "${arg}" == "--help" ]] && return 0
    done
    return 1
}

cli_has_database_path() {
    local arg
    for arg in "$@"; do
        [[ "${arg}" == --database_path=* ]] && return 0
    done

    local expect_value=0
    for arg in "$@"; do
        if [[ ${expect_value} -eq 1 ]]; then
            return 0
        fi
        [[ "${arg}" == "--database_path" ]] && expect_value=1
    done

    return 1
}

read_cli_database_path() {
    local arg
    for arg in "$@"; do
        case "${arg}" in
            --database_path=*)
                printf '%s\n' "${arg#*=}"
                return 0
                ;;
        esac
    done

    local expect_value=0
    for arg in "$@"; do
        if [[ ${expect_value} -eq 1 ]]; then
            printf '%s\n' "${arg}"
            return 0
        fi
        [[ "${arg}" == "--database_path" ]] && expect_value=1
    done

    return 1
}

read_config_database_path() {
    [[ -f "${CONFIG_FILE}" ]] || return 1
    awk -F"'" '/database_path/ {print $2; exit}' "${CONFIG_FILE}"
}

write_config_database_path() {
    local db_dir="$1"
    cat > "${CONFIG_FILE}" <<EOF
params {
    database_path = '${db_dir}'
}
EOF
}

prompt_database_path() {
    if [[ -t 0 ]]; then
        printf 'Database path [%s]: ' "${DEFAULT_DB_DIR}" >&2
        read -r user_path
        printf '%s\n' "${user_path:-${DEFAULT_DB_DIR}}"
        return
    fi

    printf 'Database path not configured, using default: %s\n' "${DEFAULT_DB_DIR}" >&2
    printf '%s\n' "${DEFAULT_DB_DIR}"
}

resolve_database_path() {
    local cli_path=""
    local config_path=""
    local resolved_path=""

    if cli_path=$(read_cli_database_path "$@" 2>/dev/null); then
        printf '%s\n' "${cli_path}"
        return
    fi

    if config_path=$(read_config_database_path 2>/dev/null); then
        printf '%s\n' "${config_path}"
        return
    fi

    resolved_path=$(prompt_database_path)
    write_config_database_path "${resolved_path}"
    printf '%s\n' "${resolved_path}"
}

database_ready() {
    local db_dir="$1"
    local suffix

    for suffix in "${DB_SUFFIXES[@]}"; do
        [[ -f "${db_dir}/${DB_PREFIX}.${suffix}" ]] || return 1
    done

    return 0
}

download_database() {
    local db_dir="$1"
    local suffix=""
    local file=""
    local url=""

    mkdir -p "${db_dir}"
    for suffix in "${DB_SUFFIXES[@]}"; do
        file="${DB_PREFIX}.${suffix}"
        url="${DB_BASE_URL}/${file}"
        printf 'Downloading %s\n' "${file}" >&2
        wget -q --show-progress -O "${db_dir}/${file}" "${url}"
    done
}

ensure_database() {
    local db_dir="$1"

    if database_ready "${db_dir}"; then
        printf 'Database ready at %s\n' "${db_dir}" >&2
        return
    fi

    printf 'Database files not found in %s\n' "${db_dir}" >&2
    download_database "${db_dir}"
    printf 'Database ready at %s\n' "${db_dir}" >&2
}

setup_database() {
    local db_dir=""
    db_dir=$(resolve_database_path "$@")
    ensure_database "${db_dir}"
}

run_pipeline() {
    local db_dir=""

    if has_help_flag "$@"; then
        nextflow run "${PROJECT_DIR}/main.nf" "$@"
        return
    fi

    db_dir=$(resolve_database_path "$@")
    ensure_database "${db_dir}"

    if cli_has_database_path "$@"; then
        nextflow run "${PROJECT_DIR}/main.nf" "$@"
        return
    fi

    nextflow run "${PROJECT_DIR}/main.nf" --database_path "${db_dir}" "$@"
}

run_smoke() {
    local smoke_dir=""
    smoke_dir=$(mktemp -d "${TMPDIR:-/tmp}/ssuextract-smoke.XXXXXX")
    trap "rm -rf '${smoke_dir}'" EXIT
    cp "${SMOKE_FASTA}" "${smoke_dir}/"
    run_pipeline --querydir "${smoke_dir}" --outdir "${PROJECT_DIR}/results/smoke" --threads_per_job 1 "$@"
}

main "$@"
