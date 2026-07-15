#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
PROJECT_DIR=$(cd -- "${SCRIPT_DIR}/.." && pwd)
CONFIG_FILE="${PROJECT_DIR}/config/local.config"
DEFAULT_DB_DIR="${PROJECT_DIR}/resources/database"
DEFAULT_DB_PROFILE="curated"
DATABASE_MANAGER="${PROJECT_DIR}/scripts/database_manager.py"
SMOKE_FASTA="${PROJECT_DIR}/data/example/LKH565_P11_Ci.fna"

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

has_information_flag() {
    local arg
    for arg in "$@"; do
        [[ "${arg}" == "--help" || "${arg}" == "--version" ]] && return 0
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

cli_has_database_profile() {
    local arg
    for arg in "$@"; do
        [[ "${arg}" == --database_profile=* ]] && return 0
    done

    local expect_value=0
    for arg in "$@"; do
        if [[ ${expect_value} -eq 1 ]]; then
            return 0
        fi
        [[ "${arg}" == "--database_profile" ]] && expect_value=1
    done
    return 1
}

read_cli_database_profile() {
    local arg
    for arg in "$@"; do
        case "${arg}" in
            --database_profile=*)
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
        [[ "${arg}" == "--database_profile" ]] && expect_value=1
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
    if grep -Eq '^[[:space:]]*database_path[[:space:]]*=' "${CONFIG_FILE}"; then
        awk -F"'" '/^[[:space:]]*database_path[[:space:]]*=/ {print $2; exit}' "${CONFIG_FILE}"
    else
        sed -n '1p' "${CONFIG_FILE}"
    fi
}

write_config_database_path() {
    local db_dir="$1"
    local temp_file="${CONFIG_FILE}.tmp.$$"
    printf '%s\n' "${db_dir}" > "${temp_file}"
    mv "${temp_file}" "${CONFIG_FILE}"
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
    local profile="$2"
    if python3 "${DATABASE_MANAGER}" validate \
        --root "${db_dir}" \
        --profile "${profile}" \
        >/dev/null 2>&1; then
        return 0
    fi

    # The resolver is needed only for the deprecated unprofiled curated layout.
    python3 "${DATABASE_MANAGER}" resolve \
        --root "${db_dir}" \
        --profile "${profile}" \
        --model RF00177 \
        >/dev/null 2>&1
}

ensure_database() {
    local db_dir="$1"
    local profile="$2"

    if database_ready "${db_dir}" "${profile}"; then
        printf 'Database profile %s ready at %s\n' "${profile}" "${db_dir}" >&2
        return
    fi

    printf 'Database profile %s not found in %s\n' "${profile}" "${db_dir}" >&2
    python3 "${DATABASE_MANAGER}" install \
        --root "${db_dir}" \
        --profile "${profile}"
    if ! database_ready "${db_dir}" "${profile}"; then
        printf 'Installed database profile failed validation: %s/%s\n' \
            "${db_dir}" "${profile}" >&2
        return 1
    fi
    printf 'Database profile %s ready at %s\n' "${profile}" "${db_dir}" >&2
}

setup_database() {
    local db_dir=""
    local profile="${DEFAULT_DB_PROFILE}"
    db_dir=$(resolve_database_path "$@")
    profile=$(read_cli_database_profile "$@" 2>/dev/null || printf '%s\n' "${profile}")
    ensure_database "${db_dir}" "${profile}"
}

run_pipeline() {
    local db_dir=""
    local profile="${DEFAULT_DB_PROFILE}"
    local nextflow_args=()

    if has_information_flag "$@"; then
        nextflow run "${PROJECT_DIR}/main.nf" "$@"
        return
    fi

    db_dir=$(resolve_database_path "$@")
    profile=$(read_cli_database_profile "$@" 2>/dev/null || printf '%s\n' "${profile}")
    ensure_database "${db_dir}" "${profile}"

    if cli_has_database_path "$@" && cli_has_database_profile "$@"; then
        nextflow run "${PROJECT_DIR}/main.nf" "$@"
        return
    fi

    nextflow_args=(nextflow run "${PROJECT_DIR}/main.nf")
    if ! cli_has_database_path "$@"; then
        nextflow_args+=(--database_path "${db_dir}")
    fi
    if ! cli_has_database_profile "$@"; then
        nextflow_args+=(--database_profile "${profile}")
    fi
    "${nextflow_args[@]}" "$@"
}

run_smoke() {
    local smoke_dir=""
    smoke_dir=$(mktemp -d "${TMPDIR:-/tmp}/ssuextract-smoke.XXXXXX")
    trap "rm -rf '${smoke_dir}'" EXIT
    cp "${SMOKE_FASTA}" "${smoke_dir}/"
    run_pipeline --querydir "${smoke_dir}" --outdir "${PROJECT_DIR}/results/smoke" --threads_per_job 1 "$@"
}

if [[ "${BASH_SOURCE[0]}" == "$0" ]]; then
    main "$@"
fi
