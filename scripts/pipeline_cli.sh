#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
PROJECT_DIR=$(cd -- "${SCRIPT_DIR}/.." && pwd)
CONFIG_FILE="${PROJECT_DIR}/config/local.config"
DEFAULT_DB_DIR="${PROJECT_DIR}/resources/database"
DEFAULT_DB_PROFILE="curated"
DATABASE_MANAGER="${PROJECT_DIR}/scripts/database_manager.py"
SMOKE_FASTA="${PROJECT_DIR}/data/example/LKH565_P11_Ci.fna"
PYTHON=${PYTHON:-python3}

main() {
    local command="${1:-run}"
    shift || true

    case "${command}" in
        setup) setup_database "$@" ;;
        run) run_pipeline "$@" ;;
        smoke) run_smoke "$@" ;;
        *)
            echo "Usage: scripts/pipeline_cli.sh {setup|run|smoke} [arguments]" >&2
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

run_nextflow() {
    if [[ -n "${TERM:-}" ]] && ! tput colors >/dev/null 2>&1; then
        TERM=xterm-256color command nextflow "$@"
        return
    fi
    command nextflow "$@"
}

cli_has_option() {
    local option="$1"
    shift
    local arg
    for arg in "$@"; do
        [[ "${arg}" == "${option}" || "${arg}" == "${option}="* ]] && return 0
    done
    return 1
}

cli_has_database_path() {
    cli_has_option --database_path "$@"
}

cli_has_database_profile() {
    cli_has_option --database_profile "$@"
}

read_cli_option() {
    local option="$1"
    shift
    local arg
    local expect_value=0
    for arg in "$@"; do
        if [[ ${expect_value} -eq 1 ]]; then
            printf '%s\n' "${arg}"
            return 0
        fi
        case "${arg}" in
            "${option}="*)
                printf '%s\n' "${arg#*=}"
                return 0
                ;;
            "${option}") expect_value=1 ;;
        esac
    done
    return 1
}

read_cli_database_profile() {
    read_cli_option --database_profile "$@"
}

read_cli_database_path() {
    read_cli_option --database_path "$@"
}

read_config_database_path() {
    [[ -f "${CONFIG_FILE}" ]] || return 1
    if grep -Eq '^database_path=' "${CONFIG_FILE}"; then
        sed -n 's/^database_path=//p' "${CONFIG_FILE}" | sed -n '1p'
    elif grep -Eq '^[[:space:]]*database_path[[:space:]]*=' "${CONFIG_FILE}"; then
        awk -F"'" '/^[[:space:]]*database_path[[:space:]]*=/ {print $2; exit}' "${CONFIG_FILE}"
    else
        sed -n '1p' "${CONFIG_FILE}"
    fi
}

read_config_database_profile() {
    [[ -f "${CONFIG_FILE}" ]] || return 1
    if grep -Eq '^database_profile=' "${CONFIG_FILE}"; then
        sed -n 's/^database_profile=//p' "${CONFIG_FILE}" | sed -n '1p'
    elif grep -Eq '^[[:space:]]*database_profile[[:space:]]*=' "${CONFIG_FILE}"; then
        awk -F"'" '/^[[:space:]]*database_profile[[:space:]]*=/ {print $2; exit}' "${CONFIG_FILE}"
    else
        return 1
    fi
}

write_database_config() {
    local db_dir="$1"
    local profile="$2"
    local temp_file="${CONFIG_FILE}.tmp.$$"
    if [[ -z "${db_dir}" || "${db_dir}" == *$'\n'* || "${db_dir}" == *$'\r'* ]]; then
        printf 'Database path must be non-empty, single-line text.\n' >&2
        return 1
    fi
    if [[ ! "${profile}" =~ ^[A-Za-z0-9][A-Za-z0-9._-]*$ ]]; then
        printf 'Database profile is not a safe identifier: %s\n' "${profile}" >&2
        return 1
    fi
    printf 'database_path=%s\ndatabase_profile=%s\n' "${db_dir}" "${profile}" > "${temp_file}"
    mv "${temp_file}" "${CONFIG_FILE}"
}

profile_from_selection() {
    local selection="$1"
    local default_profile="$2"
    local profiles="$3"
    local selected=""

    if [[ -z "${selection}" ]]; then
        selected="${default_profile}"
    elif [[ "${selection}" =~ ^[0-9]+$ ]]; then
        selected=$(printf '%s\n' "${profiles}" | sed -n "${selection}p" | cut -f1)
    else
        selected="${selection}"
    fi
    if printf '%s\n' "${profiles}" | cut -f1 | grep -Fxq -- "${selected}"; then
        printf '%s\n' "${selected}"
        return 0
    fi
    return 1
}

prompt_database_profile() {
    local default_profile="$1"
    local profiles=""
    local default_number=""
    local profile_count=""
    local selection=""
    local selected=""

    if ! profiles=$("${PYTHON}" "${DATABASE_MANAGER}" profiles); then
        printf 'Could not list database profiles.\n' >&2
        return 1
    fi
    if [[ -z "${profiles}" ]]; then
        printf 'No database profiles are configured.\n' >&2
        return 1
    fi
    default_number=$(
        printf '%s\n' "${profiles}" |
            awk -F '\t' -v profile="${default_profile}" '$1 == profile {print NR; exit}'
    )
    profile_count=$(printf '%s\n' "${profiles}" | awk 'END {print NR}')
    [[ -n "${default_number}" ]] || default_number=1
    printf 'Available database profiles:\n' >&2
    printf '%s\n' "${profiles}" |
        awk -F '\t' '{printf "  %d) %s v%s (%s) - %s\n", NR, $1, $2, $3, $4}' >&2
    while true; do
        printf 'Database profile (1-%s) [default: %s (%s)]: ' \
            "${profile_count}" \
            "${default_number}" \
            "${default_profile}" >&2
        read -r selection
        if selected=$(profile_from_selection "${selection}" "${default_profile}" "${profiles}"); then
            printf '%s\n' "${selected}"
            return 0
        fi
        printf 'Choose a listed number or profile name.\n' >&2
    done
}

prompt_database_path() {
    local default_path="$1"
    local user_path=""
    printf 'Database path [%s]: ' "${default_path}" >&2
    read -r user_path
    printf '%s\n' "${user_path:-${default_path}}"
}

resolve_database_profile() {
    local configured=""
    if read_cli_database_profile "$@"; then
        return
    fi
    if configured=$(read_config_database_profile 2>/dev/null); then
        printf '%s\n' "${configured}"
        return
    fi
    if [[ -t 0 ]]; then
        prompt_database_profile "${DEFAULT_DB_PROFILE}"
        return
    fi
    printf 'Database profile not configured, using %s.\n' "${DEFAULT_DB_PROFILE}" >&2
    printf '%s\n' "${DEFAULT_DB_PROFILE}"
}

resolve_setup_database_profile() {
    local configured="${DEFAULT_DB_PROFILE}"
    if read_cli_database_profile "$@"; then
        return
    fi
    configured=$(read_config_database_profile 2>/dev/null || printf '%s\n' "${configured}")
    if [[ -t 0 ]]; then
        prompt_database_profile "${configured}"
    else
        printf 'Database profile not specified, using %s.\n' "${configured}" >&2
        printf '%s\n' "${configured}"
    fi
}

resolve_database_path() {
    local configured=""
    if read_cli_database_path "$@"; then
        return
    fi
    if configured=$(read_config_database_path 2>/dev/null); then
        printf '%s\n' "${configured}"
        return
    fi
    if [[ -t 0 ]]; then
        prompt_database_path "${DEFAULT_DB_DIR}"
        return
    fi
    printf 'Database path not configured, using %s.\n' "${DEFAULT_DB_DIR}" >&2
    printf '%s\n' "${DEFAULT_DB_DIR}"
}

database_ready() {
    local db_dir="$1"
    local profile="$2"
    if "${PYTHON}" "${DATABASE_MANAGER}" validate \
        --root "${db_dir}" \
        --profile "${profile}" \
        >/dev/null 2>&1; then
        return 0
    fi

    # The resolver is needed only for the deprecated unprofiled curated layout.
    "${PYTHON}" "${DATABASE_MANAGER}" resolve \
        --root "${db_dir}" \
        --profile "${profile}" \
        --model RF00177 \
        >/dev/null 2>&1
}

database_version() {
    local db_dir="$1"
    local profile="$2"
    "${PYTHON}" "${DATABASE_MANAGER}" version \
        --root "${db_dir}" \
        --profile "${profile}" 2>/dev/null
}

install_database() {
    local db_dir="$1"
    local profile="$2"
    local install_args=(
        "${PYTHON}" "${DATABASE_MANAGER}" install
        --root "${db_dir}"
        --profile "${profile}"
        --latest
    )
    if [[ -e "${db_dir}/${profile}" ]]; then
        printf 'Existing profile %s failed validation and will be replaced safely.\n' \
            "${db_dir}/${profile}" >&2
        install_args+=(--force)
    fi
    "${install_args[@]}"
}

ensure_database() {
    local db_dir="$1"
    local profile="$2"
    local version=""

    if ! database_ready "${db_dir}" "${profile}"; then
        printf 'Database profile %s is not ready in %s; starting installation.\n' \
            "${profile}" "${db_dir}" >&2
        install_database "${db_dir}" "${profile}"
        if ! database_ready "${db_dir}" "${profile}"; then
            printf 'Installed database profile failed validation: %s/%s\n' \
                "${db_dir}" "${profile}" >&2
            return 1
        fi
    fi

    if version=$(database_version "${db_dir}" "${profile}"); then
        printf 'Database profile %s v%s ready at %s\n' \
            "${profile}" "${version}" "${db_dir}" >&2
    else
        printf 'Legacy database ready at %s; reinstall the curated profile to track versions.\n' \
            "${db_dir}" >&2
    fi
}

check_database_update() {
    local db_dir="$1"
    local profile="$2"
    if ! database_version "${db_dir}" "${profile}" >/dev/null; then
        return 0
    fi
    printf 'Checking Zenodo for database updates...\n' >&2
    "${PYTHON}" "${DATABASE_MANAGER}" check-update \
        --root "${db_dir}" \
        --profile "${profile}" >&2
}

setup_requests_update() {
    local arg
    for arg in "$@"; do
        [[ "${arg}" == "--update" ]] && return 0
    done
    return 1
}

setup_database() {
    local db_dir=""
    local profile=""
    local result=""
    local status=""
    local installed=""
    local latest=""
    local reason=""
    local answer=""

    profile=$(resolve_setup_database_profile "$@")
    db_dir=$(resolve_database_path "$@")
    write_database_config "${db_dir}" "${profile}"

    if ! database_ready "${db_dir}" "${profile}"; then
        ensure_database "${db_dir}" "${profile}"
        return
    fi
    ensure_database "${db_dir}" "${profile}"
    if ! database_version "${db_dir}" "${profile}" >/dev/null; then
        return
    fi

    printf 'Checking Zenodo for database updates...\n' >&2
    result=$("${PYTHON}" "${DATABASE_MANAGER}" check-update \
        --root "${db_dir}" \
        --profile "${profile}" \
        --format tsv)
    IFS=$'\t' read -r status installed latest reason <<< "${result}"
    case "${status}" in
        update_available)
            printf 'Database update available for %s: v%s -> v%s.\n' \
                "${profile}" "${installed}" "${latest}" >&2
            if setup_requests_update "$@"; then
                answer=y
            elif [[ -t 0 ]]; then
                printf 'Install the update now? [y/N]: ' >&2
                read -r answer
            fi
            if [[ "${answer}" =~ ^[Yy]$ ]]; then
                "${PYTHON}" "${DATABASE_MANAGER}" install \
                    --root "${db_dir}" \
                    --profile "${profile}" \
                    --latest \
                    --force
            else
                printf "Run 'pixi run setup --database_profile %s --update' to update.\n" \
                    "${profile}" >&2
            fi
            ;;
        current)
            printf 'Database profile %s is current at v%s.\n' \
                "${profile}" "${installed}" >&2
            ;;
        installed_newer)
            printf 'Installed database v%s is newer than Zenodo v%s.\n' \
                "${installed}" "${latest}" >&2
            ;;
        unavailable)
            printf 'Database update check unavailable: %s\n' "${reason}" >&2
            ;;
        *)
            printf 'Unexpected database update status: %s\n' "${status}" >&2
            return 1
            ;;
    esac
}

run_pipeline() {
    local db_dir=""
    local profile=""
    local nextflow_args=()

    if [[ "$#" -eq 1 && "$1" == "--version" ]]; then
        "${PYTHON}" "${PROJECT_DIR}/scripts/check_version.py"
        return
    fi
    if has_information_flag "$@"; then
        run_nextflow run "${PROJECT_DIR}/main.nf" "$@"
        return
    fi

    profile=$(resolve_database_profile "$@")
    db_dir=$(resolve_database_path "$@")
    if ! cli_has_database_path "$@" && ! cli_has_database_profile "$@"; then
        write_database_config "${db_dir}" "${profile}"
    fi
    ensure_database "${db_dir}" "${profile}"
    check_database_update "${db_dir}" "${profile}"

    if cli_has_database_path "$@" && cli_has_database_profile "$@"; then
        run_nextflow run "${PROJECT_DIR}/main.nf" "$@"
        return
    fi

    nextflow_args=(run_nextflow run "${PROJECT_DIR}/main.nf")
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
