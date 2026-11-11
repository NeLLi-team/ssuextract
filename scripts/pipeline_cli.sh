#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
PROJECT_DIR=$(cd -- "${SCRIPT_DIR}/.." && pwd)
CONFIG_FILE="${PROJECT_DIR}/config/local.config"
DEFAULT_DB_DIR="${PROJECT_DIR}/resources/database"
DEFAULT_DB_PROFILE="curated"
DATABASE_MANAGER="${PROJECT_DIR}/scripts/database_manager.py"
SMOKE_FASTAS=(
    "${PROJECT_DIR}/data/example/LKH462_P08_Rh.fna"
    "${PROJECT_DIR}/data/example/LKH565_P11_Ci.fna"
)
EXAMPLE_EXPECTATIONS="${PROJECT_DIR}/data/example/expected_annotations.tsv"
EXAMPLE_VALIDATOR="${PROJECT_DIR}/scripts/validate_example_output.py"
EXAMPLE_THREADS_PER_JOB=2
PYTHON=${PYTHON:-python3}
# Exit statuses of `database_manager.py validate` for an absent or damaged profile.
DB_NOT_INSTALLED=3
DB_DAMAGED=4
# Set to 1 by ensure_database when it installed or replaced the profile.
DATABASE_INSTALLED=0
# The selection used by the last successful run, also used by the example validator.
RUN_DATABASE_PATH=""
RUN_DATABASE_PROFILE=""

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

validate_database_options() {
    while [[ "$#" -gt 0 ]]; do
        case "$1" in
            --database_profile | --database_path)
                if [[ "$#" -lt 2 || -z "$2" || "$2" == -* ]]; then
                    printf '%s requires a value that does not start with "-".\n' "$1" >&2
                    return 2
                fi
                shift 2
                ;;
            --database_profile=* | --database_path=*)
                if [[ -z "${1#*=}" || "${1#*=}" == -* ]]; then
                    printf '%s requires a value that does not start with "-".\n' \
                        "${1%%=*}" >&2
                    return 2
                fi
                shift
                ;;
            *) shift ;;
        esac
    done
}

require_known_profile() {
    local profile="$1"
    local names=""
    local profiles=""
    if ! profiles=$("${PYTHON}" "${DATABASE_MANAGER}" profiles); then
        printf 'Could not list database profiles.\n' >&2
        return 1
    fi
    names=$(printf '%s\n' "${profiles}" | cut -f1)
    if printf '%s\n' "${names}" | grep -Fxq -- "${profile}"; then
        return 0
    fi
    printf 'Unknown database profile: %s\nValid profiles: %s\n' "${profile}" \
        "$(printf '%s\n' "${names}" | awk 'NR > 1 {printf ", "} {printf "%s", $0}')" >&2
    return 2
}

read_config_value() {
    local key="$1"
    local value=""
    [[ -f "${CONFIG_FILE}" ]] || return 1
    value=$(sed -n "s/^${key}=//p" "${CONFIG_FILE}" | sed -n '1p')
    [[ -n "${value}" ]] || return 1
    printf '%s\n' "${value}"
}

read_config_database_path() {
    read_config_value database_path
}

read_config_database_profile() {
    read_config_value database_profile
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
    if [[ -f "${CONFIG_FILE}" && \
        "$(<"${CONFIG_FILE}")" == "database_path=${db_dir}"$'\n'"database_profile=${profile}" ]]; then
        return 0
    fi
    # Callers use `|| return`, which turns errexit off here: check each step.
    if ! printf 'database_path=%s\ndatabase_profile=%s\n' "${db_dir}" "${profile}" \
        > "${temp_file}" || ! mv "${temp_file}" "${CONFIG_FILE}"; then
        rm -f "${temp_file}"
        printf 'Could not save the database selection to %s.\n' "${CONFIG_FILE}" >&2
        return 1
    fi
}

profile_from_selection() {
    local selection="$1"
    local default_profile="$2"
    local profiles="$3"
    local selected=""

    if [[ -z "${selection}" ]]; then
        selected="${default_profile}"
    elif [[ "${selection}" =~ ^0*[1-9][0-9]*$ ]]; then
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
    user_path="${user_path:-${default_path}}"
    case "${user_path}" in
        "~") user_path="${HOME}" ;;
        "~/"*) user_path="${HOME}/${user_path#\~/}" ;;
    esac
    [[ "${user_path}" == /* ]] || user_path="${PWD}/${user_path}"
    printf '%s\n' "${user_path}"
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
    shift 2
    "${PYTHON}" "${DATABASE_MANAGER}" validate \
        --root "${db_dir}" \
        --profile "${profile}" \
        "$@" \
        >/dev/null
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
    local replace="${3:-0}"
    local install_args=(
        "${PYTHON}" "${DATABASE_MANAGER}" install
        --root "${db_dir}"
        --profile "${profile}"
        --latest
    )
    if [[ "${replace}" -eq 1 ]]; then
        install_args+=(--force)
    fi
    ${install_args[@]+"${install_args[@]}"}
}

report_database_ready() {
    local db_dir="$1"
    local profile="$2"
    local version=""
    version=$(database_version "${db_dir}" "${profile}") || return 1
    printf 'Database profile %s v%s ready at %s\n' \
        "${profile}" "${version}" "${db_dir}" >&2
}

# Installs an absent profile. Replaces a damaged profile unless the user gave
# --database_path. Any other validation failure is an error.
ensure_database() {
    local db_dir="$1"
    local profile="$2"
    local explicit_path="${3:-0}"
    local validation_args=("${@:4}")
    local status=0

    DATABASE_INSTALLED=0
    database_ready "${db_dir}" "${profile}" ${validation_args[@]+"${validation_args[@]}"} || status=$?
    case "${status}" in
        0) ;;
        "${DB_NOT_INSTALLED}")
            printf 'Database profile %s is not installed in %s; starting installation.\n' \
                "${profile}" "${db_dir}" >&2
            install_database "${db_dir}" "${profile}" 0 || return 1
            DATABASE_INSTALLED=1
            ;;
        "${DB_DAMAGED}")
            if [[ "${explicit_path}" -eq 1 ]]; then
                printf 'Database profile %s in %s is damaged. It is not replaced because --database_path was given.\n' \
                    "${profile}" "${db_dir}" >&2
                printf 'To reinstall it, remove %s/%s and run setup again.\n' \
                    "${db_dir}" "${profile}" >&2
                return 1
            fi
            printf 'Database profile %s in %s is damaged and will be replaced safely.\n' \
                "${profile}" "${db_dir}" >&2
            install_database "${db_dir}" "${profile}" 1 || return 1
            DATABASE_INSTALLED=1
            ;;
        *) return 1 ;;
    esac

    if [[ "${DATABASE_INSTALLED}" -eq 1 ]] && \
        ! database_ready "${db_dir}" "${profile}" ${validation_args[@]+"${validation_args[@]}"}; then
        printf 'Installed database profile failed validation: %s/%s\n' \
            "${db_dir}" "${profile}" >&2
        return 1
    fi
    report_database_ready "${db_dir}" "${profile}"
}

check_database_update() {
    local db_dir="$1"
    local profile="$2"
    local allow_prompt="${3:-0}"
    local force_update="${4:-0}"
    local include_database_path="${5:-0}"
    local update_args=("${@:6}")
    local result=""
    local status=""
    local installed=""
    local latest=""
    local reason=""
    local answer=""

    result=$("${PYTHON}" "${DATABASE_MANAGER}" check-update \
        --root "${db_dir}" \
        --profile "${profile}" \
        --format tsv ${update_args[@]+"${update_args[@]}"}) || return 1
    IFS=$'\t' read -r status installed latest reason <<< "${result}"
    case "${status}" in
        update_available)
            printf 'Database update available for %s: v%s -> v%s.\n' \
                "${profile}" "${installed}" "${latest}" >&2
            if [[ "${force_update}" -eq 1 ]]; then
                answer=y
            elif can_prompt_for_database_update "${allow_prompt}"; then
                printf 'Update database now? [y/N]: ' >&2
                if ! read -r answer; then
                    answer=""
                fi
            fi
            if [[ "${answer}" =~ ^[Yy]$ ]]; then
                if ! install_database_update "${db_dir}" "${profile}" "${latest}"; then
                    return 1
                fi
            else
                print_database_update_command \
                    "${db_dir}" "${profile}" "${include_database_path}"
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
            if [[ "${force_update}" -eq 1 ]]; then
                return 1
            fi
            ;;
        *)
            printf 'Unexpected database update status: %s\n' "${status}" >&2
            return 1
            ;;
    esac
}

can_prompt_for_database_update() {
    local allow_prompt="$1"
    [[ "${allow_prompt}" -eq 1 && -t 0 ]]
}

install_database_update() {
    local db_dir="$1"
    local profile="$2"
    local version="$3"
    if ! "${PYTHON}" "${DATABASE_MANAGER}" install \
        --root "${db_dir}" \
        --profile "${profile}" \
        --latest \
        --require-version "${version}" \
        --force; then
        printf 'Database update failed; profile %s remains installed and the pipeline was not started.\n' \
            "${profile}" >&2
        return 1
    fi
    database_ready "${db_dir}" "${profile}" || return 1
    report_database_ready "${db_dir}" "${profile}"
}

print_database_update_command() {
    local db_dir="$1"
    local profile="$2"
    local include_database_path="$3"
    local quoted_db_dir=""

    if [[ "${include_database_path}" -eq 1 ]]; then
        printf -v quoted_db_dir '%q' "${db_dir}"
        printf 'Run: pixi run setup --database_path %s --database_profile %s --update\n' \
            "${quoted_db_dir}" "${profile}" >&2
    else
        printf 'Run: pixi run setup --database_profile %s --update\n' \
            "${profile}" >&2
    fi
}

setup_requests_update() {
    local arg
    for arg in "$@"; do
        [[ "${arg}" == "--update" ]] && return 0
    done
    return 1
}

print_setup_usage() {
    cat <<'USAGE'
Usage: pixi run setup [--database_profile NAME] [--database_path DIR] [--update]

Install a database profile, or with --update check for and install a newer
release. In a terminal, setup asks for a profile or path that you do not give.
USAGE
}

validate_setup_arguments() {
    validate_database_options "$@" || return
    while [[ "$#" -gt 0 ]]; do
        case "$1" in
            --update | --database_profile=* | --database_path=*) shift ;;
            --database_profile | --database_path) shift 2 ;;
            *)
                printf 'Unknown setup argument: %s\n' "$1" >&2
                print_setup_usage >&2
                return 2
                ;;
        esac
    done
}

setup_has_help() {
    local arg
    for arg in "$@"; do
        [[ "${arg}" == "-h" || "${arg}" == "--help" ]] && return 0
    done
    return 1
}

setup_database() {
    local db_dir=""
    local profile=""
    local force_update=0
    local explicit_path=0

    if setup_has_help "$@"; then
        print_setup_usage
        return 0
    fi
    validate_setup_arguments "$@" || return
    profile=$(resolve_setup_database_profile "$@") || return
    require_known_profile "${profile}" || return
    db_dir=$(resolve_database_path "$@") || return
    if cli_has_database_path "$@"; then
        explicit_path=1
    fi
    ensure_database "${db_dir}" "${profile}" "${explicit_path}" || return
    write_database_config "${db_dir}" "${profile}" || return
    if setup_requests_update "$@"; then
        force_update=1
    fi
    check_database_update "${db_dir}" "${profile}" 1 "${force_update}"
}

run_pipeline() {
    local db_dir=""
    local profile=""
    local allow_update_prompt=1
    local include_database_path=0
    local nextflow_args=()
    local normalized_args=()

    while [[ "$#" -gt 0 ]]; do
        case "$1" in
            -q)
                if [[ "$#" -lt 2 ]]; then
                    printf '%s\n' '-q requires a FASTA file or directory path.' >&2
                    return 2
                fi
                if [[ -z "$2" || "$2" == -* ]]; then
                    printf '%s\n' '-q requires a FASTA file or directory path.' >&2
                    return 2
                fi
                normalized_args+=(--query "$2")
                shift 2
                ;;
            -q=*)
                if [[ -z "${1#*=}" ]]; then
                    printf '%s\n' '-q requires a FASTA file or directory path.' >&2
                    return 2
                fi
                normalized_args+=(--query "${1#*=}")
                shift
                ;;
            *)
                normalized_args+=("$1")
                shift
                ;;
        esac
    done
    set -- ${normalized_args[@]+"${normalized_args[@]}"}

    if [[ "$#" -eq 1 && "$1" == "--version" ]]; then
        "${PYTHON}" "${PROJECT_DIR}/scripts/check_version.py"
        return
    fi
    if has_information_flag "$@"; then
        run_nextflow run "${PROJECT_DIR}/main.nf" "$@"
        return
    fi

    validate_database_options "$@" || return
    profile=$(resolve_database_profile "$@") || return
    require_known_profile "${profile}" || return
    db_dir=$(resolve_database_path "$@") || return
    if cli_has_database_path "$@"; then
        allow_update_prompt=0
        include_database_path=1
    fi
    ensure_database "${db_dir}" "${profile}" "${include_database_path}" --quick || return
    if ! cli_has_database_path "$@" && ! cli_has_database_profile "$@"; then
        write_database_config "${db_dir}" "${profile}" || return
    fi
    if [[ "${SSUEXTRACT_SKIP_UPDATE_CHECK:-0}" != 1 ]]; then
        if ! check_database_update \
            "${db_dir}" "${profile}" "${allow_update_prompt}" 0 "${include_database_path}" --cached; then
            return 1
        fi
    fi
    RUN_DATABASE_PATH="${db_dir}"
    RUN_DATABASE_PROFILE="${profile}"

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
    ${nextflow_args[@]+"${nextflow_args[@]}"} "$@"
}

run_smoke() {
    local database_version=""
    local db_dir=""
    local profile=""
    local smoke_dir=""
    local smoke_args=()
    if has_information_flag "$@"; then
        run_pipeline "$@"
        return
    fi
    validate_database_options "$@" || return
    smoke_dir=$(mktemp -d "${TMPDIR:-/tmp}/ssuextract-smoke.XXXXXX") || return
    trap "rm -rf '${smoke_dir}'" EXIT
    cp ${SMOKE_FASTAS[@]+"${SMOKE_FASTAS[@]}"} "${smoke_dir}/" || return
    smoke_args=(
        --query "${smoke_dir}"
        --outdir "${PROJECT_DIR}/results/smoke"
        --threads_per_job "${EXAMPLE_THREADS_PER_JOB}"
    )
    # Reuse the selection resolved by run_pipeline, including unsaved overrides.
    run_pipeline ${smoke_args[@]+"${smoke_args[@]}"} "$@" || return
    profile="${RUN_DATABASE_PROFILE}"
    db_dir="${RUN_DATABASE_PATH}"
    database_version=$("${PYTHON}" "${DATABASE_MANAGER}" version \
        --root "${db_dir}" --profile "${profile}") || return
    "${PYTHON}" "${EXAMPLE_VALIDATOR}" \
        --summary "${PROJECT_DIR}/results/smoke/cmsearch_summary.tsv" \
        --expectations "${EXAMPLE_EXPECTATIONS}" \
        --profile "${profile}" \
        --database-version "${database_version}"
}

if [[ "${BASH_SOURCE[0]}" == "$0" ]]; then
    main "$@"
fi
