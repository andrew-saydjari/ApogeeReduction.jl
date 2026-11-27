#!/usr/bin/env bash
# Render the systemd user units in airflow/systemd/*.in (services + the
# watchdog timer) into OUTDIR.
#
# It only WRITES INTO OUTDIR: it never installs, reloads or restarts
# anything. Render into a scratch directory, diff against what is installed,
# then copy by hand (airflow/INSTALL.md step 4).
#
# Usage:
#   render_systemd_units.sh --sshpass PATH --pass-file PATH [options] OUTDIR
#
# Options (defaults are the production deployment on ccalin051):
#   --airflow-home DIR    AIRFLOW_HOME, the unit's WorkingDirectory
#                         [/mnt/home/sdssv/airflow]
#   --env-file FILE       EnvironmentFile (secrets; never committed)
#                         [<airflow-home>/airflow_env.sh]
#   --airflow-bin FILE    airflow executable
#                         [/mnt/home/sdssv/uv_env/airflow_env/bin/airflow]
#   --sshpass FILE        sshpass executable (REQUIRED, account-specific)
#   --pass-file FILE      file holding the Utah password (REQUIRED; the path
#                         only — its content never goes anywhere near git)
#   --tunnel-host ALIAS   ssh alias of the Utah host [mwm]
#   --db-forward SPEC     ssh -L spec for the almanac DB
#                         [63333:operations.sdss.org:5432]
#   --repo DIR            ApogeeReduction.jl clone the watchdog script runs
#                         from [/mnt/home/sdssv/gitcode/ApogeeReduction.jl]
#   --python FILE         python for the watchdog (stdlib only)
#                         [<dir of --airflow-bin>/python]
set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
src="$here/../systemd"

airflow_home=/mnt/home/sdssv/airflow
env_file=""
airflow_bin=/mnt/home/sdssv/uv_env/airflow_env/bin/airflow
sshpass=""
pass_file=""
tunnel_host=mwm
db_forward=63333:operations.sdss.org:5432
repo=/mnt/home/sdssv/gitcode/ApogeeReduction.jl
python=""
outdir=""

usage() { awk 'NR>1 && /^#/ {sub(/^# ?/, ""); print; next} NR>1 {exit}' \
            "${BASH_SOURCE[0]}"; exit "${1:-2}"; }

while [ $# -gt 0 ]; do
  case "$1" in
    --airflow-home) airflow_home="$2"; shift 2;;
    --env-file) env_file="$2"; shift 2;;
    --airflow-bin) airflow_bin="$2"; shift 2;;
    --sshpass) sshpass="$2"; shift 2;;
    --pass-file) pass_file="$2"; shift 2;;
    --tunnel-host) tunnel_host="$2"; shift 2;;
    --db-forward) db_forward="$2"; shift 2;;
    --repo) repo="$2"; shift 2;;
    --python) python="$2"; shift 2;;
    -h|--help) usage 0;;
    -*) echo "unknown option: $1" >&2; usage;;
    *) [ -z "$outdir" ] || usage; outdir="$1"; shift;;
  esac
done

[ -n "$outdir" ] || { echo "OUTDIR is required" >&2; usage; }
[ -n "$sshpass" ] || { echo "--sshpass is required" >&2; usage; }
[ -n "$pass_file" ] || { echo "--pass-file is required" >&2; usage; }
env_file="${env_file:-$airflow_home/airflow_env.sh}"
python="${python:-$(dirname "$airflow_bin")/python}"

for v in "$airflow_home" "$env_file" "$airflow_bin" "$sshpass" "$pass_file" \
         "$repo" "$python"; do
  case "$v" in
    /*) ;;
    *) echo "paths must be absolute (systemd does not expand them): $v" >&2
       exit 2;;
  esac
done
for v in "$airflow_home" "$env_file" "$airflow_bin" "$sshpass" "$pass_file" \
         "$tunnel_host" "$db_forward" "$repo" "$python"; do
  case "$v" in
    *[[:space:]\|\&\\%]*|"")
      echo "refusing a value with whitespace or one of | & \\ %: '$v'" >&2
      exit 2;;
  esac
done
[[ "$db_forward" =~ ^[0-9]+:[A-Za-z0-9.-]+:[0-9]+$ ]] || {
  echo "--db-forward must look like PORT:HOST:PORT, got '$db_forward'" >&2
  exit 2; }
db_port="${db_forward%%:*}"

mkdir -p "$outdir"
for in in "$src"/*.in; do
  out="$outdir/$(basename "${in%.in}")"
  sed -e "s|@AIRFLOW_HOME@|$airflow_home|g" \
      -e "s|@AIRFLOW_ENV_FILE@|$env_file|g" \
      -e "s|@AIRFLOW_BIN@|$airflow_bin|g" \
      -e "s|@SSHPASS@|$sshpass|g" \
      -e "s|@PASS_FILE@|$pass_file|g" \
      -e "s|@TUNNEL_HOST@|$tunnel_host|g" \
      -e "s|@DB_FORWARD@|$db_forward|g" \
      -e "s|@DB_PORT@|$db_port|g" \
      -e "s|@REPO@|$repo|g" \
      -e "s|@PYTHON@|$python|g" \
      "$in" > "$out"
  if grep -n '@[A-Z_]*@' "$out" >&2; then
    echo "unfilled placeholder in $out" >&2
    exit 1
  fi
  echo "rendered $out"
done
