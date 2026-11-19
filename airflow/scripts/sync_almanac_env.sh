#!/bin/bash
# Re-sync the production almanac env with the almanac clone that update.repo
# just pulled, and FAIL LOUDLY if that does not work.
#
# Usage: sync_almanac_env.sh <almanac env python> <almanac clone>
#
# The env holds an editable install of the clone, so code changes need no
# action, but a NEW dependency (e.g. healpy for `almanac add metadata
# --extended`) is only installed when the install is re-run. This re-runs
#     uv pip install -p <python> -e <clone>
# which adds/upgrades only what the clone's pyproject.toml now requires (it
# keeps installed versions that still satisfy it, and does not read uv.lock),
# then checks that the env can import almanac's CLI.
#
# Replaces the inline `uv pip install ... --quiet 2>&1 | tail -1 || echo
# WARNING` of update.repo, whose `|| echo` could never fire: the exit status
# of a pipeline is that of its LAST command (`tail`), so a failed install was
# reported as success and the env silently stayed stale.
#
# Exit status: 0 = env in sync and importable; 1 = install or import failed
# (the task fails, is retried, and notifies; nothing else is changed).

set -u
py="${1:?usage: sync_almanac_env.sh <almanac env python> <almanac clone>}"
dir="${2:?usage: sync_almanac_env.sh <almanac env python> <almanac clone>}"
uv_bin="${UV:-uv}"

if ! command -v "$uv_bin" > /dev/null 2>&1; then
  echo "FATAL: '$uv_bin' not found on PATH; cannot re-sync the almanac env"
  exit 1
fi

echo "=== almanac env re-sync: $uv_bin pip install -p $py -e $dir"
out=$("$uv_bin" pip install -p "$py" -e "$dir" 2>&1)
rc=$?
# uv lists what it changed ("+ healpy==1.20.1", "~ sdss-almanac ...").
printf '%s\n' "$out" | tail -n 25
if [ $rc -ne 0 ]; then
  echo "FATAL: almanac env re-sync failed (uv exit $rc). The clone was pulled"
  echo "but its dependencies may be missing from $py; fix and re-run update.repo."
  exit 1
fi

if ! imp=$("$py" -c "import almanac.cli, almanac; print(almanac.__version__, almanac.__file__)" 2>&1); then
  printf '%s\n' "$imp" | tail -n 25
  echo "FATAL: the almanac env cannot import almanac.cli after the re-sync"
  exit 1
fi
echo "almanac env OK: $imp"
