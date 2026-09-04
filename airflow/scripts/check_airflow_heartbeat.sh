#!/bin/bash
# Cross-machine watchdog for the APOGEE Airflow instance (card O3/N3).
#
# Run this from a DIFFERENT machine than the one hosting Airflow (e.g. from
# rusty scrontab while Airflow lives on ccalin051), so a dead node is caught
# too. It checks:
#   1. the heartbeat file the apogee_heartbeat DAG touches every 10 min
#   2. (SLA substitute; sla_miss_callback was removed in Airflow 3.0) that
#      the daily metrics table gained a row within the last SLA_HOURS
# and posts to Slack when something is stale. A sentinel file rate-limits
# alerts to one per ALERT_COOLDOWN_HOURS per condition.
#
# Suggested scrontab line (AKS to install; see airflow/INSTALL.md):
#   */30 * * * * /mnt/home/sdssv/gitcode/ApogeeReduction.jl/airflow/scripts/check_airflow_heartbeat.sh
#
# Env (or edit defaults below): SLACK_TOKEN required for posting;
# AR_SLACK_CHANNEL optional (defaults to the dev channel).

set -u

HEARTBEAT_FILE="${AR_HEARTBEAT_FILE:-/mnt/ceph/users/sdssv/work/daily/metrics/airflow_heartbeat.txt}"
METRICS_FILE="${AR_METRICS_FILE:-/mnt/ceph/users/sdssv/work/daily/metrics/daily_metrics.csv}"
CHANNEL="${AR_SLACK_CHANNEL:-C07KQ7BJY5P}"   # dev; prod = C08B7FKMP16
HEARTBEAT_STALE_MIN="${HEARTBEAT_STALE_MIN:-45}"    # 10-min beat + slack
SLA_HOURS="${SLA_HOURS:-30}"                         # daily row expected
ALERT_COOLDOWN_HOURS="${ALERT_COOLDOWN_HOURS:-6}"
STATE_DIR="${AR_WATCHDOG_STATE:-$HOME/.ar_watchdog}"

mkdir -p "$STATE_DIR"

post_slack() {
    local key="$1" text="$2"
    local sentinel="$STATE_DIR/last_alert_$key"
    if [ -f "$sentinel" ]; then
        local age=$(( ($(date +%s) - $(stat -c %Y "$sentinel")) / 3600 ))
        if [ "$age" -lt "$ALERT_COOLDOWN_HOURS" ]; then
            echo "suppressed (cooldown, ${age}h < ${ALERT_COOLDOWN_HOURS}h): $text"
            return 0
        fi
    fi
    if [ -z "${SLACK_TOKEN:-}" ]; then
        echo "NO SLACK_TOKEN — would post: $text"
        return 1
    fi
    curl -sS -X POST https://slack.com/api/chat.postMessage \
        -H "Authorization: Bearer ${SLACK_TOKEN}" \
        -H "Content-Type: application/json; charset=utf-8" \
        -d "{\"channel\": \"${CHANNEL}\", \"text\": \"${text}\"}" \
        > /dev/null && touch "$sentinel"
    echo "alerted: $text"
}

now=$(date +%s)
status=0

# 1. heartbeat freshness
if [ ! -f "$HEARTBEAT_FILE" ]; then
    post_slack heartbeat ":rotating_light: APOGEE Airflow watchdog: heartbeat file missing ($HEARTBEAT_FILE) — is airflow standalone running on ccalin051? (systemctl --user status airflow-standalone)"
    status=1
else
    age_min=$(( (now - $(stat -c %Y "$HEARTBEAT_FILE")) / 60 ))
    if [ "$age_min" -gt "$HEARTBEAT_STALE_MIN" ]; then
        post_slack heartbeat ":rotating_light: APOGEE Airflow watchdog: heartbeat is ${age_min} min old (limit ${HEARTBEAT_STALE_MIN}) — scheduler down or apogee_heartbeat DAG paused."
        status=1
    else
        echo "heartbeat OK (${age_min} min old)"
        rm -f "$STATE_DIR/last_alert_heartbeat"
    fi
fi

# 2. daily-metrics freshness (SLA-miss substitute)
if [ -f "$METRICS_FILE" ]; then
    age_h=$(( (now - $(stat -c %Y "$METRICS_FILE")) / 3600 ))
    if [ "$age_h" -gt "$SLA_HOURS" ]; then
        post_slack sla ":warning: APOGEE Airflow watchdog: no new daily metrics row in ${age_h}h (SLA ${SLA_HOURS}h) — the daily DAGs may be failing, skipping, or paused."
        status=1
    else
        echo "daily metrics OK (last row ${age_h}h old)"
        rm -f "$STATE_DIR/last_alert_sla"
    fi
else
    echo "metrics file not present yet ($METRICS_FILE) — skipping SLA check"
fi

exit $status
