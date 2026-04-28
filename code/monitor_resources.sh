#!/bin/bash
# monitor_resources.sh
# Monitors RAM and CPU usage of an R process at a fixed interval.
# Writes to a CSV log that survives if the monitored process crashes.
#
# Usage:
#   bash monitor_resources.sh <pid> <log_file> [interval_seconds]
#
# Examples:
#   bash monitor_resources.sh 12345 ~/resource_log.csv 5
#   bash monitor_resources.sh 12345 ~/resource_log.csv       # defaults to 5s interval
#
# To start your R script and monitor it in one go:
#   Rscript my_analysis.R &
#   R_PID=$!
#   bash monitor_resources.sh $R_PID ~/resource_log.csv 5
#
# To read the log in R afterwards:
#   log <- readr::read_csv("~/resource_log.csv")
#   ggplot(log, aes(x = as.POSIXct(timestamp), y = rss_mb)) + geom_line()

# ── Arguments ─────────────────────────────────────────────────────────────────

PID=$1
LOG=$2
INTERVAL=${3:-5}

if [ -z "$PID" ] || [ -z "$LOG" ]; then
  echo "Usage: bash monitor_resources.sh <pid> <log_file> [interval_seconds]"
  exit 1
fi

if ! kill -0 "$PID" 2>/dev/null; then
  echo "Error: no process found with PID $PID"
  exit 1
fi

# ── System RAM helper ─────────────────────────────────────────────────────────
# Tries 'free' first (Linux); falls back to vm_stat (macOS)

get_system_ram_mb() {
  if command -v free &>/dev/null; then
    free -m | awk '/^Mem:/ { printf "%s,%s,%s", $2, $3, $4 }'
  elif command -v vm_stat &>/dev/null; then
    # vm_stat reports pages; page size is typically 4096 bytes
    PAGE=$(vm_stat | awk '/page size/ { print $8 }')
    PAGE=${PAGE:-4096}
    vm_stat | awk -v page="$PAGE" '
      /Pages active/    { active = $3 }
      /Pages wired/     { wired  = $4 }
      /Pages free/      { free   = $3 }
      END {
        used_mb  = int((active + wired) * page / 1048576)
        free_mb  = int(free * page / 1048576)
        total_mb = used_mb + free_mb
        printf "%s,%s,%s", total_mb, used_mb, free_mb
      }
    '
  else
    echo "NA,NA,NA"
  fi
}

# ── Header ────────────────────────────────────────────────────────────────────

echo "timestamp,pid,rss_mb,vsz_mb,cpu_pct,mem_pct,sys_total_mb,sys_used_mb,sys_free_mb" > "$LOG"

echo "monitoring PID $PID every ${INTERVAL}s — writing to $LOG"
echo "press Ctrl+C to stop monitoring manually"

# ── Main loop ─────────────────────────────────────────────────────────────────

while kill -0 "$PID" 2>/dev/null; do
  TS=$(date '+%Y-%m-%d %H:%M:%S')

  # ps columns: rss (kb), vsz (kb), %cpu, %mem
  STATS=$(ps -p "$PID" -o rss=,vsz=,pcpu=,pmem= 2>/dev/null)

  if [ -n "$STATS" ]; then
    RSS=$(echo $STATS | awk '{ printf "%.1f", $1/1024 }')
    VSZ=$(echo $STATS | awk '{ printf "%.1f", $2/1024 }')
    CPU=$(echo $STATS | awk '{ print $3 }')
    MEM=$(echo $STATS | awk '{ print $4 }')
  else
    # process exists but ps returned nothing — log NAs rather than skipping
    RSS="NA"; VSZ="NA"; CPU="NA"; MEM="NA"
  fi

  SYS=$(get_system_ram_mb)

  echo "$TS,$PID,$RSS,$VSZ,$CPU,$MEM,$SYS" >> "$LOG"

  sleep "$INTERVAL"
done

# ── Process ended ─────────────────────────────────────────────────────────────

END_TS=$(date '+%Y-%m-%d %H:%M:%S')
echo "$END_TS,$PID,NA,NA,NA,NA,NA,NA,NA" >> "$LOG"
echo ""
echo "[$END_TS] PID $PID is no longer running. Log written to $LOG"