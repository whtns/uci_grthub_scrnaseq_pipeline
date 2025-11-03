#!/usr/bin/env bash
set -euo pipefail

# Sync local web_summaries directory to Nextcloud using rclone
# Usage: sync_web_summaries_to_nextcloud.sh [--remote REMOTE] [--local LOCAL_DIR] [--dry-run]

REMOTE_DEFAULT="grthub_nextcloud"
# remote base path template (will append project parent and current dir)
REMOTE_BASE_PATH="/projects"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

PARENT_DIR_NAME="$(basename "$(dirname "$REPO_ROOT")")"
CURRENT_DIR_NAME="$(basename "$REPO_ROOT")"

# Default local web_summaries locations to try (in order)
try_local_paths=("$REPO_ROOT/results/web_summaries")

REMOTE_PATH="${REMOTE_BASE_PATH}/${PARENT_DIR_NAME}/${CURRENT_DIR_NAME}/results/web_summaries"
REMOTE="${REMOTE_DEFAULT}:${REMOTE_PATH}"

DRY_RUN=0
LOCAL_DIR=""

usage() {
  cat <<EOF
Usage: $0 [--remote REMOTE] [--local LOCAL_DIR] [--dry-run] [--help]

Examples:
  $0                       # sync detected local web_summaries to ${REMOTE}
  $0 --dry-run             # show what would be synced
  $0 --remote other:foo    # sync to a different rclone remote:path
  $0 --local /path/to/dir   # specify a local web_summaries dir
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --remote)
      REMOTE="$2"
      shift 2
      ;;
    --local)
      LOCAL_DIR="$2"
      shift 2
      ;;
    --dry-run)
      DRY_RUN=1
      shift
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "Unknown argument: $1" >&2
      usage
      exit 2
      ;;
  esac
done

if [[ -z "$LOCAL_DIR" ]]; then
  for p in "${try_local_paths[@]}"; do
    if [[ -d "$p" ]]; then
      LOCAL_DIR="$p"
      break
    fi
  done
fi

# Pre-sync: collect any web_summary.html files under output/ into results/web_summaries
# This creates a consolidated directory with files named <sample>_web_summary.html
RESULTS_DIR="$REPO_ROOT/results/web_summaries"
mkdir -p "$RESULTS_DIR"

found=0
while IFS= read -r -d '' file; do
  # Try to extract sample name from path like .../cellranger/{sample}/outs/web_summary.html
  sample=""
  if [[ "$file" =~ ${REPO_ROOT}/cellranger/([^/]+)/outs/web_summary.html ]]; then
    sample="${BASH_REMATCH[1]}"
    else
    # Fallbacks: prefer the grandparent (sample) when parent is 'outs', otherwise use parent
    parent="$(basename "$(dirname "$file")")"
    grandparent="$(basename "$(dirname "$(dirname "$file")")")"
    if [[ "$parent" == "outs" && "$grandparent" != "" && "$grandparent" != "." ]]; then
      sample="$grandparent"
    elif [[ "$parent" != "" && "$parent" != "." ]]; then
      sample="$parent"
    elif [[ "$grandparent" != "" && "$grandparent" != "." ]]; then
      sample="$grandparent"
    else
      sample="$(basename "$file" .html)"
    fi
  fi
  dest="$RESULTS_DIR/${sample}_web_summary.html"
  echo "Collecting: $file -> $dest"
  # Optimize copy: fast checks first (size + mtime), then cmp fallback
  if [[ -f "$dest" ]]; then
    src_size=$(stat -c%s "$file" 2>/dev/null || echo "-1")
    dst_size=$(stat -c%s "$dest" 2>/dev/null || echo "-2")
    src_mtime=$(stat -c%Y "$file" 2>/dev/null || echo "0")
    dst_mtime=$(stat -c%Y "$dest" 2>/dev/null || echo "0")

    if [[ "$src_size" == "$dst_size" && "$src_mtime" == "$dst_mtime" ]]; then
      echo "  Skipping copy: size and mtime match: $dest"
    else
      # If sizes differ, definitely copy. If sizes equal but mtimes differ, do a quick cmp.
      if [[ "$src_size" != "$dst_size" ]]; then
        echo "  Copying (size differs): $dest"
        cp -p "$file" "$dest"
      else
        # same size, different mtime: compare contents
        if cmp -s "$file" "$dest"; then
          echo "  Contents identical but mtime differs; updating mtime on destination: $dest"
          # Preserve content but update mtime to match source so sync tools see correct timestamp
          touch -r "$file" "$dest"
        else
          echo "  Copying (same size but content differs): $dest"
          cp -p "$file" "$dest"
        fi
      fi
    fi
  else
    cp -p "$file" "$dest"
  fi
  found=1
done < <(find "$REPO_ROOT/output" -type f -name 'web_summary.html' -print0 2>/dev/null || true)

if [[ $found -eq 1 && -z "$LOCAL_DIR" ]]; then
  LOCAL_DIR="$RESULTS_DIR"
fi

if [[ -z "$LOCAL_DIR" ]]; then
  echo "Could not find a local web_summaries directory. Tried:" >&2
  for p in "${try_local_paths[@]}"; do echo "  - $p" >&2; done
  echo "You can provide one with --local /path/to/web_summaries" >&2
  exit 3
fi

LOG_DIR="$REPO_ROOT/logs"
mkdir -p "$LOG_DIR"
TIMESTAMP="$(date +%Y%m%d_%H%M%S)"
LOG_FILE="$LOG_DIR/sync_web_summaries_${TIMESTAMP}.log"

echo "Local dir: $LOCAL_DIR"
echo "Remote target: $REMOTE"
echo "Log file: $LOG_FILE"

RCLONE_CMD=(rclone sync "$LOCAL_DIR" "$REMOTE" --progress --create-empty-src-dirs)
if [[ $DRY_RUN -eq 1 ]]; then
  RCLONE_CMD+=(--dry-run)
fi

echo "Running: ${RCLONE_CMD[*]}" | tee "$LOG_FILE"
"${RCLONE_CMD[@]}" 2>&1 | tee -a "$LOG_FILE"

EXIT_CODE=${PIPESTATUS[0]:-0}
if [[ $EXIT_CODE -ne 0 ]]; then
  echo "rclone exited with code $EXIT_CODE" | tee -a "$LOG_FILE"
  exit $EXIT_CODE
fi

echo "Sync completed successfully." | tee -a "$LOG_FILE"
