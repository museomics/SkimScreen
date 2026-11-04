#!/bin/bash
#SBATCH --job-name=kraken2_screen
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --cpus-per-task=8
#SBATCH --mem=250G
#SBATCH --partition=himem
set -euo pipefail
# ── Args ──────────────────────────────────────────────────────────────────────
DATABASE="${1:?Usage: $0 <db_dir> <out_dir> <manifest.csv>}"
BASE_OUT="${2:?Usage: $0 <db_dir> <out_dir> <manifest.csv>}"
MANIFEST="${3:?Usage: $0 <db_dir> <out_dir> <manifest.csv>}"
# ── Environment ───────────────────────────────────────────────────────────────
export PATH="/mnt/apps/users/dparsons/conda/bin:$PATH"
source /mnt/apps/users/$USER/conda/etc/profile.d/conda.sh
conda activate read_qc

# ── Validation: database ──────────────────────────────────────────────────────
[[ -d "$DATABASE" ]] || { echo "ERROR: Database directory not found: $DATABASE"; exit 1; }
[[ -f "$MANIFEST" ]] || { echo "ERROR: Manifest not found: $MANIFEST"; exit 1; }
# Only these three files are needed by `k2 classify`
DB_FILES=("${DATABASE}/opts.k2d" "${DATABASE}/taxo.k2d" "${DATABASE}/hash.k2d")
for f in "${DB_FILES[@]}"; do
    [[ -f "$f" ]] || { echo "ERROR: Database file not found: $f"; exit 1; }
done
THREADS="${SLURM_CPUS_PER_TASK:-8}"
# ── Validation: manifest header (columns found by name) ───────────────────────
# Comma-separated. Required columns: ID, forward (or R1), reverse (or R2).
# Any other columns are ignored, and column order does not matter.
IFS= read -r header < "$MANIFEST" || { echo "ERROR: Manifest is empty: $MANIFEST"; exit 1; }
header="${header%$'\r'}"                 # tolerate Windows line endings
header="${header#$'\xEF\xBB\xBF'}"       # tolerate a UTF-8 BOM
IFS=',' read -r -a cols <<< "$header"
ID_COL=""; R1_COL=""; R2_COL=""
for i in "${!cols[@]}"; do
    name="${cols[$i]//\"/}"
    name="${name//[[:space:]]/}"
    case "$name" in
        ID)         ID_COL=$i ;;
        forward|R1) R1_COL=$i ;;
        reverse|R2) R2_COL=$i ;;
    esac
done
missing_cols=()
[[ -n "$ID_COL" ]] || missing_cols+=("ID")
[[ -n "$R1_COL" ]] || missing_cols+=("forward (or R1)")
[[ -n "$R2_COL" ]] || missing_cols+=("reverse (or R2)")
if (( ${#missing_cols[@]} > 0 )); then
    echo "ERROR: Manifest header is missing column(s): ${missing_cols[*]}"
    echo "  Header found: ${header}"
    exit 1
fi
# ── Validation: every sample, before staging the database ─────────────────────
SAMPLE_IDS=(); SAMPLE_R1=(); SAMPLE_R2=()
failed_samples=()
declare -A seen_ids=()
total=0
line_no=1
echo "Checking manifest: ${MANIFEST}"
{
    read -r _   # skip header
    while IFS= read -r line || [[ -n "$line" ]]; do
        line_no=$((line_no + 1))
        line="${line%$'\r'}"
        [[ -z "${line//[[:space:],]/}" ]] && continue   # blank line
        IFS=',' read -r -a fields <<< "$line"
        ID="${fields[$ID_COL]:-}"; ID="${ID//\"/}"
        R1="${fields[$R1_COL]:-}"; R1="${R1//\"/}"
        R2="${fields[$R2_COL]:-}"; R2="${R2//\"/}"
        total=$((total + 1))
        if [[ -z "$ID" ]]; then
            echo "  ERROR: line ${line_no}: no ID — skipping"
            failed_samples+=("line ${line_no} (no ID)")
            continue
        fi
        if [[ -n "${seen_ids[$ID]:-}" ]]; then
            echo "  ERROR: line ${line_no}: duplicate ID ${ID} — skipping"
            failed_samples+=("${ID} (duplicate ID, line ${line_no})")
            continue
        fi
        seen_ids[$ID]=1
        if [[ ! -f "$R1" ]]; then
            echo "  ERROR: ${ID}: R1 not found: ${R1:-<empty>} — skipping"
            failed_samples+=("${ID} (R1 not found)")
            continue
        fi
        if [[ ! -f "$R2" ]]; then
            echo "  ERROR: ${ID}: R2 not found: ${R2:-<empty>} — skipping"
            failed_samples+=("${ID} (R2 not found)")
            continue
        fi
        SAMPLE_IDS+=("$ID"); SAMPLE_R1+=("$R1"); SAMPLE_R2+=("$R2")
    done
} < "$MANIFEST"
n_run=${#SAMPLE_IDS[@]}
echo "Manifest check: ${n_run}/${total} samples ready to run"
if (( n_run == 0 )); then
    echo "ERROR: No runnable samples — not staging the database"
    exit 1
fi
# ── Stage database to node-local SSD ──────────────────────────────────────────
# The DB is memory-mapped from local disk rather than loaded into RAM, so
# --mem does not need to exceed hash.k2d. Copying once per job avoids
# random reads against BeeGFS and per-sample reloads.
LOCAL_ROOT="/tmp"
LOCAL_DB=""

cleanup() {
    if [[ -n "$LOCAL_DB" && -d "$LOCAL_DB" ]]; then
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Removing local database copy: ${LOCAL_DB}"
        rm -rf -- "$LOCAL_DB"
    fi
}
trap cleanup EXIT
trap 'exit 143' TERM   # e.g. scancel / time limit -> EXIT trap still runs
trap 'exit 130' INT

REQUIRED_BYTES=$(du -cb "${DB_FILES[@]}" | tail -1 | cut -f1)
AVAIL_BYTES=$(df --output=avail -B1 "$LOCAL_ROOT" | tail -1 | tr -d ' ')
echo "Node: $(hostname)"
echo "Database size: $((REQUIRED_BYTES / 1024**3)) GiB; free on ${LOCAL_ROOT}: $((AVAIL_BYTES / 1024**3)) GiB"
if (( AVAIL_BYTES < REQUIRED_BYTES * 105 / 100 )); then
    echo "ERROR: Not enough local disk on $(hostname):${LOCAL_ROOT} to stage the database"
    exit 1
fi

LOCAL_DB=$(mktemp -d "${LOCAL_ROOT}/k2db_${SLURM_JOB_ID:-local}_XXXX")
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Copying database to ${LOCAL_DB} ..."
SECONDS=0
for f in "${DB_FILES[@]}"; do
    cp -- "$f" "$LOCAL_DB/"
done
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Database staged in $((SECONDS / 60))m $((SECONDS % 60))s"
# ── Run ───────────────────────────────────────────────────────────────────────
passed=0
for i in "${!SAMPLE_IDS[@]}"; do
    ID="${SAMPLE_IDS[$i]}"
    R1="${SAMPLE_R1[$i]}"
    R2="${SAMPLE_R2[$i]}"
    OUT_DIR="${BASE_OUT}/${ID}"
    echo "────────────────────────────────────────────────────────────"
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Processing sample $((i + 1))/${n_run}: ${ID}"
    echo "  R1: ${R1}"
    echo "  R2: ${R2}"
    echo "  Output: ${OUT_DIR}"
    mkdir -p "$OUT_DIR"
    if k2 classify \
        --db "$LOCAL_DB" \
        --memory-mapping \
        --threads "$THREADS" \
        --output           "${OUT_DIR}/kraken_output.txt" \
        --report           "${OUT_DIR}/kraken_report.txt" \
        --classified-out   "${OUT_DIR}/classified#.fastq" \
        --unclassified-out "${OUT_DIR}/unclassified#.fastq" \
        --use-names \
        --paired \
        "$R1" "$R2"; then
        echo "  [$(date '+%Y-%m-%d %H:%M:%S')] SUCCESS: ${ID}"
        passed=$((passed + 1))
    else
        echo "  [$(date '+%Y-%m-%d %H:%M:%S')] FAILED: ${ID} (k2 classify exited non-zero)"
        failed_samples+=("${ID} (k2 classify failed)")
    fi
done
# ── Summary ───────────────────────────────────────────────────────────────────
echo "════════════════════════════════════════════════════════════"
echo "SUMMARY: ${passed}/${total} samples completed successfully"
if [[ ${#failed_samples[@]} -gt 0 ]]; then
    echo "FAILED samples:"
    for s in "${failed_samples[@]}"; do
        echo "  - ${s}"
    done
    exit 1
fi
echo "All samples completed successfully."
