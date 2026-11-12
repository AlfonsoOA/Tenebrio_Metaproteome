#!/bin/bash
set -euo pipefail

# ========= 📂 0. DIRECTORY STRUCTURE (THE ANCHOR) =========
# Root of the annotation workspace (defaults to the user home)
BASE_DIR="${HOME}/FASTA_Annotations"
INPUT_DIR="${BASE_DIR}/1_FASTA_Files"
CLEAN_DIR="${BASE_DIR}/2_FASTA_Clean"
IPS_DIR="${BASE_DIR}/3_InterProScan_Annotation"
EGG_DIR="${BASE_DIR}/4_EggNOG_Annotation"
RESULTS_DIR="${BASE_DIR}/5_Annotation_results"

# --- Optional: destination folder to copy the annotation results back to ---
# Set this to the output folder of the R analysis script to synchronise results.
# Leave empty to skip the copy step.
WIN_DEST=""

# Internal Paths
EGG_DB_DATA="${EGG_DIR}/data"
IPS_EXE="${IPS_DIR}/interproscan/interproscan.sh"
EGG_EXE="${EGG_DIR}/emapper.py"

# Create folders ONLY if they don't exist
mkdir -p "$INPUT_DIR" "$CLEAN_DIR" "$IPS_DIR" "$EGG_DIR" "$RESULTS_DIR"

# ========= 🛠️ 1. PRE-FLIGHT CHECKS =========

# Check for Biopython (Required for EggNOG)
if ! python3 -c "import Bio" &>/dev/null; then
    echo "📦 Biopython missing. Installing via APT..."
    sudo apt update && sudo apt install python3-biopython -y
else
    echo "✅ Biopython is ready."
fi

# ========= 🚀 2. MAIN LOOP =========

shopt -s nullglob
FASTAS=("$INPUT_DIR"/*.faa "$INPUT_DIR"/*.fasta)

if [ ${#FASTAS[@]} -eq 0 ]; then
    echo "❌ No files found in $INPUT_DIR. Drop your .faa files there to start!"
    exit 1
fi

for f in "${FASTAS[@]}"; do
    RAW_NAME=$(basename "$f")
    ID="${RAW_NAME%.*}"
    
    echo -e "\n"
    echo "########################################################"
    echo "  🧬 PROCESSING SAMPLE: $ID"
    echo "########################################################"

    # A. CLEANING MODULE
    CLEAN_FASTA="${CLEAN_DIR}/${ID}_cleaned.faa"
    if [[ ! -f "$CLEAN_FASTA" ]]; then
        echo "🧼 Cleaning FASTA (removing '*' stop codons)..."
        sed 's/\*//g' "$f" > "$CLEAN_FASTA"
    else
        echo "✅ Clean FASTA already exists for $ID."
    fi

    # B. EggNOG-MAPPER MODULE (With Anchor Check)
    EGG_FINAL="${RESULTS_DIR}/EggNOG_${ID}.emapper.annotations"
    if [[ ! -f "$EGG_FINAL" ]]; then
        # Check if DB is already in the anchor folder
        if [[ -f "${EGG_DB_DATA}/eggnog.db" ]]; then
            echo "✅ EggNOG DB found in anchor folder. Skipping download."
            echo "🥚 Running EggNOG-mapper..."
            python3 "$EGG_EXE" \
                -i "$CLEAN_FASTA" \
                -o "EggNOG_${ID}" \
                --output_dir "$RESULTS_DIR" \
                --cpu 7 -m diamond --data_dir "$EGG_DB_DATA" --override
        else
            echo "📡 EggNOG DB not found in $EGG_DB_DATA."
            echo "💡 Please move your downloaded 'data' folder to $EGG_DIR before running."
        fi
    else
        echo "✅ EggNOG results already exist for $ID. Skipping."
    fi

    # C. INTERPROSCAN MODULE (With Anchor Check)
    IPS_FINAL="${RESULTS_DIR}/InterProScan_${ID}.tsv"
    if [[ ! -f "$IPS_FINAL" ]]; then
        if [[ -f "$IPS_EXE" ]]; then
            echo "✅ InterProScan found in anchor folder."
            echo "🔍 Running InterProScan Analysis..."
            "$IPS_EXE" -i "$CLEAN_FASTA" -f tsv -cpu 7 --goterms --pathways --iprlookup -b "${RESULTS_DIR}/InterProScan_${ID}"
            # InterProScan appends the format extension to the -b basename,
            # so the output should already be at $IPS_FINAL. Verify it.
            if [[ ! -f "$IPS_FINAL" ]]; then
                echo "⚠️  Expected output not found: $IPS_FINAL"
            fi
        else
            echo "📡 InterProScan not found in $IPS_DIR."
            echo "💡 Please install/move InterProScan to $IPS_DIR before running."
        fi
    else
        echo "✅ InterProScan results already exist for $ID. Skipping."
    fi

    echo "✨ Sample $ID finished successfully!"
done

# ========= 🚛 3. WIN-SYNCHRONIZATION (THE BRIDGE) =========
echo -e "\n"
echo "########################################################"
echo "  🚛 SYNCHRONISING RESULTS"
echo "########################################################"

# Copy only the annotation outputs to the R analysis results folder
if [[ -n "$WIN_DEST" && -d "$WIN_DEST" ]]; then
    cp "${RESULTS_DIR}"/*.emapper.annotations "$WIN_DEST/" 2>/dev/null || true
    cp "${RESULTS_DIR}"/*.tsv "$WIN_DEST/" 2>/dev/null || true
    echo "✅ Results synchronized to: $WIN_DEST"
else
    echo "ℹ️  No destination folder set (WIN_DEST). Skipping synchronisation."
fi
echo -e "\n🏁 [PIPELINE COMPLETE] Check your results in $RESULTS_DIR"