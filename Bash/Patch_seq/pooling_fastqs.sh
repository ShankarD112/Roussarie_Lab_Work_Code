#----- variables you already have --------------------------------------------
RAW_DIR="path"
LINK_DIR="path"
SHEET="path"
mkdir -p "$LINK_DIR"

#----- 1. build an associative array of sample->barcode from the sheet -------
declare -A BC
while IFS=, read -r samp bc; do
  [[ $samp == "Sample" ]] && continue   # skip header
  BC[$samp]=$bc
done < "$SHEET"

#----- 2. loop over every hashed folder and make symlinks --------------------
shopt -s nullglob
for d in "$RAW_DIR"/AM*_ds.*; do
    samp=$(basename "$d" | awk -F_ '{print $1}')   # e.g. AM92
    tag=${samp#AM}                                 # e.g. 92  (used to keep files unique)
    for fq in "$d"/${samp}_S*_L*_R[12]_001.fastq.gz; do
        # original file name parts
        lane=$(echo "$fq" | sed -n 's/.*_L\([0-9]\{3\}\)_.*/\1/p')
        read=$(echo "$fq" | sed -n 's/.*_R\([12]\)_001.fastq.gz/\1/p')

        # new link name: AM_<tag> keeps uniqueness, rest identical
        ln -s "$fq" "$LINK_DIR/AM${tag}_L${lane}_R${read}.fastq.gz"
    done
done
