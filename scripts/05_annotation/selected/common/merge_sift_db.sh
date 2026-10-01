#!/usr/bin/env bash
set -euo pipefail
: "${GENOME_DIR:?}" "${ORG_VERSION:?}" "${DB_DIR:?}"
[[ "$GENOME_DIR" = /* && "$DB_DIR" = /* ]] || { echo 'Use absolute paths' >&2; exit 1; }
list=${CHR_LIST:-$GENOME_DIR/configs/chr_allow.list}
mapfile -t chrs < <(awk 'NF && $1 !~ /^#/ {print $1}' "$list")
(( ${#chrs[@]} > 0 )) || exit 1
for chr in "${chrs[@]}"; do [[ -f "$GENOME_DIR/build_$chr/.SUCCESS" ]] || { echo "Incomplete: $chr" >&2; exit 1; }; done
dirs=(fasta subst SIFT_alignments SIFT_predictions singleRecords_with_scores logs "$ORG_VERSION")
for dir in "${dirs[@]}"; do mkdir -p "$DB_DIR/$dir"; done
for chr in "${chrs[@]}"; do
  build=$GENOME_DIR/build_$chr
  for dir in "${dirs[@]}"; do
    [[ -d "$build/$dir" ]] || continue
    while IFS= read -r -d '' file; do
      target=$DB_DIR/$dir/$(basename "$file")
      if [[ -L "$target" && "$(readlink "$target")" == "$file" ]]; then continue; fi
      if [[ -e "$target" || -L "$target" ]]; then echo "Duplicate filename: $target" >&2; exit 1; fi
      ln -s "$file" "$target"
    done < <(find "$build/$dir" -type f -print0)
  done
done
# The source Ped/Zmay merge appended all_prot.fasta. Rebuild it once to keep
# reruns idempotent; Atau may not have this file.
tmp=$(mktemp "$DB_DIR/all_prot.fasta.XXXXXX")
found=0
for chr in "${chrs[@]}"; do
  file=$GENOME_DIR/build_$chr/all_prot.fasta
  if [[ -s "$file" ]]; then cat "$file" >> "$tmp"; found=1; fi
done
if (( found )); then mv "$tmp" "$DB_DIR/all_prot.fasta"; else rm -f "$tmp"; fi
echo "Merged database: $DB_DIR"
