#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "$0")" && pwd)"
cd "$ROOT"

export PYTHONHASHSEED=0
export SOURCE_DATE_EPOCH=1787011200

python3 -m py_compile \
  globalisation_holonomy.py \
  verify_A4_globalisation.py \
  verify_cross_rank_globalisation.py

python3 verify_A4_globalisation.py \
  --json A4_globalisation_certificate.json \
  --tex A4_globalisation_results.tex \
  | tee A4_globalisation_output.txt

python3 verify_cross_rank_globalisation.py \
  --json cross_rank_globalisation.json \
  --tex cross_rank_results.tex \
  | tee cross_rank_globalisation_output.txt

latexmk -C globalizable_hall_quantisations.tex >/dev/null 2>&1 || true
latexmk -pdf -interaction=nonstopmode -halt-on-error \
  globalizable_hall_quantisations.tex

if command -v latexpand >/dev/null 2>&1; then
  latexpand globalizable_hall_quantisations.tex \
    > globalizable_hall_quantisations_single.tex
else
  cp globalizable_hall_quantisations.tex \
    globalizable_hall_quantisations_single.tex
fi

FILES=(
  globalizable_hall_quantisations.tex
  globalizable_hall_quantisations_single.tex
  globalizable_hall_quantisations.pdf
  globalisation_holonomy.py
  verify_A4_globalisation.py
  verify_cross_rank_globalisation.py
  A4_globalisation_certificate.json
  cross_rank_globalisation.json
  A4_globalisation_results.tex
  cross_rank_results.tex
  A4_globalisation_output.txt
  cross_rank_globalisation_output.txt
  README_zh.md
  PROOF_STATUS_zh.md
  LITERATURE_NOTES_zh.md
  references.tex
  reproduce.sh
)
while IFS= read -r -d '' file; do
  FILES+=("${file#./}")
done < <(find sections -type f -name '*.tex' -print0 | sort -z)

sha256sum "${FILES[@]}" > MANIFEST_SHA256.txt

printf '\nReproduction completed.\n'
printf 'PDF: %s\n' "$ROOT/globalizable_hall_quantisations.pdf"
printf 'Manifest: %s\n' "$ROOT/MANIFEST_SHA256.txt"
