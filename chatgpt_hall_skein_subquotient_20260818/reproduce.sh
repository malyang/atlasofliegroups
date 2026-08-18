#!/usr/bin/env bash
set -euo pipefail

export PYTHONHASHSEED=0
export SOURCE_DATE_EPOCH=${SOURCE_DATE_EPOCH:-1787011200}

mkdir -p generated

python -m py_compile \
  hall_skein.py \
  verify_A4_hall_skein.py \
  verify_typeA_hall_skein.py

python verify_A4_hall_skein.py \
  --json generated/A4_hall_skein_certificate.json \
  --tex generated/A4_hall_skein_results.tex \
  | tee generated/A4_hall_skein_output.txt

python verify_typeA_hall_skein.py \
  --json generated/typeA_hall_skein_certificate.json \
  --tex generated/typeA_hall_skein_results.tex \
  | tee generated/typeA_hall_skein_output.txt

latexmk -pdf -file-line-error -interaction=nonstopmode -halt-on-error \
  hall_skein_subquotient.tex

LOG=hall_skein_subquotient.log
PDF=hall_skein_subquotient.pdf

test -s "$PDF"
if grep -Eq 'LaTeX Warning: (Reference|Citation).*undefined|There were undefined references|multiply defined' "$LOG"; then
  echo "Unresolved reference or citation found" >&2
  exit 1
fi
if grep -Eq 'Overfull \\hbox|Overfull \\vbox' "$LOG"; then
  echo "Overfull box found" >&2
  exit 1
fi

{
  echo "Python syntax: pass"
  echo "A4 polygon verifier: pass"
  echo "Cross-rank polygon verifier: pass"
  echo "LaTeX compilation: pass"
  echo "Undefined references: 0"
  echo "Overfull boxes: 0"
} > QA_REPORT.txt

if command -v pdfinfo >/dev/null 2>&1; then
  pdfinfo "$PDF" > generated/pdfinfo.txt
  PAGES=$(awk '/^Pages:/ {print $2}' generated/pdfinfo.txt)
  echo "PDF pages: $PAGES" >> QA_REPORT.txt
fi
if command -v pdffonts >/dev/null 2>&1; then
  pdffonts "$PDF" > generated/pdffonts.txt
fi
if command -v gs >/dev/null 2>&1; then
  gs -q -dNOPAUSE -dBATCH -sDEVICE=nullpage "$PDF"
  echo "Ghostscript parse: pass" >> QA_REPORT.txt
fi

printf 'Reproduction complete: %s\n' "$PDF"
