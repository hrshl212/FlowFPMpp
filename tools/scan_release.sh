#!/usr/bin/env bash
# Pre-publication scan of a repository tree: secrets, site paths and accounts,
# personal identifiers, binaries, and large files. Exits nonzero on findings.
# Site-specific names (accounts, user names, institutions) are not stored here;
# pass them as an extended regex in SCAN_EXTRA_PATTERNS.
set -uo pipefail
dir=${1:?usage: tools/scan_release.sh DIR}
status=0
report() { echo "FINDING: $1"; status=1; }
cd "$dir"
g() { grep -rInE --exclude=scan_release.sh --exclude-dir=.git "$@" .; }
if g '(password|passwd|secret|token|api[_-]?key|BEGIN [A-Z ]*PRIVATE KEY)'; then report "possible secret"; fi
if g '/scratch|/home/|/data/apps|/users/|@[A-Za-z0-9.-]+\.(edu|com|org)'; then report "absolute site path or e-mail address"; fi
if [ -n "${SCAN_EXTRA_PATTERNS:-}" ] && g "$SCAN_EXTRA_PATTERNS"; then report "site-specific pattern"; fi
# (--include must come first: grep keeps unmatched files when the first rule is --exclude.)
if grep -rInE --include='*.f90' --include=Makefile --include='*.example' --include='*.yml' 'private_inputs' .; then
  report "private data reference in code or scripts"; fi
big=$(find . -type f -not -path './.git/*' -size +1M); [ -n "$big" ] && { echo "$big"; report "file larger than 1 MB"; }
# PNG figures are allowed under docs/images; everything else must be text.
bin=$(find . -type f -not -path './.git/*' -not -path './docs/images/*.png' -exec file {} + | grep -viE 'text|empty' || true); [ -n "$bin" ] && { echo "$bin"; report "non-text file"; }
[ $status = 0 ] && echo "scan clean: $(find . -type f | wc -l) files"
exit $status
