#!/usr/bin/env bash
# Reconstruct the FeGenie.py used in this study.
#
# We do not redistribute FeGenie itself. This script fetches the exact upstream
# revision and applies the patch shipped alongside it.
#
#   FeGenie_asrun.py   - the version behind the published results (Fig. 9),
#                        described in scripts/local/fegenie/README.md
#
# Usage:  ./get_fegenie.sh [outdir]
# Then copy the file over MagicLamp's genies/FeGenie.py.

set -e

outdir="${1:-.}"
commit=fc7b6f7e3b99897bab744e9cc0363a6f6c46c459   # MagicLamp genies/FeGenie.py, 2022-11-25
url=https://raw.githubusercontent.com/Arkadiy-Garber/MagicLamp/${commit}/genies/FeGenie.py

here=$(cd "$(dirname "$0")" && pwd)

curl -sL "$url" -o "${outdir}/FeGenie_upstream.py"

# Verify the fetch, so a proxy error page never gets patched by mistake
echo "39d285dff9d84deebf1ae6f0c61c4323d649e43d977d5e763717f3de7831415b  ${outdir}/FeGenie_upstream.py" \
    | sha256sum -c -

cp "${outdir}/FeGenie_upstream.py" "${outdir}/FeGenie_asrun.py"
patch -p0 "${outdir}/FeGenie_asrun.py" < "${here}/FeGenie_asrun.patch"

python3 -m py_compile "${outdir}/FeGenie_asrun.py"
echo "OK: ${outdir}/FeGenie_asrun.py"
