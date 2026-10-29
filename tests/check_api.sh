#!/bin/sh
# check_api.sh - detect API drift: every public symbol referenced/declared
# in include/*.h must be defined in the built static library, unless listed
# in the allowlist (Phase 2/3 removals and not-yet-implemented declarations).
#
# Usage: sh tests/check_api.sh <library.a> <allowlist.txt>

set -eu

lib=${1:-build/libbignum.a}
allow=${2:-tests/api_allow.txt}

if [ ! -f "$lib" ]; then
    echo "check-api: library not found: $lib" >&2
    exit 2
fi

tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT

# LTO objects need the plugin-aware nm wrapper when available.
NM="${CC:-cc}-nm"
command -v "$NM" >/dev/null 2>&1 || NM=nm
"$NM" --defined-only "$lib" | awk 'NF >= 3 {print $NF}' | sort -u \
    > "$tmp/defined"

# Domain-prefixed identifiers written call-style (`name(`) in public headers:
# covers real declarations plus Doxygen cross-references (`name()`), which
# keeps docs and API in sync. Prose like "u64s (" does not match.
grep -hoE '\b(rz|rmatq|rmat|rvec|rpol|rq|ec|u64|u128)[a-z0-9_]*\(' \
    include/*.h | sed 's/($//' | sort -u > "$tmp/declared"

: > "$tmp/allow"
if [ -f "$allow" ]; then
    sed 's/#.*//' "$allow" | tr -d ' \t\r' | grep -v '^$' | sort -u \
        > "$tmp/allow"
fi

comm -23 "$tmp/declared" "$tmp/defined" > "$tmp/missing_all"
comm -23 "$tmp/missing_all" "$tmp/allow" > "$tmp/missing"

if [ -s "$tmp/missing" ]; then
    echo "check-api: declared but undefined public symbols:" >&2
    sed 's/^/  /' "$tmp/missing" >&2
    echo "check-api: implement them or add to $allow" >&2
    exit 1
fi

echo "check-api: OK ($(wc -l < "$tmp/declared" | tr -d ' ') public names checked, \
$(wc -l < "$tmp/missing_all" | tr -d ' ') allowlisted)"
