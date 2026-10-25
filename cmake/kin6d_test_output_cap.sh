#!/usr/bin/env bash
# Usage: kin6d_test_output_cap.sh LIMIT_BYTES COMMAND [ARGS...]
#
# CTest test launcher. CTest stores the complete output of every test in
# Testing/Temporary/LastTest.log, and its CTEST_CUSTOM_MAXIMUM_*_OUTPUT_SIZE
# limits do not apply to that file. A debug-build runtime-warning flood once
# wrote 128-143 GB there (decision 54). This launcher keeps the first
# LIMIT_BYTES of the merged stdout/stderr and the last 200 lines, and reports
# how much was omitted. COMMAND replaces this shell (exec), so CTest sees its
# own exit status, signals and timeouts unchanged.
limit=$1
shift
exec "$@" > >(LC_ALL=C awk -v limit="$limit" -v keep=200 '
{
    len = length($0) + 1
    if (kept + len <= limit) {
        print
        kept += len
        next
    }
    omitted += len
    n++
    tail[n % keep] = substr($0, 1, 4096)
}
END {
    if (n == 0) exit
    first = n > keep ? n - keep + 1 : 1
    printf "[kin6d test output cap: %d bytes (%d lines) beyond the first %d bytes; the last %d lines follow]\n", omitted, n, kept, n - first + 1
    for (i = first; i <= n; i++) print tail[i % keep]
}') 2>&1
