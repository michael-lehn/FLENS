#!/bin/sh
#
#  Evaluate the output files of the Netlib LAPACK test suite.
#
#  usage: check-results.sh EXPECTED OUTFILE...
#
#  The test programs always exit with status 0, so the results are taken
#  from their output.  A test passes if it ran to completion ("End of
#  tests"), no test ratio exceeded the threshold and no routine reported an
#  illegal argument.
#
#  EXPECTED lists known deviations (lines "file|message|reason").  A test
#  whose deviations match the listed ones exactly is reported as EXPECTED;
#  every other deviation is a failure.  "message" is either a line of the
#  output (with whitespace collapsed, prefixed by the number of occurrences)
#  or INCOMPLETE for a test that did not run to completion.
#
#  Exit status: 0 if there were no unexpected results, 1 otherwise.
#

expected=$1
shift

pass=0
known=0
failed=0
total=0

for out in "$@"; do
    name=$(basename "$out")

    if [ -s "$out" ] && grep -aq 'End of tests' "$out"; then
        actual=$(grep -aE 'failed to pass the threshold|illegal value' "$out" \
                 | sed -e 's/  */ /g' -e 's/^ //' -e 's/ $//' \
                 | sort | uniq -c | sed -e 's/^ *//')
    else
        actual=INCOMPLETE
    fi
    listed=$(grep -v '^#' "$expected" 2>/dev/null \
             | awk -F'|' -v f="$name" '$1==f { print $2 }' | sort)
    reason=$(grep -v '^#' "$expected" 2>/dev/null \
             | awk -F'|' -v f="$name" '$1==f { print $3; exit }')

    n=$(grep -aoE '\( *[0-9]+ tests run\)' "$out" 2>/dev/null \
        | grep -oE '[0-9]+' | awk '{ s += $1 } END { print s+0 }')
    total=$((total + n))

    if [ -z "$actual" ]; then
        if [ -n "$listed" ]; then
            echo "PASS      $name ($n tests; listed deviation no longer occurs)"
        else
            echo "PASS      $name ($n tests)"
        fi
        pass=$((pass + 1))
    elif [ "$actual" = "$(echo "$listed" | sort)" ]; then
        echo "EXPECTED  $name ($n tests): $reason"
        echo "$actual" | sed -e 's/^/            /'
        known=$((known + 1))
    else
        echo "FAIL      $name ($n tests)"
        echo "$actual" | sed -e 's/^/            /'
        failed=$((failed + 1))
    fi
done

echo "Netlib LAPACK tests: $pass passed, $known expected deviations," \
     "$failed failed; $total tests run in completed blocks"
test $failed -eq 0
