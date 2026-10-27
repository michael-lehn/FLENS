#! /bin/sh
#
#  check-exports.sh LIBRARY LIST
#
#  Checks that the Fortran-callable routines exported by LIBRARY (all
#  defined global symbols of the form name_) are exactly the routines listed
#  in LIST.  Lines of LIST:
#
#    name           routine name_ must be exported
#    !name reason   routine name_ is not implemented and must not be exported
#
#  Empty lines and lines starting with '#' are ignored.
#

if [ $# -ne 2 ]; then
    echo "usage: $0 LIBRARY LIST" >&2
    exit 2
fi

lib=$1
list=$2
NM=${NM:-nm}

tmp=${TMPDIR:-/tmp}/check-exports.$$
trap 'rm -f $tmp.*' 0 1 2 15

#  Mach-O prefixes C symbols with '_'
strip=cat
if [ "$(uname -s)" = Darwin ]; then
    strip="sed s/^_//"
fi

$NM -g "$lib" 2>/dev/null \
    | awk 'NF>=3 && $2 ~ /^[TDSC]$/ { print $3 }' \
    | $strip \
    | grep -E '^_*[a-z][a-z0-9_]*_$' \
    | sort -u > $tmp.exported

grep -vE '^[[:space:]]*(#|$)' "$list" | grep -v '^!' \
    | awk '{ print $1 "_" }' | sort -u > $tmp.supported
grep '^!' "$list" | sed 's/^!//' \
    | awk '{ print $1 "_" }' | sort -u > $tmp.unsupported

if [ ! -s $tmp.exported ]; then
    echo "check-exports: no routines found in $lib" >&2
    exit 1
fi

status=0
for name in $(comm -13 $tmp.supported $tmp.exported); do
    if grep -qx "$name" $tmp.unsupported; then
        echo "check-exports: $lib exports unsupported routine $name"
    else
        echo "check-exports: $lib exports $name, which is not listed in $list"
    fi
    status=1
done
for name in $(comm -23 $tmp.supported $tmp.exported); do
    echo "check-exports: $lib does not export $name"
    status=1
done

if [ $status -eq 0 ]; then
    echo "$lib exports $(wc -l < $tmp.exported | tr -d ' ') routines" \
         "as listed in $list" \
         "($(wc -l < $tmp.unsupported | tr -d ' ') listed as not implemented)"
fi
exit $status
