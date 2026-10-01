#!/bin/sh
# Builds PHANTAST_Live-1.1.jar against the ImageJ library of a Fiji installation.
# Needs a Java compiler (javac, version 8 or newer) on the PATH.
# Usage:  ./build.sh /path/to/Fiji.app
set -e

FIJI="$1"
if [ -z "$FIJI" ]; then
    echo "Usage: ./build.sh /path/to/Fiji.app"
    exit 1
fi

IJJAR=$(ls "$FIJI"/jars/ij-*.jar 2>/dev/null | tail -n 1)
if [ -z "$IJJAR" ]; then
    echo "Could not find ij-*.jar in $FIJI/jars"
    exit 1
fi

HERE=$(cd "$(dirname "$0")" && pwd)
OUT="$HERE/build"
rm -rf "$OUT"
mkdir -p "$OUT"

javac --release 8 -Xlint:-options -cp "$IJJAR" -d "$OUT" "$HERE/src/PHANTAST_Live.java"
cp "$HERE/plugins.config" "$HERE/src/PHANTAST_Live.java" "$OUT/"
[ -f "$HERE/../LICENSE" ] && cp "$HERE/../LICENSE" "$OUT/LICENSE.txt"
jar cf "$HERE/PHANTAST_Live-1.1.jar" -C "$OUT" .

echo "Built $HERE/PHANTAST_Live-1.1.jar"
echo "Copy it into $FIJI/plugins and restart Fiji."
