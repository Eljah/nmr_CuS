#!/usr/bin/env sh
set -eu
cd -- "$(dirname -- "$0")"

# Optional overrides point to executables, not to their parent directories.
nqr_java=${SUPERCELL_JAVA:-}
if [ -z "$nqr_java" ]; then
    if command -v java >/dev/null 2>&1; then
        nqr_java=$(command -v java)
    elif [ -x /usr/lib/jvm/java-21-openjdk-amd64/bin/java ]; then
        nqr_java=/usr/lib/jvm/java-21-openjdk-amd64/bin/java
    else
        echo 'Java 17+ is required; set SUPERCELL_JAVA to its java executable.' >&2
        exit 1
    fi
fi

mkdir -p target/classes
nqr_javac=${SUPERCELL_JAVAC:-}
if [ -z "$nqr_javac" ] && command -v javac >/dev/null 2>&1; then
    nqr_javac=$(command -v javac)
fi
if [ -n "$nqr_javac" ]; then
    "$nqr_javac" --release 17 -encoding UTF-8 -d target/classes \
        src/main/java/org/eljah/nmr/SupercellNqr.java \
        src/main/java/org/eljah/nmr/TemperatureSeries.java
else
    # This workspace provides jdk.compiler, but no javac launcher or ct.sym.
    # The source uses only Java 17 APIs; -source/-target emits class version 61.
    "$nqr_java" -m jdk.compiler/com.sun.tools.javac.Main -source 17 -target 17 \
        -encoding UTF-8 -d target/classes src/main/java/org/eljah/nmr/SupercellNqr.java \
        src/main/java/org/eljah/nmr/TemperatureSeries.java
fi

"$nqr_java" -cp target/classes org.eljah.nmr.SupercellNqr "$@"
# The default run also processes the archived temperature series, when present.
# Explicit arguments remain arguments to SupercellNqr; TemperatureSeries also
# has its own CLI, callable with java -cp target/classes org.eljah.nmr.TemperatureSeries.
if [ "$#" -eq 0 ] && [ -f inputs/temperature_series.csv ]; then
    exec "$nqr_java" -cp target/classes org.eljah.nmr.TemperatureSeries
fi
