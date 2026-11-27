#!/usr/bin/env sh
set -eu
cd -- "$(dirname -- "$0")"

# Overrides are useful when Java is installed outside PATH.
java_bin=${DIGENITE_JAVA:-}
if [ -z "$java_bin" ]; then
    if command -v java >/dev/null 2>&1; then
        java_bin=$(command -v java)
    elif [ -x /usr/lib/jvm/java-21-openjdk-amd64/bin/java ]; then
        java_bin=/usr/lib/jvm/java-21-openjdk-amd64/bin/java
    else
        echo 'Java 17+ is required; set DIGENITE_JAVA to its java executable.' >&2
        exit 1
    fi
fi

mkdir -p target/classes
javac_bin=${DIGENITE_JAVAC:-}
if [ -z "$javac_bin" ] && command -v javac >/dev/null 2>&1; then
    javac_bin=$(command -v javac)
fi
if [ -n "$javac_bin" ]; then
    "$javac_bin" --release 17 -d target/classes src/main/java/org/eljah/nmr/DigeniteEfg.java
else
    # The supplied runtime has jdk.compiler but no javac executable or ct.sym.
    # -source/-target 17 produces version-61 class files; the source only uses
    # Java 17 APIs. A normal JDK/Maven build uses --release 17 instead.
    "$java_bin" -m jdk.compiler/com.sun.tools.javac.Main -source 17 -target 17 \
        -d target/classes src/main/java/org/eljah/nmr/DigeniteEfg.java
fi

exec "$java_bin" -Xmx512m -cp target/classes org.eljah.nmr.DigeniteEfg \
    --config=model.properties "$@"
