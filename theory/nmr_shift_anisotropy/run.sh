#!/usr/bin/env sh
set -eu
cd "$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)"
mkdir -p target/classes
printf '%s\n' src/main/java/org/eljah/nmr/*.java > target/sources.list
if [ -n "${JAVA_HOME:-}" ]; then
  NMR_JAVA="$JAVA_HOME/bin/java"
  NMR_JAVAC="$JAVA_HOME/bin/javac"
else
  NMR_JAVA=java
  NMR_JAVAC=javac
fi
if command -v "$NMR_JAVAC" >/dev/null 2>&1; then
  "$NMR_JAVAC" -encoding UTF-8 --release 17 -d target/classes @target/sources.list
else
  "$NMR_JAVA" -m jdk.compiler/com.sun.tools.javac.Main -encoding UTF-8 --release 17 -d target/classes @target/sources.list
fi
exec "$NMR_JAVA" -Xmx1500m -cp target/classes org.eljah.nmr.NmrShiftAnisotropy "$@"
