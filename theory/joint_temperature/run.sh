#!/usr/bin/env sh
set -eu
cd -- "$(dirname -- "$0")"
joint_java=${JOINT_JAVA:-}
if [ -z "$joint_java" ]; then
  if command -v java >/dev/null 2>&1; then joint_java=$(command -v java)
  elif [ -x /usr/lib/jvm/java-21-openjdk-amd64/bin/java ]; then joint_java=/usr/lib/jvm/java-21-openjdk-amd64/bin/java
  else echo 'Java 17+ is required; set JOINT_JAVA to its java executable.' >&2; exit 1; fi
fi
mkdir -p target/classes
if command -v javac >/dev/null 2>&1; then
  javac --release 17 -encoding UTF-8 -d target/classes src/main/java/org/eljah/nmr/joint/*.java
else
  # The execution environment contains jdk.compiler, but no javac launcher/ct.sym.
  # The code uses Java 17 APIs and emits class-file version 61.
  "$joint_java" -m jdk.compiler/com.sun.tools.javac.Main -source 17 -target 17 -encoding UTF-8 -d target/classes src/main/java/org/eljah/nmr/joint/*.java
fi
exec "$joint_java" -cp target/classes org.eljah.nmr.joint.JointTemperature "$@"
