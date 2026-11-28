#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
if [[ ! -f target/combined-review-1.0.0.jar ]]; then
  mvn -B package
fi
exec java -jar target/combined-review-1.0.0.jar "$@"
