#!/bin/bash
# Upload the committed code (git HEAD, including data/raw) to S3.
# Usage: aws/pack_code.sh <bucket>
set -euo pipefail
BUCKET=$1
cd "$(git rev-parse --show-toplevel)"
if ! git diff --quiet HEAD; then echo "Uncommitted changes; commit first." >&2; exit 1; fi
git archive --format=tar.gz -o /tmp/pcps-sim.tar.gz HEAD
aws s3 cp /tmp/pcps-sim.tar.gz s3://${BUCKET}/pcps-sim/code/pcps-sim.tar.gz
echo "Uploaded $(git rev-parse --short HEAD) to s3://${BUCKET}/pcps-sim/code/"
