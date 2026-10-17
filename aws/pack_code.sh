#!/bin/bash
# Upload the committed code (git HEAD) plus the untracked Portik et al. 2023 tree
# (data/raw/README.md) to S3.
# Usage: aws/pack_code.sh <bucket>
set -euo pipefail
BUCKET=$1
cd "$(git rev-parse --show-toplevel)"
if ! git diff --quiet HEAD; then echo "Uncommitted changes; commit first." >&2; exit 1; fi
TREE=data/raw/portik2023_timetree.tre
echo "d871f3643575709e846f167ce297834a63afe21743296902d1548db4b4691e90  ${TREE}" | shasum -a 256 -c -
git archive --format=tar -o /tmp/pcps-sim.tar HEAD
tar -rf /tmp/pcps-sim.tar "${TREE}"
gzip -f /tmp/pcps-sim.tar
aws s3 cp /tmp/pcps-sim.tar.gz s3://${BUCKET}/pcps-sim/code/pcps-sim.tar.gz
echo "Uploaded $(git rev-parse --short HEAD) to s3://${BUCKET}/pcps-sim/code/"
