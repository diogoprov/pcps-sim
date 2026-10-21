#!/bin/bash
# Keep a spot run alive: while the final results.csv of the profile is not in
# S3, relaunch the run (resuming from the targets store) whenever its instance
# ends, trying several 192-vCPU instance types. Prints one line per event.
# Usage: AWS_PROFILE=... aws/supervise.sh <bucket> <profile> <region> <instance_id> <max_hours> [max_relaunches]
set -uo pipefail
BUCKET=$1; PROFILE=$2; REGION=$3; ID=$4; HOURS=$5; MAXR=${6:-5}
cd "$(git rev-parse --show-toplevel)"
done_run() { aws s3 ls "s3://${BUCKET}/pcps-sim/results/${PROFILE}/${PROFILE}/results.csv" >/dev/null 2>&1; }
state() { aws ec2 describe-instances --region "$REGION" --instance-ids "$1" \
  --query 'Reservations[].Instances[].State.Name' --output text 2>/dev/null || echo gone; }
n=0
while true; do
  if done_run; then echo "${PROFILE}: finished ($(date))"; exit 0; fi
  s=$(state "$ID")
  case "$s" in
    running|pending) sleep 300; continue;;
  esac
  [ "$n" -ge "$MAXR" ] && { echo "${PROFILE}: gave up after ${n} relaunches ($(date))"; exit 1; }
  n=$((n + 1))
  for T in c7a.48xlarge c7i.48xlarge m7a.48xlarge c6a.48xlarge m7i.48xlarge; do
    out=$(AWS_REGION=$REGION AWS_DEFAULT_REGION=$REGION aws/launch_spot.sh "$BUCKET" "$PROFILE" "$T" 190 pcps-sim-ec2 "$HOURS" 2>&1 | tail -1)
    case "$out" in i-*) ID=$out; echo "${PROFILE}: instance ended (${s}); relaunch ${n} on ${T} ${ID} ($(date))"; break;; esac
  done
  sleep 300
done
