#!/bin/bash
# Launch one spot instance that runs the pipeline and terminates at the end.
# Usage: aws/launch_spot.sh <bucket> <profile> <instance_type> <workers> <instance_profile> [max_hours]
#   e.g. aws/launch_spot.sh my-bucket bench c7a.2xlarge 8 pcps-sim-ec2 1
#        aws/launch_spot.sh my-bucket full c7a.48xlarge 190 pcps-sim-ec2 10
# max_hours: hard limit after which the instance terminates (default 12).
# <instance_profile>: IAM instance profile with read/write access to the bucket.
set -euo pipefail
BUCKET=$1; PROFILE=$2; TYPE=$3; WORKERS=$4; ROLE=$5; MAX_HOURS=${6:-12}
cd "$(dirname "$0")"
sed -e "s/__BUCKET__/${BUCKET}/" -e "s/__PROFILE__/${PROFILE}/" -e "s/__WORKERS__/${WORKERS}/" \
  -e "s/__MAX_HOURS__/${MAX_HOURS}/" \
  user_data.sh > /tmp/pcps_user_data.sh

# latest Amazon Linux 2023 AMI (x86_64), looked up without SSM permissions
AMI=$(aws ec2 describe-images --owners amazon \
  --filters "Name=name,Values=al2023-ami-2023.*-kernel-*-x86_64" "Name=state,Values=available" \
  --query 'sort_by(Images, &CreationDate)[-1].ImageId' --output text)
echo "AMI ${AMI}" >&2

aws ec2 run-instances \
  --image-id "${AMI}" \
  --instance-type "${TYPE}" \
  --instance-market-options '{"MarketType":"spot","SpotOptions":{"SpotInstanceType":"one-time","InstanceInterruptionBehavior":"terminate"}}' \
  --instance-initiated-shutdown-behavior terminate \
  --iam-instance-profile Name="${ROLE}" \
  --block-device-mappings '[{"DeviceName":"/dev/xvda","Ebs":{"VolumeSize":60,"VolumeType":"gp3","DeleteOnTermination":true}}]' \
  --user-data file:///tmp/pcps_user_data.sh \
  --tag-specifications "ResourceType=instance,Tags=[{Key=Name,Value=pcps-sim-${PROFILE}}]" \
  --query 'Instances[0].InstanceId' --output text
