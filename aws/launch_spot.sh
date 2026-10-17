#!/bin/bash
# Launch one spot instance that runs the pipeline and terminates at the end.
# Usage: aws/launch_spot.sh <bucket> <profile> <instance_type> <workers> <instance_profile>
#   e.g. aws/launch_spot.sh my-bucket bench c7a.2xlarge 8 pcps-sim-ec2
#        aws/launch_spot.sh my-bucket full c7a.48xlarge 190 pcps-sim-ec2
# <instance_profile>: IAM instance profile with read/write access to the bucket.
set -euo pipefail
BUCKET=$1; PROFILE=$2; TYPE=$3; WORKERS=$4; ROLE=$5
cd "$(dirname "$0")"
sed -e "s/__BUCKET__/${BUCKET}/" -e "s/__PROFILE__/${PROFILE}/" -e "s/__WORKERS__/${WORKERS}/" \
  user_data.sh > /tmp/pcps_user_data.sh

aws ec2 run-instances \
  --image-id resolve:ssm:/aws/service/ami-amazon-linux-latest/al2023-ami-kernel-default-x86_64 \
  --instance-type "${TYPE}" \
  --instance-market-options '{"MarketType":"spot","SpotOptions":{"SpotInstanceType":"one-time","InstanceInterruptionBehavior":"terminate"}}' \
  --instance-initiated-shutdown-behavior terminate \
  --iam-instance-profile Name="${ROLE}" \
  --block-device-mappings '[{"DeviceName":"/dev/xvda","Ebs":{"VolumeSize":60,"VolumeType":"gp3","DeleteOnTermination":true}}]' \
  --user-data file:///tmp/pcps_user_data.sh \
  --tag-specifications "ResourceType=instance,Tags=[{Key=Name,Value=pcps-sim-${PROFILE}}]" \
  --query 'Instances[0].InstanceId' --output text
