#!/bin/bash
# EC2 user data (Amazon Linux 2023, x86_64). Placeholders __BUCKET__, __PROFILE__
# and __WORKERS__ are filled in by aws/launch_spot.sh.
# Resumes from the targets store saved in S3, checkpoints it every 10 min and on
# a spot interruption notice, and terminates the instance when the run ends.
set -uxo pipefail
BUCKET=__BUCKET__
PROFILE=__PROFILE__
WORKERS=__WORKERS__
STORE=_targets_${PROFILE}
WORK=/opt/pcps
S3=s3://${BUCKET}/pcps-sim

dnf install -y docker
systemctl start docker
mkdir -p ${WORK}
aws s3 cp ${S3}/code/pcps-sim.tar.gz /tmp/pcps-sim.tar.gz
tar xzf /tmp/pcps-sim.tar.gz -C ${WORK}
mkdir -p ${WORK}/${STORE} ${WORK}/data/derived
aws s3 sync ${S3}/store/${STORE} ${WORK}/${STORE} --only-show-errors   # resume

checkpoint() {
  aws s3 sync ${WORK}/${STORE} ${S3}/store/${STORE} --only-show-errors
  aws s3 sync ${WORK}/data/derived ${S3}/results/${PROFILE} --only-show-errors
  aws s3 cp ${WORK}/run.log ${S3}/logs/${PROFILE}_$(date +%Y%m%dT%H%M%S).log --only-show-errors || true
}

# periodic checkpoint
( while true; do sleep 600; checkpoint; done ) &

# spot interruption notice (2-minute warning): checkpoint once and stop polling
(
  while true; do
    TOKEN=$(curl -s -X PUT http://169.254.169.254/latest/api/token \
      -H "X-aws-ec2-metadata-token-ttl-seconds: 300")
    CODE=$(curl -s -o /dev/null -w "%{http_code}" -H "X-aws-ec2-metadata-token: ${TOKEN}" \
      http://169.254.169.254/latest/meta-data/spot/instance-action)
    if [ "${CODE}" = "200" ]; then checkpoint; break; fi
    sleep 5
  done
) &

cd ${WORK}
docker build -t pcps-sim . > ${WORK}/build.log 2>&1
docker run --rm \
  -v ${WORK}/${STORE}:/project/${STORE} \
  -v ${WORK}/data/derived:/project/data/derived \
  pcps-sim Rscript analysis/run_pipeline.R ${PROFILE} ${WORKERS} > ${WORK}/run.log 2>&1
echo "exit status $? at $(date)" >> ${WORK}/run.log
cat ${WORK}/build.log >> ${WORK}/run.log
checkpoint
shutdown -h now
