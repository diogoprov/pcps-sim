# Running the simulation pipeline on AWS EC2 (spot)

The pipeline (`_targets.R`, profiles in `R/pipeline.R`) runs inside a Docker image
(`Dockerfile`, `rocker/r-ver:4.6.1` + `renv.lock`) on a single spot instance. The
`targets` store is synced to S3 every 10 min and on a spot interruption notice,
so a relaunch resumes from the last checkpoint (only the branches that were
running are lost). The instance terminates itself at the end.

## One-time setup

1. Install and configure the AWS CLI (`aws configure`), region e.g. `us-east-1`.
2. Create a bucket, e.g. `aws s3 mb s3://<bucket>`.
3. Create an IAM role for EC2 with read/write access to that bucket and an
   instance profile with the same name (e.g. `pcps-sim-ec2`).
4. Check the spot vCPU quota for the instance family
   (Service Quotas > EC2 > "All Standard ... Spot Instance Requests"); a
   c7a.48xlarge needs 192 vCPUs.

## Each run

```bash
aws/pack_code.sh <bucket>                      # uploads git HEAD (commit first)
aws/launch_spot.sh <bucket> bench c7a.2xlarge 8 pcps-sim-ec2 1   # ~15 min, < US$ 1
# check s3://<bucket>/pcps-sim/results/bench/ and the log, then:
aws/launch_spot.sh <bucket> full c7a.48xlarge 190 pcps-sim-ec2 10
```

- Progress: `aws s3 ls s3://<bucket>/pcps-sim/logs/` (log copied at each checkpoint).
- Results: `aws s3 sync s3://<bucket>/pcps-sim/results/full data/derived/` and
  the store: `aws s3 sync s3://<bucket>/pcps-sim/store/_targets_full _targets_full`
  (then `targets::tar_read(results, store = "_targets_full")`).
- If the instance is interrupted, run `launch_spot.sh` again with the same profile.
- After the run: empty the bucket prefix if no longer needed (storage is charged).

## Cost control

- Time per replicate measured locally: ~2.6 / 3.1 / 7.9 s for N = 50 / 100 / 200
  (D_P tests + PCC, 999 permutations). Use the `bench` profile to measure it on
  the target instance before the full run.
- Prices change; check the current spot price for the instance type and region.
