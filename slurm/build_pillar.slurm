#!/bin/bash
#SBATCH --job-name=tobe-build
#SBATCH --account=anakano_429
#SBATCH --partition=main
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=2:00:00
#SBATCH --output=%x-%j.out

# Build a pillar on a compute node (large R is too heavy for the login node). Run from the run directory:
#   mkdir -p /scratch1/$USER/tobermorite/split_R60 && cd $_
#   sbatch --job-name=R60-build /home1/mm18756/tobermorite/slurm/build_pillar.slurm --R 60 --L 240 --vac 40
# Arguments are passed to clamped/build_pillar.py. Writes pillar.data, pillar0.xyz and build.log (R0 is printed there).

module purge
module load gcc/13.3.0 python/3.11.9          # numpy from ~/.local (pip --user)
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export TMPDIR=/scratch1/$USER/tmp
mkdir -p "$TMPDIR"

set -o pipefail                                # a failed build must fail the job (not hidden by tee)
time python /home1/mm18756/tobermorite/clamped/build_pillar.py "$@" 2>&1 | tee build.log
