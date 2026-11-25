#!/bin/bash
#SBATCH --job-name=tobe-pillar
#SBATCH --account=anakano_429
#SBATCH --partition=main
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2G
#SBATCH --time=8:00:00
#SBATCH --output=%x-%j.out

# ReaxFF tobermorite-11A pillar. Run from the run directory (containing pillar.data):
#   cd /scratch1/mm18756/tobermolite/run1
#   sbatch --job-name=tobe-equil    /home1/mm18756/tobermorite/slurm/run_pillar.slurm equil
#   sbatch --job-name=tobe-compress /home1/mm18756/tobermorite/slurm/run_pillar.slurm compress
# Finite pillar (clamped/build_pillar.py --L ...), stages equil_finite / compress_finite:
#   sbatch run_pillar.slurm equil_finite
#   sbatch run_pillar.slurm compress_finite -var ends clamped -var R0 <R0>   (or -var ends wall)
stage=${1:?usage: sbatch run_pillar.slurm equil|compress}

module purge
module load gcc/13.3.0 openmpi/5.0.5
LMP=/home1/$USER/tobermorite/lammps/bin/lmp   # LAMMPS 30Sep2026 stable, REAXFF+OPENMP, CPU-only
export OMP_NUM_THREADS=1

export TMPDIR=/scratch1/$USER/tmp
mkdir -p "$TMPDIR"

# in.damp_finite overwrites thermal.restart: keep the undamped one
[[ $stage == damp* ]] && cp thermal.restart thermal_nvt.restart

srun $LMP -in /home1/mm18756/tobermorite/inputs/in.$stage -log log.$stage "${@:2}"   # extra args, e.g. -var nnpt 400

# backed-up copy of the small outputs; trajectories stay on scratch until reviewed
dest=/home1/$USER/tobermorite/results/$(basename "$PWD")   # per run directory
mkdir -p "$dest"
cp log.$stage "$dest"/
[[ $stage == equil* ]]    && cp relaxed.data "$dest"/
[[ $stage == compress* ]] && cp species.out compressed.data "$dest"/
exit 0
