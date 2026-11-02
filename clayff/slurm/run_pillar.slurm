#!/bin/bash
#SBATCH --job-name=clayff-pillar
#SBATCH --account=anakano_81
#SBATCH --partition=main
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2G
#SBATCH --time=1:00:00
#SBATCH --output=%x-%j.out

# ClayFF tobermorite-11A pillar. Run from the run directory (containing pillar_clayff.data):
#   cd /scratch1/anakano/tobermolite/clayff/run1
#   sbatch --job-name=clayff-equil    /home1/anakano/tobermolite/clayff/slurm/run_pillar.slurm equil
#   sbatch --job-name=clayff-compress /home1/anakano/tobermolite/clayff/slurm/run_pillar.slurm compress
stage=${1:?usage: sbatch run_pillar.slurm equil|compress}

module purge
module load gcc/13.3.0 openmpi/5.0.5
LMP=/home1/$USER/tobermolite/clayff/lammps/bin/lmp   # LAMMPS 30Sep2026, MOLECULE+KSPACE+OPENMP
export OMP_NUM_THREADS=1

export TMPDIR=/scratch1/$USER/tmp
mkdir -p "$TMPDIR"

srun $LMP -in /home1/anakano/tobermolite/clayff/inputs/in.$stage -log log.$stage

# backed-up copy of the small outputs; trajectories stay on scratch until reviewed
dest=/home1/$USER/tobermolite/clayff/results/run1
mkdir -p "$dest"
cp log.$stage "$dest"/
[ "$stage" = equil ]    && cp relaxed.data "$dest"/
[ "$stage" = compress ] && cp compressed.data "$dest"/
exit 0
