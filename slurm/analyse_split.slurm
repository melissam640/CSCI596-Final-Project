#!/bin/bash
#SBATCH --job-name=split-analyse
#SBATCH --account=anakano_429
#SBATCH --partition=main
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=3
#SBATCH --mem=48G
#SBATCH --time=4:00:00
#SBATCH --output=%x-%j.out

# SPLIT_PLAN step e on a compute node (225k-atom dumps are too heavy for the login node).
# Run from the run directory, e.g.
#   cd /scratch1/$USER/tobermorite/split_R60
#   sbatch --job-name=R60-analyse /home1/mm18756/tobermorite/slurm/analyse_split.slurm split_R60_clamped
# Arg 1 = label (also names results/stress_strain_<label>.dat, which damage_map.py overlays).
# Optional env: STRAINS_S (stress maps, %), STRAINS_D (damage maps, %), CORE (damage core radius, A).
# Runs stress_map.py, damage_map.py and gauge_strain.py in parallel (one core each).

LABEL=${1:?label}
STRAINS_S=${STRAINS_S:-5,10,14,16,17,18,20,25}
STRAINS_D=${STRAINS_D:-10,14,16,17,18,20,25}
CORE=${CORE:-30}
REPO=/home1/mm18756/tobermorite
RUN=${LABEL%_clamped}; RUN=${RUN%_wall}

module purge
module load gcc/13.3.0 python/3.11.9          # numpy/matplotlib from ~/.local (pip --user)
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MPLBACKEND=Agg
export TMPDIR=/scratch1/$USER/tmp
mkdir -p "$TMPDIR" "$REPO/results/stress_$RUN" "$REPO/results/damage_$RUN"

set -o pipefail
( time python $REPO/scripts/stress_map.py $REPO/results/stress_$RUN $LABEL=stress_compress.lammpstrj \
      --strains $STRAINS_S ) > stress_map.log 2>&1 &
( time python $REPO/scripts/damage_map.py $REPO/results/damage_$RUN $LABEL=traj_compress.lammpstrj \
      --strains $STRAINS_D --strain-per-frame 0.001 --core $CORE ) > damage_map.log 2>&1 &
( time python $REPO/scripts/gauge_strain.py traj_compress.lammpstrj $REPO/results/gauge_$LABEL.dat ) \
      > gauge.log 2>&1 &
fail=0
for p in $(jobs -p); do wait $p || fail=1; done
tail -n 30 stress_map.log damage_map.log gauge.log
exit $fail
