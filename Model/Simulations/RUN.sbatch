#!/bin/bash

#SBATCH --job-name=RUN     # create a short name for your job
#SBATCH --output=RUN.%a.out # stdout file
#SBATCH --error=RUN.%a.err  # stderr file
#SBATCH --nodes=1                # node count
#SBATCH --ntasks-per-node=1      # number of tasks per node

#SBATCH --requeue
#SBATCH --array=0-4

##SBATCH --export=NONE
#SBATCH --time=16:00:00

##--- modify as needed ---##
#SBATCH --cluster=c4
#SBATCH --account=gfdl_o
#SBATCH --partition=batch

##to load your elmer enviroment
source ~/elmer_env
## -----------------------##

echo "SLURM_ARRAYID="$SLURM_ARRAYID
echo "SLURM_ARRAY_JOB_ID"=$SLURM_ARRAY_JOB_ID
echo "SLURM_ARRAY_TASK_ID"=$SLURM_ARRAY_TASK_ID

##melange fraction for new damage
MF=(0.0 0.0 0.0 1.0 1.0)

##maximum rift-flank-contact fraction
CF=(1.0 1.0 1.0 0.0 0.0)

##Maximum principal vertically-averaged D_max for initialized melange section of rift near GIR
ORD=(0.995 0.995 0.98 0.86 0.995)

##Melange fraction for initialized melange section of rift near GIR
ORF=(0.0 0.4 0.4 1.0 1.0)

##maximum principal vertically-averaged damage
DAM=(0.995 0.995 0.995 0.86 0.995)

##accumulated strain in rift-opening direction at which rift flanks are no longer in contact
DSM=(0.04 0.04 0.04 0.0 0.0)

##initial distance from tip over which rift flank contact is linearly decreased from 100% to 0%
DSD=(30e3 30e3 30e3 0.0 0.0)

##damage stress threshold (MPa)
ST=(0.154 0.22 0.165 0.21 0.7)

##factor to increase restart and output frequency
ISF=(1.0 1.0 1.0 1.0 10.0)


i=$SLURM_ARRAY_TASK_ID

TARGETDIR=./M${MF[i]}_C${CF[i]}_IMD${ORD[i]}F${ORF[i]}_D${DAM[i]}_DSM${DSM[i]}D${DSD[i]}_ST${ST[i]}_
if [ ! -d "$TARGETDIR" ]; then
    mkdir $TARGETDIR
fi

cd $TARGETDIR

OUTDIR=./results
if [ ! -d "$OUTDIR" ]; then
    mkdir $OUTDIR
fi

sed "s/<mf>/${MF[i]}/g;s/<cf>/${CF[i]}/g;s/<ord>/${ORD[i]}/g;s/<orf>/${ORF[i]}/g;s/<dam>/${DAM[i]}/g;s/<dsm>/${DSM[i]}/g;s/<dsd>/${DSD[i]}/g;s/<st>/${ST[i]}/g;s/<isf>/${ISF[i]}/g" ./../RUN.sif > RUN.sif

echo RUN.sif > ELMERSOLVER_STARTINFO

timeout 15.8h srun -n1 ElmerSolver_mpi

if [[ $? -eq 124 ]]; then
  echo "job $SLURM_JOBID reached timeout, re-queueing"
  echo
  scontrol requeue $SLURM_JOBID
fi
