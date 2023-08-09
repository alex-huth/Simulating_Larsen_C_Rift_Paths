The following slurm scripts will run the experiments from Figure 8. Modify these scripts where indicated for other computing systems.

Each script runs a job array of experients, where each experiment differs in its rift-flank conditions and damage parameters.

First, run INIT.sbatch (`sbatch INIT.sbatch`). This will read in the data and save the initial state to a restart file for each experiment.

Then, run the forward simulations from these initial states:  `sbatch RUN.sbatch`.

As the forward simulations run, their updated states are saved to their respective restart files at least every 0.1 simulation years (by default). If a job times out, it will requeue and the simulation corresponding to the job will restart from the last state saved in its restart file.

