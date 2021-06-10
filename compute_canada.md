# Getting started with Compute Canada

1. Make an account at https://ccdb.computecanada.ca; To use some clusters like Niagara, you need to specially request access from the account.

2. Connect to the specific cluster using SSH in the terminal of your local machine. For example, if you want to access Niagara, you would do

```
ssh -Y username@niagara.computecanada.ca	
```
	
The username and password (for which you will be prompted) are those of your Compute Canada account.

3. You will probably want to copy files from your machine to your cluster filesystem. You do this using

```
scp pathToFile username@niagara.computecanada.ca:destinationPath
```
	
Here `pathToFile` is the path to the file you want to copy and `destinationPath` is where you want the file to be copied to. Copying files from the cluster to your machine is similar:

```
scp username@niagara.computecanada.ca:filePath destinationPath
```
	
4. Usually, you will have to copy mainly three files to the cluster: the job script, your uncompiled C++ code, and a makefile (not strictly neccessary).

5. In compiling your code, you will want to include a special flag that will allow you to take advantage of the cluster's optimized CPUs. If you are using gcc to compile, including the flag `-march=native` should be enough (don't forget the `-O3` flag also, which is useful even on your local machine).

6. A simple job script (.sh file) looks something like this:

```
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=40
#SBATCH --time=24:00:00
#SBATCH --job-name a5000_20000000_X8=1
#SBATCH --mail-user= email@gmail.com
#SBATCH --mail-type=ALL
cd $SLURM_SUBMIT_DIR
module load NiaEnv/2018a
module load gcc/7.3.0
module load boost/1.66.0
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
make
./PRecCppCl_stack 5000 20000000 $SLURM_CPUS_PER_TASK 1
```

This is a sample Niagara job script using 40 CPUs on a single node requesting 24 hours of wall-time. It arranges automated emails to notify you when the job begins, ends, or is cancelled/fails for some reason. It uses `OpenMp` for multi-threading utilizing 40 threads. The `module load ...` lines load the compiler and the Boost library. You might want to use more recent versions of these libraries. You can check which versions are available like this:

```
module spider gcc
```
The `make` invokes the makefile -- alternatively, you can instead just put compilation instruction directly there.

7. You can submit jobs using
```
sbatch jobScriptFile.sh
```
When you submit a job, you receive a job id, which is useful to know more about your job status.

8. The following commands are useful to obtain info about your job:

```
squeue -u username        # (about your queued jobs)

squeue -j JOBID    	  # (status of job with JOBID)

scancel -i JOBID          # (cancel a job)

sacct      		  # (info about recent jobs, even running/finished ones)

jobperf JOBID	          # (gives you performance stats for the job, once it's running)
```

9. If possible, it's useful to save intermediate data -- the cluster nodes fail surprisingly often.

10. More useful info: https://docs.computecanada.ca/wiki/Niagara_Quickstart (see sidebar for other clusters).
