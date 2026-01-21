#!/bin/bash
##ENVIRONMENT SETTINGS; CHANGE WITH CAUTION
#SBATCH --export=NONE                #Do not propagate environment
#SBATCH --get-user-env=L             #Replicate login environment

##NECESSARY JOB SPECIFICATIONS
#SBATCH --job-name=ESI_Final_Prop
#SBATCH --time=30:00:00
#SBATCH --nodes=1                   
#SBATCH --ntasks-per-node=40  
#SBATCH --mem=300G
#SBATCH --output=OUT_tm_model.%j
#SBATCH --account=132772335884

##OPTIONAL JOB SPECIFICATIONS
#SBATCH --mail-type=ALL 
#SBATCH --mail-user=brentvela@tamu.edu



module load GCC/9.3.0
module load Python/3.8.2
module load OpenBLAS/0.3.9
module load Thermo-Calc/2023.1
cd $SCRATCH
source /scratch/user/brentvela/venv/tcpython-2023.1/bin/activate
cd $SLURM_SUBMIT_DIR
python et_melt_pool_script.py
