#!/bin/bash
##ENVIRONMENT SETTINGS; CHANGE WITH CAUTION
#SBATCH --export=NONE                #Do not propagate environment
#SBATCH --get-user-env=L             #Replicate login environment

##NECESSARY JOB SPECIFICATIONS
#SBATCH --job-name=et_melt_pool
#SBATCH --time=30:00:00
#SBATCH --nodes=1                   
#SBATCH --ntasks-per-node=40  
#SBATCH --mem=300G
#SBATCH --output=OUT_tm_model.%j
#SBATCH --account=YOUR_ACCOUNT

##OPTIONAL JOB SPECIFICATIONS
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=YOUR_EMAIL



module load GCC/9.3.0
module load Python/3.8.2
module load OpenBLAS/0.3.9
module load Thermo-Calc/2023.1
cd ${SCRATCH:-$SLURM_SUBMIT_DIR}
source ${VENV_PATH:-$HOME/venv/bin/activate}
cd $SLURM_SUBMIT_DIR
python et_melt_pool_script.py
