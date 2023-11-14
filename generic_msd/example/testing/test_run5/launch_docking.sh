#!/bin/bash
#SBATCH --job-name=test_run5_launch_docking
#SBATCH --distribution=cyclic:cyclic
#SBATCH --ntasks=1
#SBATCH --output=test_run5_launch_docking.log
#SBATCH --partition=cleanup_queue
#SBATCH --time=00-04:00


bash gather_output.sh
cd results
python3 /nas/longleaf/home/leaverfa/ortho_nuc/h3h4/twomutpairs_merge_RINC/pyscripts/dock_jobs_run.py --pdb-complexes complex_sets.list --docking_flags_files /nas/longleaf/home/leaverfa/ortho_nuc/h3h4/twomutpairs_merge_RINC/input_files/post_processing/options_full_relax.flags  
ln -s ../djv_submit.sh
python3 /nas/longleaf/home/leaverfa/pyscripts/submit_dependent_script.py dock/dock_submission.log djv_submit.sh 

