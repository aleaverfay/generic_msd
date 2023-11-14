#!/bin/bash
#SBATCH --job-name=test_run5_dock_jobs_view
#SBATCH --distribution=cyclic:cyclic
#SBATCH --ntasks=1
#SBATCH --output=test_run5_djv.log
#SBATCH --partition=cleanup_queue
#SBATCH --time=00-04:00


python3 /nas/longleaf/home/leaverfa/ortho_nuc/h3h4/twomutpairs_merge_RINC/pyscripts/dock_jobs_view.py -l complex_sets.list -o after_docking_dGbind.txt --docking_flags_files /nas/longleaf/home/leaverfa/ortho_nuc/h3h4/twomutpairs_merge_RINC/input_files/post_processing/options_full_relax.flags  

