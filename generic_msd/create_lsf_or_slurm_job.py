from enum import Enum
import socket
import os
import typing
from .server_identification import ServerIdentifier, KnownComputers

class SchedulerType(Enum):
    LSF_SCHEDULER = 1
    SLURM_SCHEDULER = 2
    DEFAULT = LSF_SCHEDULER


def scheduler_type_for_server(si: ServerIdentifier):
    server = si.what_computer()
    assert server == KnownComputers.KILLDEVIL or server == KnownComputers.DOGWOOD
    return (
        SchedulerType.LSF_SCHEDULER
        if server == KnownComputers.KILLDEVIL
        else SchedulerType.SLURM_SCHEDULER
    )

class SubmissionOptions:
    def __init__(self):
        self.scheduler = SchedulerType.LSF_SCHEDULER
        self.num_nodes = 0
        self.queue = "auto"
        self.job_name = ""
        self.logfilename = ""
        self.timelimit = (1, 0, 0)
        self.submission_script_fname = ""
        self.mpi_job = False


def command_and_submission_script_for_job(
    job_options: SubmissionOptions, command_line: str
) -> typing.Tuple[str, str]:
    """Create either LSF or SLURM command for launching a job

    If the scheduler is SLURM, then create a job submission file
    that holds the command for launching the input job
    and return the command for launching
    """

    if job_options.scheduler == SchedulerType.LSF_SCHEDULER:
        subcmd_list = ["bsub -n "]
        subcmd_list.append(str(job_options.num_nodes))
        subcmd_list.append(" -q " + job_options.queue)
        subcmd_list.append(" -o " + job_options.logfilename)
        if job_options.mpi_job:
            subcmd_list.append(" -a mvapich mpirun " + command_line)
        else:
            subcmd_list.append(command_line)
        return ("".join(subcmd_list), None)
    else:
        assert job_options.scheduler == SchedulerType.SLURM_SCHEDULER

        subcmd = "sbatch " + job_options.submission_script_fname
        script_lines = ["#!/bin/bash"]
        script_lines.append("#SBATCH --job-name=%s" % job_options.job_name)
        script_lines.append("#SBATCH --distribution=cyclic")
        script_lines.append("#SBATCH --ntasks=%d" % job_options.num_nodes)
        script_lines.append("#SBATCH --output=%s" % job_options.logfilename)
        if (
            job_options.queue == "auto"
        ):  # Figure out what queue to use based on the number of nodes
            # that have been requested for this job.
            if job_options.num_nodes < 45:
                script_lines.append("#SBATCH --partition=debug_queue")
            elif job_options.num_nodes < 529:
                script_lines.append("#SBATCH --partition=528_queue")
            else:
                script_lines.append("#SBATCH --partition=2112_queue")
        else:
            script_lines.append("#SBATCH --partition=%s" % job_options.queue)
        # script_lines.append( "#SBATCH --time=%02d:%02d:%02d" % job_options.timelimit )
        script_lines.append("\n")
        if job_options.mpi_job:
            script_lines.append("srun -n $SLURM_NPROCS --mpi=pmi2 %s\n" % command_line)
        else:
            script_lines.append(command_line + "\n")
        with open(job_options.submission_script_fname, "w") as fid:
            fid.writelines("\n".join(script_lines))
        return (subcmd, job_options.submission_script_fname)