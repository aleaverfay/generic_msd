from enum import Enum
import socket
import os
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
) -> str:
    """Create either LSF or SLURM command for launching a job

    Create a file (possibly a job submission file)
    that holds the command(s) that should be run
    and return the command for submitting the job.

    NOTE: This command may include an input redirect operator
    ("<") so if you pass this command to another command,
    make sure to wrap it in quotes!
    """

    if job_options.scheduler == SchedulerType.LSF_SCHEDULER:
        queue = job_options.queue if job_options.queue != "auto" else "week"
        script_lines = ["#!/bin/bash\n"]
        script_lines.append("#BSUB -n %d\n" % job_options.num_nodes)
        script_lines.append("#BSUB -q %s\n" % queue)
        script_lines.append("#BSUB -o %s\n" % job_options.logfilename)
        #subcmd_list = ["bsub -n"]
        #subcmd_list.append(str(job_options.num_nodes))
        #subcmd_list.append("-q " + job_options.queue)
        #subcmd_list.append("-o " + job_options.logfilename)
        if job_options.mpi_job:
            script_lines.append("#BSUB -a mvapich\n")
            script_lines.append("mpirun %s\n" % command_line)
            #subcmd_list.append("-a mvapich mpirun " + command_line)
        else:
            script_lines.append(command_line)
            #subcmd_list.append(command_line)
        with open(job_options.submission_script_fname, "w") as fid:
            fid.writelines(script_lines)
        return "bsub < %s" % job_options.submission_script_fname

    else:
        assert job_options.scheduler == SchedulerType.SLURM_SCHEDULER


        script_lines = ["#!/bin/bash"]
        # as of 8/16, this line causes problems script_lines.append("#SBATCH --tasks-per-node=44")
        script_lines.append("#SBATCH --job-name=%s" % job_options.job_name)
        script_lines.append("#SBATCH --distribution=cyclic:cyclic")
        script_lines.append("#SBATCH --ntasks=%d" % job_options.num_nodes)
        script_lines.append("#SBATCH --output=%s" % job_options.logfilename)
        if (
            job_options.queue == "auto"
        ):  # Figure out what queue to use based on the number of nodes
            # that have been requested for this job.
            if job_options.num_nodes < 45:
                script_lines.append("#SBATCH --partition=cleanup_queue")
                job_options.timelimit = (0,4,0)
            elif job_options.num_nodes < 529:
                script_lines.append("#SBATCH --partition=528_queue")
            else:
                script_lines.append("#SBATCH --partition=2112_queue")
        elif job_options.queue == "debug" or job_options.queue == "debug_queue":
            script_lines.append("#SBATCH --partition=cleanup_queue")
            job_options.timelimit = (0,4,0)
        else:
            script_lines.append("#SBATCH --partition=%s" % job_options.queue)
        # TO DO:
        # Logic to trim the time limit for the various dogwood queues
        script_lines.append( "#SBATCH --time=%02d-%02d:%02d" % job_options.timelimit )
        script_lines.append("\n")
        if job_options.mpi_job:
            script_lines.append("$MPI_HOME/bin/mpirun  %s\n" % command_line)
        else:
            script_lines.append(command_line + "\n")
        with open(job_options.submission_script_fname, "w") as fid:
            fid.writelines("\n".join(script_lines))
        return "sbatch " + job_options.submission_script_fname

