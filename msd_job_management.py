import blargs
from enum import Enum
import os
from server_identification import ServerIdentifier
from server_identification import KnownComputers as comp
import typing
import math


def leading_zero_string(num: int, maxrange: int) -> str:
    return ("%%0%sd" % (int(math.floor(math.log10(maxrange))),)) % (num,)


def mkdir(dirname):
    if not os.path.exists(dirname):
        os.makedirs(dirname)
        return True
    return False


class SchedulerType(Enum):
    LSF_SCHEDULER = 1
    SLURM_SCHEDULER = 2
    DEFAULT = LSF_SCHEDULER


def scheduler_type_for_server(si):
    server = si.what_computer()
    assert server == comp.KILLDEVIL or server == comp.DOGWOOD
    return (
        SchedulerType.LSF_SCHEDULER
        if server == comp.KILLDEVIL
        else SchedulerType.SLURM_SCHEDULER
    )


class JobExecutionOptions:
    def __init__(self, opts, si: ServerIdentifier):
        self.num_cpu = opts.num_cpu
        self.num_states_per_cpu = opts.num_states_per_cpu
        self.launch = hasattr(opts, "launch")
        self.queue = opts.queue
        self.si = si

    @staticmethod
    def add_options(blargs_parser: blargs.Parser, si: ServerIdentifier):
        blargs_parser.int("num_cpu").shorthand("n").default(-1)
        blargs_parser.int("num_states_per_cpu").default(1)
        blargs_parser.flag("launch")
        server = si.what_computer()
        if server == comp.KILLDEVIL:
            blargs_parser.str("queue").shorthand("q").default("week")
        elif server == comp.DOGWOOD:
            blargs_parser.str("queue").shorthand("q").default("auto")


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
        subcmd_list.append("-q " + job_options.queue)
        subcmd_list.append("-o " + job_options.logfilename)
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


def rosetta_dir(si: ServerIdentifier):
    server = si.what_computer()
    assert server == comp.KILLDEVIL or server == comp.DOGWOOD
    return (
        "/nas02/home/l/e/leaverfa/GIT/main/"
        if server == comp.KILLDEVIL
        else "/nas/longleaf/home/leaverfa/GIT/main/"
    )


def bin_dir(si: ServerIdentifier):
    return rosetta_dir(si) + "source/bin/"


def db_path(si: ServerIdentifier):
    return rosetta_dir(si) + "database/"


def mpi_msd_exe(si: ServerIdentifier):
    server = si.what_computer()
    assert server == comp.KILLDEVIL or server == comp.DOGWOOD

    return (
        bin_dir(si) + "mpi_msd.mpi.linuxgccrelease"
        if server == comp.KILLDEVIL
        else bin_dir(si) + "mpi_msd.mpiserialization.linuxgccrelease"
    )


def pyscripts_path(si: ServerIdentifier):
    return (
        "/nas02/home/l/e/leaverfa/pyscripts"
        if si.what_computer() == comp.KILLDEVIL
        else "/nas/longleaf/home/leaverfa/pyscripts"
    )


class MSDJobManager:
    def __init__(self, msd_job, opts, si):
        self.options = JobExecutionOptions(opts, si)
        self.msd_job = msd_job
        self.launch_script = []
        self.gather_script = []

        self.gather_script.append("mkdir results\n")
        self.gather_script.append("\n")
        self.si = si
        self.base_dir = opts.base_dir

    @property
    def job_name(self):
        return self.msd_job.job_name()

    def prepare_job(self):
        """Main entry method for creating (and perhaps launching) a set of MSD jobs"""

        if not mkdir(self.job_name):
            raise ValueError(
                "ERROR: job directory", self.msd_job.job_name(), "already exists"
            )

        os.chdir(self.job_name)
        for subjob in self.msd_job.subjobs():
            if not mkdir(subjob):
                raise ValueError("Error creating subdirectory", subjob)
            os.chdir(subjob)
            self.create_symlinks(subjob)
            self.write_fitness_file(subjob)
            self.create_submission_command(subjob)
            os.chdir("..")

        self.prepare_postprocessing_step()
        self.write_submission_and_gather_files()
        if self.options.launch:
            self.launch()
        os.chdir("..")

    def create_symlinks(self, subdir):
        filename_pairs = self.msd_job.files_to_symlink(subdir)
        for src_fname, dest_fname in filename_pairs:
            if not os.path.isfile(src_fname):
                raise ValueError(
                    "Could not locate file "
                    + src_fname
                    + " when creating symlinks for sub-job "
                    + subdir
                    + " for job "
                    + self.jobname
                    + "."
                )
            else:
                os.system(" ".join(("ln -s", src_fname, dest_fname)))

    def write_fitness_file(self, subdir):
        fitness_lines = self.msd_job.fitness_lines(subdir)
        with open("fitness.daf", "w") as fid:
            fid.writelines(fitness_lines)

    def create_submission_command(self, subdir):
        nprocs = self.nprocs_for_job()
        queue_name = self.options.queue
        popsize = self.msd_job.popsize()
        ngen = self.msd_job.ngen()

        job_options = SubmissionOptions()
        job_options.scheduler = scheduler_type_for_server(self.si)
        job_options.num_nodes = nprocs
        job_options.queue = queue_name
        job_options.job_name = subdir
        job_options.logfilename = subdir + ".log"
        job_options.submission_script_fname = "submit.sh"
        job_options.mpi_job = True

        command_line = (
            mpi_msd_exe(self.si)
            + " -database "
            + db_path(self.si)
            + " -entity_resfile entity.resfile "
            + "-fitness_file fitness.daf -ms::pop_size "
            + str(popsize)
            + " -ms::generations "
            + str(ngen)
            + " -ms::numresults "
            + str(self.msd_job.ntop_msd_results_to_dock)
        )
        for flagsfile in self.msd_job.msd_flags(subdir):
            command_line += " @" + flagsfile
        if self.msd_job.seeded(subdir):
            command_line += " -seed_sequences " + " ".join(self.msd_job.seeds())
        if self.msd_job.fill_gen1_from_seeds():
            command_line += " -fill_gen1_from_seed_sequences"

        submission_command, command_filename = command_and_submission_script_for_job(
            job_options, command_line
        )
        self.launch_script.append(submission_command + " >> ../msd_submission.log\n")

        self.gather_script.append("cd " + subdir + "\n")
        species = self.msd_job.states_to_save()
        complexes = self.msd_job.complexes_to_postprocess()
        n_to_postprocess = self.msd_job.n_results_to_postprocess()
        for spec in species:
            for result_ind in range(1, n_to_postprocess + 1):
                lzri = leading_zero_string(result_ind, n_to_postprocess)
                self.gather_script.append(
                    "for i in `ls msd_output_%d_%s*`; do if [ ! -h $i ];"
                    " then cp $i ../results/%s_%s_%s.pdb; fi; done\n"
                    % (result_ind, spec, subdir, lzri, spec)
                )
        self.gather_script.append("cd ..\n")

        for result_ind in range(1, n_to_postprocess + 1):
            lzri = leading_zero_string(result_ind, n_to_postprocess)
            prefix = subdir + "_" + lzri + "_"
            self.gather_script.append(
                "echo "
                + prefix
                + (".pdb" + prefix).join(complexes)
                + ".pdb >> results/complex_sets.list\n"
            )
        self.gather_script.append("\n")

    def nprocs_for_job(self):
        if self.options.num_cpu == -1:
            nstates = self.msd_job.state_version.nstates_total()
            return int(math.ceil(nstates / self.options.num_states_per_cpu))
        else:
            return self.options.num_cpu

    def prepare_postprocessing_step(self):
        launch_docking_script = []
        launch_docking_script.append("bash gather_output.sh\n")
        launch_docking_script.append("cd results\n")

        # part of me really wants to have the logic for this docking script wrapped
        # up in the same envelope as the logic for this script. It's weird for the
        # options for the second script to "bleed through" to the options for this script
        # and I'd like to avoid duplicating any "what's the base_dir" logic.
        dock_jobs_command = [
            "python3",
            self.base_dir + "pyscripts/dock_jobs_run.py",
            "--pdb-triples",
            "complex_sets.list",
            self.msd_job.post_processing_opts.to_command_line(),
            "\n",
        ]

        launch_docking_script.append(" ".join(dock_jobs_command))

        # now prepare the dock_jobs_view command and append it to the launch_docking script
        dock_jobs_view_job_opts = SubmissionOptions()
        dock_jobs_view_job_opts.scheduler = scheduler_type_for_server(self.si)
        dock_jobs_view_job_opts.num_nodes = 1
        dock_jobs_view_job_opts.queue = "debug_queue"
        dock_jobs_view_job_opts.job_name = self.msd_job.job_name() + "_dock_jobs_view"
        dock_jobs_view_job_opts.logfilename = self.msd_job.job_name() + "_djv.log"
        dock_jobs_view_job_opts.submission_script_fname = "djv_submit.sh"

        dock_jobs_view_command_line = (
            "python3 "
            + self.base_dir
            + "pyscripts/dock_jobs_view.py -l complex_sets.list -o after_docking_dGbind.txt"
        )
        dock_jobs_view_submission_command, dock_jobs_view_command_filename = command_and_submission_script_for_job(
            dock_jobs_view_job_opts, dock_jobs_view_command_line
        )

        if dock_jobs_view_command_filename is None:
            with open("djv_submit.sh", "w") as fid:
                fid.writelines(dock_jobs_view_submission_command)
        dock_jobs_view_command = [
            "python3",
            pyscripts_path(self.si) + "submit_dependent_script.py",
            "dock/dock_submission.log djv_submit.sh",
            self.msd_job.post_processing_opts.to_command_line(),
            "\n",
        ]

        # if options.relax :
        #    dock_jobs_view_command += " --relax"
        #    # if options.relax_protocol != "" :
        #    #     dock_jobs_view_command += " --relax-protocol " + options.relax_protocol

        launch_docking_script.append(" ".join(dock_jobs_view_command))

        launch_docking_job_opts = SubmissionOptions()
        launch_docking_job_opts.scheduler = scheduler_type_for_server(self.si)
        launch_docking_job_opts.num_nodes = 1
        launch_docking_job_opts.queue = "debug_queue"
        launch_docking_job_opts.job_name = self.msd_job.job_name() + "_launch_docking"
        launch_docking_job_opts.logfilename = (
            self.msd_job.job_name() + "_launch_docking.log"
        )
        launch_docking_job_opts.submission_script_fname = "launch_docking.sh"

        launch_docking_submission_command, launch_docking_command_filename = command_and_submission_script_for_job(
            launch_docking_job_opts, "".join(launch_docking_script)
        )

        if launch_docking_command_filename is None:
            open("launch_docking.sh", "w").writelines(launch_docking_submission_command)
        open("prepare_for_docking.sh", "w").writelines(
            " ".join(
                [
                    "python3",
                    pyscripts_path(self.si) + "submit_dependent_script.py",
                    "msd_submission.log launch_docking.sh\n",
                ]
            )
        )

    def write_submission_and_gather_files(self):
        with open("submit_all_jobs.sh", "w") as fid:
            fid.writelines(self.launch_script)
        with open("gather_output.sh", "w") as fid:
            fid.writelines(self.gather_script)

    def launch(self):
        pass
