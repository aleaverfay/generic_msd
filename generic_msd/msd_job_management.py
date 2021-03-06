import blargs
from enum import Enum
import os
from .server_identification import ServerIdentifier, KnownComputers
from .create_lsf_or_slurm_job import (SubmissionOptions, command_and_submission_script_for_job, SchedulerType, scheduler_type_for_server)
from .msd_interface_design import add_required_opts_to_blargs_parser
import typing
import math
import sys


def leading_zero_string(num: int, maxrange: int) -> str:
    return ("%%0%sd" % (int(math.floor(math.log10(maxrange))),)) % (num,)


def mkdir(dirname):
    if not os.path.exists(dirname):
        os.makedirs(dirname)
        return True
    return False


class JobExecutionOptions:
    def __init__(self, opts, si: ServerIdentifier):
        self.num_cpu = opts.num_cpu
        self.num_states_per_cpu = opts.num_states_per_cpu
        self.launch = opts.launch
        self.queue = opts.queue
        self.si = si

    @classmethod
    def add_options(cls, blargs_parser: blargs.Parser, si: ServerIdentifier):
        cls.add_required_options(blargs_parser, si)
        cls.add_non_required_options(blargs_parser, si)

    @classmethod
    def add_required_options(cls, blargs_parser: blargs.Parser, si: ServerIdentifier):
        add_required_opts_to_blargs_parser(blargs_parser, cls.required_opts())

    @classmethod
    def required_opts(cls):
        return []

    @classmethod
    def add_non_required_options(cls, blargs_parser: blargs.Parser, si: ServerIdentifier):
        blargs_parser.int("num_cpu").shorthand("n").default(-1)
        blargs_parser.int("num_states_per_cpu").default(1)
        blargs_parser.flag("launch")
        server = si.what_computer()
        if server == KnownComputers.KILLDEVIL:
            blargs_parser.str("queue").shorthand("q").default("week")
        elif server == KnownComputers.DOGWOOD:
            blargs_parser.str("queue").shorthand("q").default("auto")
        else:
            # for debugging on other machines
            blargs_parser.str("queue").shorthand("q").default("auto")

    def to_command_line(self):
        args = []
        if self.num_cpu != -1 :
            args.append("--num_cpu")
            args.append(str(self.num_cpu))
        if self.num_states_per_cpu != 1:
            args.append("--num_states_per_cpu")
            args.append(str(self.num_states_per_cpu))
        if self.launch:
            args.append("--launch")
        return " ".join(args)


def rosetta_dir(si: ServerIdentifier):
    server = si.what_computer()
    assert (
        server == KnownComputers.KILLDEVIL or
        server == KnownComputers.DOGWOOD or
        server == KnownComputers.LONGLEAF
    )
        
    return (
        "/nas02/home/l/e/leaverfa/GIT/main/"
        if server == KnownComputers.KILLDEVIL
        else "/nas/longleaf/home/leaverfa/GIT/main/"
    )


def bin_dir(si: ServerIdentifier):
    return rosetta_dir(si) + "source/bin/"


def db_path(si: ServerIdentifier):
    return rosetta_dir(si) + "database/"


def mpi_msd_exe(si: ServerIdentifier):
    server = si.what_computer()
    assert (
        server == KnownComputers.KILLDEVIL or
        server == KnownComputers.DOGWOOD or
        server == KnownComputers.LONGLEAF
    )

    return (
        bin_dir(si) + "mpi_msd.mpiserialization.linuxgccrelease"
        if server == KnownComputers.KILLDEVIL
        else bin_dir(si) + "mpi_msd.mpiserialization.linuxgccrelease"
    )

def rosetta_scripts_jd3_exe(si: ServerIdentifier):
    #server = si.what_computer()
    #assert server == KnownComputers.KILLDEVIL or server == KnownComputers.DOGWOOD
    return bin_dir(si) + "rosetta_scripts_jd3.mpiserialization.linuxgccrelease"



def pyscripts_path(si: ServerIdentifier):
    return (
        "/nas02/home/l/e/leaverfa/pyscripts/"
        if si.what_computer() == KnownComputers.KILLDEVIL
        else "/nas/longleaf/home/leaverfa/pyscripts/"
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
        self.save_creation_command()
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
        job_options.server = self.si.what_computer()
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

        submission_command = command_and_submission_script_for_job(
            job_options, command_line
        )
        self.launch_script.append("cd " + subdir + "\n")
        self.launch_script.append(submission_command + " >> ../msd_submission.log\n")
        self.launch_script.append("cd ..\n")
        self.launch_script.append("\n")

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
                + (".pdb " + prefix).join(complexes)
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
            os.path.join(self.base_dir, "pyscripts/dock_jobs_run.py"),
            "--pdb-complexes",
            "complex_sets.list",
            self.msd_job.post_processing_opts.to_command_line(),
            self.msd_job.extra_post_processing_options(),
            "\n",
        ]

        launch_docking_script.append(" ".join(dock_jobs_command))
        # we need the djv_submit script to live in the results directory, but that directory
        # doesn't exist yet and won't exist until after the gather_output script runs
        launch_docking_script.append("ln -s ../djv_submit.sh\n")

        # now prepare the dock_jobs_view command and append it to the launch_docking script
        dock_jobs_view_job_opts = SubmissionOptions()
        dock_jobs_view_job_opts.server = self.si.what_computer()
        dock_jobs_view_job_opts.scheduler = scheduler_type_for_server(self.si)
        dock_jobs_view_job_opts.num_nodes = 1
        dock_jobs_view_job_opts.queue = "debug_queue"
        dock_jobs_view_job_opts.job_name = self.msd_job.job_name() + "_dock_jobs_view"
        dock_jobs_view_job_opts.logfilename = self.msd_job.job_name() + "_djv.log"
        dock_jobs_view_job_opts.submission_script_fname = "djv_submit.sh"

        dock_jobs_view_command_line = " ".join([
            "python3",
            os.path.join(self.base_dir, "pyscripts/dock_jobs_view.py"),
            "-l complex_sets.list -o after_docking_dGbind.txt",
            self.msd_job.post_processing_opts.to_command_line(),
            self.msd_job.extra_post_processing_options(),
            "\n"
            ])
        command_and_submission_script_for_job(
           dock_jobs_view_job_opts, dock_jobs_view_command_line
        )

        dock_jobs_view_submission_command = [
            "python3",
            os.path.join(pyscripts_path(self.si), "submit_dependent_script.py"),
            "dock/dock_submission.log",
            dock_jobs_view_job_opts.submission_script_fname,
            "\n",
        ]

        # if options.relax :
        #    dock_jobs_view_command += " --relax"
        #    # if options.relax_protocol != "" :
        #    #     dock_jobs_view_command += " --relax-protocol " + options.relax_protocol

        launch_docking_script.append(" ".join(dock_jobs_view_submission_command))

        launch_docking_job_opts = SubmissionOptions()
        launch_docking_job_opts.server = self.si.what_computer()
        launch_docking_job_opts.scheduler = scheduler_type_for_server(self.si)
        launch_docking_job_opts.num_nodes = 1
        launch_docking_job_opts.queue = "debug"
        launch_docking_job_opts.job_name = self.msd_job.job_name() + "_launch_docking"
        launch_docking_job_opts.logfilename = (
            self.msd_job.job_name() + "_launch_docking.log"
        )
        launch_docking_job_opts.submission_script_fname = "launch_docking.sh"
        
        command_and_submission_script_for_job(
            launch_docking_job_opts, "".join(launch_docking_script)
        )
            
        with open("prepare_for_docking.sh", "w") as fid:
            fid.writelines(
                " ".join(
                    [
                        "python3",
                        os.path.join(pyscripts_path(self.si), "submit_dependent_script.py"),
                        "msd_submission.log launch_docking.sh\n",
                    ]
                )
            )

    def write_submission_and_gather_files(self):
        with open("submit_all_jobs.sh", "w") as fid:
            fid.writelines(self.launch_script)
        with open("gather_output.sh", "w") as fid:
            fid.writelines(self.gather_script)

    def save_creation_command(self):
        with open("creation_command.txt", "w") as fid:
            fid.writelines( " ".join(sys.argv) + "\n")
        
    def launch(self):
        os.system("bash submit_all_jobs.sh")
        os.system("bash prepare_for_docking.sh")
