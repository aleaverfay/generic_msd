from generic_msd.tests.dummy_desdef import dummy_design_species
from generic_msd.tests.dummy_state_version import dummy_state_version
from generic_msd.tests.dummy_interface_job import H3H4InterfaceMSDJob
from generic_msd.tests.mergebb_design_species import MergeH3H4DesignSpecies
from generic_msd.tests.mergebb_msd_job import H3H4MergeBBInterfaceMSDJob
from generic_msd.opt_holder import OptHolder
from generic_msd.msd_interface_design import (
    MSDIntDesJobOptions,
    DesignDefinitionOpts,
    StateVersionOpts,
    PostProcessingOpts,
)
from generic_msd.msd_job_management import MSDJobManager, JobExecutionOptions
from generic_msd.server_identification import KnownComputers, ServerIdentifier
import blargs
import os


def recursively_rm_directory(dirname):
    fnames = os.listdir(dirname)
    for fname in fnames:
        fname_path = dirname + "/" + fname
        if os.path.isdir(fname_path):
            recursively_rm_directory(fname_path)
        else:
            os.remove(fname_path)
    os.rmdir(dirname)


class KilldevilServerIdentifier(ServerIdentifier):
    def what_computer(self) -> KnownComputers:
        return KnownComputers.KILLDEVIL


class DogwoodServerIdentifier(ServerIdentifier):
    def what_computer(self) -> KnownComputers:
        return KnownComputers.DOGWOOD


def test_setup_msd_job_killdevil():
    opts = OptHolder()
    currpath = os.path.dirname(os.path.abspath(__file__))
    basedir = currpath + "/dummy/"
    opts["base_dir"] = basedir

    p = blargs.Parser(opts)
    DesignDefinitionOpts.add_options(p)
    StateVersionOpts.add_options(p)
    MSDIntDesJobOptions.add_options(p)
    PostProcessingOpts.add_options(p)

    si = KilldevilServerIdentifier()
    JobExecutionOptions.add_options(p, si)

    p.process_command_line(
        [
            "--des_def",
            "desdef1_only_hphobes",
            "--state_version",
            "frwt_v1mock_dd1",
            "--job_name",
            "test_job1_killdevil",
            "--daf",
            basedir + "input_files/fitness_functions/example_func.txt",
            "--w_dGdiff_bonus_weights_file",
            basedir + "input_files/scan_values/scan_1_2.txt",
            "--entfunc_weights_file",
            basedir + "input_files/scan_values/scan_1_2_3.txt",
        ]
    )

    msd_job = H3H4InterfaceMSDJob(opts)

    if os.path.isdir("test_job1_killdevil"):
        recursively_rm_directory("test_job1_killdevil")
    msd_manager = MSDJobManager(msd_job, opts, si)
    msd_manager.prepare_job()

    # Assertions
    assert os.path.isdir("test_job1_killdevil/")
    focus_dir = "test_job1_killdevil/test_job1_killdevil_1.0w_dGdiff_1.0Ent/"
    assert os.path.isdir(focus_dir)
    for dgwt in ["1.0w", "2.0w"]:
        for entwt in ["1.0Ent", "2.0Ent", "3.0Ent"]:
            assert os.path.isdir(
                focus_dir.replace("1.0w", dgwt).replace("1.0Ent", entwt)
            )
    assert os.path.isfile(focus_dir + "fitness.daf")

    #recursively_rm_directory("test_job1_killdevil")


def test_setup_msd_job_dogwood():
    opts = OptHolder()
    currpath = os.path.dirname(os.path.abspath(__file__))
    basedir = currpath + "/dummy/"
    opts["base_dir"] = basedir

    p = blargs.Parser(opts)
    DesignDefinitionOpts.add_options(p)
    StateVersionOpts.add_options(p)
    MSDIntDesJobOptions.add_options(p)
    PostProcessingOpts.add_options(p)

    si = DogwoodServerIdentifier()
    JobExecutionOptions.add_options(p, si)

    p.process_command_line(
        [
            "--des_def",
            "desdef1_only_hphobes",
            "--state_version",
            "frwt_v1mock_dd1",
            "--job_name",
            "test_job1_dogwood",
            "--daf",
            basedir + "input_files/fitness_functions/example_func.txt",
            "--w_dGdiff_bonus_weights_file",
            basedir + "input_files/scan_values/scan_1_2.txt",
            "--entfunc_weights_file",
            basedir + "input_files/scan_values/scan_1_2_3.txt",
        ]
    )

    msd_job = H3H4InterfaceMSDJob(opts)

    if os.path.isdir("test_job1_dogwood"):
        recursively_rm_directory("test_job1_dogwood")
    msd_manager = MSDJobManager(msd_job, opts, si)
    msd_manager.prepare_job()

    # Assertions
    assert os.path.isdir("test_job1_dogwood/")
    focus_dir = "test_job1_dogwood/test_job1_dogwood_1.0w_dGdiff_1.0Ent/"
    assert os.path.isdir(focus_dir)
    for dgwt in ["1.0w", "2.0w"]:
        for entwt in ["1.0Ent", "2.0Ent", "3.0Ent"]:
            assert os.path.isdir(
                focus_dir.replace("1.0w", dgwt).replace("1.0Ent", entwt)
            )
    assert os.path.isfile(focus_dir + "fitness.daf")
    assert os.path.isfile(focus_dir + "submit.sh")

    submit_lines = None
    with open(focus_dir + "submit.sh") as fid:
        submit_lines = fid.readlines()
    # should be a slurm submission file
    assert (
        submit_lines[1] == "#SBATCH --job-name=test_job1_dogwood_1.0w_dGdiff_1.0Ent\n"
    )

    #recursively_rm_directory("test_job1_dogwood")

def test_setup_mergebb_msd_job_killdevil():
    opts = OptHolder()
    currpath = os.path.dirname(os.path.abspath(__file__))
    basedir = os.path.join(currpath,"merge_bb_inputs")
    opts["base_dir"] = basedir

    p = blargs.Parser(opts)
    DesignDefinitionOpts.add_options(p)
    StateVersionOpts.add_options(p)
    MSDIntDesJobOptions.add_options(p)
    PostProcessingOpts.add_options(p)

    si = KilldevilServerIdentifier()
    JobExecutionOptions.add_options(p, si)

    p.process_command_line(
        [
            "--des_def",
            "dd1",
            "--state_version",
            "frwt_v1mock_dd1",
            "--job_name",
            "test_job2_killdevil",
            "--daf",
            os.path.join(basedir, "input_files/fitness_functions/example_func.txt"),
            "--w_dGdiff_bonus_weights_file",
            os.path.join(basedir, "input_files/scan_values/scan_1_2.txt"),
            "--entfunc_weights_file",
            os.path.join(basedir, "input_files/scan_values/scan_1_2_3.txt"),
        ]
    )

    msd_job = H3H4MergeBBInterfaceMSDJob(opts)

    if os.path.isdir("test_job2_killdevil"):
        recursively_rm_directory("test_job2_killdevil")
    msd_manager = MSDJobManager(msd_job, opts, si)
    msd_manager.prepare_job()

    # Assertions
    assert os.path.isdir("test_job2_killdevil/")
    focus_dir = "test_job2_killdevil/test_job2_killdevil_1.0w_dGdiff_1.0Ent/"
    assert os.path.isdir(focus_dir)
    for dgwt in ["1.0w", "2.0w"]:
        for entwt in ["1.0Ent", "2.0Ent", "3.0Ent"]:
            assert os.path.isdir(
                focus_dir.replace("1.0w", dgwt).replace("1.0Ent", entwt)
            )
    assert os.path.isfile(focus_dir + "fitness.daf")

    #recursively_rm_directory("test_job1_killdevil")
