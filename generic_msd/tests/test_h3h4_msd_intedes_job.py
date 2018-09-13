from msd_interface_design import (
    MSDIntDesJobOptions,
    DesignDefinitionOpts,
    StateVersionOpts,
    PostProcessingOpts,
)
from tests.dummy_desdef import dummy_design_species
from tests.dummy_state_version import dummy_state_version
from tests.dummy_interface_job import H3H4InterfaceMSDJob
from opt_holder import OptHolder
import blargs
import os


def test_create_h3h4_msd_job():
    opts = OptHolder()
    currpath = os.path.dirname(os.path.abspath(__file__))
    basedir = currpath + "/dummy/"
    opts["base_dir"] = basedir

    p = blargs.Parser(opts)
    DesignDefinitionOpts.add_options(p)
    StateVersionOpts.add_options(p)
    MSDIntDesJobOptions.add_options(p)
    PostProcessingOpts.add_options(p)

    p.process_command_line(
        [
            "--des_def",
            "desdef1_only_hphobes",
            "--state_version",
            "frwt_v1mock_dd1",
            "--job_name",
            "test_job1",
            "--daf",
            basedir + "input_files/fitness_functions/example_func.txt",
        ]
    )

    msd_job = H3H4InterfaceMSDJob(opts)
    linkable_files = msd_job.files_to_symlink(None)

    for p in linkable_files:
        assert os.path.isfile(p[0])


def test_create_h3h4_fitness_lines():
    opts = OptHolder()
    currpath = os.path.dirname(os.path.abspath(__file__))
    basedir = currpath + "/dummy/"
    opts["base_dir"] = basedir

    p = blargs.Parser(opts)
    DesignDefinitionOpts.add_options(p)
    StateVersionOpts.add_options(p)
    MSDIntDesJobOptions.add_options(p)
    PostProcessingOpts.add_options(p)

    p.process_command_line(
        [
            "--des_def",
            "desdef1_only_hphobes",
            "--state_version",
            "frwt_v1mock_dd1",
            "--job_name",
            "test_job1",
            "--daf",
            basedir + "input_files/fitness_functions/example_func.txt",
        ]
    )

    msd_job = H3H4InterfaceMSDJob(opts)
    msd_job.files_to_symlink(None)
    hypothetical_job = "testjob_1.5w_dGdiff_3.25Ent"
    fitness_func_lines = msd_job.fitness_lines(hypothetical_job)

    assert (
        fitness_func_lines[-1]
        == "FITNESS best_MH3_MH4 + 1.500000 * best_dGbind + 3.250000 * entfunc\n"
    )
