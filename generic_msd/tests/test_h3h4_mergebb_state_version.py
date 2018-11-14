from generic_msd.msd_interface_design import (
    DesignSpecies,
    DesignDefinitionOpts,
    MergeBBDesDefFnames,
    MergeBBStateVersion,
    StateVersionOpts,
)
from generic_msd.tests.mergebb_design_species import MergeH3H4DesignSpecies

import os


class mock_opts:
    pass

def test_create_state_version():
    currpath = os.path.dirname(os.path.abspath(__file__))
    basedir = os.path.join(currpath, "merge_bb_inputs")

    opts = mock_opts()
    setattr(opts, "base_dir", basedir)
    setattr(opts, "state_version", "frwt_v1mock_dd1")
    state_ver_opts = StateVersionOpts(opts)
    stateverdir = os.path.join(basedir,"input_files","state_versions","frwt_v1mock_dd1")

    # let's delete all of the .states files in this directory
    # if they are already present so that we can regenerate them
    # and make sure that the regenerated files are what we want
    states_files =  [fname for fname in os.listdir(stateverdir) if fname.endswith(".states")]
    for fname in states_files:
        #print("deleting", os.path.join(stateverdir,fname))
        os.system("rm " + os.path.join(stateverdir,fname))

    spec = MergeH3H4DesignSpecies()
    states = MergeBBStateVersion(state_ver_opts, spec)

    assert states.nstates_total() == 14

    states_files =  [fname for fname in os.listdir(stateverdir) if fname.endswith(".states")]
    assert len(states_files) == 5
    with open(os.path.join(stateverdir,"MH3_WTH4.states")) as fid:
        lines = fid.readlines()
    assert len(lines) == 4
    assert lines[0] == "1KX5_chAB_0331.pdb MH3_WTH4.corr MH3_WTH4.2resfile\n"

    for fname in states_files:
        #print("deleting", os.path.join(stateverdir,fname))
        os.system("rm " + os.path.join(stateverdir,fname))

