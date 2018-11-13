from generic_msd.msd_interface_design import (
    IsolateBBDesignSpecies,
    DesignDefinitionOpts,
    IsolateBBDesDefFnames,
    StateVersion,
    StateVersionOpts,
)
import generic_msd.tests.separate_h3h4_interface
from generic_msd.tests.dummy_desdef import dummy_design_species
from generic_msd.tests.dummy_state_version import dummy_state_version
import os



class empty_class:
    pass


def test_create_state_version_opts():
    currpath = os.path.dirname(os.path.abspath(__file__))
    basedir = currpath + "/dummy/"

    dummy_opts = empty_class()
    setattr(dummy_opts, "base_dir", basedir)
    setattr(dummy_opts, "state_version", "frwt_v1mock_dd1")
    state_ver_opts = StateVersionOpts(dummy_opts)
    assert state_ver_opts.state_version_dir == "frwt_v1mock_dd1"
    assert state_ver_opts.base_dir == basedir


def test_create_state_version():
    currpath = os.path.dirname(os.path.abspath(__file__))
    basedir = currpath + "/dummy/"

    dummy_opts = empty_class()
    setattr(dummy_opts, "base_dir", basedir)
    setattr(dummy_opts, "state_version", "frwt_v1mock_dd1")
    setattr(dummy_opts, "des_def", "desdef1_only_hphobes")

    desdef_opts = DesignDefinitionOpts(dummy_opts)
    state_ver_opts = StateVersionOpts(dummy_opts)

    spec = dummy_design_species()
    fnames = IsolateBBDesDefFnames(desdef_opts, spec)
    state_ver = StateVersion(state_ver_opts, spec)

    assert hasattr(state_ver, "backbone_names")
    assert hasattr(state_ver, "negbackbone_names")
    assert hasattr(state_ver, "design_species")
    assert state_ver.design_species is spec


def test_state_version_determine_pdbs():
    currpath = os.path.dirname(os.path.abspath(__file__))
    basedir = currpath + "/dummy/"

    dummy_opts = empty_class()
    setattr(dummy_opts, "base_dir", basedir)
    setattr(dummy_opts, "state_version", "frwt_v1mock_dd1")
    setattr(dummy_opts, "des_def", "desdef1_only_hphobes")

    desdef_opts = DesignDefinitionOpts(dummy_opts)
    state_ver_opts = StateVersionOpts(dummy_opts)

    spec = dummy_design_species()
    fnames = IsolateBBDesDefFnames(desdef_opts, spec)
    state_ver = dummy_state_version(state_ver_opts, spec)

    state_ver.determine_pdbs()
    state_ver.backbone_names == ["rwt_0982", "rwt_0331"]
    gold_pdbs = {
        "1KX5_chAB_0982_sep.pdb",
        "1KX5_chAB_0331.pdb",
        "1KX5_chAB_0331_sep.pdb",
        "dock_test_frwt_v0_dd1_ff1a_run2_1dG_1.0Ent_6_MH3_WTH4.pdb",
        "dock_test_frwt_v0_dd1_ff1a_run2_1dG_1.0Ent_2_WTH3_MH4.pdb",
        "1KX5_chAB_0982.pdb",
        "dock_test_frwt_v0_dd1_ff1a_run2_1dG_1.0Ent_3_MH3_WTH4.pdb",
        "dock_test_frwt_v0_dd1_ff1a_run2_1dG_1.0Ent_6_WTH3_MH4.pdb",
    }
    assert state_ver.pdbs() == gold_pdbs
    gold_pos_states = {
        "1KX5_chAB_0982.pdb": "rwt_0982",
        "1KX5_chAB_0331.pdb": "rwt_0331",
    }
    assert state_ver.pos_states == gold_pos_states

    gold_neg_states_wtA = {
        "dock_test_frwt_v0_dd1_ff1a_run2_1dG_1.0Ent_2_WTH3_MH4.pdb": "wtA1",
        "dock_test_frwt_v0_dd1_ff1a_run2_1dG_1.0Ent_6_WTH3_MH4.pdb": "wtA2",
    }
    assert state_ver.neg_states_wtA == gold_neg_states_wtA
    gold_neg_states_wtB = {
        "dock_test_frwt_v0_dd1_ff1a_run2_1dG_1.0Ent_3_MH3_WTH4.pdb": "wtB1",
        "dock_test_frwt_v0_dd1_ff1a_run2_1dG_1.0Ent_6_MH3_WTH4.pdb": "wtB2",
    }
    assert state_ver.neg_states_wtB == gold_neg_states_wtB


def test_state_version_generate_states_file():
    currpath = os.path.dirname(os.path.abspath(__file__))
    basedir = currpath + "/dummy/"

    dummy_opts = empty_class()
    setattr(dummy_opts, "base_dir", basedir)
    setattr(dummy_opts, "state_version", "frwt_v1mock_dd2")
    setattr(dummy_opts, "des_def", "desdef1_only_hphobes")

    desdef_opts = DesignDefinitionOpts(dummy_opts)
    state_ver_opts = StateVersionOpts(dummy_opts)

    spec = dummy_design_species()
    fnames = IsolateBBDesDefFnames(desdef_opts, spec)
    state_ver = dummy_state_version(state_ver_opts, spec)

    assert (
        state_ver.state_version_dir
        == os.path.join(basedir,"input_files/state_versions/frwt_v1mock_dd2")
    )
    svdir = state_ver.state_version_dir

    bbs = ["rwt_0982", "wtA1", "wtA2", "wtB1", "wtB2"]
    wtAbbs = ["rwt_0982", "wtA1", "wtA2"]
    wtBbbs = ["rwt_0982", "wtB1", "wtB2"]

    # Cleanup from previous runs?
    for species in spec.species():
        for bb in bbs:
            os.system("rm -f " + svdir + species + "_for_" + bb + ".states")

    state_ver.determine_pdbs()
    state_ver.create_state_file_lists()

    assert os.path.isfile(os.path.join(svdir, "MH3_MH4_for_rwt_0982.states"))
    for bb in wtAbbs:
        assert os.path.isfile(os.path.join(svdir, "WTH3_MH4_for_" + bb + ".states"))
    for bb in wtBbbs:
        assert os.path.isfile(os.path.join(svdir, "MH3_WTH4_for_" + bb + ".states"))

    with open(os.path.join(svdir, "MH3_MH4_for_rwt_0982.states")) as fid:
        lines = fid.readlines()
        gold_lines = ["1KX5_chAB_0982.pdb MH3_MH4.corr MH3_MH4.2resfile\n"]
        assert lines == gold_lines

    with open(os.path.join(svdir, "MH3_MH4_for_rwt_0982.states")) as fid:
        lines = fid.readlines()
        gold_lines = ["1KX5_chAB_0982.pdb MH3_MH4.corr MH3_MH4.2resfile\n"]
        assert lines == gold_lines

    for bb in wtAbbs:
        with open(os.path.join(svdir, "WTH3_MH4_for_" + bb + ".states")) as fid:
            lines = fid.readlines()
            gold_lines = [
                ("%s WTH3_MH4.corr WTH3_MH4.2resfile\n" % pdb)
                for pdb, pdbbb in state_ver.pos_states.items()
                if pdbbb == bb
            ]
            gold_lines.extend(
                [
                    ("%s WTH3_MH4.corr WTH3_MH4.2resfile\n" % pdb)
                    for pdb, pdbbb in state_ver.neg_states_wtA.items()
                    if pdbbb == bb
                ]
            )
            assert lines == gold_lines
    for bb in wtAbbs:
        with open(os.path.join(svdir, "WTH3_p_MH4_for_" + bb + ".states")) as fid:
            lines = fid.readlines()
            gold_lines = [
                ("%s WTH3_p_MH4.corr WTH3_p_MH4.2resfile\n" % pdb)
                for pdb, pdbbb in state_ver.sep_states_wtA.items()
                if pdbbb == bb
            ]
            assert lines == gold_lines

    for bb in wtBbbs:
        with open(os.path.join(svdir, "MH3_WTH4_for_" + bb + ".states")) as fid:
            lines = fid.readlines()
            gold_lines = [
                ("%s MH3_WTH4.corr MH3_WTH4.2resfile\n" % pdb)
                for pdb, pdbbb in state_ver.pos_states.items()
                if pdbbb == bb
            ]
            gold_lines.extend(
                [
                    ("%s MH3_WTH4.corr MH3_WTH4.2resfile\n" % pdb)
                    for pdb, pdbbb in state_ver.neg_states_wtB.items()
                    if pdbbb == bb
                ]
            )
            assert lines == gold_lines
    for bb in wtBbbs:
        with open(os.path.join(svdir, "MH3_p_WTH4_for_" + bb + ".states")) as fid:
            lines = fid.readlines()
            gold_lines = [
                ("%s MH3_p_WTH4.corr MH3_p_WTH4.2resfile\n" % pdb)
                for pdb, pdbbb in state_ver.sep_states_wtB.items()
                if pdbbb == bb
            ]
            assert lines == gold_lines

    # Cleanup
    for species in spec.species():
        for bb in bbs:
            os.system("rm -f " + svdir + species + "_for_" + bb + ".states")


def test_state_version_pdbs_function():

    currpath = os.path.dirname(os.path.abspath(__file__))
    basedir = currpath + "/dummy/"

    dummy_opts = empty_class()
    setattr(dummy_opts, "base_dir", basedir)
    setattr(dummy_opts, "state_version", "frwt_v1mock_dd1")
    setattr(dummy_opts, "des_def", "desdef1_only_hphobes")

    desdef_opts = DesignDefinitionOpts(dummy_opts)
    state_ver_opts = StateVersionOpts(dummy_opts)

    spec = dummy_design_species()
    fnames = IsolateBBDesDefFnames(desdef_opts, spec)
    state_ver = dummy_state_version(state_ver_opts, spec)

    gold_pdbs = set(
        [
            "dock_test_frwt_v0_dd1_ff1a_run2_1dG_1.0Ent_3_MH3_WTH4.pdb",
            "1KX5_chAB_0982_sep.pdb",
            "1KX5_chAB_0331.pdb",
            "1KX5_chAB_0982.pdb",
            "1KX5_chAB_0331_sep.pdb",
            "dock_test_frwt_v0_dd1_ff1a_run2_1dG_1.0Ent_6_WTH3_MH4.pdb",
            "dock_test_frwt_v0_dd1_ff1a_run2_1dG_1.0Ent_2_WTH3_MH4.pdb",
            "dock_test_frwt_v0_dd1_ff1a_run2_1dG_1.0Ent_6_MH3_WTH4.pdb",
        ]
    )
    assert state_ver.pdbs() == gold_pdbs


def test_state_version_nstates_total():

    currpath = os.path.dirname(os.path.abspath(__file__))
    basedir = currpath + "/dummy/"

    dummy_opts = empty_class()
    setattr(dummy_opts, "base_dir", basedir)
    setattr(dummy_opts, "state_version", "frwt_v1mock_dd1")
    setattr(dummy_opts, "des_def", "desdef1_only_hphobes")

    desdef_opts = DesignDefinitionOpts(dummy_opts)
    state_ver_opts = StateVersionOpts(dummy_opts)

    spec = dummy_design_species()
    fnames = IsolateBBDesDefFnames(desdef_opts, spec)
    state_ver = dummy_state_version(state_ver_opts, spec)
    state_ver.determine_pdbs()

    assert state_ver.nstates_total() == 20
