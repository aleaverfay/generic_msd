from generic_msd.msd_interface_design import DesignSpecies, DesignDefinitionOpts, IsolateBBDesDefFnames
from generic_msd.tests.dummy_desdef import dummy_design_species

import os


def test_create_design_species():
    spec = dummy_design_species()
    assert spec is not None


class empty_class:
    pass


def test_create_desdef_options():
    opts_mock = empty_class()
    setattr(opts_mock, "des_def", "desdef1_only_hphobes")
    setattr(opts_mock, "base_dir", "test123")
    desdef_opts = DesignDefinitionOpts(opts_mock)
    assert hasattr(desdef_opts, "des_def")
    assert hasattr(desdef_opts, "base_dir")


def test_desdef_fnames():
    currpath = os.path.dirname(os.path.abspath(__file__))
    basedir = currpath + "/dummy/"

    opts_mock = empty_class()
    setattr(opts_mock, "base_dir", basedir)
    setattr(opts_mock, "des_def", "desdef1_only_hphobes")
    desdef_opts1 = DesignDefinitionOpts(opts_mock)
    setattr(opts_mock, "des_def", "desdef1_only_hphobes/")
    desdef_opts2 = DesignDefinitionOpts(opts_mock)

    spec = dummy_design_species()

    fnames1 = IsolateBBDesDefFnames(desdef_opts1, spec)
    assert (
        fnames1.desdef_dir
        == basedir + "input_files/design_definitions/desdef1_only_hphobes"
    )
    fnames2 = IsolateBBDesDefFnames(desdef_opts2, spec)
    assert (
        fnames2.desdef_dir
        == basedir + "input_files/design_definitions/desdef1_only_hphobes/"
    )

    corr_gold = {
        "MH3_MH4": "MH3_MH4.corr",
        "MH3_WTH4": "MH3_WTH4.corr",
        "WTH3_MH4": "WTH3_MH4.corr",
        "MH3_p_MH4": "MH3_MH4.corr",
        "MH3_p_WTH4": "MH3_WTH4.corr",
        "WTH3_p_MH4": "WTH3_MH4.corr",
    }

    twores_gold = {
        "MH3_MH4": "design_both.2resfile",
        "MH3_WTH4": "design_A.2resfile",
        "WTH3_MH4": "design_B.2resfile",
        "MH3_p_MH4": "design_both_sep.2resfile",
        "MH3_p_WTH4": "design_A_sep.2resfile",
        "WTH3_p_MH4": "design_B_sep.2resfile",
    }
    assert fnames1.corr == corr_gold
    assert fnames1.secresfiles == twores_gold
    assert fnames1.entfunc == "h3h4_native.entfunc"
