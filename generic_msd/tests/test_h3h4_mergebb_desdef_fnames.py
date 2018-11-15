from generic_msd.msd_interface_design import MergeBBDesignSpecies, MergeBBDesDefFnames, DesignDefinitionOpts, DesignSpecies
from generic_msd.tests.mergebb_design_species import MergeH3H4DesignSpecies
import os



class mock_opts:
    pass

def test_mergebb_design_species_read_yaml_design_fnames():
    opts = mock_opts()
    setattr(opts, "des_def", "dd1")
    basedir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "merge_bb_inputs")
    setattr(opts, "base_dir", basedir)
    desdef_opts = DesignDefinitionOpts(opts)

    spec = MergeH3H4DesignSpecies()
    fnames = MergeBBDesDefFnames(desdef_opts, spec)

    assert len(fnames.corr) == 5
    assert len(fnames.secresfiles) == 5
    assert fnames.corr["MH3_MH4"] == "MH3_MH4.corr"
    assert fnames.secresfiles["MH3_MH4"] == "MH3_MH4.2resfile"
    assert fnames.entfunc == "h3h4_native.entfunc"
    
