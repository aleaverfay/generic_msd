from generic_msd.msd_interface_design import MergeBBDesignSpecies, MergeBBDesDefFnames, DesignDefinitionOpts, DesignSpecies
import os


class MergeH3H4DesignSpecies(MergeBBDesignSpecies):
    def species(self):
        return [ "MH3_MH4", "MH3_WTH4", "WTH3_MH4", "MH3", "MH4", "WTH3", "WTH4" ]
    def is_complex(self, spec):
        return "_" in spec
    def is_monomer(self, spec):
        return not self.is_complex(spec)

class MergeH3H4DesDefFnames(MergeBBDesDefFnames):
    def __init__(self, opts: DesignDefinitionOpts, design_species: DesignSpecies):
        super(MergeH3H4DesDefFnames, self).__init__(opts, design_species)


class mock_opts:
    pass

def test_read_yaml_design_fnames():
    opts = mock_opts()
    setattr(opts, "des_def", "dd1")
    basedir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "merge_bb_inputs")
    setattr(opts, "base_dir", basedir)
    desdef_opts = DesignDefinitionOpts(opts)

    spec = MergeH3H4DesignSpecies()
    fnames = MergeH3H4DesDefFnames(desdef_opts, spec)

    
