from generic_msd.msd_interface_design import MergeBBDesignSpecies, MergeBBDesDefFnames, DesignDefinitionOpts, DesignSpecies

class MergeH3H4DesignSpecies(MergeBBDesignSpecies):
    def species(self):
        return [ "MH3_MH4", "MH3_WTH4", "WTH3_MH4", "MH3", "MH4" ]
    def is_complex(self, spec):
        return "_" in spec
    def is_monomer(self, spec):
        return not self.is_complex(spec)
