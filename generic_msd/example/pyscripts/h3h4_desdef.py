from generic_msd.msd_interface_design import MergeBBDesignSpecies
import toolz.functoolz


class H3H4TwoMutDesignSpecies(MergeBBDesignSpecies):
    @toolz.functoolz.memoize
    def species(self):
        return [
            "AH3_AH4",
            "BH3_BH4",
            "AH3_BH4",
            "BH3_AH4",
            "AH3_WTH4",
            "WTH3_AH4",
            "BH3_WTH4",
            "WTH3_BH4",
            "AH3m",
            "BH3m",
            "AH4m",
            "BH4m",
        ]

    def is_complex(self, spec):
        return "_" in spec
    def is_monomer(self, spec):
        return not self.is_complex(spec)
