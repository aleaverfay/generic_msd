from generic_msd.msd_interface_design import IsolateBBDesignSpecies
import toolz.functoolz


class dummy_design_species(IsolateBBDesignSpecies):
    @toolz.functoolz.memoize
    def species(self):
        return [
            "MH3_MH4",
            "MH3_WTH4",
            "WTH3_MH4",
            "MH3_p_MH4",
            "MH3_p_WTH4",
            "WTH3_p_MH4",
        ]

    def is_negative_species(self, spec):
        return "WT" in spec

    def is_positive_species(self, spec):
        return not self.is_negative_species(spec)

    def is_separated_species(self, spec):
        return "_p_" in spec

    def is_chAwt(self, spec):
        return "WTH3" in spec

    def is_chBwt(self, spec):
        return "WTH4" in spec

    def both_chains_mutated(self, spec):
        return "MH3" in spec and "MH4" in spec

    def sep_species_for_complex(self, spec):
        assert not self.is_separated_species(spec)
        assert spec in self.species()
        cols = spec.split("_")
        return cols[0] + "_p_" + cols[1]
