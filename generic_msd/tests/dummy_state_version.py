from generic_msd.msd_interface_design import (
    DesignSpecies,
    DesignDefinitionOpts,
    DesDefFnames,
    StateVersion,
    StateVersionOpts,
)


class dummy_state_version(StateVersion):
    def __init__(self, opts: StateVersionOpts, design_species: DesignSpecies):
        super().__init__(opts, design_species)
        self.all_pdbs = set([])
        self.pos_states = {}
        self.sep_states_both_mut = {}
        self.sep_states_wtA = {}
        self.sep_states_wtB = {}
        self.neg_states_wtA = {}
        self.neg_states_wtB = {}

    def pdbs(self):
        if len(self.all_pdbs) == 0:
            self.determine_pdbs()
        return self.all_pdbs

    def valid_mut_statuses(self):
        return set(["chAwt", "chBwt"])

    def pos_states(self):
        return self.pos_states

    def add_positive_complex(self, bbname, complex_name):
        self.pos_states[complex_name] = bbname
        self.all_pdbs.add(complex_name)

    def add_positive_sep(self, bbname, sep_pdb):
        self.sep_states_both_mut[sep_pdb] = bbname
        self.sep_states_wtA[sep_pdb] = bbname
        self.sep_states_wtB[sep_pdb] = bbname
        self.all_pdbs.add(sep_pdb)

    def add_negative_complex(self, mut_status, bbname, complex_pdb):
        assert mut_status == "chAwt" or mut_status == "chBwt"
        if mut_status == "chAwt":
            self.neg_states_wtA[complex_pdb] = bbname
        else:
            self.neg_states_wtB[complex_pdb] = bbname
        self.all_pdbs.add(complex_pdb)

    def add_negative_sep(self, mut_status, bbname, sep_pdb):
        assert mut_status == "chAwt" or mut_status == "chBwt"
        if mut_status == "chAwt":
            self.sep_states_wtA[sep_pdb] = bbname
        else:
            self.sep_states_wtB[sep_pdb] = bbname
        self.all_pdbs.add(sep_pdb)

    def separate_pdb(self, complex_pdb, sep_pdb):
        # separate_h3h4_interface.separate_h3h4_int(
        #    self.state_version_dir + complex_pdb,
        #    self.state_version_dir + sep_pdb )
        if not hasattr(self, "separations"):
            self.separations = []
        self.separations.append((complex_pdb, sep_pdb))

    def pdbs_for_target_bb(_, bb_for_pdbs, bbname):
        return [pdb for pdb, pdbbb in bb_for_pdbs.items() if pdbbb == bbname]

    def pdbs_for_species_for_bb(self, spec, bbname):
        pdbs = []

        if self.design_species.is_separated_species(spec):
            if self.design_species.is_chAwt(spec):
                pdbs.extend(self.pdbs_for_target_bb(self.sep_states_wtA, bbname))
            elif self.design_species.is_chBwt(spec):
                pdbs.extend(self.pdbs_for_target_bb(self.sep_states_wtB, bbname))
            else:
                pdbs.extend(self.pdbs_for_target_bb(self.sep_states_both_mut, bbname))
        else:
            # we're looking at a complex
            pdbs.extend(self.pdbs_for_target_bb(self.pos_states, bbname))
            if self.design_species.is_negative_species(spec):
                if self.design_species.is_chAwt(spec):
                    pdbs.extend(self.pdbs_for_target_bb(self.neg_states_wtA, bbname))
                else:
                    assert self.design_species.is_chBwt(spec)
                    pdbs.extend(self.pdbs_for_target_bb(self.neg_states_wtB, bbname))

        return pdbs

    def nstates_total(self):
        return (
            len(self.pos_states) * 3
            + len(self.sep_states_both_mut)
            + len(self.sep_states_wtA)
            + len(self.sep_states_wtB)
            + len(self.neg_states_wtA)
            + len(self.neg_states_wtB)
        )
