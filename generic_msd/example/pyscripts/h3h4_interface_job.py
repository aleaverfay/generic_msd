from generic_msd.msd_interface_design import (
    MergeBBInterfaceMSDJob,
    MergeBBDesignSpecies,
    MergeBBDesDefFnames,
    MergeBBStateVersion,
    MSDIntDesJobOptions,
    DesignDefinitionOpts,
    StateVersionOpts,
    PostProcessingOpts,
)
from h3h4_desdef import H3H4TwoMutDesignSpecies
# from h3h4_state_version import H3H4TwoMutStateVersion ??
from generic_msd.opt_holder import OptHolder
import blargs


class H3H4TwoMutInterfaceMSDJob(MergeBBInterfaceMSDJob):
    @staticmethod
    def add_options(p: blargs.Parser):
        DesignDefinitionOpts.add_options(p)
        StateVersionOpts.add_options(p)
        MSDIntDesJobOptions.add_options(p)
        PostProcessingOpts.add_options(p)

    def __init__(self, all_opts):
        self.all_opts = all_opts
        msd_opts = MSDIntDesJobOptions(all_opts)
        super().__init__(msd_opts)

    ###################
    # Factory methods #
    ###################
    def create_design_species(self) -> MergeBBDesignSpecies:
        return H3H4TwoMutDesignSpecies()

    def create_desdef_fnames(self, design_species: MergeBBDesignSpecies) -> MergeBBDesDefFnames:
        dd_opts = DesignDefinitionOpts(self.all_opts)
        return MergeBBDesDefFnames(dd_opts, design_species)

    def create_state_version(self, design_species: MergeBBDesignSpecies) -> MergeBBStateVersion:
        sv_opts = StateVersionOpts(self.all_opts)
        return MergeBBStateVersion(sv_opts, design_species)

    def create_post_processing_options(self) -> PostProcessingOpts:
        return PostProcessingOpts(self.all_opts)
