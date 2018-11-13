from generic_msd.msd_interface_design import (
    MSDIntDesJobOptions,
    IsolateBBInterfaceMSDJob,
    IsolateBBDesignSpecies,
    IsolateBBDesDefFnames,
    IsolateBBStateVersion,
    DesignDefinitionOpts,
    StateVersionOpts,
    PostProcessingOpts,
)
from generic_msd.tests.dummy_desdef import dummy_design_species
from generic_msd.tests.dummy_state_version import dummy_state_version
from generic_msd.opt_holder import OptHolder
import blargs


class H3H4InterfaceMSDJob(IsolateBBInterfaceMSDJob):
    @staticmethod
    def add_options(p: blargs.Parser):
        DesignDefinitionOpts.add_options(p)
        StateVersionOpts.add_options(p)
        MSDIntDesJobOptions.add_options(p)
        PostProcessingOpts.add_options(p)

    def __init__(self, all_opts):
        self.all_opts = all_opts
        msd_opts = MSDIntDesJobOptions(all_opts)
        super(H3H4InterfaceMSDJob, self).__init__(msd_opts)

    def files_to_symlink(self, subdir):
        return super().files_to_symlink(subdir)

    def states_to_save(self):
        ds = self.design_species
        return [x for x in ds.species() if not ds.is_separated_species(x)]

    def complexes_to_postprocess(self):
        return self.states_to_save()

    ###################
    # Factory methods #
    ###################
    def create_design_species(self) -> IsolateBBDesignSpecies:
        return dummy_design_species()

    def create_desdef_fnames(self, design_species: IsolateBBDesignSpecies) -> IsolateBBDesDefFnames:
        dd_opts = DesignDefinitionOpts(self.all_opts)
        return IsolateBBDesDefFnames(dd_opts, design_species)

    def create_state_version(self, design_species: IsolateBBDesignSpecies) -> IsolateBBStateVersion:
        sv_opts = StateVersionOpts(self.all_opts)
        return dummy_state_version(sv_opts, design_species)

    def create_post_processing_options(self) -> PostProcessingOpts:
        return PostProcessingOpts(self.all_opts)
