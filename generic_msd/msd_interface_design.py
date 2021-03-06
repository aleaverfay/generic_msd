import blargs
import os
import math
import traceback
import yaml
import itertools

# The goal of the functionality provided in these classes is to make the work
# requried to write the script for preparing and launching a set of multistate
# design jobs relatively easy.
# 
# The "setup_jobs.py" script will look something like this:
# 
#     import (some things)
# 
#     si = ServerIdentifier()
#     opts = OptHolder()
#     with blargs.Parser(opts) as p:
#         DerivedMSDJob.add_options(p)
#         JobExecutionOptions.add_options(p)
#
#     msd_job = DerivedMSDJob(opts)
#     msd_job_manager = MSDJobManager(msd_job, opts, si)
#     msd_job_manager.prepare_job()
#
# and that's it.
#
# The programmer will define a class DerivedMSDJob, derived from either
# the IsolateBBInterfaceMSDJob class or the MergeBBInterfaceMSDJob, and
# with the DerivedMSDJob, a DerivedDesignSpecies class, derived from
# either the IsolateBBDesignSpecies or the MergeBBDesignSpecies classes,
# and a DerivedStateVersion class, derived from either the 
# IsolateBBStateVersion or MergeBBStateVersion classes. These three
# classes will specify how the design job should proceed in broad
# strokes, where their base classes will do the heavy lifting.


def resolve_abs_path(fname):
    """Given a file name, resolve the absolute path for that file name
    from the current working directory"""
    if len(fname) > 0 and fname[0] == "/":
        return fname
    else:
        print(os.getcwd(), fname)
        return os.path.join(os.getcwd(), fname)

def add_required_opts_to_blargs_parser(
    blargs_parser: blargs.Parser,
    opt_args
):
    """Add a set of required options to the blargs Parser
    where the functions and their arguments are specified
    in opt_args.
    
    Opt args should be a list of tuples; each tuple is
      1 the type of option, as a string
      2 the name of the option, as a string
      3 a list of tuples; each element is
          1 the name of an additional function to call, as a string
          2 a tuple of the arguments to pass to that function
    """
    add_opts_to_blargs_parser(blargs_parser, opt_args, True)

def add_opts_to_blargs_parser(
    blargs_parser: blargs.Parser,
    opt_args,
    required: bool
):
    """Add a set of options to the blargs Parser
    where the functions and their arguments are specified
    in opt_args.
    
    Opt args should be a list of tuples; each tuple is
      1 the type of option, as a string
      2 the name of the option, as a string
      3 a list of tuples; each element is
          1 the name of an additional function to call, as a string
          2 a tuple of the arguments to pass to that function
    """
    for opt in opt_args:
        m = getattr(blargs_parser, opt[0])
        blarg_opt = m(opt[1])
        for exmethod, args in opt[2]:
            m = getattr(blarg_opt, exmethod)
            m(*args)
        if required:
            blarg_opt.required()


class StateVersionOpts:
    """The set of options that are needed to define the state version used in the design job"""

    @classmethod
    def add_options(cls, blargs_parser: blargs.Parser):
        cls.add_required_options(blargs_parser)
        cls.add_non_required_options(blargs_parser)
        
    @classmethod
    def add_required_options(cls, blargs_parser: blargs.Parser):
        add_required_opts_to_blargs_parser(blargs_parser, cls.required_opts())

    @classmethod
    def required_opts(cls):
        return [("str","state_version", [("shorthand",("s",))])]

    @classmethod
    def add_non_required_options(cls, blargs_parser: blargs.Parser):
        pass

    def __init__(self, opts):
        self.state_version_dir = opts.state_version
        self.base_dir = opts.base_dir

    def to_command_line(self):
        args = []
        args.append("--state_version")
        args.append(self.state_version_dir)
        return " ".join(args)


class DesignDefinitionOpts:
    """The set of options that are needed to define the design definition used in the design job"""

    def __init__(self, opts):
        self.des_def = opts.des_def
        self.base_dir = opts.base_dir

    @classmethod
    def add_options(cls, blargs_parser: blargs.Parser):
        cls.add_required_options(blargs_parser)
        cls.add_non_required_options(blargs_parser)

    @classmethod
    def add_required_options(cls, blargs_parser: blargs.Parser):
        add_required_opts_to_blargs_parser(blargs_parser, cls.required_opts())

    @classmethod
    def required_opts(cls):
        return [("str","des_def",[("shorthand",("d",))])]

    @classmethod
    def add_non_required_options(cls, blargs_parser: blargs.Parser):
        pass

    def to_command_line(self):
        args = []
        args.append("--des_def")
        args.append(self.des_def)
        return " ".join(args)


class DesignSpecies:
    """This is a class that represents the set of species that are represented
    in the design simulation. It's a very simple class, but fundamental for the
    behavior of a lot of other classes."""

    def species(self):
        """Return a list of strings representing all the species in the system"""
        raise NotImplementedError()


class DesDefFnames:
    """The DesDefFNames class loads the information for each species about the
    correspondence file and secondary resfile that should be used for it. It loads
    the name of the entity resfile that will be used. It
    also stores the entity function that should be used (if any). This is
    the base class; the derived classes do most of the heavy lifting."""

    def __init__(self, opts: DesignDefinitionOpts, design_species: DesignSpecies):
        self.desdef_dir = os.path.join(
            opts.base_dir,
            "input_files/design_definitions/",
            opts.des_def,
        )
        self.corr = {}
        self.secresfiles = {}
        self.entfunc = ""
        self.species = design_species
        self.n_entities = self._read_nentities_from_entity_resfile()

    def _read_nentities_from_entity_resfile(self):
        with open(os.path.join(self.desdef_dir, "entity.resfile")) as fid:
            lines = fid.readlines()
        return int(lines[0].strip())


class StateVersion:
    def __init__(self, opts: StateVersionOpts, design_species: DesignSpecies):
        self.state_version_dir = os.path.join(
            opts.base_dir, "input_files/state_versions/", opts.state_version_dir
        )
        self.design_species = design_species

    def nstates_total(self):
        """Return the (integer) number of individual states that will be
        simultaneously optimized by the multistate design algorithm. This
        number will be used by the MSDJobManager to determine the number
        of CPUs to request for running the each msd job"""
        raise NotImplementedError()

    def pdbs(self):
        """Return the complete list of all PDBs (without their path) that are to
        be used in the design simulation. These will be sim-linked by the
        MSDJobManager into each directory where the MSD jobs will run"""
        raise NotImplementedError()

    def pose_energy_vector_lists(self):
        """Return the list of files that will be used in a POSE_ENERGY_VECTOR
        command of a fitness file for the relevant conformations that will
        not be used in design"""
        raise NotImplementedError()

# class DesignDefinition:
#    """Base class from which particular design defintion classes will derive"""
#
#    pass


class MSDIntDesJobOptions:
    def __init__(self, cl_opts):
        self.base_dir = cl_opts.base_dir
        self.ntop_msd_results_to_dock = cl_opts.ntop_msd_results_to_dock
        self.job_name = cl_opts.job_name
        self.flags_files = cl_opts.flags_files
        self.w_dGdiff_bonus_weights_file = cl_opts.w_dGdiff_bonus_weights_file
        self.daf = cl_opts.daf
        self.entfunc_weights_file = cl_opts.entfunc_weights_file
        self.preserve_DAF = cl_opts.preserve_DAF
        self.positive_states = cl_opts.positive_states
        self.ngen_scale = cl_opts.ngen_scale
        self.seed_sequences = cl_opts.seed_sequences
        self.pdb_seed_pairs = cl_opts.pdb_seed_pairs
        self.single_round = cl_opts.single_round
        self.fill_gen1_from_seeds = cl_opts.fill_gen1_from_seeds
        self.pop_size = cl_opts.pop_size

    @classmethod
    def add_options(cls, blargs_parser: blargs.Parser):
        cls.add_non_required_options(blargs_parser)
        cls.add_required_options(blargs_parser)

        
    @classmethod
    def add_non_required_options(cls, blargs_parser: blargs.Parser):
        p = blargs_parser
        p.int("ntop_msd_results_to_dock").default(1)
        p.multiword("flags_files").shorthand("f").default("").cast(lambda x: x.split())
        p.str("w_dGdiff_bonus_weights_file").shorthand("w")
        p.str("entfunc_weights_file").shorthand("e")
        p.flag("preserve_DAF").shorthand("p")
        p.str("positive_states").default("")
        p.float("ngen_scale").default(15)
        ss = p.multiword("seed_sequences").cast(lambda x: x.split())
        p.multiword("pdb_seed_pairs").cast(lambda x: x.split()).conflicts(ss)
        sr = p.flag("single_round")
        p.flag("fill_gen1_from_seeds").conflicts(sr)
        p.int("pop_size").default(100)

    @classmethod
    def add_required_options(cls, blargs_parser: blargs.Parser):
        add_required_opts_to_blargs_parser(blargs_parser, cls.required_opts())

    @classmethod
    def required_opts(cls):
        return [("str", "job_name", [("shorthand", ("j",))]),
                ("str", "daf", [("shorthand", ("a",))])]

    def to_command_line(self):
        args = []
        args.append("--job_name")
        args.append(self.job_name)
        args.append("--daf")
        args.append(self.daf)
        args.append("--ntop_msd_results_to_dock")
        args.append(str(self.ntop_msd_results_to_dock))
        if self.flags_files:
            args.append("--flags_files")
            args.extend(self.flags_files)
        if self.w_dGdiff_bonus_weights_file:
            args.append("--w_dGdiff_bonus_weights_file")
            args.append(self.w_dGdiff_bonus_weights_file)
        if self.entfunc_weights_file:
            args.append("--entfunc_weights_file")
            args.append(self.entfunc_weights_file)
        if self.preserve_DAF:
            args.append("--preserve_DAF")
        if self.positive_states != "":
            args.append("--positive_states")
            args.append(self.positive_states)
        args.append("--ngen_scale")
        args.append(str(self.ngen_scale))
        if self.seed_sequences != 0:
            args.append("--seed_sequences")
            args.extend(self.seed_sequences)
        if self.pdb_seed_pairs:
            args.append("--pdb_seed_pairs")
            args.extend(self.pdb_seed_pairs)
        if self.single_round:
            args.append("--single_round")
        if self.fill_gen1_from_seeds:
            args.append("--fill_gen1_from_seeds")
        args.append("--pop_size")
        args.append(str(self.pop_size))

        return " ".join(args)


class PostProcessingOpts:
    def __init__(self, cl_opts):
        self.base_dir = cl_opts.base_dir
        self.docking_flags_files = [resolve_abs_path(x) for x in cl_opts.docking_flags_files]
        self.relax = cl_opts.relax
        self.docking_n_cpu = cl_opts.docking_n_cpu
        self.docking_queue = cl_opts.docking_queue

    @classmethod
    def add_options(cls, blargs_parser: blargs.Parser):
        cls.add_required_options(blargs_parser)
        cls.add_non_required_options(blargs_parser)

    @classmethod
    def add_required_options(cls, blargs_parser: blargs.Parser):
        add_required_opts_to_blargs_parser(blargs_parser, cls.required_opts())

    @classmethod
    def required_opts(cls):
        return []

    @classmethod
    def add_non_required_options(cls, blargs_parser: blargs.Parser):
        p = blargs_parser
        p.multiword("docking_flags_files").default("").cast(lambda x: x.split())
        p.flag("relax")
        p.int("docking_n_cpu")
        p.str("docking_queue")

    def to_command_line(self):
        args = []
        if self.docking_flags_files:
            args.append("--docking_flags_files")
            args.extend(self.docking_flags_files)
        if self.relax:
            args.append("--relax")
        if self.docking_n_cpu:
            args.append("--docking_n_cpu")
            args.append(str(self.docking_n_cpu))
        if self.docking_queue:
            args.append("--docking_queue")
            args.append(self.docking_queue)
        return " ".join(args)


class InterfaceMSDJob:

    ##########################################
    # interface presented to the MSDJobManager
    ##########################################

    def job_name(self):
        return self.job_name_

    def subjobs(self):
        return self.default_subjobs()

    def files_to_symlink(self, subdir):
        raise NotImplementedError()

    def popsize(self):
        if self.single_round:
            return len(self.seeds())
        else:
            return self.pop_size_

    def ngen(self):
        """Default implementation scales the number of generations by the number of
        designable positions, which is given by the design definition"""
        if self.single_round:
            return 1
        else:
            return int(math.ceil(self.desdef_fnames.n_entities * self.ngen_scale))

    def msd_flags(self, subdir):
        return self.flags_files

    def seeded(self, subdir):
        return (
            (self.seed_sequences and len(self.seed_sequences) != 0)
            or (self.pdb_seed_pairs and len(self.pdb_seed_pairs) != 0)
            or self.fill_gen1_from_seeds_
        )
    def seeds(self):
        if self.seed_sequences:
            return self.seed_sequences
        # oh, man, there's no code here for dealing with PDB seed pairs.
        # there needs to be!
        raise NotImplementedError()

    def fill_gen1_from_seeds(self):
        return self.fill_gen1_from_seeds_

    def states_to_save(self):
        raise NotImplementedError()

    def complexes_to_postprocess(self):
        raise NotImplementedError()

    def n_results_to_postprocess(self):
        return self.ntop_msd_results_to_dock

    def fitness_lines(self, subdir):
        raise NotImplementedError()

    def extra_post_processing_options(self):
        """Return a string containign additional options to be
        passed to dock_jobs_run.py"""
        return ""

    #######################################################
    # Interface between base class and concrete derived class
    #######################################################

    def create_design_species(self) -> DesignSpecies:
        """Factory method for a derived DesignSpecies class"""
        raise NotImplementedError()

    def create_desdef_fnames(self, design_species: DesignSpecies) -> DesDefFnames:
        """Factory method for a (possibly?) derived DesDefFNames class"""
        raise NotImplementedError()

    def create_state_version(self, design_species: DesignSpecies) -> StateVersion:
        """Factory method for a derived StateVersion class"""
        raise NotImplementedError()

    def create_post_processing_options(self) -> PostProcessingOpts:
        raise NotImplementedError()

    ################

    def __init__(self, msd_opts: MSDIntDesJobOptions):
        self.design_species = self.create_design_species()
        self.desdef_fnames = self.create_desdef_fnames(self.design_species)
        self.state_version = self.create_state_version(self.design_species)
        self.post_processing_opts = self.create_post_processing_options()

        self.ntop_msd_results_to_dock = msd_opts.ntop_msd_results_to_dock
        self.base_dir = msd_opts.base_dir
        self.job_name_ = msd_opts.job_name
        self.flags_files = [resolve_abs_path(x) for x in msd_opts.flags_files]
        self.w_dGdiff_bonus_weights_file = resolve_abs_path(
            msd_opts.w_dGdiff_bonus_weights_file
        ) if msd_opts.w_dGdiff_bonus_weights_file else ""
        self.daf = resolve_abs_path(msd_opts.daf)
        self.entfunc_weights_file = resolve_abs_path(msd_opts.entfunc_weights_file) if msd_opts.entfunc_weights_file else ""
        self.preserve_DAF = msd_opts.preserve_DAF
        self.positive_states = msd_opts.positive_states
        self.ngen_scale = msd_opts.ngen_scale
        self.seed_sequences = msd_opts.seed_sequences
        self.pdb_seed_pairs = msd_opts.pdb_seed_pairs
        self.single_round = msd_opts.single_round
        self.fill_gen1_from_seeds_ = msd_opts.fill_gen1_from_seeds
        self.pop_size_ = msd_opts.pop_size

    def entfunc_weights_from_file(self):
        if self.entfunc_weights_file != "":
            with open(self.entfunc_weights_file) as fid:
                lines = fid.readlines()
            entfunc_weights = [float(line.strip()) for line in lines]
        else:
            entfunc_weights = [1.0]
        return entfunc_weights

    def dGdiff_bonus_weights_from_file(self):
        if self.w_dGdiff_bonus_weights_file != "":
            with open(self.w_dGdiff_bonus_weights_file) as fid:
                lines = fid.readlines()
            dGdiff_bonus_weights = [float(line.strip()) for line in lines]
        else:
            dGdiff_bonus_weights = [
                float(x) for x in ("5", "10", "15", "20", "25", "30")
            ]
        return dGdiff_bonus_weights

    def default_subjobs(self):
        entfunc_weights = self.entfunc_weights_from_file()
        dGdiff_bonus_weights = self.dGdiff_bonus_weights_from_file()

        names = [
            "{0}_{1:.1f}w_dGdiff_{2:.1f}Ent".format(self.job_name_, w_dGdiff, w_entfunc)
            for w_dGdiff in dGdiff_bonus_weights
            for w_entfunc in entfunc_weights
        ]

        return names

    def dG_bonus_weight_from_subdir_name(self, subdir):
        cols = subdir.split("_")
        return float(cols[-3][:-1])

    def entfunc_weight_from_subdir_name(self, subdir):
        cols = subdir.split("_")
        return float(cols[-1][:-3])

    def replace_wscan_and_add_entfunc_to_fitness_func(self, subdir, orig_lines):
        """Replace the "%(WSCAN)s" string and add in an ENTITY_FUNCTION line
        ahead of the FITNESS line if the design definition files includes and
        entity function."""

        dG_bonus_weight = self.dG_bonus_weight_from_subdir_name(subdir)
        w_ent_func = self.entfunc_weight_from_subdir_name(subdir)
        replace_dict = {"WSCAN": "%f" % dG_bonus_weight}

        newlines = []
        for line in orig_lines:
            if line[0] != "#":
                # append the entity function line to the input daf
                if line[0:7] == "FITNESS" and self.desdef_fnames.entfunc:
                    efunc_fname = self.desdef_fnames.entfunc
                    if efunc_fname.find("/") != -1:
                        efunc_fname = efunc_fname.rpartition("/")[2]
                    newlines.append("ENTITY_FUNCTION entfunc " + efunc_fname + "\n")
                    newlines.append("\n")
                    line = line[:-1] + " + %f * entfunc\n" % w_ent_func
                newlines.append(line % replace_dict)
            else:
                newlines.append(line)
        return newlines


############################################################################################
############################################################################################
############################################################################################
### The idea of the "IsolateBB" classes is that they are going to wrap up a set of MSD jobs
### that only compare energies between two states if they come from the "same backbone"
### i.e. the internal DOFs of their backbones are identical. This will allow you to compare
### the energies between the complex and separated-complex, e.g., to determine a binding
### energy.
###
### For each backbone, there is a "complex" and a "separated complex" pair. The different
### species will be thread across these pairs to allow a binding energy for that species
### on that backbone. It is possible to have many different rigid-body docked orientations
### for a particular backbone, so that a single state for the separated complex may be used
### (it would be expensive/wasteful to have a state for each of the complexes that differ
### only by their rigid-body orientation). If there is only a single backbone (e.g. the
### wild type backbone), that is simply handled as a special case of there being multiple
### backbones.
###
### The fitness functions will often look to the complex with the largest binding energy.
###
### The IsolateBBInterfaceMSDJob class will generate a large portion of the fitness-
### function file, reading from the state version class to list all of the different
### backbones for a species, and then to define dG_bind values for each backbone/species
### combo, putting all of the dGbinds into a vector variable (vdGbind_<complex>).


class IsolateBBDesignSpecies(DesignSpecies):
    """Derived IsolateBBDesignSpecies classes identify what species are being modeled
    and whether those species are complexes or separated complexes and
    whether those species are positive or negative species"""

    ###############################
    #### Functions to override ####
    ###############################
    def species(self):
        """Return a list of strings representing all the species in the system"""
        raise NotImplementedError()

    def is_negative_species(self, spec):
        """Return whether a particular species is a negative species"""
        raise NotImplementedError()

    def is_positive_species(self, spec):
        """Return whether a particular species is a positive species"""
        raise NotImplementedError()

    def is_separated_species(self, spec):
        """Return whether a particular species represents a separated pair"""
        raise NotImplementedError()

    def sep_species_for_complex(self, spec):
        """Return the name of the separated species corresponding to a complex"""
        raise NotImplementedError()

    #################################
    # Simple helpers based on above #
    # functions #####################
    #################################

    def complexes(self):
        """Return the list of species that represent complexes"""
        return [spec for spec in self.species() if not self.is_separated_species(spec)]


class IsolateBBDesDefFnames(DesDefFnames):
    """The DesDefFNames class loads the information for each species about the
    correspondence file and secondary resfile that should be used for it. It
    also stores the entity function that should be used (if any)."""

    def __init__(self, opts: DesignDefinitionOpts, design_species: DesignSpecies):
        super(IsolateBBDesDefFnames, self).__init__(opts, design_species)

        fname = os.path.join(self.desdef_dir, "definition_files.txt")
        with open(fname) as fid:
            lines = fid.readlines()
        self._read_from_lines(fname, lines)
        self.n_entities = self._read_nentities_from_entity_resfile()

    def _read_from_lines(self, desdef_fname, lines):
        """Load the correspondence file names and secondary-resfile file names
        for the species defined in definition_files.txt. These files will be
        referred to in the .states files as <spec>.corr and <spec>.2res, and
        the aliasing will be handled by symlinks"""

        assert len(self.corr) == 0
        assert len(self.secresfiles) == 0
        for line in lines:
            if line == "\n":
                continue
            cols = line.split()
            header = cols[0]

            if header in self.species.species():
                i = 1
                while i < len(cols):
                    if cols[i] == "corr=":
                        self.corr[header] = cols[i + 1]
                        i += 2
                    elif cols[i] == "2res=":
                        self.secresfiles[header] = cols[i + 1]
                        i += 2
                    else:
                        raise ValueError(
                            "ERROR reading",
                            desdef_fname,
                            "at column",
                            i,
                            "on line:\n",
                            line.strip("\n"),
                        )

            elif header == "entfunc":
                self.entfunc = cols[1]
            else:
                print("allowed species:", self.species.species())
                raise ValueError(
                    "ERROR reading",
                    desdef_fname,
                    "unrecognized command",
                    header,
                    "on line:\n",
                    line.strip("\n"),
                )


class IsolateBBStateVersion(StateVersion):
    """Derived IsolateBBStateVersion classes must track which backbone should be used for
    which species so that the .states files can be constructed. It will rely
    on the IsolateBBDesignSpecies class to identify which are positive and which are
    negative states. The base class will read three files in the state-version
    directory that state what PDB should be used in which context: the
    "pos_backbones.list" file, the "neg_complexes.list" file and the
    "neg_backbones.list" file.

    The "pos_backbones.list" file will list all the backbones and PDBs that
    should be used by both positive and negative species. It is important
    that for the PDBs listed in this file that the wild-type sequence be
    present in every position that is un-designed in any species. E.g., if
    a complex between chains A and B is given here, and there is some species
    that looks at the mutant-A--WT-B interaction, then all positions on chain B
    must be wild type. If there is also some species that looks at the
    WT-A--mutant-B interaction, then chain A must also be wild type.

    The "neg_backbones.list" file will list new backbone conformations that
    should be used only for negative species. Each backbone may appear only
    a single time. It must give an identifier for each backbone, its
    "mut_status," that will tell the derived class which negative species
    this backbone should be used for. This "mut_status" identifier carries
    no meaning to the StateVersion base class and is merely passed along
    to the derived class. The StateVersion base class will, however, verify
    that the provided mut_status identifier is valid using the set returned
    by the derived class's "valid_mut_statuses" method. When giving a new
    negative backbone, both a complex and a separated complex PDB must be
    given.

    The "neg_complexes.list" file will list all the backbones and PDBs that
    should be used only for negative species, and must state that any PDB
    within represents another conformation of a backbone first identified
    in either the "pos_backbones.list" or "neg_backbones.list" files, and
    must also state WHICH mutant species this negaive conformation should
    be used for using a "mut_status" id, (described above). Only the
    complex PDB should be listed for PDBs in this file -- the separated
    complex PDB that it corresponds to will be the one for the same backbone
    listed in either the "pos_backbones.list" or "neg_backbones.list" files.
    """

    def __init__(self, opts: StateVersionOpts, design_species: IsolateBBDesignSpecies):
        super(IsolateBBStateVersion, self).__init__(opts, design_species)
        self.backbone_names = set([])
        self.negbackbone_names = set([])
        self._pdbs_determined = False

    ###############################################
    # Functions invoked by the base class which the
    # derived class must provide
    ###############################################
    def valid_mut_statuses(self):
        """Return a set of valid mut_status identifiers (strings) that can
        appear in the neg_complexes.list and neg_backbones.list files.
        These identifiers are what the derived class will use to decide
        which species the backbones should be used with, since not all
        backbones can necessarily be used with all species"""
        raise NotImplementedError()

    def pdbs(self):
        """Return the complete list of all PDBs (without their path) that are to
        be used in the design simulation. These will be sim-linked by the
        MSDJobManager into each directory where the MSD jobs will run"""
        raise NotImplementedError()

    def pose_energy_vector_lists(self):
        """Base class provides an effective no-op"""
        return []


    def pdbs_for_species_for_bb(self, spec, bbname):
        """Return the list of PDBs for a particular species that all have the same backbone.
        This is the main function that the base class needs to invoke in order to create
        the .states files. Once the .states files have been created, though, it is not
        invoked again later. The derived class will be told what PDB files fall in what
        category of use by the base class through the methods below:
            * add_positive_complex
            * add_positive_sep
            * add_negative_complex
            * add_negative_sep
        along with the "mut status" identifier provided in the state definition input
        file. The derived class must then use this information to later report which
        PDB file should be used with which species when this function is called.
        """
        raise NotImplementedError()

    def add_positive_complex(self, bbname, complex_pdb):
        """Store a new positive complex, identified by its backbone name. The
        derived class must decide which species should use this complex
        for later invocations of pdbs_for_species_for_bb.
        """
        raise NotImplementedError()

    def add_positive_sep(self, bbname, sep_pdb):
        """Store a new positive separated complex, identified by its backbone name.
        The derived class must decide which species should use this separated complex
        for later invocations of pdbs_for_species_for_bb.
        """
        raise NotImplementedError()

    def add_negative_complex(self, mut_status, bbname, complex_pdb):
        """Store a negative complex, identified by its mut status and backbone name.
        The derived class must decide which species should use this complex for
        later invocations of pdbs_for_species_for_bb.
        """
        raise NotImplementedError()

    def add_negative_sep(self, mut_status, bbname, sep_pdb):
        """Store a negative separated complex, identified by its mut status
        and backbone name. The derived class must decide which species should
        use this complex for later invocations of pdbs_for_species_for_bb
        """
        raise NotImplementedError()

    def separate_pdb(self, complex_pdb, sep_pdb):
        """Write a new PDB file to disk representing the separated structure of a
        particular complex. The derived class must read in the complex pdb
        and then separate the chains appropriately."""
        raise NotImplementedError()

    def state_file_name(self, spec, bbname):
        """Return the name of the state file for a particular species
        with a particular backbone. Derived classes are welcome to override
        this method."""
        return spec + "_for_" + bbname + ".states"

    def nstates_total(self):
        """Return the (integer) number of individual states that will be
        simultaneously optimized by the multistate design algorithm. This
        number will be used by the MSDJobManager to determine the number
        of CPUs to request for running the each msd job"""
        raise NotImplementedError()

    ####################################################################
    # Functions that this class implements to create the state files
    # These functions are invoked by the IsolateBBInterfaceMSDJob
    # class
    ####################################################################
    def determine_pdbs(self):
        if self._pdbs_determined:
            return
        self.determine_pdbs_from_pos_backbones()
        self.determine_pdbs_from_neg_backbones()
        self.determine_pdbs_from_neg_complexes()
        self._pdbs_determined = True

    def create_state_file_lists(self):
        for spec in self.design_species.species():
            for bbname in self.backbone_names:
                self.create_state_file_list_for_spec_and_bb(spec, bbname)
            if self.design_species.is_positive_species(spec):
                continue
            # negative backbones also need state.list files for the negbackbone set
            for bbname in self.negbackbone_names:
                # print "Creating state.list file for negative backbone bbname"
                self.create_state_file_list_for_spec_and_bb(spec, bbname)

    #####################################################################
    # Functions that do the heavy lifting for this class
    #####################################################################

    def create_state_file_list_for_spec_and_bb(self, spec, bbname):
        """This function is tasked with creating the .states file for a combination
        of species and backbone names. It will read the contents of the existing
        .states file (if any) and generate the contents of a new .states file,
        and if they are different, then it will write the new contents to disk"""

        # print "create_state_file_list_for_spec_and_bb ", spec, bbname
        dirname = self.state_version_dir
        old_lines = None
        if os.path.isfile(os.path.join(dirname, self.state_file_name(spec, bbname))):
            old_lines = open(
                os.path.join(dirname, self.state_file_name(spec, bbname))
            ).readlines()
        else:
            old_lines = []  # perhaps the file doesn't exist because it shouldn't exist
        lines = self.state_file_lines_for_spec_and_bb(spec, bbname)
        if lines != old_lines and len(lines) > 0:
            print(
                "Writing list file for", self.state_file_name(spec, bbname), "species"
            )
            open(
                os.path.join(dirname, self.state_file_name(spec, bbname)), "w"
            ).writelines(lines)

    def state_file_lines_for_spec_and_bb(self, spec, bbname):
        """This function will be given a species and then asked to
        construct the contents of the state-file for that species. It
        relies on the pdbs_for_species_for_bb function"""

        pdbs = self.pdbs_for_species_for_bb(spec, bbname)
        return [(pdb + " " + spec + ".corr " + spec + ".2resfile\n") for pdb in pdbs]

    ###########
    def check_input_pdb_files(self, complex_name, separated, source, line):
        """Verify that both the complex pdb and the separated-complex pdb are
        present on disk. If the complex pdb is absent, we cannot continue.
        If the separated pdb is absent, and if the name for the separated pdb
        is the name of the complex pdb ending with "_sep", then we'll ask
        the derived class to construct that pdb using the "separate_pdb"
        function, which should write the new pdb to disk"""
        if not os.path.isfile(os.path.join(self.state_version_dir, complex_name)):
            raise FileNotFoundError(
                "Could not find file "
                + os.path.join(self.state_version_dir, complex_name)
                + " which was given in "
                + source
                + "On the line: "
                + line
            )
        if not os.path.isfile(os.path.join(self.state_version_dir, separated)):
            if (
                len(separated) > 8
                and separated[-4:] == ".pdb"
                and separated[-8:] == "_sep.pdb"
                and len(complex_name) > 4
                and separated[:-8] == complex_name[:-4]
            ):
                self.separate_pdb(complex_name, separated)
            else:
                raise FileNotFoundError(
                    "Could not find file "
                    + os.path.join(self.state_version_dir, separated)
                    + " which was given in "
                    + source
                    + "On the line: "
                    + line
                )

    def determine_pdbs_from_pos_backbones(self):
        dir = self.state_version_dir
        pos_backbones = os.path.join(dir, "pos_backbones.list")
        lines = [x.strip() for x in open(pos_backbones).readlines()]
        for line in lines:
            if len(line) == 0:
                continue
            if line[0] == "#":
                continue
            cols = line.split()
            bbname = cols[0]
            # if bbname in self.backbone_names:
            #     raise ValueError(
            #         "Backbone with name: "
            #         + bbname
            #         + " was already given; but appears again on the line:"
            #         + line)
            assert cols[1] == "complex="
            complex_name = cols[2]
            assert cols[3] == "separated="
            separated = cols[4]

            self.check_input_pdb_files(
                complex_name, separated, "pos_backbones.list", line
            )

            self.backbone_names.add(bbname)
            self.add_positive_complex(bbname, complex_name)
            self.add_positive_sep(bbname, separated)

    def determine_pdbs_from_neg_backbones(self):
        """TO DO: Document format of the "neg_backbones.list" file"""
        dir = self.state_version_dir
        valid_statuses = self.valid_mut_statuses()
        neg_backbones = os.path.join(dir, "neg_backbones.list")
        if not os.path.isfile(neg_backbones):
            return
        lines = [x.strip() for x in open(neg_backbones).readlines()]
        for line in lines:
            if len(line) == 0:
                continue
            if line[0] == "#":
                continue
            cols = line.split()
            mut_status = cols[0]
            if mut_status not in valid_statuses:
                raise ValueError(
                    'Mut_status (first column) given on line: "'
                    + line.strip("\n")
                    + '" of file '
                    + neg_backbones
                    + " is invalid\n"
                    + "Valid mut statuses are:"
                    + ", ".join(valid_statuses)
                )
            bbname = cols[1]
            if bbname in self.backbone_names:
                raise ValueError(
                    "Backbone with name: "
                    + bbname
                    + " was already given as a postive backbone; but appears again on the line:"
                    + line.strip("\n")
                )
            if bbname in self.negbackbone_names:
                raise ValueError(
                    "Negative backbone with name: "
                    + bbname
                    + " was already given; but appears again on the line:"
                    + line.strip("\n")
                )
                assert bbname not in self.negbackbone_names
            assert cols[2] == "complex="
            complex_name = cols[3]
            assert cols[4] == "separated="
            separated = cols[5]

            self.check_input_pdb_files(
                complex_name, separated, "neg_backbones.list", line
            )

            self.negbackbone_names.add(bbname)
            self.add_negative_complex(mut_status, bbname, complex_name)
            self.add_negative_sep(mut_status, bbname, separated)

    def determine_pdbs_from_neg_complexes(self):
        """TO DO: Document the format of the "neg_complexes.list" file"""
        dir = self.state_version_dir
        valid_statuses = self.valid_mut_statuses()
        neg_complexes = os.path.join(dir, "neg_complexes.list")
        lines = [x.strip() for x in open(neg_complexes).readlines()]
        for line in lines:
            if len(line) == 0:
                continue
            if line[0] == "#":
                continue
            cols = line.split()
            mut_status = cols[0]
            if mut_status not in valid_statuses:
                raise ValueError(
                    'Mut_status (first column) given on line: "'
                    + line.strip("\n")
                    + '" of file '
                    + neg_complexes
                    + " is invalid\n"
                    + "Valid mut statuses are:"
                    + ", ".join(valid_statuses)
                )
            bbname = cols[1]
            assert bbname in self.backbone_names or bbname in self.negbackbone_names
            complex_name = cols[2]
            assert os.path.isfile(os.path.join(dir, complex_name))
            self.add_negative_complex(mut_status, bbname, complex_name)


class IsolateBBInterfaceMSDJob(InterfaceMSDJob):
    """The IsolateBBInterfaceMSDJob class defines the set of jobs
    that will be run for the MSDJobManager. It's most difficult task
    is to construct the fitness file for each of the jobs.

    The fitness file it constructs is a combination of 1) lines it
    generates based on the state version, and 2) lines provided by
    the user that describe how the energies for the states should
    be combined. See the fitness_lines function for more detail.

    Classes that derive from IsolateBBInterfaceMSDJob need only
    implement the following functions:
        * create_design_species
        * create_desdef_fnames
        * create_state_version
        * create_post_processing_options
        * add_options (static; used by setup_jobs.py)
        * 
    """
    

    def __init__(self, msd_opts: MSDIntDesJobOptions):
        super(IsolateBBInterfaceMSDJob, self).__init__(msd_opts)
        assert isinstance(self.design_species, IsolateBBDesignSpecies)
        assert isinstance(self.desdef_fnames, IsolateBBDesDefFnames)
        assert isinstance(self.state_version, IsolateBBStateVersion)

    def files_to_symlink(self, subdir):
        return self.default_symlinkable_file_list()

    def default_symlinkable_file_list(self):
        """Setup the list of symlinked files that are defined by the state
        version, the design definition, and the flags files

        returns a list of tuples:
        0: original file path+name, and
        1: new file name ("." for same name)
        """
        input_files = []
        self.state_version.determine_pdbs()
        self.state_version.create_state_file_lists()

        stateverdir = self.state_version.state_version_dir
        desdefdir = self.desdef_fnames.desdef_dir

        for pdb in self.state_version.pdbs():
            input_files.append((os.path.join(stateverdir, pdb), "."))
        for spec in self.state_version.design_species.species():
            for bb in self.state_version.backbone_names:
                if not os.path.isfile(
                    os.path.join(
                        stateverdir, self.state_version.state_file_name(spec, bb)
                    )
                ):
                    continue
                input_files.append(
                    (
                        os.path.join(
                            stateverdir, self.state_version.state_file_name(spec, bb)
                        ),
                        ".",
                    )
                )
            if self.design_species.is_negative_species(spec):
                for bb in self.state_version.negbackbone_names:
                    if not os.path.isfile(
                        os.path.join(
                            stateverdir, self.state_version.state_file_name(spec, bb)
                        )
                    ):
                        continue
                    input_files.append(
                        (
                            os.path.join(
                                stateverdir,
                                self.state_version.state_file_name(spec, bb),
                            ),
                            ".",
                        )
                    )
        for spec in self.design_species.species():
            input_files.append(
                (os.path.join(desdefdir, self.desdef_fnames.corr[spec]), spec + ".corr")
            )
            input_files.append(
                (
                    os.path.join(desdefdir, self.desdef_fnames.secresfiles[spec]),
                    spec + ".2resfile",
                )
            )
        input_files.append((os.path.join(desdefdir, "entity.resfile"), "."))
        if self.desdef_fnames.entfunc != "":
            input_files.append(
                (os.path.join(desdefdir, self.desdef_fnames.entfunc), ".")
            )

        for fname in self.flags_files:
            if not os.path.isfile(fname):
                raise ValueError("Requested flags file", fname, "not found!")

        input_files.extend([(fname, ".") for fname in self.flags_files])
        return input_files

    def is_spec_and_bb_combo_valid(self, spec, bb):
        """Derived class is welcome to implement a faster version; hacky
        version just asks "does the state file for this species and backbone
        exist on the file system?" which is not a terribly efficient way
        of answering the question, and also requires that this be invoked
        only after the StateVersion has create all .states files"""
        state_file_name = self.state_version.state_file_name(spec, bb)
        return os.path.isfile(
            os.path.join(self.state_version.state_version_dir, state_file_name)
        )

    def fitness_lines(self, subdir):
        """Default fitness function definition that will, from the state version,
        construct the variables "v(spec)" for each of the species and
        vdGbind_(comp) for each of the complexes. It will paste the user-provided
        fitness function at the bottom of the state vector definitions that
        this function provides, and then it will edit the FITNESS function
        to replace the word "%(WSCAN)s" with the dG_bonus_weight for this subdirectory
        and also insert an ENTITY_FUNCTION line above the fitness line if
        the design definition has one, and will insert the result of the
        entity function into the FITNESS line, weighted with the entity-function
        weight that this subdir requests.

        This function is invoked by the MSDJobManager directly. Feel free to
        override it in the derived class.
        """
        # def fitness_lines(orig_lines, desdefnames, statedef, opts, weight_bonus, w_ent_func):

        orig_lines = None
        with open(self.daf) as fid:
            orig_lines = fid.readlines()

        newlines = []
        # header
        # declare the state vectors

        newlines.append("# species: " + " ".join(self.design_species.species()) + "\n")
        newlines.append(
            "# bbnames: "
            + " ".join(self.state_version.backbone_names)
            + (
                ""
                if not self.state_version.negbackbone_names
                else " ".join(self.state_version.negbackbone_names)
            )
            + "\n"
        )

        # Create STATE_VECTOR variables for all states
        for spec in self.design_species.species():
            for bb in self.state_version.backbone_names:
                if self.is_spec_and_bb_combo_valid(spec, bb):
                    line = (
                        "STATE_VECTOR "
                        + spec
                        + "_"
                        + bb
                        + " "
                        + self.state_version.state_file_name(spec, bb)
                        + "\n"
                    )
                    newlines.append(line)
            if self.design_species.is_negative_species(spec):
                for bb in self.state_version.negbackbone_names:
                    if self.is_spec_and_bb_combo_valid(spec, bb):
                        line = (
                            "STATE_VECTOR "
                            + spec
                            + "_"
                            + bb
                            + " "
                            + self.state_version.state_file_name(spec, bb)
                            + "\n"
                        )
                        newlines.append(line)
            newlines.append("\n")

        # declare sub expressions for the best energy of each species/bb pair
        for spec in self.design_species.species():
            for bb in self.state_version.backbone_names:
                if self.is_spec_and_bb_combo_valid(spec, bb):
                    line = (
                        "SCALAR_EXPRESSION best_"
                        + spec
                        + "_"
                        + bb
                        + " = vmin( "
                        + spec
                        + "_"
                        + bb
                        + ")\n"
                    )
                    newlines.append(line)
            if self.design_species.is_negative_species(spec):
                for bb in self.state_version.negbackbone_names:
                    if self.is_spec_and_bb_combo_valid(spec, bb):
                        line = (
                            "SCALAR_EXPRESSION best_"
                            + spec
                            + "_"
                            + bb
                            + " = vmin( "
                            + spec
                            + "_"
                            + bb
                            + ")\n"
                        )
                        newlines.append(line)
            newlines.append("\n")

        # now declare vector variables as follows:
        # best energy per species vectors v(spec)
        # dGbind_(comp) - one per complex
        newlines.append(
            "# best energies for a single species on each of its available backbones\n"
        )
        for spec in self.design_species.species():
            line = "VECTOR_VARIABLE v" + spec + " = "
            for bb in self.state_version.backbone_names:
                state_file_name = self.state_version.state_file_name(spec, bb)
                if not os.path.isfile(
                    os.path.join(self.state_version.state_version_dir, state_file_name)
                ):
                    continue
                line += "best_" + spec + "_" + bb + " "
            if self.design_species.is_negative_species(spec):
                for bb in self.state_version.negbackbone_names:
                    state_file_name = self.state_version.state_file_name(spec, bb)
                    if not os.path.isfile(
                        os.path.join(
                            self.state_version.state_version_dir, state_file_name
                        )
                    ):
                        continue
                    line += "best_" + spec + "_" + bb + " "
            newlines.append(line[:-1] + "\n")
        newlines.append("\n")
        for comp in self.design_species.complexes():
            sep = self.design_species.sep_species_for_complex(comp)
            newlines.append(
                "VECTOR_EXPRESSION FOR c IN v"
                + comp
                + " , s IN v"
                + sep
                + " :  vdGbind_"
                + comp
                + " = ( c - s ) * ite( lt( c - s, 0 ), 1, 0 ) \n"
            )
        newlines.append("\n")

        # note: stale # now the rest of the fitness function has access to
        # note: stale # the following variables and should be constructed
        # note: stale # from them:
        # note: stale # vMH3_MH4, vMH3_WTH4, vWTH3_MH4
        # note: stale # vMH3_p_MH4, vMH3_p_WTH4, vWTH3_p_MH4
        # note: stale # vdGbind_MH3_MH4, vdGBind_MH3_WTH4, vdGBind_WTH3_MH4

        remainder = self.replace_wscan_and_add_entfunc_to_fitness_func(subdir, orig_lines)

        return newlines + remainder


###################################################################################
###################################################################################
###################################################################################


class MergeBBDesignSpecies(DesignSpecies):
    """Derived MergeBBDesignSpecies classes identify a set of complexes and a set
    of monomers."""

    def species(self):
        raise NotImplementedError()

    def is_complex(self, spec):
        raise NotImplementedError()

    def is_monomer(self, spec):
        raise NotImplementedError()


class MergeBBDesDefFnames(DesDefFnames):
    def __init__(self, opts: DesignDefinitionOpts, design_species: DesignSpecies):
        super(MergeBBDesDefFnames, self).__init__(opts, design_species)

        fname = os.path.join(self.desdef_dir, "definition_files.yaml")
        self._read_from_file(fname)

    def _read_from_file(self, fname):
        with open(fname) as fid:
            raw = yaml.load(fid, Loader=yaml.Loader)
        for spec in raw["species"]:
            spec_name = spec["name"]
            if spec_name not in self.species.species():
                print("spec name not found:", spec_name)
            assert spec_name in self.species.species()
            spec_corr = spec["corr"]
            spec_2res = spec["2res"]
            print("species:", spec_name, spec_corr, spec_2res)
            assert os.path.isfile(os.path.join(self.desdef_dir, spec_corr))
            assert os.path.isfile(os.path.join(self.desdef_dir, spec_2res))
            self.corr[spec_name] = spec_corr
            self.secresfiles[spec_name] = spec_2res
        if "entfunc" in raw:
            self.entfunc = raw["entfunc"]
            assert os.path.isfile(os.path.join(self.desdef_dir,self.entfunc))
        assert os.path.isfile(os.path.join(self.desdef_dir,"entity.resfile"))


class MergeBBStateVersion(StateVersion):
    def __init__(self, opts: StateVersionOpts, design_species: MergeBBDesignSpecies):
        super(MergeBBStateVersion, self).__init__(opts, design_species)
        self._determined_pdbs = False
        self.count_n_states = 0
        self.pev_lists = []

    def pdbs(self):
        if not self._determined_pdbs:
            self.determine_pdbs()
        return self.all_pdbs

    def nstates_total(self):
        if not self._determined_pdbs:
            self.determine_pdbs()
        return self.count_n_states

    def pose_energy_vector_lists(self):
        return self.pev_lists

    def determine_pdbs(self):
        # read the "states.yaml" file
        # this will give the species for each PDB file
        # there may be more than one PDB file per species
        self.pdbs_for_spec = {}
        self.all_pdbs = set([])
        fname = os.path.join(self.state_version_dir, "states.yaml")
        with open(fname) as fid:
            raw = yaml.load(fid, Loader=yaml.Loader)
        for entry in raw["pdbs"]:
            self.count_n_states += 1
            spec = entry["species"]
            pdb = entry["pdb"]
            assert spec in self.design_species.species()
            assert os.path.isfile(os.path.join(self.state_version_dir, pdb))
            if spec not in self.pdbs_for_spec:
                self.pdbs_for_spec[spec] = []
            self.pdbs_for_spec[spec].append(pdb)
            self.all_pdbs.add(pdb)

        for entry in raw["pose_energy_vector_lists"]:
            fname = os.path.join(self.state_version_dir, entry)
            with open(fname) as fid:
                lines = fid.readlines()
                for line in lines:
                    if len(line) == 0:
                        continue
                    if line[0] == "#":
                        continue
                    pdb = line.strip()
                    self.all_pdbs.add(pdb)
            self.pev_lists.append(fname)

        # now, let's write the .states files if we haven't done so already
        for spec in self.design_species.species():
            states_fname = os.path.join(self.state_version_dir, spec + ".states")
            if not os.path.isfile(states_fname):
                lines = []
                for pdb in self.pdbs_for_spec[spec]:
                    lines.append("%s %s %s\n" % (pdb, spec + ".corr", spec + ".2resfile"))
                with open(states_fname, "w") as fid:
                    fid.writelines(lines)

        self.yaml_contents = raw
        self._determined_pdbs = True


class MergeBBInterfaceMSDJob(InterfaceMSDJob):
    def __init__(self, msd_opts: MSDIntDesJobOptions):
        super(MergeBBInterfaceMSDJob, self).__init__(msd_opts)
        assert isinstance(self.design_species, MergeBBDesignSpecies)
        assert isinstance(self.desdef_fnames, MergeBBDesDefFnames)
        assert isinstance(self.state_version, MergeBBStateVersion)

    def files_to_symlink(self, subdir):
        species = self.design_species.species()
        svdir = self.state_version.state_version_dir
        pdbs = [os.path.join(svdir, pdb) for pdb in self.state_version.pdbs()]
        state_files = [os.path.join(svdir, spec + ".states") for spec in species]
        pevs = [os.path.join(svdir, fname) for fname in self.state_version.pose_energy_vector_lists()]

        desdefdir = self.desdef_fnames.desdef_dir
        corr_pairs = [
            (os.path.join(desdefdir, self.desdef_fnames.corr[spec]), spec + ".corr")
            for spec in species
        ]
        secres_pairs = [
            (
                os.path.join(desdefdir, self.desdef_fnames.secresfiles[spec]),
                spec + ".2resfile",
            )
            for spec in species
        ]

        pdb_state_flag_and_pev_pairs = [
            (x, ".") for x in itertools.chain(pdbs, state_files, self.flags_files, pevs)
        ]

        extras = []
        extras.append(
            (os.path.join(self.desdef_fnames.desdef_dir, "entity.resfile"), ".")
        )
        if self.desdef_fnames.entfunc != "":
            extras.append(
                (
                    os.path.join(
                        self.desdef_fnames.desdef_dir, self.desdef_fnames.entfunc
                    ),
                    "."
                )
            )

        return list(
            itertools.chain(
                pdb_state_flag_and_pev_pairs, corr_pairs, secres_pairs, extras
            )
        )

    def states_to_save(self):
        """Instruct the MSD job manager to save all species by default"""
        return self.design_species.species()

    def complexes_to_postprocess(self):
        """Post process all complexes"""
        return [
            spec
            for spec in self.design_species.species()
            if self.design_species.is_complex(spec)
        ]

    def fitness_lines(self, subdir):
        """Replace the WSCAN and add in an ENTITY_FUNCTION line ahead of the FITNESS line
        if the design definition files includes and entity function.

        The MergeBBInterfaceMSDJob does a lot less fitness-function-boilerplate creation than
        the IsolateBBInterfaceMSDJob does; the input fitness funcion should be basically
        complete (and tailored specifically for the species listed in the MergeBBDesignSpecies)
        needing only the weight and entity function addition."""

        with open(self.daf) as fid:
            orig_lines = fid.readlines()
        return self.replace_wscan_and_add_entfunc_to_fitness_func(subdir, orig_lines)
