# from optparse import OptionParser
import blargs
import os
import sys
import math
import socket
import separate_h3h4_interface
import dup_h2o_and_separate
import pdb_structure
import amino_acids as aa

from create_lsf_or_slurm_job import *


class FauxOptions:
    pass


# read the command line
def initialize_options():
    opt_dict = {}
    with blargs.Parser(opt_dict) as p:
        p.int("num_cpu").shorthand("n").default(-1)
        p.int("num_states_per_cpu").shorthand("m").default(-1)
        p.int("ntop_msd_results_to_dock").default(1)
        p.str("des_def").shorthand("d").required()
        p.str("state_def").shorthand("s").required()
        p.str("job_name").shorthand("j").required()
        p.multiword("flags_files").shorthand("f").default("").cast(lambda x: x.split())
        p.multiword("docking_flags_files").shorthand("o").default("").cast(
            lambda x: x.split()
        )
        p.str("w_dGdiff_bonus_weights_file").shorthand("w")  # weights file
        p.str("daf").shorthand("a")
        p.str("entfunc_weights_file").shorthand("e")
        p.flag("preserve_DAF").shorthand("p")
        p.flag("launch")
        p.flag("relax")
        # frelax = p.flag( "relax" )
        # p.flag( "relax_w_mpi" ).requires( frelax )
        # p.multiword( "relax_flags_files" ).cast( lambda x : x.split() ).requires( frelax ) # APL NOTE: use docking_flags_files instead
        # p.str( "relax_protocol" ).requires( frelax )
        p.str("positive_states").default("")
        p.str("queue").shorthand("q").default("week")
        p.int("ngen_scale").default(15)
        ss = p.multiword("seed_sequences").cast(lambda x: x.split())
        ps = p.multiword("pdb_seed_pairs").cast(lambda x: x.split()).conflicts(ss)
        sr = p.flag("single_round")
        p.flag("fill_gen1_from_seeds").conflicts(sr)
        p.int("pop_size").default(100)
        p.flag(
            "wet"
        )  # in generating the separated PDB file, should waters be duplicated?

    if "seed_sequences" not in opt_dict or opt_dict["seed_sequences"] is None:
        opt_dict["seed_sequences"] = []
    if "pdb_seed_pairs" not in opt_dict or opt_dict["pdb_seed_pairs"] is None:
        opt_dict["pdb_seed_pairs"] = []
    mylocals = locals()
    opts = FauxOptions()
    for i in opt_dict:
        # print "option: ", i
        # if i == "p" or i == "mylocals" or i == "help" or i == "ps" or i == "sr" or i == "ss" : continue
        setattr(opts, i, opt_dict[i])
        # print i, opt_dict[i]

    return opts


def base_dir():
    hostname = socket.gethostname()[:9]
    pwd_cols = os.getcwd().split("/")
    # print pwd_cols
    if on_killdevil():
        if "leaverfa" in pwd_cols:
            return "/nas02/home/l/e/leaverfa/ortho_nuc/h3h4/setup_scheme2/"
        elif "catson4" in pwd_cols:
            return "/nas02/home/c/a/catson4/ortho_nuc/h3h4/setup_scheme2/"
    elif on_dogwood:
        if "leaverfa" in pwd_cols:
            return "/nas/longleaf/home/leaverfa/ortho_nuc/h3h4/setup_scheme2/"
        elif "catson4" in pwd_cols:
            return "/nas/longleaf/home/catson4/ortho_nuc/h3h4/setup_scheme2/"
    elif hostname == "wiggins":
        return "/home/andrew/rosetta/orthogonal_nucleosome/h3h4/killdevil_setup/setup_scheme2/"
    else:
        # assume we're on crito
        return "/Users/andrew/rosetta/ortho_nuc/h3h4/killdevil_setup/setup_scheme2/"


def state_ver_dir(opt):
    return base_dir() + "input_files/state_versions/" + opt.state_def + "/"


def design_def_dir(opt):
    return base_dir() + "input_files/design_definitions/" + opt.des_def + "/"


class OrthoNucSpecies:
    def __init__(self):
        self.species = [
            "MH3_MH4",
            "MH3_WTH4",
            "WTH3_MH4",
            "MH3_p_MH4",
            "MH3_p_WTH4",
            "WTH3_p_MH4",
        ]
        self.complexes = ["MH3_MH4", "MH3_WTH4", "WTH3_MH4"]
        self.dual_mut_complexes = ["MH3_MH4"]
        self.wt_binding_complexes = ["MH3_WTH4", "WTH3_MH4"]
        self.pos_complexes = ["MH3_MH4"]
        self.pos_species_set = set(["MH3_MH4", "MH3_p_MH4"])
        self.neg_species_set = set(["MH3_WTH4", "WTH3_MH4", "MH3_p_WTH4", "WTH3_p_MH4"])

    def negative_species(self, spec):
        assert spec in self.species
        return spec in self.neg_species_set

    def positive_species(self, spec):
        assert spec in self.species
        return spec in self.pos_species_set

    def separated_species(self, spec):
        cols = spec.split("_")
        return len(cols) > 2 and cols[1] == "p"

    def chAwt(self, species):
        return species.split("_")[0][:-2] == "WT"

    def chBwt(self, species):
        if species in self.complexes:
            return species.split("_")[1][:-2] == "WT"
        else:
            return species.split("_")[2][:-2] == "WT"


class DesDefFnames:
    def __init__(self):
        self.corr = {}
        self.secresfiles = {}
        self.entfunc = ""
        self.species = OrthoNucSpecies()

    def read_from_lines(self, desdef_fname, lines):
        for line in lines:
            if line == "\n":
                continue
            cols = line.split()
            header = cols[0]

            if header in self.species.species:
                i = 1
                while i < len(cols):
                    if cols[i] == "corr=":
                        self.corr[header] = cols[i + 1]
                        i += 2
                    elif cols[i] == "2res=":
                        self.secresfiles[header] = cols[i + 1]
                        i += 2
                    else:
                        print(
                            "ERROR reading",
                            desdef_fname,
                            "at column",
                            i,
                            "on line:",
                            line,
                            end=" ",
                        )
                        sys.exit(1)

            elif header == "entfunc":
                self.entfunc = cols[1]
            else:
                print(
                    "ERROR reading",
                    desdef_fname,
                    "unrecognized command",
                    header,
                    "on line:",
                    line,
                    end=" ",
                )
                sys.exit(1)


class StateDefinitions:
    def __init__(self):
        self.species = OrthoNucSpecies()
        self.pdbs = set([])
        self.pos_states = {}  # map from pdb of positive state to bb name
        self.sep_states_both_mut = {}  # map from pdb of separated pdbs to bb name where both chains are mutating
        self.sep_states_wtA = {}  # map from pdb of separated pdbs to bb name where chain A is WT
        self.sep_states_wtB = {}  # map from pdb of separated pdbs to bb name where chain B is WT
        # self.neg_states_both_mut = {} # map from pdb of negative states to bb name where both chains are mutating
        self.neg_states_wtA = {}  # map from pdb of negative states to bb name where chain A is WT
        self.neg_states_wtB = {}  # map from pdb of negative states to bb name where chain B is WT
        self.backbone_names = []
        self.negbackbone_names = []  # extra backbones used only for the negative states

    def check_input_pdb_files(self, dir, complex_name, separated, source, line, opt):
        if not os.path.isfile(dir + "/" + complex_name):
            print(
                "Could not find file "
                + dir
                + "/"
                + complex_name
                + " which was given in "
                + source
            )
            print("On the line: " + line)
            assert os.path.isfile(dir + "/" + complex_name)
        if not os.path.isfile(dir + "/" + separated):
            print(
                "Could not find file "
                + dir
                + "/"
                + separated
                + " which was given in "
                + source
            )
            if (
                len(separated) > 8
                and separated[-4:] == ".pdb"
                and separated[-8:] == "_sep.pdb"
                and len(complex_name) > 4
                and separated[:-8] == complex_name[:-4]
            ):
                print("Writing separated pdb file", dir + separated)
                # complex_name_pdb = pdb_structure.pdbstructure_from_file( dir + "/" + complex_name )
                if opt.wet:
                    dup_h2o_and_separate.dup_h2o_and_sep_h3h4(
                        dir + complex_name, dir + separated
                    )
                else:
                    separate_h3h4_interface.separate_h3h4_int(
                        dir + complex_name, dir + separated
                    )
            else:
                print("On the line: " + line)
                assert os.path.isfile(dir + "/" + separated)

    def determine_pdbs_from_pos_backbones(self, opt):
        dir = state_ver_dir(opt)
        pos_backbones = dir + "pos_backbones.list"
        lines = [x.strip() for x in open(pos_backbones).readlines()]
        for line in lines:
            if len(line) == 0:
                continue
            if line[0] == "#":
                continue
            cols = line.split()
            bbname = cols[0]
            if bbname in self.backbone_names:
                print(
                    "Backbone with name: "
                    + bbname
                    + " was already given; but appears again on the line:"
                )
                print(line)
                assert bbname not in self.backbone_names
            assert cols[1] == "complex="
            complex_name = cols[2]
            assert cols[3] == "separated="
            separated = cols[4]

            self.check_input_pdb_files(
                dir, complex_name, separated, "pos_backbones.list", line, opt
            )

            self.backbone_names.append(bbname)
            self.pos_states[complex_name] = bbname
            self.sep_states_both_mut[separated] = bbname
            self.sep_states_wtA[separated] = bbname
            self.sep_states_wtB[separated] = bbname
            self.pdbs.add(complex_name)
            self.pdbs.add(separated)

    def determine_pdbs_from_neg_backbones(self, opt):
        dir = state_ver_dir(opt)
        neg_backbones = dir + "neg_backbones.list"
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
            assert mut_status == "chAwt" or mut_status == "chBwt"
            bbname = cols[1]
            if bbname in self.backbone_names:
                print(
                    "Backbone with name: "
                    + bbname
                    + " was already given as a postive backbone; but appears again on the line:"
                )
                print(line)
                assert bbname not in self.backbone_names
            if bbname in self.negbackbone_names:
                print(
                    "Negative backbone with name: "
                    + bbname
                    + " was already given; but appears again on the line:"
                )
                print(line)
                assert bbname not in self.negbackbone_names
            assert cols[2] == "complex="
            complex_name = cols[3]
            assert cols[4] == "separated="
            separated = cols[5]

            self.check_input_pdb_files(
                dir, complex_name, separated, "neg_backbones.list", line, opt
            )

            self.negbackbone_names.append(bbname)
            if mut_status == "chAwt":
                self.neg_states_wtA[complex_name] = bbname
                self.sep_states_wtA[separated] = bbname
            else:  # mut_status == "chBwt" :
                self.neg_states_wtB[complex_name] = bbname
                self.sep_states_wtB[separated] = bbname

            self.pdbs.add(complex_name)
            self.pdbs.add(separated)

    def determine_pdbs_from_neg_complexes(self, opt):
        dir = state_ver_dir(opt)
        neg_complexes = dir + "neg_complexes.list"
        lines = [x.strip() for x in open(neg_complexes).readlines()]
        for line in lines:
            if len(line) == 0:
                continue
            if line[0] == "#":
                continue
            cols = line.split()
            mut_status = cols[0]
            # print mut_status
            assert mut_status == "chAwt" or mut_status == "chBwt"
            bbname = cols[1]
            assert bbname in self.backbone_names or bbname in self.negbackbone_names
            complex_name = cols[2]
            assert os.path.isfile(dir + "/" + complex_name)
            if mut_status == "chAwt":
                self.neg_states_wtA[complex_name] = bbname
            else:  # mut_status == "chBwt" :
                self.neg_states_wtB[complex_name] = bbname
            self.pdbs.add(complex_name)

    def determine_pdbs(self, opt):
        self.determine_pdbs_from_pos_backbones(opt)
        self.determine_pdbs_from_neg_backbones(opt)
        self.determine_pdbs_from_neg_complexes(opt)

    # Figure out for a particular species what the set of PDBs are that it should use
    def pdbs_for_species_for_bb(self, spec, bbname):
        pdbs = []
        if self.species.separated_species(spec):
            print("species ", spec, " is a separated species")
            if self.species.chAwt(spec):
                for pdb in self.sep_states_wtA:
                    if self.sep_states_wtA[pdb] == bbname:
                        pdbs.append(pdb)
            elif self.species.chBwt(spec):
                for pdb in self.sep_states_wtB:
                    if self.sep_states_wtB[pdb] == bbname:
                        pdbs.append(pdb)
            else:
                for pdb in self.sep_states_both_mut:
                    if self.sep_states_both_mut[pdb] == bbname:
                        pdbs.append(pdb)
        else:
            print("species ", spec, " is a complex")
            neg = self.species.negative_species(spec)
            for pdb in self.pos_states:
                if self.pos_states[pdb] == bbname:
                    pdbs.append(pdb)
            if neg:
                if self.species.chAwt(spec):
                    for pdb in self.neg_states_wtA:
                        if self.neg_states_wtA[pdb] == bbname:
                            pdbs.append(pdb)
                else:
                    assert self.species.chBwt(spec)
                    for pdb in self.neg_states_wtB:
                        if self.neg_states_wtB[pdb] == bbname:
                            pdbs.append(pdb)
                # else :
                #    pass # we shouldn't get here, right?
                #    #for pdb in self.neg_states_both_mut :
                #    #    if self.neg_states_both_mut[ pdb ] == bbname : pdbs.append( pdb )

        return pdbs

    # This function will be given a species and then asked to construct the state-file
    # for that species. It relies on the pdbs_for_species_for_bb function
    def create_state_file_list_for_spec_and_bb(self, opt, spec, bbname):
        # print "create_state_file_list_for_spec_and_bb ", spec, bbname
        dirname = state_ver_dir(opt)
        old_lines = None
        if os.path.isfile(dirname + self.state_file_name(spec, bbname)):
            old_lines = open(dirname + self.state_file_name(spec, bbname)).readlines()
        else:
            old_lines = []  # perhaps the file doesn't exist because it shouldn't exist
        pdbs = self.pdbs_for_species_for_bb(spec, bbname)
        lines = []
        for pdb in pdbs:
            print(pdb)
            lines.append(pdb + " " + spec + ".corr " + spec + ".2resfile\n")
        if lines != old_lines and len(lines) > 0:
            print(
                "Writing list file for", self.state_file_name(spec, bbname), "species"
            )
            open(dirname + self.state_file_name(spec, bbname), "w").writelines(lines)

    def create_state_file_lists(self, opt):
        for spec in self.species.species:
            for bbname in self.backbone_names:
                self.create_state_file_list_for_spec_and_bb(opt, spec, bbname)
            if self.species.positive_species(spec):
                continue
            # negative backbones also need state.list files for the negbackbone set
            for bbname in self.negbackbone_names:
                # print "Creating state.list file for negative backbone bbname"
                self.create_state_file_list_for_spec_and_bb(opt, spec, bbname)

    def state_file_name(self, spec, bbname):
        return spec + "_for_" + bbname + ".states"


def read_positive_species_corr_files(design_def, opts):
    desdefdir = design_def_dir(opts)
    corrs = []
    for spec in design_def.species.pos_complexes:
        lines = open(desdefdir + "/" + design_def.corr[spec]).readlines()
        corr = {}
        for line in lines:
            if line[0] == "#":
                continue
            cols = line.split()
            corr[int(cols[0])] = (cols[1], cols[2])  # resid, chainid
        corrs.append(corr)
    return corrs


def seed_sequence_from_pdb_pair(corr_pair, pdbs):
    # print corr_pair
    print("corr_pair", corr_pair)
    seed_seq = ["X"] * (len(corr_pair[0]) * 2)
    for i in range(2):
        pdbstruct = pdb_structure.pdbstructure_from_file(pdbs[i])
        for j in range(len(corr_pair[i])):
            print(i, j)
            res, ch = corr_pair[i][2 * j + i + 1]
            seed_seq[2 * j + i] = aa.longer_names[pdbstruct.residue(ch, res).resname]
    ss = "".join(seed_seq)
    print(pdbs[0], "and", pdbs[1], ":", ss)
    return ss


def nprocs_for_job(state_def, opts):
    if opts.num_cpu < 0:
        nstates = 0
        # actually go and open the .list files that were generated from the (original) contents of
        # the state-definition directory
        for spec in state_def.species.species:
            for bb in state_def.backbone_names:
                dirname = state_ver_dir(opts)
                state_file_name = dirname + "/" + state_def.state_file_name(spec, bb)
                if not os.path.isfile(state_file_name):
                    continue
                lines = open(state_file_name).readlines()
                for line in lines:
                    if line != "\n":
                        nstates += 1
            if state_def.species.positive_species(spec):
                continue
            for bb in state_def.negbackbone_names:
                dirname = state_ver_dir(opts)
                state_file_name = dirname + "/" + state_def.state_file_name(spec, bb)
                if not os.path.isfile(state_file_name):
                    continue
                lines = open(state_file_name).readlines()
                for line in lines:
                    if line != "\n":
                        nstates += 1
        if opts.num_states_per_cpu < 0:
            # count the number of states in the state_def files
            return nstates
        else:
            return math.ceil(nstates / opts.num_states_per_cpu)
    return opts.num_cpu


def nentities_for_job(design_def, opts):
    lines = open(design_def_dir(opts) + "entity.resfile").readlines()
    return int(lines[0].strip())


def ngen_for_job(design_def, opts):
    return opts.ngen_scale * nentities_for_job(design_def, opts)


def validate_seed_sequences(design_def, opts):
    if len(opts.seed_sequences) == 0:
        return
    nents = nentities_for_job(design_def, opts)
    for ss in opts.seed_sequences:
        if len(ss) != nents:
            print(
                "Seed sequence",
                ss,
                "has the wrong length for the given design definition;",
                len(ss),
                "vs",
                nents,
            )
            assert len(ss) == nents


def extra_flags_fnames(opts):
    extraflags_files = opts.flags_files
    fnames = []
    for extraflags_file in extraflags_files:
        if extraflags_file.find("/") != -1:
            extraflags_file = extraflags_file.rpartition("/")[2]
            fnames.append(extraflags_file)
    return fnames


def designdir_name(opts, design_def, w_dGdiff, w_entfunc):
    name = opts.job_name + "_" + w_dGdiff + "dG"
    if design_def.entfunc != "":
        name += "_" + str(w_entfunc) + "Ent"
    return name


def create_job(
    state_def, design_def, dir_name, opts, w_dG, w_ent, launch_script, gather_script
):
    create_symlinks(opts, state_def, design_def)
    if opts.daf == "":
        orig_lines = default_daf()
    else:
        daf_fname = opts.daf
        if daf_fname[0] != "/":
            daf_fname = "../../" + daf_fname
        orig_lines = open(daf_fname).readlines()

    if opts.preserve_DAF:
        newdaf = orig_lines
    else:
        newdaf = fitness_lines(
            orig_lines, design_def, state_def, opts, float(w_dG), w_ent
        )
    nprocs = nprocs_for_job(state_def, opts)
    queue_name = opts.queue

    with open("fitness.daf", "w") as fid:
        fid.writelines(newdaf)

    # set both popsize and ngen to 1 if the user says to only run a single round of multistate design
    popsize = (
        str(max(opts.pop_size, len(opts.seed_sequences)))
        if not opts.single_round
        else str(len(opts.seed_sequences))
    )
    ngen = str(ngen_for_job(design_def, opts)) if not opts.single_round else "1"

    job_options = SubmissionOptions()
    job_options.scheduler = scheduler_type_for_server()
    job_options.num_nodes = nprocs
    job_options.queue = queue_name
    job_options.job_name = dir_name
    job_options.logfilename = dir_name + ".log"
    job_options.submission_script_fname = "submit.sh"
    job_options.mpi_job = True

    command_line = (
        mpi_msd_exe()
        + " -database "
        + db_path()
        + " -entity_resfile entity.resfile "
        + "-fitness_file fitness.daf -ms::pop_size "
        + popsize
        + " -ms::generations "
        + ngen
        + " -ms::numresults "
        + str(opts.ntop_msd_results_to_dock)
    )
    if not opts.flags_files is None:
        extraflags_files = extra_flags_fnames(opts)
        for extraflags_file in extraflags_files:
            command_line += " @../" + extraflags_file
    if len(opts.seed_sequences) != 0:
        command_line += " -seed_sequences " + " ".join(opts.seed_sequences)
    if opts.fill_gen1_from_seeds:
        command_line += " -fill_gen1_from_seed_sequences"

    submission_command, command_filename = command_and_submission_script_for_job(
        job_options, command_line
    )
    launch_script.append(submission_command + " >> ../msd_submission.log\n")
    # if command_filename : # no need to move the file -- we're already in the correct directory
    #    os.rename( command_filename, dir_name + "/" + command_filename )

    gather_script.append("cd " + dir_name + "\n")
    for spec in design_def.species.complexes:
        for result_ind in range(1, opts.ntop_msd_results_to_dock + 1):
            lzri = leading_zero_string(result_ind, opts.ntop_msd_results_to_dock)
            gather_script.append(
                "for i in `ls msd_output_%s_%s*`; do if [ ! -h $i ]; then cp $i ../results/%s_%s_%s.pdb; fi; done\n"
                % (str(result_ind), spec, dir_name, lzri, spec)
            )
    gather_script.append("cd ..\n")
    for result_ind in range(1, opts.ntop_msd_results_to_dock + 1):
        prefix = (
            " "
            + dir_name
            + "_"
            + leading_zero_string(result_ind, opts.ntop_msd_results_to_dock)
            + "_"
        )
        gather_script.append(
            "echo"
            + prefix
            + (".pdb" + prefix).join(design_def.species.complexes)
            + ".pdb >> results/job_triples.list\n"
        )
    gather_script.append("\n")


# link all of the files in the statedef and designdef
def create_symlinks(opts, statedefs, desdefnames):
    for pdb in statedefs.pdbs:
        os.system("ln -s " + state_ver_dir(opts) + pdb)
    for spec in statedefs.species.species:
        for bb in statedefs.backbone_names:
            if not os.path.isfile(
                state_ver_dir(opts) + statedefs.state_file_name(spec, bb)
            ):
                continue
            os.system(
                "ln -s " + state_ver_dir(opts) + statedefs.state_file_name(spec, bb)
            )
        if statedefs.species.negative_species(spec):
            for bb in statedefs.negbackbone_names:
                if not os.path.isfile(
                    state_ver_dir(opts) + statedefs.state_file_name(spec, bb)
                ):
                    continue
                os.system(
                    "ln -s " + state_ver_dir(opts) + statedefs.state_file_name(spec, bb)
                )
    for spec in desdefnames.species.species:
        os.system(
            "ln -s "
            + design_def_dir(opts)
            + desdefnames.corr[spec]
            + " "
            + spec
            + ".corr"
        )
        os.system(
            "ln -s "
            + design_def_dir(opts)
            + desdefnames.secresfiles[spec]
            + " "
            + spec
            + ".2resfile"
        )
    os.system("ln -s " + design_def_dir(opts) + "entity.resfile")
    if desdefnames.entfunc != "":
        os.system("ln -s " + design_def_dir(opts) + desdefnames.entfunc)


# What do I need to setup a bunch of runs?
# I need a fitness function
# I need a bunch of symlinks to:
#   a) pdb files
#   b) corr files
#   c) state files
#   d) secondary resfiles
#   e) the entity func file
# I need to write a submit script
# I need to write an output gatherng script


# Replace a value in the original fitness function file
# WSCAN: a weight that can be scanned by this script and
def fitness_lines(orig_lines, desdefnames, statedef, opts, weight_bonus, w_ent_func):
    newlines = []
    replace_dict = {"WSCAN": "%f" % weight_bonus}
    # header
    # declare the state vectors

    newlines.append("# species: " + " ".join(statedef.species.species) + "\n")
    newlines.append(
        "# bbnames: "
        + " ".join(statedef.backbone_names)
        + (
            ""
            if not statedef.negbackbone_names
            else " ".join(statedef.negbackbone_names)
        )
        + "\n"
    )
    for spec in statedef.species.species:
        for bb in statedef.backbone_names:
            state_file_name = statedef.state_file_name(spec, bb)
            if not os.path.isfile(state_ver_dir(opts) + state_file_name):
                continue
            line = "STATE_VECTOR " + spec + "_" + bb + " " + state_file_name + "\n"
            newlines.append(line)
        if statedef.species.negative_species(spec):
            for bb in statedef.negbackbone_names:
                state_file_name = statedef.state_file_name(spec, bb)
                if not os.path.isfile(state_ver_dir(opts) + state_file_name):
                    continue
                line = "STATE_VECTOR " + spec + "_" + bb + " " + state_file_name + "\n"
                newlines.append(line)
        newlines.append("\n")

    # declare sub expressions for the best energy of each species/bb pair
    for spec in statedef.species.species:
        for bb in statedef.backbone_names:
            state_file_name = statedef.state_file_name(spec, bb)
            if not os.path.isfile(state_ver_dir(opts) + state_file_name):
                continue
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
        if statedef.species.negative_species(spec):
            for bb in statedef.negbackbone_names:
                state_file_name = statedef.state_file_name(spec, bb)
                if not os.path.isfile(state_ver_dir(opts) + state_file_name):
                    continue
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
    # best energy per species vectors vAA, vApA, etc
    # dGbind - one per complex
    newlines.append(
        "# best energies for a single species on each of its available backbones\n"
    )
    for spec in statedef.species.species:
        line = "VECTOR_VARIABLE v" + spec + " = "
        for bb in statedef.backbone_names:
            state_file_name = statedef.state_file_name(spec, bb)
            if not os.path.isfile(state_ver_dir(opts) + state_file_name):
                continue
            line += "best_" + spec + "_" + bb + " "
        if statedef.species.negative_species(spec):
            for bb in statedef.negbackbone_names:
                state_file_name = statedef.state_file_name(spec, bb)
                if not os.path.isfile(state_ver_dir(opts) + state_file_name):
                    continue
                line += "best_" + spec + "_" + bb + " "
        newlines.append(line[:-1] + "\n")
    newlines.append("\n")
    for comp in statedef.species.complexes:
        comps = comp.split("_")
        assert len(comps) == 2
        sep = comps[0] + "_p_" + comps[1]
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

    # now the rest of the fitness function has access to the following variables and should be constructed
    # from them:
    # vMH3_MH4, vMH3_WTH4, vWTH3_MH4
    # vMH3_p_MH4, vMH3_p_WTH4, vWTH3_p_MH4
    # vdGbind_MH3_MH4, vdGBind_MH3_WTH4, vdGBind_WTH3_MH4

    for line in orig_lines:
        if line[0] != "#":
            # append the entity function line to the input daf
            if line[0:7] == "FITNESS" and desdefnames.entfunc:
                efunc_fname = desdefnames.entfunc
                if efunc_fname.find("/") != -1:
                    efunc_fname = efunc_fname.rpartition("/")[2]
                print("appending entity function", efunc_fname)
                newlines.append("ENTITY_FUNCTION entfunc " + efunc_fname + "\n")
                newlines.append("\n")
                line = line[:-1] + " + %f * entfunc\n" % w_ent_func
                print("FITNESS line: ", line)
            newlines.append(line % replace_dict)
        else:
            newlines.append(line)
    return newlines


def dGdiff_bonus_weights(opts):
    if opts.w_dGdiff_bonus_weights_file != "":
        lines = open(opts.w_dGdiff_bonus_weights_file).readlines()
        weights = []
        for line in lines:
            weights.append(line.strip())
        return weights

    default_dGdiff_bonus_weights = ["5", "10", "15", "20", "25", "30"]
    return default_dGdiff_bonus_weights


def entfunc_weights(opts):
    if opts.entfunc_weights_file != "":
        lines = open(opts.entfunc_weights_file).readlines()
        weights = []
        for line in lines:
            weights.append(float(line.strip()))
        return weights
    return [1.0]


def leading_zero_string(num, maxrange):
    nzeros_tot = int(math.floor(math.log10(maxrange))) + 1
    ndigits_num = int(math.floor(math.log10(num))) + 1
    ostr = ""
    for i in range(ndigits_num, nzeros_tot):
        ostr += "0"
    ostr += str(num)
    return ostr


# MAIN FUNCTION
if __name__ == "__main__":
    options = initialize_options()

    design_def = DesDefFnames()
    design_def.read_from_lines(
        options.des_def,
        open(design_def_dir(options) + "definition_files.txt").readlines(),
    )

    if options.pdb_seed_pairs:
        corr = read_positive_species_corr_files(design_def, options)
        print(len(corr))
        seed_seqs = []
        for i in range(len(options.pdb_seed_pairs) / 2):
            seed_seqs.append(
                seed_sequence_from_pdb_pair(
                    corr,
                    (options.pdb_seed_pairs[2 * i], options.pdb_seed_pairs[2 * i + 1]),
                )
            )
        options.seed_sequences = seed_seqs

    # verify that the seed sequences match the design definition's number of entities
    validate_seed_sequences(design_def, options)

    state_def = StateDefinitions()
    state_def.determine_pdbs(options)
    state_def.create_state_file_lists(options)

    exflag_lines = None
    if options.flags_files:
        exflag_lines = {}
        for flags_file in options.flags_files:
            exflag_lines[flags_file] = open(flags_file).readlines()

    dGbonus_weights = dGdiff_bonus_weights(options)
    for w_dGdiff in dGbonus_weights:
        print("Using dGbonus weight ", w_dGdiff)

    ent_func_weights = [1.0]
    if design_def.entfunc != "":
        ent_func_weights = entfunc_weights(options)
    for w_entfunc in ent_func_weights:
        print("Using entity function weight", w_entfunc)

    ok = mkdir(options.job_name)
    if not ok:
        print("ERROR: job directory", options.job_name, "already exists")
        sys.exit(1)
    os.chdir(options.job_name)

    # write a copy of each of the extra flags files to this base directory
    if exflag_lines:
        for fname, exflag_lines in list(exflag_lines.items()):
            print("copying options file:", fname.rpartition("/")[2].strip())
            open(fname.rpartition("/")[2].strip(), "w").writelines(exflag_lines)

    # now iterate across all dGdiff_bonus weights and across all entfunc_weights
    # and create jobs for each combo
    launch_script = []
    gather_script = []
    gather_script.append("mkdir results\n")
    gather_script.append("\n")
    dir_counts = {}
    for w_dGdiff in dGbonus_weights:
        for w_entfunc in ent_func_weights:
            desdir_name = designdir_name(options, design_def, w_dGdiff, w_entfunc)
            if desdir_name not in dir_counts:
                dir_counts[desdir_name] = 1
            else:
                dir_counts[desdir_name] += 1
                desdir_name = (
                    desdir_name
                    + "_"
                    + leading_zero_string(dir_counts[desdir_name], len(dGbonus_weights))
                )
            ok = mkdir(desdir_name)
            os.chdir(desdir_name)
            launch_script.append("cd " + desdir_name + "\n")
            create_job(
                state_def,
                design_def,
                desdir_name,
                options,
                w_dGdiff,
                w_entfunc,
                launch_script,
                gather_script,
            )
            launch_script.append("cd ..\n")
            launch_script.append("\n")
            os.chdir("..")

    gather_script.append("cd results\n")
    gather_script.append(
        "tar -cjf " + options.job_name + ".tar.bz2 *.pdb job_triples.list\n"
    )
    gather_script.append("cd ..\n")

    with open("submit_all_jobs.sh", "w") as fid:
        fid.writelines(launch_script)
    with open("gather_output.sh", "w") as fid:
        fid.writelines(gather_script)
    launch_docking_script = []
    launch_docking_script.append("bash gather_output.sh\n")
    launch_docking_script.append("cd results\n")

    dock_jobs_command = (
        "python3 "
        + base_dir()
        + "pyscripts/dock_jobs_run.py --pdb-triples job_triples.list"
    )
    if options.docking_flags_files:
        dock_jobs_command += " --flag-files"
        for fname in options.docking_flags_files:
            dock_jobs_command += " " + fname
    if options.relax:
        dock_jobs_command += " --relax"
    dock_jobs_command += "\n"

    launch_docking_script.append(dock_jobs_command)

    # now prepare the dock_jobs_view command and append it to the launch_docking script
    dock_jobs_view_job_opts = SubmissionOptions()
    dock_jobs_view_job_opts.scheduler = scheduler_type_for_server()
    dock_jobs_view_job_opts.num_nodes = 1
    dock_jobs_view_job_opts.queue = "debug_queue"
    dock_jobs_view_job_opts.job_name = options.job_name + "_dock_jobs_view"
    dock_jobs_view_job_opts.logfilename = options.job_name + "_djv.log"
    dock_jobs_view_job_opts.submission_script_fname = "djv_submit.sh"

    dock_jobs_view_command_line = (
        "python3 "
        + base_dir()
        + "pyscripts/dock_jobs_view.py -l job_triples.list -o after_docking_dGbind.txt"
    )
    dock_jobs_view_submission_command, dock_jobs_view_command_filename = command_and_submission_script_for_job(
        dock_jobs_view_job_opts, dock_jobs_view_command_line
    )

    if dock_jobs_view_command_filename == None:
        with open("djv_submit.sh", "w") as fid:
            fid.writelines(dock_jobs_view_submission_command)
    dock_jobs_view_command = "python3 /nas/longleaf/home/leaverfa/bin/submit_dependent_script.py dock/dock_submission.log djv_submit.sh"

    # if options.relax :
    #    dock_jobs_view_command += " --relax"
    #    # if options.relax_protocol != "" :
    #    #     dock_jobs_view_command += " --relax-protocol " + options.relax_protocol
    dock_jobs_view_command += "\n"
    launch_docking_script.append(dock_jobs_view_command)

    launch_docking_job_opts = SubmissionOptions()
    launch_docking_job_opts.scheduler = scheduler_type_for_server()
    launch_docking_job_opts.num_nodes = 1
    launch_docking_job_opts.queue = "debug_queue"
    launch_docking_job_opts.job_name = options.job_name + "_launch_docking"
    launch_docking_job_opts.logfilename = options.job_name + "_launch_docking.log"
    launch_docking_job_opts.submission_script_fname = "launch_docking.sh"

    launch_docking_submission_command, launch_docking_command_filename = command_and_submission_script_for_job(
        launch_docking_job_opts, "".join(launch_docking_script)
    )

    if launch_docking_command_filename is None:
        open("launch_docking.sh", "w").writelines(launch_docking_submission_command)
    open("prepare_for_docking.sh", "w").writelines(
        [
            "python3 ~/pyscripts/submit_dependent_script.py msd_submission.log launch_docking.sh\n"
        ]
    )

    # if ( options.skip_docking ) :
    #    # forget the launch_docking script above, we're going directly into the cart_relax_run script
    #    launch_relax = [];
    #    launch_relax.append( "bash gather_output.scr\n" )
    #    launch_relax.append( "cd results\n" )
    #    relax_command = "python " + base_dir() + "pyscripts/cart_relax_run.py --launch"
    #    if options.relax_w_mpi : #use 10 processors or just one?
    #        relax_command += " --mpi"
    #    if options.relax_flags_files :
    #        relax_command += " --flags " + " ".join( options.relax_flags_files )
    #    # if options.relax_protocol :
    #    #     relax_command += " --protocol " + options.relax_protocol
    #    relax_command += "\n"
    #    launch_relax.append( relax_command )
    #    open( "launch_relax.scr", "w" ).writelines( launch_relax )
    #    open( "prepare_for_relax.scr", "w" ).writelines( ["python ~/pyscripts/submit_dependent_script.py msd_submission.log bash launch_relax.scr\n"] )

    command = ""
    for arg in sys.argv:
        command += arg + " "
    command = command[:-1] + "\n"
    open("creation_command.txt", "w").writelines(command)

    if options.launch:
        print("launching")
        print("ls?")
        print(os.getcwd())
        print(os.listdir("."))
        os.system("bash submit_all_jobs.sh")
        os.system("bash prepare_for_docking.sh")
