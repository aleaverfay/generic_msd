import blargs
from compare_h3h4_sequence import h3h4_pdb
import pdb_structure
import amino_acids as aas
import sys

from yaml import load, dump
try:
    from yaml import CLoader as Loader, CDumper as Dumper
except ImportError:
    from yaml import Loader, Dumper



def write_corr_files(design_positions):
    # create the .corr files
    # There are two "species" -- AH3/AH4 and BH3/BH4
    # The species are interleaved with odd positions
    # representing the A species and even positions
    # representing the B species.
    # "design positions" is a list of tuples, with
    # tup[0] == the species (either A or B)
    # tup[1] == the chain, either A or B
    # tup[2] == the residue index

    def write_corr_file( des_positions, filter, fname_noext ):
        lines = []
        for i, despos in enumerate(des_positions):
            if filter(despos):
                lines.append( "%d %d %s\n" % (i + 1, despos[2], despos[1]))
        open(fname_noext + ".corr", "w").writelines(lines)

    write_corr_file(design_positions, lambda x: x[0] == "A", "AH3_AH4")
    write_corr_file(design_positions, lambda x: x[0] == "B", "BH3_BH4")
    write_corr_file(design_positions, lambda x: (x[0] == "A" and x[1] == "A") or (x[0] == "B" and x[1] == "B"), "AH3_BH4")
    write_corr_file(design_positions, lambda x: (x[0] == "B" and x[1] == "A") or (x[0] == "A" and x[1] == "B"), "BH3_AH4")
    write_corr_file(design_positions, lambda x: x[0] == "A" and x[1] == "A", "AH3")
    write_corr_file(design_positions, lambda x: x[0] == "B" and x[1] == "A", "BH3")
    write_corr_file(design_positions, lambda x: x[0] == "A" and x[1] == "B", "AH4")
    write_corr_file(design_positions, lambda x: x[0] == "B" and x[1] == "B", "BH4")

def write_entity_resfile(design_positions, fixed_pos, forbidden_aas, allowed_aas):
    # ok, now the entitiy resfile, respecting the fixed positions and the forbidden mutations as indicated
    allowed_map = {}
    # for i in xrange(len(design_positions)) :
    #    allowed_map.append( list (aas.amino_acids) ) # clone the amino_acid list
    #    allowed_map[ i ].remove( "C" ) # don't design cysteine

    for spec, ch, res, mut in fixed_pos:
        allowed_map[design_positions.index((spec, ch, res))] = [mut]
    for spec, ch, res, allowed in allowed_aas:
        allowed_map[design_positions.index((spec, ch, res))] = allowed
    for spec, ch, res, aalist in forbidden_aas:
        entity_id = design_positions.index((spec, ch, res))
        not_forbidden = list(aas.amino_acids)  # clone the amino_acid list
        not_forbidden.remove("C")  # don't design cysteine
        for mut in aalist:
            if mut in not_forbidden:
                not_forbidden.remove(mut)
        allowed_map[entity_id] = not_forbidden

    lines = []
    lines.append("%d\n" % (len(design_positions)))
    lines.append("ALLAAxc EX 1 EX 2\n")
    lines.append("start\n")

    for entity_id in allowed_map:
        lines.append(
            "%d A PIKAA %s EX 1 EX 2\n"
            % (entity_id + 1, "".join(allowed_map[entity_id]))
        )
    open("entity.resfile", "w").writelines(lines)


def write_secondary_resfiles(
    design_positions, secondary_positions, optimize_water_placement
):
    def write_secresfile(despos_filter, secres_filter, fname_noex):
        lines = []
        lines.append("NATRO\n")
        lines.append("start\n")
        for despos in design_positions:
            if despos_filter(despos):
                lines.append("%d %s NATAA USE_INPUT_SC EX 1 EX 2\n" % (despos[2], despos[1]))
        for secpos in secondary_positions:
            if secres_filter(secpos):
                lines.append("%d %s NATAA USE_INPUT_SC EX 1\n" % (secpos[1], secpos[0]))
        # if optimize_water_placement:
        #     lines.append("* W NC WAT NC VRT1\n")
        open(fname_noex + ".2resfile", "w").writelines(lines)

    def take_all(_):
        return True
    def take_none(_):
        return False
    def take_chA(x):
        return x[-2] == "A"
    def take_chB(x):
        return x[-2] == "B"

    write_secresfile(take_none, take_all, "design_both")
    write_secresfile(take_chB, take_all, "desA_wtB")
    write_secresfile(take_chA, take_all, "wtA_desB")
    write_secresfile(take_none, take_chA, "desA")
    write_secresfile(take_none, take_chB, "desB")

def write_design_constraints(design_positions, reference_pdbs_A, reference_pdbs_B):
    # design constraints
    if nmut_limit:
        reference_structs_A = []
        reference_structs_B = []
        if not reference_pdbs_A:
            reference_pdbs_A.append(h3h4_pdb())
        if not reference_pdbs_B:
            reference_pdbs_B.append(h3h4_pdb())

        for refpdb in reference_pdbs_A:
            print("AH3 AH4 reference pdb", refpdb)
            reference_structs_A.append(
                pdb_structure.pdbstructure_from_file(refpdb)
            )
        reference_structs_A = [pdb_structure.pdbstructure_from_file(refpdb) for refpdb in reference_pdbs_A]
        reference_structs_B = [pdb_structure.pdbstructure_from_file(refpdb) for refpdb in reference_pdbs_B]

        lines = []
        for i, refstruct in enumerate(reference_structs_A):
            # top line -- state the amino-acid string for this set of designs
            lines.append(
                "# original AH3_AH4 sequence for this design definition in "
                + reference_pdbs_A[i].rpartition("/")[2]
                + " is given by\n"
            )
            line = "# "
            for despos in design_positions:
                if despos[0] == "A":
                    origaa = aas.longer_names[
                        refstruct.residue(despos[1], str(despos[2])).resname
                    ]
                    line += origaa
            lines.append(line + "\n")

        for i, refstruct in enumerate(reference_structs_B):
            # top line -- state the amino-acid string for this set of designs
            lines.append(
                "# original BH3_BH4 sequence for this design definition in "
                + reference_pdbs_B[i].rpartition("/")[2]
                + " is given by\n"
            )
            line = "# "
            for despos in design_positions:
                if despos[0] == "B":
                    origaa = aas.longer_names[
                        refstruct.residue(despos[1], str(despos[2])).resname
                    ]
                    line += origaa
            lines.append(line + "\n")

        lines.append("\n")
        # declare the eeXnat variables
        for i, idespos in enumerate(design_positions):
            origaa = set([])
            refstructs = reference_structs_A if i % 2 == 0 else reference_structs_B
            for refstruct in refstructs:
                origaa.add(
                    aas.longer_names[
                        refstruct.residue(idespos[1], str(idespos[2])).resname
                    ]
                )
            lines.append(
                "# entity element ee_%d corresponds to %s%s%d: %s \n"
                % (i + 1, idespos[0], idespos[1], idespos[2], ",".join(origaa))
            )
            lines.append(
                "SET_CONDITION ee%dnat = ee_%d in { %s }\n"
                % (i + 1, i + 1, ", ".join(origaa))
            )
            lines.append("\n")

        # count the number of mutations for the designed H3/H3 pairs
        line = "SUB_EXPRESSION nnatA = "
        line +=  " + ".join(["ee%dnat" % (i + 1) for i, dp in enumerate(design_positions) if dp[0] == "A"])
        line += "\n"
        lines.append(line)
        line = "SUB_EXPRESSION nnatB = "
        line +=  " + ".join(["ee%dnat" % (i + 1) for i, dp in enumerate(design_positions) if dp[0] == "B"])
        line += "\n"
        lines.append(line)
        line = "SUB_EXPRESSION nmutA = %d - nnatA\n" % (len(design_positions) / 2)
        lines.append(line)
        lines.append("\n")
        line = "SUB_EXPRESSION nmutB = %d - nnatB\n" % (len(design_positions) / 2)
        lines.append(line)
        lines.append("\n")

        line = "SUB_EXPRESSION mut_penaltyA = ite( gt( nmutA, %d ), nmutA - %d, 0 )\n" % (
            nmut_limit,
            nmut_limit,
        )
        lines.append(line)
        lines.append("\n")

        line = "SUB_EXPRESSION mut_penaltyB = ite( gt( nmutB, %d ), nmutB - %d, 0 )\n" % (
            nmut_limit,
            nmut_limit,
        )
        lines.append(line)
        lines.append("\n")

        lines.append("SCORE %f * ( mut_penaltyA + mut_penaltyB)\n" % mut_penalty_scale)
        open("h3h4_native.entfunc", "w").writelines(lines)


def write_definition_file():

    # # now, the definition_files.txt file
    # lines = []
    # lines.append("AH3_AH4    corr= AH3_AH4.corr   2res= design_both.2resfile\n")
    # lines.append("BH3_BH4    corr= BH3_BH4.corr   2res= design_both.2resfile\n")
    # lines.append("AH3_BH4    corr= AH3_BH4.corr   2res= design_both.2resfile\n")
    # lines.append("BH3_AH4    corr= BH3_AH4.corr   2res= design_both.2resfile\n")
    # lines.append("AH3_WTH4   corr= AH3.corr       2res= desA_wtB.2resfile\n")
    # lines.append("BH3_WTH4   corr= BH3.corr       2res= desA_wtB.2resfile\n")
    # lines.append("WTH3_AH4   corr= AH4.corr       2res= wtA_desB.2resfile\n")
    # lines.append("WTH3_BH4   corr= BH4.corr       2res= wtA_desB.2resfile\n")
    # lines.append("AH3        corr= AH3.corr       2res= desA.2resfile\n")
    # lines.append("BH3        corr= BH3.corr       2res= desA.2resfile\n")
    # lines.append("AH4        corr= AH4.corr       2res= desB.2resfile\n")
    # lines.append("BH4        corr= BH4.corr       2res= desB.2resfile\n")
    # 
    # if nmut_limit:
    #     lines.append("entfunc h3h4_native.entfunc\n")

    yaml_contents = {}
    yaml_contents["species"] = [
        {"name":"AH3_AH4",  "corr": "AH3_AH4.corr", "2res": "design_both.2resfile"},
        {"name":"BH3_BH4",  "corr": "BH3_BH4.corr", "2res": "design_both.2resfile"},
        {"name":"AH3_BH4",  "corr": "AH3_BH4.corr", "2res": "design_both.2resfile"},
        {"name":"BH3_AH4",  "corr": "BH3_AH4.corr", "2res": "design_both.2resfile"},
        {"name":"AH3_WTH4", "corr": "AH3.corr",     "2res": "desA_wtB.2resfile"},
        {"name":"BH3_WTH4", "corr": "BH3.corr",     "2res": "desA_wtB.2resfile"},
        {"name":"WTH3_AH4", "corr": "AH4.corr",     "2res": "wtA_desB.2resfile"},
        {"name":"WTH3_BH4", "corr": "BH4.corr",     "2res": "wtA_desB.2resfile"},
        {"name":"AH3m",     "corr": "AH3.corr",     "2res": "desA.2resfile"},
        {"name":"BH3m",     "corr": "BH3.corr",     "2res": "desA.2resfile"},
        {"name":"AH4m",     "corr": "AH4.corr",     "2res": "desB.2resfile"},
        {"name":"BH4m",     "corr": "BH4.corr",     "2res": "desB.2resfile"},
    ]
    if nmut_limit:
        yaml_contents["entfunc"] = "h3h4_native.entfunc"

    output = dump(yaml_contents, Dumper=Dumper)
    open("definition_files.yaml", "w").writelines(output)


def write_script_arguments_to_file(outfile_name="creation_command.txt"):
    open(outfile_name, "w").writelines([" ".join(sys.argv) + "\n"])


def error_check_aalist_for_species_and_position(
    aas_for_positions, flagname, design_positions
):
    # error checking: allowed_aas
    for spec, ch, pos, aalist in aas_for_positions:
        if (spec, ch, pos) not in design_positions:
            print(
                "ERROR, could not find the position, "
                + spec
                + ch
                + str(pos)
                + ", for the flag "
                + flagname
                + ", in the list of designable positions"
            )
            assert (spec, ch, pos) in design_positions
        for allowed_aa in aalist:
            if allowed_aa not in aas.amino_acids:
                print(
                    "ERROR, the aa "
                    + allowed_aa
                    + " given with the flag "
                    + flagname
                    + " with values "
                    + spec
                    + ch
                    + str(pos)
                    + " is not a valid amino acid"
                )
                assert allowed_aa in aas.amino_acids


def error_check_exclusive(allowed_aas, forbidden_aas):
    allowed_pos = []
    for spec, ch, pos, aalist in allowed_aas:
        allowed_pos.append((spec, ch, pos))
    for spec, ch, pos, aalist in forbidden_aas:
        if (spec, ch, pos) in allowed_aas:
            print(
                "ERROR: a position is listed with both allowed and forbidden aas.  Only one flag per position may be used"
            )
            assert (spec, ch, pos) not in allowed_aas


if __name__ == "__main__":
    with blargs.Parser(locals()) as p:
        p.require_one(
            p.multiword("design_positions")
            .shorthand("d")
            .cast(lambda x: [(val[0], int(val[1:])) for val in x.split()]),
            p.str("design_position_file"),
        )
        p.str("all_repackable_positions_file").shorthand("r").required()
        # e.g. --fixed_pos AB116R = species, chain, residue ID, fixed amino acid
        p.multiword("fixed_pos").default([]).cast(
            lambda x: [(val[0], val[1], int(val[2:-1]), val[-1]) for val in x]
        )
        p.float("mut_penalty_scale").default(2.0).requires(p.int("nmut_limit"))
        # e.g --forbidden_aas BD114_HQ = species, chain, residue ID, "_", forbidden amino acids
        p.multiword("forbidden_aas").cast(
            lambda x: [
                (val.split("_")[0][0], val.split("_")[0][1], int(val.split("_")[0][2:]), val.split("_")[1])
                for val in x.split()
            ]
        )
        # e.g. --allowed_aas BB188_FW = species, chain, residue ID, "_", allowed amino acids
        p.multiword("allowed_aas").cast(
            lambda x: [
                (val.split("_")[0][0], val.split("_")[0][1], int(val.split("_")[0][2:]), val.split("_")[1])
                for val in x.split()
            ]
        )
        p.multiword("A_reference_pdbs").default("").cast(
            lambda x: [val for val in x.split()]
        )
        p.multiword("B_reference_pdbs").default("").cast(
            lambda x: [val for val in x.split()]
        )
        p.flag(
            "optimize_water_placement"
        )  # assume water chains are "W" for the complex, and "W" and "V" for the separated chains

    if design_position_file:
        design_positions = []
        dplines = open(design_position_file).readlines()
        for dpline in dplines:
            cols = dpline.split()
            for col in cols[1:]:
                design_positions.append((cols[0], int(col)))

    if not fixed_pos:
        fixed_pos = []
    if not forbidden_aas:
        forbidden_aas = []
    if not allowed_aas:
        allowed_aas = []

    for spec, ch, pos, mut in fixed_pos:
        print("fixed pos", spec, ch, pos, mut)
        assert ch == "A" or ch == "B"

    # error checking: design positions
    design_positions_sorted = sorted(
        design_positions, key=lambda x: (0 if x[0] == "A" else 1000) + x[1]
    )
    design_positions = []
    for pos in design_positions_sorted:
        if pos[0] != "A" and pos[0] != "B":
            print(
                "Design position",
                pos[0],
                pos[1],
                "requires the chain to be either 'A' or 'B'",
            )
            sys.exit(1)
        print("designable positions:", pos[0], pos[1])
        design_positions.append(
            ("A", pos[0], pos[1])
        )
        design_positions.append(
            ("B", pos[0], pos[1])
        )

    # error checking: fixed positions
    for spec, ch, pos, mut in fixed_pos:
        if (spec, ch, pos) not in design_positions:
            print(
                "ERROR, could not find the fixed position, "
                + spec
                + ch
                + str(pos)
                + ", in the list of designable positions"
            )
            assert (spec, ch, pos) in design_positions

    error_check_aalist_for_species_and_position(
        allowed_aas, "allowed_aas", design_positions
    )
    error_check_aalist_for_species_and_position(
        forbidden_aas, "forbidden_aas", design_positions
    )
    error_check_exclusive(allowed_aas, forbidden_aas)

    secondary_positions = []
    lines = open(all_repackable_positions_file).readlines()
    for line in lines:
        cols = line.split()
        for col in cols[1:]:
            val = ("A", cols[0], int(col))
            if val not in design_positions:
                secondary_positions.append((val[1], val[2]))
    for secpos in secondary_positions:
        print("secondary position:", secpos[0], secpos[1])

    write_corr_files(design_positions)
    write_entity_resfile(design_positions, fixed_pos, forbidden_aas, allowed_aas)
    write_secondary_resfiles(
        design_positions, secondary_positions, optimize_water_placement
    )
    write_design_constraints(design_positions, A_reference_pdbs, B_reference_pdbs)
    write_definition_file()
    write_script_arguments_to_file()
