import pdb_structure
import vector3d
import blargs
import copy


def translation_vector_for_h3h4(pdb):
    ca1 = pdb.residue("A", "96").atom(" CA ").xyz
    ca2 = pdb.residue("B", "58").atom(" CA ").xyz

    trans = (ca2 - ca1).normalize().scale(50)
    return trans


def dup_h2o_and_sep_h3h4(input_pdbfilename, output_pdbfilename):
    # 1. take the input pdb structure, duplicate its waters (in chain W) to chain V
    # 2. send chain L and chain V off by 30 angstroms

    pdb = pdb_structure.pdbstructure_from_file(input_pdbfilename)
    chW = pdb.chainmap["W"]
    chV = pdb_structure.Chain()
    chV.chain_name = "V"

    for res in chW.residues:
        if res.resname != "WAT" and res.resname != "HOH":
            continue
        res.chain = None  # temporary
        copyres = copy.deepcopy(res)
        res.chain = chW
        chV.add_residue(copyres)
    pdb.add_chain(chV)

    trans = translation_vector_for_fab(pdb)

    for chid in ["B", "V"]:
        for res in pdb.chainmap[chid].residues:
            for atom in res.atoms:
                atom.xyz += trans
    open(output_pdbfilename, "w").writelines(pdb.pdb_lines())


if __name__ == "__main__":
    with blargs.Parser(locals()) as p:
        p.str("input").required()
        p.str("output").required()

    dup_h2o_and_sep_fab(input, output)
