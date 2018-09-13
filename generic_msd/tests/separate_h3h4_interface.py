import blargs
import pdb_structure
from dup_h2o_and_separate import translation_vector_for_h3h4


def separate_h3h4_int(input_pdbname, output_pdbname):
    pose = pdb_structure.pdbstructure_from_file(input_pdbname)
    trans = translation_vector_for_h3h4(pose)
    for res in pose.chainmap["B"].residues:
        for atom in res.atoms:
            atom.xyz += trans
    open(output_pdbname, "w").writelines(pose.pdb_lines())


if __name__ == "__main__":
    with blargs.Parser(locals()) as p:
        p.str("input").required()
        p.str("output").required()

    separate_h3h4_int(input, output)
