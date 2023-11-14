from pdb_structure import pdbstructure_from_file
from setup_h3h4_job import base_dir
import blargs
import amino_acids


def h3h4_pdb():
    # print base_dir() + "input_files/starting_structures/1kx5/1KX5_chAB.pdb"
    return base_dir() + "input_files/starting_structures/1kx5/1KX5_chAB.pdb"


def compare_h3h4_seq_against_design(design, compact, chAonly=False, chBonly=False):
    return compare_h3h4_original_against_design(
        h3h4_pdb(), design, compact, chAonly, chBonly
    )

def find_mutations_for_chain(orig_struct, design_struct, chain):
    orig_ch = orig_struct.chainmap[chain]
    des_ch = design_struct.chainmap[chain]
    muts = []
    for resid, orig_res in orig_ch.resmap.items():
        if resid not in des_ch.resmap:
            continue
        des_res = des_ch.resmap[resid]
        if (orig_res.resname != des_res.resname):
            mut = (orig_res.resname, resid, des_res.resname)
            muts.append(mut)
    return muts

def compare_h3h4_original_against_design(
    orig, design, compact, chAonly=False, chBonly=False
):
    orig_struct = pdbstructure_from_file(orig)
    design_struct = pdbstructure_from_file(design)
    chAmuts = find_mutations_for_chain(orig_struct, design_struct, "A")
    chBmuts = find_mutations_for_chain(orig_struct, design_struct, "B")
    # for ch in orig_struct.chainmap:
    #     if ch != "A" and ch != "B":
    #         continue  # skip water chains
    #     assert ch in design_struct.chainmap
    #     for res in orig_struct.chainmap[ch].resmap:
    #         if res not in design_struct.chainmap[ch].resmap:
    #             continue
    #         if (
    #             orig_struct.residue(ch, res).resname
    #             != design_struct.residue(ch, res).resname
    #         ):
    #             mut = (
    #                 orig_struct.residue(ch, res).resname,
    #                 res,
    #                 design_struct.residue(ch, res).resname,
    #             )
    #             if ch == "A":
    #                 chAmuts.append(mut)
    #             else:
    #                 chBmuts.append(mut)
    chAmuts = sorted(chAmuts, key=lambda mut: int(mut[1]))
    chBmuts = sorted(chBmuts, key=lambda mut: int(mut[1]))
    line = ""
    if compact:
        if not chBonly:
            line += "chA "
            for mut in chAmuts:
                line += "%c%s%c " % (
                    amino_acids.longer_names[mut[0]],
                    mut[1],
                    amino_acids.longer_names[mut[2]],
                )
        if not chAonly:
            line += "chB "
            for mut in chBmuts:
                line += "%c%s%c " % (
                    amino_acids.longer_names[mut[0]],
                    mut[1],
                    amino_acids.longer_names[mut[2]],
                )
        line = line[:-1]
    else:
        if not chBonly:
            for mut in chAmuts:
                line += "chA: %s %s %s\n" % mut
        if not chBonly:
            for mut in chBmuts:
                line += "chB: %s %s %s\n" % mut
    return line


if __name__ == "__main__":
    with blargs.Parser(locals()) as p:
        p.str("design").shorthand("d").required()
        p.str("reference_pdb").default(h3h4_pdb())
        p.flag("compact").shorthand("c")
        chao = p.flag("chH3only").shorthand("a")
        p.flag("chH4only").shorthand("b").conflicts(chao)

    line = compare_h3h4_original_against_design(
        reference_pdb, design, compact, chH3only, chH4only
    )
    print(line, end=" ")
