import pdb_structure
import blargs

def place_design_into_xtal_context(ctxt_pdb, design, output):
    ctxt = pdb_structure.pdbstructure_from_file(ctxt_pdb)
    des = pdb_structure.pdbstructure_from_file(design)

    newstruct = pdb_structure.PDBStructure()
    for chain in ctxt.chains:
        ch = chain
        if chain.chain_name == "A" or chain.chain_name == "B":
            ch = des.chainmap[chain.chain_name]
        newstruct.add_chain(ch)
    with open(output,"w") as fid:
        fid.writelines(newstruct.pdb_lines())

if __name__ == "__main__":
    with blargs.Parser(locals()) as p:
        p.str("ctxt_pdb").required()
        p.str("design").required()
        p.str("output").required()
    place_design_into_xtal_context(ctxt_pdb, design, output)
