import itertools
import amino_acids as aas
import pdb_structure
import blargs
from compare_h3h4_sequence import h3h4_pdb, find_mutations_for_chain
from pyrosetta import *
import pyrosetta.rosetta.core.select.residue_selector as rs
import pyrosetta.rosetta.core.pack.task.operation as taskops

def revert_mutations_for_struct(pdb):
    orig_struct = pdb_structure.pdbstructure_from_file(h3h4_pdb())
    design_struct = pdb_structure.pdbstructure_from_file(pdb)
    chAmuts = find_mutations_for_chain(orig_struct, design_struct, "A")
    chBmuts = find_mutations_for_chain(orig_struct, design_struct, "B")

    # ok, so?
    # let's create a packer task and repack in the neighborhood of
    # the mutated residues
    des_pose = pose_from_file(pdb)
    des_resA = [des_pose.pdb_info().pdb2pose("A",int(mut[1])) for mut in chAmuts]
    des_resB = [des_pose.pdb_info().pdb2pose("B",int(mut[1])) for mut in chBmuts]

    sfxn = get_score_function()
    sfxn(des_pose)
    resinds = pyrosetta.rosetta.utility.vector1_unsigned_long()
    for ind in itertools.chain(des_resA, des_resB):
        resinds.append(ind)
    des_res_selector = rs.ResidueIndexSelector(resinds)
    neighbors = rs.NeighborhoodResidueSelector(des_res_selector, 10.0, False)
    neighbors_or_mutated = rs.OrResidueSelector(neighbors, des_res_selector)

    not_neighbors = rs.NotResidueSelector(neighbors_or_mutated)

    disable_packing = taskops.PreventRepackingRLT()
    disable_neighbors_op = taskops.OperateOnResidueSubset(disable_packing, not_neighbors)
    disable_design = taskops.RestrictToRepackingRLT()
    repack_neighbors_op = taskops.OperateOnResidueSubset(disable_design, neighbors)
    ex12_neighbors = taskops.ExtraRotamersGenericRLT()
    ex12_neighbors.ex1(True)
    ex12_neighbors.ex2(True)
    ex12_neighbors_op = taskops.OperateOnResidueSubset(ex12_neighbors, neighbors)

    task = pyrosetta.rosetta.core.pack.task.TaskFactory.create_packer_task(des_pose)
    disable_neighbors_op.apply(des_pose, task)
    repack_neighbors_op.apply(des_pose, task)
    ex12_neighbors_op.apply(des_pose, task)

    allmuts = chAmuts + chBmuts
    allres  = des_resA + des_resB
    print("allmuts", allmuts)
    print("allres", allres)

    for i, mut in enumerate(allmuts):
        keep_aas = pyrosetta.rosetta.utility.vector1_bool(20)
        for j in range(1,21):
            keep_aas[j] = False
        onelet = aas.longer_names[mut[0]]
        pos = pyrosetta.rosetta.core.chemical.aa_from_one_or_three(mut[0])
        keep_aas[pos] = True
        task.nonconst_residue_task(allres[i]).restrict_absent_canonical_aas(keep_aas)

    task.or_include_current(True)

    pyrosetta.rosetta.core.pack.pack_rotamers(des_pose, sfxn, task)
    outname = pdb.split(".pdb")[0] + "_wt.pdb"
    des_pose.dump_pdb(outname)

    return outname



def write_chain_as_pdb(complex_struct, chain, fname, suffix):
    onechain = pdb_structure.PDBStructure()
    onechain.add_chain(complex_struct.chainmap[chain])
    onechain_fname = fname.partition(".pdb")[0] + "_" + suffix + ".pdb"
    with open(onechain_fname,"w") as fid:
       fid.writelines(onechain.pdb_lines())

if __name__ == "__main__":
    with blargs.Parser(locals()) as p:
        p.multiword("input_pdbs").cast(lambda x: x.split()).required()

    init()

    for pdb in input_pdbs:
        wt_pdb = revert_mutations_for_struct(pdb)
        complex_struct = pdb_structure.pdbstructure_from_file(wt_pdb)
        write_chain_as_pdb(complex_struct, "A", wt_pdb, "h3")
        write_chain_as_pdb(complex_struct, "B", wt_pdb, "h4")
