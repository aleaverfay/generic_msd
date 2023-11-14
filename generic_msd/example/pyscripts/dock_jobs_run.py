import os
import sys
import subprocess
from optparse import OptionParser
from generic_msd.create_lsf_or_slurm_job import *
from generic_msd.opt_holder import OptHolder
from generic_msd.msd_interface_design import PostProcessingOpts
from generic_msd.server_identification import ServerIdentifier
from generic_msd.msd_job_management import rosetta_scripts_jd3_exe
from setup_h3h4_job import base_dir
from replace_design_in_h3h4ctxt import place_design_into_xtal_context

import blargs


def job_skelleton():
    skelleton = [
        ' <Job nstruct="{njobs:d}">\n',
        "  <Input>\n",
        '   <PDB filename="{filename:s}"/>\n',
        "  </Input>\n",
        "  <Output>\n",
        '   <PDB path="{outpath:s}"/>\n',
        "  </Output>\n",
        "  <SecondaryOutput>\n",
        '   <ScoreFile path="{outpath:s}" filename="score.sc"/>\n',
        "  </SecondaryOutput>\n",
        " </Job>\n",
    ]
    return "".join(skelleton)


def job_definition_file_for_pdbs(pdblist, njobs):
    skelly = job_skelleton()
    xml_lines = []
    xml_lines.append("<JobDefinitionFile>\n")
    replacedict = {"njobs": njobs}
    for pdb in pdblist:
        assert pdb[-4:] == ".pdb"
        replacedict["filename"] = pdb
        replacedict["outpath"] = "dock_" + pdb[:-4]
        xml_lines.append(skelly.format(**replacedict))
    xml_lines.append("</JobDefinitionFile>\n")
    return xml_lines


if __name__ == "__main__":

    #parser = initialize_options_parser()
    #(options, args) = parser.parse_args()

    options = OptHolder()
    with blargs.Parser(options) as p:
        PostProcessingOpts.add_options(p)
        p.str("pdb-complexes").required()

    pdb_complexes = open(options.pdb_complexes).readlines()

    si = ServerIdentifier()

    ctxt_pdb = os.path.join(base_dir(), "input_files/starting_structures/1kx5/1KX5_h3h4ctxt.pdb")
    pdbs = []
    new_pdb_complexes = []
    for pdb_tuple in pdb_complexes:
        cols = pdb_tuple.split()
        line = []
        for pdb in cols:
            print(pdb)
            if pdb != "\n":
                newpdbname = pdb.replace(".pdb", "_ctxt.pdb")
                place_design_into_xtal_context(ctxt_pdb, pdb, newpdbname)
                pdbs.append(newpdbname)
                line.append(newpdbname)
        new_pdb_complexes.append(" ".join(line) + "\n")
    filename_parts = os.path.splitext(options.pdb_complexes)
    new_complex_filename = filename_parts[0] + "_in_ctxt" + filename_parts[1]
    with open(new_complex_filename, "w") as fid:
        fid.writelines(new_pdb_complexes)

    # design_directory = sys.argv[1] #exists already, has the designs
    dock_directory = "dock"
    os.mkdir(dock_directory)
    os.chdir(dock_directory)
    for pdb in pdbs:
        os.system("ln -s ../" + pdb + " .")
        os.mkdir("dock_" + pdb[:-4])

    jobdef_xml = job_definition_file_for_pdbs(pdbs, 20)
    with open("jobdef.xml", "w") as fid:
        fid.writelines(jobdef_xml)

    # one head node to dispatch jobs
    # one archive node for ever 50 compute nodes
    # try and run twenty docking jobs
    narchives = (len(pdbs) - 1) // 50 + 1
    ncpus = 2 * len(pdbs) + 1 + narchives
    ncpus = max(ncpus, 45) # use at least 45 CPUs -- min number for dogwood

    if ncpus > 200:
        # eh, just round up to get into the 2112 queue.        
        ncpus = 13 * 44
        narchives = (ncpus-1) // 50 + 1

    open("creation_command.py", "w").writelines([" ".join(sys.argv) + "\n"])

    dock_submit_opts = SubmissionOptions()
    dock_submit_opts.scheduler = scheduler_type_for_server(si)
    dock_submit_opts.num_nodes = ncpus
    dock_submit_opts.queue = "auto"
    dock_submit_opts.job_name = "docking"
    dock_submit_opts.logfilename = "docking.log"
    dock_submit_opts.submission_script_fname = "docking.sh"
    dock_submit_opts.mpi_job = True

    dock_command = "".join(
        [
            rosetta_scripts_jd3_exe(si),
            " -parser:protocol ",
            os.path.join(base_dir(), "input_files/post_processing/relax/fastrelax_in_nuc_ctxt_and_intana.xml"),
            " -job_definition_file jobdef.xml",
            " -n_archive_nodes ",
            str(narchives),
            " -mpi_tracer_to_file docking.log",
        ]
    )

    submission_command = command_and_submission_script_for_job(
        dock_submit_opts, dock_command
    )

    print(submission_command)
    os.system(submission_command + " >> dock_submission.log")
