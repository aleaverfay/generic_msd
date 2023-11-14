import os
import sys
import math
from read_IntAnalyzer_log import *
from optparse import OptionParser
from compare_h3h4_sequence import compare_h3h4_seq_against_design
from generic_msd.msd_interface_design import PostProcessingOpts
from generic_msd.opt_holder import OptHolder
import blargs


class IAResult:
    def __init__(self):
        self.fname = ""
        self.job_index = 0
        self.dGbind = 0
        self.dSASA = 0
        self.eDens = 0
        self.pack = 0
        self.uns = 0
        self.total_score = 0

def scorefile_column_ranges_from_lines( line ):
    """Compute the column ranges
    
    The column ranges are taken as the beginning of the space characters
    ahead of the first non-whitespace character to the last character
    of the column header. The last column, however, is considered to
    go all the way to the end of the line"""
    
    col_ranges = []
    found_non_ws = False
    start = 0
    for i in range(len(line)):
        c = line[i]
        if " " == c:
            if found_non_ws:
                col_ranges.append((start,i))
                found_non_ws = False
                start = i+1
        elif "\n" == c and found_non_ws and i == len(line)-1:
            col_ranges.append((start,-1))
        else:
            found_non_ws = True
    if line[-1] != "\n":
        col_ranges.append((start,-1))
    else:
        start,stop = col_ranges.pop()
        col_ranges.append((start,-1))
    print("col ranges:", col_ranges)
    return col_ranges


def read_iaresults_from_score_file(sc_fname, ia_prefix=""):
    with open(sc_fname) as fid:
        lines = fid.readlines()
    score_line = lines[1]
    assert score_line.startswith("SCORE:")
    col_ranges = scorefile_column_ranges_from_lines(lines[1])
    headers = {lines[1][x[0]:x[1]].strip(): x for x in col_ranges}
    results = []
    for i in range(2,len(lines)):

        if lines[i].startswith("SEQUENCE"):
            continue

        # unpack the score line as pairs of descriptions and values
        line_params = {h: lines[i][r[0]:r[1]] for h,r in headers.items()}
        if line_params["total_score"] == "total_score":
            continue

        res = IAResult()
        res.fname = line_params["description"]
        print("result:", res.fname)
        pref = "" if ia_prefix == "" else (ia_prefix + "_")
        res.job_index = line_params["description"][:-4]
        res.dGbind = float(line_params[pref + "dG_separated"])
        res.dSASA = float(line_params[pref + "dSASA_int"])
        res.uns = float(line_params[pref + "delta_unsatHbonds"])
        res.total_score = float(line_params["total_score"])
        results.append(res)
    return results

def avg_binding_energy_from_top10(results1, results2):
    top10_by_total = sorted(range(len(results1), key=lambda ind: results1[ind].total_score)[0:10]

    return sum((results1[ind].dGbind + results2[ind].dGbind) for ind in top10_by_total) / 10


if __name__ == "__main__":
    #parser = initialize_options_parser()
    #(options, args) = parser.parse_args()
    options = OptHolder()
    with blargs.Parser(options) as p:
        PostProcessingOpts.add_options(p)
        p.str("pdb-triples").shorthand("l").required()
        p.str("output-file").shorthand("o").required()
        p.flag("repeat-run").shorthand("r")
        p.flag("no_delete").shorthand("n")

    filename_parts = os.path.splitext(options.pdb_triples) # "triples" = 8-tuples
    new_complex_filename = filename_parts[0] + "_in_ctxt" + filename_parts[1]
    triples_list = open(new_complex_filename).readlines()

    if not options.repeat_run:
        os.mkdir("best_docked")
    os.chdir("dock")
    #os.system("rm out.*") # cleanup delayed for debugging
    lines = []
    for triple in triples_list:
        if triple == "\n":
            continue
        cols = triple.split()
        ia_dat = []
        for col in cols:
            dirname = "dock_" + col[0:-4]
            os.chdir(dirname)
            #avgtop10 = find_avgtop10_dGbind_from_fname(dirname + ".logIA")
            if not os.path.isfile("score.sc"):
                os.system("cat score.sc.* > score.sc")
            iaresults_h3h4_vs_nuc = read_iaresults_from_score_file("score.sc", "h3h4_vs_nuc")
            iaresults_h3_vs_h4 = read_iaresults_from_score_file("score.sc", "h3_vs_h4")
            avgtop10 = avg_binding_energy_from_top10(iaresults_h3h4_vs_nuc, iaresults_h3_vs_h4)
            best_result = min(iaresults, key=lambda x: x.total_score)
            print(dirname, avgtop10, best_result)
            if not options.repeat_run:
                os.system(
                    "cp "
                    + best_result.fname
                    + ".pdb"
                    + " ../../best_docked/"
                    + dirname
                    + ".pdb"
                )
            ia_dat.append(avgtop10)
            os.chdir("..")
        worst_on_target_dGbind = max(ia_dat[0], ia_dat[1])
        line = "%s, %f, %f, %f, %f, %f, %f, %f, %f ( %f, %f ) ( %f, %f, %f, %f ) " % (
            cols[0],
            ia_dat[0],  # on target mutant pairs
            ia_dat[1],
            ia_dat[2],  # off target mutant / mutant pairs
            ia_dat[3],
            ia_dat[4],  # off target mutant / wt pairs
            ia_dat[5],
            ia_dat[6],  # off target wt / mutant pairs
            ia_dat[7],
            worst_on_target_dGbind - ia_dat[2], # dGbind on-target/off-target
            worst_on_target_dGbind - ia_dat[3],
            ia_dat[ 0 ] - ia_dat[4], ia_dat[1] - ia_dat[5], # dGbind mut/wt
            ia_dat[ 0 ] - ia_dat[6], ia_dat[1] - ia_dat[7] ) # dGbind wt/mut
        line += " AH3_AH4: " + compare_h3h4_seq_against_design( cols[0], True )
        line += " BH3_BH4: " + compare_h3h4_seq_against_design( cols[1], True )

        line += "\n"
        lines.append(line)

    os.chdir("..")
    os.chdir("best_docked")
    if options.output_file:
        open(options.output_file, "w").writelines(lines)
    else:
        for line in lines:
            print(line, end=" ")
    os.chdir("..")
    if not options.no_delete:
        print("TEMP NOT Deleting docked PDBs")
        # TEMP! os.system("find dock | grep -e '.pdb$' | xargs rm")
    else:
        print("Not deleting docked PDBs")

    # if options.relax :
    #     os.chdir( "best_docked" )
    #     relax_command = "python " + setup_jobs.base_dir() + "/pyscripts/cart_relax_run.py --launch"
    #     # if relax_protocol != "" :
    #     #     relax_command += " --protocol " + relax_protocol
    #     os.system( relax_command )
