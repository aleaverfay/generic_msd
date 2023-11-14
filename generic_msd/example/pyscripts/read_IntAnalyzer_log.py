import os
import sys


class IA_Result:
    def __init__(self):
        self.fname = ""
        self.job_index = 0
        self.dGbind = 0
        self.dSASA = 0
        self.eDens = 0
        self.pack = 0
        self.uns = 0
        # print "instantiating new IA_Result"

    def read_resultline(self, line):
        cols = line.split()
        first = False
        if self.fname == "":
            self.fname = cols[0][:-1] + ".pdb"
            # print "setting fname:", self.fname, "cols[0]", cols[0], cols[0][:-1]
            first = True
        if not first and cols[0][:-1] != self.fname[:-4]:
            return False
        if (
            len(cols) >= 6
            and " ".join(cols[1:5]) == "SEPARATED INTERFACE ENERGY DIFFERENCE:"
        ):
            # print "setting dGbind"
            self.dGbind = float(cols[5])
        if len(cols) >= 5 and " ".join(cols[1:4]) == "INTERFACE DELTA SASA:":
            # print "setting dSASA"
            self.dSASA = float(cols[4])
        if (
            len(cols) >= 7
            and " ".join(cols[1:6])
            == "INTERFACE DELTA SASA/SEPARATED INTERFACE ENERGY:"
        ):
            # print "setting eDens"
            self.eDens = float(cols[6])
        if len(cols) >= 4 and " ".join(cols[1:4]) == "INTERFACE PACK STAT:":
            # print "seting pack"
            self.pack = float(cols[4])
        if len(cols) >= 4 and " ".join(cols[1:4]) == "DELTA UNSTAT HBONDS:":
            # print "setting uns"
            self.uns = float(cols[4])

        return True

    def read_resultline_from_mpi_silent_logfile(self, line):
        cols = line.split()
        if cols[0] != ":":
            return False

        if (
            len(cols) >= 6
            and " ".join(cols[2:6]) == "SEPARATED INTERFACE ENERGY DIFFERENCE:"
        ):
            # print "setting dGbind"
            self.dGbind = float(cols[6])
        if len(cols) >= 5 and " ".join(cols[2:5]) == "INTERFACE DELTA SASA:":
            self.dSASA = float(cols[5])
            # print "setting dSASA", self.dSASA
        if (
            len(cols) >= 8
            and " ".join(cols[2:7])
            == "INTERFACE DELTA SASA/SEPARATED INTERFACE ENERGY:"
        ):
            # print "setting eDens"
            self.eDens = float(cols[7])
        if len(cols) >= 4 and " ".join(cols[2:5]) == "INTERFACE PACK STAT:":
            # print "seting pack"
            self.pack = float(cols[5])
        if len(cols) >= 4 and " ".join(cols[2:5]) == "DELTA UNSTAT HBONDS:":
            # print "setting uns"
            self.uns = float(cols[5])

        return True


def read_IA_file(fname):
    results = []
    lines = open(fname).readlines()
    new_IAresult1 = False
    new_IAresult2 = False
    dat = None
    for line in lines:
        # print "--", new_IAresult1, "--", line,
        if (
            line.find("protocols.analysis.InterfaceAnalyzerMover.interface_selection")
            != -1
        ):
            new_IAresult1 = True
            # print "found new_IAresult"
            continue
        elif new_IAresult1:
            new_IAresult1 = False
            new_IAresult2 = True
        elif new_IAresult2:
            new_IAresult2 = False
            dat = IA_Result()
            dat.read_resultline(line)
        elif dat:
            good = dat.read_resultline(line)
            if not good:
                # print "pushing back completed IA_Result", dat.fname
                results.append(dat)
                dat = None
    if dat:
        results.append(dat)
    return results


def read_IA_log_from_silent_file(fname):
    results = []
    lines = open(fname).readlines()
    ready1 = False
    ready2 = False
    dat = None
    for line in lines:
        # print "--", ready1, ready2, "--", line,
        if line.find("Received job id") != -1:
            dat = IA_Result()
            dat.job_index = int(line.split()[8])
            # print "read job index ", dat.job_index, "on line", line.strip()
            if dat.job_index == 0:
                dat = None
            continue
        if (
            line.find("protocols.analysis.InterfaceAnalyzerMover.interface_selection")
            != -1
        ):
            ready1 = True
        elif ready1:
            ready1 = False
            ready2 = True
        elif ready2:
            # print line
            good = dat.read_resultline_from_mpi_silent_logfile(line)
            if not good:
                # print "pushing back completed IA_Result", dat.fname
                results.append(dat)
                dat = None
                ready2 = False
    if dat:
        results.append(dat)
    return results


def find_best_interface_by_dGbind(result_list):
    best = None
    for result in result_list:
        if not best or best.dGbind > result.dGbind:
            best = result
    # print best.fname, best.dGbind
    return best


def find_best_dGbind_from_fname(fname):
    results = read_IA_file(fname)
    best = find_best_interface_by_dGbind(results)
    if not best:
        print("Trouble encountered with dGbind from file ", fname)
    return best


# take the lowest energy structure by dGbind, but return the average of the best 10
def find_avgtop10_dGbind_from_fname(fname):
    results = read_IA_file(fname)
    if len(results) == 0:
        print("Failed to read any interface-analyzer results from", fname)
        assert len(results) > 0
    sresults = sorted(results, lambda x, y: x.dGbind < y.dGbind)
    best = sresults[0]
    best.dGbind = sum([x.dGbind for x in sresults[:10]]) / 10
    return best
