from h3h4_interface_job import H3H4TwoMutInterfaceMSDJob
from generic_msd.opt_holder import OptHolder
from generic_msd.msd_job_management import MSDJobManager, JobExecutionOptions
from generic_msd.server_identification import KnownComputers, ServerIdentifier


import blargs
import os

def base_dir():
    """This function can ony be invoked by functions in setup_scheme3"""
    currpath = os.path.dirname(os.path.abspath(__file__))
    return currpath.rpartition("/")[0] + "/"


#def recursively_rm_directory(dirname):
#    fnames = os.listdir(dirname)
#    for fname in fnames:
#        fname_path = dirname + "/" + fname
#        if os.path.isdir(fname_path):
#            recursively_rm_directory(fname_path)
#        else:
#            os.remove(fname_path)
#    os.rmdir(dirname)

if __name__ == "__main__":
    opts = OptHolder()
    print("basedir", base_dir())
    opts["base_dir"] = base_dir()
    si = ServerIdentifier()

    with blargs.Parser(opts) as p:
        H3H4TwoMutInterfaceMSDJob.add_options(p)
        JobExecutionOptions.add_options(p, si)
        print( "parser", p._options)

    msd_job = H3H4TwoMutInterfaceMSDJob(opts)

    msd_manager = MSDJobManager(msd_job, opts, si)
    msd_manager.prepare_job()

