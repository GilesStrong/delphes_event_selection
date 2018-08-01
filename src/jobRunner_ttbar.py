from __future__ import division

import glob
import optparse
import os

workDir = '/lstore/cms/giles/DiTauSkims/ttbar/'
dataDir = '/gstore/t3cms/store/ser/gstrong/powheg_ttbar'
softDir = '/lstore/cms/giles/delphes_event_selection/src'

def makeJOFile(inputFile, uid, opts):
    outputFile = workDir + "ttbar_" + str(uid)
    cmd = "./delphes_event_selection "
    cmd += "-i " + inputFile
    cmd += " -o " + outputFile
    cmd += " -d " + str(opts.debug[-1])
    joName = "analysis_" + str(uid) + ".job"
    joFile = open(joName, "w")
    joFile.write("echo Beginning\ job\n")
    joFile.write("module load gcc-4.8.3\n")
    joFile.write("module load python-2.7.11\n")
    joFile.write("export PATH=/lstore/cms/giles/programs/bin:$PATH\n")
    joFile.write("export LD_LIBRARY_PATH=/lstore/cms/giles/programs/lib64/:/lstore/cms/giles/programs/lib/:/lstore/cms/giles/programs/lib/root/:/lstore/cms/giles/programs/delphes/:$LD_LIBRARY_PATH\n")
    joFile.write("source /lstore/cms/giles/programs/bin/thisroot.sh\n")
    joFile.write("cd " + softDir + "\n")
    joFile.write("echo Paths\ set\n")
    joFile.write(cmd + "\n")
    joFile.close()
    sub = "qsub " + joName
    print "Submitting: " + sub
    if not opt.dryrun:
        os.system(sub)

if __name__ == "__main__":
    parser = optparse.OptionParser(usage = __doc__)
    parser.add_option("-d", "--dryrun", dest = "dryrun", action = "store", default = False, help = "Dry run?, default False")
    parser.add_option("-q", "--queue", dest = "queue", action = "store", default = '8nh', help = "Queue for submission. Default 8nh")
    opts, args = parser.parse_args()

    os.system("rm -r " + workDir)
    os.system("mkdir " + workDir)

    for uid, sample in enumerate(glob.glob(dataDir + "*ll*to*.root")):
        submitJob(sample, str(uid), opts)
