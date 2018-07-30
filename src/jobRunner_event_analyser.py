from __future__ import division

import glob
import optparse
import os

workDir = '/afs/cern.ch/work/g/gstrong/private/DiTauSkims/'
dataDir =  '/eos/cms/store/cmst3/group/exovv/clange/Xtautau/'
fileName = 'Xtautau_delphes_events.root'
softDir = '$HOME/delphes_event_selection/src'

def submitJob(inputFiles, mass, queue, dryrun, uid, truth):
    job = 'export HOME=/afs/cern.ch/user/g/gstrong/\n'
    job += 'source /afs/cern.ch/sw/lcg/external/gcc/4.9.3/x86_64-slc6/setup.sh\n'
    job += 'source /afs/cern.ch/sw/lcg/app/releases/ROOT/6.06.00/x86_64-slc6-gcc49-opt/root/bin/thisroot.sh\n'
    job += 'export LD_LIBRARY_PATH=$HOME/programs/Delphes-3.3.2/:$HOME/programs/lib/:$HOME/programs/lib64/:$LD_LIBRARY_PATH\n'
    job += 'export PATH=$HOME/miniconda2/bin:$HOME/programs/bin:$PATH\n'

    job += 'mkdir ' + workDir + mass + '\n'

    job += softDir + './delphes_event_selection'
    job += ' -i ' + inputFile
    job += ' -o ' + workDir + mass + '/' + uid
    job += ' -t ' + truth

    jobName = mass + '_' + uid + '.job'
    with open(jobName, 'w') as f:
        f.write(job)

    sub = "bsub -q " + queue + " " + jobName
    print "Submitting: " + sub

    os.system("chmod 744 " + jobName)
    if not dryrun:
        os.system(sub)

if __name__ == "__main__":
    parser = optparse.OptionParser(usage = __doc__)
    parser.add_option("-d", "--dryrun", dest = "dryrun", action = "store", default = False, help = "Dry run?, default False")
    parser.add_option("-q", "--queue", dest = "queue", action = "store", default = '8nh', help = "Queue for submission. Default 8nh")
    parser.add_option("-t", "--truth", dest = "queue", action = "store", default = False, help = "Produce gen info? Deafult False")
    opts, args = parser.parse_args()

    os.system("rm -r " + workDir)
    os.system("mkdir " + workDir)

    masses = glob.glob(dataDir + '*')
    print len(masses), "mass points found"
    
    for masspoint in masses:
        mass = masspoint[masspoint.rfind('/')+1:]
        samples = glob.glob(masspoint + '/*')
        print len(samples), "samples found for mass point", mass
        for uid, sample in enumerate(samples):
            submitJob(sample, mass, opts.queue, opts.dryrun, str(uid), opts.truth)
