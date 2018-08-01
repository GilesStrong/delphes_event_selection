from __future__ import division

import glob
import optparse
import os

workDir = '/afs/cern.ch/work/g/gstrong/private/DiTauSkims/'
dataDir =  '/eos/user/g/gstrong/cms/n-tuples/DiTau_MassRegression/'

if __name__ == "__main__":
    os.system("rm -r " + dataDir + 'signal')
    os.system("mkdir " + dataDir + 'signal')

    masses = glob.glob(workDir + '*')
    print len(masses), "mass points found"
    
    for masspoint in masses:
        mass = masspoint[masspoint.rfind('/')+1:]
        samples = glob.glob(masspoint + '/*')
        print len(samples), "samples found for mass point", mass
        os.system('hadd ' dataDir + 'signal/' + mass + '.root ' + " ".join(map(str, samples)))
