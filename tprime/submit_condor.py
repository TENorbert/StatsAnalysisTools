#!/usr/bin/env python

################################################
#
# Submit graviton combination to Condor
#
#
# Gena Kukartsev, April 2011
#
# July 2011: adapted for t' semileptonic analysis
#
################################################

import sys
import os

#-------> parameters ----------------------------------------------------
#
# channel: ejets, mujets, combined
#
#
#    mode: observed,
#          expected,
#          mass_
#

simulate = True

channel = 'ejets'

mode = 'observed'

method = 'cls'

prefix = channel         # this may be changed

mass_min  =     400
mass_max  =     500.1
mass_inc  =      50.0
n_iter    =   20000
n_burn_in =     100

exp_ntoys_per_job = 125           # these params are now used for observed as well
exp_ntoys_per_mass_point = 5000

# CLs stuff
npoints = 20
poimin = 0.0
poimax = 4.0  # does not matter now: estimating with PLR anyway
npoints_per_job = 20

#
#------------------------------------------------------------------------
                            
legend = '[submit_condor]:'

# create workdir
from datetime import datetime
d = datetime.now()
dt = d.strftime("%d-%b-%Y_%H-%M-%S")
_dir = prefix+'_'+dt
os.system('rm *.so *.d')
os.system('mkdir '+_dir)

# create the tar file
os.system('root -l -b -q -n hf_tprime.C+')
#os.system('tar -czhvf '+_dir+'/tarfile.tgz cmslpc_standalone_setup.sh twobody_C.so run_diphoton_comb.C ws_*.root data361cmssw_diphoton.root')
os.system('tar -czhvf '+_dir+'/tarfile.tgz cmslpc_standalone_setup.sh hf_tprime_C.so run_limit.C results_*/*root')
os.chdir(_dir)

prefix = prefix + '_' + mode

cmd_file_name = prefix + '.cmd'

cmd_file = open(cmd_file_name, "w")

stdout_save = sys.stdout
sys.stdout = cmd_file

print '''
# -*- sh -*- # for font lock mode
# variable definitions
- env = source ./tardir/cmslpc_standalone_setup.sh
- tag = 
- output = outputFile=
- tagmode = none
- tarfile = tarfile.tgz
- untardir = tardir

'''



_npoints = int((mass_max-mass_min)/mass_inc)+1
_njobs = _npoints

_i = 0
while _i < _njobs:
    _peak = mass_min + _i*mass_inc
    
    # CLS specific: split scan by two points per job
    n_scan_jobs = 1
    if method == 'cls':
        n_scan_jobs = int( (npoints+npoints_per_job-1)/npoints_per_job )
        cls_step = (poimax-poimin)/float(npoints - 1)
    else:
        cls_step = 0.0

    _j = 0
    while _j < n_scan_jobs:

        #print legend,
        print prefix+'_'+str(_peak)+'_'+str(_j)+'_$(JID).log',

        # commented out because ntoys is now used for observed too,
        # to run several times and everage
        #if mode == 'observed' or mode == 'mass':
        #    exp_ntoys_per_job = 1

        #    print 'root -l -b -q -n tardir/run_diphoton_comb.C(' + \
        print 'root -l -b -q -n tardir/run_limit.C(' + \
              '"'+channel+'",' + \
              '"'+mode+'",' + \
              '"'+method+'",' + \
              str(_peak)+',' + \
              '"'+prefix+'_'+str(_peak)+'_'+str(_j)+'_$(JID)",' + \
              str(exp_ntoys_per_job)+',' + \
              str(npoints)+',' + \
              str(poimin+float(npoints_per_job*_j)*cls_step)+',' + \
              str(poimin+(float(npoints_per_job*_j)+float(npoints_per_job)-1.0)*cls_step)+',' + \
              str(n_iter)+',' + \
              str(n_burn_in)+',' + \
              '"tardir/"' + \
              ')'
        
        print ''
        
        _j += 1
    
    _i += 1
    
cmd_file . close()
sys.stdout = stdout_save

if mode == 'expected':
    _nsubmit = int(exp_ntoys_per_mass_point/exp_ntoys_per_job)
    if exp_ntoys_per_mass_point % exp_ntoys_per_job != 0:
        _nsubmit +=1
else:
    #_nsubmit = 1
    _nsubmit = int(exp_ntoys_per_mass_point/exp_ntoys_per_job)
    if exp_ntoys_per_mass_point % exp_ntoys_per_job != 0:
        _nsubmit +=1


for i in range(0, _nsubmit):
    if simulate:
        os.system('echo ../runManySections.py --submitCondor '+cmd_file_name)
    else:
        os.system('../runManySections.py --submitCondor '+cmd_file_name)

print
#print legend, 'data:                ', data
print legend, 'mode:                ', mode
print legend, 'method:              ', method
print legend, 'mass points:         ', _npoints
print legend, 'channel:             ', channel
if mode == 'observed':
    print legend, 'number of jobs:      ', _njobs*n_scan_jobs*_nsubmit
if mode == 'expected':
    print legend, 'toys per job:                  ', exp_ntoys_per_job
    print legend, 'toys per mass point requested: ', exp_ntoys_per_mass_point
    print legend, 'toys per mass point submitted: ', _nsubmit*exp_ntoys_per_job
    print legend, 'jobs per mass point submitted: ', _nsubmit
    print legend, 'TOTAL jobs submitted:          ', _npoints*_nsubmit
print
_whoami = os.environ['USER']
print legend, 'on LPC, monitor job progress with condor_q -global -submitter', _whoami
print    

