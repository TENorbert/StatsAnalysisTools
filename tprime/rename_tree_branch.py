#!/usr/bin/env python
##########################################################################################
######
###### simple manipulations with ROOT trees
######
######  - rename branch in a tree
######  - integrate background fit
###### 
######
###### Author: Gena Kukartsev, June 22, 2011
######
##########################################################################################

from __future__ import division

legend = '[rename_tree_branch]:'

import sys, ROOT

#
# command line parser
#
from optparse import OptionParser
add_help_option = "./rename_tree_branch.py [options]"

parser = OptionParser(add_help_option)

parser.add_option("-i", "--input-file", dest="input_file", default=None,
                  help="Input file name")

parser.add_option("-o", "--output-file", dest="output_file", default=None,
                  help="Output file name")

parser.add_option("-n", "--name", dest="name", default=None,
                  help="Name of an object, usage depends on context")

parser.add_option("--rename-branches", dest="rename_branches", default=None,
                  help="Specify branch renaming scheme, e.g, oldmass=newmass:oldenergy=newenergy")

parser.add_option("--background-integral", dest="background_integral", default=None,
                  help="Integrate the background fit")

print legend, 'parsing command line options...',
(options, args) = parser.parse_args()
print 'done'

# rename branches, save to a new file
if options.rename_branches:

    # Get old file, old tree and set top branch address
    infile = ROOT.TFile(options.input_file, 'read')
    intree = infile.Get(options.name)
    intree.SetBranchStatus("*",1)

    # Create a new file + a clone of old tree header. Do not copy events
    outfile = ROOT.TFile(options.output_file,"recreate")
    _temp_tree = intree.CopyTree("1==1")
    outtree = intree.CloneTree(0)
    outtree.CopyEntries(intree)
   
    branch_pairs = options.rename_branches.split(':')
    for pair in branch_pairs:
        name = pair.split('=')
        print 'renaming', name[0], 'to', name[1]
        outtree.GetBranch(name[0]).SetName(name[1])

    outtree.Print()
    outfile.Write()



if options.background_integral:
    from ROOT import TF1
    f = TF1("bkgShape","exp([0]+[1]*x)*TMath::Power(x,[2])",200,2000);
    f.SetParameter(0, 20.71)
    f.SetParameter(1, -0.002268)
    f.SetParameter(2, -3.717)
    print legend, 'integral =', f.Integral(200,2000)
    
