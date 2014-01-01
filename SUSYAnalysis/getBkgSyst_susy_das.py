#!/usr/bin/env python

import commands
import os
import sys
import optparse
from math import exp, sqrt
from susy_das_inputs import *
from susy_das_functions import *

parser = optparse.OptionParser("usage: %prog [options]\
<input directory> \n")


options, args = parser.parse_args()

# Nominal function is "f1", which 1/ST^p:
nomfunc = "f1"
start = sbrange[0]

# Determine systematic uncertainty for 90 steps of 10 GeV each:
step = 10.
nsteps = 90


modList  = ["SS2", "SS3","SS4","SS5", "SS6", "SS7", "SS8","SS9","SS10", "SS11", "SS12", "SS13", "SS14", "SS15", "SS16", "SS17", "SS18"]

######################################################################
# Estimate systematic uncertainty related to statistical uncertainty
# on background parameter from fit to 2+3-jet data
######################################################################
for mod in modList:
    # Get cut value
    cut    = optcut[mod]
    # Get nominal parameter value from par dictionary of susy_das_inputs.py
    parval = par[nomfunc]
    # Integrate background PDF from cut -> infinity (~10000), for nominal parameter value
    frac0  = getPDFFrac(cut, 10000., start, 10000., nomfunc, [parval])
    # To run interactively:
    # python
    # >>> from susy_das_inputs import *
    # >>> from susy_das_functions import *
    # >>>
    # >>> frac0 = getPDFFrac(800., 10000., 600., 10000., "f1", [5.019])
    # >>> print frac0
    # 
    

    # DAS_EX4.1
    # Similar to the above integration, integrate the background PDF from cut -> infinity,
    # for nominal parameter value - uncertainty
    # Get the uncertainty from the "par" dictionary of susy_das_inputs.py:
    # frac_varied =

    # DAS_EX4.2
    # Compute the systematic uncertainty from nominal frac0 and varied frac_varied.
    # A 10% relative uncertainty should be written 1.10.  Numbers must be positive and >1
    # syst_temp =

    # Print the results for loading into susy_das_inputs.py
    print "systs[\"p0\", \"%s\"] = %4.2f" % (mod, syst_temp)
    
######################################################################
# Estimate systematic uncertainty related to choice of nominal fit
# function (1/x^p) by comparing to exponential e^(p*x)
######################################################################

for mod in modList:
    # Get cut value
    cut    = optcut[mod]
    # Integrate background PDF from cut -> infinity, for nominal function
    frac0  = getPDFFrac(cut, 10000., start, 10000., nomfunc, [parval])

    # DAS_EX4.3
    # This time vary the function from nomfunc -> "f2" and and the paramter value
    # appropriate for "f2" from susy_das_inputs.py.  The integrate:
    # frac_varied = 

    # DAS_EX4.4
    # Compute the systematic uncertainty from nominal frac0 and varied frac_varied:
    # A 10% relative uncertainty should be written 1.10.  Numbers must be positive and >1
    # syst_temp =

    # Print the results for loading into susy_das_inputs.py
    print "systs[\"func\", \"%s\"] = %4.2f" % (mod, syst_temp)
