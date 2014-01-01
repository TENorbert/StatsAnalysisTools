#!/usr/bin/env python

import commands
import os
import sys
import optparse
from math import exp, sqrt
from susy_das_inputs import *
from susy_das_functions import *

print "Getting options"

parser = optparse.OptionParser("usage: %prog [options]\
<input directory> \n")

parser.add_option ('--compareLO_NLO', action="store_true",
                   dest="compareLO_NLO", default=False,
                   help="Compare LO and NLO using 2011 cuts")

options, args = parser.parse_args()
compareLO_NLO = options.compareLO_NLO

# Need to use 2011 cuts to compare LO and NLO,
# since we have sub-process efficiencies
# only for 2011 cuts
if compareLO_NLO:
    optcut = optcut_2011

# edit
ffdir    = "susy_das_inputfiles/"
ffname   = "_geq4jets_das.txt"
modList  = ["SS2", "SS3","SS4","SS5","SS6","SS7", "SS8","SS9","SS10", "SS11", "SS12", "SS13", "SS14", "SS15", "SS16", "SS17", "SS18"]

# Read in data file and put lines of file in dictionary for later looping
lines = {}
ffname_data   = "_geq4jets_das_2col.txt"
lines["data"] = open(ffdir+"data"+ffname_data, "r").readlines()

##################################################################
# Count events for all mass points, three ST ranges,
# 1) ST > optimized ST threshold
# 2) ST in sideband
# 3) ST > sideband
# and
# four njet categories 2, 3, 4, >=5
################################################################
typeList = ["opt", "sb", "gtsb"] 
njetList = [2,3,4,5]

# initial event counters for sideband and all events for 4 and >=5 jet bins:
nevt = {}
nevt["sb", 4] = 0
nevt["sb", 5] = 0
nevt["tot", 4] = 0
nevt["tot", 5] = 0

# Loop over mass points
for mod in modList:
    # initial event counts for signal MC and data
    for njet in njetList:
        for type in typeList:
            nevt[mod, njet, type] = 0
            nevt["data", mod, njet, type] = 0

    # Read in signal MC file as a list of lines
    lines[mod] = open(ffdir+mod+ffname, "r").readlines()
    # Loop over the lines in signal MC
    for line in lines[mod]:
        # get variables from ntuple: number of jets and ST
        njet = int(line.split()[0])
        st   = float(line.split()[1])
        # Count events with ST > optimized cut
        if st>optcut[mod]:  nevt[mod, njet, "opt"] += 1
            
        # Count events with ST in sideband and > sideband:
        if st>sbrange[0] and st<sbrange[1]: nevt[mod, njet, "sb"] += 1
        elif st > sbrange[1]: nevt[mod, njet, "gtsb"] += 1
            
    # Loop over data:
    for line in lines["data"]:
        njet = int(line.split()[0])
        st   = float(line.split()[1])

        # Count data events passing optimized cut for each mass point
        if st>optcut[mod]:
            nevt["data", mod, njet, "opt"] += 1

        # Count data events in sideband and with ST>600 GeV (all events).
        # Since these counts do not depend on the
        # optimized ST cut, count just once during first mass point (SS2):
        if mod == "SS2":
            nevt["tot",njet] += 1
            if st>sbrange[0] and st<sbrange[1]:
                nevt["sb", njet] += 1



####################################################################
# After counting compute efficiency and number of expected events
####################################################################
eff = {}
nexp = {}
for mod in modList:  # for each signal model (mass point)
    for njet in njetList: # for 2,3,4, >=5 jets
        for type in typeList: # for sideband region, optimal region, and ST>700 GeV
            # DAS_EX5.1
            # Compute signal efficiency from event count for each mod, njet, and type using nevt dictionary
            # and number of events generated (nsiggen in susy_das_inputs.py):
            # eff[mod, njet, type] =

            # DAS_EX5.2
            # Compute statistical uncertainty on efficiency from efficiency and
            # number of events generated (nsiggen in susy_das_inputs.py) assuming binomial statistics:
            # eff["err", mod, njet, type] =

            # DAS_EX5.3
            # Compute the number of signal events expected  with uncertaintiy
            # from efficiency, signal cross section, and luminosity:
            # nexp[       mod, njet, type] = 
            # nexp["err", mod, njet, type] = 

    
###############################################
# Print results
###############################################

# Efficiency dictionary to add to susy_das_utils
print ""
print ""
print "SIGNAL EFFICIENCY (LO)"
print "======================"

for mod in modList:
    for njet in [4,5]:
        print "effsig[\"%s\", %i] = %5.4f" % (mod, njet, eff[mod, njet, "opt"])


# Number of observed events to add to susy_das_utils
print ""
print ""
print "# OBSERVED SIDEBAND EVENTS"
print "=========================="
for njet in [4,5]:
    print "nobs[\"sb\", %i] = %4.0f" % (njet, nevt["sb", njet])

print ""
print "# OBSERVED TOTAL EVENTS"
print "=========================="
for njet in [4,5]:
    print "nobs[\"tot\", %i] = %4.0f" % (njet, nevt["tot", njet])

print ""
print "# OBSERVED EVENTS WITH ST > OPTIMIZED CUT"
print "========================================="
for mod in modList:
    for njet in [4,5]:
        print 'nobs["%s", %i] = %i' % (mod, njet, nevt["data", mod, njet, "opt"])


if compareLO_NLO:

###############################################
# Convert LO efficiency to "NLO effiency"
###############################################

# At NLO the production sub-processes (SPC) contribute
# in different proportion than at LO.
# Assuming that the SPC efficiency is the same at LO
# and NLO, we can modify the LO efficiency from signal MC
# to account for the LO/NLO difference in SPC contributions.

# The sub-process list comes from susy_das_inputs:
# spcList = ["gg", "sg", "ss", "sb", "oo"]

# Production sub-processes are combinations of
# g = gluon
# s = squark
# b = antisquark
# oo = processes involving slepton or gaugino

    # number of SPC
    npc = len(spcList)
    # loop over signal models (masses) and # of jets
    for mod in modList:
        for njet in [4,5]:
            enlo = 0. # initial NLO efficiency
            for spc in spcList: # loop over all SPC
                # Skip SPC if it is NOT one of 4 dominant SPC; i.e. if it involves slepton or gaugino
                if spc != "oo":
                    # Check that LO cross section fraction for SPC is > 0; if not, ignore.
                    if xsfrac[mod, "lo"][spcid[spc]] > 0.:  # spcid just turns "spc" string into integer
                        lo_eff        = spcnevt[mod, njet, spc]/nsiggen
                        nlo_xsec_frac = xsfrac[mod, "nlo"][spcid[spc]]
                        lo_xsec_frac  = xsfrac[mod, "lo"][spcid[spc]]
                        # DAS_EX5_BONUS
                        # Compute reweighted NLO efficiency with sum over SPC using
                        # 1) LO PYTHIA efficiency for SPC (lo_eff)
                        # 2) Fraction of total cross section for SPC at NLO (nlo_xsec_frac)
                        # 3) Fraction of total cross section for SPC at LO (lo_xsec_frac)
                        # NLO efficiency = SUM_over_SPC of (PYTHIA efficiency for SPC) * (XSEC fraction at NLO for SPC)/(XSEC fraction at LO for SPC)
            print 'effsig["%s", %i] = %5.4f #2011' % (mod, njet, enlo)

            # Fill efficiency:
            eff[mod, njet, "nlo"] = enlo
            # Compute error assuming binomial statistics:
            eff["err",mod, njet, "nlo"] = sqrt(enlo*(1.-enlo)/nsiggen)


############################
# Print Tables
############################

if compareLO_NLO:
    print ""
    print ""
    print "COMPARE NLO AND LO EFFICIENCIES"
    print "==============================="

    print "Model : eff(LO, 4-jet) : eff(LO, 5-jet) : eff(NLO, 4-jet) : eff(NLO, 5-jet)"

    for mod in modList:
        print "%4s & $%3.1f\\pm%3.1f$ &  $%3.1f\\pm%3.1f$ & $%3.1f\\pm%3.1f$ &  $%3.1f\\pm%3.1f$\\\\" % (mod,
                                                                                                         100.*eff[mod,4, "opt"], 100.*eff["err",mod,4,"opt"],
                                                                                                         100.*eff[mod,5, "opt"], 100.*eff["err",mod,5,"opt"],
                                                                                                         100.*eff[mod,4, "nlo"], 100.*eff["err",mod,4,"nlo"],
                                                                                                         100.*eff[mod,5, "nlo"], 100.*eff["err",mod,5,"nlo"],
                                                                                                         ) 

#############################################
# Print numbers of events with ST > 700
#############################################
#print "Printing # Events with ST > 700 GeV"
#for mod in modList:
#    print "%4s" % mod,
#    for type in ["sb", "700"]:
#        for njet in njetList:
#            if nexp[mod, njet, type] != 0.:
#                print "& $%5.1f\\pm%5.1f$ " % (nexp[mod, njet, type], nexp["err", mod, njet, type]),
#            elif njet == 2 or njet == 3:
#                print "& -- ",
#            else:
#                print "& --              ",
#    print "\\\\\n",
