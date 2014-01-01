#!/usr/bin/env python

import commands
import os
import sys
import optparse
from math import exp, sqrt
print "Importing inputs"
from susy_das_inputs import *
print "Importing functions"
from susy_das_functions import *
print "Done importing functions"
#import ROOT

print "Getting options"

parser = optparse.OptionParser("usage: %prog [options]\
<input directory> \n")

parser.add_option ('--nosyst', action="store_false",
                   dest="doSyst", default=True,
                   help="Float bkg param")
parser.add_option ('--onlyPrint', action="store_true",
                   dest="onlyPrint", default= False,
                   help="Float bkg param")
parser.add_option ('--mod', type='string',
                   dest="modstr", default='das',
                   help="modifier for file names")


options, args = parser.parse_args()
doSyst = options.doSyst
onlyPrint = options.onlyPrint
modstr = options.modstr

template = "susy_das_datacard_template.txt"
# list of signal models (mass points):
modList  = ["SS2", "SS3", "SS4", "SS5", "SS6", "SS7", "SS8","SS9","SS10", "SS11", "SS12", "SS13", "SS14", "SS15", "SS16", "SS17", "SS18"]
nom_func = "f1" # Nominal bkg function
vals = {}
lines = {}

# List of needed inputs:
# OBS         = number of observed events
# RATE_SIG    = expected number of signal events
# RATE_BKG    = expected number of background events
# BKG_FIT_ERR = uncertainty on background due to
# SIGEFF_ERR  = uncertainty on signal efficiency
valList = ["OBS", "RATE_SIG", "RATE_BKG", "BKG_FIT_ERR", "SIGEFF_ERR"]

# Get Normalization
# NBKG_NJET = (NSB4 + NSB5)/SBFRAC * FRAC_NJET * TAILFRAC

#######################################################
# Compute factors for extrapolating from sideband
# to signal regions for each mass point 
#######################################################
norm     = {}
extrap   = {}

# DAS_EX6.1
# Compute the fraction of background expected in sideband using the getPDFFrac()
# function from the susy_das_functions.py
# sbFRAC = getPDFFrac(...)

nSB_4 = nobs["sb", 4] # number of 4-jet events in sideband
nSB_5 = nobs["sb", 5] # number of 5-or-more-jet events in sideband
nSB             = nSB_4 + nSB_5  # total number of events in sideband
nSB_relerr      = sqrt(nSB)/nSB  # relative Gaussian,  statistical uncertainty

# Compute fraction of sideband events with >=5 jets
# (Useful to compute the fraction separately from the total SB numbers because
# the fraction can be taken from MC if desired.)
frac5        = float(nSB_5)/nSB # fraction of sideband events with >=5 jets
frac5_err    = sqrt((1.-frac5)*frac5/nSB) # binomial uncertainty on fraction

if frac5 > 0.:
    frac5_relerr = frac5_err/frac5 # relative error
else:
    frac5_relerr = 0.
    
# DAS_EX6.2
# Compute factor for converting the total number of
# events in the sideband (4-jet + >=5-jet) into
# the normalization of the 4-jet and >=5-jet background shapes
# norm_factor = {}
# norm_factor[4] = 
# norm_factor[5] =
# norm_factor[4, "err"] = # store XX% relative error as 1.XX
# norm_factor[5, "err"] = # store XX% relative error as 1.XX

extrapSB = {}

# For each mass point, compute the extrapolation factor for
# number of events with ST > optimized requirement:
for mod in modList:    
    for njet in [4, 5]:
        # DAS_EX6.3
        # Compute the fraction of background expected in tail with ST > optimized cut
        # using the getPDFFrac() function and the optcut dictionary from the susy_das_functions.py.
        # tailfrac = 

        
        # DAS_EX6.4
        # Compute the extrapolation factor that we will multiply by the total number of events
        # in the ST sideband to get the number of events expected with ST > optimized cut
        # extrapSB[mod, njet] =

# Compute total normalization uncertainty as quadrature sum of
# 1) relative uncertainty of extrapolation factor
# 2) relative uncertainty of number of events in sideband
norm[4, "err"] = 1.+sqrt( frac5_relerr*frac5_relerr + nSB_relerr*nSB_relerr )
norm[5, "err"] = 1.+sqrt( frac5_relerr*frac5_relerr + nSB_relerr*nSB_relerr )


##################################################
# Compute inputs to limits
##################################################

# Compute largest diff b/w func forms
for mod in modList:    
    for njet in [4, 5]:
        
        ##############################################################
        # Number of observed and expected signal and background events
        ##############################################################

        # Number of observed events from susy_das_inputs.py
        vals["OBS"        , njet, mod] = nobs[mod, njet]
        # Number of expected background events = (number of observed sideband evts) * (extrapolation factor)
        vals["RATE_BKG"   , njet, mod] = nSB*extrapSB[mod, njet]
        # Number of expected signal events = (integrated luminosity) * (signal cross section) * (signal efficiency)
        vals["RATE_SIG"   , njet, mod] = lumi_all*effsig[mod, njet]*sigxsecs[mod]

        ################################################################
        # Systematic uncertainty on expected signal 
        ################################################################
        vals["LUMI_ERR"   , njet, mod] = lumi_err_all # given in susy_das_inputs.py
        vals["JES_ERR"    , njet, mod] = jes_err_all  # given in susy_das_inputs.py
        vals["PDF_ERR"    , njet, mod] = pdf_err_all  # given in susy_das_inputs.py

        # DAS_EX6.5 : Compute binomial uncertainty on signal efficiency (eff)
        # and store it in the vals dictionary.  Remember that a 10% uncertainty should
        # be stored as 1.10
        eff    = effsig[mod, njet]
        # vals["SIGEFF_ERR" , njet, mod] = 

        ################################################################
        # Systematic uncertainty on expected background
        ################################################################
        vals["BKG_NORM_ERR", njet, mod] = norm[njet, "err"]  # given in susy_das_inputs.py
        vals["BKG_FNC_ERR", njet, mod]  = systs["func", mod] # given in susy_das_inputs.py
        vals["BKG_P0_ERR" , njet, mod]  = systs["p0"  , mod] # given in susy_das_inputs.py

        ################################################################
        # Combine systematic uncertainties
        ################################################################
        # total background error
        vals["TOT_BKG_ERR"     , njet, mod] = getMinusOneQuadSum([norm[njet, "err"], systs["func", mod], systs["p0", mod]])
        # background statistical error only
        vals["TOT_BKG_STAT_ERR", njet, mod] = getMinusOneQuadSum([norm[njet, "err"], systs["p0"  , mod]])
        # backgrond systematic error only
        vals["TOT_BKG_SYST_ERR", njet, mod] = getMinusOneQuadSum([systs["func", mod]])
        # total signal error
        vals["TOT_SIG_ERR", njet, mod]      = getMinusOneQuadSum([lumi_err_all, jes_err_all, pdf_err_all, vals["SIGEFF_ERR" , njet, mod]])


#############################
# Make datacards
#############################
if not onlyPrint:
    print "Making datacards"
    for mod in modList:
        datacard = "datacard_"+mod+"_susy_"+modstr+".txt"
            

        # Replace STRINGS in template with correct values
        commandList = []
        commandList.append("cp "+template+" "+datacard)
        for njet in [4, 5]:
            commandList.append("replace OBS_%i         %i    -- %s" % (njet, int(vals["OBS"    , njet, mod]),datacard))
            commandList.append("replace RATE_SIG_%i    %4.4f -- %s" % (njet, vals["RATE_SIG"   , njet, mod] ,datacard))
            commandList.append("replace NSB            %i    -- %s" % (int(nSB)                            ,datacard))
            commandList.append("replace EXTRAP_%i      %4.3f -- %s" % (njet, extrapSB[mod, njet],datacard))
            commandList.append("replace RATE_BKG_%i    %4.5f -- %s" % (njet, vals["RATE_BKG"   , njet, mod] ,datacard))

            # SIGNAL UNCERTAINTY
            commandList.append("replace SIGEFF_ERR_%i  %4.4f -- %s" % (njet, vals["SIGEFF_ERR" , njet, mod],datacard))
            commandList.append("replace LUMI_ERR_%i    %4.2f -- %s" % (njet, vals["LUMI_ERR"   , njet, mod],datacard))
            commandList.append("replace JES_ERR_%i     %4.2f -- %s" % (njet, vals["JES_ERR"    , njet, mod],datacard))
            commandList.append("replace PDF_ERR_%i     %4.2f -- %s" % (njet, vals["PDF_ERR"    , njet, mod],datacard))

            # BACKGROUND UNCERTAINTY
            commandList.append("replace BKG_NORM_ERR_%i %4.2f -- %s" % (njet, vals["BKG_NORM_ERR", njet, mod],datacard))
            commandList.append("replace BKG_FNC_ERR_%i %4.2f -- %s"  % (njet, vals["BKG_FNC_ERR" , njet, mod],datacard))
            commandList.append("replace BKG_P0_ERR_%i  %4.2f -- %s"  % (njet, vals["BKG_P0_ERR"  , njet, mod],datacard))

            
        for cmd in commandList:
            os.system(cmd)

#############################
# Print results table
#############################

resInfo  = {}

# Put results info into convenient dictionary
# and combine results for 4- and >=5-jet categories.

for mod in modList:
    resInfo["nobs", 4         , mod]  = int(vals["OBS", 4, mod ])
    resInfo["nobs", 5         , mod]  = int(vals["OBS", 5, mod ])
    resInfo["nsig", 4         , mod]  = vals["RATE_SIG", 4, mod]
    resInfo["nsig", 5         , mod]  = vals["RATE_SIG", 5, mod]
    resInfo["nbkg", 4         , mod]  = vals["RATE_BKG", 4, mod]
    resInfo["nbkg", 5         , mod]  = vals["RATE_BKG", 5, mod]
    
    resInfo["nsig_err", 4     , mod]  = resInfo["nsig", 4, mod] * (vals["TOT_SIG_ERR", 4, mod]-1.     )
    resInfo["nsig_err", 5     , mod]  = resInfo["nsig", 5, mod] * (vals["TOT_SIG_ERR", 5, mod]-1.     )
    resInfo["nbkg_err", 4     , mod]  = resInfo["nbkg", 4, mod] * (vals["TOT_BKG_ERR", 4, mod]-1.     )
    resInfo["nbkg_err", 5     , mod]  = resInfo["nbkg", 5, mod] * (vals["TOT_BKG_ERR", 5, mod]-1.     )
    resInfo["nbkg_syst_err", 4, mod]  = resInfo["nbkg", 4, mod] * (vals["TOT_BKG_SYST_ERR", 4, mod]-1.)
    resInfo["nbkg_syst_err", 5, mod]  = resInfo["nbkg", 5, mod] * (vals["TOT_BKG_SYST_ERR", 5, mod]-1.)
    resInfo["nbkg_stat_err", 4, mod]  = resInfo["nbkg", 4, mod] * (vals["TOT_BKG_STAT_ERR", 4, mod]-1.)
    resInfo["nbkg_stat_err", 5, mod]  = resInfo["nbkg", 5, mod] * (vals["TOT_BKG_STAT_ERR", 5, mod]-1.)

    resInfo["nobs", mod]          = resInfo["nobs", 4, mod] + resInfo["nobs", 5, mod]
    resInfo["nbkg", mod]          = resInfo["nbkg", 4, mod] + resInfo["nbkg", 5, mod]
    resInfo["nsig", mod]          = resInfo["nsig", 4, mod] + resInfo["nsig", 5, mod]
    resInfo["nsig_err", mod]      = sqrt(resInfo["nsig_err", 4, mod] * resInfo["nsig_err",4, mod] + resInfo["nsig_err",5, mod] * resInfo["nsig_err",5, mod])
    resInfo["nbkg_err", mod]      = sqrt(resInfo["nbkg_err", 4, mod] * resInfo["nbkg_err",4, mod] + resInfo["nbkg_err",5, mod] * resInfo["nbkg_err",5, mod])
    resInfo["nbkg_stat_err", mod] = sqrt(resInfo["nbkg_stat_err", 4, mod] * resInfo["nbkg_stat_err",4, mod] + resInfo["nbkg_stat_err",5, mod] * resInfo["nbkg_stat_err",5, mod])
    resInfo["nbkg_syst_err", mod] = sqrt(resInfo["nbkg_syst_err", 4, mod] * resInfo["nbkg_syst_err",4, mod] + resInfo["nbkg_syst_err",5, mod] * resInfo["nbkg_syst_err",5, mod])

    # Get mass value for model / mass point
    mass = int(sigmass[mod])

####################################################      
# Print results tables for 4- and >=5-jet categories
####################################################

for njet in [4,5]:
    print " "
    print "%i-Jet Results Table:" % njet
    print "============================"
    print "Model& Mass & ST Cut  & # Bkg Expected           & # Sig Expected    & # Obs"
    print "-----  ----   ------    -----------------------    -----------------   -----"

    for mod in modList:
        mass = int(sigmass[mod])
        if mod in ["SS14", "SS15", "SS16","SS17","SS18"]:
            print "%4s & %4i & %4i  & $%6.1f\\pm%6.1f\\pm%6.1f$ & $%7.3f\\pm%7.3f$ & %2i \\\\" % (mod,
                                                                                                  mass,
                                                                                                  optcut[mod],
                                                                                                  resInfo["nbkg", njet, mod],
                                                                                                  resInfo["nbkg_stat_err", njet, mod],
                                                                                                  resInfo["nbkg_syst_err", njet, mod],
                                                                                                  resInfo["nsig", njet, mod],
                                                                                                  resInfo["nsig_err", njet, mod],
                                                                                                  resInfo["nobs", njet, mod])
        else:
            print "%4s & %4i & %4i  & $%6.1f\\pm%6.1f\\pm%6.1f$ & $%6.1f\\pm%6.1f$ & %2i \\\\" % (mod, mass,
                                                                                                  optcut[mod],
                                                                                                  resInfo["nbkg", njet, mod],
                                                                                                  resInfo["nbkg_stat_err", njet, mod],
                                                                                                  resInfo["nbkg_syst_err", njet, mod],
                                                                                                  resInfo["nsig", njet, mod],
                                                                                                  resInfo["nsig_err", njet, mod],
                                                                                                  resInfo["nobs", njet, mod])
##############################################            
# Print systematics table
##############################################
print " "
print "Systematics Table:"
print "============================"
print "            4-jet                 ||             >=5-jet               "
print "MC Stat | Norm   Function   Shape || MC Stat | Norm   Function   Shape "
print "-----------------------------------------------------------------------"
for mod in modList:    
    sigeff_4  = 100.*(vals["SIGEFF_ERR"   , 4, mod]-1.)
    lumi_4    = 100.*(vals["LUMI_ERR"     , 4, mod]-1.)
    jes_4     = 100.*(vals["JES_ERR"      , 4, mod]-1.)
    pdf_4     = 100.*(vals["PDF_ERR"      , 4, mod]-1.)
    bkgfit_4  = 100.*(vals["BKG_NORM_ERR" , 4, mod]-1.)
    bkgfunc_4 = 100.*(vals["BKG_FNC_ERR"  , 4, mod]-1.)
    bkgp0_4   = 100.*(vals["BKG_P0_ERR"   , 4, mod]-1.)

    sigeff_5  = 100.*(vals["SIGEFF_ERR"   , 5, mod]-1.)
    lumi_5    = 100.*(vals["LUMI_ERR"     , 5, mod]-1.)
    jes_5     = 100.*(vals["JES_ERR"      , 5, mod]-1.)
    pdf_5     = 100.*(vals["PDF_ERR"      , 5, mod]-1.)
    bkgfit_5  = 100.*(vals["BKG_NORM_ERR" , 5, mod]-1.)
    bkgfunc_5 = 100.*(vals["BKG_FNC_ERR"  , 5, mod]-1.)
    bkgp0_5   = 100.*(vals["BKG_P0_ERR"   , 5, mod]-1.)

    if bkgfunc_4 < 0.5: bkgfunc_4 = 0.6

    # PRint with JES
    #    print "%s & %3.0f & %3.0f & %3.0f & %3.0f & %3.0f & %3.0f & %3.0f & %3.0f & %3.0f & %3.0f \\\\" % (mod, 
    #                                                                                                       sigeff_4, jes_4, bkgfit_4, bkgfunc_4, bkgp0_4,
    #                                                                                                       sigeff_5, jes_5, bkgfit_5, bkgfunc_5, bkgp0_5)
    # PRint without JES:
    print "%5s & %4.0f & %4.0f & %4.0f & %4.0f & %4.0f & %4.0f & %4.0f & %4.0f \\\\" % (mod, 
                                                                                       sigeff_4, bkgfit_4, bkgfunc_4, bkgp0_4,
                                                                                       sigeff_5, bkgfit_5, bkgfunc_5, bkgp0_5)

###################################################################
# Print dictionary for comparing obs and bkg (p-value computation)
###################################################################

#print " "
#print "For Input to P-value:"
#print "===================="
#print "mod  obs4  bkg4  bkg4_err obs5  bkg5  bkg5_err"
#print "----------------------------------------------"
#for mod in modList:
#    print 'bkgInfo["%3s"] = [%2i,  %6.3f,  %6.3f,  %2i,  %6.3f,  %6.3f]' % (mod, resInfo["nobs",4, mod], resInfo["nbkg",4, mod], resInfo["nbkg_err", 4, mod],
#                                                                            resInfo["nobs",5, mod], resInfo["nbkg",5, mod], resInfo["nbkg_err", 5, mod])


#################################################################
# Print results table for combined 4 + >=5 jet event samples
################################################################

#for mod in modList:
#    if mod in ["SS16","SS17","SS18"]:
#        # (use different precision for SS16-18):
#        print "%4s & %4i & %4i  & $%6.1f\\pm%6.1f\\pm%6.1f$ & $%6.2f\\pm%6.2f$ & %2i \\\\" % (mod,
#                                                                                              mass,
#                                                                                              optcut[mod],
#                                                                                              resInfo["nbkg",mod],
#                                                                                              resInfo["nbkg_stat_err",mod],
#                                                                                              resInfo["nbkg_syst_err",mod],
#                                                                                              resInfo["nsig",mod],
#                                                                                              resInfo["nsig_err",mod],
#                                                                                              resInfo["nobs",mod])
#    else:
#        print "%4s & %4i & %4i  & $%6.1f\\pm%6.1f\\pm%6.1f$ & $%6.1f\\pm%6.1f$ & %2i \\\\" % (mod, mass,
#                                                                                              optcut[mod],
#                                                                                              resInfo["nbkg",mod],
#                                                                                              resInfo["nbkg_stat_err",mod],
#                                                                                              resInfo["nbkg_syst_err",mod],
#                                                                                              resInfo["nsig",mod],
#                                                                                              resInfo["nsig_err",mod],
#                                                                                              resInfo["nobs",mod])
#


#grPBI = {}
#for mod in modList:
#    npes = 100000
#    grPBI[mod, 4] = ROOT.pbi(resInfo["nobs",4, mod], resInfo["nbkg",4, mod], resInfo["nbkg_err", 4, mod])
#    grPBI[mod, 5] = ROOT.pbi(resInfo["nobs",5, mod], resInfo["nbkg",5, mod], resInfo["nbkg_err", 5, mod])
#    grPBI[mod]    = 2.*grPBI[mod, 4]*grPBI[mod, 5]
#    print "%s & %4.3f & %4.3f & %4.3f \\\\" % (mod, grPBI[mod, 4],grPBI[mod, 5],grPBI[mod])
                                                   

