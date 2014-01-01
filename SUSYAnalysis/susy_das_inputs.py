#!/usr/bin/env python
from math import sqrt

# Number of events generated for each
# signal mass point

# DAS_EX1.1
# define variable nsiggen (the number of signal events generated) to be 10000.
nsiggen = 10000
# DAS_EX1.2
# define variable lumi_all (integrated luminosity) to be 5.03 (/fb):
lumi_all = 5.03
# DAS_EX1.3
# Define python tuple sbrange = (lower, upper) with ST sideband range.
sbrange = ( 600, 1500)
# DAS_EX1.4
# Define dictionary: nobs = {} for storing number of observed events
# in the sideband region for 4-jet category and >=5-jet category:
nobs = {}
nobs["sb", 4]  = 1492
nobs["sb", 5]  = 455

# Store the total number of events and signal regions for each mass point
# and 4- and >=5-jet events
nobs["tot", 4] = 2085
nobs["tot", 5] = 4315

# DAS_EX2
# Add background parameter values and errors
# from fit to 2-jet+3-jet sample
# for functions f1, f2, f3.  f3 has two
# parameters:
# f1 is the nominal power law functions
# f2 is the alternative exponential
# f3 is the alternative 2nd order power law function

# par = {}
# par["f1"]   =
# par["f2"]   =
# par["f3",0] =
# par["f3",1] =
# par["f2","err"]
# par["f1","err"]
# par["f3",0, "err"]
# par["f3",1, "err"]


#####################################
# Optimized ST cut values
# (from optimizeSTcuts_susy_das.py)
#####################################

# DAS_EX3
# optcut = {}
# optcut["SS2"]  =
# ...
# optcut["SS18"]  =

###########################################
# Systematic Uncertainties
###########################################

# Take the following uncertainties on expected signal
# as given.  We will not compute these.
jes_err_all  = 1.05  # jet energy scale uncertainty
lumi_err_all = 1.022 # integrated luminosity uncertainty 
pdf_err_all  = 1.05  # parton distribution function uncertainty

# DAS_EX4.1
# Determine single number for uncertainty related to normalization of background.
# norm_err = 


# DAS_EX4:
# Add other background systematic uncertainties
systs = {}
# systs["p0", "SS2" ] = 
# ...
# systs["p0", "SS18" ] = 
# systs["func", "SS2" ] = 
# ...
# systs["func", "SS18" ] = 

# Theoretical uncertainty on the predicted cross section
# from PDF, scale, and alpha_S variation (dominated by PDF)
# is provided and comes from the NLO calculation:
systs["SS2" ,"pdf"] = 0.092
systs["SS3" ,"pdf"] = 0.092
systs["SS4" ,"pdf"] = 0.090
systs["SS5" ,"pdf"] = 0.098
systs["SS6" ,"pdf"] = 0.113
systs["SS7" ,"pdf"] = 0.126
systs["SS8" ,"pdf"] = 0.144
systs["SS9" ,"pdf"] = 0.161
systs["SS10","pdf"] = 0.184
systs["SS11","pdf"] = 0.211
systs["SS12","pdf"] = 0.243
systs["SS13","pdf"] = 0.288
systs["SS14","pdf"] = 0.341
systs["SS15","pdf"] = 0.400
systs["SS16","pdf"] = 0.464
systs["SS17","pdf"] = 0.523
systs["SS18","pdf"] = 0.575

#####################################
# Signal efficiency for optimized ST cut
# for each mass point and 4- and >=5-jet events
# (from getSigEffs_susy_das.py)
#########################################
# DAS_EX5
# Add signal efficiency numbers
# 
#effsig = {}
#effsig["SS2", 4] =
#effsig["SS2", 5] =
#...
#effsig["SS18", 4] =
#effsig["SS18", 5] =

#####################################################
# Number of events observed for each Njet category
# and mass point from getSigEffs_susy_das.py
##################################################
#DAS_EX5
# Add numbers of observed events
# 
#nobs = {}
#nobs["SS2", 4] =
#nobs["SS2", 5] =
#...
#nobs["SS18", 4] =
#nobs["SS18", 5] =


##############################
# Paste limit results here
##############################
# DAS_EX8
# Paste the observed, expected, expected-1sigma, expected-2sigma,
# expected+1sigma, expected+2sigma numbers here in the limit_info dictionary.
# limit_info = {} #OBS      #EXP     #-1sig   #-2sig  #+1sig   #+2sig
# limit_info["SS2"]  = [xx, xx, xx, xx, xx, xx]
# ...
# limit_info["SS18"]  = [xx, xx, xx, xx, xx, xx]

#####################################
# NLO signal cross sections
#####################################

# DAS_EX1.5
# Fill in the signal cross sections from
# from http://uaf-2.t2.ucsd.edu/~spadhi/slha/xsection/2012/stealth_susy/combined_tot_mg1500.txt

# sigxsecs = {}
# sigxsecs["SS2" ] = 4884.960000 # Mass = 400
# sigxsecs["SS3" ] = # Mass = 500
# ...
# sigxsecs["SS18"] = # M = 2000 

###################################
# STRING -> MASS map
###################################
sigmass = {}
sigmass["SS1"]  = 300. 
sigmass["SS2"]  = 400. 
sigmass["SS3"]  = 500. 
sigmass["SS4"]  = 600. 
sigmass["SS5"]  = 700. 
sigmass["SS6"]  = 800. 
sigmass["SS7"]  = 900. 
sigmass["SS8"]  = 1000.
sigmass["SS9"]  = 1100.
sigmass["SS10"] = 1200.
sigmass["SS11"] = 1300.
sigmass["SS12"] = 1400.
sigmass["SS13"] = 1500.
sigmass["SS14"] = 1600.
sigmass["SS15"] = 1700.
sigmass["SS16"] = 1800.
sigmass["SS17"] = 1900.
sigmass["SS18"] = 2000.

#################################################
# Add/confirm background parameter values here
#################################################

# Production processes:
# g = gluon
# s = squark
# b = antisquark
# oo = processes involving slepton or gaugino

spcList = ["gg", "sg", "ss", "sb", "oo"]

spcid = {}
spcid["gg"] = 0
spcid["sg"] = 1
spcid["ss"] = 2
spcid["sb"] = 3

# Cross section fraction
xsfrac = {}                #gg     #sg      #ss     #sb
xsfrac["SS2" , "lo"] = [0.0000,  0.0068,  0.4769,  0.4857]
xsfrac["SS3" , "lo"] = [0.0000,  0.0112,  0.5981,  0.3646]  
xsfrac["SS4" , "lo"] = [0.0000,  0.0159,  0.6814,  0.2827]
xsfrac["SS5" , "lo"] = [0.0001,  0.0215,  0.7359,  0.2234]
xsfrac["SS6" , "lo"] = [0.0000,  0.0281,  0.7781,  0.1781]       
xsfrac["SS7" , "lo"] = [0.0000,  0.0407,  0.8010,  0.1428]
xsfrac["SS8" , "lo"] = [0.0000,  0.0441,  0.8191,  0.1221]
xsfrac["SS9" , "lo"] = [0.0004,  0.0780,  0.8193,  0.0891]
xsfrac["SS10", "lo"] = [0.0012,  0.0857,  0.8249,  0.0749]
xsfrac["SS11", "lo"] = [0.0025,  0.1061,  0.8169,  0.0542]
xsfrac["SS12", "lo"] = [0.0022,  0.1564,  0.7816,  0.0458]
xsfrac["SS13", "lo"] = [0.0058,  0.2215,  0.7159,  0.0412]
xsfrac["SS14", "lo"] = [0.0136,  0.2611,  0.6728,  0.0234]
xsfrac["SS15", "lo"] = [0.0261,  0.3435,  0.5836,  0.0142]
xsfrac["SS16", "lo"] = [0.0780,  0.3634,  0.4811,  0.0098]
xsfrac["SS17", "lo"] = [0.0947,  0.4476,  0.3808,  0.0082]
xsfrac["SS18", "lo"] = [0.1741,  0.4525,  0.2521,  0.0035]
 
xsfrac["SS2" , "nlo"] = [0.0000,  0.0180,  0.4720,  0.5140]
xsfrac["SS3" , "nlo"] = [0.0000,  0.0180,  0.4720,  0.5140]
xsfrac["SS4" , "nlo"] = [0.0000,  0.0280,  0.5610,  0.4050]
xsfrac["SS5" , "nlo"] = [0.0000,  0.0420,  0.6330,  0.3250]
xsfrac["SS6" , "nlo"] = [0.0010,  0.0590,  0.6860,  0.2620]
xsfrac["SS7" , "nlo"] = [0.0010,  0.0800,  0.7160,  0.2090]
xsfrac["SS8" , "nlo"] = [0.0020,  0.1060,  0.7310,  0.1680]
xsfrac["SS9" , "nlo"] = [0.0040,  0.1400,  0.7320,  0.1350]
xsfrac["SS10", "nlo"] = [0.0070,  0.1830,  0.7150,  0.1070]
xsfrac["SS11", "nlo"] = [0.0130,  0.2370,  0.6820,  0.0860]
xsfrac["SS12", "nlo"] = [0.0240,  0.3040,  0.6270,  0.0680]
xsfrac["SS13", "nlo"] = [0.0440,  0.3820,  0.5450,  0.0510]
xsfrac["SS14", "nlo"] = [0.0790,  0.4620,  0.4450,  0.0360]
xsfrac["SS15", "nlo"] = [0.1330,  0.5250,  0.3370,  0.0240]
xsfrac["SS16", "nlo"] = [0.2110,  0.5590,  0.2300,  0.0150]
xsfrac["SS17", "nlo"] = [0.3120,  0.5490,  0.1410,  0.0080]
xsfrac["SS18", "nlo"] = [0.4290,  0.4990,  0.0780,  0.0040]

spcnevt = {}
spcnevt["SS2", 2, "gg"] = 0
spcnevt["SS2", 2, "sg"] = 0
spcnevt["SS2", 2, "ss"] = 1
spcnevt["SS2", 2, "sb"] = 1
spcnevt["SS2", 2, "oo"] = 0
spcnevt["SS2", 3, "gg"] = 0
spcnevt["SS2", 3, "sg"] = 0
spcnevt["SS2", 3, "ss"] = 11
spcnevt["SS2", 3, "sb"] = 8
spcnevt["SS2", 3, "oo"] = 0
spcnevt["SS2", 4, "gg"] = 0
spcnevt["SS2", 4, "sg"] = 0
spcnevt["SS2", 4, "ss"] = 102
spcnevt["SS2", 4, "sb"] = 74
spcnevt["SS2", 4, "oo"] = 3
spcnevt["SS2", 5, "gg"] = 0
spcnevt["SS2", 5, "sg"] = 19
spcnevt["SS2", 5, "ss"] = 678
spcnevt["SS2", 5, "sb"] = 707
spcnevt["SS2", 5, "oo"] = 14
spcnevt["SS3", 2, "gg"] = 0
spcnevt["SS3", 2, "sg"] = 0
spcnevt["SS3", 2, "ss"] = 1
spcnevt["SS3", 2, "sb"] = 0
spcnevt["SS3", 2, "oo"] = 0
spcnevt["SS3", 3, "gg"] = 0
spcnevt["SS3", 3, "sg"] = 1
spcnevt["SS3", 3, "ss"] = 27
spcnevt["SS3", 3, "sb"] = 17
spcnevt["SS3", 3, "oo"] = 1
spcnevt["SS3", 4, "gg"] = 0
spcnevt["SS3", 4, "sg"] = 2
spcnevt["SS3", 4, "ss"] = 227
spcnevt["SS3", 4, "sb"] = 109
spcnevt["SS3", 4, "oo"] = 7
spcnevt["SS3", 5, "gg"] = 0
spcnevt["SS3", 5, "sg"] = 36
spcnevt["SS3", 5, "ss"] = 1564
spcnevt["SS3", 5, "sb"] = 957
spcnevt["SS3", 5, "oo"] = 37
spcnevt["SS4", 2, "gg"] = 0
spcnevt["SS4", 2, "sg"] = 0
spcnevt["SS4", 2, "ss"] = 1
spcnevt["SS4", 2, "sb"] = 1
spcnevt["SS4", 2, "oo"] = 1
spcnevt["SS4", 3, "gg"] = 0
spcnevt["SS4", 3, "sg"] = 0
spcnevt["SS4", 3, "ss"] = 39
spcnevt["SS4", 3, "sb"] = 16
spcnevt["SS4", 3, "oo"] = 1
spcnevt["SS4", 4, "gg"] = 0
spcnevt["SS4", 4, "sg"] = 0
spcnevt["SS4", 4, "ss"] = 317
spcnevt["SS4", 4, "sb"] = 82
spcnevt["SS4", 4, "oo"] = 3
spcnevt["SS4", 5, "gg"] = 0
spcnevt["SS4", 5, "sg"] = 64
spcnevt["SS4", 5, "ss"] = 2053
spcnevt["SS4", 5, "sb"] = 834
spcnevt["SS4", 5, "oo"] = 24
spcnevt["SS5", 2, "gg"] = 0
spcnevt["SS5", 2, "sg"] = 0
spcnevt["SS5", 2, "ss"] = 2
spcnevt["SS5", 2, "sb"] = 1
spcnevt["SS5", 2, "oo"] = 0
spcnevt["SS5", 3, "gg"] = 0
spcnevt["SS5", 3, "sg"] = 0
spcnevt["SS5", 3, "ss"] = 48
spcnevt["SS5", 3, "sb"] = 11
spcnevt["SS5", 3, "oo"] = 0
spcnevt["SS5", 4, "gg"] = 0
spcnevt["SS5", 4, "sg"] = 1
spcnevt["SS5", 4, "ss"] = 302
spcnevt["SS5", 4, "sb"] = 81
spcnevt["SS5", 4, "oo"] = 7
spcnevt["SS5", 5, "gg"] = 0
spcnevt["SS5", 5, "sg"] = 76
spcnevt["SS5", 5, "ss"] = 2345
spcnevt["SS5", 5, "sb"] = 689
spcnevt["SS5", 5, "oo"] = 18
spcnevt["SS6", 2, "gg"] = 0
spcnevt["SS6", 2, "sg"] = 0
spcnevt["SS6", 2, "ss"] = 4
spcnevt["SS6", 2, "sb"] = 0
spcnevt["SS6", 2, "oo"] = 0
spcnevt["SS6", 3, "gg"] = 0
spcnevt["SS6", 3, "sg"] = 0
spcnevt["SS6", 3, "ss"] = 43
spcnevt["SS6", 3, "sb"] = 6
spcnevt["SS6", 3, "oo"] = 1
spcnevt["SS6", 4, "gg"] = 0
spcnevt["SS6", 4, "sg"] = 2
spcnevt["SS6", 4, "ss"] = 345
spcnevt["SS6", 4, "sb"] = 61
spcnevt["SS6", 4, "oo"] = 6
spcnevt["SS6", 5, "gg"] = 0
spcnevt["SS6", 5, "sg"] = 129
spcnevt["SS6", 5, "ss"] = 2498
spcnevt["SS6", 5, "sb"] = 591
spcnevt["SS6", 5, "oo"] = 19
spcnevt["SS7", 2, "gg"] = 0
spcnevt["SS7", 2, "sg"] = 0
spcnevt["SS7", 2, "ss"] = 0
spcnevt["SS7", 2, "sb"] = 0
spcnevt["SS7", 2, "oo"] = 0
spcnevt["SS7", 3, "gg"] = 0
spcnevt["SS7", 3, "sg"] = 0
spcnevt["SS7", 3, "ss"] = 44
spcnevt["SS7", 3, "sb"] = 5
spcnevt["SS7", 3, "oo"] = 0
spcnevt["SS7", 4, "gg"] = 0
spcnevt["SS7", 4, "sg"] = 3
spcnevt["SS7", 4, "ss"] = 362
spcnevt["SS7", 4, "sb"] = 50
spcnevt["SS7", 4, "oo"] = 5
spcnevt["SS7", 5, "gg"] = 2
spcnevt["SS7", 5, "sg"] = 157
spcnevt["SS7", 5, "ss"] = 2696
spcnevt["SS7", 5, "sb"] = 440
spcnevt["SS7", 5, "oo"] = 23
spcnevt["SS8", 2, "gg"] = 0
spcnevt["SS8", 2, "sg"] = 0
spcnevt["SS8", 2, "ss"] = 1
spcnevt["SS8", 2, "sb"] = 0
spcnevt["SS8", 2, "oo"] = 0
spcnevt["SS8", 3, "gg"] = 0
spcnevt["SS8", 3, "sg"] = 2
spcnevt["SS8", 3, "ss"] = 52
spcnevt["SS8", 3, "sb"] = 5
spcnevt["SS8", 3, "oo"] = 1
spcnevt["SS8", 4, "gg"] = 0
spcnevt["SS8", 4, "sg"] = 4
spcnevt["SS8", 4, "ss"] = 401
spcnevt["SS8", 4, "sb"] = 39
spcnevt["SS8", 4, "oo"] = 2
spcnevt["SS8", 5, "gg"] = 0
spcnevt["SS8", 5, "sg"] = 183
spcnevt["SS8", 5, "ss"] = 2681
spcnevt["SS8", 5, "sb"] = 360
spcnevt["SS8", 5, "oo"] = 23
spcnevt["SS9", 2, "gg"] = 0
spcnevt["SS9", 2, "sg"] = 0
spcnevt["SS9", 2, "ss"] = 1
spcnevt["SS9", 2, "sb"] = 0
spcnevt["SS9", 2, "oo"] = 0
spcnevt["SS9", 3, "gg"] = 0
spcnevt["SS9", 3, "sg"] = 1
spcnevt["SS9", 3, "ss"] = 53
spcnevt["SS9", 3, "sb"] = 3
spcnevt["SS9", 3, "oo"] = 3
spcnevt["SS9", 4, "gg"] = 0
spcnevt["SS9", 4, "sg"] = 1
spcnevt["SS9", 4, "ss"] = 460
spcnevt["SS9", 4, "sb"] = 34
spcnevt["SS9", 4, "oo"] = 1
spcnevt["SS9", 5, "gg"] = 3
spcnevt["SS9", 5, "sg"] = 282
spcnevt["SS9", 5, "ss"] = 2803
spcnevt["SS9", 5, "sb"] = 308
spcnevt["SS9", 5, "oo"] = 26
spcnevt["SS10", 2, "gg"] = 0
spcnevt["SS10", 2, "sg"] = 0
spcnevt["SS10", 2, "ss"] = 2
spcnevt["SS10", 2, "sb"] = 0
spcnevt["SS10", 2, "oo"] = 0
spcnevt["SS10", 3, "gg"] = 0
spcnevt["SS10", 3, "sg"] = 0
spcnevt["SS10", 3, "ss"] = 63
spcnevt["SS10", 3, "sb"] = 2
spcnevt["SS10", 3, "oo"] = 1
spcnevt["SS10", 4, "gg"] = 0
spcnevt["SS10", 4, "sg"] = 8
spcnevt["SS10", 4, "ss"] = 420
spcnevt["SS10", 4, "sb"] = 41
spcnevt["SS10", 4, "oo"] = 5
spcnevt["SS10", 5, "gg"] = 4
spcnevt["SS10", 5, "sg"] = 343
spcnevt["SS10", 5, "ss"] = 2636
spcnevt["SS10", 5, "sb"] = 234
spcnevt["SS10", 5, "oo"] = 23
spcnevt["SS11", 2, "gg"] = 0
spcnevt["SS11", 2, "sg"] = 0
spcnevt["SS11", 2, "ss"] = 0
spcnevt["SS11", 2, "sb"] = 0
spcnevt["SS11", 2, "oo"] = 0
spcnevt["SS11", 3, "gg"] = 0
spcnevt["SS11", 3, "sg"] = 0
spcnevt["SS11", 3, "ss"] = 52
spcnevt["SS11", 3, "sb"] = 2
spcnevt["SS11", 3, "oo"] = 1
spcnevt["SS11", 4, "gg"] = 0
spcnevt["SS11", 4, "sg"] = 11
spcnevt["SS11", 4, "ss"] = 412
spcnevt["SS11", 4, "sb"] = 28
spcnevt["SS11", 4, "oo"] = 4
spcnevt["SS11", 5, "gg"] = 6
spcnevt["SS11", 5, "sg"] = 392
spcnevt["SS11", 5, "ss"] = 2443
spcnevt["SS11", 5, "sb"] = 160
spcnevt["SS11", 5, "oo"] = 23
spcnevt["SS12", 2, "gg"] = 0
spcnevt["SS12", 2, "sg"] = 0
spcnevt["SS12", 2, "ss"] = 2
spcnevt["SS12", 2, "sb"] = 0
spcnevt["SS12", 2, "oo"] = 0
spcnevt["SS12", 3, "gg"] = 0
spcnevt["SS12", 3, "sg"] = 0
spcnevt["SS12", 3, "ss"] = 59
spcnevt["SS12", 3, "sb"] = 0
spcnevt["SS12", 3, "oo"] = 1
spcnevt["SS12", 4, "gg"] = 0
spcnevt["SS12", 4, "sg"] = 16
spcnevt["SS12", 4, "ss"] = 468
spcnevt["SS12", 4, "sb"] = 19
spcnevt["SS12", 4, "oo"] = 6
spcnevt["SS12", 5, "gg"] = 10
spcnevt["SS12", 5, "sg"] = 491
spcnevt["SS12", 5, "ss"] = 2437
spcnevt["SS12", 5, "sb"] = 124
spcnevt["SS12", 5, "oo"] = 35
spcnevt["SS13", 2, "gg"] = 0
spcnevt["SS13", 2, "sg"] = 0
spcnevt["SS13", 2, "ss"] = 1
spcnevt["SS13", 2, "sb"] = 0
spcnevt["SS13", 2, "oo"] = 0
spcnevt["SS13", 3, "gg"] = 0
spcnevt["SS13", 3, "sg"] = 2
spcnevt["SS13", 3, "ss"] = 50
spcnevt["SS13", 3, "sb"] = 1
spcnevt["SS13", 3, "oo"] = 2
spcnevt["SS13", 4, "gg"] = 1
spcnevt["SS13", 4, "sg"] = 15
spcnevt["SS13", 4, "ss"] = 444
spcnevt["SS13", 4, "sb"] = 13
spcnevt["SS13", 4, "oo"] = 6
spcnevt["SS13", 5, "gg"] = 23
spcnevt["SS13", 5, "sg"] = 640
spcnevt["SS13", 5, "ss"] = 2277
spcnevt["SS13", 5, "sb"] = 98
spcnevt["SS13", 5, "oo"] = 39
spcnevt["SS14", 2, "gg"] = 0
spcnevt["SS14", 2, "sg"] = 0
spcnevt["SS14", 2, "ss"] = 0
spcnevt["SS14", 2, "sb"] = 0
spcnevt["SS14", 2, "oo"] = 0
spcnevt["SS14", 3, "gg"] = 0
spcnevt["SS14", 3, "sg"] = 1
spcnevt["SS14", 3, "ss"] = 4
spcnevt["SS14", 3, "sb"] = 0
spcnevt["SS14", 3, "oo"] = 1
spcnevt["SS14", 4, "gg"] = 0
spcnevt["SS14", 4, "sg"] = 4
spcnevt["SS14", 4, "ss"] = 43
spcnevt["SS14", 4, "sb"] = 3
spcnevt["SS14", 4, "oo"] = 4
spcnevt["SS14", 5, "gg"] = 37
spcnevt["SS14", 5, "sg"] = 652
spcnevt["SS14", 5, "ss"] = 1697
spcnevt["SS14", 5, "sb"] = 92
spcnevt["SS14", 5, "oo"] = 24
spcnevt["SS15", 2, "gg"] = 0
spcnevt["SS15", 2, "sg"] = 0
spcnevt["SS15", 2, "ss"] = 0
spcnevt["SS15", 2, "sb"] = 0
spcnevt["SS15", 2, "oo"] = 0
spcnevt["SS15", 3, "gg"] = 0
spcnevt["SS15", 3, "sg"] = 0
spcnevt["SS15", 3, "ss"] = 1
spcnevt["SS15", 3, "sb"] = 1
spcnevt["SS15", 3, "oo"] = 3
spcnevt["SS15", 4, "gg"] = 0
spcnevt["SS15", 4, "sg"] = 2
spcnevt["SS15", 4, "ss"] = 12
spcnevt["SS15", 4, "sb"] = 4
spcnevt["SS15", 4, "oo"] = 7
spcnevt["SS15", 5, "gg"] = 60
spcnevt["SS15", 5, "sg"] = 761
spcnevt["SS15", 5, "ss"] = 1553
spcnevt["SS15", 5, "sb"] = 65
spcnevt["SS15", 5, "oo"] = 47
spcnevt["SS16", 2, "gg"] = 0
spcnevt["SS16", 2, "sg"] = 0
spcnevt["SS16", 2, "ss"] = 0
spcnevt["SS16", 2, "sb"] = 0
spcnevt["SS16", 2, "oo"] = 2
spcnevt["SS16", 3, "gg"] = 0
spcnevt["SS16", 3, "sg"] = 0
spcnevt["SS16", 3, "ss"] = 0
spcnevt["SS16", 3, "sb"] = 1
spcnevt["SS16", 3, "oo"] = 8
spcnevt["SS16", 4, "gg"] = 1
spcnevt["SS16", 4, "sg"] = 2
spcnevt["SS16", 4, "ss"] = 5
spcnevt["SS16", 4, "sb"] = 2
spcnevt["SS16", 4, "oo"] = 10
spcnevt["SS16", 5, "gg"] = 117
spcnevt["SS16", 5, "sg"] = 956
spcnevt["SS16", 5, "ss"] = 1219
spcnevt["SS16", 5, "sb"] = 55
spcnevt["SS16", 5, "oo"] = 72
spcnevt["SS17", 2, "gg"] = 0
spcnevt["SS17", 2, "sg"] = 0
spcnevt["SS17", 2, "ss"] = 0
spcnevt["SS17", 2, "sb"] = 0
spcnevt["SS17", 2, "oo"] = 2
spcnevt["SS17", 3, "gg"] = 1
spcnevt["SS17", 3, "sg"] = 0
spcnevt["SS17", 3, "ss"] = 1
spcnevt["SS17", 3, "sb"] = 1
spcnevt["SS17", 3, "oo"] = 9
spcnevt["SS17", 4, "gg"] = 0
spcnevt["SS17", 4, "sg"] = 4
spcnevt["SS17", 4, "ss"] = 3
spcnevt["SS17", 4, "sb"] = 12
spcnevt["SS17", 4, "oo"] = 9
spcnevt["SS17", 5, "gg"] = 230
spcnevt["SS17", 5, "sg"] = 1143
spcnevt["SS17", 5, "ss"] = 882
spcnevt["SS17", 5, "sb"] = 99
spcnevt["SS17", 5, "oo"] = 98
spcnevt["SS18", 2, "gg"] = 0
spcnevt["SS18", 2, "sg"] = 0
spcnevt["SS18", 2, "ss"] = 0
spcnevt["SS18", 2, "sb"] = 0
spcnevt["SS18", 2, "oo"] = 2
spcnevt["SS18", 3, "gg"] = 0
spcnevt["SS18", 3, "sg"] = 2
spcnevt["SS18", 3, "ss"] = 0
spcnevt["SS18", 3, "sb"] = 1
spcnevt["SS18", 3, "oo"] = 11
spcnevt["SS18", 4, "gg"] = 4
spcnevt["SS18", 4, "sg"] = 4
spcnevt["SS18", 4, "ss"] = 1
spcnevt["SS18", 4, "sb"] = 16
spcnevt["SS18", 4, "oo"] = 35
spcnevt["SS18", 5, "gg"] = 397
spcnevt["SS18", 5, "sg"] = 1198
spcnevt["SS18", 5, "ss"] = 590
spcnevt["SS18", 5, "sb"] = 96
spcnevt["SS18", 5, "oo"] = 112


