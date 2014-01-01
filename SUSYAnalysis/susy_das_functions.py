#!/usr/bin/env python
from susy_das_inputs import *

#####################################
# Analytical integrals for bkg shapes
#####################################
def funcInt(x, func, parlist):
    import sys
    from math import exp
    if func == "f1":
        parval = parlist[0]
        if parval < 1.:
            print "PROBLEM WITH INT. EXITING."
            print "par = ", parval
            
            sys.exit()
        return -1./(x/7000.)**(parval-1.)
    elif func == "f2":
        parval = parlist[0]
        if parval > 0.:
            print "PROBLEM WITH INT. EXITING."
            print "par = ", parval
            sys.exit()
        return exp(parval*x/7000.)/parval

#####################################
# Compute fraction of PDF defined on
# range (min, max) that is in range (lo, hi)
# based on parameter values in parlist.
#####################################

def getPDFFrac(lo, hi, min, max, func, parlist):
    if func == "f1" or func == "f2":
        num = funcInt(hi , func, parlist) - funcInt(lo , func, parlist)
        den = funcInt(max, func, parlist) - funcInt(min, func, parlist)
    else:
        print "Integral not defined.  Exiting."
        sys.exit()

    if den > 0.:
        return num/den
    else:
        print "getPDFFrac num/den = %4.3f / %4.3f.  Exiting." % (num, den)
        sys.exit()


#####################################
# Combine uncertainties in quadrature
#####################################
def getMinusOneQuadSum(errlist):
    # Input is list of relative errors [1.22, 1.05, 1.35]
    # Output is total error with also with the "plus 1"
    from math import sqrt
    min1_errlist = []
    for err in errlist:
        min1_errlist.append(1.-err)
    sq_sum = 0.
    for err in min1_errlist:
        sq_sum += err*err
    return 1.+sqrt(sq_sum)

#######################################
# Function for getting error vs. ST for drawing
# basic data vs. bkg plot
#######################################
#stInfo = [lower end of range, upper end of range, ST step size (use 5 GeV)]

def getErrVST(stInfo):
    from math import sqrt
    # Assuming stAll
    # assuming nomfunc
    nomfunc = "f1"
    ststart = stInfo[0]
    ststop  = stInfo[1]
    ststep  = stInfo[2]
    errPts = [ststart+1.0]
    for ipt in range(100000):
        ival = (ststart+ststep)+ipt*ststep
        if ival > (ststop-ststep): break 
        errPts.append(ival)
    errPts.append(ststop-1.)

    errIn = {}
    errInPar = {}
    errInFnc = {}
    for istp in errPts:
        parval = par[nomfunc]
        parval_err = par[nomfunc, "err"]
        # Get error due to stat uncertainty on parameter value of nominal function:
        frac0  = getPDFFrac(istp, 10000., ststart, 10000., nomfunc, [parval])
        frachi = getPDFFrac(istp, 10000., ststart, 10000., nomfunc, [parval-parval_err])
        epar = frachi/frac0-1.
        if epar < 0.10: epar = 0.10
        errInPar[istp] = epar
        
        # Get error due to differences with f0:
        fracf5 = getPDFFrac(istp, 10000., ststart, 10000., nomfunc, [parval])
        fracf0 = getPDFFrac(istp, 10000., ststart, 10000., "f2", [par["f2"]])
        unc = (fracf5-fracf0)/fracf5
        if unc<0.05: unc = 0.05
        errInFnc[istp] = unc

        # last entry is total normalization error:
        totE = sqrt(epar*epar + unc*unc + 0.15*0.15)
        errIn[istp] = totE

    return (errPts, errIn, errInPar, errInFnc)

#########################
# Format plots
#########################
def setPadPasMargin(pad, rightMargin=0.05):
  pad.SetFrameFillStyle(1001)
  pad.SetTicks()
  pad.SetTopMargin(0)
  pad.SetFillColor(0)
  #  pad.Draw()
  #  pad.cd()
  leftMargin   = 0.16
  topMargin    = 0.1
  bottomMargin = 0.15
  pad.SetLeftMargin(leftMargin)
  pad.SetRightMargin(rightMargin)
  pad.SetTopMargin(topMargin)
  pad.SetBottomMargin(bottomMargin)
