#!/usr/bin/env python

print "Importing modules"
import os
import optparse
from susy_das_inputs import *
from susy_das_functions import *
import ROOT
from math import exp, sqrt
import sys  # For exiting program

try:
    # ROOT stuff
    from ROOT import TFile, TTree, TH1F, TH2F, TH3F, gROOT
    # RooFit stuff
    from ROOT import RooAddPdf, RooArgList, RooArgSet, RooDataHist, RooFit, RooHistPdf, RooRealVar
except Exception, e:
    print e
    print ("Use a python that has PyROOT installed.")
    sys.exit(0)


#######################
# Get options
#######################

print "Getting options"

parser = optparse.OptionParser("usage: %prog [options]\
<input directory> \n")

parser.add_option ('--infile', dest='infile', type='string',
                   default = '',
                   help="Input data")
parser.add_option ('--o', dest='outd', type='string',
                   default = '-1',
                   help="directory for output gifs and eps")
parser.add_option ('--cut', dest='stCut', type='float',
                   default = 600.,
                   help="st cut")
parser.add_option ('--doZBI', action="store_true",
                   dest="doZBI", default=False,
                   help="Do ZBI in addition to S/sqrt(B)")

options, args = parser.parse_args()
infile = options.infile
outd = options.outd
stCut = options.stCut
doZBI = options.doZBI

if outd == "-1" :
    outdir = "./"
elif outd == "make":
    outdir = "/afs/fnal.gov/files/home/room2/jhirsch/public_html/stealthsusy/das_plots/"
else :
    outdir = outd+"/"
if not os.path.isdir(outdir) : os.system("mkdir "+outdir)

# Function for computing the Z_Bi variable from
# obs  = expected number of total events (signal+background)
# bkg  = expected number of background events
# bkge = uncertainty on the expected number of background events

def getZBI(obs, bkg, bkge):
    if not bkg > 0:
        print "Bkg must be >0 for ZBI.  Exiting."
        sys.exit()
    if obs==0: return 0
    tau = bkg/(bkge*bkge);
    n_off = tau*bkg;
    P_BI = ROOT.TMath.BetaIncomplete(1./(1.+tau), obs, n_off+1)
    if P_BI == 0.:
        P_BI = 1e-300
    Z_BI = sqrt(2.0)*ROOT.TMath.ErfcInverse(2.*P_BI)
    if Z_BI>=obs:  Z_BI=obs
    return Z_BI

# Compute the error related to changing the background parameter
# The result is a number great than 1.  That is, 15% error will be 1.15
def getBkgError(cut, func, par, parErr):
    frac0  = getPDFFrac(cut, 10000., stCut, 10000., func, [par])
    frachi = getPDFFrac(cut, 10000., stCut, 10000., func, [par-parErr])
    if frac0 != 0.:
        return frachi/frac0


# Set up root
print "Setting ROOT options"
ROOT.gROOT.SetBatch()
ROOT.gROOT.SetStyle("Plain")
ROOT.gStyle.SetOptStat(1111)
ROOT.gStyle.SetOptTitle(0)
ROOT.gStyle.SetPalette(1)
ROOT.gStyle.SetNdivisions(405,"x");
ROOT.gStyle.SetEndErrorSize(0.)
ROOT.gStyle.SetErrorX(0.001)

# Define some variables
func = "f1"
stVar = "st_all"
indir    = "susy_inputs_das/mc_"
filebase = "_das.root"
intvars = ["jets_n"]
fltvars = ["st_all"]
vars = list(intvars)
vars.extend(fltvars)
nsiggen = 10000.
step  = 50.  # Optimization step size (GeV)
nstep = 48  # Number of optimization steps
#step  = 10.  # Optimization step size (GeV)
#nstep = 240  # Number of optimization steps
lumi = 5.
eff  = {}
cut  = {}
sob  = {}
sosb = {}
zbi  = {}
gr   = {}
grscale = {}

for ss in range(2, 19):
    grscale["SS"+str(ss)] = 1.
    grscale["SS"+str(ss),"zbi"] = 1.
    
# these variables come from susy_das_inputs.py
parval    = par[func]
parvalerr = par[func, "err"]


# List of SS points to run over
ssList = [2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18]
#ssList = [5,6]
#ssList = [3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18]

######################################
# Loop over stealth SUSY mass points
######################################
gr["cut", "sob"]  = ROOT.TGraph(len(ssList))
gr["cut", "sosb"] = ROOT.TGraph(len(ssList))
gr["cut", "zbi"]  = ROOT.TGraph(len(ssList))

# Loop over the mass points
for ss in ssList:
    mod = "SS"+str(ss)
    gr[mod, "sob"]  = ROOT.TGraph(nstep)
    gr[mod, "sosb"] = ROOT.TGraph(nstep)
    gr[mod, "zbi"]  = ROOT.TGraph(nstep)


    # Open signal ntuple
    infile   = indir+"SS"+str(ss)+filebase
    tfile = ROOT.TFile(infile)
    ntp = tfile.Get("tree")
    entries = ntp.GetEntries()

    # initialize efficiency (sig) and extrapolation (bkg) dictionaries
    for ipt in range(0, nstep):
        ist = int(stCut+ipt*step)
        eff[mod, "sig", ist] = 0
        eff[mod, "bkg", ist] = 0

    # Run over the signal ntuple
    for ievt in range(entries):
        ntp.GetEvent(ievt)
        vals = {}
        # Loop over the variables to retreive from the ntuple
        for var in vars:
            # Get the value of each variable for each event
            vals[var] = ntp.GetLeaf(var).GetValue()

        if vals["jets_n"] < 4 or vals[stVar] < stCut: continue

        # DAS_EX3.1 
        # Loop over ST cut values to be optimized to determine
        # signal efficiency for each.

        # For the loop, you can use:
        # 1) the number of steps ("nstep")
        # 2) the size of each step in GeV ("step")
        # 3) the starting point ("stCut"), which is 600 GeV in this case
        # 4) the python range object
        
        # Within the cut loop, if event passes ST cut being tested,
        # update the appropriate signal efficiency
        # object from the "eff" dictionary (eff[mod, "sig", st_cut])
        # Hint1: Increment the efficiency by 1./nsiggen, where nsiggen is the number of generated events
        # Hint2: Make st_cut the integer value of the cut for convenient indexing
        


    # DAS_EX3.2
    # Loop over ST cut values to be optimized to determine
    # background "efficiency" for each.   Do this outside
    # the signal ntuple loop.
    
    # For the loop, you can use:
    # 1) the number of steps ("nstep")
    # 2) the size of each step in GeV ("step")
    # 3) the starting point ("stCut"), which is 600 GeV in this case
    # 4) the python range object
    
    # Within the cut loop, the "efficiency" should be defined
    # as the fraction of the background PDF above the ST cut as
    # (PDF integral from ST cut to 10000 GeV) / (PDF integral from ST cut to 600 GeV)
    # Review and use the getPDFFrac() function from susy_das_functions.py

    # Update the appropriate bkg "efficiency"
    # object from the "eff" dictionary (eff[mod, "bkg", st_cut])
    # Hint1: Make st_cut the integer value of the cut for convenient indexing
    

    zbi[mod]  = -999.
    sosb[mod] = -999.
    sob[mod]  = -999.
    cut[mod, "zbi" ]  = 100.
    cut[mod, "sob" ]  = 100.
    cut[mod, "sosb"]  = 100.

    print "Performing optimization for ", mod
    for ipt in range(0, nstep):
        ist = int(stCut+ipt*step)
        fst = float(stCut+ipt*step)

        
        # DAS_EX3.3
        # Compute expected number of signal events from
        # 1) integrated luminosity (lumi_all in susy_das_inputs.py)
        # 2) efficiency from (eff dictionary filled above in 3.1)
        # 3) signal cross section (sigxsecs dictionary from susy_das_inputs.py)

        # DAS_EX3.4
        # Compute expected background from
        # 1) the number of observed events in normalization region of the >=5-jet category from susy_das_inputs.py
        # 2) the fraction of the PDF in the sideband -- Review and use the getPDFFrac() function from susy_das_functions.py
        # 3) the "efficiency" fraction stored in eff[mod, "bkg",ist]

        # DAS_EX3.5
        # Compute s/sqrt(b), s/sqrt(s+b)

        # Get the uncertainty on the background expectation related to the
        # uncertainty on the background parameter.
        bkgetail = getBkgError(fst, func, parval, parvalerr)

        # Get total uncertainty
        # Input is list of relative errors [1.22, 1.05, 1.35]
        # Output is total relative error.  15% error is 1.15
        # rel_bkge = getMinusOneQuadSum([bkgetail, 1.+norm_err, systs["func", mod]])
        rel_bkge = 1.15

        # Compute background uncertainty in number of events
        # from relative error
        bkge = bkg*rel_bkge-bkg

        # DAS_EX4.1
        # Compute Z_Bi variable
        zbi = 1.
        
        # Put results into a graph, the scale "grscale" is so that all
        # plots can be visible on the same graph
        gr[mod, "sob"] .SetPoint(ipt, fst, isob )
        gr[mod, "sosb"].SetPoint(ipt, fst, isosb)
        gr[mod, "zbi"] .SetPoint(ipt, fst, izbi )
        
        if isob > sob[mod]:
            sob[mod] = isob
            cut[mod, "sob"] = ist
            cut[mod, "sigeff", "sob"] = eff[mod, "sig", ist]
            cut[mod, "bkgeff", "sob"] = eff[mod, "bkg", ist]
        if isosb > sosb[mod]:
            sosb[mod] = isosb
            cut[mod, "sosb"] = ist
            cut[mod, "sigeff", "sosb"] = eff[mod, "sig", ist]
            cut[mod, "bkgeff", "sosb"] = eff[mod, "bkg", ist]
        if izbi > zbi[mod]:
            zbi[mod] = izbi
            cut[mod,"zbi"] = ist
            cut[mod, "sigeff", "zbi"] = eff[mod, "sig", ist]
            cut[mod, "bkgeff", "zbi"] = eff[mod, "bkg", ist]

# Z_BI computation fails when nsig >> nbkg for SS2,3,4.
# The p-value from the BetaIncomplete function, which is really small (1e-400)
# is reported as zero.  Since we expect Z_BI == S/sqrt(B) in this range, it doesn't matter.
# Take Z_BI optimization to be the same as S/sqrt(B) for simplicity.

for ss in [2,3,4]:
    mod = "SS"+str(ss)
    key = (mod, "sob")
    if key in cut:
        cut[mod,"zbi"] = cut[mod, "sob"] 
        cut[mod, "sigeff", "zbi"] = cut[mod, "sigeff", "sob"]
        cut[mod, "bkgeff", "zbi"] = cut[mod, "bkgeff", "sob"]
    

ipt = 0
if doZBI:
    print "%4s & %4s & %4s & %4s \\\\" % ("Mass Point", "S/sqrt(B)", "S/sqrt(S+B)", "Z_Bi")
else:
    print "%4s & %4s & %4s \\\\" % ("Mass Point", "S/sqrt(B)", "S/sqrt(S+B)")

for ss in ssList:
    mod =  "SS"+str(ss)
    if doZBI:
        print "%4s & %4i & %4i & %4i \\\\" % (mod, cut[mod,"sob"], cut[mod,"sosb"], cut[mod,"zbi"])
    else:
        print "%4s & %4i & %4i \\\\" % (mod, cut[mod,"sob"], cut[mod,"sosb"])
        
    for fom in ["sob", "sosb", "zbi"]:
        gr["cut", fom] .SetPoint(ipt, sigmass[mod], cut[mod,fom])
    ipt += 1
    


print "SOB"
for ss in ssList:
    mod =  "SS"+str(ss) 
    print "optcut[\"%s\"] = %4i. #sob"    % (mod, cut[mod,"sob"])


if doZBI:
    print "ZBI"
    for ss in ssList:
        mod =  "SS"+str(ss) 
        print "optcut[\"%s\"] = %4i. #zbi"    % (mod, cut[mod,"zbi"])




#######################
# Make plots
#######################

cname = "optCut_v_mass"
canv = ROOT.TCanvas(cname,cname,400,424)
pad=canv.GetPad(0)
setPadPasMargin(pad)


gr["cut", "sob"].GetYaxis().SetLabelSize(0.055)
gr["cut", "sob"].GetYaxis().SetTitleSize(0.055)
gr["cut", "sob"].GetYaxis().SetTitleOffset(1.46)
gr["cut", "sob"].GetXaxis().SetLabelSize(0.055)
gr["cut", "sob"].GetXaxis().SetLabelSize(0.055)
gr["cut", "sob"].GetXaxis().SetTitleSize(0.055)
gr["cut", "sob"].GetXaxis().SetNdivisions(506)
gr["cut", "sob"].GetXaxis().SetTitleOffset(1.25)
gr["cut", "sob"].GetXaxis().SetLabelFont(62)
gr["cut", "sob"].GetYaxis().SetLabelFont(62)
gr["cut", "sob"].GetXaxis().SetTitleFont(62)
gr["cut", "sob"].GetYaxis().SetTitleFont(62)
gr["cut", "sob"].GetXaxis().SetNdivisions(504,1);
gr["cut", "sob"].GetXaxis().SetTitle("Squark Mass (GeV)")
gr["cut", "sob"].GetYaxis().SetTitle("Optimal S_{T} Threshold (GeV)")
gr["cut", "sob"].GetXaxis().SetRangeUser(stCut, stCut+nstep*step)

gr["cut", "zbi"] .SetLineColor(1)
gr["cut", "sob"] .SetLineColor(2)
gr["cut", "sosb"].SetLineColor(4)

gr["cut", "zbi"] .SetLineWidth(2)
gr["cut", "sob"] .SetLineWidth(2)
gr["cut", "sosb"].SetLineWidth(2)

gr["cut", "sob"] .Draw("al")
gr["cut", "sosb"].Draw("l")
if doZBI:
    gr["cut", "zbi"] .Draw("l")

textsize = 0.04; xstart = 0.6; ystart = 0.2; ystartleg = 0.0
legend = ROOT.TLegend(xstart, ystart, xstart+0.3, ystart+0.2)
legend.SetFillColor(0)
legend.SetTextSize(textsize)
legend.SetColumnSeparation(0.0)
legend.SetEntrySeparation(0.1)
legend.SetMargin(0.2)
legend.AddEntry(gr["cut", "sob"] , "S/#sqrt{B}"   ,"l")
legend.AddEntry(gr["cut", "sosb"], "S/#sqrt{S+B}" ,"l")
if doZBI:
    legend.AddEntry(gr["cut", "zbi"] , "Z_{Bi}"       ,"l")
legend.Draw()
for end in [".pdf",".gif"]:
    canv.SaveAs(outdir+cname+end)



cname = "stCut_opt_sob"
canv = ROOT.TCanvas(cname,cname,400,424)
pad=canv.GetPad(0)
setPadPasMargin(pad)
gr["SS2", "sob"].GetYaxis().SetLabelSize(0.055)
gr["SS2", "sob"].GetYaxis().SetTitleSize(0.055)
gr["SS2", "sob"].GetYaxis().SetTitleOffset(1.46)
gr["SS2", "sob"].GetXaxis().SetLabelSize(0.055)
gr["SS2", "sob"].GetXaxis().SetLabelSize(0.055)
gr["SS2", "sob"].GetXaxis().SetTitleSize(0.055)
gr["SS2", "sob"].GetXaxis().SetNdivisions(506)
gr["SS2", "sob"].GetXaxis().SetTitleOffset(1.15)
gr["SS2", "sob"].GetXaxis().SetLabelFont(62)
gr["SS2", "sob"].GetYaxis().SetLabelFont(62)
gr["SS2", "sob"].GetXaxis().SetTitleFont(62)
gr["SS2", "sob"].GetYaxis().SetTitleFont(62)
gr["SS2", "sob"].GetXaxis().SetNdivisions(504,1);
gr["SS2", "sob"].GetXaxis().SetTitle("S_{T} Threshold (GeV)")
gr["SS2", "sob"].GetYaxis().SetTitle("S / #sqrt{B}")
gr["SS2", "sob"].GetXaxis().SetRangeUser(stCut, stCut+nstep*step)

gr["SS2", "sob"].SetLineColor(1)
gr["SS5", "sob"].SetLineColor(2)
gr["SS10", "sob"].SetLineColor(4)
gr["SS15", "sob"].SetLineColor(6)

gr["SS2", "sob"].SetLineWidth(2)
gr["SS5", "sob"].SetLineWidth(2)
gr["SS10", "sob"].SetLineWidth(2)
gr["SS15", "sob"].SetLineWidth(2)
gr["SS2" , "sob"] .Draw("al")
gr["SS5" , "sob"] .Draw("l")
gr["SS10", "sob"] .Draw("l")
gr["SS15", "sob"] .Draw("l")
textsize = 0.04; xstart = 0.6; ystart = 0.6; ystartleg = 0.0
legend = ROOT.TLegend(xstart, ystart, xstart+0.3, ystart+0.2)
legend.SetFillColor(0)
legend.SetTextSize(textsize)
legend.SetColumnSeparation(0.0)
legend.SetEntrySeparation(0.1)
legend.SetMargin(0.2)
legend.AddEntry(gr["SS2", "sob"]  , "SS15" ,"l")
legend.AddEntry(gr["SS5", "sob"]  , "SS16" ,"l")
legend.AddEntry(gr["SS10", "sob"]  , "SS17" ,"l")
legend.AddEntry(gr["SS15", "sob"]  , "SS18" ,"l")
legend.Draw()
for end in [".pdf",".gif"]:
    canv.SaveAs(outdir+cname+end)


