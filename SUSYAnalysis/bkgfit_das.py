#!/usr/bin/env python

print "Importing modules"
import os
import optparse
from susy_das_functions import *
from ROOT import  gROOT,gStyle, RooDataSet, RooArgList,RooGenericPdf,RooAbsPdf
from math import exp, sqrt

import sys  # For exiting program
try:
    # ROOT stuff
    from ROOT import TFile, TTree, TH1F, TH2F, TH3F, gROOT, TCanvas, TLegend, TPave,TLatex
    # RooFit stuff
    from ROOT import RooAddPdf, RooArgList, RooArgSet, RooDataHist, RooFit, RooHistPdf, RooRealVar, RooDataSet, RooArgList
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

parser.add_option ('--nj', dest='nj', type='string',
                   default = '2',
                   help="Flat file to use (number of jets)")
parser.add_option ('--fitRange', dest='fitrange', type='float', nargs = 2,
                   default = '0',
                   help="fit range")
parser.add_option ('--o', dest='outd', type='string',
                   default = '-1',
                   help="directory for output gifs and eps")

options, args = parser.parse_args()
nj       = options.nj
outd     = options.outd
fitrange = options.fitrange

# Define dictionaries
files  = {}
tfiles = {}
hists  = {}
rhists = {}

# Options for which dataset to fit
if nj == "2":
    files["dat"] = 'susy_das_inputfiles/data_eq2jets_das.txt'
elif nj == "3":
    files["dat"] = 'susy_das_inputfiles/data_eq3jets_das.txt'
elif nj == "23":
    files["dat"] = 'susy_das_inputfiles/data_eq23jets_das.txt'
elif nj == "4":
    files["dat"] = 'susy_das_inputfiles/data_geq4jets_das.txt'

# Setup output directory for plots
if outd == "-1" :
    outdir = "./"
else :
    outdir = outd+"/"
if not os.path.isdir(outdir) : os.system("mkdir "+outdir)


# Setup ROOT
print "Setting ROOT options"
gROOT.SetBatch()
gROOT.SetStyle("Plain")
gStyle.SetOptStat(11111111)
gStyle.SetOptTitle(0)
gStyle.SetPalette(1)
gStyle.SetNdivisions(405,"x");
gStyle.SetEndErrorSize(0.)
gStyle.SetErrorX(0.001)

#########################################################
# Make fitFuncs dictionary containing all information
# for defining three fit functions
#########################################################

fitFuncs = {}

# Power law
######################
fitFuncs["f1","col"]  = 1
fitFuncs["f1","npar"] = 1
fitFuncs["f1","name"]  = "1 / x^{p0}"
fitFuncs["f1","str"]  = "1. / ( st/7000.)^(f1_p0)"
fitFuncs["f1","par"]  = ("f1_p0", 5.0, 0., 10.)

# fitFuncs[function , njet, parameter index, "par"] = (start, lower limt, upper limit)
fitFuncs["f1", 2, 0, "par"]  = (5.0, 0., 10.)
fitFuncs["f1", 3, 0, "par"]  = (5.0, 0., 10.)
fitFuncs["f1", 4, 0, "par"]  = (5.0, 0., 10.)

# Exponential
######################
fitFuncs["f2","col"]  = 4
fitFuncs["f2","npar"] = 2
fitFuncs["f2","name"]  = "e^{ p1*x }"
fitFuncs["f2","str"]  = "exp (f2_p0 + f2_p1*st/7000.)"

# fitFuncs[function , njet, parameter index, "par"] = (start, lower limt, upper limit)
fitFuncs["f2", 2, 0, "par"]  = (  0.0 ,   10., 10.)
fitFuncs["f2", 2, 1, "par"]  = ( -6.0 , -100.,  0.)

# 2nd order power law
######################
fitFuncs["f3","col"] = 2
fitFuncs["f3","npar"] = 2
fitFuncs["f3","name"]  = "1 / x^{p0 + p1*log(x)}"
fitFuncs["f3","str"]  = "1. / (st/7000.)^(f3_p0+f3_p1*log(st/7000.))"

# fitFuncs[function , njet, parameter index, "par"] = (start, lower limt, upper limit)
fitFuncs["f3", 2,  0, "par"]  = ( 30.0,   10., 50.)
fitFuncs["f3", 2,  1, "par"]  = (  3.0,   0., 40.)         

# Make ranges dictionary for storing fit range.
ranges = {}

if fitrange == 0:
    ranges["fit","lo"]  = 600.
    ranges["fit","hi"]  = 1900.
else:
    ranges["fit","lo"]  = fitrange[0]
    ranges["fit","hi"]  = fitrange[1]

nsiggen = 10000.
anamult = 5
ptcut = 20

anaRange = [2]
catList = ["dat"]

datevts = {}
sigeffs = {}

############################################
# Definitions
############################################

binWidth = 50. #GeV

# DAS_EX_2.1
# Define a RooRealVar object for the ST variable
st = RooRealVar("st", "st", ranges["fit","lo"], ranges["fit","hi"])


############################################
# Get data and models
############################################

# 
print "Reading Files"
# DAS_EX_2.2
# Read data from file into RooDataSet using the RooDataSet.read() method.
#rdata = RooDataSet.add(files["dat"], RooArgList(st))
rdata = RooDataSet.read(files["dat"], RooArgList(st))

###################################
# Set up fit PDFs with information
# from fitFuncs dictionary
###################################

# Declare a few dictionaries
pdfs = {}; pars = {}; aset = {}

for func in ["f2", "f1", "f3"]:

    # DAS_EX_2.3
    # Define and initialize (with your ST RooRealVar) a RooArgSet for use in the pdf
    # We need one RooArgSet for each function, so store in python dictionary indexed
    # by "func" string:
    aset[func] = RooArgSet(st)

    # Define fit parameters
    for ipar in range(fitFuncs[func,"npar"]):
        njet = 2  #Use same parameter definition for all jet multiplicities
        # start, lo, and hi are the starting value, lower limit, and upper limit for each parameter
        start = fitFuncs[func, njet, ipar, "par"][0]
        lo    = fitFuncs[func, njet, ipar, "par"][1]
        hi    = fitFuncs[func, njet, ipar, "par"][2]

        # DAS_EX_2.4
        # Make RooRealVar for each parameter
        # For some functions there are more than one parameter, so we store each
        # parameter in a dictionary indexed by "func" string and "ipar" integer:
        # Make sure the RooRealVar name and title are unique
        pars[func, ipar] = RooRealVar(func+"_p"+str(ipar), func+"_p"+str(ipar), start, lo, hi)

        # DAS_EX_2.5
        # Add each param to the RooArgSet aset[func]
        aset[func].add(pars[func,ipar])
        

    # DAS_EX_2.6
    # Define RooGenericPdf for each "func" from function string in the fitFuncs dictionary
    # and RooArgSet.  Make sure PDF name, title are uniqe

    # Get function string from fitFuncs dictionary
    funcStr = fitFuncs[func,"str"]

    # Make sure to use this name and title when defining the PDF
    name = "bkgpdf_"+func
    title = "bkgpdf_"+func

    # Define pdf
    # pdfs[func] = 
    pdfs[func] = RooGenericPdf(name, title, funcStr, RooArgList(aset[func]))
    
######################
# Do Fit
######################

roores = {}
for func in ["f2", "f1", "f3"]:
    # DAS_EX_2.7
    # Using the fitTo() method of the PDF, fit the data with each function
    # Pass the option RooFit.Save(True) to the fitTo method, and save each result in the roores dictionary:
    roores[func] = pdfs[func].fitTo(rdata, RooFit.Save(True))

#########################
# Make Plots
#########################

# Define canvas
cname = "bkgfit_"+str(nj)+"jet"
canv  = TCanvas(cname,cname,400,424)

# Define RooPlot
nbins = int((ranges["fit","hi"]-ranges["fit","lo"])/binWidth)
rooplot = st.frame(ranges["fit","lo"] ,  ranges["fit","hi"], nbins)
# Add data to plot
rdata.plotOn(rooplot, RooFit.Invisible())

# Add fit functions to plot
for func in ["f1", "f2", "f3"]:
    pdfs[func] .plotOn(rooplot, RooFit.LineColor(fitFuncs[func,"col"]))
    rdata.plotOn(rooplot)

# Draw RooPlot
rooplot.Draw()

############################
# Make legend/text and save plot
############################

# Get objects from rooplot for adding in legend
for func in ["f2", "f1", "f3"]:
    pdfs[func, "obj"] = rooplot.findObject("bkgpdf_"+func+"_Norm[st]")
    
# Make latex text and legend
width = 0.4
textsize = 0.04; xstart = 0.4; ystart = 0.75

latex = TLatex()
latex.SetNDC()
latex.SetTextAlign(12)
latex.SetTextSize(textsize)
latex.DrawLatex(xstart, 0.85, "CMS Preliminary")
latex.DrawLatex(xstart, 0.80, "4.96 fb^{-1}, #sqrt{s}=7TeV")

legend = TLegend(xstart, ystart-7.7*textsize, xstart+width, ystart)
legend.SetFillColor(0)
legend.SetTextSize(textsize)
legend.SetColumnSeparation(0.0)
legend.SetEntrySeparation(0.1)
legend.SetMargin(0.2)

# Add data to legend
if nj == "2":
    legend.AddEntry(rooplot.findObject("h_dataset"), "2-jet Data", "p")
elif nj == "3":
    legend.AddEntry(rooplot.findObject("h_dataset"), "3-jet Data", "p")
elif nj == "4":
    legend.AddEntry(rooplot.findObject("h_dataset"), "#geq4-jet Data", "p")
elif nj == "23":
    legend.AddEntry(rooplot.findObject("h_dataset"), "2+3-jet Data", "p")

# Add functions to legend
for func in ["f1", "f2", "f3"]:
    legend.AddEntry(pdfs[func,"obj"], fitFuncs[func,"name"], "l")

# Draw legend
legend.Draw()

# Save plot
for end in [".pdf",".gif",".png",]:
    canv.SaveAs(outdir+cname+end)
        

