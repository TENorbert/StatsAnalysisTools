#!/usr/bin/env python

################################################
#
# This script for making limit plots and such
#
#   - plot limits from ASCII file
#   
# Usage:
#        plot by number: ./plots.py --plot-number 2 (specified in tprime_plots.py)
#
# Gena Kukartsev, 2010-2012
# Mike Luk, 2011
#
#
################################################

from ROOT import TGraph
from ROOT import TMultiGraph
from ROOT import TCanvas
from ROOT import TLatex
from array import array
from ROOT import RooWorkspace
from ROOT import RooArgSet
from ROOT import RooFit
from ROOT import gROOT
from ROOT import TFile
from ROOT import TGraphAsymmErrors
from ROOT import TF1
from ROOT import RooNumIntConfig

from plotter import *

import ROOT
import math,random
from ROOT import TMath
from ROOT import TH1F
from ROOT import RooDataHist
from ROOT import RooArgList
from ROOT import RooRealVar
from ROOT import RooVoigtian

from operator import mod
import os

_legend = '[plots:]'

#
# command line parser
#
from optparse import OptionParser
add_help_option = "./metZplots.py -p PLOT_SET [other options]"

parser = OptionParser(add_help_option)

parser.add_option("-p", "--plot-set", dest="plot_set", default=None,
                  help="Which plot(s) to make", metavar="PLOTSET")

parser.add_option("-f", "--plot-format", dest="plot_format", default='png',
                  help="Plot file format", metavar="PLOTFORMAT")

parser.add_option("-l", "--lumi", dest="lumi", default=0.8,
                  help="Integrated luminosity", metavar="LUMI")
parser.add_option("--legend", dest="legend", default='channel',
                  help="Legend for the plot", metavar="LUMI")

parser.add_option("-d", "--dir", dest="dir", default=None,
                  help="Input directory", metavar="DIR")
parser.add_option("-e", "--dir2", dest="dir2", default="nofile",
                  help="2nd Input directory for comparison", metavar="DIR")

parser.add_option("-t", "--test-statistic", dest="teststat", default=None,
                  help="Value of the test statistic for look-elsewhere effect")

parser.add_option("-m", "--mass-pts",dest="masspts",default=6,
                  help="Information of limit plot - insert number of pts to consider")

parser.add_option("-c", "--plot-comp",dest="clscomp", default=False,
                  help="output 2 different limit plots ratio")

# options for asymcls
parser.add_option("-u", "--mass-min",dest="mass_min", default=400,
                  help="min of mass range")
parser.add_option("-o", "--mass-max",dest="mass_max", default=600.1,
                  help="max of mass range")
parser.add_option("-i", "--mass-inc",dest="mass_inc", default=50,
                  help="mass increment to use")
parser.add_option("-q", "--n-toys",dest="ntoys", default=20000,
                  help="how many toys for each mass point")
parser.add_option("-s", "--n-points",dest="npoints", default=200,
                  help="how many points to scan in range")
parser.add_option("-x","--poi-min",dest="poimin",default=0.,
                  help="point of interest minimum of scan")
parser.add_option("-y","--poi-max",dest="poimax",default=10.,
                  help="point of interest maximum of scan")
parser.add_option("-r","--channel",dest="channel",default="ejets",
                  help="decay channel ejets,mujets or comb")

# options for toys
parser.add_option("--plot-toys", dest="plot_toys",
                  default=False,
                  action="store_true",
                  help="Make toy-related plots")

# plots by number
parser.add_option("--plot-number", dest="plot_number", default=None,
                  action="append",
                  help="Make plot with given number as specified in tprime_plots.py")

print _legend, 'parsing command line options...',
(options, args) = parser.parse_args()
print 'done'

format = options.plot_format
lumi   = options.lumi

do_limits      = False
do_mass_limits = False
do_lee         = False
do_robust      = False
do_fit         = False
do_cls         = False
do_asymcls     = False
do_cls_comp    = False
do_plot_info   = False
do_plot_comp   = False
#tprime plots
do_tprime_limits = False

# define which plots to make based on params
if options.plot_set == 'zprime_limits':
    do_limits = True

if options.plot_set == 'mass':
    do_mass_limits = True

if options.plot_set == 'lee':
    do_lee = True

if options.plot_set == 'robust':
    do_robust = True

if options.plot_set == 'tprime_limits':
    do_tprime_limits = True

if options.plot_set == 'fit':
    do_fit = True

if options.plot_set == 'cls':
    do_cls = True

if options.plot_set == 'asymcls':
    do_asymcls = True

#if options.dir2 != 'nofile':
#    do_cls_comp= True
#    do_cls      =False

if options.clscomp:
    do_cls_comp = True
    do_plot_info= True
    do_cls     = False
    
# end of command line parser

from limit_plots import * 
from mass_limits import * 
from look_elsewhere import * 
#from zprime_robust import * 


############################################################
#
# Main
#
#

gROOT.SetStyle('Plain')

if do_limits:
    limit_plots(format)

if do_mass_limits:
    mass_limits(format)

if do_lee:
    ROOT.gStyle.SetOptStat(0)
    look_elsewhere(float(options.teststat), options.dir, format)

if do_robust:
    robust_plots(format)

if do_tprime_limits:
    from tprime_limits import *
    #tprime_limits(format)
    tprime_limits_mcmc(format)
    #tprime_cls_limits(format)

if do_cls_comp:
    # make plot comparison
    #       ./plots.py -d directory1 -e directory2 -m (optional) number of pass points to consider --plot-comp
    # outputs limit plots with both compared images
    from cls_limits import *
    cls_limits(options.dir,format)
    cls_limits(options.dir2,format)
    tprime_cls_limits_comp(options.dir,options.dir2,options.lumi, options.legend,format)
    

if do_cls:
    # - takes condor output dir as input
    # - read input root files from condor dir
    # - merge the results for each mass point and get CLs numbers
    # - make the xsec limit plot as function of mass
    #
    # Usage:
    #        ./plots.py -p cls -d <condor_dir_name> -l <lumi> --legend <#mu+jets> -f format
    #
    from cls_limits import *
    cls_limits(options.dir,format)
    tprime_cls_limits(options.dir,options.lumi, options.legend,format)

if do_asymcls:
    mode    = 'observed'    
    from datetime import datetime
    d = datetime.now()
    dt = d.strftime("%d-%b-%Y_%H-%M-%S")
    _dir = options.channel+'_'+dt
    os.system('rm *.so *.d')
    os.system('mkdir '+_dir)

    _npoints = int((options.mass_max-options.mass_min)/options.mass_inc)+1    
    _i = 0
    
    os.system('root -l -b -q -n hf_tprime.C+')
    while _i < _npoints:
        _peak = options.mass_min + _i*options.mass_inc
        os.system("root -l -q "+"'"+'run_limit2.C+("'+options.channel+'",'+'"'+mode+'",'+'"asymcls",'+str(_peak)+','+'"asymcls_'+str(_peak)+'.ascii",'+str(options.ntoys)+','+str(options.npoints)+','+str(options.poimin)+','+str(options.poimax)+',0,0,'+'""'+')'+"'")
        os.system("mv tprime_limit_cls_"+options.channel+"_"+mode+"_asymcls_"+str(_peak)+".ascii.ascii "+_dir)
        _i += 1
    os.system('./merge_ascii.py '+_dir+'/*ascii.ascii > '+_dir+'/for_limit_plot.ascii')
    from cls_limits import *
    tprime_cls_limits(_dir,options.lumi, options.legend,format)
        
if do_fit:
    fit(format)

if do_plot_info:
    # -Usage:
    #       ./plots.py -d directory1 -m (optional defaul is 6) number of pass points to consider 
    # prints to screen the numbers from limit plot
    cls_limits(options.dir,format)
    plot_info(options.dir,int(options.masspts))

if do_plot_comp:
    from cls_limits import *
    #       ./plots.py -d directory1 -e directory2 -m (optional) number of pass points to consider
    # prints to screen the difference in the limit plots of the 2 directories
    cls_limits(options.dir,format)
    cls_limits(options.dir2,format)
    plot_comp(options.dir,options.dir2,int(options.masspts))

if options.plot_toys:
    #
    # Toy-related plots
    #
    # Usage:
    #        ./plots.py --plot-toys
    #
    from tprime_plots import *
    tprime_plots(format)

    
if options.plot_number:
    #
    # Plots by number (from tprime_plots.py)
    #
    # Usage:
    #        ./plots.py --plot-number 2 --plot-number 5
    #
    from tprime_plots import *
    tprime_plots(format, options.plot_number,options.dir,options.dir2,options.channel)

    
else:
    print _legend, 'nothing to do, exiting...'
