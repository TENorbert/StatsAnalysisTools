#!/usr/bin/env python
####################################################################
#
# add_hist.py
#
# Create and add histograms to a root file
#
# Usage:
#        add_new_hist(filename,histname,
#                     list_of_values,list_of_errors)
#
#
# Author: Gena Kukartsev, September 2011
#
####################################################################
from __future__ import division

from ROOT import TFile
from ROOT import TH1F

def add_new_hist(filename, histname,
                 values, errors):

    _legend = '[add_hist.add_new_hist]:'

    _file = TFile(filename, 'UPDATE')
    _file.cd()

    _nbins = len(values)

    _hist = TH1F(histname, histname, _nbins, 0, _nbins)

    print _legend, values, errors

    i = 0
    for v in values:
        _hist.SetBinContent(i+1,float(v))
        _hist.SetBinError(i+1,float(errors[i]))
        i += 1

    _hist.Print()
    _hist.Write()

    _file.Close()
