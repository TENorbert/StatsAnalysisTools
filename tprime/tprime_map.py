#!/usr/bin/env python
#########################################################################
#
# tprime_map.py
#
# Checks and manipulations with bin maps
#
# Usage: 
#        ./tprime.py --bin-map dummy -f input.binmap --channel '#mu+jets'
#
#
# Author: Gena Kukartsev, January 2012
#
#########################################################################
from __future__ import division

import ROOT
import cms_prel

def CountBins(options):
    #
    # count bins in the map
    #

    legend = '[tprime_map.CountBins]:'

    _nbins = 36
    _nmaxbins = 36
    _textsize = 0.065

    cms_prel.SetPrlStyle()
    ROOT.gStyle.SetOptLogy(0)
    ROOT.gStyle.SetPadBottomMargin(0.15)
    ROOT.gStyle.SetPadRightMargin(0.15)
    ROOT.gStyle.SetPadLeftMargin(0.15)

    hist = ROOT.TH2F('sb_map', '', 45,100,1000,90,200,2000)
    
    if options.in_file == None:
        print legend, 'no input map file specified, exiting'
        return 0
    else:
        infile = options.in_file[0]

    rows = []
    bins = []
    with open(infile) as mapfile:
        for line in mapfile:
            if len(line.strip())==0:
                continue

            row = line.strip().split(':')

            _row = int(row[0].strip())
            rows.append( _row )
            for bin in row[1].split():
                _bin = int(bin.strip())
                bins.append( _bin )

                # fill in the visualization hist
                #hist.SetBinContent(_bin, _row)
                # reverse
                hist.SetBinContent(_bin, _nbins-_row+1)

    print legend, len(rows), 'bins in the rebinned template'
    print legend, len(bins), 'bins in the original template'

    c = ROOT.TCanvas('c','c')
    c.SetFrameLineWidth(2)

    #ROOT.gStyle.SetOptStat(0)
    #hist.GetXaxis().SetNdivisions(405)
    #hist.GetYaxis().SetNdivisions(405)
    #hist.GetXaxis().SetLabelSize(0.04)
    #hist.GetYaxis().SetLabelSize(0.04)
    #hist.SetContour( len(rows) )
    hist.SetContour( _nmaxbins )
    #ROOT.gStyle.SetPalette( 1 )
    hist.SetMaximum( _nmaxbins )
    hist.SetMinimum(1)
    hist.Draw('COLZ')

    cms_prel.XLabel('M_{fit} [GeV]', 0.90, 0.03, text_size = _textsize)
    cms_prel.YLabel('H_{T} [GeV]', 0.04, 0.90, text_size = _textsize)
    title = 'Bin map: '
    if options.channel:
        title +=options.channel
    cms_prel.Title('t\'#bar{t\'} (550 GeV)', 0.90, 0.92, text_size = _textsize)
    cms_prel.CmsSim(0.15, 0.92, text_size = _textsize)
    cms_prel.Cms7Tev(0.45, 0.92, text_size = _textsize)
    cms_prel.YLabel('s/b Rank', x=0.99, y=0.95, text_size = _textsize)
    cms_prel.Legend(title, 0.20, 0.60, text_size = _textsize)

    #c.SaveAs('binmap.pdf')
    c.SaveAs('binmap.C')
    raw_input('press enter to continue')
