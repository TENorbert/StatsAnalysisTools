#!/usr/bin/env python
#########################################################################
#
# hist_plot.py
#
# Plot histograms
#
# Usage: ./root2.py --hist -i infile.txt -o outfile
#
#
# Author: Gena Kukartsev, December 2011
#
#########################################################################
from __future__ import division

import ROOT
import cms_prel

def Plot(in_fname, out_fname,
         xmin = None, xmax = None, nbin = None,
         title = 'hist',
         xlabel = 'x',
         ylabel = 'number of entries',
         legend = None,
         kstest = True):

    _num = []

    with open(in_fname, 'r') as infile:
        for line in infile:
            _words = line.strip().split()
            if len(_words)>0:
                if _words[0][0]=='#':
                    continue
                _num.append(float(_words[0]))

    if xmin:
        _xmin=float(xmin)
    else:
        _xmin = min(_num)
        
    if xmax:
        _xmax=float(xmax)
    else:
        _xmax = max(_num)

    if nbin:
        _nbin=int(nbin)
    else:
        _nbin = max(5, int(len(_num)/100))

    _hist = ROOT.TH1F(title, '', _nbin, _xmin, _xmax)

    for value in _num:
        print 'test:', value
        _hist.Fill(value)
        
    _c = ROOT.TCanvas()
    _c.SetFrameLineWidth(2)
    _c.SetRightMargin(0.05)
    #_c.SetFrameFillColor(ROOT.kGray-6)
    _c.SetTickx(1)
    _c.SetTicky(1)
    
    ROOT.gStyle.SetOptStat(0)
    _hist.SetMinimum(0)
    _hist.SetLineStyle(1)
    _hist.SetLineColor(ROOT.kBlack)
    _hist.SetLineWidth(2)
    _hist.SetFillColor(ROOT.kYellow)
    _hist.SetFillStyle(1001)
    
    _hist.GetXaxis().SetNdivisions(405)
    _hist.GetYaxis().SetNdivisions(405)
    _hist.GetXaxis().SetLabelSize(0.04)
    _hist.GetYaxis().SetLabelSize(0.04)
    
    _hist.Draw()
    
    _fit = _hist.Fit('pol0','S')
    #print _fit.Chi2()
    #print _fit.Prob()
    
                                                                                
    if kstest:
        _kshist= ROOT.TH1F(title, '', _nbin, _xmin, _xmax)
        _histtemp= ROOT.TH1F(title, '', _nbin, _xmin, _xmax)

        _totcount=0
        for ibin in range(1,_nbin+1):
            _totcount = _totcount + _hist.GetBinContent(ibin)

        _binavg = int(_totcount/_nbin)
        if _binavg != _totcount/_nbin:
            print "trucating to where the number of pvalue entries is a multiple of the number of bins (for KS test vs flat)"
            from math import floor
            _totcount = floor(_totcount)
            
        for ibin in range(0,_nbin):
            for iter in range(0,_binavg):
                _kshist.Fill(ibin/_nbin)

        for value in _num[0:int(_totcount)]:
            _histtemp.Fill(value)
                
        _c2 = ROOT.TCanvas()
        _c2.SetFrameLineWidth(2)
        _c2.SetRightMargin(0.05)
        #_c2.SetFrameFillColor(ROOT.kGray-6)
        _c2.SetTickx(1)
        _c2.SetTicky(1)
        
        ROOT.gStyle.SetOptStat(0)
        _kshist.SetMinimum(0)
        _kshist.SetLineStyle(1)
        _kshist.SetLineColor(ROOT.kBlack)
        _kshist.SetLineWidth(2)
        _kshist.SetFillColor(ROOT.kYellow)
        _kshist.SetFillStyle(1001)
        
        _kshist.GetXaxis().SetNdivisions(405)
        _kshist.GetYaxis().SetNdivisions(405)
        _kshist.GetXaxis().SetLabelSize(0.04)
        _kshist.GetYaxis().SetLabelSize(0.04)

        _kshist.Draw() 
        _c2.SaveAs("ks_flathisto.png")

        print "------------------------------------------------------------------------"
        print "KS test stats                       = ", _histtemp.KolmogorovTest(_kshist)
        print "test done against flat dist, y mean = ",_binavg
        print "number of entries                   = ",_totcount
        print "-------------------------------------------------------------------------"
        
           

    # draw title
    if title:
        cms_prel.Title(title)

    # draw axis labels
    cms_prel.XLabel(xlabel)
    cms_prel.YLabel(ylabel)

    # legend
    if legend:
        xLegend = 0.50
        yLegend = 0.75
        legendWidth = 0.40
        legendHeight = 0.08
        tlegend = ROOT.TLegend(xLegend, yLegend,
                               xLegend + legendWidth,
                               yLegend + legendHeight)
        #legend.SetShadowColor(0)
        #legend.SetFillColor(0)
        #legend.SetLineColor(0)
        #legend.SetFillStyle(0)
        #legend.SetBorderSize(0)
        tlegend.SetTextSize(0.04)
        tlegend.SetTextFont(42)
        tlegend.AddEntry(_hist, legend, 'f')
        tlegend.Draw()

    cms_prel.CmsPrel()
        
    #raw_input('hit enter to continue...')

    _c.SaveAs(out_fname)
