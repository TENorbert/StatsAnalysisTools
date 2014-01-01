#!/usr/bin/env python

################################################
#
# This script calculates significance with LEE
# and plots the p-value plot
#
#   - input data from ASCII files
#   
# Usage:
#       ./plots.py -p lee -t '2.692381' -f zprime_significance_714ee_746mumu_22-Jun-2011_20-49-11/*ascii
#         
# Gena Kukartsev, 2011
#
################################################

import ROOT
import fileinput
import math
from plotter import *
from ROOT import TMath
from ROOT import TCanvas
from ROOT import TH1F
from ROOT import TLine


def compare_columns(a, b):
    # sort on ascending index 1, then 0
    return cmp(a[1], b[1]) or cmp(a[0], b[0])



def look_elsewhere(observed_test_stat, dir, plot_format = 'png'):
    sys.argv.pop(1)
    sys.argv.pop(1)
    sys.argv.pop(1)
    sys.argv.pop(1)
    sys.argv.pop(1)

    intLumiEe = 1.1;
    intLumiMumu = 1.1;

    if observed_test_stat == None:
        print 'Observed test statistic value not specified (use -t value)'

    _limits =[]

    for line in fileinput.input():
        #print legend, line.strip()
        
        _words = line.strip().split()
        
        if _words[0][0] == '#':
            continue

        _number = []
        for word in _words:
            _number.append(float(word))
            
        _limits.append(_number)

    _limits.sort(compare_columns)

    _h = TH1F("-2 ln L/L_{0}", "", 50, 0, 10)

    # print sorted result
    #_p = 3.734019
    #_p = 3.308118 # 3 events
    #_p = 3.798105  # 4 events
    #_p = 3.55  # 4 events + flat
    #_p = 3.219322  # 5 events fake 600/pb strange mumu bg
    #_p = 3.512337  # 5 events 350+370/pb strange mumu bg
    #_p = 2.313430  # 5 events fake 600/pb
    #_p = 2.645160  # 5 events 350+370/pb strange mumu bg
    _p = observed_test_stat

    _count = 0
    for l in _limits:
        #for n in l:
            #print n, '   ',
            
        #print

        if float(l[1])<_p:
            _count += 1

        _h.Fill(l[1])

    _nup = len(_limits) - _count
    _pval = _nup/float(len(_limits))
    print 'test_stat:', _p
    print 'p-value:', _pval
    _z = math.sqrt(2.0)*TMath.ErfInverse(1.0-2.0*_pval)
    print 'significance:', _z

    c = TCanvas("c", "c", 600, 400)
    c.SetFrameLineWidth(2)
    c.SetRightMargin(0.05)
    c.cd()

    _h.SetLineWidth(2)
    _h.SetLineColor(ROOT.kBlue+2)
    _h.SetFillStyle(1001)
    _h.SetFillColor(ROOT.kBlue-9)
    _h.Draw("")

    _l = TLine(_p,0,_p,_h.GetMaximum() * 0.9)
    _l.SetLineWidth(2)
    _l.SetLineColor(ROOT.kRed)
    _l.Draw("same")
    latex = TLatex()
    latex.SetNDC()
    latex.SetTextSize(0.045)
    latex.SetTextAlign(5) # align left
    #latex.DrawLatex(0.45,0.5,"#font[62]{M = 954.3 GeV}")
    latex.DrawLatex(0.45,0.45,"#font[62]{p-value = "+'%.2f'%_pval+"}")
    latex.DrawLatex(0.45,0.4,"#font[62]{Z = "+'%.1f'%_z+"}")
    #latex.DrawLatex(0.45,0.5,"#font[62]{M = 962.2 GeV}")
    #latex.DrawLatex(0.45,0.45,"#font[62]{p-value = 0.0068}")
    #latex.DrawLatex(0.45,0.4,"#font[62]{Z = 2.5}")

    # x axis
    latex.SetTextAlign(31) # align left
    latex.DrawLatex(0.95,0.01, 'test statistic')

    # y axis
    latex.SetTextAlign(31) # align left
    latex.SetTextAngle(90)
    latex.DrawLatex(0.04,0.95, 'number of pseud-experiments')

    # legend
    xLegend = 0.65
    yLegend = 0.75
    legendWidth = 0.3
    legendHeight = 0.15
    tLegend = TLegend(xLegend, yLegend,
                      xLegend + legendWidth, yLegend + legendHeight)
    tLegend.AddEntry(_h, "pseudo-experiments", 'l')
    tLegend.AddEntry(_l, "observed data", 'l')
    tLegend.Draw()

    # CMS preliminary
    latex.SetTextAngle(0)
    #latex.DrawLatex(0.9,0.929,"#font[62]{CMS preliminary, #int #font[42]{L}dt = "+'%.1ffb^{-1} (#mu^{+}#mu^{-})'%intLumiMumu+"}")
    #latex.DrawLatex(0.9,0.929,"#font[62]{CMS preliminary, #int #font[42]{L}dt = "+'%.1ffb^{-1} (ee)'%intLumiEe+"}")
    latex.DrawLatex(0.9,0.929,"#font[62]{CMS preliminary, #int #font[42]{L}dt = "+'%.1ffb^{-1} (#mu^{+}#mu^{-}), %.1fpb^{-1} (ee)'%(intLumiMumu,intLumiEe)+"}")

    #latex.DrawLatex(0.9,0.929,"#font[62]{CMS, #int #font[42]{L}dt = "+'%.1ffb^{-1} (#mu^{+}#mu^{-})'%intLumiMumu+"}")
    #latex.DrawLatex(0.9,0.929,"#font[62]{CMS, #int #font[42]{L}dt = "+'%.1ffb^{-1} (ee)'%intLumiEe+"}")
    #latex.DrawLatex(0.9,0.929,"#font[62]{CMS, #int #font[42]{L}dt = "+'%.1ffb^{-1} (#mu^{+}#mu^{-}), %.1fpb^{-1} (ee)'%(intLumiMumu,intLumiEe)+"}")
    
#    _p = 2.67
#    _count = 0
#    for l in _limits:
#
#        if float(l[1])<_p:
#            _count += 1
#
#        _h.Fill(l[1])
#
#    _nup = len(_limits) - _count
#    _pval = _nup/float(len(_limits))
#    print 'test_stat:', _p
#    print 'p-value:', _pval
#    print 'significance:', math.sqrt(2.0)*TMath.ErfInverse(1.0-2.0*_pval)
#    _l2 = TLine(_p,0,_p,1000)
#    _l2.Draw("same")
#    latex.DrawLatex(0.15,0.3,"#font[62]{M = 942.2 GeV}")
#    latex.DrawLatex(0.15,0.25,"#font[62]{p-value = 0.05}")
    
    c.SaveAs("lee.pdf")
    
