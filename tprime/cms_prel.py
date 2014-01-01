#!/usr/bin/env python
#########################################################################
#
# cms_prel.py
#
# Write CMS preliminary and other standard disclaimers on plots
#
# Usage:
#        import cms_prel
#        cms_prel.cmsPrel()
#
#
# Author: Gena Kukartsev, December 2011
#
#########################################################################
from __future__ import division

import ROOT


def SetPasStyle():

    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetNdivisions(405, "X")
    ROOT.gStyle.SetNdivisions(405, "Y")
    ROOT.gStyle.SetLabelSize(0.045, "X")
    ROOT.gStyle.SetLabelSize(0.045, "Y")
    ROOT.gStyle.SetLabelSize(0.045, "Z")
    ROOT.gStyle.SetLabelFont(22, "X")
    ROOT.gStyle.SetLabelFont(22, "Y")
    ROOT.gStyle.SetLabelFont(22, "Z")
    ROOT.gStyle.SetPalette( 1 )
    
    return



def SetPrlStyle():

    _leftMargin = 0.2;
    _rightMargin = 0.045;
    _topMargin = 0.09;
    _bottomMargin = 0.18;
    _labelSize = 0.07;
    _titleSize = 0.07;
    _textSize  = 0.065;
    _latextextSize = 0.08;
    _xtitoffset=1.2;
    _ytitoffset=1.4;
    _logy = 1
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetCanvasBorderMode(0)
    ROOT.gStyle.SetCanvasColor(ROOT.kWhite)
    ROOT.gStyle.SetPadBorderMode(0)
    ROOT.gStyle.SetPadColor(ROOT.kWhite)
    ROOT.gStyle.SetTitleXOffset(_xtitoffset)
    ROOT.gStyle.SetTitleYOffset(_ytitoffset)
    ROOT.gStyle.SetTitleSize(_titleSize, "XYZ")
    ROOT.gStyle.SetLabelSize(_labelSize, "XYZ")
    ROOT.gStyle.SetTitleFont(22,"XYZ")
    ROOT.gStyle.SetLabelFont(22,"XYZ")
    #
    ROOT.gStyle.SetPadRightMargin(_rightMargin)
    ROOT.gStyle.SetPadBottomMargin(_bottomMargin)
    ROOT.gStyle.SetOptLogy(_logy)
    #
    ROOT.gStyle.SetNdivisions(405, "X")
    ROOT.gStyle.SetNdivisions(405, "Y")
    ROOT.gStyle.SetPalette( 1 )
    
    return



def Title(title = '', x=0.93, y=0.93, text_size = 0.05):
    latex = ROOT.TLatex()
    latex.SetNDC()
    latex.SetTextSize(text_size)
    latex.SetTextFont(22)

    latex.SetTextAlign(31) # align right
    #latex.SetTextAlign(5) # align left
    latex.DrawLatex(x, y, '#font[22]{'+title+'}')

    return



def XLabel(label = '', x=0.95, y=0.02, text_size = 0.05):
    latex = ROOT.TLatex()
    latex.SetNDC()
    latex.SetTextSize(text_size)
    latex.SetTextFont(22)

    latex.SetTextAlign(31) # align right
    #latex.SetTextAlign(5) # align left
    latex.DrawLatex(x, y, label)

    return



def YLabel(label = '', x=0.03, y=0.95, text_size = 0.05):
    latex = ROOT.TLatex()
    latex.SetNDC()
    latex.SetTextSize(text_size)
    latex.SetTextFont(22)

    latex.SetTextAlign(31) # align right
    #latex.SetTextAlign(5) # align left
    latex.SetTextAngle(90)
    latex.DrawLatex(x, y, label)

    return



def CmsPrel(x=0.10, y=0.93, text_size = 0.05):
    latex = ROOT.TLatex()
    latex.SetNDC()
    latex.SetTextSize(text_size)
    latex.SetTextFont(22)

    #latex.SetTextAlign(31) # align right
    latex.SetTextAlign(5) # align left
    latex.DrawLatex(x, y, "#font[22]{CMS preliminary}")

    return



def CmsSim(x=0.10, y=0.93, text_size = 0.05):
    latex = ROOT.TLatex()
    latex.SetNDC()
    latex.SetTextSize(text_size)
    latex.SetTextFont(22)

    #latex.SetTextAlign(31) # align right
    latex.SetTextAlign(5) # align left
    latex.DrawLatex(x, y, "#font[22]{CMS simulation}")

    return



def Cms(x=0.10, y=0.93, text_size = 0.05):
    latex = ROOT.TLatex()
    latex.SetNDC()
    latex.SetTextSize(text_size)
    latex.SetTextFont(22)

    #latex.SetTextAlign(31) # align right
    latex.SetTextAlign(5) # align left
    latex.DrawLatex(x, y, "#font[22]{CMS}")

    return



def Cms7Tev(x=0.11, y=0.85, text_size = 0.05):
    latex = ROOT.TLatex()
    latex.SetNDC()
    latex.SetTextSize(text_size)
    latex.SetTextFont(22)

    #latex.SetTextAlign(31) # align right
    latex.SetTextAlign(5) # align left
    latex.DrawLatex(x, y, "#font[22]{#sqrt{s}=7 TeV}")

    return



def Legend(legend, x=0.15, y=0.22, text_size = 0.05):
    latex = ROOT.TLatex()
    latex.SetNDC()
    latex.SetTextSize(text_size)
    latex.SetTextFont(22)

    #latex.SetTextAlign(31) # align right
    latex.SetTextAlign(5) # align left
    latex.DrawLatex(x, y, '#font[22]{'+legend+'}')

    return


