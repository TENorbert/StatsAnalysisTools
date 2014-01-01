#!/usr/bin/env python
##########################################################################
# hist_factory.py
#
# Module for fitting model to data
#
# Usage:
#
#     ./tprime.py --fit
#
# Author: Gena Kukartsev, January 2012
#
#########################################################################

import os
import ROOT


def Fit():

    legend = '[tprime_fit.Fit]:'

    # workspace file for doing the fit
    #wsFileName = 'comb_reb20jan11/combined_test_model.root' #'combined_reb/combined_test_model.root'
    wsFileName = 'mujets_reb20jan11/mujets_test_model.root' 
    massorht = 'mass'
    folder   = 'mass20jan11'
    # workspace object name in the file
    emu    = 'mujets'
    wsName = emu
    #    wsName = emu+'jets'
    
    # workspace file for plotting (mass, mass)
    wsPlotFileName = emu+'_'+folder+'/'+wsName+'_test_model.root'
    
    # data needs to be read in separately as the errors aren't propagated properly or not printed properly in the model?
    datafile = emu+'_'+folder+'/hf_input_'+wsName+'_test.root'

    
    wsFile = ROOT.TFile(wsFileName, 'r')
    
    #    ws = wsFile.Get('combined') #combined if fitting to combined.
    ws = wsFile.Get(wsName) 
    
    
    #get model config from workspace
    modelConfig = ws.genobj('ModelConfig')
    
    ws.Print()
    modelConfig.Print()
    
    obsData = ws.data('obsData')
    
    #get full model pdf
    model = modelConfig.GetPdf()
    model.Print()
    
    #do fit of the full model to data
    fitResult = model.fitTo(obsData)
    
    #observable
    _obs = ws.var('obs_x_'+wsName)
    
    #three main event categories, all systematics included
    _signal = ws.function('L_x_signal_'+wsName+'_overallSyst_x_HistSyst')
    _top = ws.function('L_x_top_'+wsName+'_overallSyst_x_HistSyst')
    _ewk = ws.function('L_x_ewk_'+wsName+'_overallSyst_x_HistSyst')
    
    print 'signal integral test:', _signal.createIntegral(ROOT.RooArgSet(_obs)).getVal()
    print 'top integral test:', _top.createIntegral(ROOT.RooArgSet(_obs)).getVal()
    print 'ewk integral test:', _ewk.createIntegral(ROOT.RooArgSet(_obs)).getVal()
    
    
    # now get the workspace that we want to plot
    wsPlotFile = ROOT.TFile(wsPlotFileName, 'r')
    wsPlot = wsPlotFile.Get(wsName)
    modelConfigPlot = wsPlot.genobj('ModelConfig')

    # parameters from the fit
    xsecFit = ws.var('xsec').getVal()
    print 'xsec fitted value is :',xsecFit
    nuisFit = modelConfig.GetNuisanceParameters()

    # set parameters for the plot model
    # to be the same as in the fit

    # POI
    wsPlot.var('xsec').setVal(xsecFit)

    # all nuisance parameters
    _iter = nuisFit.createIterator()
    _key = _iter.Next()
    #_iter.Reset()
    while _key:
        _name = _key.GetName()
        print _name
        if ("eff_mu" in _name) and ( wsName == "ejets"):
            _key = _iter.Next()
            _name = _key.GetName()
        if ("eff_e" in _name ) and (wsName == "mujets"):
            _key = _iter.Next()
            _name = _iter.GetName()
        print _name
        wsPlot.var(_name).setVal( ws.var(_name).getVal() )
        print _name, wsPlot.var(_name).getVal(),"<--->",  ws.var(_name).getVal()
        _key = _iter.Next()


    # get plotting components
    _obs = wsPlot.var('obs_x_'+wsName)
    _frame = _obs.frame()
    obsData = wsPlot.data('obsData')

    obsData.plotOn(_frame,
                   ROOT.RooFit.MarkerSize(0.4),
                   ROOT.RooFit.MarkerStyle(ROOT.kFullCircle)#,
                   #ROOT.RooFit.XErrorSize(0)
                   #ROOT.RooFit.DrawOption("")
                   )    
    # full model
    model = modelConfigPlot.GetPdf()

    # signal
    _signalPlot = wsPlot.function('L_x_signal_'+wsName+'_overallSyst_x_HistSyst')

    # top background
    _topPlot = wsPlot.function('L_x_top_'+wsName+'_overallSyst_x_HistSyst')

    #ewk background
    _ewkPlot = wsPlot.function('L_x_ewk_'+wsName+'_overallSyst_x_HistSyst')

    print 'signal integral test2 - of plotted dists :', _signalPlot.createIntegral(ROOT.RooArgSet(_obs)).getVal()
    print 'top integral test2 - of plotted dists :', _topPlot.createIntegral(ROOT.RooArgSet(_obs)).getVal()
    print 'ewk integral test2 - of plotted dists :', _ewkPlot.createIntegral(ROOT.RooArgSet(_obs)).getVal()

    """
    model.plotOn(_frame,
                 ROOT.RooFit.LineColor(ROOT.kBlack),
                 ROOT.RooFit.FillColor(8),
                 ROOT.RooFit.LineStyle(1),
                 ROOT.RooFit.LineWidth(1),
                 ROOT.RooFit.FillStyle(1001),
                 ROOT.RooFit.DrawOption("L && F"),
                 ROOT.RooFit.MoveToBack()
                 )

    _topPlot.plotOn(_frame,
                ROOT.RooFit.LineColor(8),
                ROOT.RooFit.LineWidth(2),
                ROOT.RooFit.FillColor(8),
                ROOT.RooFit.DrawOption("F"),
                ROOT.RooFit.FillStyle(1001),
                ROOT.RooFit.MoveToBack()
                )

    _ewkPlot.plotOn(_frame,
                ROOT.RooFit.LineColor(ROOT.kBlack),
                ROOT.RooFit.LineWidth(2),
                ROOT.RooFit.FillColor(ROOT.kBlue),
                ROOT.RooFit.DrawOption("FL"),
                ROOT.RooFit.FillStyle(1001))

    _signalPlot.plotOn(_frame,
                ROOT.RooFit.LineColor(ROOT.kBlack),
                ROOT.RooFit.LineWidth(2),
                ROOT.RooFit.LineStyle(2))#,
                   #                ROOT.RooFit.Normalization(1,0) )

    _frame.GetXaxis().SetTitle("")
    _frame.GetYaxis().SetTitle("Events")
    _frame.SetTitle("CMS preliminary #int L = 4.7 fb^{-1}, e+jets 575GeV")
    _frame.Draw()
    raw_input("")
    """
    
    #Mike added plotting of histograms
    c1 = ROOT.TCanvas("c1","My Canvas", 0,0,600,400)
    c1.Clear()
    c1.SetFillColor(0)
    c1.SetLeftMargin(0.1)
    c1.SetRightMargin(0.04)
    c1.SetTopMargin(0.09)
    c1.SetBottomMargin(0.13)
    c1.SetFrameBorderMode(0)

    if wsName =="ejets":
        lumi = "4.683"
    else:
        lumi = "4.601"
    text = ROOT.TLatex(0.6,.94,"CMS preliminary #int L = "+lumi+" fb^{-1}")
    text.SetTextSize(0.04)
    text.SetNDC()
    
              
    topHisto = _topPlot.createHistogram("topHisto",_obs)
    signalHisto = _signalPlot.createHistogram("signalHisto",_obs)
    signalHisto2 = _signalPlot.createHistogram("signalHisto",_obs)
    ewkHisto  = _ewkPlot.createHistogram("ewkHisto",_obs)
    #    dataHisto = obsData.createHistogram("dataHisto",obs)

    topHisto.SetFillColor(8)
    ewkHisto.SetFillColor(ROOT.kBlue)
    signalHisto.SetFillColor(ROOT.kRed)
    signalHisto.SetLineStyle(ROOT.kDashed)
    SMHisto = ROOT.THStack("SMHisto","")

    if (massorht == "mass") or (massorht == "ht"):
        signal_scaling = _signal.createIntegral(ROOT.RooArgSet(_obs)).getVal()/_signalPlot.createIntegral(ROOT.RooArgSet(_obs)).getVal()
        top_scaling = _top.createIntegral(ROOT.RooArgSet(_obs)).getVal()/_topPlot.createIntegral(ROOT.RooArgSet(_obs)).getVal()
        ewk_scaling = _ewk.createIntegral(ROOT.RooArgSet(_obs)).getVal()/_ewkPlot.createIntegral(ROOT.RooArgSet(_obs)).getVal()
        topHisto.Scale(top_scaling)
        signalHisto.Scale(signal_scaling)
        signalHisto2.Scale(signal_scaling)
        ewkHisto.Scale(ewk_scaling)
        print "HACK rescaled by (sig,top,ewk) : ", signal_scaling,top_scaling,ewk_scaling
    
    SMHisto.Add(ewkHisto)
    SMHisto.Add(topHisto)
    SMHisto.Add(signalHisto)
    fData = ROOT.TFile(datafile, 'r')
    hData = fData.Get('DATA')
    if hData.GetMaximum() > SMHisto.GetMaximum:
        SMHisto.SetMaximum(hData.GetMaximum()*0.98)

    ROOT.gPad.SetLogy()    
    SMHisto.Draw()
    #hData.Draw("SAME")
    text.Draw("SAME")
    if massorht == "mass":
        SMHisto.GetXaxis().SetTitle("Mass [GeV]")
        SMHisto.GetYaxis().SetTitle("Events [/20GeV]")
    elif massorht == "ht":
        SMHisto.GetXaxis().SetTitle("H_{T} [GeV] Fitted Objects")
        SMHisto.GetYaxis().SetTitle("Events [/40GeV]")
    else:
        SMHisto.GetYaxis().SetTitle("Events")


    hData.GetYaxis().SetTitleOffset(1.)
    hData.GetXaxis().SetTitleSize(0.05)
    hData.GetYaxis().SetTitleSize(0.05)

    
    _legend = ROOT.TLegend(0.63,0.6,0.85,0.85)

    if wsName =="mujets":
        _legend.SetHeader("#mu+jets")
    else:
        _legend.SetHeader("e+jets")

    



    hData.SetLineColor(ROOT.kBlack)
    hData.SetMarkerStyle(20)
    hData.SetMarkerSize(0.5)
    
    #_legend.AddEntry(hData,"Data ","LPE")
    _legend.AddEntry(ewkHisto,"Ewk ","F")
    _legend.AddEntry(topHisto,"t#bar{t} ","F")
    _legend.AddEntry(signalHisto,"t'#bar{t'} (575GeV)","F")
    _legend.SetTextSize(0.035)
    _legend.SetFillColor(0)
    _legend.SetLineColor(0)
    

    signalHisto2.SetFillColor(0)
    signalHisto2.Scale(200)
    signalHisto2.SetLineColor(ROOT.kBlack)
    signalHisto2.SetLineStyle(ROOT.kDashed)
    _legend.AddEntry(signalHisto2,"x200 t'#bar{t'} (575GeV)","f")
    
    signalHisto2.Draw("SAME")
    _legend.Draw("SAME")
    
    #hData = ROOT.TH1D('hdata', 'hdata',
    #                  signalHisto.GetNbinsX(), signalHisto.GetBinLowEdge(1), signalHisto.GetBinLowEdge(signalHisto.GetNbinsX()+1)) #needs to be bin width of histo
    #hData = obsData.fillHistogram(hData,ROOT.RooArgList(_obs) )
    #hData.Sumw2()
    #hData.Draw("SAME")
    ROOT.c1.SaveAs(wsName+"_"+massorht+"__fitted_20jan12_logy.pdf")


    #create pull histogram
    c2 = ROOT.TCanvas("c2","My Canvas", 0,0,600,400)
    c2.Clear()
    c2.SetFillColor(0)
    c2.SetLeftMargin(0.1)
    c2.SetRightMargin(0.04)
    c2.SetTopMargin(0.09)
    c2.SetBottomMargin(0.13)
    c2.SetFrameBorderMode(0)
    c2.cd()

    pullHisto = ROOT.TH1D("pullhisto","pull : (data-MC)/MC",1000,-2,2.)
    ratioHisto = ROOT.TH1D("ratiohisto","ratio : data/MC",hData.GetNbinsX(),0,hData.GetBinLowEdge(hData.GetNbinsX()))

    for bins in range(0,hData.GetNbinsX()):
        _mcbg = ewkHisto.GetBinContent(bins) + topHisto.GetBinContent(bins) + signalHisto.GetBinContent(bins)
        if _mcbg != 0:
            pullHisto.Fill((hData.GetBinContent(bins) - _mcbg)/_mcbg)
            ratioHisto.Fill(hData.GetBinCenter(bins),hData.GetBinContent(bins)/_mcbg)
            
        
    pullHisto.Draw()
    ROOT.c2.SaveAs(wsName+"_"+massorht+"__pull_20jan12.pdf")

    c3 = ROOT.TCanvas("c3","My Canvas", 0,0,600,400)
    c3.Clear()
    c3.SetFillColor(0)
    c3.SetLeftMargin(0.1)
    c3.SetRightMargin(0.04)
    c3.SetTopMargin(0.09)
    c3.SetBottomMargin(0.13)
    c3.SetFrameBorderMode(0)
    c3.cd()

    ratioHisto.SetMinimum(0.001)
    ROOT.gPad.SetLogy()
    ratioHisto.Draw()
    ROOT.c3.SaveAs(wsName+"_"+massorht+"__ratio_20jan12.pdf")

    raw_input('hit enter to continue...')
    


