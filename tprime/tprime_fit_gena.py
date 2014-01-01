#!/usr/bin/env python
#########################################################################
#
# hist_factory.py
#
# Module for fitting model to data
#
# Usage:
#     
#     do fit
#     ./tprime.py --fit fit 
#     
#     make fit plots of templates
#     ./tprime.py --fit plots
#
# Author: Gena Kukartsev, January 2012
#
#########################################################################

import os
import math
import copy
import ROOT

import cms_prel as cstyle


def FitSimple(wsFileName, wsName, mass):
    #
    # Fit the HistFactory model to data
    # Save fit results to a file
    # "mass" is for book keeping only
    #

    legend = '[tprime_fit.FitSimple]:'

    # workspace file for doing the fit
    #wsFileName = 'comb_test/mujets_test_model.root'
    #wsFileName = 'comb_test/combined_test_model.root'
    #wsFileName = 'comb_test/combined_test_model.root'
    
    # workspace object name in the file
    #wsName = 'ejets'
    #wsName = 'mujets'
    #wsName = 'combined'

    if wsFileName:
        wsFile = ROOT.TFile(wsFileName, 'r')
    else:
        print legend, 'workspace file not specified, exiting'
        return

    if wsName and wsFile:
        ws = wsFile.Get(wsName)
    else:
        print legend, 'workspace file not found or workspace name incorrect, exiting'
        return
    
    #get model config from workspace
    modelConfig = ws.genobj('ModelConfig')
    
    ws.Print()
    modelConfig.Print()
    
    obsData = ws.data('obsData')
    
    #get full model pdf
    model = modelConfig.GetPdf()
    model.Print()

    # optionally set and fix some parameters
    #
    # values from mu fit
    #ws.var('Lumi').setVal(0.999561)
    #ws.var('alpha_jes').setVal(0.495713)
    #ws.var('alpha_lepton_eff_e').setVal(0.28623)
    #ws.var('alpha_lepton_eff_mu').setVal(-0.314082)
    #ws.var('alpha_match').setVal(0.203847)
    #ws.var('alpha_norm_ewk_syst').setVal(-0.609317)
    #ws.var('alpha_norm_top_syst').setVal(0.0101145)
    #ws.var('xsec').setVal(0.0755228)

    
    # tweak the minimizer
    minName  = 'Minuit'
    algoName = 'migradimproved'
    #ROOT::Math::MinimizerOptions::SetDefaultPrecision(0.000001)
    #ROOT::Math::MinimizerOptions::SetDefaultPrintLevel(1)
    #ROOT.Math.MinimizerOptions.SetDefaultMinimizer(minName, algoName)
    #ROOT.Math.MinimizerOptions.SetDefaultTolerance(0.1);
    #ROOT.Math.MinimizerOptions.SetDefaultMaxIterations(3);
    print legend, 'minimizer type:', ROOT.Math.MinimizerOptions.DefaultMinimizerType()
    print legend, 'minimizer algo:', ROOT.Math.MinimizerOptions.DefaultMinimizerAlgo()
    print legend, 'minimizer precision:', ROOT.Math.MinimizerOptions.DefaultPrecision()
    print legend, 'minimizer tolerance:', ROOT.Math.MinimizerOptions.DefaultTolerance()

    #do fit of the full model to data
    #fitResult = model.fitTo(obsData,
    #                        ROOT.RooFit.Hesse(ROOT.kFALSE),
    #                        ROOT.RooFit.Minimizer(minName, 'Simplex'))
    fitResult = model.fitTo(obsData,
                            ROOT.RooFit.Hesse(ROOT.kFALSE),
                            ROOT.RooFit.Minimizer(minName, algoName))
    
    #observable
    _obs = ws.var('obs_x_'+wsName)

    # dictionary with two strings to return:
    #  - parameter names
    #  - parameter values
    sResult = {}
    sNames = ''
    sValues = ''
    
    # parameters from the fit
    xsecFit = ws.var('xsec').getVal()
    print 'xsec fitted value is :',xsecFit
    nuisFit = modelConfig.GetNuisanceParameters()
    sNames +='xsec   '
    sValues += str(xsecFit)+'   '

    # all nuisance parameters
    _iter = nuisFit.createIterator()
    _key = _iter.Next()
    #iter.Reset()
    while _key:
        _name = _key.GetName()
        _value = ws.var(_name).getVal()
        print _name, _value
        sNames += _name+'   '
        sValues += str(_value)+'   '
        _key = _iter.Next()

    sResult['names'] = sNames
    sResult['values'] = sValues

    print sResult['names']
    print sResult['values']

    return sResult



def MakeFitPlots():

    legend = '[tprime_fit.MakeFitPlots]:'

    _mass = 550 # for legend only
    _ymin = 2
    _sigXSec = 0.171 # 550
    _sigScale = 50.0
    _textsize = 0.065

    _xleg  = 0.55
    _xleg1 = 0.75
    _yleg  = 0.60
    _yleg1 = 0.90
    
    color1 = ROOT.TColor(231,0.4,1,0.4,"")
    color2 = ROOT.TColor(232,0,0.7,0,"")
    col1 = 231
    col2 = 232

    cstyle.SetPrlStyle()

    # workspace file for doing the fit
    wsFileName = 'templates_test/ejets_test_model.root'
    #wsFileName = 'templates_test/mujets_test_model.root'
    #wsFileName = 'templates_test/combined_test_model.root'
    #wsFileName = 'templ_comb_08may2012v1/combined_test_model_475.root'
    
    # workspace file for plotting (mass, ht)
    wsPlotFileName = wsFileName
    #wsPlotFileName = 'templ_mass_08may2012v1/combined_test_model_475.root'
    #wsPlotFileName = 'templ_ht_08may2012v1/combined_test_model_475.root'

    # workspace object name in the file
    #wsName = 'mujets'
    wsName = 'ejets'
    #wsName = 'combined'

    channel = 'ejets' # needed for combination. Used to be channel=wsName
    #channel = 'mujets' # needed for combination. Used to be channel=wsName
    massorht = 'rebinned' #specify 'mass','ht' or 'rebinned' which is after ganging

    _bg_only_model = True
    
    wsFile = ROOT.TFile(wsFileName, 'r')
    
    ws = wsFile.Get(wsName)
    
    #get model config from workspace
    modelConfig = ws.genobj('ModelConfig')
    
    ws.Print()
    modelConfig.Print()
    
    obsData = ws.data('obsData')
    
    #get full model pdf
    model = modelConfig.GetPdf()
    model.Print()
    
    if _bg_only_model:
        #do fit of the BG-only model to data
        modelConfig.GetParametersOfInterest().first().setVal(0.0)
        modelConfig.GetParametersOfInterest().first().setConstant(1)

    fitResult = model.fitTo(obsData)
    
    #observable
    #_obs = ws.var('obs_x_'+wsName)
    _obs = ws.var('obs_x_'+channel)
    
    #three main event categories, all systematics included
    _signal = ws.function('L_x_signal_'+channel+'_overallSyst_x_HistSyst')
    _top = ws.function('L_x_top_'+channel+'_overallSyst_x_HistSyst')
    _ewk = ws.function('L_x_ewk_'+channel+'_overallSyst_x_HistSyst')
    
    _signal_integral = _signal.createIntegral(ROOT.RooArgSet(_obs)).getVal()
    _top_integral = _top.createIntegral(ROOT.RooArgSet(_obs)).getVal()
    _ewk_integral =  _ewk.createIntegral(ROOT.RooArgSet(_obs)).getVal()
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

    # get expected theoretical signal
    #modelConfig.GetParametersOfInterest().first().setVal(_sigXSec)
    ws.var('xsec').setVal(_sigXSec)
    _signal_expected = _signal.createIntegral(ROOT.RooArgSet(_obs)).getVal()
    modelConfig.GetParametersOfInterest().first().setVal(xsecFit)

    # set parameters for the plot model
    # to be the same as in the fit

    # POI
    if _bg_only_model:
        wsPlot.var('xsec').setVal(_sigXSec)
    else:
        wsPlot.var('xsec').setVal(xsecFit)

    # all nuisance parameters
    _iter = nuisFit.createIterator()
    _key = _iter.Next()
    #iter.Reset()
    while _key:
        _name = _key.GetName()
        wsPlot.var(_name).setVal( ws.var(_name).getVal() )
        print _name, wsPlot.var(_name).getVal()
        _key = _iter.Next()


    # get plotting components
    _obs = wsPlot.var('obs_x_'+channel)
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
    _signal = wsPlot.function('L_x_signal_'+channel+'_overallSyst_x_HistSyst')

    # top background
    _top = wsPlot.function('L_x_top_'+channel+'_overallSyst_x_HistSyst')

    #ewk background
    _ewk = wsPlot.function('L_x_ewk_'+channel+'_overallSyst_x_HistSyst')

    
#    model.plotOn(_frame,
#                 ROOT.RooFit.LineColor(ROOT.kBlack),
#                 ROOT.RooFit.FillColor(8),
#                 ROOT.RooFit.LineStyle(1),
#                 ROOT.RooFit.LineWidth(1),
#                 ROOT.RooFit.FillStyle(1001),
#                 ROOT.RooFit.DrawOption("L && F"),
#                 ROOT.RooFit.MoveToBack()
#                 )

    _top.plotOn(_frame,
                ROOT.RooFit.LineColor(8),
                ROOT.RooFit.LineWidth(2),
                ROOT.RooFit.FillColor(8),
                ROOT.RooFit.DrawOption("F"),
                ROOT.RooFit.FillStyle(1001),
                ROOT.RooFit.MoveToBack()
                )

    _ewk.plotOn(_frame,
                ROOT.RooFit.LineColor(ROOT.kBlack),
                ROOT.RooFit.LineWidth(2),
                ROOT.RooFit.FillColor(ROOT.kBlue),
                ROOT.RooFit.DrawOption("FL"),
                ROOT.RooFit.FillStyle(1001))

    _signal.plotOn(_frame,
                ROOT.RooFit.LineColor(ROOT.kBlack),
                ROOT.RooFit.LineWidth(2),
                ROOT.RooFit.LineStyle(2),
                ROOT.RooFit.Normalization(1,0) )


    #Mike added plotting of histograms
    c1 = ROOT.TCanvas("c1","My Canvas", 0,0,600,400)
    c1.Clear()
    c1.SetFillColor(0)
    c1.SetLeftMargin(0.1)
    c1.SetRightMargin(0.04)
    c1.SetTopMargin(0.09)
    c1.SetBottomMargin(0.13)
    c1.SetFrameBorderMode(0)


    # set arbitrary signal rate for plotting
    modelConfigPlot.GetParametersOfInterest().first().setVal(_sigXSec)
              
    topHisto = _top.createHistogram("topHisto",_obs)
    signalHisto = _signal.createHistogram("signalHisto",_obs)
    print legend, 'DEBUG signal integral', signalHisto.Integral()
    ewkHisto  = _ewk.createHistogram("ewkHisto",_obs)
    #    dataHisto = obsData.createHistogram("dataHisto",_obs)

    topHisto.SetLineColor(ROOT.kBlack)
    ewkHisto.SetLineColor(ROOT.kBlack)
    #signalHisto.SetFillColor(ROOT.kRed)
    signalHisto.SetLineStyle(ROOT.kDashed)
    signalHisto.SetLineWidth(2)
    signalHisto.SetLineColor(ROOT.kBlack)
    SMHisto = ROOT.THStack("SMHisto","")

    _nbins = signalHisto.GetNbinsX()

    topTempl = copy.deepcopy(topHisto)
    ewkTempl = copy.deepcopy(ewkHisto)
    signalTempl = copy.deepcopy(signalHisto)

    # reverse bins: i -> N-i+1
    for bin in range(1, _nbins+1):
        if (massorht == 'rebinned'):
            topTempl.SetBinContent(bin, topHisto.GetBinContent(_nbins-bin+1))
            ewkTempl.SetBinContent(bin, ewkHisto.GetBinContent(_nbins-bin+1))
            signalTempl.SetBinContent(bin, signalHisto.GetBinContent(_nbins-bin+1))

    topTempl.SetFillColor(col1)
    ewkTempl.SetFillColor(col2)
    if _bg_only_model==False:
        signalTempl.SetFillColor(ROOT.kRed)

    # rescale to correct norm
    topTempl.Scale(_top_integral/topTempl.Integral())
    ewkTempl.Scale(_ewk_integral/ewkTempl.Integral())
    if _bg_only_model==False:
        signalTempl.Scale(_signal_integral/signalTempl.Integral())
    else:
        signalTempl.Scale(_signal_expected/signalTempl.Integral())
    
    SMHisto.Add(ewkTempl)
    SMHisto.Add(topTempl)
    if _bg_only_model==False:
        SMHisto.Add(signalTempl)

    _ymax = SMHisto.GetMaximum()

    #c1.DrawFrame(0,_ymin,_nbins,_ymax*5)
    c1.DrawFrame(signalTempl.GetXaxis().GetXmin(),_ymin,signalTempl.GetXaxis().GetXmax(),_ymax*5)

    SMHisto.Draw('same')
    #text.Draw("SAME")
    if massorht == "mass":
        SMHisto.GetXaxis().SetTitle("Mass [GeV]")
        SMHisto.GetYaxis().SetTitle("Events [/20GeV]")
    elif massorht == "ht":
        SMHisto.GetXaxis().SetTitle("H_{T} [GeV] Fitted Objects")
        SMHisto.GetYaxis().SetTitle("Events [/40GeV]")
    SMHisto.GetYaxis().SetTitleOffset(1.)
    SMHisto.GetXaxis().SetTitleSize(0.05)
    SMHisto.GetYaxis().SetTitleSize(0.05)

    
    
    if (massorht == 'rebinned'):
        hData = ROOT.TH1F('hdata', 'hdata',
                          _nbins, 0, _nbins)
    else:
        #hData = copy.deepcopy(signalHisto)
        hData = ROOT.TH1F('hdata', 'hdata',
                          _nbins, signalTempl.GetXaxis().GetXmin(), signalTempl.GetXaxis().GetXmax())
        hData.Clear()

    if wsName=='combined':
        hData = obsData.fillHistogram(hData,ROOT.RooArgList(_obs),'channelCat==channelCat::'+channel )
    else:
        hData = obsData.fillHistogram(hData,ROOT.RooArgList(_obs) )
        
    _obs.Print()
    obsData.Print()
    print legend, 'data integral', hData.Integral()
    print legend, 'top integral', topTempl.Integral()
    print legend, 'EWK integral', ewkTempl.Integral()
    print legend, 'signal integral', signalTempl.Integral()

    hData.Sumw2() # for whatever reason this works incorrectly
    #print 'bin 36', hData.GetBinContent(36), '+/-', hData.GetBinError(36)

    hDataRev = copy.deepcopy(hData)

    # reverse bin order
    for bin in range(1, _nbins+1):
        if (massorht == 'rebinned'):
            hDataRev.SetBinContent(bin, hData.GetBinContent(_nbins-bin+1))
            # fix for incorrectly computed errors
            hDataRev.SetBinError(bin, math.sqrt(hData.GetBinContent(_nbins-bin+1)))
        else:
            # fix for incorrectly computed errors
            hDataRev.SetBinError(bin, math.sqrt(hData.GetBinContent(bin)))

    hDataRev.SetMarkerStyle(20)
    hDataRev.SetMarkerSize(0.7)
    hDataRev.Draw("SAME:e")

    # plot signal template, magnified
    if _bg_only_model:
        signalTempl.Scale(_sigScale)
    signalTempl.Draw('same')

    c1.RedrawAxis()
           
    _legend = ROOT.TLegend(_xleg, _yleg, _xleg1, _yleg1)

    #if wsName =="mujets":
        #_legend.SetHeader("#mu+jets")
    #else:
        #_legend.SetHeader("e+jets")
        
    _legend.AddEntry(hDataRev,"Data ","p")
    _legend.AddEntry(topTempl,"t#bar{t} ","F")
    _legend.AddEntry(ewkTempl,"Other Bkg ","F")
    if _bg_only_model:
        #_legend.AddEntry(signalTempl,"t'#bar{t'} (550 GeV) #times 50","lp")
        #_legend.AddEntry(signalTempl,'t\'#bar{t\'} #times'+str(int(_sigScale)),"lp")
        _legend.AddEntry(signalTempl,'t\'#bar{t\'} ('+str(_mass)+' GeV) #times'+str(int(_sigScale)),"lp")
    else:
        _legend.AddEntry(signalTempl, 't\'#bar{t\'}', "f")
    _legend.SetTextSize(_textsize)
    _legend.SetTextAlign(11)
    _legend.SetFillColor(0)
    _legend.SetLineColor(0)
    _legend.Draw("SAME")

    cstyle.Cms(0.10, 0.92, text_size = _textsize)
    cstyle.Cms7Tev(0.30, 0.92, text_size = _textsize)
    if massorht=='rebinned':
        cstyle.XLabel('s/b Rank', text_size = _textsize)
    elif massorht=='mass':
        cstyle.XLabel('m_{t\'}', text_size = _textsize)
    elif massorht=='ht':
        cstyle.XLabel('H_{T}', text_size = _textsize)
    cstyle.YLabel('Events/Bin', text_size = _textsize)
    if channel=='ejets':
        cstyle.Legend('e+jets', text_size = _textsize)
        cstyle.Title('L=5.0fb^{-1}', text_size = _textsize)
    else:
        cstyle.Legend('#mu+jets', text_size = _textsize)
        cstyle.Title('L=4.9fb^{-1}', text_size = _textsize)

    if _bg_only_model:
        ROOT.c1.SaveAs(wsName+'_'+channel+'_'+massorht+'_bg_only__fitted.pdf')
    else:
        ROOT.c1.SaveAs(wsName+'_'+channel+'_'+massorht+'_fitted.pdf')

    raw_input('hit enter to continue...')
    


#-----> main module interface
def Fit(options):
    if options.fit == 'plots':
        MakeFitPlots()
    elif options.fit == 'fit':
        FitSimple(options.in_file[0],
                  options.name,
                  options.mass)
