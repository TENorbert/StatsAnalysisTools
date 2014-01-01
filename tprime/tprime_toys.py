#!/usr/bin/env python
#########################################################################
#
# tprime_toys.py
#
# Module for toy study of the tprime HistFactory model, checking for bias
#
# Usage:
#     
#     ./tprime.py --toys limit
#     ./tprime.py --toys pvalue
#     ./tprime.py --toys templates
#     ./tprime.py --toys fit
#
# Author: Gena Kukartsev, October 2011
#
#
# Notes: two bugs to report:
#          - 5.32.00-rc2 reads back workspace with errors
#          - ModelConfig cannot set WS, at least in PyROOT
#
#########################################################################

#--- configuration: set what you want to run here ----------------

#action = 'test'
#action = 'limit_2d_toys'
#action = 'generate_2d_toys'
#action = 'toys'
action = 'fit_pdfs'
#action = 'generate_1d_templates'
#action = 'generate_2d_templates'

# test
#pdfNames = ['sig_mass_pdf', 'bg_mass_pdf', 'sig_ht_sum_pdf', 'bg_ht_sum_pdf']
#varNames = ['mass', 'ht']
#varNBins = {'mass':45, 'ht':90}

# generate_1d_templates
#inPdfFileName   = 'data/toys/pdfs_ejets_smooth_1d_21nov2011.root'
#workspaceName   = 'toyWs'
#outTemplateFile = 'template_1d_generic.root'
#pdfNames  = ['bg_mass_pdf',   'sig_mass_pdf', 'bg_mass_pdf', 'bg_mass_pdf',         'sig_mass_pdf',       'bg_mass_pdf',       'bg_mass_pdf',           'sig_mass_pdf',         'bg_mass_pdf',         'bg_mass_pdf']
#histNames = [       'data', 'tprime450_mass',    'top_mass',    'ewk_mass', 'tprime450_mass_jesup',    'top_mass_jesup',    'ewk_mass_jesup', 'tprime450_mass_jesdown',    'top_mass_jesdown',    'ewk_mass_jesdown']
#vars      = [       'mass',           'mass',        'mass',        'mass',                 'mass',              'mass',              'mass',                   'mass',                'mass',                'mass']
#histLow   = [          100,              100,           100,           100,                    100,                 100,                 100,                      100,                   100,                   100]
#histHigh  = [         1000,             1000,          1000,          1000,                   1000,                1000,                1000,                     1000,                  1000,                  1000]
#histNBins = [           45,               45,            45,            45,                     45,                  45,                  45,                       45,                    45,                    45]
#nEvents   = [         4231,           100000,        100000,        100000,                 100000,              100000,              100000,                   100000,                100000,                100000]
#yields    = [           -1,            203.2,        3455.1,         776.0,                  207.7,              3767.4,               878.3,                    198.4,                3152.9,                 705.1]

# generate_2d_templates
#inPdfFileName   = 'data/toys/pdfs_ejets_smooth_1d_21nov2011.root'
#workspaceName   = 'toyWs'
#outTemplateFile = 'template_2d_generic.root'
#pdfNames   = [  'bg_mass_pdf',          'sig_mass_pdf',     'bg_mass_pdf',     'bg_mass_pdf',                'sig_mass_pdf',           'bg_mass_pdf',           'bg_mass_pdf',                  'sig_mass_pdf',             'bg_mass_pdf',             'bg_mass_pdf']
#pdfNamesY  = ['bg_ht_sum_pdf',        'sig_ht_sum_pdf',   'bg_ht_sum_pdf',   'bg_ht_sum_pdf',              'sig_ht_sum_pdf',         'bg_ht_sum_pdf',         'bg_ht_sum_pdf',                'sig_ht_sum_pdf',           'bg_ht_sum_pdf',           'bg_ht_sum_pdf']
#histNames  = [         'Data','TPrime450_ht35:fitMass','Top_ht35:fitMass','Ewk_ht35:fitMass','TPrime450_ht35:fitMass_JESup','Top_ht35:fitMass_JESup','Ewk_ht35:fitMass_JESup','TPrime450_ht35:fitMass_JESdown','Top_ht35:fitMass_JESdown','Ewk_ht35:fitMass_JESdown']
#vars       = [         'mass',                  'mass',            'mass',            'mass',                        'mass',                  'mass',                  'mass',                          'mass',                    'mass',                    'mass']
#varsY      = [           'ht',                    'ht',              'ht',              'ht',                          'ht',                    'ht',                    'ht',                            'ht',                      'ht',                      'ht']
#histNBins  = [             45,                      45,                45,                45,                            45,                      45,                      45,                              45,                        45,                        45]
#histNBinsY = [             90,                      90,                90,                90,                            90,                      90,                      90,                              90,                        90,                        90]
##nEvents    = [           4231,                 1000000,           1000000,           1000000,                       1000000,                 1000000,                 1000000,                         1000000,                   1000000,                   1000000]
#nEvents    = [           4231,                    15000,             50000,             10000,                         15000,                   50000,                   10000,                           15000,                     50000,                    10000]
#yields     = [             -1,                   203.2,            3455.1,             776.0,                         207.7,                  3767.4,                   878.3,                           198.4,                    3152.9,                     705.1]

# fit_pdfs
outFileName     = 'ws_pdf.root'
tprimeMass        = 450.0
templates2d       = 'data/toys/ejets_3.56_second_pseudo_2d_15nov2011v1.root'

# toys
##templates2d       = 'data/toys/ejets_3.56_second_pseudo_2d_15nov2011v1.root'
##templates2d       = 'data/toys/ejets_3.56_first_pseudo_2d_15nov2011v1.root'
##binMap            = 'data/toys/ejets_3.56_first_pseudo_merged_15nov2011v1.root.binmap'
#templates2d       = 'data/toys/templates_toy_smooth_2d_22nov2011v2.root'
#binMap            = 'data/toys/ejets_toy_smooth_merged_22nov2011v1.root.binmap'
#sbPdfName         = 'model_ejets'
#nSbToys           = 1
#tprimeMass        = 450.0
#tprimeXsec        =   0.662
##workspaceFileName = 'results_test/tprime_ejets_tprimeCrossSection_model_first_15nov2011v1.root'
#workspaceFileName = 'results_test/ejets_toy_smooth_2d_28nov2011v1_model.root'
#workspaceName     = 'ejets'
#obsName           = 'obs_x_ejets'
#poiName           = 'xsec'
#pdfFileName       = 'data/toys/pdfs_ejets_smooth_1d_21nov2011.root'
#pdfWsName         = 'toyWs'
#pdfNames = ['sig_mass_pdf', 'bg_mass_pdf', 'sig_ht_sum_pdf', 'bg_ht_sum_pdf']
#varNames = ['mass', 'ht']
#varNBins = {'mass':45, 'ht':90}

# generate_2d_toys
#toyMcFile = 'toys.root'
#pdfFileName       = 'data/toys/pdfs_ejets_smooth_1d_21nov2011.root'
#pdfWsName         = 'toyWs'
#pdfNames = ['sig_mass_pdf', 'bg_mass_pdf', 'sig_ht_sum_pdf', 'bg_ht_sum_pdf']
#varNames = ['mass', 'ht']
#varNBins = {'mass':45, 'ht':90}
#nSbToys           = 100
#sbPdfName         = 'model_ejets'
#tprimeMass        = 450.0
#tprimeXsec        =   0.662
#templates2d       = 'data/toys/templates_toy_smooth_2d_22nov2011v2.root'

# limit_2d_toys
#toyMcFile = 'toys.root'
#workspaceFileName = 'results_test/ejets_toy_smooth_2d_30nov2011v1_model.root'
#workspaceName     = 'ejets'
#obsName           = 'obs_x_ejets'
#poiName           = 'xsec'
#binMap            = 'data/toys/ejets_toy_smooth_merged_22nov2011v1.root.binmap'

#useInverter = True
useInverter = False

#--- end of configurtion -----------------------------------------



import sys
import copy
import string
import ROOT

from load_data import TprimeData
#from combine_bins import *
import combine_bins as cb

if useInverter:
    # load StandardHypoTestInvDemo
    ROOT.gROOT.ProcessLine('.L StandardHypoTestInvDemo.5.32.C+')
    invCalc = ROOT.RooStats.HypoTestInvTool()



workspace = None

class Toys:

    def __init__(self, data):

        self.legend = '[tprime_toys.datastore]:'
        self.data = data
        self.workspace = None
        self.random = ROOT.TRandom3(0)
        print self.legend, 'initialized'
        self.massPdf = None
        self.pdfs = {}
        self.vars = {} # cache some useful vars



    def loadBinMap(self, filename):
        #
        # load 2d->1d map from an ascii file created by combineBins routine
        #
        legend = '[Toys.loadBinMap]:'

        self.binMap = {}
        
        with open(filename) as mFile:

            lines=mFile.readlines()

            bin = 1
            for line in lines:
                _columns = line.strip().split(':')
                
                if _columns[0][0] == '#':
                    continue

                self.binMap[bin] = list(map(int,_columns[1].split()))
                bin += 1

        print legend, 'map for', len(self.binMap), 'bins loaded'

        return self.binMap




    def drawHist(self, name, style=''):
        #print data.hist_names
        self.data.hists[name].Draw(style)
        raw_input('press <enter> to continue...')



    def generateToy(self):
        # for now return the template
        return self.data.hists[name]

    
        
    def loadWorkspace(self, wsName, fileName):
        infile = ROOT.TFile(fileName, 'read')
        self.workspace = infile.Get(wsName)
        #self.workspace.Print()



    def LoadPdfs(self, inFileName, wsName, pdfNames, varNames):
        #
        # Load RooFit pdfs from file
        #

        legend = '[Toys.LoadPdfs]:'
        
        inFile = ROOT.TFile(inFileName, 'read')
        ws = inFile.Get(wsName)

        for name in pdfNames:
            self.pdfs[name] = ws.pdf(name).cloneTree()

        for name in varNames:
            self.vars[name] = ws.var(name).cloneTree()

        inFile.Close()

        print legend, self.pdfs

        return self.pdfs



    def Generate2D(self, varNameX, pdfNameX, varNameY, pdfNameY, nEvents):
        #
        # Load workspace, generate 2d histogram,
        # return
        #

        pdfs = self.pdfs
        vars = self.vars

        pdfX = pdfs[pdfNameX]
        pdfY = pdfs[pdfNameY]

        _pdf = ROOT.RooProdPdf('tmpPdf', 'tmpPdf', pdfX, pdfY)
        _spec = _pdf.prepareMultiGen(ROOT.RooArgSet( vars[varNameX], vars[varNameY] ),
                                     ROOT.RooFit.Name('ds_tmp'),
                                     ROOT.RooFit.NumEvents(nEvents),
                                     ROOT.RooFit.Extended(ROOT.kFALSE),
                                     ROOT.RooFit.Verbose(ROOT.kFALSE))
        
        _data = _pdf.generate(_spec)

        _hg = _data.createHistogram(vars[varNameX], vars[varNameY],
                                    varNBins[varNameX], varNBins[varNameY])
        _hg.Sumw2()

        # we want TH2D not TH2F
#        _h = ROOT.TH2D(name, name,
#                       _hg.GetXaxis().GetNbins(),
#                       _hg.GetXaxis().GetXmin(),
#                       _hg.GetXaxis().GetXmax(),
#                       _hg.GetYaxis().GetNbins(),
#                       _hg.GetYaxis().GetXmin(),
#                       _hg.GetYaxis().GetXmax() )
#        _h.Sumw2()
#        _h.Add(_hg)
#        _h.SetName(name)
#        _h.SetTitle(name)
            
        return _hg


    def GetToyData(self, muSig, muBg, mode=None):
        #
        # Generate and return a toy dataset
        # according to sig and bg yield Poisson means
        #
        legend = '[GetToyData]:'

        data = None

        nsig = self.random.Poisson(muSig)
        nbg = self.random.Poisson(muBg)

        print legend, 'nsig and nbg:', nsig, nbg

        # for testing, just return observed data
        #if self.workspace == None:
        #    return None

        if mode == 'uncor_2d_from_pdfs':
            print legend, 'generating 2D toys from 1d PDFs'

            toyData = self.Generate2D('mass', 'sig_mass_pdf',
                                      'ht', 'sig_ht_sum_pdf', nsig)

            toyData.Add(self.Generate2D('mass', 'bg_mass_pdf',
                                        'ht', 'bg_ht_sum_pdf', nbg))

#        else:
#            
#            nbinx = self.bgTemplate.GetNbinsX() 
#            xmin = self.bgTemplate.GetXaxis().GetBinLowEdge(1)
#            xmax = xmin + nbinx*self.bgTemplate.GetBinWidth(1)
#            nbiny = self.bgTemplate.GetNbinsY() 
#            ymin = self.bgTemplate.GetYaxis().GetBinLowEdge(1)
#            ymax = ymin + nbiny*self.bgTemplate.GetBinWidth(1)
#            _minv = ROOT.Double(0.0)
#            _ht = ROOT.Double(0.0)
#            toyData = ROOT.TH2F('toy', 'toy',
#                                nbinx, xmin, xmax,
#                                nbiny, ymin, ymax)
#            toyData.Sumw2()
#
#            # generate background toy entries
#            for i in range(0,int(nbg)):
#                self.bgTemplate.GetRandom2(_minv, _ht)
#                toyData.Fill(_minv, _ht)
#
#            # generate signal toy entries
#            for i in range(0,int(nsig)):
#                self.sigTemplate.GetRandom2(_minv, _ht)
#                toyData.Fill(_minv, _ht)

        #toyData.Draw('LEGO')
        #raw_input('hit enter...')

        return toyData



    def CombineBins(self, hist):
        #
        # Use a 2D->1D map and convert a 2D hist
        # into 1D
        #
        # Format of the map: list of bin lists that correspond to each combined bin
        #
        legend = '[Toys.CombineBins]:'
        
        #hist.Print()

        hist1d = ROOT.TH1F('combHist', 'combHist',
                           len(self.binMap),0, len(self.binMap))

        for i in self.binMap:

            for bin in self.binMap[i]:

                hist1d.AddBinContent(i,hist.GetBinContent(bin))

        #print legend, hist1d.Integral(), 'event in the combined histogram'
        #hist1d.Draw()
        #raw_input('press <enter> to continue...')

        #self.workspace.Print()
        
        data = ROOT.RooDataHist('toy', 'toy',
                                ROOT.RooArgList(self.workspace.var(obsName)),
                                hist1d)
        return data
    
        
        
    def GetTestResult(self, data, pdf, mode='fit'):
        #
        # Run a test on generated toy data
        # Let's start with getting a pull...
        #
        legend = '[GetTestResult]:'

        poi = {}

        # get workspace
        ws = self.workspace

        if mode=='fit':
            # set some constants
            #ws.var('Lumi').setConstant(ROOT.kTRUE)
            #ws.var('alpha_btag_syst').setConstant(ROOT.kTRUE)
            #ws.var('alpha_lepton_eff_e').setConstant(ROOT.kTRUE)
            #ws.var('alpha_norm_ewk_syst').setConstant(ROOT.kTRUE)
            #ws.var('alpha_norm_top_syst').setConstant(ROOT.kTRUE)
        
            #print legend
            #data.Print()
            #pdf.Print()

            # Profile the model
            #nll = pdf.createNLL(data)
            #profile = nll.createProfile(ROOT.RooArgSet())
            # this will do fit and set nuisance parameters to profiled values
            #profile.getVal()

            # or just fit
            pdf.fitTo(data, ROOT.RooFit.Save())

            poi['value'] = self.workspace.var(poiName).getVal()
            poi['error'] = self.workspace.var('xsec').getError()

            #print legend, 'poi =', poi

        elif mode=='run_inverter':
            print legend, 'running inverter'
            #invCalc.RunInverter(ws, data,
            #                    'ModelConfig', '',
            #                    2,3,
            #                    True, 10, 0, 0.5,
            #                    1,False,'')
            ROOT.RooStats.StandardHypoTestInvDemo(workspaceFileName,
                                                  workspaceName,
                                                  'ModelConfig', '',
                                                  'obsData',
                                                  2, 3,
                                                  True, 10, 0, 0.5,
                                                  1, False, '')
            poi['value'] = 1.0
            poi['error'] = 1.0

        # debug: plot fit result
        #frame = self.workspace.var(obsName).frame()
        #self.workspace.Print()
        #self.workspace.pdf('ejets_model').Print()
        #pdf.Print()
        #data.plotOn(frame)
        #pdf.plotOn(frame)
        #pdf.plotOn(frame, ROOT.RooFit.Components(self.workspace.pdf('binWidth_obs_x_ejets_1')), ROOT.RooFit.LineColor(ROOT.kRed))
        #frame.Draw()
        #raw_input('press <enter> to continue...')

        return poi



    def GetTheoryXsec(self, mass):
        #
        # Return t' theory cross section based on a few mass points
        # Interpolates linearly between the points
        #
        legend = '[GetTheoryXsec]:'

        xsec = None

        # list of tuples with point coordinates
        points=[(350.0,3.19947),
                (375,0,2.0996),
                (400,0,1.4055),
                (425,0,957418),
                (450.0,0.66227),
                (475,0,0.464403),
                (500,0,0.329632),
                (525,0,0.236506),
                (550,0,0.171393),
                (575,0,0.125289),
                (600,0,0.0922756),
                (625,0,0.0684641),
                (650,0,0.0511335)
                ]
        
        npoints = len(points)
        #_keys = sorted(points.iterkeys())

        x1 = None
        x2 = None
        y1 = None
        y2 = None

        print legend, points

        if mass<points[0][0]:
            x1 = points[0][0]
            x2 = points[1][0]
            y1 = points[0][1]
            y2 = points[1][1]
           
        elif mass>points[npoints-1][0]:
            x1 = points[npoints-2][0]
            x2 = points[npoints-1][0]
            y1 = points[npoints-2][1]
            y2 = points[npoints-1][1]

        else:
            for i in range(0,len(points)):
                if mass < points[i][0]:
                    x1 = points[i-1][0]
                    x2 = points[i][0]
                    y1 = points[i-1][1]
                    y2 = points[i][1]
                    break

        xsec = (y1-y2)/(x1-x2) * mass + (y2*x1-y1*x2)/(x1-x2)

        print legend, x1,y1,x2,y2
        
        return xsec



    def getTemplate(self, type, mass=None):

        hist = None

        #print self.data.hists
        if type == 'signal':
            #name = 'TPrime'+str(int(mass))+'_HtvsMfit'
            #name = 'TPrime'+str(int(mass))+'_ht35:fitMass'
            name = 'tprime'+str(int(mass))+'_ht35:fitMass'
            hist = self.data.hists[name]

        if type == 'background':
            #name = 'TTjets_HtvsMfit'
            #name = 'Top_ht35:fitMass'
            
            names = ['Top_ht35:fitMass', 'Ewk_ht35:fitMass']
            firstEntry = True
            for name in names:
                if firstEntry:
                    hist = copy.deepcopy(self.data.hists[name])
                    firstEntry = False
                else:
                    hist.Add(copy.deepcopy(self.data.hists[name]))
            
        if type == 'top':
            
            names = ['Top_ht35:fitMass']
            firstEntry = True
            for name in names:
                if firstEntry:
                    hist = copy.deepcopy(self.data.hists[name])
                    firstEntry = False
                else:
                    hist.Add(copy.deepcopy(self.data.hists[name]))
            
        if type == 'ewk':
            
            names = ['Ewk_ht35:fitMass']
            firstEntry = True
            for name in names:
                if firstEntry:
                    hist = copy.deepcopy(self.data.hists[name])
                    firstEntry = False
                else:
                    hist.Add(copy.deepcopy(self.data.hists[name]))
            
        return hist

    
    
    def SaveAllHists(self, filename, hist_list):
        outFile = ROOT.TFile(filename, "recreate")
        outFile.cd()
        print hist_list
        for hist in hist_list:
            hist.Write()
        outFile.Close()



    def Limit2dToys(self):
        #
        # Compute limits for provided 2D toy MC data histos

        legend = '[Toys.Generate2dToys]:'

        print legend, 'starting...'

#        ROOT.RooStats.StandardHypoTestInvDemo(workspaceFileName,
#                                              'ejets',
#                                              'ModelConfig', '',
#                                              'obsData',
#                                              2, 3, True, 11, 0, 0.5)
        invCalc.LoadWorkspace(workspaceFileName, 'ejets')
        invCalc.GetWorkspace().Print()



    def Generate2dToys(self):
        #
        # Generate 2D toy MC and save to file

        legend = '[Toys.Generate2dToys]:'
        
        # user set parameters
        #   assign with global parameters
        nToys = nSbToys
        pdfName = sbPdfName
        mass = tprimeMass
        xsec = tprimeXsec

        # load pdfs (if needed)
        self.LoadPdfs(pdfFileName, pdfWsName, pdfNames, varNames)

        # theoretical cross section
        xsecExp = self.GetTheoryXsec(mass)
        print legend, xsecExp

        # load 2D hists, to be used as templates/pseudodata
        # should be statistically independent from the model
        self.data.load_all_hists(templates2d)

        self.sigTemplate = self.getTemplate('signal', mass)
        self.bgTemplate = self.getTemplate('background')

        # expected signal and bg yields
        nSigExp = xsec/xsecExp * self.sigTemplate.Integral()*0.97
        nBgExp = self.bgTemplate.Integral()*1.067
        print legend, 'expected signal and background yields:', nSigExp, nBgExp
        
        # run toy loop
        print legend, 'beginning toy loop'
        toy_data_2d = []
        for n_toy in  range(1,nToys+1):
            _hist = self.GetToyData(nSigExp, nBgExp, 'uncor_2d_from_pdfs')
            _hist.SetName('toy_'+string.zfill(n_toy,6))
            _hist.SetTitle('toy_'+string.zfill(n_toy,6))
            toy_data_2d.append(_hist)

        self.SaveAllHists(toyMcFile, toy_data_2d)



    def run(self):
        #
        # Generate toy MC in a loop, fit, save pulls and such
        #
        # Exercise 1:
        #   - two statistically independent sets of 2D signal and background templates
        #     - no data
        #     - merge 2D templates into signal, Top and Ewk categories
        #     - pseudodata from EWK and Top templates, with exact expected yieds
        #       (probably is not necessary, only as a place holder for the interface)
        #   - 2D->1D mapping of the -first- set of templates
        #     - save the bin map file
        #   - produce single channel (ejets) HistFactory model (workspace)
        #   - generate 2D toys from the -second- set of templates
        #     - make 1D using the bin map
        #     - fit the model to the toys
        #     - save POI pull
        #
        
        legend = '[Toys.run]:'

        # user set parameters
        #   assign with global parameters
        nToys = nSbToys
        pdfName = sbPdfName
        mass = tprimeMass
        xsec = tprimeXsec

        # fit or limit or...
        testMode = 'fit'
        if useInverter:
            testMode = 'run_inverter'

        # bin map - global parameter
        self.loadBinMap(binMap)

        # load pdfs (if needed)
        self.LoadPdfs(pdfFileName, pdfWsName, pdfNames, varNames)

        # theoretical cross section
        xsecExp = self.GetTheoryXsec(mass)
        print legend, xsecExp

        # load 2D hists, to be used as templates/pseudodata
        # should be statistically independent from the model
        self.data.load_all_hists(templates2d)

        self.sigTemplate = self.getTemplate('signal', mass)
        self.bgTemplate = self.getTemplate('background')

        # expected signal and bg yields
        nSigExp = xsec/xsecExp * self.sigTemplate.Integral()*0.97
        nBgExp = self.bgTemplate.Integral()*1.067
        print legend, 'expected signal and background yields:', nSigExp, nBgExp
        
        # load the HistFactory model, which we want to study
        # should be statistically independent from the toy templates/pseudodata
        self.loadWorkspace(workspaceName, workspaceFileName)

        # get model pdf
        pdf = self.workspace.pdf(pdfName)
        pdf.Print()
        
        # create some diagnostics histograms
        hPull = ROOT.TH1F('pull', 'pull', 20, -10, 10)

        # run toy loop
        print legend, 'beginning toy loop'
        for n_toy in  range(1,nToys+1):

            toy_data_2d = self.GetToyData(nSigExp, nBgExp, 'uncor_2d_from_pdfs')
            toy_data_2d.Print()
            toy_data = self.CombineBins(toy_data_2d)

            poi = self.GetTestResult(toy_data, pdf, testMode)
            pull = (poi['value']-xsec)/poi['error']

            hPull.Fill(pull)

            print legend, 'poi =', poi['value'], '+/-', poi['error']
            #print legend, 'pull =', pull

        # draw hists
        hPull.Draw()
        raw_input('press <enter> to continue...')
        
        # draw the 2d background
        #self.drawHist('TTjets_HtvsMfit','LEGO')
    
        print self.legend, 'done.'



        #def SavePdf(self, ws, pdfName=None, dataName=None, doFit=True):
        #
        # Save pdf and data to 
        #



    def FitPdfs(self, mass):
        #
        # Get analytical PDF that describes tprime mass, to some degree
        #
        
        # load 2D templates
        self.data.load_all_hists(templates2d)
        hSig = self.getTemplate('signal', mass)
        hBg  = self.getTemplate('background')
        #hBg  = self.getTemplate('top')
        #hBg.Draw('LEGO')
        #hBg.Print()
        #raw_input('press <enter> to continue...')

        # project signal template onto mass
        hSigMass = hSig.ProjectionX("sig_mass", 0, -1, "e")
        hBgMass  = hBg.ProjectionX("bg_mass", 0, -1, "e")
        hSigHt   = hSig.ProjectionY("sig_ht", 0, -1, "e")
        hBgHt   = hBg.ProjectionY("bg_ht", 0, -1, "e")
        #hBgMass.Draw()
        #hBg.Print()
        #hSigHt.Draw()
        #hBgHt.Draw()
        #hSigHt.Print()
        #raw_input('press <enter> to continue...')

        # load the workspace (HistFactory)
        # self.loadWorkspace(workspaceName, workspaceFileName)
        # ws = self.workspace

        # open file for workspace and plots
        outFile = ROOT.TFile(outFileName, "RECREATE")

        # create workspace
        ws = ROOT.RooWorkspace('toyWs')
        

        # generic signal mass PDF
        massSigPdf = ws.factory('CBShape::sig_mass_pdf(mass[450,100,1000], sig_peak[450], sig_sigma[50], sig_alpha[1.0], sig_n[5])')

        # generic background mass PDF
        massBgPdf = ws.factory('Landau::bg_mass_pdf(mass, bg_mean[200], bg_sigma[40])')

        # generic signal HT PDF
        htSigGauss1Pdf = ws.factory('Gaussian::sig_ht_gauss1_pdf(ht[1000,200,2000], sig_ht_gmean1[850], sig_ht_gsig1[140])')
        htSigGauss2Pdf = ws.factory('Gaussian::sig_ht_gauss2_pdf(ht[1000,200,2000], sig_ht_gmean2[1200], sig_ht_gsig2[275])')
        htSigPdf       = ws.factory('SUM::sig_ht_sum_pdf(ht_sig_a1[0.75]*sig_ht_gauss1_pdf,sig_ht_gauss2_pdf)')

        # generic background HT PDF
        htBgGauss1Pdf = ws.factory('Gaussian::bg_ht_gauss1_pdf(ht, bg_ht_gmean1[540], bg_ht_gsig1[60])')
        htBgLandauPdf = ws.factory('Landau::bg_ht_gauss2_pdf(ht, bg_ht_gmean2[700], bg_ht_gsig2[50])')
        htBgPdf       = ws.factory('SUM::bg_ht_sum_pdf(ht_bg_a1[0.5]*bg_ht_gauss1_pdf,bg_ht_gauss2_pdf)')

        # tprime mass PDF
        massPdfCore = ws.factory('CBShape::tprime_mass_core_pdf(mass, peak[450,420,480], sigma[50,20,60], alpha[0.5,0,1], n[1,0,5])')
        massGauss1  = ws.factory('Gaussian::tprime_mass_gauss1_pdf(mass, gmean1[420,100,1000], gsig1[160,0,1000])')
        massPdfSum  = ws.factory('SUM::tprime_mass_sum_pdf(a1[1,0,1]*tprime_mass_gauss1_pdf,tprime_mass_core_pdf)')
        massPdf     = ws.factory('CEXPR::tprime_mass_pdf("tprime_mass_sum_pdf*(1.0-exp((100.0-mass)/100.0))",tprime_mass_sum_pdf,mass)')

        # background mass PDF
        massBgGauss   = ws.factory('Gaussian::allbg_mass_gauss_pdf(mass, gmean[175,150,200], gsig[10,5,50])')
        massBgLogn   = ws.factory('Lognormal::allbg_mass_logn_pdf(mass, lmean[175,150,300], kappa[3,0.5,10])')
        massBgPdfSum = ws.factory('SUM::allbg_mass_sum_pdf(a2[0.5,0.0,1.0]*allbg_mass_logn_pdf,allbg_mass_gauss_pdf)')
        massBgPdf    = ws.factory('CEXPR::allbg_mass_pdf("allbg_mass_sum_pdf*(1.0-exp((100.0-mass)/100.0))",allbg_mass_sum_pdf,mass)')

        # tprime HT PDF
        htTprimeGauss1Pdf = ws.factory('Gaussian::tprime_ht_gauss1_pdf(ht, tprime_ht_gmean1[800,200,2000], tprime_ht_gsig1[200,0,1000])')
        htTprimeGauss2Pdf = ws.factory('Gaussian::tprime_ht_gauss2_pdf(ht, tprime_ht_gmean2[1200,200,2000], tprime_ht_gsig2[400,0,1000])')
        htTprimePdf = ws.factory('SUM::tprime_ht_sum_pdf(ht_tprime_a1[1,0,1]*tprime_ht_gauss1_pdf,tprime_ht_gauss2_pdf)')
        
        # allbg HT PDF
        htAllbgGauss1Pdf = ws.factory('Gaussian::allbg_ht_gauss1_pdf(ht, allbg_ht_gmean1[600,200,2000], allbg_ht_gsig1[200,0,1000])')
        htAllbgLandauPdf = ws.factory('Landau::allbg_ht_gauss2_pdf(ht, allbg_ht_gmean2[1200,200,2000], allbg_ht_gsig2[400,0,1000])')
        htAllbgPdf = ws.factory('SUM::allbg_ht_sum_pdf(ht_allbg_a1[1,0,1]*allbg_ht_gauss1_pdf,allbg_ht_gauss2_pdf)')
        
        # make a dataset
        # Create a binned dataset that imports contents of TH1 and associates its contents to observable 'x'
        dhTprimeMass = ROOT.RooDataHist('dhTprimeMass', 'dhTprimeMass', ROOT.RooArgList(ws.var('mass')), hSigMass, 1) ;
        dhAllbgMass = ROOT.RooDataHist('dhAllbgMass', 'dhAllbgMass', ROOT.RooArgList(ws.var('mass')), hBgMass, 1) ;
        dhTprimeHt = ROOT.RooDataHist('dhTprimeHt', 'dhTprimeHt', ROOT.RooArgList(ws.var('ht')), hSigHt, 1) ;
        dhAllbgHt = ROOT.RooDataHist('dhAllbgHt', 'dhAllbgHt', ROOT.RooArgList(ws.var('ht')), hBgHt, 1) ;

        # fit mass shape to mass hist
        massPdf.fitTo(dhTprimeMass)
        massBgPdf.fitTo(dhAllbgMass)
        htTprimePdf.fitTo(dhTprimeHt)
        htAllbgPdf.fitTo(dhAllbgHt)

        # write workspace to file
        recreate = False
        ws.Write()
        #ws.writeToFile(outFileName,recreate)


        # plot and save to file
        
        frame = ws.var('mass').frame(ROOT.RooFit.Name('sig_mass'))
        massSigPdf.plotOn(frame)
        #frame.Draw()
        frame.Write()
        #raw_input('hit enter...')

        frame = ws.var('mass').frame(ROOT.RooFit.Name('bg_mass'))
        massBgPdf.plotOn(frame)
        #frame.Draw()
        frame.Write()
        #raw_input('hit enter...')

        frame = ws.var('ht').frame(ROOT.RooFit.Name('sig_ht'))
        htSigPdf.plotOn(frame)
        #frame.Draw()
        frame.Write()
        #raw_input('hit enter...')

        frame = ws.var('ht').frame(ROOT.RooFit.Name('bg_ht'))
        htBgPdf.plotOn(frame)
        #frame.Draw()
        frame.Write()
        #raw_input('hit enter...')

        #c1 = ROOT.TCanvas()
        frame = ws.var('mass').frame(ROOT.RooFit.Name('tprime_mass'))
        dhTprimeMass.plotOn(frame)
        massPdf.plotOn(frame)
        #frame.Draw()
        frame.Write()
        #raw_input('press <enter> to continue...')

        frame = ws.var('mass').frame(ROOT.RooFit.Name('allbg_mass'))
        dhAllbgMass.plotOn(frame)
        massBgPdf.plotOn(frame)
        #frame.Draw()
        frame.Write()
        #raw_input('press <enter> to continue...')

        frame = ws.var('ht').frame(ROOT.RooFit.Name('tprime_ht'))
        dhTprimeHt.plotOn(frame)
        htTprimePdf.plotOn(frame)
        #frame.Draw()
        frame.Write()
        #raw_input('press <enter> to continue...')

        frame = ws.var('ht').frame(ROOT.RooFit.Name('allbg_ht'))
        dhAllbgHt.plotOn(frame)
        htAllbgPdf.plotOn(frame)
        #frame.Draw()
        frame.Write()
        #raw_input('press <enter> to continue...')

        self.massPdf = massPdf
        return self.massPdf



    def Generate1DTemplates(self):
        #
        # Load workspace, generate histograms,
        # Save to file
        #

        inFile = ROOT.TFile(inPdfFileName, 'read')
        ws = inFile.Get(workspaceName)

        pdfs = []
        entry = 0
        for name in pdfNames:
            pdfs.append(ws.pdf(name).cloneTree())

        inFile.Close()

        outFile = ROOT.TFile(outTemplateFile, 'recreate')

        entry = 0
        for name in histNames:
            _h = ROOT.TH1D(name, name, histNBins[entry], histLow[entry], histHigh[entry])
            _spec = pdfs[entry].prepareMultiGen(ROOT.RooArgSet(ws.var(vars[entry])),
                                                ROOT.RooFit.Name('ds_'+name),
                                                ROOT.RooFit.NumEvents(nEvents[entry]),
                                                ROOT.RooFit.Extended(ROOT.kFALSE),
                                                ROOT.RooFit.Verbose(ROOT.kFALSE))
            _data = pdfs[entry].generate(_spec)
            _data.fillHistogram(_h, ROOT.RooArgList(ws.var(vars[entry])))

            _h.Sumw2()
            if yields[entry]>0:
                _h.Scale(yields[entry]/_h.Integral())

            #ws.var(vars[entry]) )
            #_h.Draw()
            _h.Write()
            #raw_input('hit enter...')
            entry += 1
            
        outFile.Close()
        


    def Generate2DTemplates(self):
        #
        # Load workspace, generate histograms,
        # Save to file
        #

        inFile = ROOT.TFile(inPdfFileName, 'read')
        ws = inFile.Get(workspaceName)

        pdfs  = []
        pdfsY = []
        entry = 0
        for name in pdfNames:
            pdfs.append(ws.pdf(name).cloneTree())

        for name in pdfNamesY:
            pdfsY.append(ws.pdf(name).cloneTree())

        inFile.Close()

        outFile = ROOT.TFile(outTemplateFile, 'recreate')

        entry = 0
        for name in histNames:
            _pdf = ROOT.RooProdPdf('tmpPdf', 'tmpPdf', pdfs[entry], pdfsY[entry])
            _spec = _pdf.prepareMultiGen(ROOT.RooArgSet(ws.var(vars[entry]),ws.var(varsY[entry])),
                                         ROOT.RooFit.Name('ds_'+name),
                                         ROOT.RooFit.NumEvents(nEvents[entry]),
                                         ROOT.RooFit.Extended(ROOT.kFALSE),
                                         ROOT.RooFit.Verbose(ROOT.kFALSE))
            _data = _pdf.generate(_spec)
            _hg = _data.createHistogram(ws.var(vars[entry]), ws.var(varsY[entry]),
                                     histNBins[entry], histNBinsY[entry])

            # we want TH2D not TH2F
            _h = ROOT.TH2D(name, name,
                           _hg.GetXaxis().GetNbins(),
                           _hg.GetXaxis().GetXmin(),
                           _hg.GetXaxis().GetXmax(),
                           _hg.GetYaxis().GetNbins(),
                           _hg.GetYaxis().GetXmin(),
                           _hg.GetYaxis().GetXmax() )

            _h.Sumw2()
            _hg.Sumw2()

            _h.Add(_hg)

            _h.SetName(name)
            _h.SetTitle(name)
            
            if yields[entry]>0:
                _h.Scale(yields[entry]/_h.Integral())

            #ws.var(vars[entry]) )
            #_h.Draw()
            _h.Write()
            #raw_input('hit enter...')
            entry += 1
            
        outFile.Close()
        


    def Get2dHist(self, proto, ngen = -1, poisson = False):
        #
        # Generate a histogram according to the proto hist
        #  - proto hist is used as a PDF
        #  - number of entries is taken from proto
        #    - positive ngen overrides number of entries
        #  - normalization taken from proto
        #  - poisson controls whether number of entries is randomized
        #
        
        legend = '[Get2dHist]:'

        hist = None

        if ngen < 0:
            mu = proto.GetEntries()
        else:
            mu = ngen

        if poisson:
            n_entries = self.random.Poisson(mu)
        else:
            n_entries = mu

        print legend, 'entries to generate:', n_entries

        nbinx = proto.GetNbinsX() 
        xmin = proto.GetXaxis().GetBinLowEdge(1)
        xmax = xmin + nbinx*proto.GetBinWidth(1)
        nbiny = proto.GetNbinsY() 
        ymin = proto.GetYaxis().GetBinLowEdge(1)
        ymax = ymin + nbiny*proto.GetBinWidth(1)
        _minv = ROOT.Double(0.0)
        _ht = ROOT.Double(0.0)
        toyData = ROOT.TH2F(proto.GetName(), proto.GetTitle(),
                            nbinx, xmin, xmax,
                            nbiny, ymin, ymax)
        toyData.Sumw2()

        # generate toy entries
        for i in range( 0, int(n_entries) ):
            proto.GetRandom2(_minv, _ht)
            toyData.Fill(_minv, _ht)

        #toyData.Draw('LEGO')
        #raw_input('hit enter...')

        return toyData



    def GetHist(self, proto, ngen = -1, poisson = False, do_not_scale = False):
        #
        # Generate a histogram according to the proto hist
        #  - proto hist is used as a PDF
        #  - number of entries is taken from proto
        #    - positive ngen overrides number of entries
        #  - normalization taken from proto
        #  - poisson controls whether number of entries is randomized
        #
        # Supposedly, clever enough to automatically tell
        # between 1d and 2d proto hists and generate accordingly
        #
        
        legend = '[GetHist]:'

        hist = None

        if ngen < 0:
            mu = proto.GetEntries()
        else:
            mu = ngen

        if poisson:
            n_entries = self.random.Poisson(mu)
        else:
            n_entries = mu

        print legend, 'entries to generate:', n_entries

        toyData = None

        is_2d = None

        if 'TH2' in proto.ClassName():
            is_2d = True
        elif 'TH1' in proto.ClassName():
            is_2d = False
        else:
            print legend, 'prototype hist is neither TH1 nor TH2...'
            return toyData

        nbinx = proto.GetNbinsX() 
        xmin = proto.GetXaxis().GetBinLowEdge(1)
        xmax = xmin + nbinx*proto.GetBinWidth(1)
        nbiny = None

        if is_2d:
            nbiny = proto.GetNbinsY() 
            ymin = proto.GetYaxis().GetBinLowEdge(1)
            ymax = ymin + nbiny*proto.GetBinWidth(1)
            
            _minv = ROOT.Double(0.0)
            _ht = ROOT.Double(0.0)
            toyData = ROOT.TH2F(proto.GetName(), proto.GetTitle(),
                                nbinx, xmin, xmax,
                                nbiny, ymin, ymax)
        else:
            toyData = ROOT.TH1F(proto.GetName(), proto.GetTitle(),
                                nbinx, xmin, xmax)
            
        toyData.Sumw2()

        # generate toy entries
        if is_2d:
            for i in range( 0, int(n_entries) ):
                proto.GetRandom2(_minv, _ht)
                toyData.Fill(_minv, _ht)
            # get integrals with overflow
            _int_proto = proto.Integral(0,nbinx+1, 0,nbiny+1)
            _int_toydata = toyData.Integral(0,nbinx+1, 0,nbiny+1)
        else:
            toyData.FillRandom(proto, int(n_entries))
            # get integrals with overflow
            _int_proto = proto.Integral(0,nbinx+1)
            _int_toydata = toyData.Integral(0,nbinx+1)

        #toyData.Draw('LEGO')
        #raw_input('hit enter...')

        # normalize the histogram
        if not do_not_scale:
            toyData.Scale(_int_proto/_int_toydata)

        return toyData



    def GenerateToyHists(self,
                         hdata_templates,    # templates for toy templates
                         mass,
                         hdata_data = None): # templates for toy data, if given
        #
        # generate toy data and toy templates
        # input:  TprimeData object (or 2 for separate data master template)
        # output: TprimeData object
        #

        legend = 'Toys.GenerateToyHists:'

        hdata = hdata_templates
    
        rdata = TprimeData(hdata.GetLumi())
        
        for name in hdata.hist_names:
            if mass and 'tprime' in name and str(mass) not in name:
                continue
            
            #if name != 'DATA':
            if 'DATA' not in name:
                rdata.AddHist(self.GetHist(hdata.hists[name]))

        # if separate templates for toy data provided, use them
        if hdata_data:
            hdata = hdata_data

        # generate toy DATA
        print legend, 'using', hdata.object_name, 'to generate data'
        
        _ntop = None
        _newk = None
        _binx = hdata.hists['top'].GetNbinsX()
        _biny = None
        if 'TH2' in hdata.hists['top'].ClassName():
            _biny = hdata.hists['top'].GetNbinsY()
            _ntop = hdata.hists['top'].Integral(0,_binx+1, 0,_biny+1)
            _newk = hdata.hists['ewk'].Integral(0,_binx+1, 0,_biny+1)
        elif 'TH1' in hdata.hists['top'].ClassName():
            _ntop = hdata.hists['top'].Integral(0,_binx+1)
            _newk = hdata.hists['ewk'].Integral(0,_binx+1)
        else:
            print legend, 'unknown histogram class'
            return rdata

        _data = self.GetHist(hdata.hists['top'],
                             ngen=_ntop,
                             poisson=True, do_not_scale = True)
        _data.SetName('DATA')
        _data.SetTitle('pseudo-data')
        rdata.AddHist(_data)
        _data = self.GetHist(hdata.hists['ewk'],
                             ngen=_newk,
                             poisson=True, do_not_scale = True)
        _data.SetName('DATA')
        rdata.AddHist(_data)
            
        return rdata


    
def RunToy(par):
    #
    # Run the full toy chain to check the model:
    #   - master templates
    #   - toy templates
    #   - toy data
    #   - merge toy templates
    #   - rebin toy templates and data
    #   - make histfactory configs
    #   - make histfactory model
    #   - get null hypo p-value
    #
    # Usage:
    #        ./tprime.py --toys-pvalue
    #        ./tprime.py --toys limit
    #        ./tprime.py --toys pvalue
    #        ./tprime.py --toys templates
    #
    # Parameters:
    #             par['study'] = 'pvalue'(default), 'limit'
    #             mass - tprime mass. If None, all masses are processed.
    #
    #
    #  Possible full sequence:
    #     master -> merge -> rebin -> toy -> toy_rebin -> histfactory -> limit(or pvalue)
    #         _2 ->    _2 ->    _2 ->
    #

    legend = '[tprime_toys.RunToy]:'

    import os

    # safe defaults
    input_mu_jer_plus  = None
    input_mu_jer_minus = None
    input_mu_btgsf_plus  = None
    input_mu_btgsf_minus = None
    input_e_btgsf_plus   = None
    input_e_btgsf_minus  = None
    input_e_jer_plus   = None
    input_e_jer_minus  = None
    input_e_testsyst_plus   = None
    input_e_testsyst_minus  = None
    input_mu_testsyst_plus   = None
    input_mu_testsyst_minus  = None

    #----------------------------------------------------------------------->
    #
    # parameters of the study


    # type of the toy study
    _study = 'pvalue'            #default
    if 'study' in par.keys():
        _study = par['study']

        out_fname = 'output.txt'
    if _study=='limit':
        out_fname = 'limits.txt'
    if _study=='pvalue':
        out_fname = 'pvalues.txt'
    if _study=='fit':
        out_fname = 'fit_results.txt'

    # for creating templates only: prefix for the output hist names
    #e_prefix  = ''
    #mu_prefix = ''
    e_prefix  = 'ele_tprime__'
    mu_prefix = 'mu_tprime__'


    # master templates (None if some are missing)
    # if jes plus/minus specified, will be merged
    #

    input_mu_nominal     = 'data/mujets_4600/from_gueorgi/19mar2012v1/mujets_4601ipb_2D_v1_csv_jer_nom.root'
    input_mu_jes_plus    = 'data/mujets_4600/from_gueorgi/19mar2012v1/mujets_4601ipb_2D_v1_csv_jes_plus.root'
    input_mu_jes_minus   = 'data/mujets_4600/from_gueorgi/19mar2012v1/mujets_4601ipb_2D_v1_csv_jes_minus.root'
    input_mu_match_plus  = 'data/mujets_4600/from_gueorgi/19mar2012v1/mujets_4601ipb_2D_v1_csv_matchingup.root'
    input_mu_match_minus = 'data/mujets_4600/from_gueorgi/19mar2012v1/mujets_4601ipb_2D_v1_csv_matchingdown.root'
    
    input_e_nominal     = 'data/ejets_4700/ricardo_19mar2012v1/ht4jetsAfterFit_vs_fitMass.root'
    input_e_jes_plus    = 'data/ejets_4700/ricardo_19mar2012v1/ht4jetsAfterFit_vs_fitMass_JES105.root'
    input_e_jes_minus   = 'data/ejets_4700/ricardo_19mar2012v1/ht4jetsAfterFit_vs_fitMass_JES095.root'
    input_e_match_plus  = 'data/ejets_4700/ricardo_19mar2012v1/ht4jetsAfterFit_vs_fitMass_matchup.root'
    input_e_match_minus = 'data/ejets_4700/ricardo_19mar2012v1/ht4jetsAfterFit_vs_fitMass_matchdown.root'

    # inputs used for Moriond-2012
    #input_mu_nominal   = 'data/mujets_4600/from_gueorgi/27jan2011v1/csv/v1/mujets_4601ipb_2D_v1_csv_jer_nom.root'
    #input_mu_jes_plus  = 'data/mujets_4600/from_gueorgi/27jan2011v1/csv/v1/mujets_4601ipb_2D_v1_csv_jes_plus.root'
    #input_mu_jes_minus = 'data/mujets_4600/from_gueorgi/27jan2011v1/csv/v1/mujets_4601ipb_2D_v1_csv_jes_minus.root'
    #input_mu_match_plus  = 'data/mujets_4600/from_gueorgi/27jan2011v1/csv/v1/mujets_4601ipb_2D_v1_csv_matchingup.root'
    #input_mu_match_minus = 'data/mujets_4600/from_gueorgi/27jan2011v1/csv/v1/mujets_4601ipb_2D_v1_csv_matchingdown.root'
     
    #input_e_nominal    = 'data/ejets_4700/ricardo_4683ipb_2d_csv_22jan2012v2/ht4jetsAfterFit_vs_fitMass.root'
    #input_e_jes_plus   = 'data/ejets_4700/ricardo_4683ipb_2d_csv_22jan2012v2/ht4jetsAfterFit_vs_fitMass_JES105.root'
    #input_e_jes_minus  = 'data/ejets_4700/ricardo_4683ipb_2d_csv_22jan2012v2/ht4jetsAfterFit_vs_fitMass_JES095.root'
    #input_e_match_plus  = 'data/ejets_4700/ricardo_4683ipb_2d_csv_22jan2012v2/ht4jetsAfterFit_vs_fitMass_matchup.root'
    #input_e_match_minus = 'data/ejets_4700/ricardo_4683ipb_2d_csv_22jan2012v2/ht4jetsAfterFit_vs_fitMass_matchdown.root'
    
    

    #
    # davis
    # 20%
    #data/ejets_4700/ricardo_4683ipb_2d_csv_22jan2012v2/adjacentBinsMergedHistograms/
    #input_mu_nominal     = 'data/mujets_4600/from_gueorgi/27jan2011v1/csv/v1/adjacentBinsMergedHistograms/mujets_4601ipb_2D_v1_csv_jer_nom_ganged_0.2.root'
    #input_mu_jes_plus    = 'data/mujets_4600/from_gueorgi/27jan2011v1/csv/v1/adjacentBinsMergedHistograms/mujets_4601ipb_2D_v1_csv_jes_plus_ganged_0.2.root'
    #input_mu_jes_minus   = 'data/mujets_4600/from_gueorgi/27jan2011v1/csv/v1/adjacentBinsMergedHistograms/mujets_4601ipb_2D_v1_csv_jes_minus_ganged_0.2.root'
    #input_mu_match_plus  = 'data/mujets_4600/from_gueorgi/27jan2011v1/csv/v1/adjacentBinsMergedHistograms/mujets_4601ipb_2D_v1_csv_btgsf_plus_ganged_0.2.root'
    #input_mu_match_minus = 'data/mujets_4600/from_gueorgi/27jan2011v1/csv/v1/adjacentBinsMergedHistograms/mujets_4601ipb_2D_v1_csv_btgsf_minus_ganged_0.2.root'

    #input_e_nominal     = 'data/ejets_4700/ricardo_4683ipb_2d_csv_22jan2012v2/adjacentBinsMergedHistograms/ht4jetsAfterFit_vs_fitMass_ganged_0.2.root'
    #input_e_jes_plus    = 'data/ejets_4700/ricardo_4683ipb_2d_csv_22jan2012v2/adjacentBinsMergedHistograms/ht4jetsAfterFit_vs_fitMass_JES105_ganged_0.2.root'
    #input_e_jes_minus   = 'data/ejets_4700/ricardo_4683ipb_2d_csv_22jan2012v2/adjacentBinsMergedHistograms/ht4jetsAfterFit_vs_fitMass_JES095_ganged_0.2.root'
    #input_e_match_plus  = 'data/ejets_4700/ricardo_4683ipb_2d_csv_22jan2012v2/adjacentBinsMergedHistograms/ht4jetsAfterFit_vs_fitMass_matchup_ganged_0.2.root'
    #input_e_match_minus = 'data/ejets_4700/ricardo_4683ipb_2d_csv_22jan2012v2/adjacentBinsMergedHistograms/ht4jetsAfterFit_vs_fitMass_matchdown_ganged_0.2.root'


    # shape systematic under study
    #input_e_testsyst_plus = input_e_powheg_plus
    #input_e_testsyst_minus = input_e_powheg_minus
    #input_mu_testsyst_plus = input_mu_scaleup
    #input_mu_testsyst_minus = input_mu_scaledown
    

    # optional for split MC study
    master_template_e_fname2 = None # safe default
    master_template_mu_fname2 = None # safe default


    #output_dir               = 'test_comb_08may2012v1'
    output_dir               = 'comb_diag_perp_28may2012v1'
    generate_toys            = False
    ntoys                    = 1
    channel                  = 'combined'  # ejets, mujets, combined
    masses                   = range(400,650,25)
    #masses                   = [550]
    lumi_e                   = 4900
    lumi_mu                  = 5000
    maxerr                   = 0.2
    maxerr_sig               = 0.2
    rebin_method             = 'signal_mass_diag_perp'    # brown, brown_smooth, brown_2d,
                                          # copy, davis, mass, brown_mike
                                          # signal_box_1bin, signal_box_2bin,signal_box_diagonal,signal_mass_diagonal, signal_mass_diag_perp
    poimin                   = 0.0
    poimax                   = 2.0
    cls_scan_points          = 150
    #out_fname                = 'limits.txt'
    rebin_input              = False # rebin master templates
                                     #once on the fly
                                     #only needed for split MC test
    make_rebin_fit_plots     = False
    
    toy_hists_ejets_fname    = 'test_toy_hists_ejets.root'
    toy_hists_mujets_fname   = 'test_toy_hists_mujets.root'
    hfModelFileNameBase      = 'test_model.root'
    hfInputEjetsFileName     = 'hf_input_ejets_test.root'
    hfInputMujetsFileName    = 'hf_input_mujets_test.root'


    # histfactory templates (old)
    #comb_template_name = 'hf_tprime_template_v1.xml'
    #e_template_name    = 'hf_tprime_e_template_v1.xml'
    #mu_template_name   = 'hf_tprime_mu_template_v1.xml'

    # histfactory templates (nominal) jes+match
    comb_template_name = 'hf_tprime_template_v3.xml'
    e_template_name    = 'hf_tprime_e_template_v3.xml'
    mu_template_name   = 'hf_tprime_mu_template_v3.xml'

    # histfactory templates (test)
    #comb_template_name = 'hf_tprime_template_v4.xml'
    #e_template_name    = 'hf_tprime_e_template_v4.xml'
    #mu_template_name   = 'hf_tprime_mu_template_v4.xml'

    # histfactory templates for toys
    comb_template_name_toy = 'hf_tprime_toy_template.xml'
    e_template_name_toy    = 'hf_tprime_toy_e_template.xml'
    mu_template_name_toy   = 'hf_tprime_toy_mu_template.xml'


    #----------------------------------------------------------------------->
    #
    # end of setting, beginning the macro
    #

    #-----> safety for prefix
    if _study=='templates':
        print legend, 'making templates, ejets histogram name prefixes are <'+e_prefix+'>'
        print legend, 'making templates, mujets histogram name prefixes are <'+mu_prefix+'>'
    else:
        e_prefix = ''
        mu_prefix = ''

    #-----> create output directory
    output_dir = output_dir.strip().rstrip('/')
    os.system('mkdir -p '+output_dir)

    #-----> add output dir name to output file names
    toy_hists_ejets_fname  = output_dir+'/'+toy_hists_ejets_fname
    toy_hists_mujets_fname = output_dir+'/'+toy_hists_mujets_fname
    #hfInputEjetsFileName   = output_dir+'/'+hfInputEjetsFileName
    #hfInputMujetsFileName  = output_dir+'/'+hfInputMujetsFileName
    out_fname              = output_dir+'/'+out_fname

    #-----> select histfactory template set
    if generate_toys:
        comb_template_name = comb_template_name_toy
        e_template_name    = e_template_name_toy
        mu_template_name   = mu_template_name_toy


    #-----> merge mu+jets inputs

    input_file_list = []
    suffix_list = []

    if input_mu_nominal:
        input_file_list.append(input_mu_nominal)
        suffix_list.append('')
    if input_mu_jes_plus:
        input_file_list.append(input_mu_jes_plus)
        suffix_list.append('__jes__plus')
    if input_mu_jes_minus:
        input_file_list.append(input_mu_jes_minus)
        suffix_list.append('__jes__minus')
    if input_mu_match_plus:
        input_file_list.append(input_mu_match_plus)
        suffix_list.append('__match__plus')
    if input_mu_match_minus:
        input_file_list.append(input_mu_match_minus)
        suffix_list.append('__match__minus')
    if input_mu_testsyst_plus:
        input_file_list.append(input_mu_testsyst_plus)
        suffix_list.append('__testsyst__plus')
    if input_mu_testsyst_minus:
        input_file_list.append(input_mu_testsyst_minus)
        suffix_list.append('__testsyst__minus')
    if input_mu_btgsf_plus:
        input_file_list.append(input_mu_btgsf_plus)
        suffix_list.append('__btgsf__plus')
    if input_mu_btgsf_minus:
        input_file_list.append(input_mu_btgsf_minus)
        suffix_list.append('__btgsf__minus')
    if input_mu_jer_plus:
        input_file_list.append(input_mu_jer_plus)
        suffix_list.append('__jer__plus')
    if input_mu_jer_minus:
        input_file_list.append(input_mu_jer_minus)
        suffix_list.append('__jer__minus')

    in_ds_mu = TprimeData(lumi_mu)

    i_suff = 0
    for name in input_file_list:
        _suffix = suffix_list[i_suff]
        _excluded_list = []
        if len(_suffix) > 0:
            _excluded_list.append('DATA__jes__plus')
            _excluded_list.append('DATA__jes__minus')
            _excluded_list.append('DATA__match__plus')
            _excluded_list.append('DATA__match__minus')
            _excluded_list.append('DATA__testsyst__plus')
            _excluded_list.append('DATA__testsyst__minus')
            _excluded_list.append('DATA__btgsf__plus')
            _excluded_list.append('DATA__btgsf__minus')
            _excluded_list.append('DATA__jer__plus')
            _excluded_list.append('DATA__jer__minus')
        in_ds_mu.load_all_hists(name,
                                merge_ewk=True,
                                merge_top=True,
                                strip_suffix = True,
                                suffix = _suffix,
                                excluded_list = _excluded_list)
        i_suff += 1
        
    master_template_mu_fname = output_dir+'/input_merged_mu.root'
    in_ds_mu.SaveAllHists(master_template_mu_fname)


    #-----> merge e+jets inputs

    input_file_list = []
    suffix_list = []

    if input_e_nominal:
        input_file_list.append(input_e_nominal)
        suffix_list.append('')
    if input_e_jes_plus:
        input_file_list.append(input_e_jes_plus)
        suffix_list.append('__jes__plus')
    if input_e_jes_minus:
        input_file_list.append(input_e_jes_minus)
        suffix_list.append('__jes__minus')
    if input_e_match_plus:
        input_file_list.append(input_e_match_plus)
        suffix_list.append('__match__plus')
    if input_e_match_minus:
        input_file_list.append(input_e_match_minus)
        suffix_list.append('__match__minus')
    if input_e_testsyst_plus:
        input_file_list.append(input_e_testsyst_plus)
        suffix_list.append('__testsyst__plus')
    if input_e_testsyst_minus:
        input_file_list.append(input_e_testsyst_minus)
        suffix_list.append('__testsyst__minus')
    if input_e_btgsf_plus:
        input_file_list.append(input_e_btgsf_plus)
        suffix_list.append('__btgsf__plus')
    if input_e_btgsf_minus:
        input_file_list.append(input_e_btgsf_minus)
        suffix_list.append('__btgsf__minus')
    if input_e_jer_plus:
        input_file_list.append(input_e_jer_plus)
        suffix_list.append('__jer__plus')
    if input_e_jer_minus:
        input_file_list.append(input_e_jer_minus)
        suffix_list.append('__jer__minus')

    in_ds_e = TprimeData(lumi_e)

    i_suff = 0
    for name in input_file_list:
        _suffix = suffix_list[i_suff]
        _excluded_list = []
        if len(_suffix) > 0:
            _excluded_list.append('DATA__jes__plus')
            _excluded_list.append('DATA__jes__minus')
            _excluded_list.append('DATA__match__plus')
            _excluded_list.append('DATA__match__minus')
            _excluded_list.append('DATA__testsyst__plus')
            _excluded_list.append('DATA__testsyst__minus')
            _excluded_list.append('DATA__btgsf__plus')
            _excluded_list.append('DATA__btgsf__minus')
            _excluded_list.append('DATA__jer__plus')
            _excluded_list.append('DATA__jer__minus')
        in_ds_e.load_all_hists(name,
                               merge_ewk=True,
                               merge_top=True,
                               strip_suffix = True,
                               suffix = _suffix,
                               excluded_list = _excluded_list)
        i_suff += 1
        
    master_template_e_fname = output_dir+'/input_merged_e.root'
    in_ds_e.SaveAllHists(master_template_e_fname)


    #-----> rebin e+jets master templates (merged input)
    
    if rebin_input:

        _newname = output_dir+'/input_merged_rebinned_e.root'
        cb.CombineBins(masses[0], scale=1.0, maxErr=maxerr,
                       infile=master_template_e_fname,
                       outfile=_newname,
                       prefix='',
                       method=rebin_method,
                       maxErrSig=maxerr_sig)
        _binmap_filename = _newname+'.binmap'
        master_template_e_fname  = _newname
        if master_template_e_fname2:
            _newname = output_dir+'/input_merged_rebinned_e_2.root'

            cb.CombineBinsWithMap(_binmap_filename,
                                  master_template_e_fname2,
                                  _newname)
            master_template_e_fname2 = _newname


    #-----> rebin mu+jets master templates (merged input)

        _newname = output_dir+'/input_merged_rebinned_mu.root'
        cb.CombineBins(masses[0], scale=1.0, maxErr=maxerr,
                       infile=master_template_mu_fname,
                       outfile=_newname,
                       prefix='',
                       method=rebin_method,
                       maxErrSig=maxerr_sig)
        _binmap_filename = _newname+'.binmap'
        master_template_mu_fname  = _newname
        if master_template_mu_fname2:
            _newname = output_dir+'/input_merged_rebinned_mu_2.root'

            cb.CombineBinsWithMap(_binmap_filename,
                                  master_template_mu_fname2,
                                  _newname)
            master_template_mu_fname2 = _newname

        rebin_method = 'copy'


    if generate_toys == False:
        ntoys = 1

    pvalue = None

    # create output dir
    #os.system('mkdir -p '+outdir.rstrip('/'))

    # load appropriate module for the requested study
    if _study == 'pvalue':
        ROOT.gROOT.ProcessLine('.L tprime_pvalue.C+')
        from ROOT import GetPValue
    elif _study == 'limit':
        ROOT.gROOT.ProcessLine('.L StandardHypoTestInvDemo.5.32.C+')
        from ROOT import StandardHypoTestInvDemo

    # load histfactory module
    import hist_factory as hf
    hf.comb_template_name = comb_template_name
    hf.e_template_name = e_template_name
    hf.mu_template_name = mu_template_name

    # load master templates
    # ejets
    ds_e = TprimeData(lumi_e)
    ds_e.load_all_hists(master_template_e_fname, strip_suffix = False)
    ds_e2 = None
    if master_template_e_fname2:
        ds_e2 = TprimeData(lumi_e, 'alternative_ds')
        ds_e2.load_all_hists(master_template_e_fname2, strip_suffix = False)
    e_toys = Toys(master_template_e_fname)
    # mujets
    ds_mu = TprimeData(lumi_mu)
    ds_mu.load_all_hists(master_template_mu_fname, strip_suffix = False)
    mu_toys = Toys(master_template_mu_fname)


    # toy loop
    ofile = open(out_fname, 'a')
    if _study == 'pvalue':
        ofile.write('# p-values\n')
        ofile.write('# ------------------\n')
    elif _study == 'limit':
        ofile.write('# limits\n')
        ofile.write('# \n')
        ofile.write('# mass observed -2sig -1sig median_exp +1sig +2sig\n')
        ofile.write('# -------------------------------------------\n')
    ofile.close()
    pvalue = None
    limit = None
    first_mass_point = True
    for mass in masses:

        # file name prefix for the templates action
        _e_fileprefix = ''
        _mu_fileprefix = ''
        if _study == 'templates':
            _e_fileprefix = e_prefix+str(mass)+'_'
            _mu_fileprefix = mu_prefix+str(mass)+'_'
        
        for itoy in range (0,ntoys):

            ofile = open(out_fname, 'a')

            # if toys
            if generate_toys:
                # generate and save a set of toy templates and data
                # ejets
                toy_hists_e = e_toys.GenerateToyHists(ds_e, mass, ds_e2)
                toy_hists_e.SaveAllHists(toy_hists_ejets_fname)
                # mujets
                toy_hists_mu = mu_toys.GenerateToyHists(ds_mu, mass)
                toy_hists_mu.SaveAllHists(toy_hists_mujets_fname)
            else:
                toy_hists_ejets_fname = master_template_e_fname
                toy_hists_mujets_fname = master_template_mu_fname
            
            # rebin
            # ejets
            cb.CombineBins(mass, scale=1.0, maxErr=maxerr,
                           infile=toy_hists_ejets_fname,
                           outfile=output_dir+'/'+_e_fileprefix+hfInputEjetsFileName,
                           prefix=e_prefix,
                           method=rebin_method,
                           maxErrSig=maxerr_sig,
                           make_plots = make_rebin_fit_plots,
                           channel = "ejets")
            # mujets
            cb.CombineBins(mass, scale=1.0, maxErr=maxerr,
                           infile=toy_hists_mujets_fname,
                           outfile=output_dir+'/'+_mu_fileprefix+hfInputMujetsFileName,
                           prefix=mu_prefix,
                           method=rebin_method,
                           maxErrSig=maxerr_sig,
                           make_plots = make_rebin_fit_plots,
                           channel="mujets")

            
            # make histfactory model
            if _study == 'pvalue' or _study == 'limit' or _study=='fit':
                hf.run(mass,
                       hfInputEjetsFileName,
                       hfInputMujetsFileName,
                       hfModelFileNameBase,
                       output_dir)
            
            # get null hypo p-value
            if _study == 'pvalue':
                pvalue = GetPValue(output_dir+'/'+channel+'_'+hfModelFileNameBase, channel)
                ofile.write(str(pvalue)+'  \n')
                # or get limits
            elif _study == 'limit':
                print legend, 'LIMIT'

                r = StandardHypoTestInvDemo(output_dir+'/'+channel+'_'+hfModelFileNameBase,
                                            channel,
                                            'ModelConfig',
                                            '',
                                            'obsData',
                                            2,
                                            3, 
                                            True,
                                            cls_scan_points,
                                            poimin,
                                            poimax)
                obs_limit = r.UpperLimit()
                exp_limit_min2  = r.GetExpectedUpperLimit(-2)
                exp_limit_min1  = r.GetExpectedUpperLimit(-1)
                exp_limit_med   = r.GetExpectedUpperLimit(0)
                exp_limit_plus1 = r.GetExpectedUpperLimit(1)
                exp_limit_plus2 = r.GetExpectedUpperLimit(2)
                #print legend, 'observed limit:', limit
                ofile.write(str(mass)+'  '+
                            str(obs_limit)+'  '+
                            str(exp_limit_min2)+'  '+
                            str(exp_limit_min1)+'  '+
                            str(exp_limit_med)+'  '+
                            str(exp_limit_plus1)+'  '+
                            str(exp_limit_plus2)+'  '+
                            '\n')

                ofile.close()

            elif _study == 'fit':
                print legend, 'doing a fit study'

                import tprime_fit_gena as tf
                print legend, output_dir+'/'+channel+'_'+hfModelFileNameBase
                print legend, channel
                print legend, mass
                
                fitRes = tf.FitSimple(output_dir+'/'+channel+'_'+hfModelFileNameBase,
                                      channel,
                                      mass)

                if first_mass_point:
                    ofile.write('# fit parameter values\n')
                    ofile.write('# mass   '+fitRes['names']+'\n')
                    ofile.write('# ---------------------------------------------\n')

                ofile.write('  '+str(mass)+'   '+fitRes['values']+'\n')
    
        first_mass_point = False

    return pvalue


    
def run(data):
    #
    # Main function that gets invoked in this module
    # All the functionality is implemented in
    # the underlying Toys class.
    #
    
    legend = '[tprime_toys.run]:'

    toys = Toys(data)

    if action == 'test':
        toys.LoadPdfs('data/toys/pdfs_ejets_smooth_1d_21nov2011.root',
                      'toyWs',
                      pdfNames,
                      varNames)
        _h = toys.Generate2D('mass', 'sig_mass_pdf', 'ht', 'sig_ht_sum_pdf', 10000)
        _h.Draw('LEGO')
        raw_input('hit enter...')
        
        
    elif action == 'toys':
        toys.run()
    elif action == 'limit_2d_toys':
        toys.Limit2dToys()
    elif action == 'generate_2d_toys':
        toys.Generate2dToys()
    elif action == 'fit_pdfs':
        toys.FitPdfs(tprimeMass)
    elif action == 'generate_1d_templates':
        toys.Generate1DTemplates()
    elif action == 'generate_2d_templates':
        toys.Generate2DTemplates()
    


