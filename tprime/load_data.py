#!/usr/bin/env python
#########################################################################
#
# load_data.py
#
# Data container and associated procedures
#   - creates a pseudo data histogram
#
#
# Usage: import load_data from an analysis script
#
# Author: Gena Kukartsev, May 2011
#
#########################################################################

import copy
from array import array

from ROOT import TFile
from ROOT import TIter
from ROOT import RooWorkspace
from ROOT import RooArgSet
from ROOT import RooArgList
from ROOT import RooRealVar
from ROOT import RooDataHist
from ROOT import RooHistPdf
import ROOT

class TprimeData:

    def __init__(self, lumi = None, objName = 'TprimeData'):

        self._legend = '[TprimeData]:'
        
        self.lumi = lumi         # /pb
        self.hist_names=[]
        #self.hist_ids=[]
        self.hists={}
        #self.hist_id_map={}
        self.object_name = objName
        


    def GetLumi(self):

        return self.lumi


    
    def AddHist(self, hist):

        _name = hist.GetName()

        if _name not in self.hist_names:
            self.hist_names.append(_name)
            self.hists[_name] = copy.deepcopy(hist)
        else:
            self.hists[_name].Add( copy.deepcopy(hist) )

        return


        
    def Rebin1dHist(self, hist,
                    merge_factor, low, high,
                    move_overflow = False):
        #
        # Rebin histogram by merging <merge_factor> adjacent bins.
        # Move overflow and underflow bins into visible range
        #

        legend = '[TprimeData.Rebin1dHist]:'

        if merge_factor == None:
            merge_factor = 1

        if low == None:
            low = hist.GetXaxis().GetXmin()

        if high == None:
            high = hist.GetXaxis().GetXmax()

        merge_factor = int(merge_factor)
        low = float(low)
        high = float(high)

        #hist.Print()
        if 'TH1' not in hist.ClassName():
            #print legend, hist.GetName(), 'is not a TH1, cannot rebin'
            return hist

        #print legend, 'rebinning', hist.GetName()
        
        _hist = hist
        #_hist.Sumw2()

        name = hist.GetName()

        nBinsOrig = _hist.GetNbinsX()
        binWidthOrig = _hist.GetBinWidth(1) # assume uniform bin width
        binWidth = binWidthOrig*float(merge_factor)

        # check range
        if low == None:
            low = _hist.GetXaxis().GetXmin()

        if high == None:
            high = _hist.GetXaxis().GetXmax()

        #nBins = int(nBinsOrig/merge_factor)
        nBins = int((high-low)/binWidth)

        # calculate new binning
        xbins = array('d')
        for i in range(0,nBins+1):
            xbins.append(low+float(i)*binWidth)

        #print xbins
        
        # rebin a copy
        #_rebinned = _hist.Rebin(merge_factor, name)
        _rebinned = _hist.Rebin(nBins, name, xbins)
        #_rebinned.Sumw2()

        # move under- and overflow to visible range
        if move_overflow:
            _err = ROOT.Double()

            # underflow
            _overflow = _rebinned.IntegralAndError(0, 1, _err)
            _rebinned.SetBinContent(0, 0.0)
            _rebinned.SetBinError(0, 0.0)
            _rebinned.SetBinContent(1, _overflow)
            _rebinned.SetBinError(1, _err)
        
            # overflow
            #print legend, 'overflow:', _overflow
            _overflow = _rebinned.IntegralAndError(nBins, nBins+1, _err)
            _rebinned.SetBinContent(nBins+1, 0.0)
            _rebinned.SetBinError(nBins+1, 0.0)
            _rebinned.SetBinContent(nBins, _overflow)
            _rebinned.SetBinError(nBins, _err)

        return _rebinned
        
        

    def RenameHistograms(self, schema):
        #
        # rename selected histograms according to schema
        #

        legend = '[TprimeData.RenameHistograms]:'

        name_pairs = schema.split('|')
        for pair in name_pairs:
            name = pair.split('=')
            print 'trying to rename', name[0], 'to', name[1]

            if name[0] in self.hist_names:
                self.hists[name[0]].SetName(name[1])
                self.hists[name[0]].SetTitle(name[1])
            


        
    def CreatePseudoData(self, name, mapNEvents):
        #
        # Create a histogram that mimics observed data
        # using the existing templates.
        #
        # mapNEvents is a map with a number of events for
        # each template to be generated.
        #

        legend = '[TprimeData.CreatePseudoData]:'

        print legend, 'creating pseudo data according to', mapNEvents

        if name in self.hist_names:
            print legend, 'histogram', name, 'already exists, failed to generate'
            return 0

        if len(mapNEvents) == 0:
            print legend, 'no events requested, nothing will be generated'
            return 0


        for cat in mapNEvents:
            _class = self.hists[cat].ClassName()

            if 'TH1' in _class:
                _nbins = self.hists[cat].GetNbinsX()
                _bin_width = self.hists[cat].GetBinWidth(_nbins)
                _low = self.hists[cat].GetBinLowEdge(1)
                _high = self.hists[cat].GetBinLowEdge(_nbins)+_bin_width
                self.hist_names.append(name)
                self.hists[name] = ROOT.TH1D(name, name, _nbins, _low, _high)

            elif 'TH2' in _class:
                nbinx = self.hists[cat].GetNbinsX() 
                xmin = self.hists[cat].GetXaxis().GetBinLowEdge(1)
                xmax = xmin + nbinx*self.hists[cat].GetBinWidth(1)
                nbiny = self.hists[cat].GetNbinsY() 
                ymin = self.hists[cat].GetYaxis().GetBinLowEdge(1)
                ymax = ymin + nbiny*self.hists[cat].GetBinWidth(1)
                self.hist_names.append(name)
                self.hists[name] = ROOT.TH2D(name, name,
                                             nbinx, xmin, xmax,
                                             nbiny, ymin, ymax)
            break

            
        
        for cat in mapNEvents:
            for i in range(0,mapNEvents[cat]):
                if 'TH2' in _class:
                    _minv = ROOT.Double(0.0)
                    _ht = ROOT.Double(0.0)
                    self.hists[cat].GetRandom2(_minv, _ht)
                    self.hists[name].Fill(_minv, _ht)
                elif 'TH1' in _class:
                    self.hists[name].Fill(self.hists[cat].GetRandom() )

        

    def generate_bw(self, peak, width, low, high, nbins, nev):
        ws = RooWorkspace('ws', 'ws')
        ws.factory('BreitWigner::bw(mass[200, 600], peak[350, 200, 700], width[50])')
        ws.var('mass').setRange(low,high)
        ws.var('peak').setVal(peak)
        ws.var('width').setVal(width)
        _bw = ws.pdf('bw')
        ws.var('mass').setBins(nbins)
        _ds = _bw.generateBinned(RooArgSet(ws.var('mass')), nev)
        _h = _ds.createHistogram('hBw', ws.var('mass'))
        return _h

#    def generate_exp(self, expc, low, high, nbins, nev):
#        ws = RooWorkspace('ws', 'ws')
#        ws.factory('Exponential::bgExp(expc[200, 600], peak[350, 200, 700], width[50])')
#        ws.var('mass').setRange(low,high)
#        ws.var('peak').setVal(peak)
#        _bw = ws.pdf('bw')
#        ws.var('mass').setBins(nbins)
#        _ds = _bw.generateBinned(RooArgSet(ws.var('mass')), nev)
#        _h = _ds.createHistogram('hBw', ws.var('mass'))
#        return _h
#
    def generate_poly(self, b, low, high, nbins, nev):
        ws = RooWorkspace('ws', 'ws')
        ws.factory('Polynomial::bgPol(mass[200, 500], b[-0.001])')
        ws.var('mass').setRange(low,high)
        ws.var('b').setVal(b)
        _pdf = ws.pdf('bgPol')
        ws.var('mass').setBins(nbins)
        _ds = _pdf.generateBinned(RooArgSet(ws.var('mass')), nev)
        _h = _ds.createHistogram('hBgpol', ws.var('mass'))
        return _h

    def generate_hists(self, filename):
        infile = TFile(filename, "recreate")
        _sh = self.generate_poly(-0.001, 200, 500, 50, 1000)
        _sh.SetName('data')
        _sh = self.generate_bw(350, 50, 200, 500, 50, 100000)
        _sh.SetName('sig')
        _sh = self.generate_poly(-0.001, 200, 500, 50, 100000)
        _sh.SetName('bg')
        _sh = self.generate_poly(-0.0015, 200, 500, 50, 100000)
        _sh.SetName('bgNegSyst')
        _sh = self.generate_poly(-0.0005, 200, 500, 50, 100000)
        _sh.SetName('bgPosSyst')
        _sh.Draw('H')
        raw_input('press <enter> to continue...')
        infile.Write()
        infile.Close()

    def load_all_hists(self, filename,
                       prefix = "", suffix = "",
                       zero_overflow = False,
                       projection_x = None,
                       projection_y = None,
                       merge_ewk = False,
                       merge_top = False,
                       strip_suffix = False,
                       rebin_ngroup = None,
                       rebin_low = None,
                       rebin_high = None,
                       move_overflow = False,
                       scale_bin_error = None,
                       excluded_list = []
                       ):
        #
        # merge histograms, works on 1d, 2d...
        #

        legend = '[TprimeData.load_all_hists]:'

        #print legend, excluded_list

        # special prefix for 1D projections
        prefix_1d = 'hist_'

        # these histograms will be merged together under name in ewkName
        #ewk = ['Wjets', 'WW', 'WZ', 'ZZ', 'Zjets', 'SingleToptW', 'SingleTopT', 'SM', 'QCD','SinTop']
        ewk = ['SingleToptW', 'SingleTopT', 'SingleTopS', 'SingleTopbarT', 'SingleTopbarS', 'SinTop', 'Wjets', 'WW', 'WZ', 'ZZ', 'Zjets', 'SM', 'QCD', 'Ewk']
        ewkName = 'ewk'

        # these histograms will be merged together under name in topName
        top = ['TTjets', 'TOP', 'Top', 'top']
        topName = 'top'

        # these histogram names will be changed to tprimeName
        tprime = ['TPrime', 'Tprime', 'tprime', 'TPRIME']
        tprimeName = 'tprime'

        # these histogram names will be changed to dataName
        datanames = ['DATA', 'Data', 'data']
        dataName = 'DATA'

        infile = TFile(filename, "read")
        _keys = TIter(infile.GetListOfKeys())
        _key = _keys.Begin()
        _nhists = 0
        while _key:
            _name = _key.GetName()

            # fix double underscores
#            _new_name = ''
#            _was_underscore = False
#            for _let in _name:
#                if _let=='_':
#                    if _was_underscore:
#                        print 'double underscore fixed'
#                    else:
#                        _was_underscore = True
#                        _new_name += _let
#                else:
#                    _new_name += _let
#                    _was_underscore = False
#
#            _name = _new_name

            if strip_suffix:
                #_name = _name.strip().split('_')[0]
                _newname = _name.strip().split('_')[0]
                if '__jes__plus' in  _name:
                    _newname += '__jes__plus'
                if '__jes__minus' in  _name:
                    _newname += '__jes__minus'
                _name = _newname

            _name = prefix+_name+suffix

            # merge certain histos into one
            if merge_ewk:
                for hname in ewk:
                    #_name = _name.replace(hname+'_', ewkName+'_')
                    _name = _name.replace(hname, ewkName)

            # merge top templates
            if merge_top:
                for hname in top:
                    #_name = _name.replace(hname+'_', topName+'_')
                    _name = _name.replace(hname, topName)

            # uniform signal names
            for hname in tprime:
                _name = _name.replace(hname, tprimeName)

            # uniform data names
            for hname in datanames:
                _name = _name.replace(hname, dataName)

            # skip histos from the excluded list
            if _name in excluded_list:
                _key = _keys.Next()
                continue

            _hist = copy.deepcopy(infile.Get(_key.GetName()))



            # scale error in each bin (in case of multiple use of the same source)
            _nbinx = _hist.GetNbinsX()
            _nbiny = None
            if 'TH2' in _hist.ClassName():
                _nbiny = _hist.GetNbinsY()

            if scale_bin_error:
                print legend, 'histogram errors are multiplied by', float(scale_bin_error)
                for i in range(0, _nbinx+2):
                    for j in range(0, _nbiny+2):
                        _ibin = _hist.GetBin(i,j)
                        _err = _hist.GetBinError(_ibin)
                        _hist.SetBinError(_ibin, _err*float(scale_bin_error))
                
                

            # move over- and underflow to adjacent visible bins
            # fixme: unfinished
            if move_overflow:
                _binx = _hist.GetNbinsX()
                _biny = None
                if 'TH2' in _hist.ClassName():
                    _biny = _hist.GetNbinsY()
                    for i in range(0, _binx+2):
                        for j in range(0, _biny+2):
                            _ibin = GetBin(i,j)

            # project onto 1D if requested and the hist is 2D
            _class_name = _hist.ClassName()
            if 'TH2' in _class_name:
                if projection_x or projection_y:
                    
                    _tmpname = _name+'_tmp'
                    _hist.SetName(_tmpname)

                    if projection_x:
                        _hist = copy.deepcopy(_hist.ProjectionX(_name+'_proj', 0, -1, 'e'))
                    if projection_y:
                        _hist = copy.deepcopy(_hist.ProjectionY(_name+'_proj', 0, -1, 'e'))
                        
                    _hist.SetName(_name)

                    # rebin if requested
                    if rebin_ngroup:
                        self.tmphist = _hist
                        _hist = self.Rebin1dHist(self.tmphist,
                                                 rebin_ngroup, rebin_low, rebin_high,
                                                 move_overflow = True)

            # add the hist to the list and dictionary
            if _name not in self.hist_names:
                self.hist_names.append(_name)
                self.hists[_name] = copy.deepcopy(_hist)
            else:
                self.hists[_name].Add( copy.deepcopy(_hist) )

            _nhists += 1
            self.hists[_name].SetName(_name)

            # set True to remove overflow
            if zero_overflow:
                # zero overflow/underflow bins (assuming that they are already added to the adjacent visible bins)
                # print 'Will zero out overflow and underflow bins!!!'
                
                #print self.hists[_name].ClassName()
                if 'TH2' in self.hists[_name].ClassName():
                    print legend, 'TH2 detected...'
                    print legend, 'Will zero out overflow and underflow bins (no copying)!!!'
                    _binx = self.hists[_name].GetNbinsX()
                    _biny = self.hists[_name].GetNbinsY()
                    for _x in range(0,_binx+2):
                        for _y in range(0,_biny+2):
                            if _x==0 or _x==_binx+1 or _y==0 or _y==_biny+1:
                                # sanity check
                                #if self.hists[_name].GetBinContent(_x,_y) > 0.000001:
                                #    print 'Overflow:', self.hists[_name].GetBinContent(_x,_y), 'set to 0'
                                self.hists[_name].SetBinContent(_x,_y,0.0)
                                self.hists[_name].SetBinError(_x,_y,0.0)
                else:
                    print legend, 'cannot comply with request to zero out over-/underflow for 1D hist'

            _key = _keys.Next()

        # second loop over all hists
#        for name in self.hist_names:
#            # rebin
#            if rebin_ngroup:
#                self.hists[name] = copy.deepcopy(self.Rebin1dHist(self.hists[name],
#                                                                  rebin_ngroup, rebin_low, rebin_high,
#                                                                  move_overflow = True) )
        
        print 'Merged into', len(self.hists), 'histograms'



    def LoadAllHists(self, filename, prefix = "", suffix = "",
                     low = 100, high = 1000, binwidth = 20,
                     setOverflow = True):
        #
        # Same as load_all_hists() but reduces the hist range
        # Range is hardcoded at the moment 100-550
        # Also, trying to merge electroweak category
        #

        # these histograms will be merged together under name in ewkName
        ewk = ['Wjets', 'WW', 'WZ', 'ZZ', 'Zjets', 'SingleToptW', 'SingleTopT']
        ewkName = 'Ewk'
        
        infile = TFile(filename, "read")
        _keys = TIter(infile.GetListOfKeys())
        _key = _keys.Begin()
        _nhists = 0
        while _key:
            _name = _key.GetName()
            _name = str(prefix+_name+suffix)

            # merge certain electroweak histos into one
            for hname in ewk:
                _name = _name.replace(hname+'_', ewkName+'_')

            if _name not in self.hist_names:
                self.hist_names.append(_name)
                
            #_id = tuple(_name.split('_'))
            _nhists += 1
            #self.hist_ids.append(_id)
            _xbins = array('d')
            _low = low
            _high = high
            _width = binwidth
            for ind in range(0,int((_high-_low)/_width +  0.5)+1,1):
                #for binLow in range(_low, int(_high + _width+0.1), _width):
                binLow = _low+ind*_width
                _xbins.append(binLow)
            _hist = copy.deepcopy(infile.Get(_key.GetName()))
            
            #firstOverflowBin = _hist.FindFirstBinAbove(_high-1.5*_width) # including the bin before overflow
            originalBinWidth = _hist.GetBinWidth(1) # assume the input binned uniformly
            firstOverflowBin = int((_high-_low)/originalBinWidth)
            lastOverflowBin = _hist.GetNbinsX()+1 # including the actual overflow
            overflowErr = ROOT.Double()
            _overflow = _hist.IntegralAndError(firstOverflowBin, lastOverflowBin, overflowErr)
            print originalBinWidth, firstOverflowBin, lastOverflowBin, _overflow, overflowErr
            if setOverflow:
                _hist.SetBinContent(firstOverflowBin, _overflow)
                _hist.SetBinError(firstOverflowBin, overflowErr)
            _hist = _hist.Rebin(len(_xbins)-1, 'copyhist', _xbins)

            # if this hist already exists, append to it
            if _name in self.hists:
                self.hists[_name].Add(copy.deepcopy(_hist))
            else:
                self.hists[_name] = copy.deepcopy(_hist)
                self.hists[_name].SetName(_name)

            _key = _keys.Next()
        
        print 'Loaded', len(self.hists), 'histograms'



    def SaveAllHists(self, filename):
        outFile = TFile(filename, "recreate")
        outFile.cd()
        for name in self.hist_names:
            self.hists[name].Write()
        outFile.Close()




    def ConvertToSignalBox(self,
                           xmin = None, xmax = None,
                           ymin = None, ymax = None,
                           signal_only = False):
        #
        # replace all histograms with two-bin histograms:
        # in signal box - right bin
        # out of signal box - left bin
        #

        legend = '[load_data.ConvertToSignalBox]:'

        for name in self.hist_names:

            nbinx = self.hists[name].GetNbinsX()
            lowx = self.hists[name].GetXaxis().GetXmin()
            highx = self.hists[name].GetXaxis().GetXmax()
            binwidthx = (highx-lowx)/nbinx

            nbiny = self.hists[name].GetNbinsY()
            lowy = self.hists[name].GetYaxis().GetXmin()
            highy = self.hists[name].GetYaxis().GetXmax()
            binwidthy = (highy-lowy)/nbiny

            if xmin==None:
                xmin = lowx
            if xmax==None:
                xmax = highx
            if ymin==None:
                ymin = lowy
            if ymax==None:
                ymax = highy


            binx1 = int( (xmin-lowx)/binwidthx ) + 1
            binx2 = int( (xmax-lowx)/binwidthx )
            biny1 = int( (ymin-lowy)/binwidthy ) + 1
            biny2 = int( (ymax-lowy)/binwidthy )

            _name = self.hists[name].GetName()
            _title = self.hists[name].GetTitle()

            if signal_only:
                _hist = ROOT.TH1F(_name, _title, 1, 0, 1)
                _nsig = self.hists[name].Integral(binx1, binx2, biny1, biny2)
                _nbkg = self.hists[name].Integral() - _nsig
                _hist.SetBinContent(1, _nsig)
            else:
                _hist = ROOT.TH1F(_name, _title, 2, 0, 2)
                _nsig = self.hists[name].Integral(binx1, binx2, biny1, biny2)
                _nbkg = self.hists[name].Integral() - _nsig
                _hist.SetBinContent(1, _nbkg)
                _hist.SetBinContent(2, _nsig)

            self.hists[name] = _hist

            #print legend, name, binx1, binx2, biny1, biny2
            print legend, name, _nsig+_nbkg, _nsig, _nbkg

    def ConvertToSignalDiagonal(self,
                                mass = None,
                                signal_only = False):
        for name in self.hist_names:
            
            nbinx = self.hists[name].GetNbinsX()
            lowx = self.hists[name].GetXaxis().GetXmin()
            highx = self.hists[name].GetXaxis().GetXmax()
                        
            nbiny = self.hists[name].GetNbinsY()
            lowy = self.hists[name].GetYaxis().GetXmin()
            highy = self.hists[name].GetYaxis().GetXmax()
                                    
            _name = self.hists[name].GetName()
            _title = self.hists[name].GetTitle()

            _hist = ROOT.TH1F(_name, _title, 2, 0, 2)
            if signal_only:
                _hist = ROOT.TH1F(_name, _title, 1, 0, 1)
            else:
                _hist = ROOT.TH1F(_name, _title, 2, 0, 2)

            _nsig  = 0
            _nbkg  = 0 
            for _ibins in range(0,nbinx+1):
                for _ybins in range(0,nbiny+1):
                    xbincen = self.hists[name].GetXaxis().GetBinCenter(_ibins)
                    ybincen = self.hists[name].GetYaxis().GetBinCenter(_ybins)
                    y_sig_region   = 2*(float(mass)-xbincen)+600.
                    #print y_sig_region,xbincen,ybincen
                    if ybincen > y_sig_region :
                        _nsig += self.hists[name].GetBinContent(_ibins,_ybins)                            
                    else:
                        _nbkg += self.hists[name].GetBinContent(_ibins,_ybins)
            _nbkg = self.hists[name].Integral() - _nsig
            
            
            if signal_only:
                _hist.SetBinContent(1, _nsig)
            else:
                _hist.SetBinContent(1, _nbkg)
                _hist.SetBinContent(2, _nsig)
            
            self.hists[name] = _hist
                
            print self._legend, name, _nsig+_nbkg, _nsig, _nbkg
                





    def ConvertToMassDiagonal(self,
                              mass = None,
                              rebin_ngroup = 5,
                              rebin_low = 100,
                              rebin_high = 700.
                              ):
        for name in self.hist_names:
            
            nbinx = self.hists[name].GetNbinsX()
            lowx = self.hists[name].GetXaxis().GetXmin()
            highx = self.hists[name].GetXaxis().GetXmax()
                        
            nbiny = self.hists[name].GetNbinsY()
            lowy = self.hists[name].GetYaxis().GetXmin()
            highy = self.hists[name].GetYaxis().GetXmax()
                                    
            _name = self.hists[name].GetName()
            _title = self.hists[name].GetTitle()

            _hist = ROOT.TH1F(_name, _title, nbinx, lowx, highx)#, nbiny, lowy, highy)
            for _ibins in range(0,nbinx+1):
                for _ybins in range(0,nbiny+1):
                    xbincen = self.hists[name].GetXaxis().GetBinCenter(_ibins)
                    ybincen = self.hists[name].GetYaxis().GetBinCenter(_ybins)

                    y_sig_minus= xbincen+50
                    #y_sig      = xbincen+350
                    y_sig_plus = xbincen+550
                    
                    if ybincen > y_sig_minus and ybincen < y_sig_plus :
                        #_hist.Fill(xbincen,ybincen, self.hists[name].GetBinContent(_ibins,_ybins))
                        _hist.Fill(xbincen, self.hists[name].GetBinContent(_ibins,_ybins))

            

            self.tmphist = _hist
            
            _hist = self.Rebin1dHist(self.tmphist,
                                     rebin_ngroup, rebin_low, rebin_high,
                                     move_overflow = True)
                                     
            self.hists[name] = _hist
            print self._legend, name, _hist.Integral()

            
    def ConvertToMassDiagPerp(self,
                              mass = None,
                              rebin_ngroup = 5,
                              rebin_low = 100,
                              rebin_high = 700.
                              ):
        for name in self.hist_names:
            
            nbinx = self.hists[name].GetNbinsX()
            lowx = self.hists[name].GetXaxis().GetXmin()
            highx = self.hists[name].GetXaxis().GetXmax()
            
            nbiny = self.hists[name].GetNbinsY()
            lowy = self.hists[name].GetYaxis().GetXmin()
            highy = self.hists[name].GetYaxis().GetXmax()
            
            _name = self.hists[name].GetName()
            _title = self.hists[name].GetTitle()
            
            _hist = ROOT.TH1F(_name, _title, nbinx, lowx, highx)#, nbiny, lowy, highy)
            for _ibins in range(0,nbinx+1):
                for _ybins in range(0,nbiny+1):
                    xbincen = self.hists[name].GetXaxis().GetBinCenter(_ibins)
                    ybincen = self.hists[name].GetYaxis().GetBinCenter(_ybins)
                        
                    y_sig_minus   = 2*(float(mass)-xbincen)+600.
                    y_sig_plus   = 2*(float(mass)-xbincen)+1200.
                    
                    if ybincen > y_sig_minus and ybincen < y_sig_plus :
                        #_hist.Fill(xbincen,ybincen, self.hists[name].GetBinContent(_ibins,_ybins))
                        _hist.Fill(xbincen, self.hists[name].GetBinContent(_ibins,_ybins))
                        
            self.tmphist = _hist
            _hist = self.Rebin1dHist(self.tmphist,
                                     rebin_ngroup, rebin_low, rebin_high, move_overflow = True)
            
            self.hists[name] = _hist
            print self._legend, name, _hist.Integral()
                            
                            
    def get_position(self, name, pos):
        #
        # gets position value out of name
        #
        _id = tuple(name.split('_'))
        if len(_id)>=pos:
            return _id[pos-1]
        

    def processes(self):
        procs = set()
        for name in self.hist_names:
            procs.add(self.get_position(name,1))

        for p in procs:
            print p

        return procs

    def variables(self):
        _s = set()
        for name in self.hist_names:
            _s.add(self.get_position(name,2))

        for p in _s:
            print p

        return _s

    def find_filled_interval(self, hist):
        _nbins = hist.GetNbinsX()
        _len = 0
        _running_len = 0
        _low = 1
        _running_low = 1
        for i in range(1,_nbins+1):
            if hist.Integral() < 0.000001:
                continue
            _bin = hist.GetBinContent(i)/hist.Integral()*hist.GetEntries()
            if _bin > 0.9:
                _running_len += 1
            else:
                if _running_len > _len:
                    _len = _running_len
                    _low = _running_low
                _running_len = 0
                _running_low = i+1

        if _running_len > _len:
            _len = _running_len
            _low = _running_low
        return tuple([_low, _len])

    def list_process(self, process):
        self.list(process)
        
    def list(self, pattern):
        _pats = pattern.split('*')
        histnames = set()
        print 'hist_name'.rjust(40),
        print 'n_entries'.rjust(9),
        print 'integral'.rjust(12),
        print 'with overflow'.rjust(16),
        print 'n_bins'.rjust(7),
        print 'low'.rjust(6),
        print 'high'.rjust(6)
        #print 'min_bin'.rjust(8),
        #print 'min_content'.rjust(10),
        #print 'bin_range'.rjust(10)
        print '-------------------------------------------------------------------------------------------------------'
        for name in self.hist_names:
            #_p = self.get_position(name,1)
            _match = True
            for _pat in _pats:
                if name.find(_pat)==-1:
                    _match=False
            if _match:
                histnames.add(name)

                #for name in histnames:
                _entries = self.hists[name].GetEntries()
                _nbins = self.hists[name].GetNbinsX()
                _integral = self.hists[name].Integral()
                _integral_over = None
                if 'TH1' in self.hists[name].ClassName():
                    _integral_over = self.hists[name].Integral(0,_nbins+1)
                elif 'TH2' in self.hists[name].ClassName():
                    _nbinsY = self.hists[name].GetNbinsY()
                    _integral_over = self.hists[name].Integral(0,_nbins+1, 0,_nbinsY+1)
                _bin_width = self.hists[name].GetBinWidth(_nbins)
                _low = self.hists[name].GetBinLowEdge(1)
                _high = self.hists[name].GetBinLowEdge(_nbins)+_bin_width
                _min_bin = self.hists[name].GetMinimumBin()
                _min_bin_content = self.hists[name].GetBinContent(_min_bin)
                _interval = self.find_filled_interval(self.hists[name])
                _low_range = _interval[0]
                _high_range = _interval[0]+_interval[1]-1
                print name.rjust(40), repr(_entries).rjust(9),
                print '{0:.3f}'.format(_integral).rjust(12),
                print '{0:.3f}'.format(_integral_over).rjust(16),
                print repr(_nbins).rjust(7),
                print '{0:.1f}'.format(_low).rjust(6),
                print '{0:.1f}'.format(_high).rjust(6)
                #print repr(_min_bin).rjust(8),
                #print repr(_min_bin_content).rjust(10),
                #print '[{0:3d}, {1:3d}]'.format(_low_range,_high_range)

        return histnames

    def generate_hist(self, hist, varname, n):
        _var = RooRealVar(varname, varname, hist.GetMean())
        _dh = RooDataHist("dh", "dh", RooArgList(_var), hist)
        _pdf = RooHistPdf("pdf", "pdf", RooArgSet(_var), _dh)
        _ds = _pdf.generateBinned(RooArgSet(_var), n)
        _h = _ds.createHistogram('hist', _var)
        #_h.Draw()
        #raw_input('press <enter> to continue...')
        return copy.deepcopy(_h)

        

    def reduce(self, outname):
        _ofile = TFile(outname, "recreate")
        
        h_sig = copy.deepcopy(self.hists['TPrime350_fitMass_4j_1t'])
        h_sig.SetName("signal")
        h_sig.SetTitle("signal")
        _ofile.Append(h_sig)
        
        h_sig_jesup = copy.deepcopy(self.hists['TPrime350_fitMass_4j_1t_JESup'])
        h_sig_jesup.SetName("signal_jesup")
        h_sig_jesup.SetTitle("signal_jesup")
        _ofile.Append(h_sig_jesup)
        
        h_sig_jesdown = copy.deepcopy(self.hists['TPrime350_fitMass_4j_1t_JESdown'])
        h_sig_jesdown.SetName("signal_jesdown")
        h_sig_jesdown.SetTitle("signal_jesdown")
        _ofile.Append(h_sig_jesdown)
        
        h_bg1 = copy.deepcopy(self.hists['Top_fitMass_4j_1t'])
        h_bg1.SetName("bg1")
        h_bg1.SetTitle("bg1")
        _ofile.Append(h_bg1)
        
        h_bg1_jesup = copy.deepcopy(self.hists['Top_fitMass_4j_1t_JESup'])
        h_bg1_jesup.SetName("bg1_jesup")
        h_bg1_jesup.SetTitle("bg1_jesup")
        _ofile.Append(h_bg1_jesup)
        
        h_bg1_jesdown = copy.deepcopy(self.hists['Top_fitMass_4j_1t_JESdown'])
        h_bg1_jesdown.SetName("bg1_jesdown")
        h_bg1_jesdown.SetTitle("bg1_jesdown")
        _ofile.Append(h_bg1_jesdown)
        
        h_bg2 = copy.deepcopy(self.hists['Wjets_fitMass_4j_1t'])
        h_bg2.SetName("bg2")
        h_bg2.SetTitle("bg2")
        _ofile.Append(h_bg2)
        
        h_bg2_jesup = copy.deepcopy(self.hists['Wjets_fitMass_4j_1t_JESup'])
        h_bg2_jesup.SetName("bg2_jesup")
        h_bg2_jesup.SetTitle("bg2_jesup")
        _ofile.Append(h_bg2_jesup)
        
        h_bg2_jesdown = copy.deepcopy(self.hists['Wjets_fitMass_4j_1t_JESdown'])
        h_bg2_jesdown.SetName("bg2_jesdown")
        h_bg2_jesdown.SetTitle("bg2_jesdown")
        _ofile.Append(h_bg2_jesdown)
        
        #h_data = copy.deepcopy(self.hists['Data_fitMass_4j_1t'])
        h_data = self.generate_hist(h_bg1, 'mass', 52)
        h_data.Add(self.generate_hist(h_bg2, 'mass', 10))
        #h_data.Add(self.generate_hist(h_sig, 'mass', 5))
        h_data.Print()
        h_data.SetName("data")
        h_data.SetTitle("data")
        _ofile.Append(h_data)
        
        _ofile.Write()
        _ofile.Close()

        return
