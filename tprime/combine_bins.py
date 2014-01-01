#!/usr/bin/env python
#########################################################################
#
# combine_bins.py
#
# Module for mapping 2D templates on 1D
#
# Usage:
#     
#     ./tprime.py -a combine_bins --mass 450.0 --scale-factor 1.0 --error 0.05 -f input.root -o output.root
#     ./tprime.py --combine-bins signal_box_1bin --mass 450.0 --channel mujets -f input_merged_mu.root -o output.root
#
# Author: Gena Kukartsev, November 2011
#
#########################################################################

import sys
import copy
import ROOT


binMap = {}

def loadBinMap(filename):
    #
    # load 2d->1d map from an ascii file
    # (created by combineBins routine)
    #
    # copied from tprime_toys.py
    #
    
    legend = '[combine_bins.loadBinMap]:'

    binMap.clear()
    
    with open(filename) as mFile:

        lines=mFile.readlines()
        
        bin = 1
        for line in lines:
            _columns = line.strip().split(':')
            
            if _columns[0][0] == '#':
                continue
            
            binMap[bin] = list(map(int,_columns[1].split()))
            bin += 1

    print legend, 'map for', len(binMap), 'bins loaded'

    return binMap



def CombineBinsWithMapOneHist(hist):
    #
    # Use a 2D->1D map and convert a 2D hist
    # into 1D
    #
    # Format of the map: list of bin lists that correspond
    # to each combined bin
    #
    # copied from tprime_toys.py
    #
    
    legend = '[combine_bins.CombineBinsWithMapOneHist]:'
        
    hist1d = ROOT.TH1F('combHist', 'combHist',
                       len(binMap),0, len(binMap))

    for i in binMap:
        
        for bin in binMap[i]:
            
            hist1d.AddBinContent(i,hist.GetBinContent(bin))

    hist1d.SetEntries(hist.GetEntries())
    return hist1d



def CombineBinsWithMap(binMapFileName,
                       infile, outfile):
    #
    # Use a 2D->1D map and convert a 2D hist
    # into 1D
    #
    # Format of the map: list of bin lists that correspond
    # to each combined bin
    #

    legend = '[combine_bins.CombineBinsWithMapOneHist]:'

    hists = {}
    hist_names = []

    loadBinMap(binMapFileName)

    inputFile = ROOT.TFile(infile, "read")
    _keys = ROOT.TIter(inputFile.GetListOfKeys())
    _key = _keys.Begin()
    _nhists = 0
    while _key:
        _name = _key.GetName()

        _hist = copy.deepcopy(inputFile.Get(_name))

        hist_names.append(_name)

        _rebinned_hist = CombineBinsWithMapOneHist(_hist)
        hists[_name] = copy.deepcopy(_rebinned_hist)
        hists[_name].SetName(_name)
        

        _key = _keys.Next()

    inputFile.Close()

    # save rebinned hists
    outputFile = ROOT.TFile(outfile, "recreate")
    outputFile.cd()
    for name in hist_names:
        hists[name].Write()
    outputFile.Close()

        
    return



def GangBins(channel, mass, scale,
             maxErr, inFileName, outFileName, prefix):
    #
    # Gang bins Davis-style
    #

    legend = '[combine_bins.GangBins]:'

    sigHistName = 'tprime'+str(mass)
    bgHistName = ['top', 'ewk']
    extraHistName = ['DATA',
                     sigHistName+'__jes__plus', sigHistName+'__jes__minus',
                     'top__jes__plus', 'top__jes__minus',
                     'ewk__jes__plus', 'ewk__jes__minus']


    inFile = ROOT.TFile(inFileName, 'read')

    sigHist = copy.deepcopy( inFile.Get(sigHistName) )
    bgHist = None

    for name in bgHistName:
        if bgHist == None:
            bgHist = copy.deepcopy( inFile.Get(name) )
        else:
            bgHist.Add( inFile.Get(name) )

    inFile.Close()

    #sigHist.Draw('LEGO')
    #raw_input('hit enter...')
    #bgHist.Draw('LEGO')
    #raw_input('hit enter...')

    sbHist = copy.deepcopy( sigHist )
    sbHist.SetName('sb_ratio')
    sbHist.SetTitle('sb_ratio')
    sbHist.Divide(bgHist)
    sbHist.Draw('LEGO')
    raw_input('hit enter...')



def ProjectOnMass(infile, outfile, prefix = ''):
    
    from load_data import TprimeData
    _hists = TprimeData(1)

    _hists.load_all_hists(infile, projection_x = True,
                          rebin_ngroup = 5,
                          rebin_low = 100,
                          rebin_high = 700,
                          prefix=prefix)
    _hists.SaveAllHists(outfile)

    return



def ProjectOnHt(infile, outfile, prefix = ''):
    from load_data import TprimeData
    _hists = TprimeData(1)
    _hists.load_all_hists(infile, projection_y = True,
                          rebin_ngroup = 4,
                          rebin_low    = 200,
                          rebin_high   = 2000,
                          prefix=prefix)
    _hists.SaveAllHists(outfile)

    return


def MakeSignalBoxTemplates(infile, outfile, prefix = '',
                           xmin = None, xmax = None,
                           ymin = None, ymax = None,
                           signal_only = False):
    
    from load_data import TprimeData
    _hists = TprimeData(1)

    _hists.load_all_hists(infile,
                          prefix=prefix)

    _hists.ConvertToSignalBox(xmin, xmax, ymin, ymax, signal_only=signal_only)
    
    _hists.SaveAllHists(outfile)

    return

def MakeSignalBoxDiagonal(infile, outfile, prefix = '',
                          mass = None,
                          signal_only = False):

    from load_data import TprimeData
    _hists = TprimeData(1)
    
    _hists.load_all_hists(infile,
                          prefix=prefix)
    
    _hists.ConvertToSignalDiagonal(mass,signal_only=signal_only)
    
    _hists.SaveAllHists(outfile)
    
    return

def MakeSignalMassDiagonal(infile, outfile, prefix = '',
                           mass = None,
                           rebin_ngroup = None,
                           rebin_low = None,
                           rebin_high= None
                           ):

    from load_data import TprimeData
    _hists = TprimeData(1)
    
    _hists.load_all_hists(infile,
                          prefix=prefix)
    
    _hists.ConvertToMassDiagonal(mass,rebin_ngroup,rebin_low,rebin_high)
    
    _hists.SaveAllHists(outfile)
    
    return
                        
def MakeSignalMassDiagPerp(infile, outfile, prefix = '',
                           mass = None,
                           rebin_ngroup = None,
                           rebin_low = None,
                           rebin_high= None
                           ):
    
        from load_data import TprimeData
        _hists = TprimeData(1)
        
        _hists.load_all_hists(infile,
                              prefix=prefix)
        
        _hists.ConvertToMassDiagPerp(mass,rebin_ngroup,rebin_low,rebin_high)

        _hists.SaveAllHists(outfile)
        
        return
                        

def CopyHists(infile, outfile, prefix = ''):
    
    from load_data import TprimeData
    _hists = TprimeData(1)

    _hists.load_all_hists(infile, prefix)
    _hists.SaveAllHists(outfile)

    return



def CombineBins(mass, scale, maxErr,
                infile, outfile,
                prefix = '',
                method = 'brown',
                make_plots = False,
                maxErrSig = 10.0,
                channel = "ljets"):
    #
    # Main function that gets invoked in this module
    # Main mapping is loaded from combineBins.C
    #
    
    legend = '[combine_bins]:'

    # FIXME: dummy now, need to remove
    #channel = 'ljets'

    if channel==None:
        print legend, 'no channel specified, exiting'
        return

    if "brown" in method:
        if method == 'brown_smooth':
            print legend, 'combine bins with Brown smooth algo (EPS2011), using S/B'
            ROOT.gROOT.ProcessLine('.L combineBinsSmooth.C+')
            from ROOT import combineBins
            print '####################', maxErrSig, float(maxErrSig)
            print channel, int(mass), float(scale), float(maxErr), infile, outfile, prefix, make_plots, 0.2
            combineBins(channel,
                        int(mass),
                        float(scale),
                        float(maxErr),
                        infile, outfile, prefix,
                        make_plots,
                        0.2)
                        #float(maxErrSig))
        
        elif method == 'brown':
            print legend, 'combine bins with Brown algo (EPS2011), using S/B'
            ROOT.gROOT.ProcessLine('.L combineBins.C+')
            from ROOT import combineBins
            combineBins(channel, int(mass), float(scale),
                        float(maxErr), infile, outfile, prefix,
                        make_plots,
                        float(maxErrSig))

        elif method == 'brown_2d':
            print legend, 'combine bins with Brown 2D algo, using S/B'
            ROOT.gROOT.ProcessLine('.L combineBins_2d.C+')
            from ROOT import combineBins
            combineBins(channel, int(mass), float(scale),
                        float(maxErr), infile, outfile, prefix,
                        make_plots,
                        float(maxErrSig))
            
        elif method == 'brown_mike':
            ROOT.gROOT.ProcessLine('.L brownNxNGang.C+')
            print legend, 'combine bins with 4x4 + Brown algo (EPS2011), using S/B'
            from ROOT import rebin
            rebin( int(mass), 
                   float(scale),
                   float(maxErr),
                   infile,
                   outfile,
                   make_plots)
    elif method == 'davis':
        print legend, 'combine bins with Davis algo'
        ROOT.gROOT.ProcessLine('.L davisGang.C+')
        from ROOT import rebin
        rebin( int(mass),
               float(scale),
               float(maxErr),
               infile,
               outfile)

               
    elif method == 'copy':
        print legend, 'just copy all the histograms unchanged'
        CopyHists(infile, outfile, prefix)
    
    elif method == 'mass':
        print legend, 'combine bins projecting on mass axis and rebinning'
        ProjectOnMass(infile, outfile, prefix)
        from load_data import TprimeData
        _hists = TprimeData(1)
    
    elif method == 'ht':
        print legend, 'combine bins projecting on HT axis and rebinning'
        ProjectOnHt(infile, outfile, prefix)
        from load_data import TprimeData
        _hists = TprimeData(1)
    
    elif method == 'signal_box_1bin':
        print legend, 'combine bins in and out a signal box'
        MakeSignalBoxTemplates(infile, outfile, prefix,
                               xmin=float(mass)*0.9, xmax=float(mass)*1.1,
                               ymin = 700.0, ymax = 2000.0,
                               signal_only = True)
    
    elif method == 'signal_box_2bin':
        print legend, 'combine bins in and out a signal box'
        MakeSignalBoxTemplates(infile, outfile, prefix,
                               xmin=float(mass)*0.9, xmax=float(mass)*1.1,
                               ymin = 700.0, ymax = 2000.0,
                               signal_only = False)

    elif method == 'signal_box_diagonal':
        print legend, 'combine bins 1d diagonally'
        MakeSignalBoxDiagonal(infile, outfile, prefix,
                              float(mass),
                              signal_only = False)

    elif method == 'signal_mass_diagonal':
        MakeSignalMassDiagonal(infile, outfile, prefix,
                               float(mass),
                               rebin_ngroup = 5,
                               rebin_low = 100,
                               rebin_high= 700
                               )

    elif method == 'signal_mass_diag_perp':
        MakeSignalMassDiagPerp(infile, outfile, prefix,
                               float(mass),
                               rebin_ngroup = 5,
                               rebin_low = 100,
                               rebin_high= 700
                               )
    else:
        print legend, 'no valid method for combine bins specified...'
        sys.exit(1)
