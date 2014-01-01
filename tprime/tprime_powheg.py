#!/usr/bin/env python
#########################################################################
#
# tprime_powheg.py
#
# Checks and manipulations with bin maps
#
# Usage: ./tprime.py --powheg
#
#
# Author: Gena Kukartsev, February 2012
#
#########################################################################
from __future__ import division

import ROOT
import cms_prel

import copy

from load_data import TprimeData

def CreateTemplates(options):
    #
    # create two "standard" files with "plus" and "minus"
    # shape templates for POWHEG-Nominal difference
    # Since there is only one POWHEG template, the idea is
    # to interpolate in some way between it and the nominal one
    #

    legend = '[tprime_powheg.CreateTemplates]:'

    channels = ['e','mu']

    input_nom = {}
    input_pow = {}

    input_nom['mu'] = 'data/mujets_4600/from_gueorgi/27jan2011v1/csv/v1/mujets_4601ipb_2D_v1_csv_jer_nom.root'
    input_pow['mu'] = 'data/mujets_4600/from_gueorgi/27jan2011v1/csv/v1/mujets_4601ipb_2D_v1_csv_POWHEG.root'

    input_nom['e'] = 'data/ejets_4700/ricardo_4683ipb_2d_csv_22jan2012v2/ht4jetsAfterFit_vs_fitMass.root'
    input_pow['e'] = 'data/ejets_4700/ricardo_4683ipb_2d_csv_22jan2012v2/ht4jetsAfterFit_vs_fitMass_powheg.root'

    hist_name = {}
    hist_name['e']  = 'Top_ht4jetsAfterFit:fitMass'
    hist_name['mu'] = 'TTjets_HtvsMfit'

    
    for chan in channels:

        data = TprimeData()
        data_powheg = TprimeData()
        
        _suffix = ''
        
        data.load_all_hists(input_nom[chan],
                            merge_ewk=False,
                            merge_top=False,
                            strip_suffix = False,
                            suffix = _suffix
                            )
        
        data_powheg.load_all_hists(input_pow[chan],
                                   merge_ewk=False,
                                   merge_top=False,
                                   strip_suffix = False,
                                   suffix = _suffix
                                   )
        
        hPlus = copy.deepcopy(data.hists[ hist_name[chan] ])
        hMinus = copy.deepcopy(data_powheg.hists[ hist_name[chan] ])

        hPlus.SetMarkerColor(ROOT.kRed)
        
        #hMinus.Draw()
        #hPlus.Draw('same')
        
        #raw_input('press enter to continue')
    
        hNom = copy.deepcopy(hMinus)
        hNom.Add(hPlus)
        hNom.Scale(0.5)


        # save nominal
        
        data.hists.pop( hist_name[chan] )
        iTop = 0
        for name in data.hist_names:
            if name == hist_name[chan]:
                break
            iTop += 1
        data.hist_names.pop(iTop)

        data.AddHist(hNom)

        data.SaveAllHists(chan+'_powheg_nom.root')

    
        # save minus
        
        data.hists.pop( hist_name[chan] )
        iTop = 0
        for name in data.hist_names:
            if name == hist_name[chan]:
                break
            iTop += 1
        data.hist_names.pop(iTop)

        data.AddHist(hMinus)

        data.SaveAllHists(chan+'_powheg_minus.root')

    
        # save plus
        
        data.hists.pop( hist_name[chan] )
        iTop = 0
        for name in data.hist_names:
            if name == hist_name[chan]:
                break
            iTop += 1
        data.hist_names.pop(iTop)

        data.AddHist(hPlus)

        data.SaveAllHists(chan+'_powheg_plus.root')

    
