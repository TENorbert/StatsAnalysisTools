#!/usr/bin/env python
#########################################################################
#
# tprime_ws.py
#
# Workspace operations for tprime analysis
#
# Usage: from tprime_ws import * from an analysis script
#
# Author: Gena Kukartsev, May 2011
#
#########################################################################

import copy

from ROOT import TFile
from ROOT import TIter
from ROOT import RooWorkspace
from ROOT import RooArgSet
from ROOT import RooArgList
from ROOT import RooRealVar
from ROOT import RooDataHist
from ROOT import RooHistPdf



def ws_inspect(filename, name):
    legend = '[tprime: WS inspect]:'
    infile = TFile(filename, "read")
    ws = infile.Get("combined")
    
    if name == 'all':
        ws.var("obs_tprime").setVal(200)
        ws.Print("v")

    else:
        _var = ws.var("obs_tprime")
        _var.Print()
        _frame = _var.frame(100)
        _f = ws.function(name)
        _f.plotOn(_frame)
        _frame.Draw()
        
    raw_input('press <enter> to continue...')

    return
