#!/usr/bin/env python
#########################################################################
#
# quantiles.py
#
# Calculate quantiles
#
# Usage: ./root2.py --quantiles <column_number> -i infile.txt
#
#
# Author: Gena Kukartsev, December 2011
#
#########################################################################
from __future__ import division

import ROOT
import cms_prel

from array import array

def quantiles(in_fname, column, pvalue = None):

    if pvalue:
        pvalue = float(pvalue)

    legend = '[quantiles.quantiles]:'

    _num = array('d')
    with open(in_fname, 'r') as infile:
        for line in infile:
            _words = line.strip().split()
            if len(_words)>0:
                if _words[0][0]=='#':
                    continue
                _num.append(float(_words[column]))

    #print legend, _num

    _prob = array('d')
    _quantiles = array('d')
    if pvalue:
        _prob.append(pvalue)
        _quantiles.append(0)
        _nprob = 1
    else:
        _prob.append(0.021)
        _prob.append(0.159)
        _prob.append(0.5)
        _prob.append(0.841)
        _prob.append(0.979)
        _quantiles.append(0)
        _quantiles.append(0)
        _quantiles.append(0)
        _quantiles.append(0)
        _quantiles.append(0)
        _nprob = 5

    ROOT.TMath.Quantiles(len(_num), _nprob, _num, _quantiles, _prob, False)

    for i in range(0,_nprob):
        print legend, str(_prob[i])+'-quantile =', _quantiles[i]

    return _quantiles

#                _prob.append(0.159)
#                _prob.append(0.841)
#                _prob.append(0.979)
#                _quantiles.append(0)
#                _quantiles.append(0)
#                _quantiles.append(0)

