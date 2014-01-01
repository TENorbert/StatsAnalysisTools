#!/usr/bin/env python
#########################################################################
#
# root2.py
#
# Assorted data manipulations with ROOT
#
# Usage: ./root2.py --help
#
#   Plot and save histogram from a text file
#     ./root2.py --plot-hist -i input.root -o output.root
#
#
# Author: Gena Kukartsev, December 2011
#
#########################################################################
from __future__ import division



########################################
#
# banner
#
def banner():
    print '''
+--------------------------------------------------------------
|
| root2.py
|
| Assorted data manipulations with ROOT
|
| author: Gena Kukartsev, December 2011
|
+--------------------------------------------------------------
    '''
banner()

_legend = '[tprime]:'





########################################
#
# parse command line parameters
#
from optparse import OptionParser
add_help_option = "./root2.py [options]"

parser = OptionParser(add_help_option)

parser.add_option("-t", "--test", dest="test", default=False,
                  help="Test",
                  action="store_true")

parser.add_option("-i", "--in-file", dest="in_file", default=None,
                  action="append",
                  help="Input file names", metavar="INFILE")

parser.add_option("-o", "--out-file", dest="out_file", default=None,
                  help="Output file name")

parser.add_option("--hist", dest="hist",
                  default=None,
                  action="store_true",
                  help="Compute p-value distribution for a set of toy experiments")

parser.add_option("--x-min", dest="xmin",
                  default=None,
                  help="Lower histogram boundary")

parser.add_option("--x-max", dest="xmax",
                  default=None,
                  help="Higher histogram boundary")

parser.add_option("--n-bin", dest="nbin",
                  default=None,
                  help="Number of histogram bins")

parser.add_option("--title", dest="title",
                  default='',
                  help="Histogram title")

parser.add_option("--x-label", dest="xlabel",
                  default='x',
                  help="Label for X axis")

parser.add_option("--y-label", dest="ylabel",
                  default='N',
                  help="Label for Y axis")

parser.add_option("--legend", dest="legend",
                  default=None,
                  help="Legend on the plot")

parser.add_option("--kstest", dest="kstest",
                  default=False,
                  action="store_true",
                  help="do KS test against flat distribution")

parser.add_option("--quantiles", dest="quantiles",
                  default=None,
                  help="Compute quantiles in a column of numbers from file")

parser.add_option("--p-value", dest="pvalue",
                  default=None,
                  help="p-value for quantiles and such")

print _legend, 'parsing command line options...',
(options, args) = parser.parse_args()
print 'done'

if options.hist:
    if options.in_file:
        import hist_plot as hp
        hp.Plot(options.in_file[0], options.out_file,
                options.xmin, options.xmax, options.nbin,
                options.title,
                options.xlabel, options.ylabel,
                options.legend,
                options.kstest)

if options.quantiles:
    # compute quantiles in a column of numbers from file
    if options.in_file:
        import quantiles as q
        q.quantiles(options.in_file[0],
                    int(options.quantiles),
                    options.pvalue)
