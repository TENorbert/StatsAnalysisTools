#!/usr/bin/env python

################################################
#
# This script plots limits for the t't'bar
# semileptonic analysis with CMS
#
#   - input data from ASCII files
#   
#
# Gena Kukartsev, 2011
#
################################################

import ROOT
from plotter import *
import os,re

def cls_limits(dir, format = 'png'):
    # - takes condor output dir as input
    # - read input root files from condor dir
    # - merge the results for each mass point and get CLs numbers
    # - make the xsec limit plot as function of mass
    
    legend = '[cls_limits]:'

    print legend, 'harvesting', dir, 'for CLs results'

    # find and sort by mass individual CLs outputs for merging
    mass_points={}
    for root,dirs,files in os.walk(dir):
        for file in files:
            #print file
            if file[len(file)-5:len(file)] != '.root':
                continue

            # find root file that corresponds to a mass point
            m = re.search('[0-9]+.[0-9]',file)
            if m:
                mass = m.group(0)
                if mass not in mass_points:
                    mass_points[mass] = []
                mass_points[mass].append(file)

    # write CLs output files to be merged in a file by mass
    os.system('root -l -b -q -n hf_tprime.C+')
    os.system('rm '+dir.rstrip('/')+'/hf_tprime_C.so' )
    os.system('rm '+dir.rstrip('/')+'/*ascii.ascii' )
    os.system('cp hf_tprime_C.so '+dir.rstrip('/')+'/')
    os.system('cp run_limit2.C '+dir.rstrip('/')+'/')
    os.chdir(dir)
    for mass in mass_points:
        merge_file = 'merge_'+mass+'.ascii'
        with open(merge_file,'w') as file:
            for filename in mass_points[mass]:
                file.write(filename+'\n')

        # merge CLs per mass point
        os.system('root -l -b -q -n run_limit2.C\(' + \
                  '\\"combined\\",' + \
                  '\\"observed\\",' + \
                  '\\"cls\\",' + \
                  str(mass)+',' + \
                  '\\"'+merge_file+'\\",' + \
                  str(1)+',' + \
                  str(1)+',' + \
                  str(0)+',' + \
                  str(1)+',' + \
                  str(1)+',' + \
                  str(1)+',' + \
                  '\\"\\"' + \
                  '\)')

    os.system('../merge_ascii.py *ascii.ascii > for_limit_plot.ascii')

    os.chdir('../')


    

def tprime_cls_limits(dir,lumi,channel,plot_format = 'png'):
    #
    # CLs limits
    #    

    data={}
    
    ROOT.gStyle.SetHatchesLineWidth(2)
    ROOT.gStyle.SetLineStyleString(11,"60 30");

    ###############################################################################
    #
    # mu+jets
    # CLs observed and expected limits
    #    

    ROOT.gStyle.SetHatchesLineWidth(2)
    ROOT.gStyle.SetLineStyleString(11,"60 30");

    data[4] = Data( float(lumi) )

    _scale = 1.0
    _scale_theory = 1.0
    
    data[4].add_with_errors('RS_exp2', dir.rstrip('/')+'/for_limit_plot.ascii',
                '95% expected',
                1.0,
                ROOT.kYellow, 1, 1,
                ROOT.kYellow, 0.1, 8,
                fill_style = 1002,
                index_x = 0,
                index_y = 4,
                index_ey_up = 6,
                index_ey_down = 2)
    data[4].add_with_errors('RS_exp1', dir.rstrip('/')+'/for_limit_plot.ascii',
                '68% expected',
                1.0,
                ROOT.kGreen, 1, 1,
                ROOT.kGreen, 0.1, 8,
                fill_style = 1002,
                index_x = 0,
                index_y = 4,
                index_ey_up = 5,
                index_ey_down = 3)
    data[4].add_with_errors('obs_exp_median', dir.rstrip('/')+'/for_limit_plot.ascii',
                'median expected',
                1.0,
                ROOT.kBlack, 2, 3,
                ROOT.kBlack, 0.1, 8,
                fill_style = 1002,
                index_x = 0,
                index_y = 4)
    data[4].add_with_errors('obs', dir.rstrip('/')+'/for_limit_plot.ascii',
                'data',
                1.0,
                ROOT.kBlack, 1, 3,
                ROOT.kBlack, 1.2, 8,
                fill_style = 1002,
                index_x = 0,
                index_y = 1)
    data[4].add('theory', 'tprime_theory.ascii', 't\'_{THEORY}',
                1.0,
                ROOT.kBlue, 1, 3,
                ROOT.kBlue, 0.1, 8)

    makeMultiGraph(data[4], 'limit_plot.'+plot_format,
                   400, 600, 0, 4.0,
                   xlabel = "M_{t'} [GeV]", ylabel = '#sigma [pb]',
                   xLegend = .55, yLegend = .55,
                   legendWidth = 0.45, legendHeight = 0.30,
                   fillStyle = 1002,
                   drawOption = 'APL3', letter = '',
                   toplegend = '#font[42]{CL_{S}: '+channel+'}')





def tprime_cls_limits_comp(dir,dir2,lumi,channel,plot_format = 'png'):
    #
    # CLs limits
    #    

    data={}
    
    ROOT.gStyle.SetHatchesLineWidth(2)
    ROOT.gStyle.SetLineStyleString(11,"60 30");

    ###############################################################################
    #
    # mu+jets
    # CLs observed and expected limits
    #    

    ROOT.gStyle.SetHatchesLineWidth(2)
    ROOT.gStyle.SetLineStyleString(11,"60 30");

    data[4] = Data( float(lumi) )

    _scale = 1.0
    _scale_theory = 1.0
    
    data[4].add_with_errors('RS_exp2', dir.rstrip('/')+'/for_limit_plot.ascii',
                '95% expected',
                1.0,
                ROOT.kYellow, 1, 1,
                ROOT.kYellow, 0.1, 8,
                fill_style = 1002,
                index_x = 0,
                index_y = 4,
                index_ey_up = 6,
                index_ey_down = 2)
    data[4].add_with_errors('RS_exp1', dir.rstrip('/')+'/for_limit_plot.ascii',
                '68% expected',
                1.0,
                ROOT.kGreen, 1, 1,
                ROOT.kGreen, 0.1, 8,
                fill_style = 1002,
                index_x = 0,
                index_y = 4,
                index_ey_up = 5,
                index_ey_down = 3)
    data[4].add_with_errors('obs_exp_median', dir.rstrip('/')+'/for_limit_plot.ascii',
                'median expected',
                1.0,
                ROOT.kBlack, 2, 3,
                ROOT.kBlack, 0.1, 8,
                fill_style = 1002,
                index_x = 0,
                index_y = 4)
    data[4].add_with_errors('obs', dir.rstrip('/')+'/for_limit_plot.ascii',
                'data',
                1.0,
                ROOT.kBlack, 1, 3,
                ROOT.kBlack, 1.2, 8,
                fill_style = 1002,
                index_x = 0,
                index_y = 1)
   
#2nd directory
    data[4].add_with_errors('RS_exp42', dir2.rstrip('/')+'/for_limit_plot.ascii',
                '95% expected 2',
                1.0,
                ROOT.kMagenta, 1, 1,
                ROOT.kMagenta, 0.1, 8,
                fill_style = 0000,
                index_x = 0,
                index_y = 4,
                index_ey_up = 6,
                index_ey_down = 2)
    data[4].add_with_errors('RS_exp32', dir2.rstrip('/')+'/for_limit_plot.ascii',
                '68% expected 2',
                1.0,
                ROOT.kOrange, 1, 1,
                ROOT.kOrange, 1.2, 8,
                fill_style = 0000,
                index_x = 0,
                index_y = 4,
                index_ey_up = 5,
                index_ey_down = 3)
    data[4].add_with_errors('obs_exp_median2', dir2.rstrip('/')+'/for_limit_plot.ascii',
                'median expected 2',
                1.0,
                ROOT.kMagenta, 2, 3,
                ROOT.kMagenta, 1.2, 8,
                fill_style = 0000,
                index_x = 0,
                index_y = 4)
    data[4].add_with_errors('obs2', dir2.rstrip('/')+'/for_limit_plot.ascii',
                'data 2',
                1.0,
                ROOT.kBlue, 1, 3,
                ROOT.kBlue, 1.2, 8,
                fill_style = 0000,
                index_x = 0,
                index_y = 1)

#    data[4].add('theory', 'tprime_theory.ascii', 't\'_{THEORY}',
#                1.0,
#                ROOT.kBlue, 1, 3,
#                ROOT.kBlue, 0.1, 8)


    makeMultiGraph(data[4], 'limit_plot.'+plot_format,
                   400, 600, 0, 4.0,
                   xlabel = "M_{t'} [GeV]", ylabel = '#sigma [pb]',
                   xLegend = .55, yLegend = .55,
                   legendWidth = 0.45, legendHeight = 0.30,
                   fillStyle = 1002,
                   drawOption = 'APL3', letter = '',
                   toplegend = '#font[42]{CL_{S}: '+channel+'}')

