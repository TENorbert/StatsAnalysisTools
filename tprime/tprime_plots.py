#!/usr/bin/env python

################################################
#
# This script plots limits for the t't'bar
# semileptonic analysis with CMS
#
#   - input data from ASCII files
#
# 1,2,3 - nominal limit plots
# 33    - POWHEG syst for the AN
# 26,24 - expected syst plots for AN
#
# Gena Kukartsev, 2011
#
################################################

import ROOT
from plotter import *


ymax = 10
ymin = 0.07
suffix = 'diagperp_28may12'
#directory = 'comb'
#22mar2012vM'


#for uly
#directory3 = '1889_4j'
#directory5 = '1889_5j'
#directory = '4jv0_1889'
#directory4 = '5jv0_1889'

emu = 'ejets'

#directory1=  'comb_01feb2012v2'
#directory2=  'comb_22mar2012vM'




#directory3=  '16feb12_10pc' #wewk_reb14feb12'
#directory3 = 'reb14feb12_reb450x450'
plot_observed = True

def tprime_plots(plot_format = 'png',
                 plot_number = [1],
                 directory1   = 'comb',
                 directory2  = '',
                 channel     = ''
                 ):
    #
    # Make plots with specified plot numbers
    #


    #------------------------------------------------------------>
    #1
    # e+jets
    # Asymptotic CLs limit
    # brown 20 (s and b)
    #    

    if str(1) in plot_number:

        channel = 'e+jets (5.0fb^{-1}) mass/diag2'

        ROOT.gStyle.SetHatchesLineWidth(2)
        ROOT.gStyle.SetLineStyleString(11,"60 30");
        
        data = Data(4700)
        
        _scale = 1.0
        _scale_theory = 1.0

        #_dir = 'e_brown_20'
        #_dir = 'e_btag'
        #_dir = 'ejets_19jan2012v1'
        #_dir = 'ejets_19jan2012v2'
        #_dir = 'ejets_22jan2012v4'

        _dir = directory1
        #'ejets_'+directory
        #_dir = 'ejets_10feb2012v2'
        #_dir = 'ejets_21mar2012v2' # paper 2011
        #_dir = 'e_2bin_21may2012v1'

        
        data.add_with_errors('RD_exp1', _dir+'/limits.txt',
                             '#pm2#sigma expected',
                             1.0,
                             ROOT.kYellow, 1, 3,
                             ROOT.kYellow, 0.1, 8,
                             fill_style = 1002,
                             index_x = 0,
                             index_y = 4,
                             index_ey_up = 6,
                             index_ey_down = 2,
                             legend_index = 4)
        data.add_with_errors('RD_exp2', _dir+'/limits.txt',
                             '#pm1#sigma expected',
                             1.0,
                             ROOT.kGreen, 1, 3,
                             ROOT.kGreen, 0.1, 8,
                             fill_style = 1002,
                             index_x = 0,
                             index_y = 4,
                             index_ey_up = 5,
                             index_ey_down = 3,
                             legend_index = 3)
        if plot_observed:
            data.add_with_errors('obs', _dir+'/limits.txt',
                                 'observed 95% C.L.',
                                 1.0,
                                 ROOT.kBlack, 1, 3,
                                 ROOT.kBlack, 1.2, 8,
                                 fill_style = 1002,
                                 index_x = 0,
                                 index_y = 1,
                                 legend_index = 1)
            
        data.add_with_errors('obs2', _dir+'/limits.txt',
                             'expected',
                             1.0,
                             ROOT.kBlack, 2, 3,
                             ROOT.kBlack, 0.1, 2,
                             index_x = 0,
                             index_y = 4,
                             legend_index = 2)
        
        data.add_with_errors('obs3', 'ascii/tprime_theory.ascii', 't\'_{THEORY}',
                             1.0,
                             ROOT.kBlue, 1, 3,
                             ROOT.kBlue, 0.1, 8,
                             index_x = 0,
                             index_y = 1,
                             legend_index = 5)


        makeMultiGraph(data, 'ejets_limit_5.0ifb_'+suffix+'.'+plot_format,
                       400, 625, ymin, ymax,
                       xlabel = "M_{t'} [GeV]", ylabel = '#sigma [pb]',
                       xLegend = .55, yLegend = .55,
                       legendWidth = 0.45, legendHeight = 0.30,
                       fillStyle = 1002,
                       drawOption = 'APL3', letter = '',
                       toplegend = '#font[22]{CL_{S}: '+channel+'}',
                       cms_prel = False)



    #------------------------------------------------------------>
    #2
    # mu+jets
    # Asymptotic CLs limit
    # brown 20 (s and b)
    #    

    if str(2) in plot_number:

        #channel = '#mu+jets (4.9fb^{-1})'
        channel = '#mu+jets (4.9fb^{-1}) mass/diag2'
    
        ROOT.gStyle.SetHatchesLineWidth(2)
        ROOT.gStyle.SetLineStyleString(11,"60 30");
        
        data = Data(4700)
        
        _scale = 1.0
        _scale_theory = 1.0
        

        #_dir = 'mu_brown_20'
        #_dir = 'mu_btag'
        #_dir = 'mujets_19jan2012v1'
        #_dir = 'mujets_19jan2012v2'
        #_dir = 'mujets_22jan2012v4'

        _dir = directory1
        #'mujets_'+directory

        data.add_with_errors('RD_exp1', _dir+'/limits.txt',
                             '#pm2#sigma expected',
                             1.0,
                             ROOT.kYellow, 1, 3,
                             ROOT.kYellow, 0.1, 8,
                             fill_style = 1002,
                             index_x = 0,
                             index_y = 4,
                             index_ey_up = 6,
                             index_ey_down = 2,
                             legend_index = 4)
        data.add_with_errors('RD_exp2', _dir+'/limits.txt',
                             '#pm1#sigma expected',
                             1.0,
                             ROOT.kGreen, 1, 3,
                             ROOT.kGreen, 0.1, 8,
                             fill_style = 1002,
                             index_x = 0,
                             index_y = 4,
                             index_ey_up = 5,
                             index_ey_down = 3,
                             legend_index = 3)
        if plot_observed:
            data.add_with_errors('obs', _dir+'/limits.txt',
                                 'observed 95% C.L.',
                                 1.0,
                                 ROOT.kBlack, 1, 3,
                                 ROOT.kBlack, 1.2, 8,
                                 fill_style = 1002,
                                 index_x = 0,
                                 index_y = 1,
                                 legend_index = 1)
            
        data.add_with_errors('obs2', _dir+'/limits.txt',
                             'expected',
                             1.0,
                             ROOT.kBlack, 2, 3,
                             ROOT.kBlack, 0.1, 2,
                             index_x = 0,
                             index_y = 4,
                             legend_index = 2)
        
        data.add_with_errors('obs3', 'ascii/tprime_theory.ascii', 't\'_{THEORY}',
                             1.0,
                             ROOT.kBlue, 1, 3,
                             ROOT.kBlue, 0.1, 8,
                             index_x = 0,
                             index_y = 1,
                             legend_index = 5)


        makeMultiGraph(data, 'mujets_limit_4.9ifb_'+suffix+'.'+plot_format,
                       400, 625, ymin, ymax,
                       xlabel = "M_{t'} [GeV]", ylabel = '#sigma [pb]',
                       xLegend = .55, yLegend = .55,
                       legendWidth = 0.45, legendHeight = 0.30,
                       fillStyle = 1002,
                       drawOption = 'APL3', letter = '',
                       toplegend = '#font[22]{CL_{S}: '+channel+'}',
                       cms_prel = False)



    #------------------------------------------------------------>
    #3
    # combined
    # Asymptotic CLs limit
    # brown 20 (s and b)
    #    

    if str(3) in plot_number:

        #        channel = '#mu+jets (4.90fb^{-1}), e+jets (4.98fb^{-1})'
        channel = '#mu+jets (4.9fb^{-1}), e+jets (5.0fb^{-1}) mass/diag2'
    
        ROOT.gStyle.SetHatchesLineWidth(2)
        ROOT.gStyle.SetLineStyleString(11,"60 30");
        
        data = Data(4700)
        
        _scale = 1.0
        _scale_theory = 1.0
        

        #_dir = 'comb_brown_20'
        #_dir = 'comb_btag'
        #_dir = 'comb_test'
        #_dir = 'comb_19jan2012v1'
        #_dir = 'comb_19jan2012v2'
        #_dir = 'comb_22jan2012v4'

        _dir = directory1
        #'comb_'+directory

        data.add_with_errors('RD_exp1', _dir+'/limits.txt',
        #data.add_with_errors('RD_exp1', 'test_comb_brown_20/limits.txt',
                             '#pm2#sigma expected',
                             1.0,
                             ROOT.kYellow, 1, 3,
                             ROOT.kYellow, 0.1, 8,
                             fill_style = 1002,
                             index_x = 0,
                             index_y = 4,
                             index_ey_up = 6,
                             index_ey_down = 2,
                             legend_index = 4)
        data.add_with_errors('RD_exp2', _dir+'/limits.txt',
        #data.add_with_errors('RD_exp2', 'test_comb_brown_20/limits.txt',
                             '#pm1#sigma expected',
                             1.0,
                             ROOT.kGreen, 1, 3,
                             ROOT.kGreen, 0.1, 8,
                             fill_style = 1002,
                             index_x = 0,
                             index_y = 4,
                             index_ey_up = 5,
                             index_ey_down = 3,
                             legend_index = 3)

        if plot_observed:
            data.add_with_errors('obs', _dir+'/limits.txt',
            #data.add_with_errors('obs', 'test_comb_brown_20/limits.txt',
                                 'observed 95% C.L.',
                                 1.0,
                                 ROOT.kBlack, 1, 3,
                                 ROOT.kBlack, 1.2, 8,
                                 fill_style = 1002,
                                 index_x = 0,
                                 index_y = 1,
                                 legend_index = 1)
            
        data.add_with_errors('obs2', _dir+'/limits.txt',
        #data.add_with_errors('obs2', 'test_comb_brown_20/limits.txt',
                             'expected',
                             1.0,
                             ROOT.kBlack, 2, 3,
                             ROOT.kBlack, 0.1, 2,
                             index_x = 0,
                             index_y = 4,
                             legend_index = 2)
        
        data.add_with_errors('obs3', 'ascii/tprime_theory.ascii', 't\'_{THEORY}',
                             1.0,
                             ROOT.kBlue, 1, 3,
                             ROOT.kBlue, 0.1, 8,
                             index_x = 0,
                             index_y = 1,
                             legend_index = 5)

        makeMultiGraph(data, 'comb_limit_'+suffix+'.'+plot_format,
                       400, 625, ymin,ymax,
                       #400, 625, 0.0, 1.0,
                       xlabel = "M_{t'} [GeV]", ylabel = '#sigma [pb]',
                       xLegend = .55, yLegend = .55,
                       legendWidth = 0.45, legendHeight = 0.30,
                       fillStyle = 1002,
                       drawOption = 'APL3', find_intersection = True,#draw_limit_line = True,
                       letter = '',
                       toplegend = '#font[22]{CL_{S}: '+channel+'}',
                       cms_prel = False)


    


    #------------------------------------------------------------>
    #3001
    # comparison of ejets, mujets and combined
    # Asymptotic CLs limit
    # brown 20 (s and b)
    #    
    if str(444) in plot_number:

        channel = emu #''#mu+jets (4.6fb^{-1}), e+jets (4.7fb^{-1})'
        
        ROOT.gStyle.SetHatchesLineWidth(2)
        ROOT.gStyle.SetLineStyleString(11,"60 30");
        
        data = Data(4700)
        
        _scale = 1.0
        _scale_theory = 1.0
        
        if plot_observed:
            #data.add_with_errors('obs', 'ejets_19jan2012v2/limits.txt',
            data.add_with_errors('obs', directory1+'/limits.txt',
                                 'observed '+emu+'+jets - 1',
                                 1.0,
                                 ROOT.kBlack, 1, 3,
                                 ROOT.kBlack, 1.2, 8,
                                 fill_style = 1002,
                                 index_x = 0,
                                 index_y = 1)
            #data.add_with_errors('obs2', 'mujets_19jan2012v2/limits.txt',
            data.add_with_errors('obs2', directory2+'/limits.txt',
                                 'observed '+emu+'+jets - 2',
                                 1.0,
                                 ROOT.kRed, 1, 3,
                                 ROOT.kRed, 1.2, 8,
                                 fill_style = 1002,
                                 index_x = 0,
                                 index_y = 1)

            data.add_with_errors('obs4', directory1+'/limits.txt',
                                 'expected '+emu+'+jets - 1',
                                 1.0,
                                 ROOT.kBlack, 2, 3,
                                 ROOT.kBlack, 0.1, 2,
                                 index_x = 0,
                                 index_y = 4)
            data.add_with_errors('obs5', directory2+'/limits.txt',
                                 'expected '+emu+'+jets - 2',
                                 1.0,
                                 ROOT.kRed, 2, 3,
                                 ROOT.kRed, 0.1, 2,
                                 index_x = 0,
                                 index_y = 4)
            
            data.add('theory', 'ascii/tprime_theory.ascii', 't\'_{THEORY}',
                     1.0,
                     ROOT.kBlue, 1, 3,
                     ROOT.kBlue, 0.1, 8)
            makeMultiGraph(data, 'cls_limit_'+emu+'_5ifb.'+plot_format,
                           300, 700, ymin, ymax,
                           xlabel = "M_{t'} [GeV]", ylabel = '#sigma [pb]',
                           xLegend = .55, yLegend = .55,
                           legendWidth = 0.45, legendHeight = 0.30,
                           fillStyle = 1002,
                           drawOption = 'APL3', letter = '',
                           toplegend = '#font[22]{CL_{S}: '+channel+'}')
            
                           
    if str(3001) in plot_number:

        channel = '#mu+jets (4.6fb^{-1}), e+jets (4.7fb^{-1})'
    
        ROOT.gStyle.SetHatchesLineWidth(2)
        ROOT.gStyle.SetLineStyleString(11,"60 30");
        
        data = Data(4700)
        
        _scale = 1.0
        _scale_theory = 1.0
        
        if plot_observed:
            #data.add_with_errors('obs', 'ejets_19jan2012v2/limits.txt',
            data.add_with_errors('obs', 'ejets_01feb2012v1/limits.txt',
                                 'observed e+jets',
                                 1.0,
                                 ROOT.kBlue, 1, 3,
                                 ROOT.kBlue, 1.2, 8,
                                 fill_style = 1002,
                                 index_x = 0,
                                 index_y = 1)
            #data.add_with_errors('obs2', 'mujets_19jan2012v2/limits.txt',
            data.add_with_errors('obs2', 'mujets_01feb2012v1/limits.txt',
                                 'observed #mu+jets',
                                 1.0,
                                 ROOT.kRed, 1, 3,
                                 ROOT.kRed, 1.2, 8,
                                 fill_style = 1002,
                                 index_x = 0,
                                 index_y = 1)
            #data.add_with_errors('obs3', 'comb_19jan2012v2/limits.txt',
            data.add_with_errors('obs3', 'comb_01feb2012v1/limits.txt',
                                 'observed combined',
                                 1.0,
                                 ROOT.kBlack, 1, 3,
                                 ROOT.kBlack, 1.2, 8,
                                 fill_style = 1002,
                                 index_x = 0,
                                 index_y = 1)
        
        data.add_with_errors('obs4', 'ejets_01feb2012v1/limits.txt',
                             'expected e+jets',
                             1.0,
                             ROOT.kBlue, 2, 3,
                             ROOT.kBlue, 0.1, 2,
                             index_x = 0,
                             index_y = 4)
        data.add_with_errors('obs5', 'mujets_01feb2012v1/limits.txt',
                             'expected #mu+jets',
                             1.0,
                             ROOT.kRed, 2, 3,
                             ROOT.kRed, 0.1, 2,
                             index_x = 0,
                             index_y = 4)
        data.add_with_errors('obs6', 'comb_01feb2012v1/limits.txt',
                             'expected combined',
                             1.0,
                             ROOT.kBlack, 2, 3,
                             ROOT.kBlack, 0.1, 2,
                             index_x = 0,
                             index_y = 4)
        
        data.add('theory', 'ascii/tprime_theory.ascii', 't\'_{THEORY}',
                 1.0,
                 ROOT.kBlue, 1, 3,
                 ROOT.kBlue, 0.1, 8)

        makeMultiGraph(data, 'cls_limit_comp_5ifb.'+plot_format,
                       400, 625, 0.05, 5.0,
                       xlabel = "M_{t'} [GeV]", ylabel = '#sigma [pb]',
                       xLegend = .55, yLegend = .55,
                       legendWidth = 0.45, legendHeight = 0.30,
                       fillStyle = 1002,
                       drawOption = 'APL3', letter = '',
                       toplegend = '#font[22]{CL_{S}: '+channel+'}')



    #------------------------------------------------------------>
    #1001
    # e+jets
    # Asymptotic CLs limit
    # brown smooth 10 (s and b)
    #    

    if str(1001) in plot_number:

        channel = 'e+jets (4.7fb^{-1})'

        ROOT.gStyle.SetHatchesLineWidth(2)
        ROOT.gStyle.SetLineStyleString(11,"60 30");
        
        data = Data(4700)
        
        _scale = 1.0
        _scale_theory = 1.0

        _dir = 'ejets_smooth10_01feb2012v1'
        
        data.add_with_errors('RD_exp1', _dir+'/limits.txt',
                             '#pm2#sigma expected',
                             1.0,
                             ROOT.kYellow, 1, 3,
                             ROOT.kYellow, 0.1, 8,
                             fill_style = 1002,
                             index_x = 0,
                             index_y = 4,
                             index_ey_up = 6,
                             index_ey_down = 2)
        data.add_with_errors('RD_exp2', _dir+'/limits.txt',
                             '#pm1#sigma expected',
                             1.0,
                             ROOT.kGreen, 1, 3,
                             ROOT.kGreen, 0.1, 8,
                             fill_style = 1002,
                             index_x = 0,
                             index_y = 4,
                             index_ey_up = 5,
                             index_ey_down = 3)
        if plot_observed:
            data.add_with_errors('obs', _dir+'/limits.txt',
                                 'observed 95% C.L.',
                                 1.0,
                                 ROOT.kBlack, 1, 3,
                                 ROOT.kBlack, 1.2, 8,
                                 fill_style = 1002,
                                 index_x = 0,
                                 index_y = 1)
            
        data.add_with_errors('obs2', _dir+'/limits.txt',
                             'expected',
                             1.0,
                             ROOT.kBlack, 2, 3,
                             ROOT.kBlack, 0.1, 2,
                             index_x = 0,
                             index_y = 4)
        
        data.add('theory', 'ascii/tprime_theory.ascii', 't\'_{THEORY}',
                 1.0,
                 ROOT.kBlue, 1, 3,
                 ROOT.kBlue, 0.1, 8)

        makeMultiGraph(data, 'ejets_limit_smooth10_4.7ifb_'+suffix+'.'+plot_format,
                       400, 625, 0.05, 5.0,
                       xlabel = "M_{t'} [GeV]", ylabel = '#sigma [pb]',
                       xLegend = .55, yLegend = .55,
                       legendWidth = 0.45, legendHeight = 0.30,
                       fillStyle = 1002,
                       drawOption = 'APL3', letter = '',
                       toplegend = '#font[42]{CL_{S}: '+channel+'}')



    #------------------------------------------------------------>
    #1002
    # mu+jets
    # Asymptotic CLs limit
    # brown smooth 10 (s and b)
    #    

    if str(1002) in plot_number:

        channel = '#mu+jets (4.6fb^{-1})'
    
        ROOT.gStyle.SetHatchesLineWidth(2)
        ROOT.gStyle.SetLineStyleString(11,"60 30");
        
        data = Data(4700)
        
        _scale = 1.0
        _scale_theory = 1.0
        
        _dir = 'mujets_smooth10_01feb2012v1'

        data.add_with_errors('RD_exp1', _dir+'/limits.txt',
                             '#pm2#sigma expected',
                             1.0,
                             ROOT.kYellow, 1, 3,
                             ROOT.kYellow, 0.1, 8,
                             fill_style = 1002,
                             index_x = 0,
                             index_y = 4,
                             index_ey_up = 6,
                             index_ey_down = 2)
        data.add_with_errors('RD_exp2', _dir+'/limits.txt',
                             '#pm1#sigma expected',
                             1.0,
                             ROOT.kGreen, 1, 3,
                             ROOT.kGreen, 0.1, 8,
                             fill_style = 1002,
                             index_x = 0,
                             index_y = 4,
                             index_ey_up = 5,
                             index_ey_down = 3)
        if plot_observed:
            data.add_with_errors('obs', _dir+'/limits.txt',
                                 'observed 95% C.L.',
                                 1.0,
                                 ROOT.kBlack, 1, 3,
                                 ROOT.kBlack, 1.2, 8,
                                 fill_style = 1002,
                                 index_x = 0,
                                 index_y = 1)
            
        data.add_with_errors('obs2', _dir+'/limits.txt',
                             'expected',
                             1.0,
                             ROOT.kBlack, 2, 3,
                             ROOT.kBlack, 0.1, 2,
                             index_x = 0,
                             index_y = 4)
        
        data.add('theory', 'ascii/tprime_theory.ascii', 't\'_{THEORY}',
                 1.0,
                 ROOT.kBlue, 1, 3,
                 ROOT.kBlue, 0.1, 8)

        makeMultiGraph(data, 'mujets_limit_smooth10_4.6ifb_'+suffix+'.'+plot_format,
                       400, 625, 0.05, 5.0,
                       xlabel = "M_{t'} [GeV]", ylabel = '#sigma [pb]',
                       xLegend = .55, yLegend = .55,
                       legendWidth = 0.45, legendHeight = 0.30,
                       fillStyle = 1002,
                       drawOption = 'APL3', letter = '',
                       toplegend = '#font[42]{CL_{S}: '+channel+'}')



    #------------------------------------------------------------>
    #1003
    # combined
    # Asymptotic CLs limit
    # brown smooth 10 (s and b)
    #    

    if str(1003) in plot_number:

        channel = '#mu+jets (4.6fb^{-1}), e+jets (4.7fb^{-1})'
    
        ROOT.gStyle.SetHatchesLineWidth(2)
        ROOT.gStyle.SetLineStyleString(11,"60 30");
        
        data = Data(4700)
        
        _scale = 1.0
        _scale_theory = 1.0
        
        _dir = 'comb_smooth10_01feb2012v1'

        data.add_with_errors('RD_exp1', _dir+'/limits.txt',
        #data.add_with_errors('RD_exp1', 'test_comb_brown_20/limits.txt',
                             '#pm2#sigma expected',
                             1.0,
                             ROOT.kYellow, 1, 3,
                             ROOT.kYellow, 0.1, 8,
                             fill_style = 1002,
                             index_x = 0,
                             index_y = 4,
                             index_ey_up = 6,
                             index_ey_down = 2)
        data.add_with_errors('RD_exp2', _dir+'/limits.txt',
        #data.add_with_errors('RD_exp2', 'test_comb_brown_20/limits.txt',
                             '#pm1#sigma expected',
                             1.0,
                             ROOT.kGreen, 1, 3,
                             ROOT.kGreen, 0.1, 8,
                             fill_style = 1002,
                             index_x = 0,
                             index_y = 4,
                             index_ey_up = 5,
                             index_ey_down = 3)

        if plot_observed:
            data.add_with_errors('obs', _dir+'/limits.txt',
            #data.add_with_errors('obs', 'test_comb_brown_20/limits.txt',
                                 'observed 95% C.L.',
                                 1.0,
                                 ROOT.kBlack, 1, 3,
                                 ROOT.kBlack, 1.2, 8,
                                 fill_style = 1002,
                                 index_x = 0,
                                 index_y = 1)
            
        data.add_with_errors('obs2', _dir+'/limits.txt',
        #data.add_with_errors('obs2', 'test_comb_brown_20/limits.txt',
                             'expected',
                             1.0,
                             ROOT.kBlack, 2, 3,
                             ROOT.kBlack, 0.1, 2,
                             index_x = 0,
                             index_y = 4)
        
        data.add('theory', 'ascii/tprime_theory.ascii', 't\'_{THEORY}',
                 1.0,
                 ROOT.kBlue, 1, 3,
                 ROOT.kBlue, 0.1, 8)

        makeMultiGraph(data, 'comb_limit_'+suffix+'.'+plot_format,
                       400, 625, ymin, ymax,
                       xlabel = "M_{t'} [GeV]", ylabel = '#sigma [pb]',
                       xLegend = .55, yLegend = .55,
                       legendWidth = 0.45, legendHeight = 0.30,
                       fillStyle = 1002,
                       drawOption = 'APL3', letter = '',
                       toplegend = '#font[42]{CL_{S}: '+channel+'}')



    #------------------------------------------------------------>
    #4
    # mu+jets algo comparison
    # Asymptotic expected CLs limit
    #    

    if str(4) in plot_number:

        channel = '#mu+jets, expected 95% C.L. limit'
    
        ROOT.gStyle.SetHatchesLineWidth(2)
        ROOT.gStyle.SetLineStyleString(11,"60 30");
        
        data = Data(4600)
        
        _scale = 1.0
        _scale_theory = 1.0
        data.add_with_errors('obs', 'mujets_'+directory+'/limits.txt',
                             'No btag in fit - 4 Jets',
                             1.0,
                             ROOT.kRed, 2, 3,
                             ROOT.kRed, 0.1, 2, 
                             fill_style = 1001,
                             index_x = 0,
                             index_y = 4)
        data.add_with_errors('obs2', 'mujets_'+directory3+'/limits.txt',
                             'Nominal - 4 Jets',
                             1.0,
                             ROOT.kBlue+2, 1, 3,
                             ROOT.kBlue+2, 0.1, 8,
                             fill_style = 1001,
                             index_x = 0,
                             index_y = 4)
        data.add_with_errors('obs3', 'mujets_'+directory4+'/limits.txt',
                             'No btag in fit - 5 Jets',
                             1.0,
                             ROOT.kGreen+2, 1, 3,
                             ROOT.kGreen+2, 1, 8,
                             fill_style = 1001,
                             index_x = 0,
                             index_y = 4)
        data.add_with_errors('obs4', 'mujets_'+directory5+'/limits.txt',
                             'Nominal - 5 Jets',
                             1.0,
                             ROOT.kBlack, 1, 3,
                             ROOT.kBlack, 1, 8,
                             fill_style = 1001,
                             index_x = 0,
                             index_y = 4)
        
        # data.add_with_errors('obs2', 'mujets_'+directory2+'/limits.txt',
        #                     '180x200 GeV Binning',
        #                     1.0,
        #                     ROOT.kBlue+2, 1, 3,
        #                     ROOT.kBlue+2, 0.1, 8,
        #                     fill_style = 1002,
        #                     index_x = 0,
        #                     index_y = 4)
        #        data.add_with_errors('obs3', 'mujets_'+directory3+'/limits.txt',
        #                     'Nominal - 5 Jets',
        #                     1.0,
        #                     ROOT.kGreen+2, 1, 3,
        #                     ROOT.kGreen+2, 1, 8,
        #                     fill_style = 1002,
        #                     index_x = 0,
        #                     index_y = 4)

        #        data.add_with_errors('obs4', 'mujets_'+directory4+'/limits.txt',
        #                     'adjacent, 20% B',
        #                     1.0,
        #                     ROOT.kRed, 2, 3,
        #                     ROOT.kRed, 0.1, 2,
        #                     index_x = 0,
        #                     index_y = 4)
        
        data.add('theory', 'ascii/tprime_theory.ascii', 't\'_{THEORY}',
                 1.0,
                 ROOT.kBlue, 1, 3,
                 ROOT.kBlue, 0.1, 8)


        makeMultiGraph(data, 'mu_1889ifb_'+suffix+'.'+plot_format,
                       400, 800, ymin, ymax,
                       xlabel = "M_{t'} [GeV]", ylabel = '#sigma [pb]',
                       xLegend = .55, yLegend = .55,
                       legendWidth = 0.45, legendHeight = 0.30,
                       fillStyle = 1002,
                       drawOption = 'APL3', letter = '',
                       toplegend = '#font[42]{CL_{S}: '+channel+'}')



    #------------------------------------------------------------>
    #41
    # mu+jets algo comparison
    # Asymptotic expected CLs limit
    #    

    if str(41) in plot_number:

        channel = 'e+jets, expected 95% C.L. limit'
    
        ROOT.gStyle.SetHatchesLineWidth(2)
        ROOT.gStyle.SetLineStyleString(11,"60 30");
        
        data = Data(4700)
        
        _scale = 1.0
        _scale_theory = 1.0

        data.add_with_errors('obs2', 'ejets_'+directory+'/limits.txt',
                             'S/B Nominal 20%',
                             1.0,
                             ROOT.kBlack, 1, 3,
                             ROOT.kBlack, 0.1, 8,
                             fill_style = 1002,
                             index_x = 0,
                             index_y = 4)
        #        data.add_with_errors('obs2', 'ejets_'+directory2+'/limits.txt',
        #                     '180x200 GeV Binning',
        #                     1.0,
        #                     ROOT.kBlue+2, 1, 3,
        #                     ROOT.kBlue+2, 0.1, 8,
        #                     fill_style = 1002,
        #                     index_x = 0,
        #                     index_y = 4)
        data.add_with_errors('obs3', 'ejet_'+directory3+'/limits.txt',
                             'S/B Nominal 10%',
                             1.0,
                             ROOT.kGreen+2, 1, 3,
                             ROOT.kGreen+2, 0.1, 8,
                             fill_style = 1002,
                             index_x = 0,
                             index_y = 4)
        #data.add_with_errors('obs4', 'ejets_adj20_01feb2012v1/limits.txt',
        #                    'adjacent, 20% B',
        #                     1.0,
        #                     ROOT.kRed, 2, 3,
        #                     ROOT.kRed, 0.1, 2,
        #                     index_x = 0,
        #                     index_y = 4)
        
        data.add('theory', 'ascii/tprime_theory.ascii', 't\'_{THEORY}',
                 1.0,
                 ROOT.kBlue, 1, 3,
                 ROOT.kBlue, 0.1, 8)


        makeMultiGraph(data, 'e_sensitivity_5ifb_'+suffix+'.'+plot_format,
                       400, 625, 0.05, 5.0,
                       xlabel = "M_{t'} [GeV]", ylabel = '#sigma [pb]',
                       xLegend = .55, yLegend = .55,
                       legendWidth = 0.45, legendHeight = 0.30,
                       fillStyle = 1002,
                       drawOption = 'APL3', letter = '',
                       toplegend = '#font[42]{CL_{S}: '+channel+'}')



    #------------------------------------------------------------>
    #42
    # mu+jets + e+jets algo comparison
    # Asymptotic expected CLs limit
    #    

    if str(42) in plot_number:

        channel = '#mu,e+jets, expected 95% C.L. limit'
    
        ROOT.gStyle.SetHatchesLineWidth(2)
        ROOT.gStyle.SetLineStyleString(11,"60 30");
        
        data = Data()
        
        _scale = 1.0
        _scale_theory = 1.0

        data.add_with_errors('obs2', 'combined_'+directory+'/limits.txt',
                             'S/B Nominal 20%',
                             1.0,
                             ROOT.kBlack, 1, 3,
                             ROOT.kBlack, 0.1, 8,
                             fill_style = 1002,
                             index_x = 0,
                             index_y = 4)
        #        data.add_with_errors('obs2', 'comb_'+directory2+'/limits.txt',
        #                     '180x200 GeV Binning',
        #                     1.0,
        #                     ROOT.kBlue+2, 1, 3,
        #                     ROOT.kBlue+2, 0.1, 8,
        #                     fill_style = 1002,
        #                     index_x = 0,
        #                     index_y = 4)
        data.add_with_errors('obs3', 'combined_'+directory3+'/limits.txt',
                             'S/B Nominal 10%',
                             1.0,
                             ROOT.kGreen+2, 1, 3,
                             ROOT.kGreen+2, 0.1, 8,
                             fill_style = 1002,
                             index_x = 0,
                             index_y = 4)
        #        data.add_with_errors('obs4', 'comb_adj20_01feb2012v1/limits.txt',
        #                     'adjacent, 20% B',
        #                     1.0,
        #                     ROOT.kRed, 2, 3,
        #                     ROOT.kRed, 0.1, 2,
        #                     index_x = 0,
        #                     index_y = 4)
        
        data.add('theory', 'ascii/tprime_theory.ascii', 't\'_{THEORY}',
                 1.0,
                 ROOT.kBlue, 1, 3,
                 ROOT.kBlue, 0.1, 8)


        makeMultiGraph(data, 'comb_sensitivity_5ifb_'+suffix+'.'+plot_format,
                       400, 625, 0.05, 5.0,
                       xlabel = "M_{t'} [GeV]", ylabel = '#sigma [pb]',
                       xLegend = .55, yLegend = .55,
                       legendWidth = 0.45, legendHeight = 0.30,
                       fillStyle = 1002,
                       drawOption = 'APL3', letter = '',
                       toplegend = '#font[42]{CL_{S}: '+channel+'}')



    #------------------------------------------------------------>
    #5
    # mu+jets algo comparison
    # Asymptotic expected CLs limit
    # 40%
    #    

    if str(5) in plot_number:

        channel = '#mu+jets, expected 95% C.L. limit'
    
        ROOT.gStyle.SetHatchesLineWidth(2)
        ROOT.gStyle.SetLineStyleString(11,"60 30");
        
        data = Data(4600)
        
        _scale = 1.0
        _scale_theory = 1.0
        
        data.add_with_errors('obs1', 'mu_brown_40/limit_test.txt',
                             'S/B ordering, 40%',
                             1.0,
                             ROOT.kBlack, 1, 3,
                             ROOT.kBlack, 0.1, 8,
                             fill_style = 1002,
                             index_x = 0,
                             index_y = 4)
        data.add_with_errors('obs2', 'mu_brown_40b/limit_test.txt',
                             'S/B ordering, 40% bg',
                             1.0,
                             ROOT.kGreen+2, 1, 3,
                             ROOT.kGreen+2, 0.1, 8,
                             fill_style = 1002,
                             index_x = 0,
                             index_y = 4)
        data.add_with_errors('obs3', 'mu_davis_40b/limit_test.txt',
                             'adjacent bins, 40% bg',
                             1.0,
                             ROOT.kRed, 2, 3,
                             ROOT.kRed, 0.1, 2,
                             index_x = 0,
                             index_y = 4)
        
        data.add('theory', 'ascii/tprime_theory.ascii', 't\'_{THEORY}',
                 1.0,
                 ROOT.kBlue, 1, 3,
                 ROOT.kBlue, 0.1, 8)

        makeMultiGraph(data, 'mu_sensitivity_40pc_4600ipb.'+plot_format,
                       400, 600, 0.05, 5.0,
                       xlabel = "M_{t'} [GeV]", ylabel = '#sigma [pb]',
                       xLegend = .55, yLegend = .55,
                       legendWidth = 0.45, legendHeight = 0.30,
                       fillStyle = 1002,
                       drawOption = 'APL3', letter = '',
                       toplegend = '#font[42]{CL_{S}: '+channel+'}')



    #------------------------------------------------------------>
    #6
    # e+jets algo comparison
    # Asymptotic expected CLs limit
    # 20%
    #    

    if str(6) in plot_number:

        channel = 'e+jets, expected 95% C.L. limit'
    
        ROOT.gStyle.SetHatchesLineWidth(2)
        ROOT.gStyle.SetLineStyleString(11,"60 30");
        
        data = Data(4600)
        
        _scale = 1.0
        _scale_theory = 1.0
        
        data.add_with_errors('obs1', 'e_brown_20/limit_test.txt',
                             'S/B ordering, 20%',
                             1.0,
                             ROOT.kBlack, 1, 3,
                             ROOT.kBlack, 0.1, 8,
                             fill_style = 1002,
                             index_x = 0,
                             index_y = 4)
        data.add_with_errors('obs2', 'e_brown_20b/limit_test.txt',
                             'S/B ordering, 20% bg',
                             1.0,
                             ROOT.kGreen+2, 1, 3,
                             ROOT.kGreen+2, 0.1, 8,
                             fill_style = 1002,
                             index_x = 0,
                             index_y = 4)
        data.add_with_errors('obs3', 'e_davis_20b/limit_test.txt',
                             'adjacent bins, 20% bg',
                             1.0,
                             ROOT.kRed, 2, 3,
                             ROOT.kRed, 0.1, 2,
                             index_x = 0,
                             index_y = 4)
        
        data.add('theory', 'ascii/tprime_theory.ascii', 't\'_{THEORY}',
                 1.0,
                 ROOT.kBlue, 1, 3,
                 ROOT.kBlue, 0.1, 8)

        makeMultiGraph(data, 'e_sensitivity_20pc_4700ipb.'+plot_format,
                       400, 600, 0.05, 5.0,
                       xlabel = "M_{t'} [GeV]", ylabel = '#sigma [pb]',
                       xLegend = .55, yLegend = .55,
                       legendWidth = 0.45, legendHeight = 0.30,
                       fillStyle = 1002,
                       drawOption = 'APL3', letter = '',
                       toplegend = '#font[42]{CL_{S}: '+channel+'}')



    #------------------------------------------------------------>
    #7
    # e+jets algo comparison
    # Asymptotic expected CLs limit
    # 40%
    #    

    if str(7) in plot_number:

        channel = 'e+jets, expected 95% C.L. limit'
    
        ROOT.gStyle.SetHatchesLineWidth(2)
        ROOT.gStyle.SetLineStyleString(11,"60 30");
        
        data = Data(4600)
        
        _scale = 1.0
        _scale_theory = 1.0
        
        data.add_with_errors('obs1', 'e_brown_40/limit_test.txt',
                             'S/B ordering, 40%',
                             1.0,
                             ROOT.kBlack, 1, 3,
                             ROOT.kBlack, 0.1, 8,
                             fill_style = 1002,
                             index_x = 0,
                             index_y = 4)
        data.add_with_errors('obs2', 'e_brown_40b/limit_test.txt',
                             'S/B ordering, 40% bg',
                             1.0,
                             ROOT.kGreen+2, 1, 3,
                             ROOT.kGreen+2, 0.1, 8,
                             fill_style = 1002,
                             index_x = 0,
                             index_y = 4)
        data.add_with_errors('obs3', 'e_davis_40b/limit_test.txt',
                             'adjacent bins, 40% bg',
                             1.0,
                             ROOT.kRed, 2, 3,
                             ROOT.kRed, 0.1, 2,
                             index_x = 0,
                             index_y = 4)
        
        data.add('theory', 'ascii/tprime_theory.ascii', 't\'_{THEORY}',
                 1.0,
                 ROOT.kBlue, 1, 3,
                 ROOT.kBlue, 0.1, 8)

        makeMultiGraph(data, 'e_sensitivity_40pc_4700ipb.'+plot_format,
                       400, 600, 0.05, 5.0,
                       xlabel = "M_{t'} [GeV]", ylabel = '#sigma [pb]',
                       xLegend = .55, yLegend = .55,
                       legendWidth = 0.45, legendHeight = 0.30,
                       fillStyle = 1002,
                       drawOption = 'APL3', letter = '',
                       toplegend = '#font[42]{CL_{S}: '+channel+'}')



    #------------------------------------------------------------>
    #8
    # e+jets sensitivity comparison
    # Asymptotic expected CLs limit
    # 40%, 20, 10... for brown
    #    

    if str(8) in plot_number:

        channel = 'e+jets, expected 95% C.L. limit'
    
        ROOT.gStyle.SetHatchesLineWidth(2)
        ROOT.gStyle.SetLineStyleString(11,"60 30");
        
        data = Data(4700)
        
        _scale = 1.0
        _scale_theory = 1.0
        
        data.add_with_errors('obs1', 'e_brown_40/limit_test.txt',
                             'S/B ordering, 40%',
                             1.0,
                             ROOT.kRed, 1, 3,
                             ROOT.kRed, 0.1, 8,
                             fill_style = 1002,
                             index_x = 0,
                             index_y = 4)
        data.add_with_errors('obs2', 'e_brown_40b/limit_test.txt',
                             'S/B ordering, 40% bg',
                             1.0,
                             ROOT.kRed, 2, 3,
                             ROOT.kRed, 0.1, 8,
                             fill_style = 1002,
                             index_x = 0,
                             index_y = 4)
        data.add_with_errors('obs3', 'e_brown_20/limit_test.txt',
                             'S/B ordering, 20%',
                             1.0,
                             ROOT.kBlack, 1, 3,
                             ROOT.kBlack, 0.1, 8,
                             fill_style = 1002,
                             index_x = 0,
                             index_y = 4)
        data.add_with_errors('obs4', 'e_brown_20b/limit_test.txt',
                             'S/B ordering, 20% bg',
                             1.0,
                             ROOT.kBlack, 2, 3,
                             ROOT.kBlack, 0.1, 8,
                             fill_style = 1002,
                             index_x = 0,
                             index_y = 4)
        
        data.add('theory', 'ascii/tprime_theory.ascii', 't\'_{THEORY}',
                 1.0,
                 ROOT.kBlue, 1, 3,
                 ROOT.kBlue, 0.1, 8)

        makeMultiGraph(data, 'e_sensitivity_brown_4700ipb.'+plot_format,
                       400, 600, 0.05, 5.0,
                       xlabel = "M_{t'} [GeV]", ylabel = '#sigma [pb]',
                       xLegend = .55, yLegend = .55,
                       legendWidth = 0.45, legendHeight = 0.30,
                       fillStyle = 1002,
                       drawOption = 'APL3', letter = '',
                       toplegend = '#font[42]{CL_{S}: '+channel+'}')



    #------------------------------------------------------------>
    #9
    # e+jets sensitivity comparison
    # Asymptotic expected CLs limit
    # 40%, 20, 10... for brown
    #    

    if str(9) in plot_number:

        channel = '#mu+jets, expected 95% C.L. limit'
    
        ROOT.gStyle.SetHatchesLineWidth(2)
        ROOT.gStyle.SetLineStyleString(11,"60 30");
        
        data = Data(4600)
        
        _scale = 1.0
        _scale_theory = 1.0
        
        data.add_with_errors('obs1', 'mu_brown_40/limit_test.txt',
                             'S/B ordering, 40%',
                             1.0,
                             ROOT.kRed, 1, 3,
                             ROOT.kRed, 0.1, 8,
                             fill_style = 1002,
                             index_x = 0,
                             index_y = 4)
        data.add_with_errors('obs2', 'mu_brown_40b/limit_test.txt',
                             'S/B ordering, 40% bg',
                             1.0,
                             ROOT.kRed, 2, 3,
                             ROOT.kRed, 0.1, 8,
                             fill_style = 1002,
                             index_x = 0,
                             index_y = 4)
        data.add_with_errors('obs3', 'mu_brown_20/limit_test.txt',
                             'S/B ordering, 20%',
                             1.0,
                             ROOT.kBlack, 1, 3,
                             ROOT.kBlack, 0.1, 8,
                             fill_style = 1002,
                             index_x = 0,
                             index_y = 4)
        data.add_with_errors('obs4', 'mu_brown_20b/limit_test.txt',
                             'S/B ordering, 20% bg',
                             1.0,
                             ROOT.kBlack, 2, 3,
                             ROOT.kBlack, 0.1, 8,
                             fill_style = 1002,
                             index_x = 0,
                             index_y = 4)
        data.add_with_errors('obs5', 'mu_brown_10/limit_test.txt',
                             'S/B ordering, 10%',
                             1.0,
                             ROOT.kGreen+2, 1, 3,
                             ROOT.kGreen+2, 0.1, 8,
                             fill_style = 1002,
                             index_x = 0,
                             index_y = 4)
        data.add_with_errors('obs6', 'mu_brown_10b/limit_test.txt',
                             'S/B ordering, 10% bg',
                             1.0,
                             ROOT.kGreen+2, 2, 3,
                             ROOT.kGreen+2, 0.1, 8,
                             fill_style = 1002,
                             index_x = 0,
                             index_y = 4)
        
        data.add('theory', 'ascii/tprime_theory.ascii', 't\'_{THEORY}',
                 1.0,
                 ROOT.kBlue, 1, 3,
                 ROOT.kBlue, 0.1, 8)

        makeMultiGraph(data, 'mu_sensitivity_brown_4600ipb.'+plot_format,
                       400, 600, 0.05, 5.0,
                       xlabel = "M_{t'} [GeV]", ylabel = '#sigma [pb]',
                       xLegend = .55, yLegend = .55,
                       legendWidth = 0.45, legendHeight = 0.30,
                       fillStyle = 1002,
                       drawOption = 'APL3', letter = '',
                       toplegend = '#font[42]{CL_{S}: '+channel+'}')



    #------------------------------------------------------------>
    #10
    # e+jets
    # Asymptotic CLs limit
    # davis 20b
    #    

    if str(10) in plot_number:

        channel = 'e+jets'
    
        ROOT.gStyle.SetHatchesLineWidth(2)
        ROOT.gStyle.SetLineStyleString(11,"60 30");
        
        data = Data(4700)
        
        _scale = 1.0
        _scale_theory = 1.0
        
        _dir = 'ejets_adj20_01feb2012v1'

        data.add_with_errors('RD_exp1', _dir+'/limits.txt',
                             '#pm2#sigma expected',
                             1.0,
                             ROOT.kYellow, 1, 3,
                             ROOT.kYellow, 0.1, 8,
                             fill_style = 1002,
                             index_x = 0,
                             index_y = 4,
                             index_ey_up = 6,
                             index_ey_down = 2)
        data.add_with_errors('RD_exp2', _dir+'/limits.txt',
                             '#pm1#sigma expected',
                             1.0,
                             ROOT.kGreen, 1, 3,
                             ROOT.kGreen, 0.1, 8,
                             fill_style = 1002,
                             index_x = 0,
                             index_y = 4,
                             index_ey_up = 5,
                             index_ey_down = 3)
        if plot_observed:
            data.add_with_errors('obs', _dir+'/limits.txt',
                                 'observed 95% C.L.',
                                 1.0,
                                 ROOT.kBlack, 1, 3,
                                 ROOT.kBlack, 1.2, 8,
                                 fill_style = 1002,
                                 index_x = 0,
                                 index_y = 1)
        data.add_with_errors('obs2', _dir+'/limits.txt',
                             'expected',
                             1.0,
                             ROOT.kBlack, 2, 3,
                             ROOT.kBlack, 0.1, 2,
                             index_x = 0,
                             index_y = 4)
        
        data.add('theory', 'ascii/tprime_theory.ascii', 't\'_{THEORY}',
                 1.0,
                 ROOT.kBlue, 1, 3,
                 ROOT.kBlue, 0.1, 8)

        makeMultiGraph(data, 'ejets_limit_adj20_4.7ifb_'+suffix+'.'+plot_format,
                       400, 625, 0.05, 5.0,
                       xlabel = "M_{t'} [GeV]", ylabel = '#sigma [pb]',
                       xLegend = .55, yLegend = .55,
                       legendWidth = 0.45, legendHeight = 0.30,
                       fillStyle = 1002,
                       drawOption = 'APL3', letter = '',
                       toplegend = '#font[42]{CL_{S}: '+channel+'}')



    #------------------------------------------------------------>
    #11
    # mu+jets
    # Asymptotic CLs limit
    # davis 20b
    #    

    if str(11) in plot_number:

        channel = '#mu+jets'
    
        ROOT.gStyle.SetHatchesLineWidth(2)
        ROOT.gStyle.SetLineStyleString(11,"60 30");
        
        data = Data(4700)
        
        _scale = 1.0
        _scale_theory = 1.0
        
        _dir = 'mujets_adj20_01feb2012v1'

        data.add_with_errors('RD_exp1', _dir+'/limits.txt',
                             '#pm2#sigma expected',
                             1.0,
                             ROOT.kYellow, 1, 3,
                             ROOT.kYellow, 0.1, 8,
                             fill_style = 1002,
                             index_x = 0,
                             index_y = 4,
                             index_ey_up = 6,
                             index_ey_down = 2)
        data.add_with_errors('RD_exp2', _dir+'/limits.txt',
                             '#pm1#sigma expected',
                             1.0,
                             ROOT.kGreen, 1, 3,
                             ROOT.kGreen, 0.1, 8,
                             fill_style = 1002,
                             index_x = 0,
                             index_y = 4,
                             index_ey_up = 5,
                             index_ey_down = 3)
        if plot_observed:
            data.add_with_errors('obs', _dir+'/limits.txt',
                                 'observed 95% C.L.',
                                 1.0,
                                 ROOT.kBlack, 1, 3,
                                 ROOT.kBlack, 1.2, 8,
                                 fill_style = 1002,
                                 index_x = 0,
                                 index_y = 1)
        data.add_with_errors('obs2', _dir+'/limits.txt',
                             'expected',
                             1.0,
                             ROOT.kBlack, 2, 3,
                             ROOT.kBlack, 0.1, 2,
                             index_x = 0,
                             index_y = 4)
        
        data.add('theory', 'ascii/tprime_theory.ascii', 't\'_{THEORY}',
                 1.0,
                 ROOT.kBlue, 1, 3,
                 ROOT.kBlue, 0.1, 8)

        makeMultiGraph(data, 'mujets_limit_adj20_4.6ifb_'+suffix+'.'+plot_format,
                       400, 625, 0.05, 5.0,
                       xlabel = "M_{t'} [GeV]", ylabel = '#sigma [pb]',
                       xLegend = .55, yLegend = .55,
                       legendWidth = 0.45, legendHeight = 0.30,
                       fillStyle = 1002,
                       drawOption = 'APL3', letter = '',
                       toplegend = '#font[42]{CL_{S}: '+channel+'}')



    #------------------------------------------------------------>
    #1101
    # mu+jets
    # Asymptotic CLs limit
    # davis 20b
    #    

    if str(1101) in plot_number:

        channel = '#mu+jets (4.6fb^{-1}), e+jets (4.7fb^{-1})'
    
        ROOT.gStyle.SetHatchesLineWidth(2)
        ROOT.gStyle.SetLineStyleString(11,"60 30");
        
        data = Data()
        
        _scale = 1.0
        _scale_theory = 1.0
        
        _dir = 'comb_adj20_01feb2012v1'

        data.add_with_errors('RD_exp1', _dir+'/limits.txt',
                             '#pm2#sigma expected',
                             1.0,
                             ROOT.kYellow, 1, 3,
                             ROOT.kYellow, 0.1, 8,
                             fill_style = 1002,
                             index_x = 0,
                             index_y = 4,
                             index_ey_up = 6,
                             index_ey_down = 2)
        data.add_with_errors('RD_exp2', _dir+'/limits.txt',
                             '#pm1#sigma expected',
                             1.0,
                             ROOT.kGreen, 1, 3,
                             ROOT.kGreen, 0.1, 8,
                             fill_style = 1002,
                             index_x = 0,
                             index_y = 4,
                             index_ey_up = 5,
                             index_ey_down = 3)
        if plot_observed:
            data.add_with_errors('obs', _dir+'/limits.txt',
                                 'observed 95% C.L.',
                                 1.0,
                                 ROOT.kBlack, 1, 3,
                                 ROOT.kBlack, 1.2, 8,
                                 fill_style = 1002,
                                 index_x = 0,
                                 index_y = 1)
        data.add_with_errors('obs2', _dir+'/limits.txt',
                             'expected',
                             1.0,
                             ROOT.kBlack, 2, 3,
                             ROOT.kBlack, 0.1, 2,
                             index_x = 0,
                             index_y = 4)
        
        data.add('theory', 'ascii/tprime_theory.ascii', 't\'_{THEORY}',
                 1.0,
                 ROOT.kBlue, 1, 3,
                 ROOT.kBlue, 0.1, 8)

        makeMultiGraph(data, 'comb_limit_adj20_4.6ifb_'+suffix+'.'+plot_format,
                       400, 625, 0.05, 5.0,
                       xlabel = "M_{t'} [GeV]", ylabel = '#sigma [pb]',
                       xLegend = .55, yLegend = .55,
                       legendWidth = 0.45, legendHeight = 0.30,
                       fillStyle = 1002,
                       drawOption = 'APL3', letter = '',
                       toplegend = '#font[42]{CL_{S}: '+channel+'}')



    #------------------------------------------------------------>
    #1102
    # comparison of ejets, mujets and combined
    # Asymptotic CLs limit
    # davis 20 b
    #    

    if str(1102) in plot_number:

        channel = '#mu+jets (4.6fb^{-1}), e+jets (4.7fb^{-1})'
    
        ROOT.gStyle.SetHatchesLineWidth(2)
        ROOT.gStyle.SetLineStyleString(11,"60 30");
        
        data = Data()
        
        _scale = 1.0
        _scale_theory = 1.0
        
        if plot_observed:
            #data.add_with_errors('obs', 'ejets_19jan2012v2/limits.txt',
            data.add_with_errors('obs', 'ejets_adj20_01feb2012v1/limits.txt',
                                 'observed e+jets',
                                 1.0,
                                 ROOT.kBlue, 1, 3,
                                 ROOT.kBlue, 1.2, 8,
                                 fill_style = 1002,
                                 index_x = 0,
                                 index_y = 1)
            #data.add_with_errors('obs2', 'mujets_19jan2012v2/limits.txt',
            data.add_with_errors('obs2', 'mujets_adj20_01feb2012v1/limits.txt',
                                 'observed #mu+jets',
                                 1.0,
                                 ROOT.kRed, 1, 3,
                                 ROOT.kRed, 1.2, 8,
                                 fill_style = 1002,
                                 index_x = 0,
                                 index_y = 1)
            #data.add_with_errors('obs3', 'comb_19jan2012v2/limits.txt',
            data.add_with_errors('obs3', 'comb_adj20_01feb2012v1/limits.txt',
                                 'observed combined',
                                 1.0,
                                 ROOT.kBlack, 1, 3,
                                 ROOT.kBlack, 1.2, 8,
                                 fill_style = 1002,
                                 index_x = 0,
                                 index_y = 1)
        
        data.add_with_errors('obs4', 'ejets_adj20_01feb2012v1/limits.txt',
                             'expected e+jets',
                             1.0,
                             ROOT.kBlue, 2, 3,
                             ROOT.kBlue, 0.1, 2,
                             index_x = 0,
                             index_y = 4)
        data.add_with_errors('obs5', 'mujets_adj20_01feb2012v1/limits.txt',
                             'expected #mu+jets',
                             1.0,
                             ROOT.kRed, 2, 3,
                             ROOT.kRed, 0.1, 2,
                             index_x = 0,
                             index_y = 4)
        data.add_with_errors('obs6', 'comb_adj20_01feb2012v1/limits.txt',
                             'expected combined',
                             1.0,
                             ROOT.kBlack, 2, 3,
                             ROOT.kBlack, 0.1, 2,
                             index_x = 0,
                             index_y = 4)
        
        data.add('theory', 'ascii/tprime_theory.ascii', 't\'_{THEORY}',
                 1.0,
                 ROOT.kBlue, 1, 3,
                 ROOT.kBlue, 0.1, 8)

        makeMultiGraph(data, 'cls_limit_adj20_comp_5ifb.'+plot_format,
                       400, 625, 0.05, 5.0,
                       xlabel = "M_{t'} [GeV]", ylabel = '#sigma [pb]',
                       xLegend = .55, yLegend = .55,
                       legendWidth = 0.45, legendHeight = 0.30,
                       fillStyle = 1002,
                       drawOption = 'APL3', letter = '',
                       toplegend = '#font[42]{CL_{S}: '+channel+'}')




    #------------------------------------------------------------>
    #12
    # e+jets
    # Asymptotic CLs limit
    # davis 40b
    #    

    if str(12) in plot_number:

        channel = 'e+jets'
    
        ROOT.gStyle.SetHatchesLineWidth(2)
        ROOT.gStyle.SetLineStyleString(11,"60 30");
        
        data = Data(4700)
        
        _scale = 1.0
        _scale_theory = 1.0
        
        data.add_with_errors('RD_exp1', 'e_davis_40b/limit_test.txt',
                             '#pm2#sigma expected',
                             1.0,
                             ROOT.kYellow, 1, 3,
                             ROOT.kYellow, 0.1, 8,
                             fill_style = 1002,
                             index_x = 0,
                             index_y = 4,
                             index_ey_up = 6,
                             index_ey_down = 2)
        data.add_with_errors('RD_exp2', 'e_davis_40b/limit_test.txt',
                             '#pm1#sigma expected',
                             1.0,
                             ROOT.kGreen, 1, 3,
                             ROOT.kGreen, 0.1, 8,
                             fill_style = 1002,
                             index_x = 0,
                             index_y = 4,
                             index_ey_up = 5,
                             index_ey_down = 3)
        data.add_with_errors('obs', 'e_davis_40b/limit_test.txt',
                             'observed 95% C.L.',
                             1.0,
                             ROOT.kBlack, 1, 3,
                             ROOT.kBlack, 1.2, 8,
                             fill_style = 1002,
                             index_x = 0,
                             index_y = 1)
        data.add_with_errors('obs2', 'e_davis_40b/limit_test.txt',
                             'expected',
                             1.0,
                             ROOT.kBlack, 2, 3,
                             ROOT.kBlack, 0.1, 2,
                             index_x = 0,
                             index_y = 4)
        
        data.add('theory', 'ascii/tprime_theory.ascii', 't\'_{THEORY}',
                 1.0,
                 ROOT.kBlue, 1, 3,
                 ROOT.kBlue, 0.1, 8)

        makeMultiGraph(data, 'e_limit_davis_40b_4700ipb.'+plot_format,
                       400, 600, 0.05, 5.0,
                       xlabel = "M_{t'} [GeV]", ylabel = '#sigma [pb]',
                       xLegend = .55, yLegend = .55,
                       legendWidth = 0.45, legendHeight = 0.30,
                       fillStyle = 1002,
                       drawOption = 'APL3', letter = '',
                       toplegend = '#font[42]{CL_{S}: '+channel+'}')



    #------------------------------------------------------------>
    #13
    # mu+jets
    # Asymptotic CLs limit
    # davis 40b
    #    

    if str(13) in plot_number:

        channel = '#mu+jets'
    
        ROOT.gStyle.SetHatchesLineWidth(2)
        ROOT.gStyle.SetLineStyleString(11,"60 30");
        
        data = Data(4700)
        
        _scale = 1.0
        _scale_theory = 1.0
        
        data.add_with_errors('RD_exp1', 'mu_davis_40b/limit_test.txt',
                             '#pm2#sigma expected',
                             1.0,
                             ROOT.kYellow, 1, 3,
                             ROOT.kYellow, 0.1, 8,
                             fill_style = 1002,
                             index_x = 0,
                             index_y = 4,
                             index_ey_up = 6,
                             index_ey_down = 2)
        data.add_with_errors('RD_exp2', 'mu_davis_40b/limit_test.txt',
                             '#pm1#sigma expected',
                             1.0,
                             ROOT.kGreen, 1, 3,
                             ROOT.kGreen, 0.1, 8,
                             fill_style = 1002,
                             index_x = 0,
                             index_y = 4,
                             index_ey_up = 5,
                             index_ey_down = 3)
        data.add_with_errors('obs', 'mu_davis_40b/limit_test.txt',
                             'observed 95% C.L.',
                             1.0,
                             ROOT.kBlack, 1, 3,
                             ROOT.kBlack, 1.2, 8,
                             fill_style = 1002,
                             index_x = 0,
                             index_y = 1)
        data.add_with_errors('obs2', 'mu_davis_40b/limit_test.txt',
                             'expected',
                             1.0,
                             ROOT.kBlack, 2, 3,
                             ROOT.kBlack, 0.1, 2,
                             index_x = 0,
                             index_y = 4)
        
        data.add('theory', 'ascii/tprime_theory.ascii', 't\'_{THEORY}',
                 1.0,
                 ROOT.kBlue, 1, 3,
                 ROOT.kBlue, 0.1, 8)

        makeMultiGraph(data, 'mu_limit_davis_40b_4600ipb.'+plot_format,
                       400, 600, 0.05, 5.0,
                       xlabel = "M_{t'} [GeV]", ylabel = '#sigma [pb]',
                       xLegend = .55, yLegend = .55,
                       legendWidth = 0.45, legendHeight = 0.30,
                       fillStyle = 1002,
                       drawOption = 'APL3', letter = '',
                       toplegend = '#font[42]{CL_{S}: '+channel+'}')



    #------------------------------------------------------------>
    #14
    # e+jets
    # Asymptotic CLs limit
    # brown 40 (s and b)
    #    

    if str(14) in plot_number:

        channel = 'e+jets'
    
        ROOT.gStyle.SetHatchesLineWidth(2)
        ROOT.gStyle.SetLineStyleString(11,"60 30");
        
        data = Data(4700)
        
        _scale = 1.0
        _scale_theory = 1.0
        
        data.add_with_errors('RD_exp1', 'e_brown_40/limit_test.txt',
                             '#pm2#sigma expected',
                             1.0,
                             ROOT.kYellow, 1, 3,
                             ROOT.kYellow, 0.1, 8,
                             fill_style = 1002,
                             index_x = 0,
                             index_y = 4,
                             index_ey_up = 6,
                             index_ey_down = 2)
        data.add_with_errors('RD_exp2', 'e_brown_40/limit_test.txt',
                             '#pm1#sigma expected',
                             1.0,
                             ROOT.kGreen, 1, 3,
                             ROOT.kGreen, 0.1, 8,
                             fill_style = 1002,
                             index_x = 0,
                             index_y = 4,
                             index_ey_up = 5,
                             index_ey_down = 3)
        data.add_with_errors('obs', 'e_brown_40/limit_test.txt',
                             'observed 95% C.L.',
                             1.0,
                             ROOT.kBlack, 1, 3,
                             ROOT.kBlack, 1.2, 8,
                             fill_style = 1002,
                             index_x = 0,
                             index_y = 1)
        data.add_with_errors('obs2', 'e_brown_40/limit_test.txt',
                             'expected',
                             1.0,
                             ROOT.kBlack, 2, 3,
                             ROOT.kBlack, 0.1, 2,
                             index_x = 0,
                             index_y = 4)
        
        data.add('theory', 'ascii/tprime_theory.ascii', 't\'_{THEORY}',
                 1.0,
                 ROOT.kBlue, 1, 3,
                 ROOT.kBlue, 0.1, 8)

        makeMultiGraph(data, 'e_limit_brown_40_4700ipb.'+plot_format,
                       400, 600, 0.05, 5.0,
                       xlabel = "M_{t'} [GeV]", ylabel = '#sigma [pb]',
                       xLegend = .55, yLegend = .55,
                       legendWidth = 0.45, legendHeight = 0.30,
                       fillStyle = 1002,
                       drawOption = 'APL3', letter = '',
                       toplegend = '#font[42]{CL_{S}: '+channel+'}')



    #------------------------------------------------------------>
    #15
    # mu+jets
    # Asymptotic CLs limit
    # brown 40 (s and b)
    #    

    if str(15) in plot_number:

        channel = '#mu+jets'
    
        ROOT.gStyle.SetHatchesLineWidth(2)
        ROOT.gStyle.SetLineStyleString(11,"60 30");
        
        data = Data(4700)
        
        _scale = 1.0
        _scale_theory = 1.0
        
        data.add_with_errors('RD_exp1', 'mu_brown_40/limit_test.txt',
                             '#pm2#sigma expected',
                             1.0,
                             ROOT.kYellow, 1, 3,
                             ROOT.kYellow, 0.1, 8,
                             fill_style = 1002,
                             index_x = 0,
                             index_y = 4,
                             index_ey_up = 6,
                             index_ey_down = 2)
        data.add_with_errors('RD_exp2', 'mu_brown_40/limit_test.txt',
                             '#pm1#sigma expected',
                             1.0,
                             ROOT.kGreen, 1, 3,
                             ROOT.kGreen, 0.1, 8,
                             fill_style = 1002,
                             index_x = 0,
                             index_y = 4,
                             index_ey_up = 5,
                             index_ey_down = 3)
        data.add_with_errors('obs', 'mu_brown_40/limit_test.txt',
                             'observed 95% C.L.',
                             1.0,
                             ROOT.kBlack, 1, 3,
                             ROOT.kBlack, 1.2, 8,
                             fill_style = 1002,
                             index_x = 0,
                             index_y = 1)
        data.add_with_errors('obs2', 'mu_brown_40/limit_test.txt',
                             'expected',
                             1.0,
                             ROOT.kBlack, 2, 3,
                             ROOT.kBlack, 0.1, 2,
                             index_x = 0,
                             index_y = 4)
        
        data.add('theory', 'ascii/tprime_theory.ascii', 't\'_{THEORY}',
                 1.0,
                 ROOT.kBlue, 1, 3,
                 ROOT.kBlue, 0.1, 8)

        makeMultiGraph(data, 'mu_limit_brown_40_4600ipb.'+plot_format,
                       400, 600, 0.05, 5.0,
                       xlabel = "M_{t'} [GeV]", ylabel = '#sigma [pb]',
                       xLegend = .55, yLegend = .55,
                       legendWidth = 0.45, legendHeight = 0.30,
                       fillStyle = 1002,
                       drawOption = 'APL3', letter = '',
                       toplegend = '#font[42]{CL_{S}: '+channel+'}')



    #------------------------------------------------------------>
    #16
    # mu+jets
    # Asymptotic CLs limit
    # brown 10 (s and b)
    #    

    if str(16) in plot_number:

        channel = '#mu+jets'
    
        ROOT.gStyle.SetHatchesLineWidth(2)
        ROOT.gStyle.SetLineStyleString(11,"60 30");
        
        data = Data(4700)
        
        _scale = 1.0
        _scale_theory = 1.0
        
        data.add_with_errors('RD_exp1', 'mu_brown_10/limit_test.txt',
                             '#pm2#sigma expected',
                             1.0,
                             ROOT.kYellow, 1, 3,
                             ROOT.kYellow, 0.1, 8,
                             fill_style = 1002,
                             index_x = 0,
                             index_y = 4,
                             index_ey_up = 6,
                             index_ey_down = 2)
        data.add_with_errors('RD_exp2', 'mu_brown_10/limit_test.txt',
                             '#pm1#sigma expected',
                             1.0,
                             ROOT.kGreen, 1, 3,
                             ROOT.kGreen, 0.1, 8,
                             fill_style = 1002,
                             index_x = 0,
                             index_y = 4,
                             index_ey_up = 5,
                             index_ey_down = 3)
        data.add_with_errors('obs', 'mu_brown_10/limit_test.txt',
                             'observed 95% C.L.',
                             1.0,
                             ROOT.kBlack, 1, 3,
                             ROOT.kBlack, 1.2, 8,
                             fill_style = 1002,
                             index_x = 0,
                             index_y = 1)
        data.add_with_errors('obs2', 'mu_brown_10/limit_test.txt',
                             'expected',
                             1.0,
                             ROOT.kBlack, 2, 3,
                             ROOT.kBlack, 0.1, 2,
                             index_x = 0,
                             index_y = 4)
        
        data.add('theory', 'ascii/tprime_theory.ascii', 't\'_{THEORY}',
                 1.0,
                 ROOT.kBlue, 1, 3,
                 ROOT.kBlue, 0.1, 8)

        makeMultiGraph(data, 'mu_limit_brown_10_4600ipb.'+plot_format,
                       400, 600, 0.05, 5.0,
                       xlabel = "M_{t'} [GeV]", ylabel = '#sigma [pb]',
                       xLegend = .55, yLegend = .55,
                       legendWidth = 0.45, legendHeight = 0.30,
                       fillStyle = 1002,
                       drawOption = 'APL3', letter = '',
                       toplegend = '#font[42]{CL_{S}: '+channel+'}')



    #------------------------------------------------------------>
    #1000
    # e+jets
    # Asymptotic CLs
    # toys: split vs single
    # brown split s20/b20
    #    

    if str(1000) in plot_number:
        
        channel = 'e+jets'
        
        ROOT.gStyle.SetHatchesLineWidth(2)
        ROOT.gStyle.SetLineStyleString(11,"60 30");
        
        data = Data(3560)
        
        _scale = 1.0
        _scale_theory = 1.0
        
        data.add_with_errors('RD_exp1', 'toy_limit_first_s20b20.ascii',
                             '#pm2#sigma expected',
                             1.0,
                             ROOT.kYellow, 1, 3,
                             ROOT.kYellow, 1.2, 8,
                             fill_style = 1002,
                             index_x = 0,
                             index_y = 1,
                             index_ey_up = 5,
                             index_ey_down = 4)
        data.add_with_errors('RD_exp2', 'toy_limit_first_s20b20.ascii',
                             '#pm1#sigma expected',
                             1.0,
                             ROOT.kGreen, 1, 3,
                             ROOT.kGreen, 1.2, 8,
                             fill_style = 1002,
                             index_x = 0,
                             index_y = 1,
                             index_ey_up = 3,
                             index_ey_down = 2)
        data.add_with_errors('obs', 'toy_limit_first_s20b20.ascii',
                                    'unsplit',
                             1.0,
                             ROOT.kBlack, 1, 3,
                             ROOT.kBlack, 1.2, 8,
                             fill_style = 1002,
                             index_x = 0,
                             index_y = 1)
        data.add_with_errors('obs_exp_median', 'toy_limit_split_s20b20.ascii',
                             'split',
                             1.0,
                             ROOT.kRed, 1, 3,
                             ROOT.kRed, 1.2, 8,
                             fill_style = 1002,
                             index_x = 0,
                             index_y = 1)
        #    data.add('theory', 'ascii/tprime_theory.ascii', 't\'_{THEORY}',
        #                1.0,
        #                ROOT.kBlue, 1, 3,
        #                ROOT.kBlue, 0.1, 8)
        
        makeMultiGraph(data, 'limit_plot.'+plot_format,
                       400, 600, 0.05, 5.0,
                       xlabel = "M_{t'} [GeV]", ylabel = '#sigma [pb]',
                       xLegend = .55, yLegend = .55,
                       legendWidth = 0.45, legendHeight = 0.30,
                       fillStyle = 1002,
                       drawOption = 'APL3', letter = '',
                       toplegend = '#font[42]{CL_{S}: '+channel+'}')



    #------------------------------------------------------------>
    #1001
    # e+jets
    # Asymptotic CLs
    # toys: split vs single
    # brown split s10/b10
    #    

    if str(1001) in plot_number:

        channel = 'e+jets'
    
        ROOT.gStyle.SetHatchesLineWidth(2)
        ROOT.gStyle.SetLineStyleString(11,"60 30");
        
        data = Data(3560)
        
        _scale = 1.0
        _scale_theory = 1.0
        
        data.add_with_errors('RD_exp1', 'toy_limit_first_s10b10.ascii',
                             '#pm2#sigma expected',
                             1.0,
                             ROOT.kYellow, 1, 3,
                             ROOT.kYellow, 1.2, 8,
                             fill_style = 1002,
                             index_x = 0,
                             index_y = 1,
                             index_ey_up = 5,
                             index_ey_down = 4)
        data.add_with_errors('RD_exp2', 'toy_limit_first_s10b10.ascii',
                             '#pm1#sigma expected',
                             1.0,
                             ROOT.kGreen, 1, 3,
                             ROOT.kGreen, 1.2, 8,
                             fill_style = 1002,
                             index_x = 0,
                             index_y = 1,
                             index_ey_up = 3,
                             index_ey_down = 2)
        data.add_with_errors('obs', 'toy_limit_first_s10b10.ascii',
                             'unsplit',
                             1.0,
                             ROOT.kBlack, 1, 3,
                             ROOT.kBlack, 1.2, 8,
                             fill_style = 1002,
                             index_x = 0,
                             index_y = 1)
        data.add_with_errors('obs_exp_median', 'toy_limit_split_s10b10.ascii',
                             'split',
                             1.0,
                             ROOT.kRed, 1, 3,
                             ROOT.kRed, 1.2, 8,
                             fill_style = 1002,
                             index_x = 0,
                             index_y = 1)
        
        makeMultiGraph(data, 'toy_limit_split_s10b10.'+plot_format,
                       400, 600, 0.05, 5.0,
                       xlabel = "M_{t'} [GeV]", ylabel = '#sigma [pb]',
                       xLegend = .55, yLegend = .55,
                       legendWidth = 0.45, legendHeight = 0.30,
                       fillStyle = 1002,
                       drawOption = 'APL3', letter = '',
                       toplegend = '#font[42]{CL_{S}: '+channel+'}')
        

    #------------------------------------------------------------>
    #18
    # e+jets
    # Bayesian with Theta
    # brown 20 (s and b)
    #    

    if str(18) in plot_number:

        channel = '#mu+jets (4.6fb^{-1})'

        ROOT.gStyle.SetHatchesLineWidth(2)
        ROOT.gStyle.SetLineStyleString(11,"60 30");
        
        data = Data(4700)
        
        _scale = 1.0
        _scale_theory = 1.0

        _dir = 'theta/mujets_20jan2012v1'
        
        data.add_with_errors('RD_exp1', _dir+'/limits.txt',
                             '#pm2#sigma expected',
                             1.0,
                             ROOT.kYellow, 1, 3,
                             ROOT.kYellow, 0.1, 8,
                             fill_style = 1002,
                             index_x = 0,
                             index_y = 2,
                             index_ey_up = 4,
                             index_ey_down = 3)
        data.add_with_errors('RD_exp2', _dir+'/limits.txt',
                             '#pm1#sigma expected',
                             1.0,
                             ROOT.kGreen, 1, 0,
                             ROOT.kGreen, 0.1, 8,
                             fill_style = 1002,
                             index_x = 0,
                             index_y = 2,
                             index_ey_up = 6,
                             index_ey_down = 5)
        if plot_observed:
            data.add_with_errors('obs', _dir+'/limits.txt',
                                 'observed 95% C.L.',
                                 1.0,
                                 ROOT.kBlack, 1, 3,
                                 ROOT.kBlack, 1.2, 8,
                                 fill_style = 1002,
                                 index_x = 0,
                                 index_y = 1)
            
        data.add_with_errors('obs2', _dir+'/limits.txt',
                             'expected',
                             1.0,
                             ROOT.kBlack, 2, 3,
                             ROOT.kBlack, 0.1, 2,
                             index_x = 0,
                             index_y = 2)
        
        data.add('theory', 'ascii/tprime_theory.ascii', 't\'_{THEORY}',
                 1.0,
                 ROOT.kBlue, 1, 3,
                 ROOT.kBlue, 0.1, 8)

        makeMultiGraph(data, 'mujets_bayesian_5ifb_'+suffix+'.'+plot_format,
                       400, 625, 0.05, 5.0,
                       xlabel = "M_{t'} [GeV]", ylabel = '#sigma [pb]',
                       xLegend = .55, yLegend = .55,
                       legendWidth = 0.45, legendHeight = 0.30,
                       fillStyle = 1002,
                       drawOption = 'APL3', letter = '',
                       toplegend = '#font[42]{Bayesian: '+channel+'}')




    #------------------------------------------------------------>
    #19
    # e+jets
    # Bayesian with Theta
    # brown 20 (s and b)
    #    

    if str(19) in plot_number:

        channel = 'e+jets (4.7fb^{-1})'

        ROOT.gStyle.SetHatchesLineWidth(2)
        ROOT.gStyle.SetLineStyleString(11,"60 30");
        
        data = Data(4700)
        
        _scale = 1.0
        _scale_theory = 1.0

        _dir = 'theta/ejets_20jan2012v1'
        
        data.add_with_errors('RD_exp1', _dir+'/limits.txt',
                             '#pm2#sigma expected',
                             1.0,
                             ROOT.kYellow, 1, 3,
                             ROOT.kYellow, 0.1, 8,
                             fill_style = 1002,
                             index_x = 0,
                             index_y = 2,
                             index_ey_up = 4,
                             index_ey_down = 3)
        data.add_with_errors('RD_exp2', _dir+'/limits.txt',
                             '#pm1#sigma expected',
                             1.0,
                             ROOT.kGreen, 1, 0,
                             ROOT.kGreen, 0.1, 8,
                             fill_style = 1002,
                             index_x = 0,
                             index_y = 2,
                             index_ey_up = 6,
                             index_ey_down = 5)
        if plot_observed:
            data.add_with_errors('obs', _dir+'/limits.txt',
                                 'observed 95% C.L.',
                                 1.0,
                                 ROOT.kBlack, 1, 3,
                                 ROOT.kBlack, 1.2, 8,
                                 fill_style = 1002,
                                 index_x = 0,
                                 index_y = 1)
            
        data.add_with_errors('obs2', _dir+'/limits.txt',
                             'expected',
                             1.0,
                             ROOT.kBlack, 2, 3,
                             ROOT.kBlack, 0.1, 2,
                             index_x = 0,
                             index_y = 2)
        
        data.add('theory', 'ascii/tprime_theory.ascii', 't\'_{THEORY}',
                 1.0,
                 ROOT.kBlue, 1, 3,
                 ROOT.kBlue, 0.1, 8)

        makeMultiGraph(data, 'ejets_bayesian_5ifb_'+suffix+'.'+plot_format,
                       400, 625, 0.05, 5.0,
                       xlabel = "M_{t'} [GeV]", ylabel = '#sigma [pb]',
                       xLegend = .55, yLegend = .55,
                       legendWidth = 0.45, legendHeight = 0.30,
                       fillStyle = 1002,
                       drawOption = 'APL3', letter = '',
                       toplegend = '#font[42]{Bayesian: '+channel+'}')




    #------------------------------------------------------------>
    #20
    # e+jets + mu+jets
    # Bayesian with Theta
    # brown 20 (s and b)
    #    

    if str(20) in plot_number:

        channel = '#mu+jets (4.6fb^{-1}), e+jets (4.7fb^{-1})'

        ROOT.gStyle.SetHatchesLineWidth(2)
        ROOT.gStyle.SetLineStyleString(11,"60 30");
        
        data = Data(4700)
        
        _scale = 1.0
        _scale_theory = 1.0

        _dir = 'theta/comb_20jan2012v1'
        
        data.add_with_errors('RD_exp1', _dir+'/limits.txt',
                             '#pm2#sigma expected',
                             1.0,
                             ROOT.kYellow, 1, 3,
                             ROOT.kYellow, 0.1, 8,
                             fill_style = 1002,
                             index_x = 0,
                             index_y = 2,
                             index_ey_up = 4,
                             index_ey_down = 3)
        data.add_with_errors('RD_exp2', _dir+'/limits.txt',
                             '#pm1#sigma expected',
                             1.0,
                             ROOT.kGreen, 1, 0,
                             ROOT.kGreen, 0.1, 8,
                             fill_style = 1002,
                             index_x = 0,
                             index_y = 2,
                             index_ey_up = 6,
                             index_ey_down = 5)
        if plot_observed:
            data.add_with_errors('obs', _dir+'/limits.txt',
                                 'observed 95% C.L.',
                                 1.0,
                                 ROOT.kBlack, 1, 3,
                                 ROOT.kBlack, 1.2, 8,
                                 fill_style = 1002,
                                 index_x = 0,
                                 index_y = 1)
            
        data.add_with_errors('obs2', _dir+'/limits.txt',
                             'expected',
                             1.0,
                             ROOT.kBlack, 2, 3,
                             ROOT.kBlack, 0.1, 2,
                             index_x = 0,
                             index_y = 2)
        
        data.add('theory', 'ascii/tprime_theory.ascii', 't\'_{THEORY}',
                 1.0,
                 ROOT.kBlue, 1, 3,
                 ROOT.kBlue, 0.1, 8)

        makeMultiGraph(data, 'comb_bayesian_5ifb_'+suffix+'.'+plot_format,
                       400, 625, 0.05, 5.0,
                       xlabel = "M_{t'} [GeV]", ylabel = '#sigma [pb]',
                       xLegend = .55, yLegend = .55,
                       legendWidth = 0.45, legendHeight = 0.30,
                       fillStyle = 1002,
                       drawOption = 'APL3', letter = '',
                       toplegend = '#font[42]{Bayesian: '+channel+'}')




    #------------------------------------------------------------>
    #21
    # e+jets + mu+jets
    # Bayesian (theta) vs CLs comparison
    # brown 20 (s and b)
    #    

    if str(21) in plot_number:

        channel = '#mu+jets (4.6fb^{-1})'

        ROOT.gStyle.SetHatchesLineWidth(2)
        ROOT.gStyle.SetLineStyleString(11,"60 30");
        
        data = Data(4700)
        
        _scale = 1.0
        _scale_theory = 1.0

        _dir1 = 'theta/mujets_20jan2012v1'
        _dir2 = 'mujets_19jan2012v2'

        if plot_observed:
            data.add_with_errors('obs', _dir1+'/limits.txt',
                                 'observed Bayesian',
                                 1.0,
                                 ROOT.kRed, 1, 3,
                                 ROOT.kRed, 1.2, 8,
                                 fill_style = 1002,
                                 index_x = 0,
                                 index_y = 1)
            
            data.add_with_errors('obs2', _dir2+'/limits.txt',
                                 'observed CLs',
                                 1.0,
                                 ROOT.kBlack, 1, 3,
                                 ROOT.kBlack, 1.2, 8,
                                 fill_style = 1002,
                                 index_x = 0,
                                 index_y = 1)
            
        data.add_with_errors('obs3', _dir1+'/limits.txt',
                             'expected Bayesian',
                             1.0,
                             ROOT.kRed, 2, 3,
                             ROOT.kRed, 0.1, 2,
                             index_x = 0,
                             index_y = 2)
        
        data.add_with_errors('obs4', _dir2+'/limits.txt',
                             'expected CLs',
                             1.0,
                             ROOT.kBlack, 2, 3,
                             ROOT.kBlack, 0.1, 2,
                             index_x = 0,
                             index_y = 4)
        
        data.add('theory', 'ascii/tprime_theory.ascii', 't\'_{THEORY}',
                 1.0,
                 ROOT.kBlue, 1, 3,
                 ROOT.kBlue, 0.1, 8)

        makeMultiGraph(data, 'mujets_cls_bayes_5ifb_'+suffix+'.'+plot_format,
                       400, 625, 0.05, 5.0,
                       xlabel = "M_{t'} [GeV]", ylabel = '#sigma [pb]',
                       xLegend = .55, yLegend = .55,
                       legendWidth = 0.45, legendHeight = 0.30,
                       fillStyle = 1002,
                       drawOption = 'APL3', letter = '',
                       toplegend = '#font[42]{'+channel+'}')





    #------------------------------------------------------------>
    #22
    # e+jets
    # Bayesian (theta) vs CLs comparison
    # brown 20 (s and b)
    #    

    if str(22) in plot_number:

        channel = 'e+jets (4.7fb^{-1})'

        ROOT.gStyle.SetHatchesLineWidth(2)
        ROOT.gStyle.SetLineStyleString(11,"60 30");
        
        data = Data(4700)
        
        _scale = 1.0
        _scale_theory = 1.0

        _dir1 = 'theta/ejets_20jan2012v1'
        _dir2 = 'ejets_19jan2012v2'

        if plot_observed:
            data.add_with_errors('obs', _dir1+'/limits.txt',
                                 'observed Bayesian',
                                 1.0,
                                 ROOT.kRed, 1, 3,
                                 ROOT.kRed, 1.2, 8,
                                 fill_style = 1002,
                                 index_x = 0,
                                 index_y = 1)
            
            data.add_with_errors('obs2', _dir2+'/limits.txt',
                                 'observed CLs',
                                 1.0,
                                 ROOT.kBlack, 1, 3,
                                 ROOT.kBlack, 1.2, 8,
                                 fill_style = 1002,
                                 index_x = 0,
                                 index_y = 1)
            
        data.add_with_errors('obs3', _dir1+'/limits.txt',
                             'expected Bayesian',
                             1.0,
                             ROOT.kRed, 2, 3,
                             ROOT.kRed, 0.1, 2,
                             index_x = 0,
                             index_y = 2)
        
        data.add_with_errors('obs4', _dir2+'/limits.txt',
                             'expected CLs',
                             1.0,
                             ROOT.kBlack, 2, 3,
                             ROOT.kBlack, 0.1, 2,
                             index_x = 0,
                             index_y = 4)
        
        data.add('theory', 'ascii/tprime_theory.ascii', 't\'_{THEORY}',
                 1.0,
                 ROOT.kBlue, 1, 3,
                 ROOT.kBlue, 0.1, 8)

        makeMultiGraph(data, 'ejets_cls_bayes_5ifb_'+suffix+'.'+plot_format,
                       400, 625, 0.05, 5.0,
                       xlabel = "M_{t'} [GeV]", ylabel = '#sigma [pb]',
                       xLegend = .55, yLegend = .55,
                       legendWidth = 0.45, legendHeight = 0.30,
                       fillStyle = 1002,
                       drawOption = 'APL3', letter = '',
                       toplegend = '#font[42]{'+channel+'}')





    #------------------------------------------------------------>
    #23
    # e+jets + mu+jets
    # Bayesian (theta) vs CLs comparison
    # brown 20 (s and b)
    #    

    if str(23) in plot_number:

        channel = '#mu+jets (4.6fb^{-1}), e+jets (4.7fb^{-1})'

        ROOT.gStyle.SetHatchesLineWidth(2)
        ROOT.gStyle.SetLineStyleString(11,"60 30");
        
        data = Data(4700)
        
        _scale = 1.0
        _scale_theory = 1.0

        _dir1 = 'theta/comb_20jan2012v1'
        _dir2 = 'comb_19jan2012v2'

        if plot_observed:
            data.add_with_errors('obs', _dir1+'/limits.txt',
                                 'observed Bayesian',
                                 1.0,
                                 ROOT.kRed, 1, 3,
                                 ROOT.kRed, 1.2, 8,
                                 fill_style = 1002,
                                 index_x = 0,
                                 index_y = 1)
            
            data.add_with_errors('obs2', _dir2+'/limits.txt',
                                 'observed CLs',
                                 1.0,
                                 ROOT.kBlack, 1, 3,
                                 ROOT.kBlack, 1.2, 8,
                                 fill_style = 1002,
                                 index_x = 0,
                                 index_y = 1)
            
        data.add_with_errors('obs3', _dir1+'/limits.txt',
                             'expected Bayesian',
                             1.0,
                             ROOT.kRed, 2, 3,
                             ROOT.kRed, 0.1, 2,
                             index_x = 0,
                             index_y = 2)
        
        data.add_with_errors('obs4', _dir2+'/limits.txt',
                             'expected CLs',
                             1.0,
                             ROOT.kBlack, 2, 3,
                             ROOT.kBlack, 0.1, 2,
                             index_x = 0,
                             index_y = 4)
        
        data.add('theory', 'ascii/tprime_theory.ascii', 't\'_{THEORY}',
                 1.0,
                 ROOT.kBlue, 1, 3,
                 ROOT.kBlue, 0.1, 8)

        makeMultiGraph(data, 'comb_cls_bayes_5ifb_'+suffix+'.'+plot_format,
                       400, 625, 0.05, 5.0,
                       xlabel = "M_{t'} [GeV]", ylabel = '#sigma [pb]',
                       xLegend = .55, yLegend = .55,
                       legendWidth = 0.45, legendHeight = 0.30,
                       fillStyle = 1002,
                       drawOption = 'APL3', letter = '',
                       toplegend = '#font[42]{'+channel+'}')





    #------------------------------------------------------------>
    #24
    # e+jets
    # number of bins vs precision
    #    

    if str(24) in plot_number:

        channel = 'e+jets, 550 GeV'

        ROOT.gStyle.SetHatchesLineWidth(2)
        ROOT.gStyle.SetLineStyleString(11,"60 30");
        
        data = Data()
        
        _scale = 1.0
        _scale_theory = 1.0

        _dir = './'

        data.add_with_errors('obs', _dir+'/ejets_nbins_550gev.txt',
                             'N bins',
                             1.0,
                             ROOT.kBlack, 1, 3,
                             ROOT.kBlack, 1.2, 8,
                             fill_style = 1002,
                             index_x = 0,
                             index_y = 1)

        makeMultiGraph(data, 'ejets_comb_bins_'+suffix+'.'+plot_format,
                       0, 60, 0.0, 300.0,
                       xlabel = "bin uncertainty [percent]", ylabel = 'number of bins',
                       xLegend = .55, yLegend = .55,
                       legendWidth = 0.45, legendHeight = 0.30,
                       fillStyle = 1002,
                       drawOption = 'APL3', letter = '',
                       toplegend = '#font[42]{'+channel+'}')





    #------------------------------------------------------------>
    #25
    # e+jets shape systematics study
    # relative deviation from nominal (JES) vs fit mass
    # observed
    #    

    if str(25) in plot_number:

        channel = 'e+jets, 4.7 fb^{-1}'

        ROOT.gStyle.SetHatchesLineWidth(2)
        ROOT.gStyle.SetLineStyleString(11,"60 30");
        
        data = Data()
        
        _scale = 1.0

        _dir = 'e_syst/'

        data.add_function('obs', _dir+'/limits.txt',
                          'b-tag SF',
                          _scale,
                          ROOT.kBlack, 1, 3,
                          ROOT.kBlack, 1.2, 8,
                          fill_style = 1002,
                          index_x = 0,
                          index_y = 1,
                          entry = 1)

        data.add_function('obs2', _dir+'/limits.txt',
                          'JER',
                          _scale,
                          ROOT.kBlue, 1, 3,
                          ROOT.kBlue, 1.2, 8,
                          fill_style = 1002,
                          index_x = 0,
                          index_y = 1,
                          entry = 2)

        data.add_function('obs3', _dir+'/limits.txt',
                          'matching',
                          _scale,
                          ROOT.kGreen, 1, 3,
                          ROOT.kGreen, 1.2, 8,
                          fill_style = 1002,
                          index_x = 0,
                          index_y = 1,
                          entry = 3)

        data.add_function('obs4', _dir+'/limits.txt',
                          'scale',
                          _scale,
                          ROOT.kRed, 1, 3,
                          ROOT.kRed, 1.2, 8,
                          fill_style = 1002,
                          index_x = 0,
                          index_y = 1,
                          entry = 4)


        makeMultiGraph(data, 'ejets_syst_obs_'+suffix+'.'+plot_format,
                       400, 625, -1.0, 1.0,
                       xlabel = "m_{fit}", ylabel = 'fraction of expected sigma',
                       xLegend = .55, yLegend = .55,
                       legendWidth = 0.45, legendHeight = 0.30,
                       fillStyle = 1002,
                       drawOption = 'APL3', letter = '',
                       toplegend = '#font[42]{'+channel+'}')





    #------------------------------------------------------------>
    #26
    # e+jets shape systematics study
    # relative deviation from nominal (JES) vs fit mass
    # expected
    #    

    if str(26) in plot_number:

        channel = 'e+jets, 4.7 fb^{-1}'

        ROOT.gStyle.SetHatchesLineWidth(2)
        ROOT.gStyle.SetLineStyleString(11,"60 30");
        
        data = Data()
        
        _scale = 1.0

        _dir = 'e_syst/'

        data.add_function('obs', _dir+'/limits.txt',
                          'b-tag SF',
                          _scale,
                          ROOT.kBlack, 1, 3,
                          ROOT.kBlack, 1.2, 8,
                          fill_style = 1002,
                          index_x = 0,
                          index_y = 4,
                          entry = 1)

        data.add_function('obs2', _dir+'/limits.txt',
                          'JER',
                          _scale,
                          ROOT.kBlue, 1, 3,
                          ROOT.kBlue, 1.2, 8,
                          fill_style = 1002,
                          index_x = 0,
                          index_y = 4,
                          entry = 2)

        data.add_function('obs3', _dir+'/limits.txt',
                          'matching',
                          _scale,
                          ROOT.kGreen, 1, 3,
                          ROOT.kGreen, 1.2, 8,
                          fill_style = 1002,
                          index_x = 0,
                          index_y = 4,
                          entry = 3)

        data.add_function('obs4', _dir+'/limits.txt',
                          'scale',
                          _scale,
                          ROOT.kRed, 1, 3,
                          ROOT.kRed, 1.2, 8,
                          fill_style = 1002,
                          index_x = 0,
                          index_y = 4,
                          entry = 4)


        makeMultiGraph(data, 'ejets_syst_exp_'+suffix+'.'+plot_format,
                       400, 625, -0.4, 0.6,
                       xlabel = "m_{fit}", ylabel = 'fraction of expected sigma',
                       xLegend = .55, yLegend = .55,
                       legendWidth = 0.45, legendHeight = 0.30,
                       fillStyle = 1002,
                       drawOption = 'APL3', letter = '',
                       toplegend = '#font[42]{'+channel+'}')





    #------------------------------------------------------------>
    #27
    # e+jets shape systematics study
    # relative deviation from nominal (JES+match) vs fit mass
    # observed
    #    

    if str(27) in plot_number:

        channel = 'e+jets, 4.7 fb^{-1}'

        ROOT.gStyle.SetHatchesLineWidth(2)
        ROOT.gStyle.SetLineStyleString(11,"60 30");
        
        data = Data()
        
        _scale = 1.0

        _dir = 'e_syst/'

        data.add_function('obs2', _dir+'/limits.txt',
                          'JER',
                          _scale,
                          ROOT.kBlue, 1, 3,
                          ROOT.kBlue, 1.2, 8,
                          fill_style = 1002,
                          index_x = 0,
                          index_y = 1,
                          entry = 1)

        makeMultiGraph(data, 'ejets_syst_obs_jer'+suffix+'.'+plot_format,
                       400, 625, -1.0, 1.0,
                       xlabel = "m_{fit}", ylabel = 'fraction of expected sigma',
                       xLegend = .55, yLegend = .55,
                       legendWidth = 0.45, legendHeight = 0.30,
                       fillStyle = 1002,
                       drawOption = 'APL3', letter = '',
                       toplegend = '#font[42]{'+channel+'}')





    #------------------------------------------------------------>
    #28
    # mu+jets powheg systematics study
    # relative deviation from nominal (JES+match) vs fit mass
    # observed
    #    

    if str(28) in plot_number:

        channel = '#mu+jets, 4.6 fb^{-1}'

        ROOT.gStyle.SetHatchesLineWidth(2)
        ROOT.gStyle.SetLineStyleString(11,"60 30");
        
        data = Data()
        
        _scale = 1.0

        _dir = 'mu_syst/'

        data.add_function('obs1', _dir+'/limits_powheg.txt',
                          'POWHEG',
                          _scale,
                          ROOT.kBlack, 1, 3,
                          ROOT.kBlack, 1.2, 8,
                          fill_style = 1002,
                          index_x = 0,
                          index_y = 1,
                          entry = 2)

        data.add_function('obs2', _dir+'/limits_powheg.txt',
                          'POWHEG nominal',
                          _scale,
                          ROOT.kBlue, 2, 3,
                          ROOT.kBlue, 0.1, 8,
                          fill_style = 1002,
                          index_x = 0,
                          index_y = 1,
                          entry = 3)

        data.add_function('obs3', _dir+'/limits_powheg.txt',
                          'POWHEG nom+var',
                          _scale,
                          ROOT.kGreen+2, 2, 3,
                          ROOT.kGreen+2, 0.1, 8,
                          fill_style = 1002,
                          index_x = 0,
                          index_y = 1,
                          entry = 1)

        makeMultiGraph(data, 'mujets_syst_powheg_obs'+suffix+'.'+plot_format,
                       400, 625, -1.5, 1.5,
                       xlabel = "m_{fit}", ylabel = 'fraction of expected sigma',
                       xLegend = .55, yLegend = .55,
                       legendWidth = 0.45, legendHeight = 0.30,
                       fillStyle = 1002,
                       drawOption = 'APL3', letter = '',
                       toplegend = '#font[42]{'+channel+'}')





    #------------------------------------------------------------>
    #29
    # mu+jets powheg systematics study
    # relative deviation from nominal (JES+match) vs fit mass
    # new nominal: average of old nominal and powheg
    # observed
    #    

    if str(29) in plot_number:

        channel = '#mu+jets, 4.6 fb^{-1}'

        ROOT.gStyle.SetHatchesLineWidth(2)
        ROOT.gStyle.SetLineStyleString(11,"60 30");
        
        data = Data()
        
        _scale = 1.0

        _dir = 'mu_syst/'

        data.add_function('obs1', _dir+'/limits_powheg_newnom.txt',
                          'POWHEG',
                          _scale,
                          ROOT.kBlack, 1, 3,
                          ROOT.kBlack, 1.2, 8,
                          fill_style = 1002,
                          index_x = 0,
                          index_y = 1,
                          entry = 1)


        makeMultiGraph(data, 'mujets_syst_powheg_newnom_obs'+suffix+'.'+plot_format,
                       400, 625, -1.0, 1.0,
                       xlabel = "m_{fit}", ylabel = 'fraction of expected sigma',
                       xLegend = .55, yLegend = .55,
                       legendWidth = 0.45, legendHeight = 0.30,
                       fillStyle = 1002,
                       drawOption = 'APL3', letter = '',
                       toplegend = '#font[42]{'+channel+'}')




    #------------------------------------------------------------>
    #30
    # mu+jets shape systematics study
    # relative deviation from nominal (JES) vs fit mass
    # observed
    #    

    if str(30) in plot_number:

        channel = '#mu+jets, 4.6 fb^{-1}'

        ROOT.gStyle.SetHatchesLineWidth(2)
        ROOT.gStyle.SetLineStyleString(11,"60 30");
        
        data = Data()
        
        _scale = 1.0

        _dir = 'mu_syst/'

        data.add_function('obs', _dir+'/limits.txt',
                          'b-tag SF',
                          _scale,
                          ROOT.kBlack, 1, 3,
                          ROOT.kBlack, 1.2, 8,
                          fill_style = 1002,
                          index_x = 0,
                          index_y = 1,
                          entry = 2)

        data.add_function('obs2', _dir+'/limits.txt',
                          'JER',
                          _scale,
                          ROOT.kBlue, 1, 3,
                          ROOT.kBlue, 1.2, 8,
                          fill_style = 1002,
                          index_x = 0,
                          index_y = 1,
                          entry = 1)

        data.add_function('obs3', _dir+'/limits.txt',
                          'c SF',
                          _scale,
                          ROOT.kMagenta, 1, 3,
                          ROOT.kMagenta, 1.2, 8,
                          fill_style = 1002,
                          index_x = 0,
                          index_y = 1,
                          entry = 3)

        data.add_function('obs4', _dir+'/limits.txt',
                          'light SF',
                          _scale,
                          ROOT.kCyan, 1, 3,
                          ROOT.kCyan, 1.2, 8,
                          fill_style = 1002,
                          index_x = 0,
                          index_y = 1,
                          entry = 4)

        data.add_function('obs5', _dir+'/limits.txt',
                          'pile-up',
                          _scale,
                          ROOT.kYellow+2, 1, 3,
                          ROOT.kYellow+2, 1.2, 8,
                          fill_style = 1002,
                          index_x = 0,
                          index_y = 1,
                          entry = 5)

        data.add_function('obs6', _dir+'/limits.txt',
                          'UNC',
                          _scale,
                          ROOT.kRed-6, 1, 3,
                          ROOT.kRed-6, 1.2, 8,
                          fill_style = 1002,
                          index_x = 0,
                          index_y = 1,
                          entry = 6)

        data.add_function('obs7', _dir+'/limits.txt',
                          'matching',
                          _scale,
                          ROOT.kGreen, 1, 3,
                          ROOT.kGreen, 1.2, 8,
                          fill_style = 1002,
                          index_x = 0,
                          index_y = 1,
                          entry = 7)

        data.add_function('obs8', _dir+'/limits.txt',
                          'scale',
                          _scale,
                          ROOT.kRed, 1, 3,
                          ROOT.kRed, 1.2, 8,
                          fill_style = 1002,
                          index_x = 0,
                          index_y = 1,
                          entry = 8)


        makeMultiGraph(data, 'mujets_syst_obs_'+suffix+'.'+plot_format,
                       400, 625, -1.0, 1.0,
                       xlabel = "m_{fit}", ylabel = 'fraction of expected sigma',
                       xLegend = .65, yLegend = .10,
                       legendWidth = 0.45, legendHeight = 0.40,
                       fillStyle = 1002,
                       drawOption = 'APL3', letter = '',
                       toplegend = '#font[42]{'+channel+'}')



    #------------------------------------------------------------>
    #31
    # e+jets powheg systematics study
    # relative deviation from nominal (JES+match) vs fit mass
    # observed
    #    

    if str(31) in plot_number:

        channel = 'e+jets, 4.7 fb^{-1}'

        ROOT.gStyle.SetHatchesLineWidth(2)
        ROOT.gStyle.SetLineStyleString(11,"60 30");
        
        data = Data()
        
        _scale = 1.0

        _dir = 'e_syst/'

        data.add_function('obs1', _dir+'/limits_powheg.txt',
                          'POWHEG',
                          _scale,
                          ROOT.kBlack, 1, 3,
                          ROOT.kBlack, 1.2, 8,
                          fill_style = 1002,
                          index_x = 0,
                          index_y = 1,
                          entry = 1)

        data.add_function('obs2', _dir+'/limits_powheg.txt',
                          'POWHEG nominal',
                          _scale,
                          ROOT.kBlue, 2, 3,
                          ROOT.kBlue, 0.1, 8,
                          fill_style = 1002,
                          index_x = 0,
                          index_y = 1,
                          entry = 2)

        data.add_function('obs3', _dir+'/limits_powheg.txt',
                          'POWHEG nom+var',
                          _scale,
                          ROOT.kGreen+2, 2, 3,
                          ROOT.kGreen+2, 0.1, 8,
                          fill_style = 1002,
                          index_x = 0,
                          index_y = 1,
                          entry = 3)

        makeMultiGraph(data, 'ejets_syst_powheg_obs'+suffix+'.'+plot_format,
                       400, 625, -1.5, 1.5,
                       xlabel = "m_{fit}", ylabel = 'fraction of expected sigma',
                       xLegend = .55, yLegend = .55,
                       legendWidth = 0.45, legendHeight = 0.30,
                       fillStyle = 1002,
                       drawOption = 'APL3', letter = '',
                       toplegend = '#font[42]{'+channel+'}')





    #------------------------------------------------------------>
    #32
    # e+jets powheg systematics study
    # relative deviation from nominal (JES+match) vs fit mass
    # new nominal: average of old nominal and powheg
    # observed
    #    

    if str(32) in plot_number:

        channel = 'e+jets, 4.7 fb^{-1}'

        ROOT.gStyle.SetHatchesLineWidth(2)
        ROOT.gStyle.SetLineStyleString(11,"60 30");
        
        data = Data()
        
        _scale = 1.0

        _dir = 'e_syst/'

        data.add_function('obs1', _dir+'/limits_powheg_newnom.txt',
                          'POWHEG',
                          _scale,
                          ROOT.kBlack, 1, 3,
                          ROOT.kBlack, 1.2, 8,
                          fill_style = 1002,
                          index_x = 0,
                          index_y = 1,
                          entry = 1)


        makeMultiGraph(data, 'ejets_syst_powheg_newnom_obs'+suffix+'.'+plot_format,
                       400, 625, -1.0, 1.0,
                       xlabel = "m_{fit}", ylabel = 'fraction of expected sigma',
                       xLegend = .55, yLegend = .55,
                       legendWidth = 0.45, legendHeight = 0.30,
                       fillStyle = 1002,
                       drawOption = 'APL3', letter = '',
                       toplegend = '#font[42]{'+channel+'}')



    #------------------------------------------------------------>
    #33
    # e+jets and mu+jets powheg systematics study
    # treated as asymmetric shape syst
    # expected
    #    

    if str(33) in plot_number:

        toplegend = 'POWHEG uncertainty'

        ROOT.gStyle.SetHatchesLineWidth(2)
        ROOT.gStyle.SetLineStyleString(11,"60 30");
        
        data = Data()
        
        _scale = 1.0

        data.add_function('obs1', 'e_syst/limits_powheg_newnom.txt',
                          'e+jets',
                          _scale,
                          ROOT.kBlue, 1, 3,
                          ROOT.kBlue, 1.2, 8,
                          fill_style = 1002,
                          index_x = 0,
                          index_y = 4,
                          entry = 1)


        data.add_function('obs2', 'mu_syst/limits_powheg.txt',
                          '#mu+jets',
                          _scale,
                          ROOT.kGreen+2, 1, 3,
                          ROOT.kGreen+2, 1.2, 8,
                          fill_style = 1002,
                          index_x = 0,
                          index_y = 4,
                          entry = 2)


        makeMultiGraph(data, 'powheg_syst_exp'+suffix+'.'+plot_format,
                       400, 625, -0.5, 0.5,
                       xlabel = "m_{fit}", ylabel = 'fraction of expected sigma',
                       xLegend = .55, yLegend = .55,
                       legendWidth = 0.45, legendHeight = 0.30,
                       fillStyle = 1002,
                       drawOption = 'APL3', letter = '',
                       toplegend = '#font[42]{'+toplegend+'}')



    #------------------------------------------------------------>
    #34
    # mu+jets shape systematics study
    # relative deviation from nominal (JES) vs fit mass
    # expected
    #    

    if str(34) in plot_number:

        channel = '#mu+jets, 4.6 fb^{-1}'

        ROOT.gStyle.SetHatchesLineWidth(2)
        ROOT.gStyle.SetLineStyleString(11,"60 30");
        
        data = Data()
        
        _scale = 1.0

        _dir = 'mu_syst/'

        data.add_function('obs', _dir+'/limits.txt',
                          'b-tag SF',
                          _scale,
                          ROOT.kBlack, 1, 3,
                          ROOT.kBlack, 1.2, 8,
                          fill_style = 1002,
                          index_x = 0,
                          index_y = 4,
                          entry = 2)

        data.add_function('obs2', _dir+'/limits.txt',
                          'JER',
                          _scale,
                          ROOT.kBlue, 1, 3,
                          ROOT.kBlue, 1.2, 8,
                          fill_style = 1002,
                          index_x = 0,
                          index_y = 4,
                          entry = 1)

        data.add_function('obs5', _dir+'/limits.txt',
                          'pile-up',
                          _scale,
                          ROOT.kYellow+2, 1, 3,
                          ROOT.kYellow+2, 1.2, 8,
                          fill_style = 1002,
                          index_x = 0,
                          index_y = 4,
                          entry = 5)

        data.add_function('obs6', _dir+'/limits.txt',
                          'UNC',
                          _scale,
                          ROOT.kRed-6, 1, 3,
                          ROOT.kRed-6, 1.2, 8,
                          fill_style = 1002,
                          index_x = 0,
                          index_y = 4,
                          entry = 6)

        data.add_function('obs7', _dir+'/limits.txt',
                          'matching',
                          _scale,
                          ROOT.kGreen, 1, 3,
                          ROOT.kGreen, 1.2, 8,
                          fill_style = 1002,
                          index_x = 0,
                          index_y = 4,
                          entry = 7)

        data.add_function('obs8', _dir+'/limits.txt',
                          'scale',
                          _scale,
                          ROOT.kRed, 1, 3,
                          ROOT.kRed, 1.2, 8,
                          fill_style = 1002,
                          index_x = 0,
                          index_y = 4,
                          entry = 8)


        makeMultiGraph(data, 'mujets_syst_exp_'+suffix+'.'+plot_format,
                       400, 625, -0.4, 0.6,
                       xlabel = "m_{fit}", ylabel = 'fraction of expected sigma',
                       xLegend = .65, yLegend = .50,
                       legendWidth = 0.45, legendHeight = 0.40,
                       fillStyle = 1002,
                       drawOption = 'APL3', letter = '',
                       toplegend = '#font[42]{'+channel+'}')




    #------------------------------------------------------------>
    #35
    # e+jets 20% and 10% S/B 1d comparison
    # expected
    #    

    if str(35) in plot_number:

        toplegend = 'Expected 95% C.L. upper limit, e+jets, 4.7 fb^{-1}'

        ROOT.gStyle.SetHatchesLineWidth(2)
        ROOT.gStyle.SetLineStyleString(11,"60 30");
        
        data = Data()
        
        _scale = 1.0

        data.add_with_errors('obs1', 'ejets_20pc/limits.txt',
                          '20% max bin uncertainty',
                          _scale,
                          ROOT.kBlue, 1, 3,
                          ROOT.kBlue, 1.2, 8,
                          fill_style = 1002,
                          index_x = 0,
                          index_y = 4)
 

        data.add_with_errors('obs2', 'ejets_10pc/limits.txt',
                          '10% max bin uncertainty',
                          _scale,
                          ROOT.kGreen+2, 1, 3,
                          ROOT.kGreen+2, 1.2, 8,
                          fill_style = 1002,
                          index_x = 0,
                          index_y = 4)


        makeMultiGraph(data, 'ejets_10pc_vs_20pc_'+suffix+'.'+plot_format,
                       400, 625, 0.0, 0.7,
                       xlabel = "m_{fit} [GeV]", ylabel = '#sigma [pb]',
                       xLegend = .40, yLegend = .55,
                       legendWidth = 0.45, legendHeight = 0.30,
                       fillStyle = 1002,
                       drawOption = 'APL3', letter = '',
                       toplegend = '#font[42]{'+toplegend+'}')



    #------------------------------------------------------------>
    #36
    # mu+jets 20% and 10% S/B 1d comparison
    # expected
    #    

    if str(36) in plot_number:

        toplegend = 'Expected 95% C.L. upper limit, #mu+jets, 4.6 fb^{-1}'

        ROOT.gStyle.SetHatchesLineWidth(2)
        ROOT.gStyle.SetLineStyleString(11,"60 30");
        
        data = Data()
        
        _scale = 1.0

        data.add_with_errors('obs1', 'mujets_20pc/limits.txt',
                          '20% max bin uncertainty',
                          _scale,
                          ROOT.kBlue, 1, 3,
                          ROOT.kBlue, 1.2, 8,
                          fill_style = 1002,
                          index_x = 0,
                          index_y = 4)
 

        data.add_with_errors('obs2', 'mujets_10pc/limits.txt',
                          '10% max bin uncertainty',
                          _scale,
                          ROOT.kGreen+2, 1, 3,
                          ROOT.kGreen+2, 1.2, 8,
                          fill_style = 1002,
                          index_x = 0,
                          index_y = 4)


        makeMultiGraph(data, 'mujets_10pc_vs_20pc_'+suffix+'.'+plot_format,
                       400, 625, 0.0, 0.7,
                       xlabel = "m_{fit} [GeV]", ylabel = '#sigma [pb]',
                       xLegend = .40, yLegend = .55,
                       legendWidth = 0.45, legendHeight = 0.30,
                       fillStyle = 1002,
                       drawOption = 'APL3', letter = '',
                       toplegend = '#font[42]{'+toplegend+'}')



    #------------------------------------------------------------>
    #37
    # combined e,mu+jets 20% and 10% S/B 1d comparison
    # expected
    #    

    if str(37) in plot_number:

        toplegend = 'Expected 95% C.L. upper limit, e,#mu+jets, 4.6 fb^{-1}'

        ROOT.gStyle.SetHatchesLineWidth(2)
        ROOT.gStyle.SetLineStyleString(11,"60 30");
        
        data = Data()
        
        _scale = 1.0

        data.add_with_errors('obs1', 'comb_20pc/limits.txt',
                          '20% max bin uncertainty',
                          _scale,
                          ROOT.kBlue, 1, 3,
                          ROOT.kBlue, 1.2, 8,
                          fill_style = 1002,
                          index_x = 0,
                          index_y = 4)
 

        data.add_with_errors('obs2', 'comb_10pc/limits.txt',
                          '10% max bin uncertainty',
                          _scale,
                          ROOT.kGreen+2, 1, 3,
                          ROOT.kGreen+2, 1.2, 8,
                          fill_style = 1002,
                          index_x = 0,
                          index_y = 4)


        makeMultiGraph(data, 'comb_10pc_vs_20pc_'+suffix+'.'+plot_format,
                       400, 625, 0.0, 0.7,
                       xlabel = "m_{fit} [GeV]", ylabel = '#sigma [pb]',
                       xLegend = .40, yLegend = .55,
                       legendWidth = 0.45, legendHeight = 0.30,
                       fillStyle = 1002,
                       drawOption = 'APL3', letter = '',
                       toplegend = '#font[42]{'+toplegend+'}')



