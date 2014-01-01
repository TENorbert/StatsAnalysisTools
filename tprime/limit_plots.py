#!/usr/bin/env python

################################################
#
# This script plots limits for the di-lepton
# and di-photon analyses
#
#   - input data from ASCII files
#   
#
# Gena Kukartsev, 2011
#
################################################

import ROOT
from plotter import *

def fit(plot_format = 'png'):
    #
    # Fit k-factor
    #    

    data={}
    
    ROOT.gStyle.SetHatchesLineWidth(2)
    ROOT.gStyle.SetLineStyleString(11,"60 30");

    data[1] = Data(9.9)

    _scale = 1.0
    _scale_theory = 1.0
    
    data[1].add('obs', 'theory/belyaev_tab_lhc_d_cteq_20dec2010v2.txt', 'CTEQ6L k-factor',
                _scale,
                1, 1, 3,
                1, 1.3, 29,
                dofit=True, fit_min = 700, fit_max = 2500)

    makeMultiGraph(data[1], 'kfactor_fit.'+plot_format,
                   100, 3000, 1.0, 1.5,
                   xlabel = "M [GeV]", ylabel = 'k-factor',
                   xLegend = .55, yLegend = .40,
                   legendWidth = 0.45, legendHeight = 0.50,
                   fillStyle = 1002,
                   drawOption = 'APC3', letter = '',
                   toplegend = '#font[42]{CTEQ6L k-factor}')






def theory_plots(plot_format = 'png'):
    #
    # Graviton theory curves
    #    

    data={}
    
    ROOT.gStyle.SetHatchesLineWidth(2)
    ROOT.gStyle.SetLineStyleString(11,"60 30");

    data[1] = Data(36.1)
#
#    _scale = 1.0
#    _scale_theory = 1.0
#    
#    data[1].add_expected('exp', 'expected_tmp.ascii', 'Expected 95% C.L. limit',
#                         _scale,
#                         ROOT.kRed-3, 0, 0,
#                         ROOT.kRed-3, 1.3, 0,
#                         value_type = 'median',
#                         error_type = 'quantile',
#                         fill_style = 1002,
#                         extra_scale_name = 'unit',
#                         fill_color = ROOT.kYellow-7,
#                         fill2_color = ROOT.kRed-3)
#    data[1].add('RS0.1', 'rsg01_mumu.txt', 'G_{#font[132]{KK}} k/#bar{M_{#font[132]{Pl}}}=0.1',
#                1.0e+9*_scale_theory,
#                ROOT.kBlue-1, 1, 1,
#                ROOT.kBlue-1, 0.1, 8,
#                var_scale = 'graviton',
#                #fill_style = 3013,
#                fill_style = 1002,
#                fill_color = ROOT.kBlue+2)
#    data[1].add('RS0.05', 'rsg005_mumu.txt', 'G_{#font[132]{KK}} k/#bar{M_{#font[132]{Pl}}}=0.05',
#                1.0e+9*_scale_theory,
#                ROOT.kGreen-6, 1, 1,
#                ROOT.kGreen -6, 0.1, 8,
#                var_scale = 'graviton',
#                fill_style = 1002)
#    data[1].add('obs', 'observed_tmp.ascii', '95% C.L. limit',
#                _scale,
#                1, 1, 3,
#                1, 1.3, 29)
#
#    makeMultiGraph(data[1], 'gkk_obs.'+plot_format,
#                   300, 1450, 0, 0.5,
#                   xlabel = "M [GeV]", ylabel = '#sigma [pb]',
#                   xLegend = .55, yLegend = .40,
#                   legendWidth = 0.45, legendHeight = 0.50,
#                   fillStyle = 1002,
#                   drawOption = 'APC3', letter = '',
#                   toplegend = '#font[42]{#gamma#gamma+#mu^{+}#mu^{-}+ee}')
#
#




def limit_plots(plot_format = 'png'):
    #
    # Dielectron limits
    #    

    data={}
    
    ROOT.gStyle.SetHatchesLineWidth(2)
    ROOT.gStyle.SetLineStyleString(11,"60 30");

    data[2] = Data(1.1)

    _scale = 10.0
    #_scale = 10.0*1.02
    _scale_theory = 10000.0/970.0
    
    #data[2].add_expected('exp', 'dielectron_expected_714pb.ascii', 'Expected 95% C.L. limit',
    #data[2].add_expected('exp', 'dielectron_mcmc_1100ipb_exp_10jul2011.ascii', 'Expected 95% C.L. limit',
    data[2].add_expected('exp', 'dielectron_mcmc_expected_1100ipb_16jul2011.ascii', 'Expected 95% C.L. limit',
                         _scale,
                         ROOT.kRed-3, 0, 0,
                         ROOT.kRed-3, 1.3, 0,
                         value_type = 'median',
                         error_type = 'quantile',
                         fill_style = 1002,
                         extra_scale_name = 'unit',
                         fill_color = ROOT.kYellow-7,
                         fill2_color = ROOT.kRed-3)
    #data[2].add('SSM', 'theory/cs_ssm_belyaev_15jan2011v1.dat', 'Z\'_{SSM}',
    data[2].add('SSM', 'theory/zprime_ssm_cteq6l1_masswindow.txt', 'Z\'_{SSM}',
                1.0e+9*_scale_theory,
                ROOT.kGreen-6, 1, 3,
                ROOT.kGreen -6, 0.1, 8,
                var_scale = 'cteq6l1_fit')
    #data[2].add('Psi', 'theory/cs_psi_belyaev_15jan2011v1.dat', 'Z\'_{#Psi}',
    data[2].add('Psi', 'theory/zprime_psi_cteq6l1_masswindow.txt', 'Z\'_{#Psi}',
                1.0e+9*_scale_theory,
                ROOT.kBlue, 1, 3,
                ROOT.kBlue, 0.1, 8,
                var_scale = 'cteq6l1_fit')
    #data[2].add('RS0.1', 'rsg01_mumu.txt', 'G_{#font[132]{KK}} k/#bar{M_{#font[132]{Pl}}}=0.1',
    data[2].add('RS0.1', 'theory/rsg01_cteq6l1_masswindow.txt', 'G_{#font[132]{KK}} k/#bar{M_{#font[132]{Pl}}}=0.1',
                1.0e+9*_scale_theory,
                ROOT.kBlue-1, 1, 1,
                ROOT.kBlue-1, 0.1, 8,
                var_scale = 'graviton',
                #fill_style = 3013,
                fill_style = 1002,
                fill_color = ROOT.kBlue+2)
    #data[2].add('RS0.05', 'rsg005_mumu.txt', 'G_{#font[132]{KK}} k/#bar{M_{#font[132]{Pl}}}=0.05',
    data[2].add('RS0.05', 'theory/rsg005_cteq6l1_masswindow.txt', 'G_{#font[132]{KK}} k/#bar{M_{#font[132]{Pl}}}=0.05',
                1.0e+9*_scale_theory,
                ROOT.kGreen-6, 1, 1,
                ROOT.kGreen -6, 0.1, 8,
                var_scale = 'graviton',
                fill_style = 1002)
    #data[2].add_median('obs', 'dielectron_mcmc_mass_graviton_1100ipb_16jul2011.ascii', '95% C.L. limit',
    #data[2].add_median('obs', 'dielectron_mcmc_mass_cteq_1100ipb_16jul2011.ascii', '95% C.L. limit',
    data[2].add_median('obs', 'dielectron_mcmc_observed_1100ipb_16jul2011.ascii', '95% C.L. limit',
    #data[2].add('obs', 'dielectron_mass_cteq_1100ipb.ascii', '95% C.L. limit',
    #data[2].add('obs', 'dielectron_mass_graviton_1100ipb.ascii', '95% C.L. limit',
    #data[2].add('obs', 'dielectron_mcmc_1100ipb_10jul2011.ascii', '95% C.L. limit',
    #data[2].add('obs', 'dielectron_mcmc_1100ipb_obs_10jul2011.ascii', '95% C.L. limit',
    #data[2].add('obs', 'dielectron_observed_714pb.ascii', '95% C.L. limit',
    #data[2].add('obs', 'dielectron_cteq_714ipb.ascii', '95% C.L. limit',
    #data[2].add('obs', 'dielectron_graviton_714ipb.ascii', '95% C.L. limit',
                _scale,
                1, 1, 3,
                1, 0.1, 29)

    makeMultiGraph(data[2], 'zprime_dielectron_1100ipb.'+plot_format,
                   300, 2200, 0, 0.4,
                   xlabel = "M [GeV]", ylabel = 'R_{#sigma} #upoint 10^{-4}',
                   xLegend = .59, yLegend = .40,
                   legendWidth = 0.45, legendHeight = 0.50,
                   fillStyle = 1002,
                   drawOption = 'APC3', letter = '',
                   toplegend = '#font[42]{ee}',
                   #find_intersection = True,
                   #draw_limit_line = True
                   )


    #########################################################
    #
    # Dimuon limits
    #    

    #data={}

    _scale = 10.0
    #_scale = 10.0*1.048
    #_scale = 10.0*1.048*1061/942
    
    ROOT.gStyle.SetHatchesLineWidth(2)
    ROOT.gStyle.SetLineStyleString(11,"60 30");

    data[3] = Data(1.1)

    #data[3].add_expected('exp', 'dimuon_expected_746ipb.ascii', 'Expected 95% C.L. limit',
    #data[3].add_expected('exp', 'dimuon_expected_1125ipb.ascii', 'Expected 95% C.L. limit',
    #data[3].add_expected('exp', 'dimuon_mcmc_1100ipb_exp_10jul2011.ascii', 'Expected 95% C.L. limit',
    data[3].add_expected('exp', 'dimuon_mcmc_expected_1100ipb_17jul2011.ascii', 'Expected 95% C.L. limit',
                         _scale,
                         ROOT.kRed-3, 0, 0,
                         ROOT.kRed-3, 1.3, 0,
                         value_type = 'median',
                         error_type = 'quantile',
                         fill_style = 1002,
                         extra_scale_name = 'unit',
                         fill_color = ROOT.kYellow-7,
                         fill2_color = ROOT.kRed-3)
    #data[3].add('SSM', 'theory/cs_ssm_belyaev_15jan2011v1.dat', 'Z\'_{SSM}',
    data[3].add('SSM', 'theory/zprime_ssm_cteq6l1_masswindow.txt', 'Z\'_{SSM}',
                1.0e+9*_scale_theory,
                ROOT.kGreen-6, 1, 3,
                ROOT.kGreen -6, 0.1, 8,
                var_scale = 'cteq6l1_fit')
    #data[3].add('Psi', 'theory/cs_psi_belyaev_15jan2011v1.dat', 'Z\'_{#Psi}',
    data[3].add('Psi', 'theory/zprime_psi_cteq6l1_masswindow.txt', 'Z\'_{#Psi}',
                1.0e+9*_scale_theory,
                ROOT.kBlue, 1, 3,
                ROOT.kBlue, 0.1, 8,
                var_scale = 'cteq6l1_fit')
    #data[3].add('RS0.1', 'rsg01_mumu.txt', 'G_{#font[132]{KK}} k/#bar{M_{#font[132]{Pl}}}=0.1',
    data[3].add('RS0.1', 'theory/rsg01_cteq6l1_masswindow.txt', 'G_{#font[132]{KK}} k/#bar{M_{#font[132]{Pl}}}=0.1',
                1.0e+9*_scale_theory,
                ROOT.kBlue-1, 1, 1,
                ROOT.kBlue-1, 0.1, 8,
                var_scale = 'graviton',
                #fill_style = 3013,
                fill_style = 1002,
                fill_color = ROOT.kBlue+2)
    #data[3].add('RS0.05', 'rsg005_mumu.txt', 'G_{#font[132]{KK}} k/#bar{M_{#font[132]{Pl}}}=0.05',
    data[3].add('RS0.05', 'theory/rsg005_cteq6l1_masswindow.txt', 'G_{#font[132]{KK}} k/#bar{M_{#font[132]{Pl}}}=0.05',
                1.0e+9*_scale_theory,
                ROOT.kGreen-6, 1, 1,
                ROOT.kGreen -6, 0.1, 8,
                var_scale = 'graviton',
                fill_style = 1002)
    #data[3].add_median('obs', 'dimuon_mcmc_mass_graviton_1100ipb_16jul2011.ascii', '95% C.L. limit',
    #data[3].add_median('obs', 'dimuon_mcmc_mass_cteq_1100ipb_16jul2011.ascii', '95% C.L. limit',
    data[3].add_median('obs', 'dimuon_mcmc_observed_1100ipb_16jul2011.ascii', '95% C.L. limit',
    #data[3].add('obs', 'dimuon_mass_cteq_1100ipb.ascii', '95% C.L. limit',
    #data[3].add('obs', 'dimuon_mass_graviton_1100ipb.ascii', '95% C.L. limit',
    #data[3].add('obs', 'dimuon_mcmc_1100ipb_10jul2011.ascii', '95% C.L. limit',
    #data[3].add('obs', 'dimuon_observed_1125ipb.ascii', '95% C.L. limit',
    #data[3].add('obs', 'dimuon_observed_746ipb.ascii', '95% C.L. limit',
    #data[3].add('obs', 'dimuon_cteq_746ipb.ascii', '95% C.L. limit',
    #data[3].add('obs', 'dimuon_graviton_746ipb.ascii', '95% C.L. limit',
                _scale,
                1, 1, 3,
                1, 0.1, 29)

    makeMultiGraph(data[3], 'zprime_dimuon_1125ipb.'+plot_format,
                   300, 2200, 0, 0.4,
                   xlabel = "M [GeV]", ylabel = 'R_{#sigma} #upoint 10^{-4}',
                   xLegend = .59, yLegend = .40,
                   legendWidth = 0.45, legendHeight = 0.50,
                   fillStyle = 1002,
                   drawOption = 'APC3', letter = '',
                   toplegend = '#font[42]{#mu^{+}#mu^{-}}',
                   #find_intersection = True,
                   #draw_limit_line = True
                   )


    ##########################################################################
    #
    # Dilepton limits
    #    

    #data={}
    
    ROOT.gStyle.SetHatchesLineWidth(2)
    ROOT.gStyle.SetLineStyleString(11,"60 30");

    data[4] = Data(1.1)

    #data[4].add_expected('exp', 'dilepton_expected_mu746ipb_e714ipb.ascii', 'Expected 95% C.L. limit',
    #data[4].add_expected('exp', 'dilepton_mcmc_1100ipb_exp_10jul2011.ascii', 'Expected 95% C.L. limit',
    data[4].add_expected('exp', 'dilepton_mcmc_expected_1100ipb_16jul2011.ascii', 'Expected 95% C.L. limit',
                         _scale,
                         ROOT.kRed-3, 0, 0,
                         ROOT.kRed-3, 1.3, 0,
                         value_type = 'median',
                         error_type = 'quantile',
                         fill_style = 1002,
                         extra_scale_name = 'unit',
                         fill_color = ROOT.kYellow-7,
                         fill2_color = ROOT.kRed-3)
    #data[4].add('SSM', 'theory/cs_ssm_belyaev_15jan2011v1.dat', 'Z\'_{SSM}',
    data[4].add('SSM', 'theory/zprime_ssm_cteq6l1_masswindow.txt', 'Z\'_{SSM}',
                1.0e+9*_scale_theory,
                ROOT.kGreen-6, 1, 3,
                ROOT.kGreen -6, 0.1, 8,
                var_scale = 'cteq6l1_fit')
    #data[4].add('Psi', 'theory/cs_psi_belyaev_15jan2011v1.dat', 'Z\'_{#Psi}',
    data[4].add('Psi', 'theory/zprime_psi_cteq6l1_masswindow.txt', 'Z\'_{#Psi}',
                1.0e+9*_scale_theory,
                ROOT.kBlue, 1, 3,
                ROOT.kBlue, 0.1, 8,
                var_scale = 'cteq6l1_fit')
    #data[4].add('RS0.1', 'rsg01_mumu.txt', 'G_{#font[132]{KK}} k/#bar{M_{#font[132]{Pl}}}=0.1',
    data[4].add('RS0.1', 'theory/rsg01_cteq6l1_masswindow.txt', 'G_{#font[132]{KK}} k/#bar{M_{#font[132]{Pl}}}=0.1',
                1.0e+9*_scale_theory,
                ROOT.kBlue-1, 1, 1,
                ROOT.kBlue-1, 0.1, 8,
                var_scale = 'graviton',
                #fill_style = 3013,
                fill_style = 1002,
                fill_color = ROOT.kBlue+2)
    #data[4].add('RS0.05', 'rsg005_mumu.txt', 'G_{#font[132]{KK}} k/#bar{M_{#font[132]{Pl}}}=0.05',
    data[4].add('RS0.05', 'theory/rsg005_cteq6l1_masswindow.txt', 'G_{#font[132]{KK}} k/#bar{M_{#font[132]{Pl}}}=0.05',
                1.0e+9*_scale_theory,
                ROOT.kGreen-6, 1, 1,
                ROOT.kGreen -6, 0.1, 8,
                var_scale = 'graviton',
                fill_style = 1002)
    #data[4].add_median('obs', 'dilepton_mcmc_mass_graviton_1100ipb_16jul2011.ascii', '95% C.L. limit',
    #data[4].add_median('obs', 'dilepton_mcmc_mass_cteq_1100ipb_16jul2011.ascii', '95% C.L. limit',
    data[4].add_median('obs', 'dilepton_mcmc_observed_1100ipb_16jul2011.ascii', '95% C.L. limit',
    #data[4].add('obs', 'dilepton_mass_cteq_1100ipb.ascii', '95% C.L. limit',
    #data[4].add('obs', 'dilepton_mass_graviton_1100ipb.ascii', '95% C.L. limit',
    #data[4].add('obs', 'dilepton_mcmc_1100ipb_10jul2011.ascii', '95% C.L. limit',
    #data[4].add('obs', 'dilepton_mcmc_1100ipb_obs_10jul2011.ascii', '95% C.L. limit',
    #data[4].add('obs', 'dilepton_observed_mu746ipb_e714ipb.ascii', '95% C.L. limit',
    #data[4].add('obs', 'dilepton_cteq_mu746ipb_e714ipb.ascii', '95% C.L. limit',
    #data[4].add('obs', 'dilepton_graviton_mu746ipb_e714ipb.ascii', '95% C.L. limit',
                _scale,
                1, 1, 3,
                1, 0.1, 29)

    makeMultiGraph(data[4], 'zprime_dilepton_1100ipb.'+plot_format,
                   300, 2200, 0, 0.4,
                   xlabel = "M [GeV]", ylabel = 'R_{#sigma} #upoint 10^{-4}',
                   xLegend = .59, yLegend = .40,
                   legendWidth = 0.45, legendHeight = 0.50,
                   fillStyle = 1002,
                   drawOption = 'APC3', letter = '',
                   toplegend = '#font[42]{ee+#mu^{+}#mu^{-}}',
                   #find_intersection = True,
                   #draw_limit_line = True
                   )



def plot_info(dir,masspts):
    import os
    for j in range(0,masspts):
        mass = str(50*j + 350)
        if "350" in mass:
            os.system("head -2 "+dir+"*"+mass+"*.ascii.ascii")
        os.system("tail -1 "+dir+"*"+mass+"*ascii.ascii")


def plot_comp(dir,dir2,masspts):
    import os
    import glob
    for j in range(0,masspts):
        mass = str(50*j + 350)
        if "350" in mass:
            os.system("head -2 "+dir+"*"+mass+"*.ascii.ascii")

        filepath = dir+"*"+mass+"*.ascii.ascii"
        filepath2= dir2+"*"+mass+"*.ascii.ascii"
        filename = glob.glob(filepath)
        filename2= glob.glob(filepath2)
        
        file = open(filename[0],'r')
        file2= open(filename2[0],'r')
        file.readline()
        file.readline()
        file2.readline()
        file2.readline()
        
        numbers = file.readline().split()
        numbers2= file2.readline().split()

        print ('%.*s' % (3, numbers[0])),
        for n in range(1,len(numbers)):
            f = float(float(numbers[n])/float(numbers2[n])) 
            print "\t",('%.*f' % (3, f)),"\t",
        print

      
