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
    
def tprime_limits(plot_format = 'png'):
    #
    # e+jets
    # Profile likelihood ratio observed and expected limits
    #    

    data={}
    
    ROOT.gStyle.SetHatchesLineWidth(2)
    ROOT.gStyle.SetLineStyleString(11,"60 30");

    data[1] = Data(194)

    _scale = 1.0
    _scale_theory = 1.0
    
    #data[1].add_expected('exp', 'tprime_limit_plr_ejets_expected_134pb.ascii', 'Expected 95% C.L. limit',
    #data[1].add_expected('exp', 'tprime_limit_plr_ejets_expected_194ipb.ascii', 'Expected 95% C.L. limit',
    data[1].add_expected('exp', 'tprime_limit_plr_ejets_expected_2d_194ipb_01jul2011v1.ascii', 'Expected 95% C.L. limit',
                         _scale,
                         ROOT.kRed-3, 1, 1,
                         ROOT.kRed-3, 1.3, 0,
                         value_type = 'median',
                         error_type = 'quantile',
                         fill_style = 1002,
                         extra_scale_name = 'unit',
                         fill_color = ROOT.kYellow-7,
                         fill2_color = ROOT.kRed-3)
    #data[1].add('obs', 'tprime_limit_plr_ejets_observed_194ipb.ascii', '95% C.L. limit',
    #data[1].add('obs', 'tprime_limit_plr_ejets_observed_2d_194ipb_01jul2011v1.ascii', '95% C.L. limit',
    #            _scale,
    #            1, 1, 3,
    #            1, 1.3, 29)

    makeMultiGraph(data[1], 'tprime_limit_plr_ejets_2d_194ipb.'+plot_format,
                   350, 500, 0, 5.0,
                   xlabel = "M [GeV]", ylabel = '#sigma [pb]',
                   xLegend = .55, yLegend = .55,
                   legendWidth = 0.45, legendHeight = 0.30,
                   fillStyle = 1002,
                   drawOption = 'APL3', letter = '',
                   toplegend = '#font[42]{M vs HT, e+jets}')


    ###############################################################################################
    # mu+jets
    # Profile likelihood ratio observed and expected limits
    #    

    data={}
    
    ROOT.gStyle.SetHatchesLineWidth(2)
    ROOT.gStyle.SetLineStyleString(11,"60 30");

    data[2] = Data(684)

    _scale = 1.0
    _scale_theory = 1.0
    
    #data[2].add_expected('exp', 'tprime_limit_plr_mujets_expected_305pb.ascii', 'Expected 95% C.L. limit',
    data[2].add_expected('exp', 'tprime_limit_plr_mujets_expected_2d_719ipb_01jul2011v1.ascii', 'Expected 95% C.L. limit',
                         _scale,
                         ROOT.kRed-3, 0, 0,
                         ROOT.kRed-3, 1.3, 0,
                         value_type = 'median',
                         error_type = 'quantile',
                         fill_style = 1002,
                         extra_scale_name = 'unit',
                         fill_color = ROOT.kYellow-7,
                         fill2_color = ROOT.kRed-3)
    #data[2].add('obs', 'tprime_limit_plr_mujets_observed_305pb.ascii', '95% C.L. limit',
    data[2].add('obs', 'tprime_limit_plr_mujets_observed_2d_719ipb_01jul2011.ascii', '95% C.L. limit',
                _scale,
                1, 1, 3,
                1, 1.3, 29)

    makeMultiGraph(data[2], 'tprime_limit_plr_mujets_2d_684ipb.'+plot_format,
                   350, 500, 0, 5.0,
                   xlabel = "M [GeV]", ylabel = '#sigma [pb]',
                   xLegend = .55, yLegend = .55,
                   legendWidth = 0.45, legendHeight = 0.30,
                   fillStyle = 1002,
                   drawOption = 'APL3', letter = '',
                   toplegend = '#font[42]{M vs HT, #mu+jets}')


    #
    # e+jets combined with mu+jets
    # Profile likelihood ratio observed and expected limits
    #    

    data={}
    
    ROOT.gStyle.SetHatchesLineWidth(2)
    ROOT.gStyle.SetLineStyleString(11,"60 30");

    data[3] = Data(684)

    _scale = 1.0
    _scale_theory = 1.0
    
    #data[3].add_expected('exp', 'tprime_limit_plr_combined_expected_e134pb_mu305pb.ascii', 'Expected 95% C.L. limit',
    #data[3].add_expected('exp', 'tprime_limit_plr_combined_expected_e134ipb_mu305ipb.ascii', 'Expected 95% C.L. limit',
    data[3].add_expected('exp', 'tprime_limit_plr_combined_expected_2d_e194ipb_mu719ipb_01jul2011.ascii', 'Expected 95% C.L. limit',
                         _scale,
                         ROOT.kRed-3, 0, 0,
                         ROOT.kRed-3, 1.3, 0,
                         value_type = 'median',
                         error_type = 'quantile',
                         fill_style = 1002,
                         extra_scale_name = 'unit',
                         fill_color = ROOT.kYellow-7,
                         fill2_color = ROOT.kRed-3)
    #data[3].add('obs', 'tprime_limit_plr_combined_observed_e134pb_mu305pb.ascii', '95% C.L. limit',
    #data[3].add('obs', 'tprime_limit_plr_combined_observed_e194ipb_mu305ipb.ascii', '95% C.L. limit',
    data[3].add('obs', 'tprime_limit_plr_combined_observed_2d_e194ipb_mu719ipb.ascii', '95% C.L. limit',
                _scale,
                1, 1, 3,
                1, 1.3, 29)

    makeMultiGraph(data[3], 'tprime_limit_plr_combined_2d_e194ipb_mu684ipb.'+plot_format,
                   350, 500, 0, 5.0,
                   xlabel = "M [GeV]", ylabel = '#sigma [pb]',
                   xLegend = .55, yLegend = .55,
                   legendWidth = 0.45, legendHeight = 0.30,
                   fillStyle = 1002,
                   drawOption = 'APL3', letter = '',
                   toplegend = '#font[42]{M vs HT, e+jets and #mu+jets}')


    ###############################################################################
    #
    # e+jets
    # CLs observed and expected limits
    #    

    ROOT.gStyle.SetHatchesLineWidth(2)
    ROOT.gStyle.SetLineStyleString(11,"60 30");

    data[4] = Data(573)

    _scale = 1.0
    _scale_theory = 1.0
    
    data[4].add_with_errors('RS_exp2', 'tprime_limit_cls_ejets_expected_134pb.ascii',
                '95% expected',
                1.0,
                ROOT.kRed-3, 1, 1,
                ROOT.kRed-3, 0.1, 8,
                fill_style = 1002,
                index_x = 0,
                index_y = 6,
                index_ey_up = 8,
                index_ey_down = 4)
    data[4].add_with_errors('RS_exp1', 'tprime_limit_cls_ejets_expected_134pb.ascii',
                '68% expected',
                1.0,
                ROOT.kYellow-7, 1, 1,
                ROOT.kYellow-7, 0.1, 8,
                fill_style = 1002,
                index_x = 0,
                index_y = 6,
                index_ey_up = 7,
                index_ey_down = 5)
    data[4].add_with_errors('obs_exp_median', 'tprime_limit_cls_ejets_expected_134pb.ascii',
                'median expected',
                1.0,
                ROOT.kBlue-7, 2, 3,
                ROOT.kBlue-7, 0.1, 8,
                fill_style = 1002,
                index_x = 0,
                index_y = 6)
    data[4].add_with_errors('obs', 'tprime_limit_cls_ejets_expected_134pb.ascii',
                'data',
                1.0,
                ROOT.kBlack, 1, 3,
                ROOT.kBlack, 1.2, 8,
                fill_style = 1002,
                index_x = 0,
                index_y = 1)

    makeMultiGraph(data[4], 'tprime_limit_cls_ejets_134pb.'+plot_format,
                   350, 450, 0, 10.0,
                   xlabel = "M [GeV]", ylabel = '#sigma [pb]',
                   xLegend = .55, yLegend = .55,
                   legendWidth = 0.45, legendHeight = 0.30,
                   fillStyle = 1002,
                   drawOption = 'APL3', letter = '',
                   toplegend = '#font[42]{CLs: e+jets}')


    ###############################################################################
    #
    # mu+jets
    # CLs observed and expected limits
    #    

    ROOT.gStyle.SetHatchesLineWidth(2)
    ROOT.gStyle.SetLineStyleString(11,"60 30");

    data[5] = Data(305)

    _scale = 1.0
    _scale_theory = 1.0
    
    data[5].add_with_errors('RS_exp2', 'tprime_limit_cls_mujets_expected_305pb.ascii',
                '95% expected',
                1.0,
                ROOT.kRed-3, 1, 1,
                ROOT.kRed-3, 0.1, 8,
                fill_style = 1002,
                index_x = 0,
                index_y = 6,
                index_ey_up = 8,
                index_ey_down = 4)
    data[5].add_with_errors('RS_exp1', 'tprime_limit_cls_mujets_expected_305pb.ascii',
                '68% expected',
                1.0,
                ROOT.kYellow-7, 1, 1,
                ROOT.kYellow-7, 0.1, 8,
                fill_style = 1002,
                index_x = 0,
                index_y = 6,
                index_ey_up = 7,
                index_ey_down = 5)
    data[5].add_with_errors('obs_exp_median', 'tprime_limit_cls_mujets_expected_305pb.ascii',
                'median expected',
                1.0,
                ROOT.kBlue-7, 2, 3,
                ROOT.kBlue-7, 0.1, 8,
                fill_style = 1002,
                index_x = 0,
                index_y = 6)
    data[5].add_with_errors('obs', 'tprime_limit_cls_mujets_expected_305pb.ascii',
                'data',
                1.0,
                ROOT.kBlack, 1, 3,
                ROOT.kBlack, 1.2, 8,
                fill_style = 1002,
                index_x = 0,
                index_y = 1)

    makeMultiGraph(data[5], 'tprime_limit_cls_mujets_305pb.'+plot_format,
                   350, 450, 0, 10.0,
                   xlabel = "M [GeV]", ylabel = '#sigma [pb]',
                   xLegend = .55, yLegend = .55,
                   legendWidth = 0.45, legendHeight = 0.30,
                   fillStyle = 1002,
                   drawOption = 'APL3', letter = '',
                   toplegend = '#font[42]{CLs: #mu+jets}')


    ###############################################################################################
    #
    # mu+jets: comparison with Brussels on 500/pb
    # Profile likelihood ratio expected limits
    #    

    data={}
    
    ROOT.gStyle.SetHatchesLineWidth(2)
    ROOT.gStyle.SetLineStyleString(11,"60 30");

    data[6] = Data(500)

    _scale = 1.0
    _scale_theory = 1.0
    
    data[6].add_expected('exp', 'tprime_limit_plr_mujets_expected_2d_500ipb_01jul2011v1.ascii', 'Expected 95% C.L. limit',
                         _scale,
                         ROOT.kRed-3, 0, 0,
                         ROOT.kRed-3, 1.3, 0,
                         value_type = 'median',
                         error_type = 'quantile',
                         fill_style = 1002,
                         extra_scale_name = 'unit',
                         fill_color = ROOT.kYellow-7,
                         fill2_color = ROOT.kRed-3)


    makeMultiGraph(data[6], 'tprime_limit_plr_mujets_2d_500ipb.'+plot_format,
                   350, 500, 0, 5.0,
                   xlabel = "M [GeV]", ylabel = '#sigma [pb]',
                   xLegend = .55, yLegend = .55,
                   legendWidth = 0.45, legendHeight = 0.30,
                   fillStyle = 1002,
                   drawOption = 'APL3', letter = '',
                   toplegend = '#font[42]{M vs HT, #mu+jets}')


    ################################################################
    #
    # e+jets 1D
    # Profile likelihood ratio observed and expected limits
    #    

    data={}
    
    ROOT.gStyle.SetHatchesLineWidth(2)
    ROOT.gStyle.SetLineStyleString(11,"60 30");

    data[7] = Data(194)

    _scale = 1.0
    _scale_theory = 1.0
    
    #data[7].add_expected('exp', 'tprime_limit_plr_ejets_expected_134pb.ascii', 'Expected 95% C.L. limit',
    data[7].add_expected('exp', 'tprime_limit_plr_ejets_expected_194ipb.ascii', 'Expected 95% C.L. limit',
                         _scale,
                         ROOT.kRed-3, 1, 1,
                         ROOT.kRed-3, 1.3, 0,
                         value_type = 'median',
                         error_type = 'quantile',
                         fill_style = 1002,
                         extra_scale_name = 'unit',
                         fill_color = ROOT.kYellow-7,
                         fill2_color = ROOT.kRed-3)
    #data[7].add('obs', 'tprime_limit_plr_ejets_observed_194ipb.ascii', '95% C.L. limit',
    #            _scale,
    #            1, 1, 3,
    #            1, 1.3, 29)

    makeMultiGraph(data[7], 'tprime_limit_plr_ejets_194ipb.'+plot_format,
                   350, 500, 0, 5.0,
                   xlabel = "M [GeV]", ylabel = '#sigma [pb]',
                   xLegend = .55, yLegend = .55,
                   legendWidth = 0.45, legendHeight = 0.30,
                   fillStyle = 1002,
                   drawOption = 'APL3', letter = '',
                   toplegend = '#font[42]{1D, e+jets}')


    ###############################################################################################
    # mu+jets
    # Profile likelihood ratio observed and expected limits
    #    

    data={}
    
    ROOT.gStyle.SetHatchesLineWidth(2)
    ROOT.gStyle.SetLineStyleString(11,"60 30");

    data[8] = Data(305)

    _scale = 1.0
    _scale_theory = 1.0
    
    data[8].add_expected('exp', 'tprime_limit_plr_mujets_expected_305pb.ascii', 'Expected 95% C.L. limit',
                         _scale,
                         ROOT.kRed-3, 0, 0,
                         ROOT.kRed-3, 1.3, 0,
                         value_type = 'median',
                         error_type = 'quantile',
                         fill_style = 1002,
                         extra_scale_name = 'unit',
                         fill_color = ROOT.kYellow-7,
                         fill2_color = ROOT.kRed-3)
    data[8].add('obs', 'tprime_limit_plr_mujets_observed_305pb.ascii', '95% C.L. limit',
                _scale,
                1, 1, 3,
                1, 1.3, 29)

    makeMultiGraph(data[8], 'tprime_limit_plr_mujets_305ipb.'+plot_format,
                   350, 500, 0, 5.0,
                   xlabel = "M [GeV]", ylabel = '#sigma [pb]',
                   xLegend = .55, yLegend = .55,
                   legendWidth = 0.45, legendHeight = 0.30,
                   fillStyle = 1002,
                   drawOption = 'APL3', letter = '',
                   toplegend = '#font[42]{1D, #mu+jets}')


    ############################################################################
    #
    # 1D
    # e+jets combined with mu+jets
    # Profile likelihood ratio observed and expected limits
    #    

    data={}
    
    ROOT.gStyle.SetHatchesLineWidth(2)
    ROOT.gStyle.SetLineStyleString(11,"60 30");

    data[9] = Data(305)

    _scale = 1.0
    _scale_theory = 1.0
    
    #data[9].add_expected('exp', 'tprime_limit_plr_combined_expected_e134pb_mu305pb.ascii', 'Expected 95% C.L. limit',
    data[9].add_expected('exp', 'tprime_limit_plr_combined_expected_e134ipb_mu305ipb.ascii', 'Expected 95% C.L. limit',
                         _scale,
                         ROOT.kRed-3, 0, 0,
                         ROOT.kRed-3, 1.3, 0,
                         value_type = 'median',
                         error_type = 'quantile',
                         fill_style = 1002,
                         extra_scale_name = 'unit',
                         fill_color = ROOT.kYellow-7,
                         fill2_color = ROOT.kRed-3)
    #data[9].add('obs', 'tprime_limit_plr_combined_observed_e134pb_mu305pb.ascii', '95% C.L. limit',
    data[9].add('obs', 'tprime_limit_plr_combined_observed_e194ipb_mu305ipb.ascii', '95% C.L. limit',
                _scale,
                1, 1, 3,
                1, 1.3, 29)

    makeMultiGraph(data[9], 'tprime_limit_plr_combined_e134ipb_mu305ipb.'+plot_format,
                   350, 500, 0, 5.0,
                   xlabel = "M [GeV]", ylabel = '#sigma [pb]',
                   xLegend = .55, yLegend = .55,
                   legendWidth = 0.45, legendHeight = 0.30,
                   fillStyle = 1002,
                   drawOption = 'APL3', letter = '',
                   toplegend = '#font[42]{1D, e+jets and #mu+jets}')



def tprime_limits_mcmc(plot_format = 'png'):

    data={}
    
    ROOT.gStyle.SetHatchesLineWidth(2)
    ROOT.gStyle.SetLineStyleString(11,"60 30");


    ########################################################
    #
    # e+jets
    # Bayesian MCMC observed and expected limits
    #    

    data[1] = Data(573)

    _scale = 1.0
    _scale_theory = 1.0
    
    #data[1].add_expected('exp', 'tprime_limit_mcmc_ejets_expected_2d_573ipb_08jul2011v1.ascii', 'Expected 95% C.L. limit',
    data[1].add_expected('exp', 'tprime_limit_mcmc_ejets_expected_2d_573ipb_10jul2011v1.ascii', 'Expected 95% C.L. limit',
                         _scale,
                         ROOT.kRed-3   , 1, 1,
                         ROOT.kRed-3,1.3, 0,
                         value_type = 'median',
                         error_type = 'quantile',
                         fill_style = 1002,
                         extra_scale_name = 'unit',
                         #fill_color = ROOT.kYellow-7,
                         fill_color = ROOT.kRed-3,
                         fill2_color = ROOT.kRed-3   )

    #data[1].add('obs', 'tprime_limit_mcmc_ejets_observed_573ipb_08jul2011v1.ascii', '95% C.L. limit',
    data[1].add('obs', 'tprime_limit_mcmc_ejets_observed_573ipb_10jul2011v1.ascii', '95% C.L. limit',
                _scale,
                1, 1, 3,
                1, 1.3, 29)

    data[1].add('theory', 'tprime_theory.ascii', 't\'_{THEORY}',
                1.0,
                ROOT.kMagenta, 1, 3,
                ROOT.kMagenta, 0.1, 8)


    makeMultiGraph(data[1], 'tprime_limit_mcmc_ejets_2d_573ipb.'+plot_format,
                   350, 500, 0, 5.0,
                   xlabel = "M [GeV]", ylabel = '#sigma [pb]',
                   xLegend = .55, yLegend = .55,
                   legendWidth = 0.45, legendHeight = 0.30,
                   fillStyle = 1002,
                   drawOption = 'APL3', letter = '',
                   toplegend = '#font[42]{M vs HT, e+jets}')


    ########################################################
    #
    # mu+jets
    # Bayesian MCMC observed and expected limits
    #    

    data[2] = Data(684)

    _scale = 1.0
    _scale_theory = 1.0
    
    #data[2].add_expected('exp', 'tprime_limit_mcmc_mujets_expected_2d_684ipb_08jul2011v1.ascii', 'Expected 95% C.L. limit',
    data[2].add_expected('exp', 'tprime_limit_mcmc_mujets_expected_2d_684ipb_10jul2011v1.ascii', 'Expected 95% C.L. limit',
                         _scale,
                         ROOT.kRed-3, 1, 1,
                         ROOT.kRed-3, 1.3, 0,
                         value_type = 'median',
                         error_type = 'quantile',
                         fill_style = 1002,
                         extra_scale_name = 'unit',
                         #fill_color = ROOT.kYellow-7,
                         fill_color = ROOT.kRed-3,
                         fill2_color = ROOT.kRed-3)

    data[2].add('obs', 'tprime_limit_mcmc_mujets_observed_684ipb_10jul2011v1.ascii', '95% C.L. limit',
                _scale,
                1, 1, 3,
                1, 1.3, 29)

    data[2].add('theory', 'tprime_theory.ascii', 't\'_{THEORY}',
                1.0,
                ROOT.kMagenta, 1, 3,
                ROOT.kMagenta, 0.1, 8)

    makeMultiGraph(data[2], 'tprime_limit_mcmc_mujets_2d_684ipb.'+plot_format,
                   350, 500, 0, 5.0,
                   xlabel = "M [GeV]", ylabel = '#sigma [pb]',
                   xLegend = .55, yLegend = .55,
                   legendWidth = 0.45, legendHeight = 0.30,
                   fillStyle = 1002,
                   drawOption = 'APL3', letter = '',
                   toplegend = '#font[42]{M vs HT, #mu+jets}')


    ########################################################
    #
    # e+jets
    # Bayesian MCMC observed and expected limits
    #    

    data[3] = Data(684)

    _scale = 1.0
    _scale_theory = 1.0
    
    #data[3].add_expected('exp', 'tprime_limit_mcmc_combined_expected_2d_684ipb_08jul2011v1.ascii', 'Expected 95% C.L. limit',
    data[3].add_expected('exp', 'tprime_limit_mcmc_combined_expected_2d_mu684ipb_e573ipb_10jul2011v1.ascii', 'Expected 95% C.L. limit',
                         _scale,
                         ROOT.kRed-3, 1, 1,
                         ROOT.kRed-3, 1.3, 0,
                         value_type = 'median',
                         error_type = 'quantile',
                         fill_style = 1002,
                         extra_scale_name = 'unit',
                         #fill_color = ROOT.kYellow-7,
                         fill_color = ROOT.kRed-3,
                         fill2_color = ROOT.kRed-3)

    data[3].add('obs', 'tprime_limit_mcmc_combined_observed_mu684ipb_e573ipb_10jul2011v1.ascii', '95% C.L. limit',
                _scale,
                1, 1, 3,
                1, 1.3, 29)

    data[3].add('theory', 'tprime_theory.ascii', 't\'_{THEORY}',
                1.0,
                ROOT.kMagenta, 1, 3,
                ROOT.kMagenta, 0.1, 8)

    makeMultiGraph(data[3], 'tprime_limit_mcmc_combined_2d_mu684ipb_e573ipb.'+plot_format,
                   350, 500, 0, 5.0,
                   xlabel = "M [GeV]", ylabel = '#sigma [pb]',
                   xLegend = .55, yLegend = .55,
                   legendWidth = 0.45, legendHeight = 0.30,
                   fillStyle = 1002,
                   drawOption = 'APL3', letter = '',
                   toplegend = '#font[42]{M vs HT, e+jets and #mu+jets}')


def tprime_cls_limits(plot_format = 'png'):
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

    data[4] = Data(0.8)

    _scale = 1.0
    _scale_theory = 1.0
    
    data[4].add_with_errors('RS_exp2', 'mujets_cls_821ipb_14jul2011v1.ascii',
                '95% expected',
                1.0,
                ROOT.kYellow, 1, 1,
                ROOT.kYellow, 0.1, 8,
                fill_style = 1002,
                index_x = 0,
                index_y = 4,
                index_ey_up = 6,
                index_ey_down = 2)
    data[4].add_with_errors('RS_exp1', 'mujets_cls_821ipb_14jul2011v1.ascii',
                '68% expected',
                1.0,
                ROOT.kGreen, 1, 1,
                ROOT.kGreen, 0.1, 8,
                fill_style = 1002,
                index_x = 0,
                index_y = 4,
                index_ey_up = 5,
                index_ey_down = 3)
    data[4].add_with_errors('obs_exp_median', 'mujets_cls_821ipb_14jul2011v1.ascii',
                'median expected',
                1.0,
                ROOT.kBlack, 2, 3,
                ROOT.kBlack, 0.1, 8,
                fill_style = 1002,
                index_x = 0,
                index_y = 4)
    data[4].add_with_errors('obs', 'mujets_cls_821ipb_14jul2011v1.ascii',
                'data',
                1.0,
                ROOT.kBlack, 1, 3,
                ROOT.kBlack, 1.2, 8,
                fill_style = 1002,
                index_x = 0,
                index_y = 1)
    data[4].add('theory', 'tprime_theory.ascii', 't\'_{THEORY}',
                1.0,
                ROOT.kMagenta, 1, 3,
                ROOT.kMagenta, 0.1, 8)

    makeMultiGraph(data[4], 'tprime_limit_cls_mujets_821pb.'+plot_format,
                   350, 500, 0, 3.0,
                   xlabel = "M [GeV]", ylabel = '#sigma [pb]',
                   xLegend = .55, yLegend = .55,
                   legendWidth = 0.45, legendHeight = 0.30,
                   fillStyle = 1002,
                   drawOption = 'APL3', letter = '',
                   toplegend = '#font[42]{CL_{S}: #mu+jets}')


    ###############################################################################
    #
    # e+jets
    # CLs observed and expected limits
    #    

    ROOT.gStyle.SetHatchesLineWidth(2)
    ROOT.gStyle.SetLineStyleString(11,"60 30");

    data[5] = Data(0.6)

    _scale = 1.0
    _scale_theory = 1.0
    
    data[5].add_with_errors('RS_exp2', 'ejets_cls_573ipb_14jul2011v1.ascii',
                '95% expected',
                1.0,
                ROOT.kYellow, 1, 1,
                ROOT.kYellow, 0.1, 8,
                fill_style = 1002,
                index_x = 0,
                index_y = 4,
                index_ey_up = 6,
                index_ey_down = 2)
    data[5].add_with_errors('RS_exp1', 'ejets_cls_573ipb_14jul2011v1.ascii',
                '68% expected',
                1.0,
                ROOT.kGreen, 1, 1,
                ROOT.kGreen, 0.1, 8,
                fill_style = 1002,
                index_x = 0,
                index_y = 4,
                index_ey_up = 5,
                index_ey_down = 3)
    data[5].add_with_errors('obs_exp_median', 'ejets_cls_573ipb_14jul2011v1.ascii',
                'median expected',
                1.0,
                ROOT.kBlack, 2, 3,
                ROOT.kBlack, 0.1, 8,
                fill_style = 1002,
                index_x = 0,
                index_y = 4)
    data[5].add_with_errors('obs', 'ejets_cls_573ipb_14jul2011v1.ascii',
                'data',
                1.0,
                ROOT.kBlack, 1, 3,
                ROOT.kBlack, 1.2, 8,
                fill_style = 1002,
                index_x = 0,
                index_y = 1)
    data[5].add('theory', 'tprime_theory.ascii', 't\'_{THEORY}',
                1.0,
                ROOT.kMagenta, 1, 3,
                ROOT.kMagenta, 0.1, 8)

    makeMultiGraph(data[5], 'tprime_limit_cls_ejets_573pb.'+plot_format,
                   350, 500, 0, 3.0,
                   xlabel = "M [GeV]", ylabel = '#sigma [pb]',
                   xLegend = .55, yLegend = .55,
                   legendWidth = 0.45, legendHeight = 0.30,
                   fillStyle = 1002,
                   drawOption = 'APL3', letter = '',
                   toplegend = '#font[42]{CL_{S}: e+jets}')

    ###############################################################################
    #
    # mu+jets and e+jets
    # CLs observed and expected limits
    #    

    ROOT.gStyle.SetHatchesLineWidth(2)
    ROOT.gStyle.SetLineStyleString(11,"60 30");

    data[6] = Data(0.8)

    _scale = 1.0
    _scale_theory = 1.0
    
    data[6].add_with_errors('RS_exp2', 'comb_cls_mu821ipb_e573ipb_14jul2011v1.ascii',
                '95% expected',
                1.0,
                ROOT.kYellow, 1, 1,
                ROOT.kYellow, 0.1, 8,
                fill_style = 1002,
                index_x = 0,
                index_y = 4,
                index_ey_up = 6,
                index_ey_down = 2)
    data[6].add_with_errors('RS_exp1', 'comb_cls_mu821ipb_e573ipb_14jul2011v1.ascii',
                '68% expected',
                1.0,
                ROOT.kGreen, 1, 1,
                ROOT.kGreen, 0.1, 8,
                fill_style = 1002,
                index_x = 0,
                index_y = 4,
                index_ey_up = 5,
                index_ey_down = 3)
    data[6].add_with_errors('obs_exp_median', 'comb_cls_mu821ipb_e573ipb_14jul2011v1.ascii',
                'median expected',
                1.0,
                ROOT.kBlack, 2, 3,
                ROOT.kBlack, 0.1, 8,
                fill_style = 1002,
                index_x = 0,
                index_y = 4)
    data[6].add_with_errors('obs', 'comb_cls_mu821ipb_e573ipb_14jul2011v1.ascii',
                'data',
                1.0,
                ROOT.kBlack, 1, 3,
                ROOT.kBlack, 1.2, 8,
                fill_style = 1002,
                index_x = 0,
                index_y = 1)
    data[6].add('theory', 'tprime_theory.ascii', 't\'_{THEORY}',
                1.0,
                ROOT.kMagenta, 1, 3,
                ROOT.kMagenta, 0.1, 8)

    makeMultiGraph(data[6], 'tprime_limit_cls_comb_e573ipb_mu821ipb.'+plot_format,
                   350, 500, 0, 3.0,
                   xlabel = "M [GeV]", ylabel = '#sigma [pb]',
                   xLegend = .55, yLegend = .55,
                   legendWidth = 0.45, legendHeight = 0.30,
                   fillStyle = 1002,
                   drawOption = 'APL3', letter = '',
                   toplegend = '#font[42]{CL_{S}: #mu+jets and e+jets}')
