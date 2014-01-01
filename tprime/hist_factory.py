#!/usr/bin/env python
#########################################################################
#
# hist_factory.py
#
# Module for automating and running HistFactory
#
# Usage:
#     
#     ./tprime.py --hist-factory --mass 450 -f input.root
#
# Author: Gena Kukartsev, December 2011
#
#########################################################################

import os
import ROOT


comb_template_name = 'hf_tprime_template.xml'
e_template_name    = 'hf_tprime_e_template.xml'
mu_template_name   = 'hf_tprime_mu_template.xml'


class HistFactory:

    def __init__(self, combTemplateName, eTemplateName, muTemplateName):
        
        legend = '[combine_bins]:'

        self.tNames = {}
        self.tNames['comb'] = combTemplateName
        self.tNames['ejets'] = eTemplateName
        self.tNames['mujets'] = muTemplateName



    def MakeConfig(self, channelName, replacements,
                   output_dir = '.'):
        #
        # make channel config file using template file
        #
        # make main combined config if channelName == comb
        #
        # replacements - a map with replacements
        #

        output_dir = output_dir.strip().rstrip('/')
        outFileName = output_dir+'/hf_'+channelName+'_auto.xml'
        
        with open(self.tNames[channelName], 'r') as infile:
            with open(outFileName, 'w') as outfile:
                for line in infile:
                    for patt in replacements:
                        line = line.replace(patt, replacements[patt])

                    outfile.write(line)

        return



    def RunFactory(self, configFileName,
                   destFileName = 'hf_model.root',
                   output_dir = '.'):
        #
        # launch the HistFactory application
        #

        output_dir = output_dir.strip().rstrip('/')
        
        _pwd = os.getcwd()

        os.chdir(output_dir+'/')

        os.system('prepareHistFactory')

        os.system('rm -f results/tprime_ejets_tprimeCrossSection_model.root')
        os.system('rm -f results/tprime_mujets_tprimeCrossSection_model.root')
        os.system('rm -f results/tprime_combined_tprimeCrossSection_model.root')

        #os.system('hist2workspace '+output_dir+'/'+configFileName)
        #os.system('mv results/tprime_ejets_tprimeCrossSection_model.root '+output_dir+'/ejets_'+destFileName)
        #os.system('mv results/tprime_mujets_tprimeCrossSection_model.root '+output_dir+'/mujets_'+destFileName)
        #os.system('mv results/tprime_combined_tprimeCrossSection_model.root '+output_dir+'/combined_'+destFileName)

        os.system('hist2workspace '+configFileName)

        print 'DEBUG 1.1******************'
        os.system('mv results/tprime_ejets_tprimeCrossSection_model.root ejets_'+destFileName)
        os.system('mv results/tprime_mujets_tprimeCrossSection_model.root mujets_'+destFileName)
        os.system('mv results/tprime_combined_tprimeCrossSection_model.root combined_'+destFileName)
        print 'DEBUG 1.2******************'

        os.chdir(_pwd)

        return
        
        

def run( mass,
         in_file_ejets = 'data/ejets_3560/tprime_ejets_3560ipb_merged_mass_05dec2011v1.root',
         in_file_mujets = 'data/mujets_4600/tprime_mujets_3560ipb_merged_mass_05dec2011v1.root',
         destFileName = 'hf_model.root',
         output_dir = '.'
         ):
    # main public interface

    legend = '[hist_factory.run]:'

    print legend, 'starting...'

    output_dir = output_dir.strip().rstrip('/')

    hf = HistFactory(comb_template_name,
                     e_template_name,
                     mu_template_name)

    replacements = {}
    hf.MakeConfig('comb', replacements, output_dir=output_dir)

    replacements['MASS'] = str(mass)

    replacements['INPUTFILE'] = in_file_ejets
    hf.MakeConfig('ejets', replacements, output_dir=output_dir)

    replacements['INPUTFILE'] = in_file_mujets
    hf.MakeConfig('mujets', replacements, output_dir=output_dir)
    
    hf.RunFactory('hf_comb_auto.xml', destFileName, output_dir)
