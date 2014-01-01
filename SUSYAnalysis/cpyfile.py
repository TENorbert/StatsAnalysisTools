#!/usr/bin/python 
#simple Python file to copy files

import sys
import  os
import shutil

myfilelist=[ "/uscms/home/jhirsch/higgsLim_428/CMSSW_4_2_8/test/susy_das_inputfiles/data_5ifb_2011_das.root",
"/uscms/home/jhirsch/higgsLim_428/CMSSW_4_2_8/test/susy_das_inputfiles/mc_SS10_das.root",
"/uscms/home/jhirsch/higgsLim_428/CMSSW_4_2_8/test/susy_das_inputfiles/mc_SS11_das.root",
"/uscms/home/jhirsch/higgsLim_428/CMSSW_4_2_8/test/susy_das_inputfiles/mc_SS12_das.root",
"/uscms/home/jhirsch/higgsLim_428/CMSSW_4_2_8/test/susy_das_inputfiles/mc_SS13_das.root",
"/uscms/home/jhirsch/higgsLim_428/CMSSW_4_2_8/test/susy_das_inputfiles/mc_SS14_das.root",
"/uscms/home/jhirsch/higgsLim_428/CMSSW_4_2_8/test/susy_das_inputfiles/mc_SS15_das.root",
"/uscms/home/jhirsch/higgsLim_428/CMSSW_4_2_8/test/susy_das_inputfiles/mc_SS16_das.root",
"/uscms/home/jhirsch/higgsLim_428/CMSSW_4_2_8/test/susy_das_inputfiles/mc_SS17_das.root",
"/uscms/home/jhirsch/higgsLim_428/CMSSW_4_2_8/test/susy_das_inputfiles/mc_SS18_das.root",
"/uscms/home/jhirsch/higgsLim_428/CMSSW_4_2_8/test/susy_das_inputfiles/mc_SS2_das.root",
"/uscms/home/jhirsch/higgsLim_428/CMSSW_4_2_8/test/susy_das_inputfiles/mc_SS3_das.root",
"/uscms/home/jhirsch/higgsLim_428/CMSSW_4_2_8/test/susy_das_inputfiles/mc_SS4_das.root",
"/uscms/home/jhirsch/higgsLim_428/CMSSW_4_2_8/test/susy_das_inputfiles/mc_SS5_das.root",
"/uscms/home/jhirsch/higgsLim_428/CMSSW_4_2_8/test/susy_das_inputfiles/mc_SS6_das.root",
"/uscms/home/jhirsch/higgsLim_428/CMSSW_4_2_8/test/susy_das_inputfiles/mc_SS7_das.root",
"/uscms/home/jhirsch/higgsLim_428/CMSSW_4_2_8/test/susy_das_inputfiles/mc_SS8_das.root",
"/uscms/home/jhirsch/higgsLim_428/CMSSW_4_2_8/test/susy_das_inputfiles/mc_SS9_das.root"]

distdir = os.mkdir("SusyFiles")
#for file in range(len(myfilelist)):
for file in myfilelist[:]:
#  os.chdir("SusyFiles")
  shutil.copy2(file, distdir)




