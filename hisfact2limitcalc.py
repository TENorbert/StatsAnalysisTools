 #!/usr/bin/env python

  #
  # A pyROOT script demonstrating
  # an example of writing a HistFactory
  # model using python
  #
  # This example was written to match
  # the example.xml analysis in
  # $ROOTSYS/tutorials/histfactory/
  #
  # Orignally Written by George Lewis
  # Adjusted for DelayedPhoton by TEN-UMN
  #

#def main():

# try:
#           import ROOT
# except:
#           print "It seems that pyROOT isn't properly configured"
#           return


import sys
import optparse
import shutil
import os
import re
import ROOT 
from ROOT import *
#from array import array

"""
Create a HistFactory measurement from python
"""
""" Creat Hist/RootFiles for CLs Using HistFactory"""
__version__ = "1.0"

data_bg_InputFile = "./data_bg_file.root"
sg_InputFile      = "./sig_gmbs600.root"

# Create the measurement
meas = ROOT.RooStats.HistFactory.Measurement("meas", "meas")
meas.SetOutputFilePrefix("./results/Delayed_Photon_UsingPy")
meas.SetPOI("SigXsecOverSM")
meas.AddConstantParam("Lumi")
meas.AddConstantParam("alpha_syst1")
meas.SetLumi(1.0)
meas.SetLumiRelErr(0.10)
meas.SetExportOnly(False)

# Create a channel

chan = ROOT.RooStats.HistFactory.Channel("channel1")
chan.SetData("h_dataTime", data_bg_InputFile)
chan.SetStatErrorConfig(0.05, "Poisson")

# Now, create some samples

# Create the signal sample
signal = ROOT.RooStats.HistFactory.Sample("h_sgTime__ctau6000_hehb", "h_sgTime__ctau6000_hehb", sg_InputFile)
signal.AddOverallSys("syst1",  0.95, 1.05)
signal.AddNormFactor("SigXsecOverSM", 1, 0, 3)
chan.AddSample(h_sgTime__ctau6000_hehb)


# Background 1
background1 = ROOT.RooStats.HistFactory.Sample("h_bgTime", "h_bgTime", data_bg_InputFile)
# background1.ActivateStatError("background1_statUncert", data_bg_InputFile)
background1.AddOverallSys("syst2", 0.95, 1.05 )
chan.AddSample(h_bgTime)


# Background 2
# background2 = ROOT.RooStats.HistFactory.Sample("background2", "background2", data_bg_InputFile)
# background2.ActivateStatError()
# background2.AddOverallSys("syst3", 0.95, 1.05 )
# chan.AddSample(background2)


# Done with this channel
# Add it to the measurement:

meas.AddChannel(chan)

# Collect the histograms from their files,
# print some output, 
meas.CollectHistograms()
meas.PrintTree();

# One can print XML code to an
# output directory:
# meas.PrintXML("xmlFromCCode", meas.GetOutputFilePrefix());

meas.PrintXML("xmlFromPy", meas.GetOutputFilePrefix());

# Now, do the measurement
ROOT.RooStats.HistFactory.MakeModelAndMeasurementFast(meas);

# pass


# if __name__ == "__main__":
#    main()
