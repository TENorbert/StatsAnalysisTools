#include "RooStats/HistFactory/Measurement.h"

using namespace RooFit;
using namespace RooStats;

   void NeutralinoDecaySimpleModel() {

   // Create a simple 1-channel model
   // using c++ and ROOT

   // Create the measurement object
   // This is the top node of the structure
   // We do some minor configuration as well
   RooStats::HistFactory::Measurement meas("my_measurement", "my measurement");

   // Set the prefix that will appear before
   // all output for this measurement
   // We Set ExportOnly to false, meaning
   // we will fit the measurement and make 
   // plots in addition to saving the workspace
   meas.SetOutputPrefix("results/my_measurement");
   meas.SetExportOnly(False);

   // Set the name of the parameter of interest
   // Note that this parameter hasn't yet been
   // created, we are anticipating it
   meas.SetPOI("SigXsecOverSM");

   // Set the luminosity
   // There are a few conventions for this.
   // Here, we assume that all histograms have
   // already been scaled by luminosity
   // We also set a 10% uncertainty
   meas.SetLumi(1.0);
   meas.SetLumiRelErr(0.10);

   // Okay, now that we've configured the measurement,
   // we'll start building the tree.
   // We begin by creating the first channel
   RooStats::HistFactory::Channel chan("channel");

   // First, we set the 'data' for this channel
   // The data is a histogram represeting the 
   // measured distribution.  It can have 1 or many bins.
   // In this example, we assume that the data histogram
   // is already made and saved in a ROOT file.  
   // So, to 'set the data', we give this channel the
   // path to that ROOT file and the name of the data
   // histogram in that root file
   // The arguments are: SetData(HistogramName, HistogramFile)
   chan.SetData("h_dataTime", "./data_bg_file.root");


   // Now that we have a channel and have attached
   // data to it, we will start creating our Samples
   // These describe the various processes that we
   // use to model the data.
   // Here, they just consist of a signal process
   // and a single background process.
   RooStats::HistFactory::Sample signal("h_sgTime__ctau6000_hehb", "h_sgTime__ctau6000_hehb", "./sig_gmbs600.root");

   // Having created this sample, we configure it
   // First, we add the cross-section scaling
   // parameter that we call SigXsecOverSM
   // Then, we add a systematic with a 5% uncertainty
   // Finally, we add it to our channel
   signal.AddNormFactor("SigXsecOverSM", 1, 0, 3);
   signal.AddOverallSys("syst1",  0.95, 1.05);
   chan.AddSample(signal);     


   // We do a similar thing for our background
   RooStats::HistFactory::Sample background1("h_bgTime", "h_bgTime", "./data_bg_file.root");
  // background1.ActivateStatError("background1_statUncert", InputFile);
   background1.AddOverallSys("syst2", 0.95, 1.05 );
   chan.AddSample(h_bgTime);         

/*   // And we create a second background for good measure
   Sample background2("background2", "background2", "data/example.root");
   background2.ActivateStatError();
   background2.AddOverallSys("syst3", 0.95, 1.05 );
   chan.AddSample(background2);
*/

   // Now that we have fully configured our channel,
   // we add it to the main measurement
   meas.AddChannel(chan);

   // At this point, we have only given our channel
   // and measurement the input histograms as strings
   // We must now have the measurement open the files,
   // collect the histograms, copy and store them.
   // This step involves I/O 
   meas.CollectHistograms();

   // Print to the screen a text representation of the model
   // just for minor debugging
   meas.PrintTree();

   // Finally, run the measurement.
   // This is the same thing that happens when
   // one runs 'hist2workspace' on an xml files
   MakeModelAndMeasurementFast(meas);         


 }

