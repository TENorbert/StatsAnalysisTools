// This is an example of a RooStats interval calculation
// for a counting experiment
//
// ROOT macro approach

/*
RooStats exercise at CMSDAS-2011:
  Fermilab, January 25-29, 2011
  inspired by official RooStats tutorials,
  see http://root.cern.ch/root/html/tutorials/roostats/
*/

#include "TCanvas.h"
#include "TFile.h"
#include "TH1D.h"
#include "RooWorkspace.h"
#include "RooArgList.h"
#include "RooDataHist.h"
#include "RooPlot.h"
#include "RooStats/BayesianCalculator.h"
#include "RooStats/ModelConfig.h"
#include "RooStats/SimpleInterval.h"

#include "RooStats/MCMCCalculator.h"
#include "RooStats/MCMCInterval.h"
#include "RooStats/MCMCIntervalPlot.h"
#include "RooStats/ProposalHelper.h"

#include "RooFitResult.h"

using namespace RooFit;
using namespace RooStats;

void tprime( void ){
  //
  // this function implements the interval calculation for the tprime analysis
  //

  //make workspace
  RooWorkspace * wspace = new RooWorkspace("myWS");
  
  // define variables
  wspace -> factory("x[300,600]");
  wspace -> factory("s[0,0,100]");
  wspace -> factory("b[1000,0,1000]");

  // define observables
  wspace -> defineSet("observables","x");

  // define parameters of interest
  wspace -> defineSet("poi","s");

  // define nuisance parameters
  wspace -> defineSet("nuisance_parameters","b");

 // load data
  TFile * _file = new TFile("PRL_hist.root", "read");
  TH1D * _sig_templ = (TH1D *)_file->Get("topmasspt_tprime_300");
  TH1D * _bg_templ = (TH1D *)_file->Get("topmasspt_top");
  TH1D * _data = (TH1D *)_file->Get("topmasspt_data");
  RooArgList * observables = new RooArgList( *wspace->set("observables") );
  RooDataHist * sig_templ = new RooDataHist("tprime_sig", "tprime_sig", *observables, _sig_templ);
  RooDataHist * bg_templ = new RooDataHist("tprime_bg", "tprime_bg", *observables, _bg_templ);
  RooDataHist * data = new RooDataHist("tprime_data", "tprime_data", *observables, _data);
  sig_templ -> SetName("sig_templ");
  bg_templ -> SetName("bg_templ");
  data -> SetName("data");

  wspace -> import(*data);
  wspace -> import(*sig_templ);
  wspace -> import(*bg_templ);

  wspace -> factory("RooHistPdf::sig_pdf(x,sig_templ)");
  wspace -> factory("RooHistPdf::bg_pdf(x,bg_templ)");

 // make model
  wspace -> factory("SUM::model_pdf(s*sig_pdf,b*bg_pdf)");
  wspace -> factory("Lognormal::likelihood_b(b,1000,1.5)");
  wspace -> factory("PROD::model_likelihood(model_pdf, likelihood_b)");
  wspace -> factory("Uniform::prior_pdf(s)");





 // model config
  ModelConfig modelConfig("model_config");
  modelConfig . SetWorkspace(*wspace);
  modelConfig . SetPdf(*wspace->pdf("model_likelihood"));
  modelConfig . SetPriorPdf(*wspace->pdf("prior_pdf"));
  modelConfig . SetParametersOfInterest(*wspace->set("poi"));
  modelConfig . SetNuisanceParameters(*wspace->set("nuisance_parameters"));
  wspace -> import(modelConfig, "model_config");

  wspace->Print();

  // Bayesian MCMC calculation
  // optional: create a proposal function
  RooFitResult * fit = wspace->pdf("model_likelihood")->fitTo(*data,Save());
  ProposalHelper ph;
  ph.SetVariables((RooArgSet&)fit->floatParsFinal());
  ph.SetCovMatrix(fit->covarianceMatrix());
  ph.SetUpdateProposalParameters(kTRUE); // auto-create mean vars and add mappings
  ph.SetCacheSize(100);
  ProposalFunction* pf = ph.GetProposalFunction();
  
  MCMCCalculator mcmc( *data, modelConfig );
  mcmc.SetConfidenceLevel(0.95);
  mcmc.SetNumIters(10000);          // Metropolis-Hastings algorithm iterations
  mcmc.SetProposalFunction(*pf);
  mcmc.SetNumBurnInSteps(100); // first N steps to be ignored as burn-in
  mcmc.SetLeftSideTailFraction(0.0);
  mcmc.SetNumBins(100);
  
  MCMCInterval * mcInt = mcmc.GetInterval();

  RooRealVar * firstPOI = (RooRealVar*) modelConfig.GetParametersOfInterest()->first();
  double _limit = mcInt->UpperLimit(*firstPOI);
  cout << "\n95% upper limit on " <<firstPOI->GetName()<<" is : "<<
    _limit <<endl;

  // posterior plot
  TCanvas c1("c1");
  MCMCIntervalPlot plot(*mcInt);
  plot.Draw();
  c1.SaveAs("tprime.png");


}
