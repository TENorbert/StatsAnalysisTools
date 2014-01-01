//////////////////////////////////////////////////////////////////////////
//
// tprime.C  version 0.1
//
// Routine for statistical inference for a tprime analysis
//
//
// Gena Kukartsev
//
// May 2011: first version
// 

/*
The code should be compiled in ROOT (5.28.00c or higher):

root -l

.L hf_tprime.C+

Usage:
       limit = limit(param, ...);

Inputs:
       param          - Parameter


Description goes here
*/

#include <iostream>
#include <fstream>

#include "TFile.h"
#include "TCanvas.h"
#include "TRandom3.h"
#include "TUnixSystem.h"
#include "TH1D.h"

#include "RooRandom.h"
#include "RooWorkspace.h"
#include "RooArgList.h"
#include "RooDataHist.h"
#include "RooPlot.h"
#include "RooFitResult.h"

#include "RooStats/BayesianCalculator.h"
#include "RooStats/ModelConfig.h"
#include "RooStats/SimpleInterval.h"
#include "RooStats/MCMCCalculator.h"
#include "RooStats/MCMCInterval.h"
#include "RooStats/MCMCIntervalPlot.h"
#include "RooStats/ProposalHelper.h"
#include "RooStats/HistFactory/LinInterpVar.h"

using namespace RooFit;
using namespace RooStats;

// ---> user interface ---------------------------------------------

void tprime_prototype( void );
void hf_tprime(void){} // dummy

// ---> implementation below ---------------------------------------

class Tprime{
  //
  // The class combines multiple channel analyses. A workspace is created
  // with the combined model, data, model config. The class can also call
  // interval calculation routines
  //

public:

  Tprime();
  ~Tprime();
  
  Int_t LoadWorkspace(std::string filename,
		      std::string ws_name);

  RooStats::ModelConfig * GetModelConfig( void );

  MCMCInterval * GetMcmcInterval(double conf_level,
				 int n_iter,
				 int n_burn,
				 double left_side_tail_fraction,
				 int n_bins);

  MCMCInterval * GetMcmcInterval(ModelConfig mc,
				 double conf_level,
				 int n_iter,
				 int n_burn,
				 double left_side_tail_fraction,
				 int n_bins);

  void makeMcmcPosteriorPlot( std::string filename );
  double printMcmcUpperLimit( double peak, std::string filename = "" );

  Double_t GetRandom( std::string pdf, std::string var );

  void ExamineTemplateFile(std::string filename);

private:

  RooWorkspace * ws;

  // roostats calculators results
  MCMCInterval * mcInt;

  // random numbers
  TRandom3 r;
};


Tprime::Tprime(){
  //
  // Default constructor.
  // Initialize data members, Set random seed.
  //

  std::string legend = "[Tprime::Tprime]: ";
  ws = 0;
  mcInt = 0;

  // set random seed
  r.SetSeed();
  UInt_t _seed = r.GetSeed();
  UInt_t _pid = gSystem->GetPid();
  std::cout << legend << "Random seed: " << _seed << std::endl;
  std::cout << legend << "PID: " << _pid << std::endl;
  _seed = 31*_seed+_pid;
  std::cout << legend << "New random seed (31*seed+pid): " << _seed << std::endl;
  r.SetSeed(_seed);

  // set RooFit random seed (it has a private copy)
  RooRandom::randomGenerator()->SetSeed(_seed);
}


Tprime::~Tprime(){
  //
  // Destructor. Clean up.
  //

  delete ws;
  delete mcInt;
}


Int_t Tprime::LoadWorkspace(std::string filename,
			    std::string ws_name){
  //
  // load a workspace from a file
  // a copy of the workspace is kept as a class member,
  // the input file is immediately closed
  //
  
  TFile * infile = new TFile(filename.c_str(), "read");
  ws = (RooWorkspace *)infile->Get("tprime")->Clone();
  delete infile;
  
  ws->Print();

  return 0;
}


RooStats::ModelConfig * Tprime::GetModelConfig( void ){
  //
  // Return a pointer to the ModelConfig in the loaded
  // workspace, or 0 if no workspace is loaded.
  // User does NOT take ownership.
  //

  if (ws){
    RooStats::ModelConfig * _mc = (RooStats::ModelConfig *)ws->genobj("ModelConfig");
    _mc -> SetWorkspace(*ws);
    //_mc->Print();
    //_mc->GetWorkspace()->Print();
    return _mc;
  }
  else return 0;
}


MCMCInterval * Tprime::GetMcmcInterval(double conf_level,
				       int n_iter,
				       int n_burn,
				       double left_side_tail_fraction,
				       int n_bins){
  //
  // Bayesian MCMC calculation using the ModelConfig from the workspace
  // Want an efficient proposal function, so derive it from covariance
  // matrix of fit
  //

  return GetMcmcInterval(*GetModelConfig(),
			 0.95,
			 10000,
			 100,
			 0,
			 100);
}

  
MCMCInterval * Tprime::GetMcmcInterval(ModelConfig mc,
				       double conf_level,
				       int n_iter,
				       int n_burn,
				       double left_side_tail_fraction,
				       int n_bins){
  //
  // Bayesian MCMC calculation using arbitrary ModelConfig
  // Want an efficient proposal function, so derive it from covariance
  // matrix of fit
  //
  
  RooAbsData * _data = ws->data("obsData");
  //RooStats::ModelConfig * _mc = (RooStats::ModelConfig *)ws->genobj("ModelConfig");
  

  RooFitResult * fit = ws->pdf("model_tprime")->fitTo(*_data,Save());
  ProposalHelper ph;
  ph.SetVariables((RooArgSet&)fit->floatParsFinal());
  ph.SetCovMatrix(fit->covarianceMatrix());
  ph.SetUpdateProposalParameters(kTRUE); // auto-create mean vars and add mappings
  ph.SetCacheSize(100);
  ProposalFunction * pf = ph.GetProposalFunction();
  
  MCMCCalculator mcmc( *_data, mc );
  mcmc.SetConfidenceLevel(conf_level);
  mcmc.SetNumIters(n_iter);          // Metropolis-Hastings algorithm iterations
  mcmc.SetProposalFunction(*pf);
  mcmc.SetNumBurnInSteps(n_burn); // first N steps to be ignored as burn-in
  mcmc.SetLeftSideTailFraction(left_side_tail_fraction);
  mcmc.SetNumBins(n_bins);
  
  mcInt = mcmc.GetInterval();
  
  return mcInt;
}


void Tprime::makeMcmcPosteriorPlot( std::string filename ){
  //
  // Produce an MCMC posterior plot from the MCMCInterval,
  // and save to file.
  //

  TCanvas c1("c1");
  MCMCIntervalPlot plot(*mcInt);
  plot.Draw();
  c1.SaveAs(filename.c_str());
  
  return;
}



double Tprime::printMcmcUpperLimit( double peak, std::string filename ){
  //
  // print out the upper limit on the first Parameter of Interest
  //

  RooStats::ModelConfig * _mc = (RooStats::ModelConfig *)ws->genobj("ModelConfig");

  RooRealVar * firstPOI = (RooRealVar*) _mc->GetParametersOfInterest()->first();
  double _limit = mcInt->UpperLimit(*firstPOI);
  cout << "\n95% upper limit on " <<firstPOI->GetName()<<" is : "<<
    _limit <<endl;

  if (filename.size()!=0){
    
    std::ofstream aFile;

    // append to file if exists
    aFile.open(filename.c_str(), std::ios_base::app);

    char buf[1024];
    sprintf(buf, "%7.1f   %7.6f", peak, _limit);

    aFile << buf << std::endl;

    // close outfile here so it is safe even if subsequent iterations crash
    aFile.close();

  }

  return _limit;
}


Double_t Tprime::GetRandom( std::string pdf, std::string var ){
  //
  // generates a random number using a pdf in the workspace
  //
  
  // generate a dataset with one entry
  RooDataSet * _ds = ws->pdf(pdf.c_str())->generate(*ws->var(var.c_str()), 1);

  Double_t _result = ((RooRealVar *)(_ds->get(0)->first()))->getVal();
  delete _ds;

  return _result;
}




void tprime_prototype( void ){
  //
  // this function implements the interval calculation for the tprime analysis
  //
 TFile* file =new TFile("results/tprime_tprime_tprimeCrossSection_model.root");
  RooWorkspace* WSpace = (RooWorkspace*) file->Get("tprime");
  RooAbsData * _data = WSpace->data("obsData");
  RooStats::ModelConfig * _mc = (RooStats::ModelConfig *)WSpace->genobj("ModelConfig");
  TFile *histfile = new TFile("hist300.root","RECREATE");
 

  TH1D* mcslimit= new TH1D("mcslimit","mcslimit",500,0.,50.);
  //Make a histogram of the upperlimit vs the number of generated events
 
  //TH1* mcslimitvsevt_histo=(TH1F*)mcs->fitParDataSet().createHistogram("hist2");
 
  for(int i=0; i<10; i++){
  Tprime tp;

  //load workspace
  tp.LoadWorkspace("results/tprime_tprime_tprimeCrossSection_model.root", "tprime");

  tp.GetMcmcInterval(0.95,
		     10000,
		     100,
		     0,
		     100);

  std::string suffix = "prototype";

  std::string _outfile = "comb_" + suffix + ".ascii";
  tp.printMcmcUpperLimit( 300.0, _outfile );

  _outfile = "comb_mcmc_posterior" + suffix + ".pdf";
  tp.makeMcmcPosteriorPlot( _outfile );

<<<<<<< hf_tprime.C
=======
  RooRealVar * firstPOI = (RooRealVar*) _mc->GetParametersOfInterest()->first();
   MCMCCalculator mcmc( *_data, *tp.GetModelConfig() );
    mcmc.SetConfidenceLevel(0.95);
  mcmc.SetNumIters(10000);          // Metropolis-Hastings algorithm iterations
  // mcmc.SetProposalFunction(*pf);
  mcmc.SetNumBurnInSteps(100); // first N steps to be ignored as burn-in
  mcmc.SetLeftSideTailFraction(0.0);
  mcmc.SetNumBins(100);
   MCMCInterval * mcInt = mcmc.GetInterval();
  double _limit = mcInt->UpperLimit(*firstPOI);
  mcslimit->Fill(_limit);
  }
  histfile->Write();
  histfile->Close();
  
  // Bayesian MCMC calculation
  // optional: create a proposal function
  /*
  TFile* file =new TFile("results/tprime_tprime_tprimeCrossSection_model.root");
  RooWorkspace* WSpace = (RooWorkspace*) file->Get("tprime");
  RooAbsData * _data = WSpace->data("obsData");
  RooStats::ModelConfig * _mc = (RooStats::ModelConfig *)WSpace->genobj("ModelConfig");
  
  
  RooFitResult * fit = WSpace->pdf("model_likelihood")->fitTo(*_data,Save());
  ProposalHelper ph;
  ph.SetVariables((RooArgSet&)fit->floatParsFinal());
  ph.SetCovMatrix(fit->covarianceMatrix());
  
  ph.SetUpdateProposalParameters(kTRUE); // auto-create mean vars and add mappings
  ph.SetCacheSize(100);
  ProposalFunction* pf = ph.GetProposalFunction();
  
  MCMCCalculator mcmc( *_data, *tp.GetModelConfig() );
    mcmc.SetConfidenceLevel(0.95);
  mcmc.SetNumIters(10000);          // Metropolis-Hastings algorithm iterations
  mcmc.SetProposalFunction(*pf);
  mcmc.SetNumBurnInSteps(100); // first N steps to be ignored as burn-in
  mcmc.SetLeftSideTailFraction(0.0);
  mcmc.SetNumBins(100);
  
  MCMCInterval * mcInt = mcmc.GetInterval();

  RooRealVar * firstPOI = (RooRealVar*) _mc->GetParametersOfInterest()->first();
  double _limit = mcInt->UpperLimit(*firstPOI);
  cout << "\n95% upper limit on " <<firstPOI->GetName()<<" is : "<<
    _limit <<endl;

  // posterior plot
  TCanvas c1("c1");
  MCMCIntervalPlot plot(*mcInt);
  plot.Draw();
  c1.SaveAs("hf_tprime_mcmc.png");
  */

>>>>>>> 1.2
}
