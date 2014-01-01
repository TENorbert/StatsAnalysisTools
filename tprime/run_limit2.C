#include "TSystem.h"
#include "hf_tprime.C"
void run_limit2( std::string channel, // ejets, mujets, combined
		std::string mode,    // observed, expected
		std::string method,  // mcmc, plr, cls
		double peak,         // resonance mass
		std::string suffix,  // suffix for output file names
		Int_t ntoys,         // number of pseudoexperiments for expected limit
		Int_t npoints,       // number points for CLs scan
		Double_t poimin,     // lower scan range
		Double_t poimax,     // upper scan range
		Int_t mcmc_iter,     // number of MCMC iterations
		Int_t mcmc_burnin,   // number of MCMC burn in steps to be discarded
		std::string inputdir) // directory with workspace files
{

  
  gSystem->SetIncludePath( "-I./");
  gSystem->Load("hf_tprime_C.so");
  limit( channel,
	 mode,
	 method,
	 peak,
	 suffix,
	 ntoys,
	 npoints,
	 poimin,
	 poimax,
	 mcmc_iter,
	 mcmc_burnin,
	 inputdir );
  return;
}
