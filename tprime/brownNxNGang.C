#include <iostream>
#include <fstream>
#include <sstream>

#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"

//#include "combineBins_ejets.C"
#include "combineBins.C"
string convertInt(int number)
{
  stringstream ss;//create a stringstream                                                         
  ss << number;//add number to the stream                                                         
  return ss.str();//return a string with the contents of the stream                               
}

void 
rebin( int mass=450, 
       double scale=1.,
       double maxErr=0.1,
       std::string inFileName="",
       std::string outFileName="",
       bool make_plots = false ){
  //
  // Mike4x4-style bin ganging
  //


  std::cout << "Doing 2D rebinning to 1D, Mike-style\n"
	    << "mass = " << mass << std::endl
	    << "scale factor (optional) = " << scale << std::endl
	    << "max uncertainty = " << maxErr << std::endl
	    << "input file with 2D templates: " << inFileName << std::endl
	    << "output file: " << outFileName << std::endl;


  std::string channel = "ejets";
  int rebin_x = 4;
  int rebin_y = 4;
  
  TString inFileN = inFileName;
  TFile *infile = new TFile(inFileN);
  const int nhists = 10;
  TString histonames[nhists]={ //array of data and background histograms        
    "tprime"+convertInt(mass),
    "DATA", // data histogram must be first in this list                                       
    "top",
    "ewk",
    "tprime"+convertInt(mass)+"__jes__plus",
    "tprime"+convertInt(mass)+"__jes__minus",
    "top__jes__plus",
    "top__jes__minus",
    "ewk__jes__plus",
    "ewk__jes__minus"
  };
  TFile *outfile = new TFile("temp_rebinfile.root","RECREATE");
  TH2F *hists[nhists];
  for(int nh = 0 ; nh < nhists; ++nh){
    hists[nh] = (TH2F*)infile->Get(histonames[nh]);
    hists[nh] -> Rebin2D(rebin_x,rebin_y);
    hists[nh] -> Write();
  }
  outfile->Write();
  outfile->Close();
  combineBins(channel, mass, scale, maxErr,"temp_rebinfile.root", outFileName, "", make_plots);
  return;
}
