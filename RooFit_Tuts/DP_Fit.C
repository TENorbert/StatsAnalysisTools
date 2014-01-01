/*//////////// Fitting Tool to Bkg 
/////////////  Orig Author: Tambe E. Norbert
/////////////  @UMN.EDU Go GOPHERS
/////////////  Sept 02, 2013
*/

//Import Needed  Classes
#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooConstVar.h"
#include "RooChebychev.h"
#include "RooAddPdf.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "TRandom.h"
#include "RooSimultaneous.h"
#include "RooCategory.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "RooPlot.h"
#include "TFile.h"
#include "TH2D.h"
#include "TH1D.h"
#include <cmath>
#include <stdio>
#include <vector>
#include <string>

using namespace RooFit ;
using namespace std;



void DP_Fit() {

 // C r e a t e   m o d e l   f o r   p h y s i c s   s a m p l e
  // -------------------------------------------------------------

  // Create observables
  RooRealVar time("time","Time(ns)",-5,25) ;

  // Construct QCD/signal pdf
  RooRealVar mean("mean","mean",0,-5,25) ;
  RooRealVar sigma("sigma","sigma",0.3,0.1,10) ;
  RooGaussian gt("gt","gt",time,mean,sigma) ;

  // Construct HALO/background pdf
  RooRealVar a0("a0","a0",-0.1,-5,25) ;
  RooRealVar a1("a1","a1",0.004,-1,1) ;
  RooChebychev pt("pt","pt",time,RooArgSet(a0,a1)) ;

  // Construct composite pdf
  RooRealVar f("f","f",0.2,0.,1.) ;
  RooAddPdf model("model","model",RooArgList(gt,pt),f) ;

//  TFile* file = IFile();
//  TH1D* ihist = get1Dhist(file);
  TH1D* ihist = makeTH1();   // Random MC

  // Create a binned dataset that imports contents of TH1 and associates its contents to observable 'x'
  RooDataHist dhist("dhist","dhist",time,Import(*ihist));

// P l o t   a n d   f i t   a   R o o D a t a H i s t
  // ---------------------------------------------------

  // Make plot of binned dataset showing Poisson error bars (RooFit default)
  RooPlot* frame = time.frame(Title("Imported Time Hist with Poisson error bars")) ;
  dhist.plotOn(frame) ; 

  // Fit a Gaussian p.d.f to the data
  gt.fitTo(dhist) ;
  gt.plotOn(frame) ;

 // If histogram has custom error (i.e. its contents is does not originate from a Poisson process
  // but e.g. is a sum of weighted events) you can data with symmetric 'sum-of-weights' error instead
  // (same error bars as shown by ROOT)
  RooPlot* frame2 = time.frame(Title("Imported Time TH1 with internal errors")) ;
  dhist.plotOn(frame2,DataError(RooAbsData::SumW2)) ; 
  gt.plotOn(frame2) ;



 // Draw all frames on a canvas
  TCanvas* c = new TCanvas("DP_Fit","DP_Fit",800,800) ;
  c->Divide(2,0) ;
  c->cd(1) ; gPad->SetLeftMargin(0.15) ; frame->GetYaxis()->SetTitleOffset(1.4) ; frame->Draw() ;
  c->cd(2) ; gPad->SetLeftMargin(0.15) ; frame2->GetYaxis()->SetTitleOffset(1.4) ; frame2->Draw() ;
/*  c->cd(3) ; gPad->SetLeftMargin(0.15) ; frame3->GetYaxis()->SetTitleOffset(1.4) ; frame3->Draw() ;
  c->cd(4) ; gPad->SetLeftMargin(0.15) ; frame4->GetYaxis()->SetTitleOffset(1.4) ; frame4->Draw() ;
*/
}

// Get input ROOT file
TFile *IFile(){

TFile* ifile = new TFile("myinputfile.root")
cout <<"Reading file ... " << ifile << endl;
return ifile;
}

//Get histogram
TH1D* get1DHist(TFile * ifile) {

TH1D* h_d = (TH1D*) ifile->Get("Rootfilename/Cur_dir/myhist");
TH1D* hclone = h_d->Clone();
return hclone;
}

TH1* makeTH1() {
  // Create ROOT TH1 filled with a Gaussian distribution

  TH1D* hh = new TH1D("hh","hh",30,-5,25) ;
  for (int i=0 ; i<100 ; i++) {
    hh->Fill(gRandom->Gaus(0,3)) ;
  }
  return hh ;
}



