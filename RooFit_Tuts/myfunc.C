// Simple function to fit
//Aept 09, 2013
// Author 10Sr@umn.edu


#include "Riostream.h"
#include "TF1.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TMath.h"

using namespace std;

void myfunc()
{
  TF1*fx = new TF1("fx","x*cos(x)", -5.0, 25.);  
  TF1*fx1 = new TF1("fx1","x*sin(x)", -5.0, 25.);  
  TF1*fx2 = new TF1("fx2","TMath::BreitWigner(x,0,1) + TMath::Exp(-0.5*x)", -5.0, 25.);  
  TF1*fx3 = new TF1("fx3","TMath::Gaus(x,0,0.5) + TMath::Exp(-0.5*x )", -5.0, 25.0);  
TCanvas * c1 = new TCanvas("c1","c1", 800., 800.);
/*TCanvas * c2 = new TCanvas("c2","c2", 800., 800.);
TCanvas * c3 = new TCanvas("c3","c3", 800., 800.);
TCanvas * c4 = new TCanvas("c4","c4", 800., 800.);
*/
std::cout <<"begin drawing..."<< std::endl;
c1->Divide(2,2,0,0);
c1->cd(1);
fx->SetLineColor(kGreen);
fx->Draw();
c1->cd(2);
fx1->SetLineColor(kBlue);
fx1->Draw();
c1->cd(3);
fx2->SetLineColor(kRed);
fx2->Draw();
c1->cd(4);
fx3->SetLineColor(kCyan);
fx3->Draw();

}


