#include <iostream>
#include <fstream>
#include <stdio.h>
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"

#include "TCanvas.h"
#include "TMath.h"
#include "TF1.h"
#include "TF2.h"
#include "TROOT.h"
#include "TRandom.h"
//
// Recent changes:
//
// - last unrebinned bin now added to the next-to-last
// - signal precision included together with total background
// - s=0 / b=0 is now  99999 (was 0)
// - s>0 / b=0 is now 999999 (was 999)

// forward declaration
void combineBins( std::string channel,
		  int mass,
		  double scale_factor,
		  double femax,              // bin precision for background
		  std::string inFileName,
		  std::string outFileName,
		  std::string prefix = "",
		  bool make_plots = false,
		  double femax_sig = 0.20 ); // separate bin precision for signal



// sig mass and ht PDF
Double_t gaussian_mass(Double_t *x, Double_t *par) {
  return par[0]*TMath::Gaus(x[0], par[1], par[2])
    + par[3]*TMath::Gaus(x[0], par[4], par[5]);
}

//Double_t gaussian_ht(Double_t *x, Double_t *par) {
//  return par[0]*TMath::Gaus(x[0], par[1], par[2])
//    + par[3]*TMath::Gaus(x[0], par[4], par[5]);
//}
Double_t gaussian_ht(Double_t *x, Double_t *par) {
  return par[0]*(par[1]*TMath::Gaus(x[0], par[2], par[3])
		 +(1.0-par[1])*TMath::Gaus(x[0], par[4], par[5]));
}

// bg mass PDF
Double_t bgMassPdf(Double_t *x, Double_t *par) {
  return par[0]*TMath::LogNormal(x[0]/200.0, par[1], par[2], par[3])
    + par[4]*TMath::Gaus(x[0], par[5], par[6]);
}

// bg HT PDF
Double_t bgHtPdf(Double_t *x, Double_t *par) {
  return par[0]*TMath::Gaus(x[0], par[1], par[2])
    + par[3]*TMath::Gaus(x[0], par[4], par[5]);
  //+ par[3]*TMath::Landau(x[0], par[4], par[5]);
}

Double_t func2D_sig(Double_t *x,Double_t *par)
{
  Double_t p0=par[0]*(x[0]-par[1])*(x[0]-par[1])*exp(-par[2]*(x[0]-par[1]))+par[3]*exp(-0.5*(x[0]-par[4])*(x[0]-par[4])/par[5]/par[5])+par[6]*exp(-0.5*(x[0]-par[7])*(x[0]-par[7])/par[8]/par[8]);

  Double_t p1=par[9]+par[10]*x[0];
  Double_t p2=par[11]+par[12]*x[0];
  Double_t p3=par[13]+par[14]*x[0];
  Double_t p4=par[15]+par[16]*x[0];
  Double_t f = par[17]*p0*((1-p4)*TMath::Gaus(x[1],p1,p2)/p2+p4*TMath::Gaus(x[1],1.3*p1,p3)/p3);
  //  cout << x[0]<<","<<x[1]<<" p= "<<p0<<", "<<p1<<", "<<p2<<", "<<p3<<", "<<p4<<" f= "<<f<<endl;
  return f;
}

Double_t func2D(Double_t *x,Double_t *par)
{
  Double_t p0=par[0]*(x[0]-par[1])*(x[0]-par[1])*exp(-par[2]*(x[0]-par[1]))
    +par[3]*TMath::Gaus(x[0],par[4],par[5]);

  Double_t p1=par[9]+par[10]*x[0];
  Double_t p2=par[11]+par[12]*x[0];
  Double_t p3=par[13]+par[14]*x[0];
  Double_t p4=par[15]+par[16]*x[0];

  Double_t f = par[17]*p0*((1-p4)*TMath::Gaus(x[1],p1,p2)/p2+p4*TMath::Gaus(x[1],1.3*p1,p3)/p3);

  //  cout << x[0]<<","<<x[1]<<" p= "<<p0<<", "<<p1<<", "<<p2<<", "<<p3<<", "<<p4<<" f= "<<f<<endl;                                                                                                 
  return f;
}

void combineBins(int mass, double scale_factor = 1.0){
  //
  // Wrapper method for backwards compatibility
  //
  //   mass = mass of tprime quark

  // input file name
  //  std::string fname = "data/tprime_templates_27dec11/tprime_ejets_4683ipb_2d_v1_27dec11.root";
    std::string fname = "ejets_reb20jan11/input_merged_e.root";
  // output file name
  //  char _template[256] = "data/ejets_3560/ht35/tprime_%i_ejets_3560ipb_ht35_2d_merged_14nov2011_001.root";
  char _template[256] = "tprime_%i_ejets_merged_01feb12.root";
  //char _template[256] = "output.root";
  char buf[256];
  sprintf(buf, _template, mass);
  std::string oname(buf);

  // highest tolerated stat uncertainty in any bin
  float femax=0.20; //max fractional error in each bin of background histogram

  combineBins("ejets", mass, scale_factor, femax, fname, oname);

  return;
}



void 
combineBins( std::string channel, int mass, double scale_factor,
	     double femax,
	     std::string inFileName,
	     std::string outFileName,
	     std::string prefix,
	     bool make_plots,
	     double femax_sig ){

  char * fname = (char *)(inFileName.c_str());

  char * oname = (char *)(outFileName.c_str());

  char sname[100]; //name of signal histogram
  sprintf(sname,"tprime%i",mass);

  std::vector<std::string> vBname;
  vBname.push_back("DATA"); // data histogram must be first in this list
  vBname.push_back("top");
  vBname.push_back("ewk");
  vBname.push_back("tprime%i__jes__plus");
  vBname.push_back("tprime%i__jes__minus");
  vBname.push_back("top__jes__plus");
  vBname.push_back("top__jes__minus");
  vBname.push_back("ewk__jes__plus");
  vBname.push_back("ewk__jes__minus");
    vBname.push_back("tprime%i__btgsf__plus");
  vBname.push_back("tprime%i__btgsf__minus");
  vBname.push_back("top__btgsf__plus");
  vBname.push_back("top__btgsf__minus");
  vBname.push_back("ewk__btgsf__plus");
  vBname.push_back("ewk__btgsf__minus");
  vBname.push_back("tprime%i__jer__plus");
  vBname.push_back("tprime%i__jer__minus");
  vBname.push_back("top__jer__plus");
  vBname.push_back("top__jer__minus");
  vBname.push_back("ewk__jer__plus");
  vBname.push_back("ewk__jer__minus");
  

  char bname[30][100];//={ //array of data and background histograms
//    "DATA", // data histogram must be first in this list
//    "top",
//    "ewk",
//    "tprime%i__jes__plus",
//    "tprime%i__jes__minus",
//    "top__jes__plus",
//    "top__jes__minus",
//    "ewk__jes__plus",
//    "ewk__jes__minus",
//    "tprime%i__btgsf__plus",
//    "tprime%i__btgsf__minus",
//    "top__btgsf__plus",
//    "top__btgsf__minus",
//    "ewk__btgsf__plus",
//    "ewk__btgsf__minus",
//    "tprime%i__jer__plus",
//    "tprime%i__jer__minus",
//    "top__jer__plus",
//    "top__jer__minus",
//    "ewk__jer__plus",
//    "ewk__jer__minus"
//  };

  int nxx = 1,nyy = 1;//rebinning before processing 
  int nb=21; //number of histograms in list
  int n_skip=3; // starting with this index, do not consider for background normalization
  ////////float femax=0.20; //max fractional error in each bin of background histogram

  // bin precision for signal and background
  std::vector<double> vFemax;
  vFemax.push_back(femax_sig);
  vFemax.push_back(femax);

  TFile *f = TFile::Open(fname);
  if (f==NULL) {
    printf("Cannot open file '%s'\n",fname);
    return;
  }

  // prune non-existent hists
  int nb_new = 0;
  for (std::vector<std::string>::const_iterator _n = vBname.begin();
       _n != vBname.end();
       ++_n){
    sprintf(bname[nb_new],_n->c_str(), mass);      
    if ( f->Get( bname[nb_new] ) ){
      // add to list
      ++nb_new;
    }
  }
  nb = nb_new;
  for (int i=0; i<nb; ++i){
    std::cout << bname[i] << std::endl;
  }
  std::cout << "Total: " << nb << " valid histograms" << std::endl;

  //TH2D* hs; f->GetObject(sname,hs); 
  TH2D * hs = (TH2D *)f->Get(sname); 
  double sEntries = hs->GetEntries(); // number of entries in signal

  if (hs==NULL) {
    printf("Cannot find histogram '%s' in '%s'\n",sname,fname);
    return;
  }
  
  // Gena: scale signal template to proper cross section
  hs->Scale(scale_factor);
  hs->Rebin2D(nxx,nyy);

  //figure out the binning
  int nx = hs->GetNbinsX()+2;
  int ny = hs->GetNbinsY()+2;

  // cross check printout
  std::cout << "2D hist name: " << hs->GetName() << std::endl;
  std::cout << "Integral with overflow: " << hs->Integral(0,nx-1,0,ny-1) << std::endl;
  std::cout << "Integral no overflow: " << hs->Integral(1,nx-2,1,ny-2) << std::endl << std::endl;

  TH2D *hb = (TH2D*)hs->Clone(); 
  hb->SetName("hb");
  hb->Reset();
  TH2D *hX[30];
  double bEntries[30]; // array with number of entries
  for (int i=0;i<nb;i++){
    std::string sBName(bname[i]);
    // GENA: get names for signal JES histos
    if (sBName.find("TPrime")!=std::string::npos ||
	sBName.find("Tprime")!=std::string::npos ||
	sBName.find("tprime")!=std::string::npos){

      sprintf(bname[i],sBName.c_str(),mass);
      std::cout << bname[i] << std::endl;
    }

    hX[i] = 0;
    //f->GetObject(bname[i],hX[i]); 
    hX[i] = (TH2D *)f->Get(bname[i]); 
    hX[i]->Rebin2D(nxx,nyy);
    //if (!hX[i]){
    //  std::cout << "No histogram " << bname[i]
    //		<< ", skipping and reducing total number of histograms" 
    //		<< std::endl;
    //  --nb;
    //  continue;
    //}
    bEntries[i] = hX[i]->GetEntries();
    
    // GENA: scale JES signal templates to proper cross section
    if (sBName.find("TPrime")!=std::string::npos ||
	sBName.find("Tprime")!=std::string::npos ||
	sBName.find("tprime")!=std::string::npos){
      hX[i]->Scale(scale_factor);      
    }
    if (hX[i]==NULL) {
      printf("Cannot find histogram '%s' in '%s'\n",bname[i],fname);
      return;
    }
    //hX[i]->Print("base");
    std::cout << "2D hist name: " << hX[i]->GetName() << std::endl;
    std::cout << "Integral with overflow: " << hX[i]->Integral(0,nx-1,0,ny-1) << std::endl;
    std::cout << "Integral no overflow: " << hX[i]->Integral(1,nx-2,1,ny-2) << std::endl << std::endl;
    //sum all background histograms into hb; do not add the data histogram
    if (i>0 && i<n_skip) hb->Add(hX[i]); 
  }

  //figure out the binning
  //int nx = hs->GetNbinsX();
  //int ny = hs->GetNbinsY();
  int nbin=nx*ny;
  std::cout << "number of bins: x="<<nx<<", y="<<ny<<std::endl;
  
  //book some 1d histograms with the same number of bins for diagnostics
  TH1F *h1sb = new TH1F("h1sb","h1sb",nbin,0,nbin);
  TH1F *h1s = new TH1F("h1s","h1s",nbin,0,nbin);
  TH1F *h1b = new TH1F("h1b","h1b",nbin,0,nbin);
  // GENA: vector to create 2D->1D bin map
  std::vector<std::vector<int> > vMap(nbin);

  //
  // -------------- fits
  // GENA: SB:
  // do simple 1d fits

  //MIKE : 2D:
  // Create a new canvas.                                                                         
  TCanvas *cx = new TCanvas("cx","Mfit projection",10,10,700,500);
  cx->SetFillColor(0);
  cx->SetFrameFillColor(0);
  TCanvas *cy = new TCanvas("cy","HT projection",10,10,700,500);
  cy->SetFillColor(0);
  cy->SetFrameFillColor(0);

  cx->cd();
  TH1D *hx = hs->ProjectionX("hx",0,ny);
  hx->SetFillColor(38);
  char title[80];
  sprintf(title,"Mfit Projection");
  hx->SetName("hx");
  hx->SetTitle(title);

  //gaussian + x^2*exp                                                                            
  TF1 *fx = new TF1("fx","[0]*(x-[1])*(x-[1])*exp(-[2]*(x-[1]))+[3]*exp(-0.5*(x-[4])*(x-[4])/[5]/[5])+[6]*exp(-0.5*(x-[7])*(x-[7])/[8]/[8])",100,1000);
   //fx->SetParameters(100,450,100,100,400,100,100,300,300);                                       
   fx->SetParameters(0.002,120,0.01,25,450,30,10,350,50);
   hx->Fit("fx","l");
   hx->GetFunction("fx")->SetLineColor(kRed);
   hx->GetFunction("fx")->SetLineWidth(6);
   cx->Update();

   //save the parameters in p_mfit                                                                
   Double_t p_mfit[9];
   TF1 *fit = hx->GetFunction("fx");
   for (int ipar=0;ipar<9;ipar++){
     p_mfit[ipar] = fit->GetParameter(ipar);
     cout << p_mfit[ipar]<<endl;
   }

   //fit the HT projection (y) - for kicks only, this is not used again                          
   cy->cd();
   TH1D *hy = hs->ProjectionY("hy",0,nx);
   hy->SetFillColor(38);
   //   char titlebg[80];
   sprintf(title,"HT Projection");
   hy->SetName("hy");
   hy->SetTitle(title);
   TF1 *fy = new TF1("fy","[0]*[0]*exp(-0.5*(x-[1])*(x-[1])/[2]/[2])+[3]*[3]*exp(-0.5*(x-[4])*(x-[4])/[5]/[5])+[6]*[6]*exp(-0.5*(x-[7])*(x-[7])/[8]/[8])+[9]*[9]*exp(-0.5*(x-[10])*(x-[10])/[11]/[11\
])",100,1000);
   fy->SetParameters(20,466,10,4,860,245,14,628,100,60,1000);
   fy->SetParameter(11,700);
   hy->Fit("fy","l");
   hy->GetFunction("fy")->SetLineColor(kRed);
   hy->GetFunction("fy")->SetLineWidth(6);
   cy->Update();

   //fit HT slices with double Gaussian                                                            

   //book histograms of parameter values vs Mfit                                                   
   TH1D *hpar[5];
   char buf2[16];
   for (int ipar=0;ipar<5;ipar++){
     sprintf(buf2, "hpar_%d", ipar);
     //hpar[ipar] = hs->ProjectionX("",1,1);
     hpar[ipar] = hs->ProjectionX(buf2,1,1);
     hpar[ipar]->Reset();
     char hname[20];
     sprintf(hname,"Par%i",ipar);
     hpar[ipar]->SetName(hname);
   }

   TCanvas *cyy = (TCanvas*)gROOT->GetListOfCanvases()->FindObject("cyy");
   if(cyy) delete cyy->GetPrimitive("Projection");
   else   cyy = new TCanvas("cyy","Projection Canvas",700,10,490,350);
   cyy->SetGrid();
   cyy->cd();

   //double gaussian                                                                            
   TF1 *f1 = new TF1("f1","[0]*((1-[4])*exp(-0.5*(x-[1])*(x-[1])/[2]/[2])/[2]+[4]*exp(-0.5*(x-1.3*[1])*(x-1.3*[1])/[3]/[3])/[3])",0,1000);
   //gROOT->ProcessLine(".L function.C");                                                         
   //TF1 *f1 = new TF1("f1",function,0,1000,5);                                                   
   cout << "number of x bins: " << nx <<endl;

   //loop over the HT slices                                                                     
   for (int binx=1;binx<nx+1;binx++){
     TH1D *hp = hs->ProjectionY("hp",binx,binx);
     if (hp->GetEntries()>50){
       hp->Scale(1/hp->Integral());
       hp->SetFillColor(38);
       char title[80];
       sprintf(title,"Slice of binx=%d",binx);
       hp->SetName("Slice");
       hp->SetTitle(title);

       //if the error in a bin is tiny make it bigger so these bins don't drive the fit         
       for (int ibin=1;ibin<ny+1;ibin++){
	 float ebin=hp->GetBinError(ibin);
	 if (ebin>0 && ebin<0.01) hp->SetBinError(ibin,0.01);
	 float cbin=hp->GetBinContent(ibin);
	 if (cbin>0) cout <<"bin "<<ibin<<" = "<<cbin<<" +/- "<<ebin<<endl;
       }

       //first fit a gaussian to get starting values for the parameters                          
       hp->Fit("gaus","l0");
       TF1 *fit = hp->GetFunction("gaus");
       Double_t p0 = fit->GetParameter(0);
       Double_t p1 = fit->GetParameter(1);
       Double_t p2 = fit->GetParameter(2);

       cyy->Update();
       f1->SetParameters(p0,0.9*p1,p2,1.5*p2,0.25);
       f1->SetParLimits(4,0,0.5);
       f1->SetParLimits(3,70,500);
       f1->SetParLimits(2,20,1000);
       cout << "parameters: "<<p0<<","<<p1<<","<<p2<<endl;

       //     int q;                                                                           
       //     cout <<"continue?"<<endl;                                                        
       //     cin>>q;                                                              
       //then fit the real thing                                                                     
       int Q=1;
       int nint=0;
       do {
	 hp->Fit("f1","lm");
	 Double_t chi2 = f1->GetChisquare();
	 if (chi2>200 && nint<5){
	   float x1=gRandom->Gaus(1,0.2);
	   float x2=gRandom->Gaus(2,1);
	   f1->SetParameters(p0,p1,x1*p2,x2*p2,0.25);
	   nint++;
	 }else{
	   Q=0;
	 }
       }while (Q==1);
       hp->GetFunction("f1")->SetLineColor(kRed);
       hp->GetFunction("f1")->SetLineWidth(6);
       cyy->Update();

       //make a png file of the canvas                                                           
       if (make_plots){  
	 char pname[80];
	 sprintf(pname,"hist%i.png",binx);
	 cyy->Print(pname);
       }

       //stop so one can look at the fit                                                         
       //int q;                                                                                 
       //cout <<"continue?"<<endl;                                                            
       //cin>>q;                                                                                 
       // write the fit results into the hpar histograms                                             
       Double_t par,err;
       for (int ipar=0;ipar<5;ipar++){
	 par=f1->GetParameter(ipar);
	 err=f1->GetParError(ipar);
	 if(err>999)err=0;
	 hpar[ipar]->SetBinContent(binx,par);
	 hpar[ipar]->SetBinError(binx,err);
       }
       hp->Delete();
     }
   }
   //end of loop over HT slice fits 
   //now fit the parameters vs Mfit with a line                                                    
   Double_t p_ht[5][2];
   for (int ipar=0;ipar<5;ipar++){
     hpar[ipar]->Fit("pol1","FW");
     TF1 *fit = hpar[ipar]->GetFunction("pol1");
     Double_t p0 = fit->GetParameter(0);
     Double_t p1 = fit->GetParameter(1);

     TF1 *fpol1=new TF1("fpol1","[0]+[1]*x+[2]*exp(-0.5*(x-[3])*(x-[3])/[4])",0,1000);
     fpol1->SetParameters(p0,p1,10,500);
     fpol1->FixParameter(4,2500);
     hpar[ipar]->Fit("fpol1","FW");
     fpol1->FixParameter(4,400);
     hpar[ipar]->Fit("fpol1","FW");
     if (make_plots){
       char pname[80];
       sprintf(pname,"par%i.png",ipar);
       cyy->Print(pname);
     }
     TF1 *fitbg = hpar[ipar]->GetFunction("fpol1");
     p_ht[ipar][0] = fitbg->GetParameter(0);
     p_ht[ipar][1] = fitbg->GetParameter(1);
     cout << p_ht[ipar][0]<< ", "<< p_ht[ipar][1]<<endl;
     //stop so one can look at the fitbg                                                             
     //int q;                                                                                      
     //cout <<"enter 1 to continue"<<endl;                                                         
     //cin>>q;                                                                                     
   }

   cyy->Update();

   //draw the HT vs Mfitbg plot and attach the 2D function                                           
   TCanvas *c1 = new TCanvas("c1","HT vs Mfitbg plot",10,10,700,500);
   c1->SetFillColor(0);
   c1->SetFrameFillColor(0);
   c1->cd();
   hs->Draw("col");

   TF2 *f2D = new TF2("f2D",func2D_sig,100,1000,200,2000,18);

   // set the parameters                                                                         
   for (int ipar=0;ipar<9;ipar++){
     f2D->FixParameter(ipar,p_mfit[ipar]);
   }
   f2D->FixParameter(9,p_ht[1][0]);
   f2D->FixParameter(10,p_ht[1][1]);
   f2D->FixParameter(11,p_ht[2][0]);
   f2D->FixParameter(12,p_ht[2][1]);
   f2D->FixParameter(13,p_ht[3][0]);
   f2D->FixParameter(14,p_ht[3][1]);
   f2D->FixParameter(15,p_ht[4][0]);
   f2D->FixParameter(16,p_ht[4][1]);
   f2D->SetParameter(17,38.3);

   cout << "parameters set" << endl;
   hs->Fit("f2D","l0");
   cout << "performed signal fit" << endl;


   // now the same for background...
   //fit the Mfit projection (x)                                                                 
   cx->cd();
   TH1D *hxbg = hb->ProjectionX("hxbg",0,ny);
   hxbg->SetFillColor(38);
   //   char title[80];
   sprintf(title,"Mfit Projection");
   hxbg->SetName("hxbg");
   hxbg->SetTitle(title);

   //gaussian + x^2*exp                                                                           
   //TF1 *fxbg = new TF1("fxbg","[0]*(p-[1])*(p-[1])*exp(-[2]*(p-[1]))+[3]*exp(-0.5*(p-[4])*(p-[4])/[5]/[5])",100,1000);
   TF1 *fxbg = new TF1("fxbg","[0]*(x-[1])*(x-[1])*exp(-[2]*(x-[1]))+[3]*exp(-0.5*(x-[4])*(x-[4])/[5]/[5])",100,1000);
   fxbg->SetParameters(0.1,100,0.02,100,175,15);
   hxbg->Fit("fxbg","l");
   hxbg->GetFunction("fxbg")->SetLineColor(kRed);
   hxbg->GetFunction("fxbg")->SetLineWidth(6);
   cx->Update();

   //save the parameters in p_mfit                                                                
   Double_t p_mfitbg[6];
   TF1 *fitbg = hxbg->GetFunction("fxbg");
   for (int ipar=0;ipar<6;ipar++){
     p_mfitbg[ipar] = fitbg->GetParameter(ipar);
   }

   //fit the HT projection (y) - for kicks only, this is not used again                           
   /*
   cy->cd();
   TH1D *hybg = hb->ProjectionY("",0,nx);
   hybg->SetFillColor(38);
   //   char title[80];

   sprintf(title,"HT Projection");
   hybg->SetName("hybg");
   hybg->SetTitle(title);
   TF1 *fybg = new TF1("fybg","[0]*[0]*exp(-0.5*(x-[1])*(x-[1])/[2]/[2])+[3]*[3]*exp(-0.5*(x-[4])*(x-[4])/[5]/[5])+[6]*[6]*exp(-0.5*(x-[7])*(x-[7])/[8]/[8])+[9]*[9]*exp(-0.5*(x-[10])*(x-[10])/[11]/[11])",100,1000);
   fybg->SetParameters(20,466,10,4,860,245,14,628,100,60,1000);
   fybg->SetParameter(11,700);
   hybg->Fit("fybg","l");
   hybg->GetFunction("fybg")->SetLineColor(kRed);
   hybg->GetFunction("fybg")->SetLineWidth(6);
   cy->Update();
   */

   //fit HT slices with double Gaussian                                                           
   //book histograms of parameter values vs Mfit                                                 
   TH1D *hparbg[5];
   for (int ipar=0;ipar<5;ipar++){
     sprintf(buf2,"hparbg_%d",ipar);
     //hparbg[ipar] = hb->ProjectionX("",1,1);
     hparbg[ipar] = hb->ProjectionX(buf2,1,1);
     hparbg[ipar]->Reset();
     char hname[20];
     sprintf(hname,"Par%i",ipar);
     hparbg[ipar]->SetName(hname);
   }

   //   TCanvas *cyy = (TCanvas*)gROOT->GetListOfCanvases()->FindObject("cyy");
   if(cyy) delete cyy->GetPrimitive("Projection");
   else   cyy = new TCanvas("cyy","Projection Canvas",700,10,490,350);
   cyy->SetGrid();
   cyy->cd();
   //double gaussian                                                                              
   TF1 *f1bg = new TF1("f1bg","[0]*((1-[4])*exp(-0.5*(x-[1])*(x-[1])/[2]/[2])/[2]+[4]*exp(-0.5*(x-1.3*[1])*(x-1.3*[1])/[3]/[3])/[3])",0,1000);
   cout << "number of x bins: " << nx <<endl;
   //loop over the HT slices                                                                     

   for (int binx=1;binx<25;binx++){
     TH1D *hp = hb->ProjectionY("hp",binx,binx);
     hp->Scale(1/hp->Integral());
     hp->SetFillColor(38);
     char title[80];
     sprintf(title,"Slice of binx=%d",binx);
     hp->SetName("Slice");
     hp->SetTitle(title);

     //if the error in a bin is tiny make it bigger so these bins don't drive the fit           
     for (int ibin=1;ibin<ny+1;ibin++){
       float ebin=hp->GetBinError(ibin);
       if (ebin>0 && ebin<0.01) hp->SetBinError(ibin,0.01);
       float cbin=hp->GetBinContent(ibin);
       if (cbin>0) cout <<"bin "<<ibin<<" = "<<cbin<<" +/- "<<ebin<<endl;
     }

     //first fit a gaussian to get starting values for the parameters                           
     hp->Fit("gaus","l0");
     TF1 *fit = hp->GetFunction("gaus");
     Double_t p0 = fit->GetParameter(0);
     Double_t p1 = fit->GetParameter(1);
     Double_t p2 = fit->GetParameter(2);

     cyy->Update();

     f1bg->SetParameters(p0,0.9*p1,0.5*p2,1.5*p2,0.25);
     f1bg->SetParLimits(4,0,0.5);
     f1bg->SetParLimits(3,70,500);
     f1bg->SetParLimits(2,20,1000);
     cout << "parameters: "<<p0<<","<<p1<<","<<p2<<endl;
     //     int q;                                                                               
     //     cout <<"continue?"<<endl;                                                            
     //     cin>>q;                                                                               
     //then fit the real thing                                                                     
     int Q=1;
     int nint=0;
     do {
       hp->Fit("f1bg","lm");
       Double_t chi2 = f1bg->GetChisquare();
       if (chi2>300 && nint<5){
         float x1=gRandom->Gaus(1,0.2);
         float x2=gRandom->Gaus(2,1);
         f1bg->SetParameters(p0,p1,x1*p2,x2*p2,0.25);
         nint++;
       }else{
         Q=0;
       }
     }while (Q==1);
     hp->GetFunction("f1bg")->SetLineColor(kRed);
     hp->GetFunction("f1bg")->SetLineWidth(6);
     cyy->Update();
     //make a png file of the canvas                                                             
     if (make_plots){
       char pname[80];  
       sprintf(pname,"hist%i.png",binx);
       cyy->Print(pname);
     }
     // write the fit results into the hpar histograms                                           
     Double_t par,err;
     for (int ipar=0;ipar<5;ipar++){
       par=f1bg->GetParameter(ipar);
       err=f1bg->GetParError(ipar);
       if(err>999)err=0;
       hparbg[ipar]->SetBinContent(binx,par);
       hparbg[ipar]->SetBinError(binx,err);
     }
     hp->Delete();
   }
   //end of loop over HT slice fits                                                                

   //now fit the parameters vs Mfit with a line                                                    
   Double_t p_htbg[5][2];
   for (int ipar=0;ipar<5;ipar++){
     hparbg[ipar]->Fit("pol1","");
     if (make_plots){
       char pname[80];  
       sprintf(pname,"par%i.png",ipar);
       cyy->Print(pname);
     }
     TF1 *fit = hparbg[ipar]->GetFunction("pol1");
     p_htbg[ipar][0] = fit->GetParameter(0);
     p_htbg[ipar][1] = fit->GetParameter(1);
     cout << p_htbg[ipar][0]<< ", "<< p_htbg[ipar][1]<<endl;
   }
  if (make_plots){  
    cyy->Update();
  }
   //draw the HT vs Mfit plot and attach the 2D function                                           
   c1->cd();
   hb->Draw("col");
   TF2 *f2Dbg = new TF2("f2Dbg",func2D,100,1000,200,2000,18);

   // set the parameters                                                                              
   for (int ipar=0;ipar<6;ipar++){
     f2Dbg->FixParameter(ipar,p_mfitbg[ipar]);
   }
   f2Dbg->FixParameter(9,p_htbg[1][0]);
   f2Dbg->FixParameter(10,p_htbg[1][1]);
   f2Dbg->FixParameter(11,p_htbg[2][0]);
   f2Dbg->FixParameter(12,p_htbg[2][1]);
   f2Dbg->FixParameter(13,p_htbg[3][0]);
   f2Dbg->FixParameter(14,p_htbg[3][1]);
   f2Dbg->FixParameter(15,p_htbg[4][0]);
   f2Dbg->FixParameter(16,p_htbg[4][1]);
   f2Dbg->SetParameter(17,38.3);

   cout << "parameters set" << endl;

   hb->Fit("f2Dbg","l0");

   cout << "performed background fit" << endl;



   /*

  TH1D * hs_mass = hs->ProjectionX("hs_mass", 0, -1, "e");
  double _xmin = hs_mass->GetXaxis()->GetXmin();
  double _xmax = hs_mass->GetXaxis()->GetXmax();
  double _nbin = hs_mass->GetNbinsX();
  double _bin_x_width = (_xmax-_xmin)/_nbin;
  int ngroup = (int)( 20.5/((_xmax-_xmin)/_nbin) );
  hs_mass->Rebin(ngroup);

  TF1 * fs_mass = new TF1("fs_mass", gaussian_mass, _xmin, _xmax, 6);
  //TF1 * fs_mass = new TF1("fs_mass", gaussian_mass, 150.0, 700.0, 6);
  //TF1 * fs_mass = new TF1("fs_mass", gaussian_mass, 115, _xmax, 6);
  //fs_mass->SetParameters(1,400,100,1,200,200);
  fs_mass->SetParameters(1,mass,50,1,mass/2.0,100);
  fs_mass->SetParLimits(0,0,10000000);
  //fs_mass->SetParLimits(1,100,1000);
  fs_mass->SetParLimits(1,mass/2.0,1000);
  //fs_mass->SetParLimits(2,10,1000);
  fs_mass->SetParLimits(2,10,mass);
  fs_mass->SetParLimits(3,0,10000000);
  //fs_mass->SetParLimits(4,100,1000);
  fs_mass->SetParLimits(4,100,mass);
  fs_mass->SetParLimits(5,10,1000);
  hs_mass->Fit("fs_mass","R");
  if (make_plots){
    TCanvas c2;
    hs_mass->Draw();
    _fit_fname = channel+"_sig_mass.pdf";
    c2.SaveAs( _fit_fname.c_str() );
  }
  TH1D * hs_ht = hs->ProjectionY("hs_ht", 0, -1, "e");
  double _ymin = hs_ht->GetXaxis()->GetXmin();
  double _ymax = hs_ht->GetXaxis()->GetXmax();
  double _nbiny = hs_ht->GetNbinsX();
  double _bin_y_width = (_ymax-_ymin)/_nbiny;
  ngroup = (int)( 20.5/((_ymax-_ymin)/_nbiny) );
  hs_ht->Rebin(ngroup);
  TF1 * fs_ht = new TF1("fs_ht", gaussian_ht, _ymin, _ymax, 6);
  //TF1 * fs_ht = new TF1("fs_ht", gaussian_ht, 300, 1400, 6);
  //TF1 * fs_ht = new TF1("fs_ht", gaussian_ht, 300, 2000, 6);
  if ( channel.find("ejets")!=std::string::npos && abs(mass-625)<1 ){
    delete fs_ht;
    fs_ht = new TF1("fs_ht", gaussian_ht, 300, 2000, 6);
  }
  if ( channel.find("ejets")!=std::string::npos && abs(mass-600)<1 ){
    delete fs_ht;
    fs_ht = new TF1("fs_ht", gaussian_ht, 400, 2000, 6);
  }
  if ( channel.find("ejets")!=std::string::npos && abs(mass-550)<1 ){
    delete fs_ht;
    fs_ht = new TF1("fs_ht", gaussian_ht, 300, 2000, 6);
  }
  if ( channel.find("ejets")!=std::string::npos && abs(mass-475)<1 ){
    delete fs_ht;
    fs_ht = new TF1("fs_ht", gaussian_ht, 300, 2000, 6);
  }
  if ( channel.find("ejets")!=std::string::npos && abs(mass-450)<1 ){
    delete fs_ht;
    fs_ht = new TF1("fs_ht", gaussian_ht, 400, 2000, 6);
  }
  if ( channel.find("mujets")!=std::string::npos && abs(mass-450)<1 ){
    delete fs_ht;
    fs_ht = new TF1("fs_ht", gaussian_ht, 300, 2000, 6);
  }
  if ( channel.find("ejets")!=std::string::npos && abs(mass-425)<1 ){
    delete fs_ht;
    fs_ht = new TF1("fs_ht", gaussian_ht, 300, 2000, 6);
  }
  if ( channel.find("mujets")!=std::string::npos && abs(mass-425)<1 ){
    delete fs_ht;
    fs_ht = new TF1("fs_ht", gaussian_ht, 300, 2000, 6);
  }
  //
  //fs_ht->SetParameters(1,1000,200,1,500,100);
  ////fs_ht->SetParameters(8,820,140,6,1070,190);
  //fs_ht->SetParLimits(0,0,1000);
  //fs_ht->SetParLimits(1,200,2000);
  //fs_ht->SetParLimits(2,10,1000);
  //fs_ht->SetParLimits(3,0,1000);
  //fs_ht->SetParLimits(4,200,2000);
  //fs_ht->SetParLimits(5,10,1000);
  //
  fs_ht->SetParameters(8,0.5,800,500,1000,500);
  fs_ht->SetParLimits(0,0,1000);
  fs_ht->SetParLimits(1,0.1,0.9);
  fs_ht->SetParLimits(2,200,2000);
  fs_ht->SetParLimits(3,140,1000);
  fs_ht->SetParLimits(4,200,2000);
  fs_ht->SetParLimits(5,140,1000);
  //
  hs_ht->Fit("fs_ht","R");
  if (make_plots){
    TCanvas c3;
    hs_ht->Draw();
    _fit_fname = channel+"_sig_ht.pdf";
    c3.SaveAs( _fit_fname.c_str() );
  }
  TH1D * hb_mass = hb->ProjectionX("hb_mass", 0, -1, "e");
  _xmin = hb_mass->GetXaxis()->GetXmin();
  _xmax = hb_mass->GetXaxis()->GetXmax();
  _nbin = hb_mass->GetNbinsX();
  ngroup = (int)( 20.5/((_xmax-_xmin)/_nbin) );
  hb_mass->Rebin(ngroup);
  TF1 *fb_mass = new TF1("fb_mass", bgMassPdf,_xmin,_xmax,7);
  //TF1 *fb_mass = new TF1("fb_mass", bgMassPdf,_xmin,600.0,7);
  fb_mass->SetParameters(1,0.5,0.0,1,1,175,10); // LogNormal
  fb_mass->SetParLimits(1,0.3,1.0); // LogNormal
  fb_mass->SetParLimits(2,0.0,0.01); // LogNormal
  fb_mass->SetParLimits(3,0.0,5.0); // LogNormal
  fb_mass->SetParLimits(5,100.0,200.0); // Gaussian
  fb_mass->SetParLimits(6,1.0,50.0); // Gaussian
  hb_mass->Fit("fb_mass","R");
  if (make_plots){
    TCanvas c4;
    hb_mass->Draw();
    _fit_fname = channel+"_bkg_mass.pdf";
    c4.SaveAs( _fit_fname.c_str() );
  }
  TH1D * hb_ht = hb->ProjectionY("hb_ht", 0, -1, "e");
  _ymin = hb_ht->GetXaxis()->GetXmin();
  _ymax = hb_ht->GetXaxis()->GetXmax();
  _nbiny = hb_ht->GetNbinsX();
  ngroup = (int)( 20.5/((_ymax-_ymin)/_nbiny) );
  hb_ht->Rebin(ngroup);
  TF1 *fb_ht = new TF1("fb_ht", bgHtPdf,_ymin,_ymax,6);
  //TF1 *fb_ht = new TF1("fb_ht", bgHtPdf,250.0,1000.0,6);
  fb_ht->SetParameters(100,545,60,100,700,50);
  fb_ht->SetParLimits(1,200,2000);
  fb_ht->SetParLimits(2,30,1000);
  fb_ht->SetParLimits(4,200,2000);
  fb_ht->SetParLimits(5,30,1000);
  hb_ht->Fit("fb_ht","R");
  if (make_plots){
    TCanvas c5;
    hb_ht->Draw();
    _fit_fname = channel+"_bkg_ht.pdf";
    c5.SaveAs( _fit_fname.c_str() );
  }
  */

  // try to make a corresponding 2D function and normalize it
  //TF2 * _f2s = new TF2("f2s", "[0]*gaussian_mass(x)*gaussian_ht(y)", 100,1000,200,2000);
  //double _int2s = _f2s->Integral(100,1000,200,2000);
  //_f2s->SetParameter(0,1.0);
  //std::cout << " *** Integral2d = " << _int2s << std::endl;

  //--------- end of fits


  float xs,xb,eb;
  //xsb holds the s/b values for each bin
  //xx are the histogram contents 
  //(0=signal, 1=total background, 2=data, 3...nb-1=individual backgrounds)
  float xsb[30000],xx[30000][30],xe[30000][30];
  int ibin; 
  double _sum = 0.0;
  TH1F hLogSb("hLogSb", "log(s/b)", 100, -40, 60);
  double _logsb;
  for (int i=0;i<nx;i++){
    for (int j=0;j<ny;j++){
      
      ibin=hs->GetBin(i,j);

      // GENA: Will fill each bin with its original index
      vMap[ibin].push_back(ibin);

      xs=hs->GetBinContent(ibin);
      xb=hb->GetBinContent(ibin);
      // This needs to change for 2D...
      // location in x,y... is...
      /*      _xmin = hs->GetXaxis()->GetXmin();
      _xmax = hs->GetXaxis()->GetXmax(); 
      _ymin = hs->GetYaxis()->GetYmin();
      _ymax = hs->GetYaxis()->GetYmax(); 
      */
      //compute signal/background
      double _x = hs->GetXaxis()->GetBinCenter(i); //_xmin + min(_xmax, max(0.0, ((double)(i)-0.5)))*_bin_x_width;
      double _y = hs->GetYaxis()->GetBinCenter(j);
      //double _y = _ymin + min(_ymax, max(0.0, ((double)(j)-0.5)))*_bin_y_width;

      //only the SB value changes now, everything else can remain the same...
      xsb[ibin] = log(1+  f2D->Eval(_x,_y)   / f2Dbg->Eval(_x,_y)
		      );/*log(
		      1.0+ fs_mass->Eval(_x)*fs_ht->Eval(_y) 
		      / fb_mass->Eval(_x)/fb_ht->Eval(_y)
		      );*/
      _logsb = log(f2D->Eval(_x,_y)   / f2Dbg->Eval(_x,_y)  );
      //edit this bunch
      //std::cout << "DEBUG*****!!!!!  145/150: " << fs_mass->Eval(145.0)
      //	<< "  /  " << fs_mass->Eval(155.0)
      //	<< std::endl;

      // s/b: 10 00 11 01
      /*
      if (xb>0){
	if (xs>0){ // [1;1000]
	  xsb[ibin]=xs/xb+1.0;
	}
	else{      // [0;1]
	  //xsb[ibin]=1.0 - xb/(xb+1.0);
	  xsb[ibin]=((float)ibin)/((float)ibin+1.0);
	}
      }else{
        if (xs>0){ // (2000; up)
          //xsb[ibin]=2000.0+xs;
          xsb[ibin]=2000.0+(float)ibin;
        }else{     // [1000;1001]
          xsb[ibin]=ibin/10000.0+1000.0;
        }
      }
      */

      // original Ulrich:
      //if (xb>0) {
      //  xsb[ibin]=xs/xb;
      //}else{
      //  if (xs>0){
      //    xsb[ibin]=999;
      //  }else{
      //    xsb[ibin]=0;
      //  }
      //}
      
      xx[ibin][0]=xs;
      xe[ibin][0]=hs->GetBinError(ibin);
      xx[ibin][1]=xb;
      xe[ibin][1]=hb->GetBinError(ibin);
      for (int k=0;k<nb;k++){
	xx[ibin][k+2]=hX[k]->GetBinContent(ibin);
	xe[ibin][k+2]=hX[k]->GetBinError(ibin);
      }
      //if (xb>0) h1sb->SetBinContent(ibin,xs/xb);
      //hLogSb.Fill(xsb[ibin]);
      
      hLogSb.Fill(_logsb);
      h1sb->SetBinContent(ibin,xsb[ibin]);
      h1s->SetBinContent(ibin,xx[ibin][0]);
      h1s->SetBinError(ibin,xe[ibin][0]);
      h1b->SetBinContent(ibin,xx[ibin][1]);
      h1b->SetBinError(ibin,xe[ibin][1]);
      
      _sum += xx[ibin][0];
    }
  }

  std::cout << "SUM: " << _sum << std::endl;

  // clean up fits
  delete f2D;//s_mass;
  delete f2Dbg;//fb_mass;
  //  delete fs_ht;
  //delete fb_ht;
  // FIXME: 
  //delete _f2s;


  if (make_plots){  
    // plot S/B hist
    TCanvas c1;
    c1.SetLogy();
    h1sb->SetMinimum(0.00001);
    h1sb->Draw();
    std::string _fit_fname;
    _fit_fname = channel+"_sb.pdf";
    c1.SaveAs( _fit_fname.c_str() );
    TCanvas c6;
    //c6.SetLogy();
    //hLogSb->SetMinimum(0.00001);
    hLogSb.Draw();
    _fit_fname = channel+"_logsb.pdf";
    c6.SaveAs( _fit_fname.c_str() );
  }

  //sort all histogram bins in decreasing s/b
  int nswap=1;
  float xtmp;

  // GENA: for bin map
  int ibin_tmp;

  while (nswap>0) {
    nswap=0;
    for (int i=0;i<nbin-1;i++) {
      if (xsb[i]<xsb[i+1]){
	xtmp=xsb[i];
	xsb[i]=xsb[i+1];
	xsb[i+1]=xtmp;

	// GENA: for bin map
	ibin_tmp = vMap[i][0];
	vMap[i][0] = vMap[i+1][0];
	vMap[i+1][0] = ibin_tmp;

	for (int j=0;j<nb+2;j++){
	  xtmp=xx[i][j];
	  xx[i][j]=xx[i+1][j];
	  xx[i+1][j]=xtmp;

	  xtmp=xe[i][j];
	  xe[i][j]=xe[i+1][j];
	  xe[i+1][j]=xtmp;
	}
	nswap=nswap+1;
      }
    }
  }

  //these histograms have the bins ordered in decrerasing s/b for diagnostics
  TH1F *h1sb1 = new TH1F("h1sb1","h1sb1",nbin,0,nbin);
  TH1F *h1fe1 = new TH1F("h1fe1","h1fe1",nbin,0,nbin);
  TH1F *h1s1 = new TH1F("h1s1","h1s1",nbin,0,nbin);
  TH1F *h1b1 = new TH1F("h1b1","h1b1",nbin,0,nbin);
  for (int i=0;i<nbin;i++){
    h1sb1->SetBinContent(i+1,xsb[i]);
    if (xx[i][1]>0) h1fe1->SetBinContent(i+1,xe[i][1]/xx[i][1]);
    h1s1->SetBinContent(i+1,xx[i][0]);
    h1s1->SetBinError(i+1,xe[i][0]);
    h1b1->SetBinContent(i+1,xx[i][1]);
    h1b1->SetBinError(i+1,xe[i][1]);
  }

  //combine bins starting with the highest s/b until the fractional error in
  //the total backround in every bin is smaller than femax
  int ncomb=1;
  //float xtmp;
  float fe=0;
  double _max_bin_err = 1.0;
  while (ncomb>0) {
    ncomb=0;
    for (int i=0;i<nbin-1;i++){

      bool err_too_large = false;
      for (int m=0; m!=2; ++m){
	if (m==2) continue;
	if (xx[i][m]>0){
	  fe=xe[i][m]/xx[i][m]; //fractional error in signal (0) and background (1)
	}else{
	  fe=1;
	}

	if ( m<vFemax.size() ) _max_bin_err = vFemax[m];
	else _max_bin_err = vFemax[1];

	if (fe>_max_bin_err) err_too_large = true;
      }

      if (err_too_large){

	// GENA: write down bin
	for (std::vector<int>::const_iterator vi=vMap[i+1].begin();
	     vi != vMap[i+1].end(); ++vi){
	  vMap[i].push_back(*vi);
	}
	//move all successive bins up
	vMap.erase(vMap.begin()+i+1);

	for (int k=0;k<nb+2;k++){ //add the next bin
	  xx[i][k]=xx[i][k]+xx[i+1][k];
	  xe[i][k]=sqrt(xe[i][k]*xe[i][k]+xe[i+1][k]*xe[i+1][k]);
	  for (int j=i+1;j<nbin-1;j++){ //move all successive bins up
	    xx[j][k]=xx[j+1][k];
	    xe[j][k]=xe[j+1][k];
	  }
	}

	ncomb++;
	nbin=nbin-1; //decrease the total number of bins
      }

    }
    
  }
  // check the last bin and add it to the previous if uncertainty is high
  int i = nbin-1;

  bool err_too_large = false;
  for (int m=0; m!=2; ++m){
    if (m==2) continue;
    if (xx[i][m]>0){
      fe=xe[i][m]/xx[i][m]; //fractional error in background
    }else{
      fe=1;
    }

    if ( m<vFemax.size() ) _max_bin_err = vFemax[m];
    else _max_bin_err = vFemax[1];

    if (fe>_max_bin_err) err_too_large = true;
  }
  
  if (err_too_large){

    std::cout << "Fixing the last bin" << std::endl;

    i = nbin-2;
    
    // GENA: write down bin
    for (std::vector<int>::const_iterator vi=vMap[i+1].begin();
	 vi != vMap[i+1].end(); ++vi){
      vMap[i].push_back(*vi);
    }
    //move all successive bins up
    vMap.erase(vMap.begin()+i+1);
    
    for (int k=0;k<nb+2;k++){ //add the next bin
      xx[i][k]=xx[i][k]+xx[i+1][k];
      xe[i][k]=sqrt(xe[i][k]*xe[i][k]+xe[i+1][k]*xe[i+1][k]);
      for (int j=i+1;j<nbin-1;j++){ //move all successive bins up
	xx[j][k]=xx[j+1][k];
	xe[j][k]=xe[j+1][k];
      }
    } // end of last bin treatment
    
    nbin=nbin-1; //decrease the total number of bins
  }
  
  // check the last bin again
  i = nbin-1;
  if (xx[i][1]>0){
    fe=xe[i][1]/xx[i][1]; //fractional error in background
  }else{
    fe=1;
  }
  std::cout << "Last bin error: " << fe << std::endl;




  //GENA: open the map file
  std::ofstream mapFile;
  std::string mapFileName = outFileName+".binmap";
  mapFile.open(mapFileName.c_str());
  int bin_count = 0;
  for (std::vector<std::vector<int> >::const_iterator i=vMap.begin();
       i != vMap.end(); ++i){

    mapFile << " " << i-vMap.begin()+1 << ":";

    for(std::vector<int>::const_iterator j=i->begin();
	j != i->end(); ++j){
      mapFile << " " << *j;
      ++bin_count;
    }

    mapFile << std::endl;

  }
  //GENA: close the map file
  mapFile.close();

  //these are the output histograms
  std::string sname_out = prefix + std::string(sname);
  TFile *f2 = TFile::Open(oname,"recreate");
  TH1F *h1feb2 = new TH1F("h1fe2","h1fe2",nbin,0,nbin);
  //TH1F *h1s2 = new TH1F(sname,sname,nbin,0,nbin);
  TH1F *h1s2 = new TH1F(sname_out.c_str(),sname_out.c_str(),nbin,0,nbin);
  TH1F *h1b2 = new TH1F("h1b2","h1b2",nbin,0,nbin);
  TH1F *h1X2[30];
  for (int i=0;i<nb;i++){
    std::string _bname = prefix + std::string(bname[i]);
    //h1X2[i] = new TH1F(bname[i],bname[i],nbin,0,nbin);
    h1X2[i] = new TH1F(_bname.c_str(),_bname.c_str(),nbin,0,nbin);
  }
  for (int i=0;i<nbin;i++){
    h1feb2->SetBinContent(i+1,xe[i][1]/xx[i][1]);
    h1s2->SetBinContent(i+1,xx[i][0]);
    h1s2->SetBinError(i+1,xe[i][0]);
    h1b2->SetBinContent(i+1,xx[i][1]);
    h1b2->SetBinError(i+1,xe[i][1]);
    for (int j=0;j<nb;j++){
      h1X2[j]->SetBinContent(i+1,xx[i][j+2]);
      h1X2[j]->SetBinError(i+1,xe[i][j+2]);
    }
  }
  
  std::cout << "Merged 1D hist name: " << h1s2->GetName() << std::endl;
  std::cout << "Integral with overflow: " << h1s2->Integral(0,nbin+1) << std::endl;
  std::cout << "Integral no overflow: " << h1s2->Integral(1,nbin) << std::endl << std::endl;
  h1s2->SetEntries(sEntries);
  h1s2->Write();
  for (int j=0;j<nb;j++){
    std::cout << "Merged 1D hist name: " << h1X2[j]->GetName() << std::endl;
    std::cout << "Integral with overflow: " << h1X2[j]->Integral(0,nbin+1) << std::endl;
    std::cout << "Integral no overflow: " << h1X2[j]->Integral(1,nbin) << std::endl << std::endl;
    h1X2[j]->SetEntries(bEntries[j]);
    h1X2[j]->Write();
  }
  
  h1s2->Print("base");

  f2->Close();
  f->Close();

  //std::cout << "map size: " << vMap.size() << " combined bins" << std::endl;
  //std::cout << "total bins merged: " << bin_count << std::endl;
}



