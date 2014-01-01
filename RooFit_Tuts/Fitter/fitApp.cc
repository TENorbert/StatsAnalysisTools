// $Id: fitApp.cc,v 1.5 2012/05/30 13:25:12 sigamani Exp $
#include "TColor.h"
#include "TStyle.h"
#include "TMath.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TTree.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "RooRealVar.h"
#include "RooFormulaVar.h"
#include "RooDataSet.h"
#include "RooGenericPdf.h"
#include "RooCBShape.h"
#include "RooAddPdf.h"
#include "RooProdPdf.h"
#include "RooDataHist.h"
#include "RooChi2Var.h"
#include "RooMinuit.h"
#include "RooPlot.h"
#include "RooHistPdf.h"
#include "RooHist.h"
#include "RooNLLVar.h"
#include "RooBinning.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TPolyLine.h"
#include "RooGlobalFunc.h"
#include "TROOT.h"
#include "RooFitResult.h"


using namespace RooFit;
using namespace std;




std::vector<double> eff_num ;
std::vector<double> ul ;

std::vector<TFile*> fvec;
std::vector<TFile*> fvec_pj;
std::vector<TFile*> fvec_data;
std::vector<TFile*> fgmsb;
std::vector<double> efficiency;




TFile* fvec_ttbar;
TFile* fvec_Wenu;

void rootlogon();

double normGauss(double mu,double sig);
TGraph* convoluteGraph(const TGraph* ingr, const double& smear);

void computeUpperLimit(TGraph* hist, double xCL);
TGraph* fit_2D( int GMSBindex, float lumin, std::string ctau, double qcdscale, double gjetscale);

void plot_limit(std::string ctau);



void dumpTH1(const TH1* hin) {
  hin->Print("v");
  for(int i=0; i<=hin->GetNbinsX()+1; ++i)
   cout << "i: " << i << "  y: " << hin->GetBinContent(i) << endl;
  return; 
}

void protectZeroBin(TH1* hin) {
  hin->Print("v");

   for(int i = 1; i<=hin->GetNbinsX()+1; i++){
     for(int j =1; j<= hin->GetNbinsY()+1;j++){
       if(hin->GetBinContent(i, j) <=0)  { 
             hin->SetBinContent(i, j, 1.e-5);
             hin->SetBinError(i, j, 1.e-5);
       }
     }
   }
}




TH2F* hdensity(const TH2F* hin){

  TH2F* hout = (TH2F*) hin->Clone("hout");
  
  for(int i = 1; i<=hout->GetNbinsX(); i++){
    for(int j =1; j<= hout->GetNbinsY();j++){
      double area = hout->ProjectionX()->GetBinWidth(i)*hout->ProjectionY()->GetBinWidth(j);
      
      hout->SetBinContent(i,j,(double)hout->GetBinContent(i,j)/area);
      hout->SetBinError(i,j,(double)hout->GetBinError(i,j)/area);

      // std::cout << "i: " << i << " j: " << j <<
	//	" val: " << hout->GetBinContent(i,j) << " area: " << area <<std::endl; 

  }
 }  
  return hout;  
}



TH1* hRedensity(const TH1* hin){

  TH1* hout = (TH1*) hin->Clone("hout");
  
  double density = (double)(hin->GetXaxis()->GetXmax()-hin->GetXaxis()->GetXmin())/hin->GetNbinsX()  ;
  std::cout << density << " : density" << std::endl;
  for(int i = 1; i<=hout->GetNbinsX(); i++){
   
      double area = hout->GetBinWidth(i);
      
      hout->SetBinContent(i,(double)hout->GetBinContent(i)*density/area);
      hout->SetBinError(i,(double)hout->GetBinError(i)*density/area);

      // std::cout << "i: " << i << " j: " << j <<
	//	" val: " << hout->GetBinContent(i,j) << " area: " << area <<std::endl; 

  
 }  
  return hout;  
}



int main(int argc, char* argv[]) {

  rootlogon();

  std::string LAMBDA = argv[1];
  std::string  CTAU = argv[2];

  TString dirpath("/u1/sigamani/default/");



  fgmsb.push_back( (TFile*)  TFile::Open(dirpath+"PhotonJet_GMSB_Lambda-"+LAMBDA+"_CTau-"+CTAU+"_7TeV_pythia6_cff-Summer11-PU_S4_START42_V11-v1-AODSIM-UNCLEAN25_pfakt5_2ndJet10_default.root"));
  fvec.push_back( (TFile*) TFile::Open(dirpath+"PhotonJet_DATA_2011_ALL_pfakt5_CONTROL_QCD_2ndJet10_default.root"));
  fvec_pj.push_back( (TFile*) TFile::Open(dirpath+"PhotonJet_DATA_2011_ALL_pfakt5_CONTROL_G_2ndJet10_default.root"));
  fvec_data.push_back( (TFile*) TFile::Open(dirpath+"PhotonJet_DATA_2011_ALL_pfakt5_2ndJet10_default.root"));
  fvec_ttbar = (TFile*) TFile::Open(dirpath+"PhotonJet_TT_2011_pfakt5_2ndJet10_default.root");
  fvec_Wenu = (TFile*) TFile::Open(dirpath+"PhotonJet_WToENU_2011_pfakt5_2ndJet10_default.root");

  double  errsig = 0.;


  if ( (CTAU=="1") && (LAMBDA=="120") ){ efficiency.push_back(0.2155+(errsig*0.2155)) ; }
  if ( (CTAU=="1") && (LAMBDA=="140") ){ efficiency.push_back(0.2746+(errsig*0.2746)) ; }
  if ( (CTAU=="1") && (LAMBDA=="160") ){ efficiency.push_back(0.3194+(errsig*0.3194)) ; }
  if ( (CTAU=="1") && (LAMBDA=="180") ){ efficiency.push_back(0.9999+(errsig*0.9999)) ; }

  if ( (CTAU=="250") && (LAMBDA=="120") ){ efficiency.push_back(0.19816+(errsig*0.19816)) ; }
  if ( (CTAU=="250") && (LAMBDA=="140") ){ efficiency.push_back(0.26474+(errsig*0.26474)) ; }
  if ( (CTAU=="250") && (LAMBDA=="160") ){ efficiency.push_back(0.30923+(errsig*0.30923)) ; }
  if ( (CTAU=="250") && (LAMBDA=="180") ){ efficiency.push_back(0.34342+(errsig*0.34342)) ; }

  if ( (CTAU=="1") && (LAMBDA=="100") ){ efficiency.push_back(0.999+(errsig*0.999)) ; }
  if ( (CTAU=="250") && (LAMBDA=="100") ){ efficiency.push_back(0.1638+(errsig*0.1638)); }
  if ( (CTAU=="500") && (LAMBDA=="100") ){ efficiency.push_back(0.1469+(errsig*0.1469)) ; }
  if ( (CTAU=="1000") && (LAMBDA=="100") ){ efficiency.push_back(0.1105+(errsig*0.1105)) ; }
  if ( (CTAU=="2000") && (LAMBDA=="100") ){ efficiency.push_back(0.0718+(errsig*0.0718)) ; }
  if ( (CTAU=="3000") && (LAMBDA=="100") ){ efficiency.push_back(0.9999+(errsig*0.9999)) ; }
  if ( (CTAU=="4000") && (LAMBDA=="100") ){ efficiency.push_back(0.0395+(errsig*0.0395)) ; }
  if ( (CTAU=="6000") && (LAMBDA=="100") ){ efficiency.push_back(0.0284+(errsig*0.0284)) ; }

 
  double gjetscale = 0.42;
  double qcdscale = 0.58;
 
  // contains the limit graph
  TGraph *t=0;
  t = fit_2D( 0, 4553.1, CTAU, qcdscale, gjetscale ); // take values from command line
  computeUpperLimit(t,0.95);
  plot_limit("250");
} 







//===================================
void computeUpperLimit(TGraph* hist, double xCL) {
//===================================

  const int nbin = hist->GetN();
  double* gy = hist->GetY(); 
  double* gx = hist->GetX(); 
  double* dx = new double[nbin];

  cout <<"Input graph with " << nbin << " bins" << endl;

  double tot = 0.0;
  for(int i=0; i<nbin; i++) {
   if(gx[i]<0) continue;
   if(i==0) {
     dx[i] = gx[i+1]-gx[i];
   } else if(i==nbin-1) {
     dx[i] = gx[i]-gx[i-1];
   } else {
    dx[i] = 0.5*(gx[i+1]-gx[i]) + 0.5*(gx[i]-gx[i-1]);
   }
   tot += gy[i]*dx[i];
   //cout <<  "bin: "<< i << " dx: "<< dx[i] << " y: " << gy[i] << endl;
  }
  cout << "Integral of input graph with positive prior: " << tot << endl;


  std::cout<<nbin<<std::endl;
  // truncated graph
  double newx[nbin];
  double newy[nbin];

  double cumul = 0.0;
  int iCL = 0;
  for(int i=0; i<nbin; i++) {
    newx[i] = gx[i];
    //if (gx[i]<0.)  continue; // postive prior 
    if(gx[i]>0.) {
       cumul  +=  gy[i]*dx[i]/tot;
       newy[i] = gy[i];
    } else {
       newy[i] = 0.0;
    }
     //cout <<  "bin: "<< i << " y: " << gy[i] << " frac cumul: " << cumul << endl;
    if(cumul>=xCL) { 
      cout << "Reached the " << xCL << " limit at bin " << i << " with x = " << gx[i] << endl;
      iCL = i; break;
    }// upper limit
  }

  if(iCL>0 && iCL<=nbin) { 
    double errPos = newx[iCL];
    //cout << "Upper limit: " << errPos  << " x 10^-5 at " << xCL << " C.L." << endl;
    cout << "Upper limit: " << errPos  << " at " << xCL << " C.L." << endl;

    for(int i=iCL+1; i<nbin; i++) { newx[i] = gx[i]; newy[i] = 0; }
  } else {
    cout << "Warning**** Upper limit was not found! Integral(0->+inf) = " << cumul << endl;
  }

  TGraph* newh = new TGraph(nbin,newx,newy);

  TCanvas* c1 = new TCanvas("c_scan","scan of LL",800,600);

  newh->SetFillColor(kBlue);

  //hist->SetTitle("B^{0}\\rightarrow D^{0} K^{*0}, All D^{0} modes combined");
  hist->SetTitle("");
  hist->GetXaxis()->SetTitle("Branching Fraction (x 10^{-5})");
  hist->GetYaxis()->SetTitle("Likelihood");
  hist->GetYaxis()->SetTitleOffset(1.2);

  hist->Draw("AC");
  newh->Draw("ACFsame");
  c1->SaveAs("scan-upperlimit.eps");
  
  ul.push_back( newx[iCL] );
}


//--------------
TGraph* fit_2D( int GMSBindex, float lumin,  std::string ctau, double  qcdscale, double gjetscale ){
//--------------
   char treeName[100];
   sprintf(treeName,"jetTree");

   //bool isMC = false;
   
     //---------------import data--------------//

   cout << "--- dataset for DATA --- " << endl;

   TString hfitname("met_varsize_x_timePhot_varsize");

   TH2F* h2_data= (TH2F*)fvec_data[0]->Get(hfitname)->Clone("h2_data");
   protectZeroBin(h2_data);
   //h2_data->Sumw2();

   cout << "projX: " << h2_data->ProjectionX() << "\tprojY: " << h2_data->ProjectionY() << endl; 


   RooRealVar epfMet("met","PF MET",h2_data->ProjectionX()->GetXaxis()->GetXmin(), h2_data->ProjectionX()->GetXaxis()->GetXmax(), "GeV");
   RooRealVar timePhotReco("time","ECAL time",h2_data->ProjectionY()->GetXaxis()->GetXmin(), h2_data->ProjectionY()->GetXaxis()->GetXmax(),"ns");

    cout << "RealVars created" << endl;

    //TH2F* h2_density = hdensity(h2_data);
    RooDataHist h_allDATA("h_allDATA","hist allDATA",RooArgSet(epfMet,timePhotReco), Import(*h2_data,false) );
    RooHistPdf h_allDATA_PDF("h_allDATA_PDF","met X time data",RooArgSet(epfMet,timePhotReco), h_allDATA ); 
   
    std::cout<<"RooDataHist create"<<std::endl;


    TH1* methist = h_allDATA.createHistogram("methist",epfMet);
    TH1* met_pdf_hist = hRedensity(methist);

    TH1* timehist = h_allDATA.createHistogram("methist",timePhotReco);
    TH1* time_pdf_hist = hRedensity(timehist);

    //------import shapes for the pdf------------//
   
    cout << "--- dataset for W -> e nu --- " << endl;
    // Wenu PDF
    TH2F* h2_Wenu_ = (TH2F*)fvec_Wenu->Get(hfitname);
    TH2F* h2_Wenu = (TH2F*)h2_Wenu_->Clone("h2_Wenu");
    protectZeroBin(h2_Wenu);

    h2_Wenu->Sumw2();
    h2_Wenu->Scale(lumin);
  
    RooDataHist hWenu("hWenu","hist Wenu",RooArgSet(epfMet,timePhotReco), Import(*h2_Wenu,false));

   TH1* methistwenu = hWenu.createHistogram("methistwenu",epfMet);
   dumpTH1(methistwenu);

   RooHistPdf Wenu_PDF("Wenu_PDF","met X time Wenu",RooArgSet(epfMet,timePhotReco), hWenu );

   std::cout << "WEnu integral: " << h2_Wenu->Integral() <<std::endl;
   
   // QCD control sample
   cout << "--- dataset for QCD --- " << endl;

   // original histo
   TH2F* h2_QCD_ = (TH2F*)fvec[0]->Get(hfitname);
   
   
   // new histo to get rid of empty bins causing trouble with LL fit
 

   TH2F* h2_QCD = (TH2F*)h2_QCD_->Clone("h2_QCD");
   protectZeroBin(h2_QCD);
   
   
   // scale control sample and MC
   h2_QCD->Sumw2();
 

   RooDataHist hQCD("hQCD","hist QCD",RooArgSet(epfMet,timePhotReco), Import(*h2_QCD,false) );
   TH1* methistqcd = hQCD.createHistogram("methistqcd",epfMet);
   dumpTH1(methistqcd);
 
   RooHistPdf QCD_PDF("QCD_PDF","met X time QCD",RooArgSet(epfMet,timePhotReco), hQCD ); 





   
   
   // signal PDF
   cout << "---- GMSB for signal --- " << endl;
   TH2F* h2_gmsb_=(TH2F*) fgmsb[GMSBindex]->Get(hfitname);
   TH2F* h2_gmsb = (TH2F*)h2_gmsb_->Clone("h2_gmsb"); // already normalized to 1 pb-1 in analysis code 
   protectZeroBin(h2_gmsb);

   h2_gmsb->Draw("hist");
   h2_gmsb->Sumw2();
   h2_gmsb->Scale(lumin);

   RooDataHist hGMSB("hGMSB","hist GMSB",RooArgSet(epfMet,timePhotReco), Import(*h2_gmsb,false) );
   RooHistPdf signal("signal","met X time GMSB",RooArgSet(epfMet,timePhotReco), hGMSB ); 
   std::cout<<"signal pdf"<<std::endl;


   // ttbar PDF
   TH2F* h2_ttbar_ = (TH2F*)fvec_ttbar->Get(hfitname);
   TH2F* h2_ttbar = (TH2F*)h2_ttbar_->Clone("h2_ttbar");
   protectZeroBin(h2_ttbar);

   h2_ttbar->Sumw2();
   h2_ttbar->Scale(lumin);
   
   RooDataHist hTTBAR("hTTBAR","hist TTBAR",RooArgSet(epfMet,timePhotReco), Import(*h2_ttbar,false) );
   TH1* methistttbar = hTTBAR.createHistogram("methistttbar",epfMet);
   dumpTH1(methistttbar);

   RooHistPdf TTbar_PDF("TTBAR_PDF","met X time TTBAR",RooArgSet(epfMet,timePhotReco), hTTBAR );


 
   // gamma + jet
   TH2F* h2_PJ_  = (TH2F*)fvec_pj[0]->Get(hfitname);
   TH2F* h2_PJ = (TH2F*)h2_PJ_->Clone("h2_PJ");
   protectZeroBin(h2_PJ);
   h2_PJ->Sumw2();
 

   RooDataHist hPJ("hPJ","hist PJ",RooArgSet(epfMet,timePhotReco), Import(*h2_PJ,false) );
   TH1* methistpj = hPJ.createHistogram("methistpj",epfMet);
   dumpTH1(methistpj);

   RooHistPdf PJ_PDF("PJ_PDF","met X time PJ",RooArgSet(epfMet,timePhotReco), hPJ ); 

   printf("GMSB: %f   QCD: %f     PJ: %f\n",h2_gmsb->Integral(),h2_QCD->Integral(),h2_PJ->Integral());
   printf("DATA %f\n",h2_data->Integral());


   double scaleQCD = qcdscale*h2_data->Integral()/(h1_QCD->Integral());
   double scaleG = gjetscale*h2_data->Integral()/h1_PJ->Integral();


   h2_QCD->Scale(scaleQCD);
   h2_PJ->Scale(scaleG);


 
   TH2F* h2_FakeMET = (TH2F*) h2_PJ->Clone("h2_FakeMET");
   h2_FakeMET->Add(h2_QCD);


   RooDataHist hFakeMET("hFakeMET","hist FakeMET",RooArgSet(epfMet,timePhotReco), Import(*h2_FakeMET,false) );
   TH1* methistFakeMET = hFakeMET.createHistogram("methistFakeMET",epfMet);
   dumpTH1(methistFakeMET);

   RooHistPdf FakeMET_PDF("FakeMET_PDF","met X time FakeMET",RooArgSet(epfMet,timePhotReco), hFakeMET ); 


    TH2F* h2_sum = (TH2F*) h2_PJ->Clone("h2_FakeMET");
    h2_sum->Add(h2_ttbar);
    h2_sum->Add(h2_Wenu);



    RooDataHist h_allMC("h_allMC","sum of all bkg samples MC and control sample", 
                         RooArgSet(epfMet,timePhotReco), Import(*h2_sum,false) );
    RooHistPdf BKG_PDF("BKG_PDF", "pdf for all bkg", RooArgSet(epfMet, timePhotReco), h_allMC);
 

   // std::cout<<"all MC"<<std::endl;

 
   RooRealVar lumi("lumi","luminosity",0., lumin,"pb^{-1}");
   RooRealVar sigEff("sigEff","signal effciency", 0., 1.);
   RooRealVar sigAcc("sigAcc","signal acceptance",0., 1.);
   RooRealVar xsection("xsection","cross section", 0.1, 0., 10., "pb");
  
 
   RooRealVar toy_met("toy_met", " GeV/c^{2}", 5., 8., 800.);
   RooRealVar toy_time("toy_time", " ns", 0.1, -2., 15.);
   RooDataSet* toy_mc_data =  BKG_PDF.generate(RooArgSet(epfMet, timePhotReco), h2_data->Integral());
   RooDataHist toy_mc_data_hist("toy_mc_data_hist", "toy_mc_data_hist",RooArgSet(epfMet, timePhotReco), *toy_mc_data );

   lumi.setVal( lumin );
   sigEff.setVal( efficiency[GMSBindex] );
   sigAcc.setVal( 1. );

   xsection.setVal(h2_gmsb->Integral()/(lumin*efficiency[GMSBindex]));
    
   lumi.setConstant();
   sigEff.setConstant();
   sigAcc.setConstant();

  // define xsection as a function of lumi, eff and #signal events
   RooFormulaVar nsig("nsig","# signal events","@0*@1*@2*@3",RooArgList(lumi,sigAcc,sigEff,xsection) );
   RooRealVar nttbar("nttbar","nttbar",1,100000);
   RooRealVar nWenu("nWenu","nWenu",1,100000);
   RooRealVar nbkg("nbkg","nbkg",1,1000000);
   RooRealVar nQCD("nQCD","# of QCD",1,1000000);
   RooRealVar nPJ("nPJ","# gamma+jet",1,1000000);
   RooRealVar nFakeMET("nFakeMET","# FakeMET",1,1000000);


   // fraction of bkg in gamma+jet and QCD. this should be coming from MC NOT from control sample!
   //RooRealVar nfrac("nfrac", "nfrac", 0., 1.);
   //nfrac.setVal(0.335);
   //nfrac.setConstant();
   
   // bkg PDF as sum of  all backgrounds
   RooAddPdf BKGFIX("bkgfix", "bkgfix", RooArgList(FakeMET_PDF,TTbar_PDF, Wenu_PDF ), RooArgList(nFakeMET, nttbar, nWenu) );
//   RooAddPdf BKG("bkg", "bkg", RooArgList(PJ_PDF,QCD_PDF,TTbar_PDF, Wenu_PDF ), RooArgList(nPJ, nQCD, nttbar, nWenu) );

   // total expected background with each component properly normalized
   //nbkg.setVal((h2_PJ->Integral()+h2_QCD->Integral())*1.15);
   //std::cout<<h2_PJ->Integral()+h2_QCD->Integral()<<" nbkg"<<std::endl;

   nPJ.setVal(h2_PJ->Integral());
   nQCD.setVal(h2_QCD->Integral());
   nFakeMET.setVal(h2_FakeMET->Integral());


   // ttbar  and Wenu normnalization. Now using theory xsec but should switch to CMS measurement
   // 
   nttbar.setVal(h2_ttbar->Integral());
   nttbar.setConstant();

   // assume MC is correctly scaled and is correct
   nWenu.setVal(h2_Wenu->Integral());
   nWenu.setConstant();

/*
   std::cout<< "tt: " << nttbar.getVal()<< "wenu: " << nWenu.getVal()<< "gmsb: " << 
               nsig.getVal()<< "QCD " << nQCD.getVal()<< "PJ " << nPJ.getVal()<< 
     " DATA: " <<h2_data->Integral()<< " <-DATA" << std::endl;
*/

   RooAddPdf model("model","Signal + Background",
                   RooArgList(signal,FakeMET_PDF,TTbar_PDF, Wenu_PDF),
                   RooArgList(nsig,nFakeMET, nttbar, nWenu));
  
   // pointer to dataset
  RooDataHist& data = h_allDATA;
 //  RooDataHist& data =toy_mc_data_hist;  // to obtain the observed switch to ---> h_allDATA
  
   
   // LL fit
   RooNLLVar nll("nll","nll",model,data,RooFit::SumW2Error(kTRUE) ,RooFit::Extended(kTRUE));
   RooMinuit mnll(nll);
   mnll.migrad() ;
 //  mnll.hesse() ;
 //  mnll.minos() ;
   
   RooFitResult* rfit = mnll.fit("r");

   rfit->Print("v");

   std::cout << "Fitted signal yield: " << nsig.getVal()
         //    << " +/- " << bla.getError() 
             << std::endl;




   ul.push_back(xsection.getError());

   // plot of roo data hist
   RooPlot *xframe = epfMet.frame(100);
   xframe->SetXTitle("Missing energy [GeV]");
   data.plotOn(xframe,Invisible());
   
  
   model.plotOn(xframe, LineWidth(1.5));
   model.plotOn(xframe,Components(signal),LineColor(kBlue)) ;
   //BKGFIX.plotOn(xframe, LineWidth(2.5), LineColor(27), LineStyle(kDashed));
   BKGFIX.plotOn(xframe, DrawOption("F"), FillColor(16));
   model.plotOn(xframe,LineWidth(2),Components(signal),LineColor(kRed), Name("signal")) ;

   RooPlot *yframe = timePhotReco.frame(100);
   yframe->SetXTitle("ECAL time [ns]");
   data.plotOn(yframe, Invisible() );
 

   model.plotOn(yframe, LineWidth(1.5));
   model.plotOn(yframe,Components(signal),LineColor(kBlue)) ;
   //BKGFIX.plotOn(yframe, LineWidth(2.5), LineColor(27), LineStyle(kDashed));
   BKGFIX.plotOn(yframe, DrawOption("F"), FillColor(16));
   model.plotOn(yframe,LineWidth(2), Components(signal),LineColor(kRed), Name("signal")) ;
  
 
  int floatParms = (rfit->floatParsFinal()).getSize();
  RooHist* pullhisto = xframe->pullHist();
  double probchi = TMath::Prob((xframe->chiSquare(floatParms)) * (pullhisto->GetN()-floatParms), (pullhisto->GetN()-floatParms));
  std::cout << "chisq probability = " << probchi << std::endl;
  std::cout << "float params = " << floatParms << std::endl;

  int DataMarkerStyle = 20;
  double DataMarkerSize = 1.5;

   TH1F* h_leg_data= new TH1F("h_leg_data", "", 100, 0, 100);
   h_leg_data->SetMarkerStyle(DataMarkerStyle);
   h_leg_data->SetMarkerSize(DataMarkerSize);
   h_leg_data->SetMarkerColor(kBlack);

   TH1F* h_leg_model= new TH1F("h_leg_model", "", 100, 0, 100);
   h_leg_model->SetLineWidth(2);
   h_leg_model->SetLineColor(kBlue);
   
   TH1F* h_leg_signal= new TH1F("h_leg_signal", "", 100, 0, 100);
   h_leg_signal->SetLineWidth(2);
   h_leg_signal->SetLineColor(kRed);

   TH1F* h_leg_bkg= new TH1F("h_leg_bkg", "", 100, 0, 100);
   //h_leg_bkg->SetLineWidth(2.5);
   //h_leg_bkg->SetLineStyle(kDashed);
   h_leg_bkg->SetFillColor(16);

   TLegend* leg = new TLegend(0.5,0.65,0.85,0.93);
   leg->SetFillColor(0);
   leg->SetBorderSize(0);
   leg ->SetTextFont(22);
   leg->AddEntry(h_leg_data,"Data","PE");
   leg->AddEntry(h_leg_model,"Signal + Background","l");
   leg->AddEntry(h_leg_bkg,"All background","f");
   leg->AddEntry(h_leg_signal,"GMSB ("+TString(ctau)+")","l");


   TCanvas c2("c2","log scale canvas",1024,768);
   c2.SetLogy();
   //c2.SetLogx();
   xframe->Draw();
   xframe->SetMinimum( 0.0001) ;
   xframe->SetMaximum( xframe->GetMaximum()*2.) ;
   met_pdf_hist->SetLineWidth(1.5); 
   met_pdf_hist->Draw("PESAME");
   leg->Draw("same");
   c2.SaveAs("plot_fit_met_data_log.eps");


    yframe->Draw();
    yframe->SetMinimum( 0.0001) ;
    yframe->SetMaximum( yframe->GetMaximum()*2.) ;
    time_pdf_hist->SetLineWidth(1.5); 
    time_pdf_hist->Draw("PESAME");
    leg->Draw("same");
    c2.SaveAs("plot_fit_time_data_log.eps");


   
   delete(xframe);
   delete(yframe);
  
   // mni LL of the fit
   double minNLL = nll.getVal();

   
   // scan LL vs. x section
   TGraph* gr_nll = 0;
    TH1* h_nll = 0;
   cout << "==== start scan of LL vs. x section" << endl; 
   const size_t nbins = 100;
   double llx[nbins], lly[nbins];
   double xmin(xsection.getVal()-5.*xsection.getError()), xmax(xsection.getVal()+5.*xsection.getError());
   h_nll = new TH1F("nllScan","Scan of LL vs. cross section", nbins,xmin,xmax);
   
  

   cout << "Scanning " << xsection.GetName() << endl;

   //xsection.printToStream(cout);
   cout <<"Scan interval [" << xmin << ":" << xmax << "]" << endl;
   xsection.setConstant(kTRUE);
   
   //RooRealVar* k1 = (RooRealVar*) r2D->floatParsFinal().find("\\kappa"); k1->setConstant();
   //RooRealVar* pde1 = (RooRealVar*) r2D->floatParsFinal().find("p1De"); pde1->setConstant();
   
   for (size_t i=0 ; i<nbins ; i++) {
     double x = xmin + (xmax-xmin)*i/nbins;
     
     llx[i] = x;
     lly[i] = 0.;
     
     xsection.setVal(x) ;
     //cout << "i: " << i << "\tx: " << x << endl;
     mnll.migrad() ;
     Double_t nllScanVal = nll.getVal() ;

     h_nll->SetBinContent(i+1,TMath::Exp(-nllScanVal+minNLL)) ;
     lly[i] = TMath::Exp(-nllScanVal+minNLL);
      
   }
   
   gr_nll = new TGraph(nbins,llx,lly);
   char nomefile[100];
   sprintf(nomefile,"%s%i%s","graphGMSB_INDEX_",GMSBindex,".root");
   TFile* fout = TFile::Open(nomefile,"RECREATE");
   gr_nll->Write();
   h_nll->Write();
   
   fout->Write();
   fout->Close();
 
   
   return gr_nll;


}

double normGauss(double mu,double sig) {
  double x = (1./sqrt(2.*TMath::Pi()*sig*sig))*exp(-(mu/sig)*(mu/sig)/2.);
  return x;
}

TGraph* convoluteGraph(const TGraph* ingr, const double& smear) {
  double* x = ingr->GetX();
  double* y = ingr->GetY();
  const int nb = ingr->GetN();

  double newx[nb], newy[nb];
  for(Int_t i=0; i<nb; ++i) { // loop over bins of new graph
    newy[i] = 0.0;
    newx[i] = x[i];
    for(Int_t j=0; j<nb; ++j) { // contribution from all bins of orig graph
      double dx  = x[i]- x[j]; // x(new) - x(old)
      double weight = x[i]>=0 ? normGauss(dx,smear) : 0.;
      double dy = weight*y[j];
      newy[i] += dy;

      if(0==1) {
        cout << "new bin: " << newx[i] << "  oldbin:" << x[j]
             << "  dx:" << dx << "  weight: " << weight 
             << " oldy: " << y[j] << endl;
      }
    }
  }
  TGraph* newg = new TGraph(nb,newx,newy);
  return newg;

}




void plot_limit(std::string ctau ){


/*Float_t upperlimit_Lambda100[5] = {0.018,0.022,0.0438671,0.054,0.074};
Float_t upperlimit_Lambda120[5] = {0.01,0.018,0.035048,0.046,0.062};
Float_t upperlimit_Lambda140[5] = {0.006,0.01,0.0240775,0.03,0.042};
Float_t upperlimit_Lambda160[5] = {0.002,0.006,0.0177448,0.022,0.03};
Float_t upperlimit_Lambda180[5] = {0.002,0.006,0.0131662,0.014,0.018};
*/
        Float_t upperlimit1_Lambda100[5] = {0.,0.,0.,0.,0.};
        Float_t upperlimit_Lambda120[5] = {0.0471724,0.0522344,0.0600325,0.0670999,0.0778567};
        Float_t upperlimit_Lambda140[5] = {0.006,0.01,0.0255327,0.03,0.038};
        Float_t upperlimit_Lambda160[5] = {0.006,0.006,0.0187741,0.022,0.03};
        Float_t upperlimit_Lambda180[5] = {0.002,0.006,0.0138085,0.014,0.018};
        Float_t upperlimit_Lambda100[5] = {0.0230419,0.0250243,0.0288547,0.0319195,0.0393824};
        Float_t upperlimit500_Lambda100[5] = {0.0248624,0.0274653,0.0306991,0.0333265,0.0433285};
        Float_t upperlimit1000_Lambda100[5] = {0.0314197,0.034027,0.0367808,0.0392051,0.0422995};
        Float_t upperlimit2000_Lambda100[5] = {0.0462303,0.0491389,0.052771,0.0559793,0.0681559};
        Float_t upperlimit4000_Lambda100[5] = {0.0745266,0.0841812,0.0891725,0.0960313,0.100699};
        Float_t upperlimit6000_Lambda100[5] = {0.10191,0.110137,0.119908,0.132735,0.154897};


   Double_t obs_lim[5]= {0.};
   Double_t exp_lim[5]= {0.};



/*	obs_lim[0] =0.0612122; 
	obs_lim[1] =0.0460931;
	obs_lim[2] =0.0379412; 
	obs_lim[3] =0.0307229; 
	obs_lim[4] =0.0241807; 

*/

obs_lim[0] =0.0549363;
obs_lim[1] =0.0445393;
obs_lim[2] =0.0332381;
obs_lim[3] =0.0283376;
obs_lim[4] =0.0249918;


   exp_lim[0] = upperlimit_Lambda100[2]; 
   exp_lim[1] = upperlimit_Lambda120[2]; 
   exp_lim[2] = upperlimit_Lambda140[2]; 
   exp_lim[3] = upperlimit_Lambda160[2]; 
   exp_lim[4] = upperlimit_Lambda180[2]; 


   Double_t theoFac = 1.39;

   Double_t mTh[5] = {  100*theoFac,
                        120*theoFac,
                        140*theoFac,
                        160*theoFac,
                        180*theoFac};

 Double_t xsTh[5] = {   0.2357,
                        0.0860,
                        0.0368,
                        0.0181,
                        0.0092 };

 Double_t x_pdf[10] = { 100*theoFac,
                        120*theoFac,
                        140*theoFac,
                        160*theoFac,
                        180*theoFac,
                        180*theoFac,
                        160*theoFac,
                        140*theoFac,
                        120*theoFac,
                        100*theoFac };


 Double_t y_pdf[10];


 for(int i = 0; i<5; i++){
   y_pdf[i] = 0.07*xsTh[i]+xsTh[i];
 }

 int j = 3;
 for(int i = 5; i<10; i++){
   y_pdf[i] = -0.07*xsTh[j]+xsTh[j];
   j--;
 }


Double_t y_pdf_exp[10];

 for(int ii = 0; ii<5; ii++){
   y_pdf_exp[ii] = 0.2*obs_lim[ii]+obs_lim[ii];
 }
 
 int jj = 3;

 for(int ii = 5; ii<10; ii++){
   y_pdf_exp[ii] = -0.2*obs_lim[jj]+obs_lim[jj];
   jj--;
 }


  // Green and yellow error bands
 
   Double_t y_pdf_1sig[10] = {  upperlimit_Lambda100[1],
				upperlimit_Lambda120[1],
				upperlimit_Lambda140[1],
				upperlimit_Lambda160[1],
				upperlimit_Lambda180[1],
				upperlimit_Lambda180[3],
				upperlimit_Lambda160[3],
				upperlimit_Lambda140[3],
				upperlimit_Lambda120[3],
				upperlimit_Lambda100[3] 
				};


   Double_t y_pdf_2sig[10] = {  upperlimit_Lambda100[0],
				upperlimit_Lambda120[0],
				upperlimit_Lambda140[0],
				upperlimit_Lambda160[0],
				upperlimit_Lambda180[0],
				upperlimit_Lambda180[4],
				upperlimit_Lambda160[4],
				upperlimit_Lambda140[4],
				upperlimit_Lambda120[4],
				upperlimit_Lambda100[4] 
				};

 /*
   Double_t y_pdf_1sig[10] = {	0.03, 
				0.034,
				0.026, 
				0.022, 
				0.018, 
				0.066, 
				0.07, 
				0.066, 
				0.058, 
				0.046 };

   Double_t y_pdf_2sig[10] = {	0.0022, 
				0.0026, 
				0.0018, 
				0.0018, 
				0.0014, 
				0.09, 
				0.094, 
				0.086, 
				0.086, 
				0.066 };
*/


   TGraph* Onesig_graph;
   Onesig_graph = new TGraph(10., x_pdf, y_pdf_1sig );

   TGraph* Twosig_graph;
   Twosig_graph = new TGraph(10., x_pdf, y_pdf_2sig );

   Onesig_graph->SetFillColor(kGreen);
   Twosig_graph->SetFillColor(kYellow);


   TGraph* exp_lim_graph;
   exp_lim_graph  = new TGraph(5, mTh, exp_lim);
   exp_lim_graph->SetLineWidth(2.5);
   exp_lim_graph->SetLineStyle(2);

   TGraph* ul_lim_graph;
   ul_lim_graph  = new TGraph(5, mTh, obs_lim);
   ul_lim_graph->SetLineColor(kBlack);
   ul_lim_graph->SetLineWidth(2);


   TMultiGraph* mg = new TMultiGraph;
  
   mg->Add(exp_lim_graph);
   mg->Add(ul_lim_graph);


  
   TCanvas* c0 = new TCanvas("exclusion limit", "exclusion limit", 1);
   c0->cd();
   c0->SetLogy();
   // c0->SetGridx();
   //c0->SetGridy();
  

// region excluded by Tevatron limits

   Double_t x_shaded[8] = {	130,
				135,
				140,
				146,
				146, 
				140, 
				135, 
				130 };


   Double_t y_shaded[8] = {	0.0007,
				0.0007,
				0.0007,
				0.0007,
				5, 
				5, 
				5, 
				5 };

   TGraph *pl = new TGraph(8,x_shaded,y_shaded);
   pl->SetFillColor(kGray);
   pl->Draw("f");
   mg->Draw("Lsame");
   Twosig_graph->Draw("Fsame");
   Onesig_graph->Draw("Fsame");
   mg->Draw("Lsame");
 
   mg->GetXaxis()->SetTitle("M_{#tilde{#chi^{0}_{1}}} [GeV/c^{2}]");
   mg->GetYaxis()->SetTitle("#sigma ( #tilde{#chi^{0}_{1}} #rightarrow #tilde{G}#gamma) #times BR [pb]");

   std::string s_lumi;
   s_lumi = "4.55";
   std::string lint = "#intLdt= "+s_lumi+" fb^{-1}";
   TLatex l1;
   l1.SetTextAlign(12);
   l1.SetTextSize(0.035);
   l1.SetTextFont(22);
   l1.SetNDC();
   l1.DrawLatex(0.155, 0.967, "CMS Preliminary");
   l1.DrawLatex(0.75, 0.96, lint.c_str());

   
   TGraph *xsTh_vs_m = new TGraph(5, mTh, xsTh);
   xsTh_vs_m->SetLineWidth(2);
   xsTh_vs_m->SetLineColor(kRed);
   xsTh_vs_m->SetMarkerSize(1.);
   xsTh_vs_m->SetMarkerStyle(22);
   xsTh_vs_m->SetMarkerColor(kRed);
   xsTh_vs_m->Draw("Csame");
   xsTh_vs_m->SetTitle("M_{#tilde{#chi^{0}_{1}}} [GeV/c^{2}]");
   xsTh_vs_m->GetXaxis()->SetTitleOffset(1.5);
   xsTh_vs_m->SetTitle("#sigma (pb) ");
   mg->GetYaxis()->SetRangeUser(1e-3, 10.);
   mg->GetXaxis()->SetRangeUser(140., 248);



   TLegend* leg = new TLegend(0.5,0.6,0.95,0.9);
   leg->SetFillStyle(0); leg->SetBorderSize(0); 
   leg->SetFillColor(0);
 
   leg->SetHeader(("#tilde{#chi}^{0}_{1} c#tau = "+ctau+" mm").c_str());
   leg ->SetTextFont(22);
   leg->AddEntry(pl,"CDF exclusion (2.6 fb^{-1})","f");
   leg->AddEntry(xsTh_vs_m,"Theoretical LO cross-section","l");
   leg->AddEntry(ul_lim_graph, "Observed  95% CL upper limit", "L");
   leg->AddEntry(exp_lim_graph, "Expected 95% CL upper limit", "L");
   leg->AddEntry(Onesig_graph, "#pm 1 #sigma Expected", "F");
   leg->AddEntry(Twosig_graph, "#pm 2 #sigma Expected", "F");
   leg->Draw("same");

   c0->SaveAs(("exclusion_limit_"+ctau+".eps").c_str());
   //c0->SaveAs(("exclusion_limit_"+ctau+".pdf").c_str());
   //c0->SaveAs(("exclusion_limit_"+ctau+".png").c_str());

}





#include "tdrstyle.C"

void rootlogon() {
   gROOT->SetStyle("Plain");
   gStyle->SetPalette(1);
   gStyle->SetOptStat(1111111);  // Show overflow, underflow + SumOfWeights 
   gStyle->SetOptFit(111110);
   gStyle->SetOptFile(1);
   gStyle->SetMarkerStyle(20);
   gStyle->SetMarkerSize(2.);
   gStyle->SetMarkerColor(1);
   gStyle->SetTitleOffset(1.20,"Y");

    //define high def color palette
    const Int_t NRGBs = 5;
    const Int_t NCont = 255;

    Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
    Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
    Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
    Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
    TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
    gStyle->SetNumberContours(NCont);

    cout << "loading TDR style and setting as default" << endl;
    gROOT->ProcessLine(".L tdrstyle.C");
    setTDRStyle();

}

