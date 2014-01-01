// $Id: fitAppToy.cc,v 1.2 2012/02/02 15:50:47 sigamani Exp $
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
void band_calculation(std::string point);

double normGauss(double mu,double sig);
TGraph* convoluteGraph(const TGraph* ingr, const double& smear);

double computeUpperLimit(TGraph* hist, double xCL);
double sigmax(const TH1F* hin, double num);

TGraph* fit_2D( int GMSBindex, float lumin,  std::string ctau);
int Nbin;


void plot_limit_lifetime(std::string LAMBDA);
void plot_limit_2d( );



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

 }  
  return hout;  
}



int main(int argc, char* argv[]) {

  rootlogon();

  std::string  sys = argv[1];
  std::string  LAMBDA = argv[2];
  std::string  CTAU = argv[3];
  int  ntoy = atoi(argv[4]);
  Nbin=ntoy;

 // dir to data files
  TString dirpath("/u1/sigamani/fit-1-2-12/"+sys+"/");
 
 // load one file for each sample
  fgmsb.push_back( (TFile*)  TFile::Open(dirpath+"PhotonJet_GMSB_Lambda-"+LAMBDA+"_CTau-"+CTAU+"_7TeV_pythia6_cff-Summer11-PU_S4_START42_V11-v1-AODSIM-UNCLEAN21_pfakt5_2ndJet10_"+sys+".root"));
  fvec.push_back( (TFile*) TFile::Open(dirpath+"PhotonJet_DATA_2011_pfakt5_CONTROL_QCD_2ndJet10_"+sys+".root"));
  fvec_pj.push_back( (TFile*) TFile::Open(dirpath+"PhotonJet_DATA_2011_pfakt5_CONTROL_G_2ndJet10_"+sys+".root"));
  fvec_data.push_back( (TFile*) TFile::Open(dirpath+"PhotonJet_DATA_2011_pfakt5_2ndJet10_"+sys+".root"));
  fvec_ttbar = (TFile*) TFile::Open(dirpath+"PhotonJet_TT_2011_pfakt5_2ndJet10_"+sys+".root");
  fvec_Wenu = (TFile*) TFile::Open(dirpath+"PhotonJet_WToENU_2011_pfakt5_2ndJet10_"+sys+".root");

 if (CTAU=="1") { efficiency.push_back(0.204) ; } //not present
 if (CTAU=="250") { efficiency.push_back(0.191767) ; }
 if (CTAU=="500") { efficiency.push_back(0.170677) ; }

 if (CTAU=="1000" && sys=="default") { efficiency.push_back(0.131104);}
 if(CTAU=="1000" && sys=="JESM_02") { efficiency.push_back(0.129658) ; }
 if (CTAU=="1000" && sys=="JESP_02") { efficiency.push_back(0.133) ; }
if(CTAU=="1000" && sys=="JESM_04") { efficiency.push_back(0.128212) ; }
 if (CTAU=="1000" && sys=="JESP_04") { efficiency.push_back(0.134478) ; }
if(CTAU=="1000" && sys=="JESM_06") { efficiency.push_back(0.126766) ; }
 if (CTAU=="1000" && sys=="JESP_06") { efficiency.push_back(0.135924) ; } 
 if (CTAU=="1000" && sys=="JER")  { efficiency.push_back(0.13054)  ; }
 if (CTAU=="1000" && sys=="PESM") { efficiency.push_back(0.12773) ; }
 if (CTAU=="1000" && sys=="PESP") { efficiency.push_back(0.130622)  ; }

 if (CTAU=="1000" && sys=="MESM") { efficiency.push_back(0.131104) ; }
 if (CTAU=="1000" && sys=="MESP") { efficiency.push_back(0.131104) ; }
 if (CTAU=="1000" && sys=="MER")  { efficiency.push_back(0.131104) ; }
 if (CTAU=="1000" && sys=="TESM") { efficiency.push_back(0.131104) ; }
 if (CTAU=="1000" && sys=="TESP") { efficiency.push_back(0.131104) ; }


 if (CTAU=="2000") { efficiency.push_back(0.086266) ; }
 if (CTAU=="3000") { efficiency.push_back(0.220) ; }//not present
 if (CTAU=="4000") { efficiency.push_back(0.048732) ; }
 if (CTAU=="6000") { efficiency.push_back(0.0346424) ; }


 
  // contains the limit graph
  double UL;
  TH1F* h_UL = new TH1F("h_UL", "h_UL", 500., 0., 2.);

  //loop for toy Mc experiment  
   for(int j = 0;j< ntoy;j++){

     TGraph *t=0;
  
     t = fit_2D( 0, 4555.6, CTAU ); // take values from command line
     UL = computeUpperLimit(t,0.95);
     if (UL <1e-4) continue;

     h_UL->Fill(UL);
   
     }

	
    

  TCanvas* c = new TCanvas ("c", "c", 1000);

  c->cd();
  h_UL->Draw("hist");
  c->SaveAs("toy_mc_UL_CT"+TString(CTAU)+"_Lambda"+TString(LAMBDA)+".png"); // change name if change GMSB signal file
  
  TFile* f_out = new TFile("f_out_CT_" + TString(CTAU) + "_Lambda" + TString(LAMBDA) + ".root", "RECREATE"); 

  f_out->cd();
  h_UL->Write();
  f_out->Write();

  band_calculation(CTAU+"_Lambda"+LAMBDA);  
  plot_limit_lifetime("100");


  return 0;

} 

void band_calculation(std::string point){

  TFile* f_in = (TFile*) TFile::Open(("f_out_CT_"+point+".root").c_str());
  TH1F* h1 = (TH1F*)f_in->Get("h_UL");
  h1->Scale(1./h1->Integral());
  std::ofstream f_out (("limit_"+point+".txt").c_str());


  f_out << "Float_t upperlimit" << TString(point) << "[5] = {" << sigmax(h1,0.023) 
       << "," << sigmax(h1,0.13623) << "," << h1->GetMean() << "," 
       << sigmax(h1,0.81823) << "," << sigmax(h1,0.95423) << "};" << endl;
				 
}


double sigmax(const TH1F* hin, double num){

          double result;
	  for(int i = 1;i<=hin->GetNbinsX();i++){
	
	  double integral = hin->Integral(0, i);
          if ( integral > num ) break;
          result = hin->GetBinCenter(i);
	  }

 return result;

}

//===================================
double computeUpperLimit(TGraph* hist, double xCL) {
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

  double UL;
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
      UL = gx[i];
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
  //c1->SaveAs("scan-upperlimit.ps");
  //c1->SaveAs("scan-upperlimit.gif");
  //c1->SaveAs("scan-upperlimit.eps");
  
  //delete c1;
  ul.push_back( newx[iCL] );
  return UL;
}


//--------------
TGraph* fit_2D( int GMSBindex, float lumin,  std::string ctau){
//--------------
   char treeName[100];
   sprintf(treeName,"jetTree");

 
   
     //---------------import data--------------//

   cout << "--- dataset for DATA --- " << endl;

   TString hfitname("met_varsize_x_timePhot_varsize");

   TH2F* h2_data= (TH2F*)fvec_data[0]->Get(hfitname)->Clone("h2_data");
   protectZeroBin(h2_data);
  
   RooRealVar epfMet("met","PF MET",h2_data->ProjectionX()->GetXaxis()->GetXmin(), h2_data->ProjectionX()->GetXaxis()->GetXmax(), "GeV");
   RooRealVar timePhotReco("time","ECAL time",h2_data->ProjectionY()->GetXaxis()->GetXmin(), h2_data->ProjectionY()->GetXaxis()->GetXmax(),"ns");

   RooDataHist h_allDATA("h_allDATA","hist allDATA",RooArgSet(epfMet,timePhotReco), Import(*h2_data,false) );
   RooHistPdf h_allDATA_PDF("h_allDATA_PDF","met X time data",RooArgSet(epfMet,timePhotReco), h_allDATA ); 
   
   std::cout<<"RooDataHist created"<<std::endl;


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
   // dumpTH1(methistwenu);

   RooHistPdf Wenu_PDF("Wenu_PDF","met X time Wenu",RooArgSet(epfMet,timePhotReco), hWenu );

  
   
   // QCD control sample
   std:: cout << "--- dataset for QCD --- " << endl;

   // original histo
   TH2F* h2_QCD_ = (TH2F*)fvec[0]->Get(hfitname);
   
   
   // new histo to get rid of empty bins causing trouble with LL fit
 

   TH2F* h2_QCD = (TH2F*)h2_QCD_->Clone("h2_QCD");
   protectZeroBin(h2_QCD);
   
   
   // scale control sample and MC
   h2_QCD->Sumw2();
 
  

   RooDataHist hQCD("hQCD","hist QCD",RooArgSet(epfMet,timePhotReco), Import(*h2_QCD,false) );
   TH1* methistqcd = hQCD.createHistogram("methistqcd",epfMet);
   //  dumpTH1(methistqcd);
 
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
   // std::cout<<"signal pdf"<<std::endl;


   // ttbar PDF
   TH2F* h2_ttbar_ = (TH2F*)fvec_ttbar->Get(hfitname);
   TH2F* h2_ttbar = (TH2F*)h2_ttbar_->Clone("h2_ttbar");
   protectZeroBin(h2_ttbar);

   h2_ttbar->Sumw2();
   h2_ttbar->Scale(lumin);
   
   RooDataHist hTTBAR("hTTBAR","hist TTBAR",RooArgSet(epfMet,timePhotReco), Import(*h2_ttbar,false) );
   TH1* methistttbar = hTTBAR.createHistogram("methistttbar",epfMet);
   // dumpTH1(methistttbar);

   RooHistPdf TTbar_PDF("TTBAR_PDF","met X time TTBAR",RooArgSet(epfMet,timePhotReco), hTTBAR );


 
   // gamma + jet
   TH2F* h2_PJ_  = (TH2F*)fvec_pj[0]->Get(hfitname);
   TH2F* h2_PJ = (TH2F*)h2_PJ_->Clone("h2_PJ");
   protectZeroBin(h2_PJ);
   h2_PJ->Sumw2();
 


   RooDataHist hPJ("hPJ","hist PJ",RooArgSet(epfMet,timePhotReco), Import(*h2_PJ,false) );
   TH1* methistpj = hPJ.createHistogram("methistpj",epfMet);
   // dumpTH1(methistpj);

   RooHistPdf PJ_PDF("PJ_PDF","met X time PJ",RooArgSet(epfMet,timePhotReco), hPJ ); 

/*
   h2_QCD->Scale(38520./57229.);
   h2_PJ->Scale(73027./95872.2);
*/

   h2_QCD->Scale(0.7/(h2_QCD->Integral()));
   h2_PJ->Scale(0.3/(h2_PJ->Integral()));




   TH2F* h2_FakeMET = (TH2F*) h2_PJ->Clone("h2_FakeMET");
   h2_FakeMET->Add(h2_QCD);


   RooDataHist hFakeMET("hFakeMET","hist FakeMET",RooArgSet(epfMet,timePhotReco), Import(*h2_FakeMET,false) );
   TH1* methistFakeMET = hFakeMET.createHistogram("methistFakeMET",epfMet);
   dumpTH1(methistFakeMET);

   RooHistPdf FakeMET_PDF("FakeMET_PDF","met X time FakeMET",RooArgSet(epfMet,timePhotReco), hFakeMET );


   
   //sum of all bkg
   TH2F* h2_sum = (TH2F*) h2_ttbar->Clone("h2_sum");
   // h2_sum->Add(h2_QCD);
    //h2_sum->Add(h2_PJ);
    h2_sum->Add(h2_FakeMET);
    h2_sum->Add(h2_Wenu);


    RooDataHist h_allMC("h_allMC","sum of all bkg samples MC and control sample", 
                         RooArgSet(epfMet,timePhotReco), Import(*h2_sum,false) );
    RooHistPdf BKG_PDF("BKG_PDF", "pdf for all bkg", RooArgSet(epfMet, timePhotReco), h_allMC);

   
 
   RooRealVar lumi("lumi","luminosity",0., lumin,"pb^{-1}");
   RooRealVar sigEff("sigEff","signal effciency", 0., 1.);
   RooRealVar sigAcc("sigAcc","signal acceptance",0., 1.);
   RooRealVar xsection("xsection","cross section", 0.1, 0., 10, "pb");
  

   lumi.setVal( lumin );
   sigEff.setVal( efficiency[GMSBindex] );
   sigAcc.setVal( 1. );

   xsection.setVal(h2_gmsb->Integral()/(lumin*efficiency[GMSBindex]));
    
   lumi.setConstant();
   sigEff.setConstant();
   sigAcc.setConstant();

  // define #signal events as a function of lumi, eff and xsec
   RooFormulaVar nsig("nsig","# signal events","@0*@1*@2*@3",RooArgList(lumi,sigAcc,sigEff,xsection) );
   RooRealVar nttbar("nttbar","nttbar",1,1000000);
   RooRealVar nWenu("nWenu","nWenu",1,1000000);
   RooRealVar nbkg("nbkg","nbkg",1,10000000);
   RooRealVar nQCD("nQCD","# of QCD",1,10000000);
   RooRealVar nPJ("nPJ","# gamma+jet",1,10000000);
   RooRealVar nFakeMET("nFakeMET","# FakeMET",1,1000000);


   // fraction of bkg in gamma+jet and QCD. this should be coming from MC NOT from control sample!
   //RooRealVar nfrac("nfrac", "nfrac", 0., 1.);
   //nfrac.setVal(0.335);
   //nfrac.setConstant();
   
   // bkg PDF as sum of  all backgrounds
   //RooAddPdf BKG("bkg", "bkg", RooArgList(PJ_PDF,QCD_PDF),nfrac);
   RooAddPdf BKG("bkg", "bkg", RooArgList(PJ_PDF,QCD_PDF,TTbar_PDF, Wenu_PDF ), RooArgList(nPJ, nQCD, nttbar, nWenu) );
   RooAddPdf BKGFIX("bkgfix", "bkgfix", RooArgList(FakeMET_PDF,TTbar_PDF, Wenu_PDF ), RooArgList(nFakeMET, nttbar, nWenu) );

   // total expected background with each component properly normalized
   //nbkg.setVal((h2_PJ->Integral()+h2_QCD->Integral())*1.15);
   //std::cout<<h2_PJ->Integral()+h2_QCD->Integral()<<" nbkg"<<std::endl;

   //****************TOY MC*******************//

   RooRealVar toy_met("toy_met", " GeV/c^{2}", 5., 8., 800.);
   RooRealVar toy_time("toy_time", " ns", 0.1, -2., 15.);
   RooDataSet* toy_mc_data =  BKG_PDF.generate(RooArgSet(epfMet, timePhotReco), h2_data->Integral());
   RooDataHist toy_mc_data_hist("toy_mc_data_hist", "toy_mc_data_hist",RooArgSet(epfMet, timePhotReco), *toy_mc_data );


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



   RooAddPdf model("model","signal + bkg",
                   //RooArgList(signal,PJ_PDF,QCD_PDF,TTbar_PDF, Wenu_PDF),
                   //RooArgList(nsig,nPJ, nQCD, nttbar, nWenu));
                   RooArgList(signal,FakeMET_PDF,TTbar_PDF, Wenu_PDF),
                   RooArgList(nsig,nFakeMET, nttbar, nWenu));
  
   // pointer to dataset
   RooDataHist& data =toy_mc_data_hist;  // to obtain the observed switch to ---> h_allDATA
  
   
   // LL fit
   RooNLLVar nll("nll","nll",model,data,RooFit::SumW2Error(kTRUE) ,RooFit::Extended(kTRUE));
   RooMinuit mnll(nll);
   mnll.migrad() ;
   mnll.hesse() ;
   //mnll.minos() ;
   
   RooFitResult* rfit = mnll.fit("r");


   ul.push_back(xsection.getError());

//RooRealVar* par1_fitresult = (RooRealVar*) rfit->floatParsFinal()->find("par1"); 
//par1_fitresult->GetAsymErrorHi() ;
  
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
   
   for (size_t i=0 ; i<nbins ; i++) {
     double x = xmin + (xmax-xmin)*i/nbins;
     
     llx[i] = x;
     lly[i] = 0.;
     
     xsection.setVal(x) ;
     //cout << "i: " << i << "\tx: " << x << endl;
     //mnll.migrad() ;
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




void plot_limit_lifetime(std::string LAMBDA ){

	//INPUT TXT HERE

	Float_t upperlimit1_Lambda100[5] =    {0.018,0.026,0.0472708,0.058,0.078};
	Float_t upperlimit1_Lambda120[5] =    {0.018,0.026,0.0473782,0.062,0.082};
	Float_t upperlimit1_Lambda140[5] =    {0.018,0.026,0.0517175,0.066,0.086};
	Float_t upperlimit1_Lambda160[5] =    {0.014,0.022,0.0445252,0.054,0.074};
	Float_t upperlimit1_Lambda180[5] =    {0.014,0.018,0.0377363,0.046,0.062};
	Float_t upperlimit250_Lambda100[5] =  {0.022,0.03,0.0524332,0.066,0.09};
	Float_t upperlimit250_Lambda120[5] =  {0.026,0.034,0.0554107,0.07,0.094};
	Float_t upperlimit250_Lambda140[5] =  {0.018,0.026,0.0541957,0.066,0.086};
	Float_t upperlimit250_Lambda160[5] =  {0.018,0.022,0.0476099,0.058,0.086};
	Float_t upperlimit250_Lambda180[5] =  {0.014,0.018,0.0397341,0.046,0.066};
	Float_t upperlimit500_Lambda100[5] =  {0.026,0.034,0.0607699,0.078,0.114};
	Float_t upperlimit500_Lambda120[5] =  {0.026,0.034,0.060671,0.074,0.102};
	Float_t upperlimit500_Lambda140[5] =  {0.022,0.034,0.0654566,0.078,0.102};
	Float_t upperlimit500_Lambda160[5] =  {0.014,0.022,0.0509845,0.062,0.09};
	Float_t upperlimit500_Lambda180[5] =  {0.018,0.022,0.0442091,0.054,0.07};
	Float_t upperlimit1000_Lambda100[5] = {0.03,0.042,0.0725348,0.094,0.126};
	Float_t upperlimit1000_Lambda120[5] = {0.034,0.042,0.0788503,0.098,0.13};
	Float_t upperlimit1000_Lambda140[5] = {0.026,0.038,0.0771741,0.094,0.13};
	Float_t upperlimit1000_Lambda160[5] = {0.018,0.034,0.0636353,0.078,0.11};
	Float_t upperlimit1000_Lambda180[5] = {0.018,0.026,0.0550671,0.07,0.09};
	Float_t upperlimit2000_Lambda100[5] = {0.038,0.054,0.0925573,0.118,0.162};
	Float_t upperlimit2000_Lambda120[5] = {0.042,0.058,0.101117,0.13,0.17};
	Float_t upperlimit2000_Lambda140[5] = {0.042,0.054,0.0953233,0.122,0.17};
	Float_t upperlimit2000_Lambda160[5] = {0.034,0.046,0.0822221,0.106,0.138};
	Float_t upperlimit2000_Lambda180[5] = {0.03,0.038,0.0743565,0.094,0.122};
	Float_t upperlimit4000_Lambda100[5] = {0.062,0.078,0.135392,0.178,0.242};
	Float_t upperlimit4000_Lambda120[5] = {0.062,0.082,0.144527,0.182,0.25};
	Float_t upperlimit4000_Lambda140[5] = {0.054,0.07,0.123553,0.158,0.21};
	Float_t upperlimit4000_Lambda160[5] = {0.046,0.062,0.112221,0.142,0.186};
	Float_t upperlimit4000_Lambda180[5] = {0.038,0.054,0.0976872,0.126,0.158};
	Float_t upperlimit6000_Lambda100[5] = {0.086,0.11,0.18215,0.234,0.334};
	Float_t upperlimit6000_Lambda120[5] = {0.078,0.098,0.167433,0.21,0.306};
	Float_t upperlimit6000_Lambda140[5] = {0.066,0.094,0.157436,0.202,0.274};
	Float_t upperlimit6000_Lambda160[5] = {0.054,0.07,0.129418,0.17,0.226};



   Double_t obs_lim[7]= {0.};
   Double_t exp_lim[7]= {0.};

   

   obs_lim[0] = 0.0307005;
   obs_lim[1] = 0.0460321;
   obs_lim[2] = 0.0652659;
   obs_lim[3] = 0.0952963;
   obs_lim[4] = 0.13442;
   obs_lim[5] = 0.208523;
   obs_lim[6] = 0.276169;
	
   exp_lim[0] = upperlimit1_Lambda100[2];
   exp_lim[1] = upperlimit250_Lambda100[2];
   exp_lim[2] = upperlimit500_Lambda100[2];
   exp_lim[3] = upperlimit1000_Lambda100[2];
   exp_lim[4] = upperlimit2000_Lambda100[2];
   exp_lim[5] = upperlimit4000_Lambda100[2];
   exp_lim[6] = upperlimit6000_Lambda100[2];
 

   Double_t y_pdf_1sig[14] = {  upperlimit1_Lambda100[1],
                                upperlimit250_Lambda100[1],
                                upperlimit500_Lambda100[1],
                                upperlimit1000_Lambda100[1],
                                upperlimit2000_Lambda100[1],
                                upperlimit4000_Lambda100[1],
                                upperlimit6000_Lambda100[1],
                                upperlimit6000_Lambda100[3],
                                upperlimit4000_Lambda100[3],
                                upperlimit2000_Lambda100[3],
                                upperlimit1000_Lambda100[3],
                                upperlimit500_Lambda100[3],
                                upperlimit250_Lambda100[3],
                                upperlimit1_Lambda100[3],
                                };


   Double_t y_pdf_2sig[14] = {  upperlimit1_Lambda100[0],
                                upperlimit250_Lambda100[0],
                                upperlimit500_Lambda100[0],
                                upperlimit1000_Lambda100[0],
                                upperlimit2000_Lambda100[0],
                                upperlimit4000_Lambda100[0],
                                upperlimit6000_Lambda100[0],
                                upperlimit6000_Lambda100[4],
                                upperlimit4000_Lambda100[4],
                                upperlimit2000_Lambda100[4],
                                upperlimit1000_Lambda100[4],
                                upperlimit500_Lambda100[4],
                                upperlimit250_Lambda100[4],
                                upperlimit1_Lambda100[4]
                                };
   
   
   // arrays of theoretical cross section
   Double_t mTh[7]  = { 1., 250.0, 500.0, 1000.0, 2000., 4000., 6000. };
   Double_t xsTh[7] = { 0.2357,0.2357,0.2357,0.2357,0.2357,0.2357,0.2357};
   
   // axes labels for the final plot
   
   Double_t x_pdf[14]      = {1., 250., 500., 1000.,2000., 4000., 6000. ,6000., 4000., 2000.,  1000., 500., 250., 1.};
   TGraph* Onesig_graph;
   Onesig_graph = new TGraph(14., x_pdf, y_pdf_1sig );
   
   TGraph* Twosig_graph;
   Twosig_graph = new TGraph(14., x_pdf, y_pdf_2sig );
   
   
   
   Onesig_graph->SetFillColor(kGreen);
   Twosig_graph->SetFillColor(kYellow);

    

 TGraph* exp_lim_graph;
 exp_lim_graph  = new TGraph(7, mTh, exp_lim);
 exp_lim_graph->SetMarkerStyle(19);
 exp_lim_graph->SetMarkerSize(1.5);
 exp_lim_graph->SetMarkerColor(kRed);
 exp_lim_graph->SetLineColor(kBlack);
 exp_lim_graph->SetLineWidth(2.5);
 exp_lim_graph->SetLineStyle(kDashed);




 TGraph* ul_lim_graph;
 ul_lim_graph  = new TGraph(7, mTh, obs_lim);
 ul_lim_graph->SetMarkerStyle(22);
 ul_lim_graph->SetMarkerSize(1.5);
 ul_lim_graph->SetMarkerColor(kBlue);
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
  

   
   exp_lim_graph->Draw("LA");
   Twosig_graph->Draw("Fsame");
   Onesig_graph->Draw("Fsame");
   ul_lim_graph->Draw("Lsame");
   exp_lim_graph->Draw("Lsame");



   exp_lim_graph->GetXaxis()->SetTitle("c#tau_{#tilde{#chi^{0}_{1}}} (mm)");
   exp_lim_graph->GetYaxis()->SetTitle("#sigma ( #tilde{#chi^{0}_{1}} #rightarrow #tilde{G}#gamma) #times BR (pb)");
   exp_lim_graph->GetYaxis()->SetRangeUser(0.01, 10.);
   exp_lim_graph->GetXaxis()->SetRangeUser(1., 6000.);

   ul_lim_graph->Draw("same");

 // integrated luminosity
   std::string s_lumi;
   s_lumi = "4.55";
   std::string lint = "#intL = "+s_lumi+" fb^{-1}";
   TLatex l1;
   l1.SetTextAlign(12);
   l1.SetTextSize(0.035);
   l1.SetTextFont(22);
   l1.SetNDC();
   l1.DrawLatex(0.165, 0.967, "CMS Preliminary");
   l1.DrawLatex(0.75, 0.96, lint.c_str());

   




   TGraph *xsTh_vs_m = new TGraph(7, mTh, xsTh);
   xsTh_vs_m->SetLineWidth(2);
   xsTh_vs_m->SetLineColor(kRed);
  
   xsTh_vs_m->SetMarkerSize(1.);
   xsTh_vs_m->SetMarkerStyle(22);
   xsTh_vs_m->SetMarkerColor(kRed);
   xsTh_vs_m->Draw("Csame");
   
   
   TLegend* leg = new TLegend(0.5,0.6,0.95,0.9);
   leg->SetFillStyle(0); leg->SetBorderSize(0);
   leg->SetFillColor(0);

   leg->SetHeader("#tilde{#chi}^{0}_{1} #Lambda =  100 TeV");
   leg ->SetTextFont(22);
   //leg->AddEntry(pl,"CDF exclusion (2.6 fb^{-1})","f");
   leg->AddEntry(xsTh_vs_m,"Theoretical LO cross-section","l");
   leg->AddEntry(ul_lim_graph, "Observed  95% CL upper limit", "L");
   leg->AddEntry(exp_lim_graph, "Expected 95% CL upper limit", "L");
   leg->AddEntry(Onesig_graph, "#pm 1 #sigma Expected", "F");
   leg->AddEntry(Twosig_graph, "#pm 2 #sigma Expected", "F");
   leg->Draw("same");
   
   //c0->SetGridx();
   //c0->SetGridy();

   c0->SaveAs(("exclusion_limit_L"+LAMBDA+".eps").c_str());
}




void plot_limit_2d( ){

   

  TCanvas* c0 = new TCanvas("exclusion limit 2d ", "exclusion limit 2d", 1);
   c0->cd();
   //c0->SetLogy();
   



// region excluded by Tevatron limits
  Double_t x_shaded[14] = {100,110,120,130,140,146,149,149,146, 140, 130,120,110, 100};// CHANGED FOR LQ2
  Double_t y_shaded[14] = {0., 0., 0., 0., 0., 0., 0., 1., 2., 2., 2., 2., 2., 2.};// CHANGED FOR LQ2
   TGraph *pl = new TGraph(14,x_shaded,y_shaded);
   //  pl->SetFillStyle(3001);
   //pl->SetLineColor(0);
   pl->SetFillColor(kYellow +8);
   //pl->SetLineColor(kGray);   // CHANGED FOR LQ2
   pl->Draw("Af");

// region excluded by Tevatron limits
   Double_t x_shaded_cms[16] = {100,110,120,130,140,160.7,161,163,163,161,160.7,  140, 130,120,110, 100};// CHANGED FOR LQ2
   Double_t y_shaded_cms[16] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.83, 6.6, 6.6,6.6, 6.6, 6.6, 6.6, 6.6, 6.6, 6.6};// CHANGED FOR LQ2
   TGraph *pl_cms = new TGraph(16,x_shaded_cms,y_shaded_cms);
   //  pl->SetFillStyle(3001);
   //pl->SetLineColor(0);
   pl_cms->SetFillColor(kYellow +6);
   //pl->SetLineColor(kGray);   // CHANGED FOR LQ2
   pl_cms->Draw("AfL");
   pl->Draw("samefL");

   pl_cms->GetXaxis()->SetTitle("M_{#tilde{#chi^{0}_{1}}} [GeV/c^{2}]");
   pl_cms->GetYaxis()->SetTitle("Lifetime #tilde{#chi^{0}_{1}} [ns]");
   pl_cms->GetYaxis()->SetRangeUser(0., 10.);
   pl_cms->GetXaxis()->SetRangeUser(100., 180.);
// integrated luminosity
   std::string s_lumi;
   s_lumi = "4.55";
   std::string lint = "#intLdt= "+s_lumi+" pb^{-1}";
   TLatex l1;
   l1.SetTextAlign(12);
   l1.SetTextSize(0.04);
   l1.SetTextFont(22);
   l1.SetNDC();
   l1.DrawLatex(0.2,0.9,"CMS");
   l1.DrawLatex(0.2,0.75,lint.c_str());
   


   TLegend* leg = new TLegend(0.5,0.7,0.9,0.9);
   leg->SetFillColor(0);
   leg->SetBorderSize(0);
   leg ->SetTextFont(22);
   leg->AddEntry(pl,"CDF exclusion (2.6 fb^{-1})","f");
   leg->AddEntry(pl_cms,"Observed exclusion region ("+TString(lint)+")","f");
  
   // leg->Draw("same");
c0->SetGridx();
   c0->SetGridy();


   c0->SaveAs(("exclusion_limit_2d_"+s_lumi+"_with_expected.eps").c_str());
   c0->SaveAs(("exclusion_limit_2d_"+s_lumi+".C").c_str());
}



#include "tdrstyle.C"

void rootlogon() {
   gROOT->SetStyle("Plain");
   gStyle->SetPalette(1);
   gStyle->SetOptStat(1111111);  // Show overflow, underflow + SumOfWeights 
   gStyle->SetOptFit(111110);
   gStyle->SetOptFile(1);
   gStyle->SetMarkerStyle(20);
   gStyle->SetMarkerSize(0.8);
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



