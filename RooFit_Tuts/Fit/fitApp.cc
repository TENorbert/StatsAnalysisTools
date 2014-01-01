// $Id: fitApp.cc,v 1.1 2011/10/17 15:18:39 sigamani Exp $
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

double normGauss(double mu,double sig);
TGraph* convoluteGraph(const TGraph* ingr, const double& smear);

void computeUpperLimit(TGraph* hist, double xCL);
TGraph* fit_2D( int GMSBindex, float lumin,  std::string year, std::string ctau);
TGraph* fit_2D_bis( );

void plot_limit(/*Double_t* exp_lim,*/  std::string ctau);
void plot_limit_2d( );

double  exclusion(float lumin=34.4, std::string year="2010",  std::string ctau="1000");


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


 // dir to data files
 TString dirpath("/cmsrm/pc24/sigamani/final-oct/");

 // load one file for each sample
  fgmsb.push_back( (TFile*)  TFile::Open(dirpath+"PhotonJet_GMSB_Lambda-100_CTau-2000_7TeV_pythia6_cff-Summer11-PU_S4_START42_V11-v1-AODSIM_pfakt5_2ndJet10.root") );
  efficiency.push_back(0.40);
  fvec.push_back( (TFile*) TFile::Open(dirpath+"PhotonJet_DATA_2011_pfakt5_CONTROL_QCD_2ndJet10.root")); 
  fvec_pj.push_back( (TFile*) TFile::Open(dirpath+"PhotonJet_DATA_2011_pfakt5_CONTROL_G_2ndJet10.root"));
  fvec_data.push_back( (TFile*) TFile::Open(dirpath+"PhotonJet_DATA_2011_pfakt5_2ndJet10.root"));
  fvec_ttbar = (TFile*) TFile::Open(dirpath+"PhotonJet_TT_TuneZ2_7TeV-pythia6-tauola-Summer11-PU_S3-AODSIM_pfakt5_2ndJet10.root");
  fvec_Wenu = (TFile*) TFile::Open(dirpath+"PhotonJet_WToENu_TuneZ2_7TeV-pythia6-Summer11-PU_S3-AODSIM_pfakt5_2ndJet10.root");
 
  // contains the limit graph
  TGraph *t=0;
  t = fit_2D( 0, 1523., "May2011", "250" ); // take values from command line
  computeUpperLimit(t,0.95);
  //t = fit_2D_bis( ); // take values from command line
  //plot_limit("2000");
  exit(0);

} 






double exclusion(float lumin, std::string year,  std::string ctau) {

/* 
  std::cout<<"after loading files"<<std::endl;
  for(size_t i=0; i<fgmsb.size(); i++){
    t = fit_2D( i , lumin, year, ctau);
     computeUpperLimit(t, 0.95);
     TGraph *tmp = convoluteGraph(t, 20.);
     computeUpperLimit(tmp, 0.95);
  }

  std::vector <Double_t> xsec_theo;
  Double_t exp_limit[4]= {0.};
  double UL(0.) ;
  for(size_t i=0; i<fgmsb.size(); i++){
    std::cout<<"name: "<<fgmsb[i]->GetName()<<"   UPPERLIMIT: "<<ul[3*i+1]<<"   UPPERLIMIT_smer: "<<ul[3*i+2]<<"   sigma: "<<ul[2*i]<<"    eff: "<<efficiency[i]<< std::endl; 
 
    exp_limit[i] = (ul[3*i+1]);
    UL = (ul[3*i+1]);
  }
*/  
 
  //plot_limit( exp_limit, lumin, ctau );

   //ul.clear();
  //return UL;
  return 0.;
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
  //c1->SaveAs("scan-upperlimit.ps");
  //c1->SaveAs("scan-upperlimit.gif");
  //c1->SaveAs("scan-upperlimit.eps");
  
  //delete c1;
  ul.push_back( newx[iCL] );
}


//--------------
TGraph* fit_2D( int GMSBindex, float lumin,  std::string year, std::string ctau){
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

   RooPlot *xframe00 = epfMet.frame(Bins(9));
   h_allDATA_PDF.plotOn(xframe00);

   TH1* methist = h_allDATA.createHistogram("methist",epfMet);
   TH1* met_pdf_hist = hRedensity(methist);

   TH1* timehist = h_allDATA.createHistogram("methist",timePhotReco);
   TH1* time_pdf_hist = hRedensity(timehist);
   /*
   for(int i = 1; i<=methist->GetNbinsX();i++){
     std::cout << "i: " << i << " DATA: " << methist->GetBinContent(i) <<
                   " PDF DATA: " << met_pdf_hist->GetBinContent(i) << 
       " wid: " << met_pdf_hist->GetBinWidth(i) <<  std::endl;
   }
   // methist->Scale(h2_data->Integral()/methist->Integral());
   //methist->SumW2();
   
   TCanvas* c_h = new TCanvas("c_h", "work canvas", 1024, 768);
   c_h->cd();
   c_h->SetLogy();
   xframe00->SetMinimum(0.0001);
   xframe00->SetMaximum(5*10e5);
   xframe00->Draw();

    met_pdf_hist->Draw("psame");
   
   c_h->SaveAs("met.png");
   // c_h->Close();
   return 0;
   */

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
 
   if(year == "May2011") h2_QCD->Scale(18.55/33.84);
   if(year == "2010") h2_QCD->Scale(0.37);

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
 
   if(year == "May2011") h2_PJ->Scale(10.08/20.50);
   if(year == "2010") h2_PJ->Scale(0.05);
   if(year == "2010BKG") h2_PJ->Scale(lumin);

   RooDataHist hPJ("hPJ","hist PJ",RooArgSet(epfMet,timePhotReco), Import(*h2_PJ,false) );
   TH1* methistpj = hPJ.createHistogram("methistpj",epfMet);
   dumpTH1(methistpj);

   RooHistPdf PJ_PDF("PJ_PDF","met X time PJ",RooArgSet(epfMet,timePhotReco), hPJ ); 

   printf("GMSB: %f   QCD: %f     PJ: %f\n",h2_gmsb->Integral(),h2_QCD->Integral(),h2_PJ->Integral());
   printf("DATA %f\n",h2_data->Integral());
   

   TH2F* h2_sum = (TH2F*) h2_PJ->Clone("h2_sum");
    h2_sum->Add(h2_QCD);
    h2_sum->Add(h2_ttbar);
    h2_sum->Add(h2_Wenu);


    RooDataHist h_allMC("h_allMC","sum of all bkg samples MC and control sample", 
                         RooArgSet(epfMet,timePhotReco), Import(*h2_sum,false) );
    

   // std::cout<<"all MC"<<std::endl;

 
   RooRealVar lumi("lumi","luminosity",0., 1000.,"pb^{-1}");
   RooRealVar sigEff("sigEff","signal effciency", 0., 1.);
   RooRealVar sigAcc("sigAcc","signal acceptance",0., 1.);
   RooRealVar xsection("xsection","cross section", 0.1, 0., 10., "pb");
  
 

   lumi.setVal( lumin );
   sigEff.setVal( efficiency[GMSBindex] );
   sigAcc.setVal( 1. );

   xsection.setVal(h2_gmsb->Integral()/(lumin*efficiency[GMSBindex]));
    
   lumi.setConstant();
   sigEff.setConstant();
   sigAcc.setConstant();

  // define xsection as a function of lumi, eff and #signal events
   RooFormulaVar nsig("nsig","# signal events","@0*@1*@2*@3",RooArgList(lumi,sigAcc,sigEff,xsection) );
   RooRealVar nttbar("nttbar","nttbar",1,10000);
   RooRealVar nWenu("nWenu","nWenu",1,10000);
   RooRealVar nbkg("nbkg","nbkg",1,100000);
   RooRealVar nQCD("nQCD","# of QCD",1,100000);
   RooRealVar nPJ("nPJ","# gamma+jet",1,100000);


   // fraction of bkg in gamma+jet and QCD. this should be coming from MC NOT from control sample!
   //RooRealVar nfrac("nfrac", "nfrac", 0., 1.);
   //nfrac.setVal(0.335);
   //nfrac.setConstant();
   
   // bkg PDF as sum of  all backgrounds
   //RooAddPdf BKG("bkg", "bkg", RooArgList(PJ_PDF,QCD_PDF),nfrac);
   RooAddPdf BKG("bkg", "bkg", RooArgList(PJ_PDF,QCD_PDF,TTbar_PDF, Wenu_PDF ), RooArgList(nPJ, nQCD, nttbar, nWenu) );

   // total expected background with each component properly normalized
   //nbkg.setVal((h2_PJ->Integral()+h2_QCD->Integral())*1.15);
   //std::cout<<h2_PJ->Integral()+h2_QCD->Integral()<<" nbkg"<<std::endl;

   nPJ.setVal(h2_PJ->Integral());
   nQCD.setVal(h2_QCD->Integral());

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

   RooAddPdf model("model","signal + bkg",
                   RooArgList(signal,PJ_PDF,QCD_PDF,TTbar_PDF, Wenu_PDF),
                   RooArgList(nsig,nPJ, nQCD, nttbar, nWenu));
  
   // pointer to dataset
   RooDataHist& data = h_allDATA;
  
   
   // LL fit
   RooNLLVar nll("nll","nll",model,data,RooFit::SumW2Error(kTRUE) ,RooFit::Extended(kTRUE));
   RooMinuit mnll(nll);
   mnll.migrad() ;
   mnll.hesse() ;
   //mnll.minos() ;
   
   RooFitResult* rfit = mnll.fit("r");

   rfit->Print("v");

   std::cout << "Fitted signal yield: " << nsig.getVal()
             //<< " +/- " << nsig.getErr() 
             << std::endl;


   ul.push_back(xsection.getError());

   // plot of roo data hist
   RooPlot *xframe = epfMet.frame(100);
   xframe->SetXTitle("Missing energy [GeV]");
   data.plotOn(xframe, Invisible() );
   
  
  
   model.plotOn(xframe, LineWidth(1.5), LineColor(kBlue), Name("model"));
   model.plotOn(xframe,Components(signal),LineColor(kGreen), Name("signal")) ;
   BKG.plotOn(xframe, LineWidth(1.5), LineColor(kRed), LineStyle(kDashed), Name("BKG"));
   

   RooPlot *yframe = timePhotReco.frame(100);
   yframe->SetXTitle("ECAL time [ns]");
   data.plotOn(yframe, Invisible() );
 
   model.plotOn(yframe, LineWidth(1.5));
   model.plotOn(yframe,Components(signal),LineColor(kGreen)) ;
   BKG.plotOn(yframe, LineWidth(1.5), LineColor(kRed), LineStyle(kDashed));
   

   TH1F* h_leg_data= new TH1F("h_leg_data", "", 100, 0, 100);
   h_leg_data->SetMarkerStyle(20);
   h_leg_data->SetMarkerSize(1.4);
   h_leg_data->SetMarkerColor(kBlack);

   TH1F* h_leg_model= new TH1F("h_leg_model", "", 100, 0, 100);
   h_leg_model->SetLineWidth(2);
   h_leg_model->SetLineColor(kBlue);
   
   TH1F* h_leg_signal= new TH1F("h_leg_signal", "", 100, 0, 100);
   h_leg_signal->SetLineWidth(2);
   h_leg_signal->SetLineColor(kGreen);

   TH1F* h_leg_bkg= new TH1F("h_leg_bkg", "", 100, 0, 100);
   h_leg_bkg->SetLineWidth(2);
   h_leg_bkg->SetLineStyle(2);
   h_leg_bkg->SetLineColor(kRed);

   TLegend* leg = new TLegend(0.5,0.65,0.85,0.93);
   leg->SetFillColor(0);
   leg->SetBorderSize(0);
   leg ->SetTextFont(22);
   leg->AddEntry(h_leg_data,"data","PE");
   leg->AddEntry(h_leg_model,"signal + bkg","l");
   leg->AddEntry(h_leg_bkg,"Background","l");
   leg->AddEntry(h_leg_signal,("GMSB c#tau ="+ctau+" mm").c_str(),"l");


   TCanvas c2("c2","log scale canvas",1024,768);
   c2.SetLogy();
   xframe->Draw();
   xframe->SetMinimum( 0.1) ;
   xframe->SetMaximum( xframe->GetMaximum()*2.) ;
   met_pdf_hist->SetLineWidth(1.5); 
   met_pdf_hist->Draw("PESAME");
   leg->Draw("same");
   c2.SaveAs("plot_fit_met_data_log.png");
   c2.SaveAs("plot_fit_met_data_log.eps");


    yframe->Draw();
    yframe->SetMinimum( 0.1) ;
    yframe->SetMaximum( yframe->GetMaximum()*2.) ;
    time_pdf_hist->SetLineWidth(1.5); 
    time_pdf_hist->Draw("PESAME");
    leg->Draw("same");
    c2.SaveAs("plot_fit_time_data_log.png");
    c2.SaveAs("plot_fit_time_data_log.eps");


    TCanvas c3("c3","lin scale canvas",1024,768);
  
   xframe->Draw();
   xframe->SetMinimum( 0.1) ;
   leg->Draw("same");
   c3.SaveAs("plot_fit_met_data.png");


    yframe->Draw();
    yframe->SetMinimum( 0.1) ;
    leg->Draw("same");
    c3.SaveAs("plot_fit_time_data.png");
   
   //   return 0;
  
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
 
 
   /* delete h2_gmsb;
   delete h2_ttbar;
   delete h2_PJ;
   delete h2_QCD;
   delete h2_data;*/
   //h_allDATA.Delete();
   
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




void plot_limit(/*Double_t* exp_lim,*/std::string ctau ){

   Double_t ul_MC[4]= {0.};
   Double_t exp_lim[4]= {0.};
if(ctau == "250"){
    ul_MC[0] = 0.123557;
    ul_MC[1] = 0.385325;
    ul_MC[2] = 0.235194;
    ul_MC[3] = 0.169909;
    // ul_MC[4] = ;
  }
 

 // arrays of LQ masses for theoretical cross section
 Double_t mTh[4] = { 100*1.39, 120.*1.39, 140*1.39, 160*1.39/*, 180*1.39*/};
 Double_t xsTh[4] = { 0.2357, 0.0860, 0.0368, 0.0181/*, 0.0092*/};
  
Double_t x_pdf[8] = {	100*1.39, 120.*1.39, 140*1.39, 160*1.39, /*180*1.39,180*1.39,*/ 160*1.39,140*1.39,120.*1.39,100*1.39	};
 Double_t y_pdf[8];

 for(int i = 0; i<4; i++){

   y_pdf[i] = 0.07*xsTh[i]+xsTh[i];
 }
 int j = 3;
for(int i = 4; i<8; i++){

   y_pdf[i] = -0.07*xsTh[j]+xsTh[j];
   j--;
 }
 
Double_t y_pdf_exp[8];

 for(int ii = 0; ii<4; ii++){

   y_pdf_exp[ii] = 0.2*ul_MC[ii]+ul_MC[ii];
   std::cout<<y_pdf_exp[ii]<<std::endl;


 }
 int jj = 3;
for(int ii = 4; ii<8; ii++){

   y_pdf_exp[ii] = -0.2*ul_MC[jj]+ul_MC[jj];
   jj--;
   std::cout<<y_pdf_exp[ii]<<std::endl;
 }




  

  



// axes labels for the final plot



 if(ctau == "250"){
 exp_lim[0] =0.032589;
 exp_lim[1] =0.14044;
 exp_lim[2] =0.095976;
 exp_lim[3] =0.06562;
 //exp_lim[4] =
 }

if(ctau == "500"){
 exp_lim[0] = 0.;
 exp_lim[1] = 0.;
 exp_lim[2] = 0.;
 exp_lim[3] = 0.;
 //exp_lim[4] =
 }


if(ctau == "1000"){
 exp_lim[0] = 0.;
 exp_lim[1] = 0.;
 exp_lim[2] = 0.;
 exp_lim[3] = 0.;
 //exp_lim[4] =
 }


if(ctau == "2000"){
 exp_lim[0] = 0.06695;
 exp_lim[1] = 0.10828;
 exp_lim[2] = 0.10388;
 exp_lim[3] = 0.07631;
 //exp_lim[4] =
 }

 //TGraph* ul_MC_graph;
 //ul_MC_graph  = new TGraph(5, mass_neutralin, ul_MC);
 TGraph* exp_lim_graph;
 exp_lim_graph  = new TGraph(4, mTh, exp_lim);
 exp_lim_graph->SetMarkerStyle(22);
 exp_lim_graph->SetMarkerSize(2.5);
 exp_lim_graph->SetMarkerColor(kBlue);
 exp_lim_graph->SetLineColor(kBlue);
 exp_lim_graph->SetLineWidth(2.5);
 exp_lim_graph->SetLineStyle(2);




 TGraph* ul_lim_graph;
 ul_lim_graph  = new TGraph(4, mTh, ul_MC);
 ul_lim_graph->SetMarkerStyle(29);
 ul_lim_graph->SetMarkerSize(1.5);
 ul_lim_graph->SetMarkerColor(kRed);
 ul_lim_graph->SetFillColor(kGreen);
 ul_lim_graph->SetLineColor(kRed);
 ul_lim_graph->SetLineWidth(2);
 ul_lim_graph->SetLineStyle(2);

 // ul_MC_graph->SetMarkerStyle(22);
 //ul_MC_graph->SetMarkerColor(kRed);
 //ul_MC_graph->SetLineColor(kRed);
 

  TMultiGraph* mg = new TMultiGraph;
  
  mg->Add(exp_lim_graph);
  // mg->Add(ul_lim_graph);

  // char c_lumin[4];
    // sprintf(c_lumin, "%f", lumin);
    //   std::string s_lumin(c_lumin);

  
   TCanvas* c0 = new TCanvas("exclusion limit", "exclusion limit", 1);
   c0->cd();
   c0->SetLogy();
   // c0->SetGridx();
   //c0->SetGridy();
  

// region excluded by Tevatron limits
   Double_t x_shaded[8] = {130,135,140,149,149, 140, 135, 130};// CHANGED FOR LQ2
   Double_t y_shaded[8] = {0.0007,0.0007,0.0007,0.0007,5, 5, 5, 5};// CHANGED FOR LQ2
   TGraph *pl = new TGraph(8,x_shaded,y_shaded);
   //  pl->SetFillStyle(3001);
   //pl->SetLineColor(0);
   pl->SetFillColor(kGray);
   //pl->SetLineColor(kGray);   // CHANGED FOR LQ2
   pl->Draw("f");

TGraph *grshade_exp = new TGraph(8,x_pdf,y_pdf_exp);
   grshade_exp->SetFillColor(kGreen);
   //grshade_exp->SetFillStyle(2);
   //  grshade_exp->Draw("fsame");
  

   mg->Draw("PLsame");
   mg->GetXaxis()->SetTitle("M_{#tilde{#chi^{0}_{1}}} [GeV/c^{2}]");
   mg->GetYaxis()->SetTitle("#sigma ( #tilde{#chi^{0}_{1}} #rightarrow #tilde{G}#gamma)X BF [pb]");
 
// integrated luminosity
   std::string s_lumi;
   s_lumi = "858.4";
   std::string lint = "#intLdt= "+s_lumi+" pb^{-1}";
   TLatex l1;
   l1.SetTextAlign(12);
   l1.SetTextSize(0.04);
   l1.SetTextFont(22);
   l1.SetNDC();
   l1.DrawLatex(0.2,0.5,"CMS");
   l1.DrawLatex(0.2,0.4,lint.c_str());
   


if(ctau == "250"){
    ul_MC[0] = 0. ;
    ul_MC[1] = 0. ;
    ul_MC[2] = 0. ;
    ul_MC[3] = 0. ;
    // ul_MC[4] = ;
  }
 




   TGraph *grshade = new TGraph(8,x_pdf,y_pdf);
   grshade->SetFillColor(7);
   //grshade->SetFillStyle(2);
   grshade->Draw("fsame");
  

  


   TGraph *xsTh_vs_m = new TGraph(4, mTh, xsTh);
   xsTh_vs_m->SetLineWidth(2);
   xsTh_vs_m->SetLineColor(kRed);
   xsTh_vs_m->SetFillColor(7);
   xsTh_vs_m->SetMarkerSize(1.);
   xsTh_vs_m->SetMarkerStyle(22);
   xsTh_vs_m->SetMarkerColor(kRed);
   xsTh_vs_m->Draw("Csame");
   xsTh_vs_m->SetTitle("M_{#tilde{#chi^{0}_{1}}} [GeV/c^{2}]");
   xsTh_vs_m->SetTitle("#sigma (pb) ");
   
   //mg->Add(xsTh_vs_m);
   
   
   mg->GetYaxis()->SetRangeUser(0.0007, 5.);
   mg->GetXaxis()->SetRangeUser(130., 255);
   









   TLegend* leg = new TLegend(0.2,0.15,0.6,0.35);
   leg->SetFillColor(0);
   leg->SetBorderSize(1);
   leg->SetHeader(("#tilde{#chi}^{0}_{1} c#tau = "+ctau+" mm").c_str());
   leg ->SetTextFont(22);
   leg->AddEntry(pl,"CDF exclusion (2.6 fb^{-1})","f");
   leg->AddEntry(xsTh_vs_m,"Theoretical LO cross-section","lf");
   leg->AddEntry(exp_lim_graph, "Observed  95% CL upper limit", "LP");
   //  leg->AddEntry(ul_lim_graph, "Expected 95% CL upper limit", "LPF");
  
   leg->Draw("same");

   c0->SaveAs(("exclusion_limit_"+s_lumi+"_"+ctau+"_with_expected.eps").c_str());
   c0->SaveAs(("exclusion_limit_"+s_lumi+"_"+ctau+".C").c_str());
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
   s_lumi = "858.4";
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
   leg->AddEntry(pl_cms,"Observed exclusion region (858.4 pb^{-1})","f");
  
   leg->Draw("same");
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
   //std::cout << "Hey dude" << endl;

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

    //gROOT->ProcessLine(".L ~/CMS/macro/cmsPrel.C");
    gROOT->ProcessLine(".L cmsPrel.C");
    cout << "use cmsPrel() for official plots" << endl;

    //gROOT->ProcessLine(".L fit2d.cc+");


    cout << "rootlogon.C loaded" << endl;
}


TGraph* fit_2D_bis( /*int GMSBindex = 1, float lumi = 858.4,  std::string year = "May2011", std::string ctau = "250"*/){


//******define variables*********//
RooRealVar epfMet("epfMet","epfMet", 8., 800. );
RooRealVar timePhotReco("timePhotReco","timePhotReco",-2.,15.);

 Double_t lumin = 858.4;
 Double_t xsec_theo = 0.237;


  //**********************GMSB***********************//
 TFile* fgmsb =  TFile::Open("PhotonJet_GMSB_Lambda-100_CTau-250_7TeV_pythia6_cff-Summer11-PU_S4_START42_V11-v1-AODSIM-UNCLEAN2_pfakt5_2ndJet10.root");
  efficiency.push_back( 0.46 );

  //****get TH2 with variable bin size from file
TH2F* h2_met_time_gmsb = (TH2F*) fgmsb->Get("met_varsize_x_timePhot_varsize")->Clone("h2_met_time_gmsb");
TH2F* h2_met_time_dens_gmsb = (TH2F*) fgmsb->Get("met_varsize_x_timePhot_varsize")->Clone("h2_met_time_dens_gmsb");

//******divide per bin size
 int nbinsx_gmsb = h2_met_time_gmsb->GetNbinsX();
 int nbinsy_gmsb = h2_met_time_gmsb->GetNbinsY();

for(int i=1; i <=nbinsx_gmsb; ++i) {
for(int j=1; j <=nbinsy_gmsb; ++j) {

   float bin_area_gmsb = h2_met_time_gmsb->ProjectionX()->GetBinWidth(i)* h2_met_time_gmsb->ProjectionY()->GetBinWidth(j);
   float dens_gmsb = h2_met_time_gmsb->GetBinContent(i,j)/bin_area_gmsb;
   /*  cout << "i: " << i << " j: " << j
          << " content:" << h2_met_time->GetBinContent(i,j)
          << " wid: " << bin_area 
          << " dens: " << dens 
          << endl;*/
   if(h2_met_time_gmsb->GetBinContent(i,j)==0.) h2_met_time_gmsb->SetBinContent(i,j, 0.000001);
     dens_gmsb = h2_met_time_gmsb->GetBinContent(i,j)/bin_area_gmsb;
     h2_met_time_dens_gmsb->SetBinContent(i,j, dens_gmsb );
     h2_met_time_dens_gmsb->SetBinError(i,j, h2_met_time_gmsb->GetBinError(i,j)/bin_area_gmsb);
     if(h2_met_time_dens_gmsb->GetBinContent(i,j)-h2_met_time_dens_gmsb->GetBinError(i,j) <= 0.) h2_met_time_dens_gmsb->SetBinError(i,j, h2_met_time_dens_gmsb->GetBinContent(i,j) );
   
   }
   }


//****scale to lumi
 h2_met_time_dens_gmsb->Sumw2();
 h2_met_time_dens_gmsb->Scale(lumin);

 //******create a roodatahist
 TH2F* h2_ptr_gmsb = h2_met_time_dens_gmsb;
 RooDataHist met_time_gmsb("met_time_gmsb", "hist of MET and timing for gmsb", RooArgList(epfMet, timePhotReco), Import( *h2_ptr_gmsb, false ) );

 h2_ptr_gmsb->Print("V");

 RooPlot *xframe00_gmsb = epfMet.frame();
 met_time_gmsb.plotOn(xframe00_gmsb);

 TCanvas* c_gmsb = new TCanvas("c_gmsb", "", 1024,768);
 c_gmsb->cd();
 xframe00_gmsb->Draw();
 h2_ptr_gmsb->SetLineColor(kRed);
 h2_met_time_dens_gmsb->ProjectionX()->Draw("histsame");

 c_gmsb->SetLogy();
 
 c_gmsb->SaveAs("met_gmsb.png");
 c_gmsb->Close();


 //*** create a roohistpdf
RooHistPdf gmsb_pdf("gmsb_pdf","met X time GMSB",RooArgSet(epfMet,timePhotReco), met_time_gmsb ); 
   std::cout<<"gmsb pdf"<<std::endl;



//**********************TTBAR***********************//
TFile* fttbar =  TFile::Open("PhotonJet_TT_TuneZ2_7TeV-pythia6-tauola-Summer11-PU_S3-AODSIM-UNCLEAN2_pfakt5_2ndJet10.root");


  //****get TH2 with variable bin size from file
TH2F* h2_met_time_ttbar = (TH2F*) fttbar->Get("met_varsize_x_timePhot_varsize")->Clone("h2_met_time_ttbar");
TH2F* h2_met_time_dens_ttbar = (TH2F*) fttbar->Get("met_varsize_x_timePhot_varsize")->Clone("h2_met_time_dens_ttbar");

//******divide per bin size
 int nbinsx_ttbar = h2_met_time_ttbar->GetNbinsX();
 int nbinsy_ttbar = h2_met_time_ttbar->GetNbinsY();

for(int i=1; i <=nbinsx_ttbar; ++i) {
for(int j=1; j <=nbinsy_ttbar; ++j) {

   float bin_area_ttbar = h2_met_time_ttbar->ProjectionX()->GetBinWidth(i)* h2_met_time_ttbar->ProjectionY()->GetBinWidth(j);
   float dens_ttbar = h2_met_time_ttbar->GetBinContent(i,j)/bin_area_ttbar;
   /*  cout << "i: " << i << " j: " << j
          << " content:" << h2_met_time->GetBinContent(i,j)
          << " wid: " << bin_area 
          << " dens: " << dens 
          << endl;*/
   if(h2_met_time_ttbar->GetBinContent(i,j)==0.) h2_met_time_ttbar->SetBinContent(i,j, 0.000001);
     dens_ttbar = h2_met_time_ttbar->GetBinContent(i,j)/bin_area_ttbar;
     h2_met_time_dens_ttbar->SetBinContent(i,j, dens_ttbar );
     h2_met_time_dens_ttbar->SetBinError(i,j, h2_met_time_ttbar->GetBinError(i,j)/bin_area_ttbar);
     if(h2_met_time_dens_ttbar->GetBinContent(i,j)-h2_met_time_dens_ttbar->GetBinError(i,j) <= 0.) h2_met_time_dens_ttbar->SetBinError(i,j,  0.0000005);
   
   }
   }


//****scale to lumi
 h2_met_time_dens_ttbar->Sumw2();
 h2_met_time_dens_ttbar->Scale(lumin);

 //******create a roodatahist
 TH2F* h2_ptr_ttbar = h2_met_time_dens_ttbar;
 RooDataHist met_time_ttbar("met_time_ttbar", "hist of MET and timing for ttbar", RooArgList(epfMet, timePhotReco), Import( *h2_ptr_ttbar, false ) );

 h2_ptr_ttbar->Print("V");

 RooPlot *xframe00_ttbar = epfMet.frame();
 met_time_ttbar.plotOn(xframe00_ttbar);

 TCanvas* c_ttbar = new TCanvas("c_ttbar", "", 1024,768);
 c_ttbar->cd();
 xframe00_ttbar->Draw();
 h2_ptr_ttbar->SetLineColor(kRed);
 h2_met_time_dens_ttbar->ProjectionX()->Draw("histsame");

 c_ttbar->SetLogy();
 
 c_ttbar->SaveAs("met_ttbar.png");
 c_ttbar->Close();


 //*** create a roohistpdf
RooHistPdf ttbar_pdf("ttbar_pdf","met X time TTBAR",RooArgSet(epfMet,timePhotReco), met_time_ttbar ); 
   std::cout<<"ttbar pdf"<<std::endl;




//**********************WENU***********************//
TFile* fWenu =  TFile::Open("PhotonJet_WToENu_TuneZ2_7TeV-pythia6-Summer11-PU_S3_START42_V11-v2-AODSIM-UNCLEAN2_pfakt5_2ndJet10.root");


  //****get TH2 with variable bin size from file
TH2F* h2_met_time_Wenu = (TH2F*) fWenu->Get("met_varsize_x_timePhot_varsize")->Clone("h2_met_time_Wenu");
TH2F* h2_met_time_dens_Wenu = (TH2F*) fWenu->Get("met_varsize_x_timePhot_varsize")->Clone("h2_met_time_dens_Wenu");

//******divide per bin size
 int nbinsx_Wenu = h2_met_time_Wenu->GetNbinsX();
 int nbinsy_Wenu = h2_met_time_Wenu->GetNbinsY();

for(int i=1; i <=nbinsx_Wenu; ++i) {
for(int j=1; j <=nbinsy_Wenu; ++j) {

   float bin_area_Wenu = h2_met_time_Wenu->ProjectionX()->GetBinWidth(i)* h2_met_time_Wenu->ProjectionY()->GetBinWidth(j);
   float dens_Wenu = h2_met_time_Wenu->GetBinContent(i,j)/bin_area_Wenu;
   /*  cout << "i: " << i << " j: " << j
          << " content:" << h2_met_time->GetBinContent(i,j)
          << " wid: " << bin_area 
          << " dens: " << dens 
          << endl;*/
   if(h2_met_time_Wenu->GetBinContent(i,j)==0.) h2_met_time_Wenu->SetBinContent(i,j, 0.000001);
     dens_Wenu = h2_met_time_Wenu->GetBinContent(i,j)/bin_area_Wenu;
     h2_met_time_dens_Wenu->SetBinContent(i,j, dens_Wenu );
     h2_met_time_dens_Wenu->SetBinError(i,j, h2_met_time_Wenu->GetBinError(i,j)/bin_area_Wenu);
     if(h2_met_time_dens_Wenu->GetBinContent(i,j)-h2_met_time_dens_Wenu->GetBinError(i,j) <= 0.) h2_met_time_dens_Wenu->SetBinError(i,j, 0.0000005);
   
   }
   }


//****scale to lumi
 h2_met_time_dens_Wenu->Sumw2();
 h2_met_time_dens_Wenu->Scale(lumin);

 //******create a roodatahist
 TH2F* h2_ptr_Wenu = h2_met_time_dens_Wenu;
 RooDataHist met_time_Wenu("met_time_Wenu", "hist of MET and timing for Wenu", RooArgList(epfMet, timePhotReco), Import( *h2_ptr_Wenu, false ) );

 h2_ptr_Wenu->Print("V");

 RooPlot *xframe00_Wenu = epfMet.frame();
 met_time_Wenu.plotOn(xframe00_Wenu);

 TCanvas* c_Wenu = new TCanvas("c_Wenu", "", 1024,768);
 c_Wenu->cd();
 xframe00_Wenu->Draw();
 h2_ptr_Wenu->SetLineColor(kRed);
 h2_met_time_dens_Wenu->ProjectionX()->Draw("histsame");

 c_Wenu->SetLogy();
 
 c_Wenu->SaveAs("met_Wenu.png");
 c_Wenu->Close();


 //*** create a roohistpdf
RooHistPdf Wenu_pdf("Wenu_pdf","met X time WENU",RooArgSet(epfMet,timePhotReco), met_time_Wenu ); 
   std::cout<<"Wenu pdf"<<std::endl;



//**********************QCD***********************//
TFile* fQCD =  TFile::Open("PhotonJet_DATA_MayReReco_2011_UNCLEAN3_pfakt5_CONTROL_QCD_2ndJet10.root");


  //****get TH2 with variable bin size from file
TH2F* h2_met_time_QCD = (TH2F*) fQCD->Get("met_varsize_x_timePhot_varsize")->Clone("h2_met_time_QCD");
TH2F* h2_met_time_dens_QCD = (TH2F*) fQCD->Get("met_varsize_x_timePhot_varsize")->Clone("h2_met_time_dens_QCD");

//******divide per bin size
 int nbinsx_QCD = h2_met_time_QCD->GetNbinsX();
 int nbinsy_QCD = h2_met_time_QCD->GetNbinsY();

for(int i=1; i <=nbinsx_QCD; ++i) {
for(int j=1; j <=nbinsy_QCD; ++j) {

   float bin_area_QCD = h2_met_time_QCD->ProjectionX()->GetBinWidth(i)* h2_met_time_QCD->ProjectionY()->GetBinWidth(j);
   float dens_QCD = h2_met_time_QCD->GetBinContent(i,j)/bin_area_QCD;
   /*  cout << "i: " << i << " j: " << j
          << " content:" << h2_met_time->GetBinContent(i,j)
          << " wid: " << bin_area 
          << " dens: " << dens 
          << endl;*/
   if(h2_met_time_QCD->GetBinContent(i,j)==0.) h2_met_time_QCD->SetBinContent(i,j, 0.000001);
     dens_QCD = h2_met_time_QCD->GetBinContent(i,j)/bin_area_QCD;
     h2_met_time_dens_QCD->SetBinContent(i,j, dens_QCD );
     h2_met_time_dens_QCD->SetBinError(i,j, h2_met_time_QCD->GetBinError(i,j)/bin_area_QCD);
     if(h2_met_time_dens_QCD->GetBinContent(i,j)-h2_met_time_dens_QCD->GetBinError(i,j) <= 0.) h2_met_time_dens_QCD->SetBinError(i,j,  0.0000005);
   
   }
   }


//****scale to lumi
 h2_met_time_dens_QCD->Sumw2();
 h2_met_time_dens_QCD->Scale(18.54/33.84);

 //******create a roodatahist
 TH2F* h2_ptr_QCD = h2_met_time_dens_QCD;
 RooDataHist met_time_QCD("met_time_QCD", "hist of MET and timing for QCD", RooArgList(epfMet, timePhotReco), Import( *h2_ptr_QCD, false ) );

 h2_ptr_QCD->Print("V");

 RooPlot *xframe00_QCD = epfMet.frame();
 met_time_QCD.plotOn(xframe00_QCD);

 TCanvas* c_QCD = new TCanvas("c_QCD", "", 1024,768);
 c_QCD->cd();
 xframe00_QCD->Draw();
 h2_ptr_QCD->SetLineColor(kRed);
 h2_met_time_dens_QCD->ProjectionX()->Draw("histsame");

 c_QCD->SetLogy();
 
 c_QCD->SaveAs("met_QCD.png");
 c_QCD->Close();


 //*** create a roohistpdf
RooHistPdf QCD_pdf("QCD_pdf","met X time QCD",RooArgSet(epfMet,timePhotReco), met_time_QCD ); 
   std::cout<<"QCD pdf"<<std::endl;


//**********************PJ***********************//
TFile* fPJ =  TFile::Open("PhotonJet_DATA_MayReReco_2011_UNCLEAN3_pfakt5_CONTROL_G_2ndJet10.root");


  //****get TH2 with variable bin size from file
TH2F* h2_met_time_PJ = (TH2F*) fPJ->Get("met_varsize_x_timePhot_varsize")->Clone("h2_met_time_PJ");
TH2F* h2_met_time_dens_PJ = (TH2F*) fPJ->Get("met_varsize_x_timePhot_varsize")->Clone("h2_met_time_dens_PJ");

//******divide per bin size
 int nbinsx_PJ = h2_met_time_PJ->GetNbinsX();
 int nbinsy_PJ = h2_met_time_PJ->GetNbinsY();

for(int i=1; i <=nbinsx_PJ; ++i) {
for(int j=1; j <=nbinsy_PJ; ++j) {

   float bin_area_PJ = h2_met_time_PJ->ProjectionX()->GetBinWidth(i)* h2_met_time_PJ->ProjectionY()->GetBinWidth(j);
   float dens_PJ = h2_met_time_PJ->GetBinContent(i,j)/bin_area_PJ;
    cout << "i: " << i << " j: " << j
          << " content:" << h2_met_time_PJ->GetBinContent(i,j)
          << " wid: " << bin_area_PJ 
          << " dens: " << dens_PJ 
          << endl;
    if(h2_met_time_PJ->GetBinContent(i,j)==0.) h2_met_time_PJ->SetBinContent(i,j, 0.000001);
     dens_PJ = h2_met_time_PJ->GetBinContent(i,j)/bin_area_PJ;
     cout  << " dens: " << dens_PJ << endl;
     h2_met_time_dens_PJ->SetBinContent(i,j, dens_PJ );
     h2_met_time_dens_PJ->SetBinError(i,j, h2_met_time_PJ->GetBinError(i,j)/bin_area_PJ);
     if((h2_met_time_dens_PJ->GetBinContent(i,j)-h2_met_time_dens_PJ->GetBinError(i,j)) <= 0.) h2_met_time_dens_PJ->SetBinError(i,j, 0.0000005);
   
   }
   }


//****scale to lumi
 h2_met_time_dens_PJ->Sumw2();
 h2_met_time_dens_PJ->Scale(10.08/20.5);

 //******create a roodatahist
 TH2F* h2_ptr_PJ = h2_met_time_dens_PJ;
 RooDataHist met_time_PJ("met_time_PJ", "hist of MET and timing for PJ", RooArgList(epfMet, timePhotReco), Import( *h2_ptr_PJ, false ) );

 h2_ptr_PJ->Print("V");

 RooPlot *xframe00_PJ = epfMet.frame();
 met_time_PJ.plotOn(xframe00_PJ);

 TCanvas* c_PJ = new TCanvas("c_PJ", "", 1024,768);
 c_PJ->cd();
 xframe00_PJ->Draw();
 h2_ptr_PJ->SetLineColor(kRed);
 h2_met_time_dens_PJ->ProjectionX()->Draw("histsame");

 c_PJ->SetLogy();
 
 c_PJ->SaveAs("met_PJ.png");
 c_PJ->Close();


 //*** create a roohistpdf
RooHistPdf PJ_pdf("PJ_pdf","met X time PJ",RooArgSet(epfMet,timePhotReco), met_time_PJ ); 
   std::cout<<"PJ pdf"<<std::endl;



//**********************SUM***********************//

   TH2F* h2_met_time_SUM = (TH2F*) h2_met_time_PJ->Clone("h2_met_time_SUM");
   h2_met_time_SUM->Add(h2_met_time_QCD);
   h2_met_time_SUM->Add(h2_met_time_ttbar);
   h2_met_time_SUM->Add(h2_met_time_Wenu);

   TH2F* h2_met_time_dens_SUM = (TH2F*) h2_met_time_SUM->Clone("h2_met_time_dens_SUM");

//******divide per bin size
 int nbinsx_SUM = h2_met_time_SUM->GetNbinsX();
 int nbinsy_SUM = h2_met_time_SUM->GetNbinsY();


 //******create a roodatahist
 TH2F* h2_ptr_SUM = h2_met_time_dens_SUM;
 RooDataHist met_time_SUM("met_time_SUM", "hist of MET and timing for SUM", RooArgList(epfMet, timePhotReco), Import( *h2_ptr_SUM, false ) );

 h2_ptr_SUM->Print("V");

 RooPlot *xframe00_SUM = epfMet.frame();
 met_time_SUM.plotOn(xframe00_SUM);

 TCanvas* c_SUM = new TCanvas("c_SUM", "", 1024,768);
 c_SUM->cd();
 xframe00_SUM->Draw();
 h2_ptr_SUM->SetLineColor(kRed);
 h2_met_time_dens_SUM->ProjectionX()->Draw("histsame");

 c_SUM->SetLogy();
 
 c_SUM->SaveAs("met_SUM.png");
 c_SUM->Close();


 
//**********************DATA***********************//
TFile* fDATA =  TFile::Open("PhotonJet_DATA_MayReReco_2011_UNCLEAN3_pfakt5_2ndJet10.root");


  //****get TH2 with variable bin size from file
TH2F* h2_met_time_DATA = (TH2F*) fDATA->Get("met_varsize_x_timePhot_varsize")->Clone("h2_met_time_DATA");
TH2F* h2_met_time_dens_DATA = (TH2F*) fDATA->Get("met_varsize_x_timePhot_varsize")->Clone("h2_met_time_dens_DATA");

//******divide per bin size
 int nbinsx_DATA = h2_met_time_DATA->GetNbinsX();
 int nbinsy_DATA = h2_met_time_DATA->GetNbinsY();

for(int i=1; i <=nbinsx_DATA; ++i) {
for(int j=1; j <=nbinsy_DATA; ++j) {

   float bin_area_DATA = h2_met_time_DATA->ProjectionX()->GetBinWidth(i)* h2_met_time_DATA->ProjectionY()->GetBinWidth(j);
   float dens_DATA = h2_met_time_DATA->GetBinContent(i,j)/bin_area_DATA;
   /*  cout << "i: " << i << " j: " << j
          << " content:" << h2_met_time->GetBinContent(i,j)
          << " wid: " << bin_area 
          << " dens: " << dens 
          << endl;*/
    if(h2_met_time_DATA->GetBinContent(i,j)==0.) h2_met_time_DATA->SetBinContent(i,j, 0.000001);
     dens_DATA = h2_met_time_DATA->GetBinContent(i,j)/bin_area_DATA;
     h2_met_time_dens_DATA->SetBinContent(i,j, dens_DATA );
     h2_met_time_dens_DATA->SetBinError(i,j, h2_met_time_DATA->GetBinError(i,j)/bin_area_DATA);
     if(h2_met_time_dens_DATA->GetBinContent(i,j)-h2_met_time_dens_DATA->GetBinError(i,j) <= 0.) h2_met_time_dens_DATA->SetBinError(i,j, 0.0000005);
   
   }
   }


 h2_met_time_dens_DATA->Sumw2();
 

 //******create a roodatahist
 TH2F* h2_ptr_DATA = h2_met_time_dens_DATA;
 RooDataHist met_time_DATA("met_time_DATA", "hist of MET and timing for DATA", RooArgList(epfMet, timePhotReco), Import( *h2_ptr_DATA, false ) );

 h2_ptr_DATA->Print("V");

 RooPlot *xframe00_DATA = epfMet.frame();
 met_time_DATA.plotOn(xframe00_DATA);

 TCanvas* c_DATA = new TCanvas("c_DATA", "", 1024,768);
 c_DATA->cd();
 xframe00_DATA->SetMinimum( 0.0001) ;
 xframe00_DATA->Draw();
 h2_ptr_DATA->SetLineColor(kRed);
 h2_met_time_dens_DATA->ProjectionX()->Draw("histsame");

 c_DATA->SetLogy();
 
 c_DATA->SaveAs("met_DATA.png");
 c_DATA->Close();






TGraph* gr_nll = 0;
return gr_nll;

 
   RooRealVar lumi("lumi","luminosity",0., 1000.,"pb^{-1}");
   RooRealVar sigEff("sigEff","signal effciency", 0., 1.);
   RooRealVar sigAcc("sigAcc","signal acceptance",0., 1.);
   RooRealVar xsection("xsection","cross section", 0., 10., "pb");
  
 

   lumi.setVal( lumin );
   sigEff.setVal( efficiency[0] );
   sigAcc.setVal( 1. );
   xsection.setVal(xsec_theo);


   lumi.setConstant();
   sigEff.setConstant();
   sigAcc.setConstant();
  
  // define xsection as a function of lumi, eff and #signal events
   RooFormulaVar nsig("nsig","# signal events","@0*@1*@2*@3",RooArgList(lumi,sigAcc,sigEff,xsection) );
   RooRealVar nttbar("nttbar","nttbar",1,10000);
   RooRealVar nWenu("nWenu","nWenu",1,10000);
   RooRealVar nbkg("nbkg","nbkg",1,100000);


   // fraction of bkg in gamma+jet and QCD. this should be coming from MC NOT from control sample!
   RooRealVar nfrac("nfrac", "nfrac", 0., 1.);
   nfrac.setVal(0.356);
   nfrac.setConstant();
   
   // bkg PDF as sum of QCD and g+j 
   RooAddPdf BKG("bkg", "bkg", RooArgList(PJ_pdf,QCD_pdf),nfrac);
  
  
   // total expected background with each component properly normalized
   nbkg.setVal(h2_ptr_QCD->Integral()+h2_ptr_PJ->Integral());

   std::cout<<h2_ptr_PJ->Integral()+h2_ptr_QCD->Integral()<<" nbkg"<<std::endl;

   // ttbar  and Wenu normnalization. Now using theory xsec but should switch to CMS measurement
   nttbar.setVal(h2_ptr_ttbar->Integral());
   nttbar.setConstant();

   nWenu.setVal(h2_ptr_Wenu->Integral());
   //nWenu.setConstant();

   RooAddPdf model("model","signal + bkg",RooArgList(gmsb_pdf,BKG,ttbar_pdf, Wenu_pdf),RooArgList(nsig,nbkg, nttbar, nWenu));
  
   
   // pointer to dataset
   RooDataHist& data = met_time_DATA;

   /*   RooPlot *xframe0 = epfMet.frame(100);
   data.plotOn(xframe0);
  
   model.plotOn(xframe0, LineWidth(1.5));
   //  model.plotOn(xframe0,Components(gmsb_pdf),LineColor(kRed)) ;
   
   
   TCanvas c0("c2","c2",1024,768);
   c0.SetLogy();
   xframe0->Draw();
   xframe0->SetMinimum( 0.1) ;
   xframe0->SetMaximum( xframe0->GetMaximum()*2.  );
   c0.SaveAs("plot_met_data_before_fit.png");
   c0.SaveAs("plot_met_data_before_fit.eps");

   */
   // LL fit
   RooNLLVar nll("nll","nll",model,data,RooFit::SumW2Error(kTRUE) ,RooFit::Extended(kTRUE));
   RooMinuit mnll(nll);
   mnll.migrad() ;
   mnll.hesse() ;
   mnll.minos() ;
   

   std::cout<<nsig.getVal()<<"   -> Nsig"<<std::endl;

   ul.push_back(xsection.getError());

   /* RooPlot *xframe = epfMet.frame(100);
   xframe->SetXTitle("Missing energy [GeV]");
   xframe->SetYTitle("Entries/[GeV]");
   data.plotOn(xframe);
  
   model.plotOn(xframe, LineWidth(1.5), LineColor(kBlue), Name("model"));
   model.plotOn(xframe,Components(gmsb_pdf),LineColor(kGreen), Name("gmsb_pdf")) ;
   model.plotOn(xframe, Components(BKG), LineWidth(1.5), LineColor(kRed), LineStyle(kDashed), Name("BKG"));
   

   RooPlot *yframe = timePhotReco.frame(100);
   yframe->SetXTitle("ECAL time [ns]");
   yframe->SetYTitle("Entries/[ns]");
   data.plotOn(yframe );
 
   model.plotOn(yframe, LineWidth(1.5), Name("model"));
   model.plotOn(yframe,Components(gmsb_pdf),LineColor(kGreen), Name("gmsb_pdf")) ;
   model.plotOn(yframe, Components(BKG), LineWidth(1.5), LineColor(kRed), LineStyle(kDashed), Name("BKG"));
   

   TH1F* h_leg_data= new TH1F("h_leg_data", "", 100, 0, 100);
   h_leg_data->SetMarkerStyle(20);
   h_leg_data->SetMarkerColor(kBlack);

   TH1F* h_leg_model= new TH1F("h_leg_model", "", 100, 0, 100);
   h_leg_model->SetLineWidth(2);
   h_leg_model->SetLineColor(kBlue);
   
   TH1F* h_leg_signal= new TH1F("h_leg_signal", "", 100, 0, 100);
   h_leg_signal->SetLineWidth(2);
   h_leg_signal->SetLineColor(kGreen);

   TH1F* h_leg_bkg= new TH1F("h_leg_bkg", "", 100, 0, 100);
   h_leg_bkg->SetLineWidth(2);
   h_leg_bkg->SetLineStyle(2);
   h_leg_bkg->SetLineColor(kRed);

   TLegend* leg = new TLegend(0.5,0.65,0.85,0.93);
   leg->SetFillColor(0);
   leg->SetBorderSize(0);
   leg ->SetTextFont(22);
   leg->AddEntry(h_leg_data,"data","P");
   leg->AddEntry(h_leg_model,"signal + bkg","l");
   leg->AddEntry(h_leg_bkg,"Background","l");
   leg->AddEntry(h_leg_signal,"GMSB c#tau = 250 mm","l");


   TCanvas c2("c2","c2",1024,768);
   c2.SetLogy();
   xframe->Draw();
   xframe->SetMinimum( 0.1) ;
   xframe->SetMaximum( xframe->GetMaximum()*2.) ;
   leg->Draw("same");
   c2.SaveAs("plot_fit_met_data_log.png");
   c2.SaveAs("plot_fit_met_data_log.eps");


    yframe->Draw();
    yframe->SetMinimum( 0.1) ;
    yframe->SetMaximum( yframe->GetMaximum()*2.) ;
    leg->Draw("same");
    c2.SaveAs("plot_fit_time_data_log.png");
    c2.SaveAs("plot_fit_time_data_log.eps");


    TCanvas c3("c3","c3",1024,768);
  
   xframe->Draw();
   xframe->SetMinimum( 0.1) ;
   leg->Draw("same");
   c3.SaveAs("plot_fit_met_data.png");


    yframe->Draw();
    yframe->SetMinimum( 0.1) ;
    leg->Draw("same");
    c3.SaveAs("plot_fit_time_data.png");
   */
   
  
   double minNLL = nll.getVal();

   // scan LL vs. x section
   // TGraph* gr_nll = 0;
   //<return gr_nll;
    TH1* h_nll = 0;/*
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
  
   delete h2_gmsb;
   delete h2_ttbar;
   delete h2_PJ;
   delete h2_QCD;
   delete h2_data;
   h_allDATA.Delete();
		   */
   return gr_nll;
}



