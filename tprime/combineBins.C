#include <iostream>
#include <fstream>
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"

#include "TCanvas.h"
#include "TMath.h"
#include "TF1.h"
#include "TF2.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"

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
		  double femax_sig = 10.0 ); // separate bin precision for signal



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



void combineBins(int mass, double scale_factor = 1.0){
  //
  // Wrapper method for backwards compatibility
  //
  //   mass = mass of tprime quark

  // input file name
  //std::string fname = "data/ejets_3560/ht35/tprime_ejets_2D_3560ipb_08nov2011_ht35.root";
  std::string fname = "data/toys/templates_toy_smooth_2d_22nov2011v2.root";
  // output file name
  //char _template[256] = "data/ejets_3560/ht35/tprime_%i_ejets_3560ipb_ht35_2d_merged_14nov2011_001.root";
  char _template[256] = "output.root";
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

  // tweak the minimizer
  const char * minName = "Minuit2";
  const char * algoName = "Migrad";
  /*
  ROOT::Math::MinimizerOptions::SetDefaultPrecision(0.000001);
  ROOT::Math::MinimizerOptions::SetDefaultPrintLevel(1);
  */
  ROOT::Math::MinimizerOptions::SetDefaultMinimizer(minName, algoName);
  ROOT::Math::MinimizerOptions::SetDefaultTolerance(0.00001);
  std::cout << "MINIMIZER!!!! " 
	    << ROOT::Math::MinimizerOptions::DefaultMinimizerType()
	    << std::endl;
  std::cout << "MINIMIZER!!!! " 
	    << ROOT::Math::MinimizerOptions::DefaultMinimizerAlgo()
	    << std::endl;
  std::cout << "MINIMIZER!!!! " 
	    << ROOT::Math::MinimizerOptions::DefaultPrecision()
	    << std::endl;
  std::cout << "MINIMIZER!!!! " 
	    << ROOT::Math::MinimizerOptions::DefaultTolerance()
	    << std::endl;
  //std::exit(-1);


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
  //vBname.push_back("tprime%i__jer__plus");
  //vBname.push_back("tprime%i__jer__minus");
  //vBname.push_back("top__jer__plus");
  //vBname.push_back("top__jer__minus");
  //vBname.push_back("ewk__jer__plus");
  //vBname.push_back("ewk__jer__minus");
  vBname.push_back("tprime%i__match__plus");
  vBname.push_back("tprime%i__match__minus");
  vBname.push_back("top__match__plus");
  vBname.push_back("top__match__minus");
  vBname.push_back("ewk__match__plus");
  vBname.push_back("ewk__match__minus");
  vBname.push_back("tprime%i__testsyst__plus");
  vBname.push_back("tprime%i__testsyst__minus");
  vBname.push_back("top__testsyst__plus");
  vBname.push_back("top__testsyst__minus");
  vBname.push_back("ewk__testsyst__plus");
  vBname.push_back("ewk__testsyst__minus");

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
  std::string _fit_fname;
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
  fs_mass->SetParLimits(0,0,1000);
  fs_mass->SetParLimits(1,mass/2.0,1000);
  fs_mass->SetParLimits(2,10,mass);
  fs_mass->SetParLimits(3,0,1000);
  fs_mass->SetParLimits(4,100,mass);
  fs_mass->SetParLimits(5,10,1000);
  if ( channel.find("ejets")!=std::string::npos && abs(mass-400)<1 ){
    fs_mass->SetParameters(15,mass,50,15,mass-100.0,100);
    fs_mass->SetParLimits(0,0,100);
    fs_mass->SetParLimits(1,mass-100.0,mass+100.0);
    fs_mass->SetParLimits(2,10,mass/3.0);
    fs_mass->SetParLimits(3,0,100);
    fs_mass->SetParLimits(4,200,mass);
    fs_mass->SetParLimits(5,50,1000);
  }
  if ( channel.find("mujets")!=std::string::npos && abs(mass-525)<1 ){
    fs_mass->SetParameters(15,mass,50,15,mass-100.0,100);
    fs_mass->SetParLimits(0,0,100);
    fs_mass->SetParLimits(1,mass-100.0,mass+100.0);
    fs_mass->SetParLimits(2,10,mass/3.0);
    fs_mass->SetParLimits(3,0,100);
    fs_mass->SetParLimits(4,200,mass);
    fs_mass->SetParLimits(5,50,1000);
  }
  if ( channel.find("mujets")!=std::string::npos && abs(mass-450)<1 ){
    fs_mass->SetParameters(15,mass,50,15,mass-100.0,100);
    fs_mass->SetParLimits(0,0,100);
    fs_mass->SetParLimits(1,mass-100.0,mass+100.0);
    fs_mass->SetParLimits(2,10,mass/3.0);
    fs_mass->SetParLimits(3,0,100);
    fs_mass->SetParLimits(4,200,mass);
    fs_mass->SetParLimits(5,50,1000);
  }
  if ( channel.find("mujets")!=std::string::npos && abs(mass-550)<1 ){
    fs_mass->SetParameters(15,mass,50,15,mass-100.0,100);
    fs_mass->SetParLimits(0,0,100);
    fs_mass->SetParLimits(1,mass-100.0,mass+100.0);
    fs_mass->SetParLimits(2,10,mass/3.0);
    fs_mass->SetParLimits(3,0,100);
    fs_mass->SetParLimits(4,200,mass);
    fs_mass->SetParLimits(5,50,1000);
  }
  hs_mass->Fit("fs_mass","RM");
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
    fs_ht = new TF1("fs_ht", gaussian_ht, 350, 2000, 6);
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
  if ( channel.find("ejets")!=std::string::npos && abs(mass-400)<1 ){
    fs_ht->SetParLimits(3,130,1000);
    fs_ht->SetParLimits(5,100,1000);
  }
  if ( channel.find("ejets")!=std::string::npos && abs(mass-425)<1 ){
    fs_ht->SetParLimits(3,130,1000);
    fs_ht->SetParLimits(5,100,1000);
  }
  if ( channel.find("mujets")!=std::string::npos && abs(mass-425)<1 ){
    fs_ht->SetParLimits(3,130,1000);
    fs_ht->SetParLimits(5,100,1000);
  }
  if ( channel.find("mujets")!=std::string::npos && abs(mass-475)<1 ){
    fs_ht->SetParLimits(3,130,1000);
    fs_ht->SetParLimits(5,100,1000);
  }
  //
  hs_ht->Fit("fs_ht","RM");
  if (make_plots){
    TCanvas c3;
    hs_ht->Draw();
    _fit_fname = channel+"_sig_ht.pdf";
    c3.SaveAs( _fit_fname.c_str() );
  }

  // DEBUG: stop here
  //if (channel.find("mujets")!=std::string::npos) std::exit(-1);

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
  if ( channel.find("mujets")!=std::string::npos && abs(mass-425)<1 ){
    //fb_mass->SetParLimits(5,100.0,200.0); // Gaussian
    //fb_mass->SetParLimits(2,0.0,0.01); // LogNormal
  }
  hb_mass->Fit("fb_mass","RM");
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
  if ( channel.find("mujets")!=std::string::npos && abs(mass-550)<1 ){
    //delete fb_ht;
    //fb_ht = new TF1("fb_ht", bgHtPdf,_ymin,1000,6);
  }
  fb_ht->SetParameters(100,545,60,100,700,50);
  fb_ht->SetParLimits(1,200,2000);
  fb_ht->SetParLimits(2,30,1000);
  fb_ht->SetParLimits(4,200,2000);
  fb_ht->SetParLimits(5,30,1000);
  hb_ht->Fit("fb_ht","RM");
  if (make_plots){
    TCanvas c5;
    hb_ht->Draw();
    _fit_fname = channel+"_bkg_ht.pdf";
    c5.SaveAs( _fit_fname.c_str() );
  }

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
      //compute signal/background
      double _x = _xmin + min(_xmax, max(0.0, ((double)(i)-0.5)))*_bin_x_width;
      double _y = _ymin + min(_ymax, max(0.0, ((double)(j)-0.5)))*_bin_y_width;
      xsb[ibin] = log(1.0+ fs_mass->Eval(_x)*fs_ht->Eval(_y) 
		      / fb_mass->Eval(_x)/fb_ht->Eval(_y));
      _logsb = log(fs_mass->Eval(_x)*fs_ht->Eval(_y) 
		      / fb_mass->Eval(_x)/fb_ht->Eval(_y));

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
  delete fs_mass;
  delete fb_mass;
  delete fs_ht;
  delete fb_ht;
  // FIXME: 
  //delete _f2s;

  if (make_plots){  
    // plot S/B hist
    TCanvas c1;
    c1.SetLogy();
    h1sb->SetMinimum(0.00001);
    h1sb->Draw();
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
