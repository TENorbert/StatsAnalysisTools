#include <iostream>
#include <fstream>
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"


// forward declaration
void combineBins( std::string channel, int mass, double scale_factor,
		  double femax,
		  std::string inFileName,
		  std::string outFileName,
		  std::string prefix = "" );



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
	     std::string prefix ){

  char * fname = (char *)(inFileName.c_str());

  char * oname = (char *)(outFileName.c_str());

  char sname[100]; //name of signal histogram
  sprintf(sname,"tprime%i",mass);

  char bname[20][100]={ //array of data and background histograms
    "DATA", // data histogram must be first in this list
    "top",
    "ewk",
    "tprime%i__jes__plus",
    "tprime%i__jes__minus",
    "top__jes__plus",
    "top__jes__minus",
    "ewk__jes__plus",
    "ewk__jes__minus"
  };

  
  int nb=9; //number of histograms in list
  int n_skip=3; // starting with this index, do not consider for background normalization
  ////////float femax=0.20; //max fractional error in each bin of background histogram

  TFile *f = TFile::Open(fname);
  if (f==NULL) {
    printf("Cannot open file '%s'\n",fname);
    return;
  }

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
  TH2D *hX[20];
  double bEntries[20]; // array with number of entries
  for (int i=0;i<nb;i++){
    std::string sBName(bname[i]);
    // GENA: get names for signal JES histos
    if (sBName.find("TPrime")!=std::string::npos ||
	sBName.find("Tprime")!=std::string::npos ||
	sBName.find("tprime")!=std::string::npos){

      sprintf(bname[i],sBName.c_str(),mass);
      std::cout << bname[i] << std::endl;
    }

    //f->GetObject(bname[i],hX[i]); 
    hX[i] = (TH2D *)f->Get(bname[i]); 
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

  float xs,xb,eb;
  //xsb holds the s/b values for each bin
  //xx are the histogram contents 
  //(0=signal, 1=total background, 2=data, 3...nb-1=individual backgrounds)
  float xsb[30000],xx[30000][20],xe[30000][20];
  int ibin; 
  double _sum = 0.0;
  for (int i=0;i<nx;i++){
    for (int j=0;j<ny;j++){

      ibin=hs->GetBin(i,j);

      // GENA: Will fill each bin with its original index
      vMap[ibin].push_back(ibin);

      xs=hs->GetBinContent(ibin);
      xb=hb->GetBinContent(ibin);
      //compute signal/background
      if (xb>0) {
	xsb[ibin]=xs/xb;
      }else{
	if (xs>0){
	  xsb[ibin]=999;
	}else{
	  xsb[ibin]=0;
	}
      }
      xx[ibin][0]=xs;
      xe[ibin][0]=hs->GetBinError(ibin);
      xx[ibin][1]=xb;
      xe[ibin][1]=hb->GetBinError(ibin);
      for (int k=0;k<nb;k++){
	xx[ibin][k+2]=hX[k]->GetBinContent(ibin);
	xe[ibin][k+2]=hX[k]->GetBinError(ibin);
      }
      if (xb>0) h1sb->SetBinContent(ibin,xs/xb);
      h1s->SetBinContent(ibin,xx[ibin][0]);
      h1s->SetBinError(ibin,xe[ibin][0]);
      h1b->SetBinContent(ibin,xx[ibin][1]);
      h1b->SetBinError(ibin,xe[ibin][1]);
      
      _sum += xx[ibin][0];
    }
  }

  std::cout << "SUM: " << _sum << std::endl;
  
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
  while (ncomb>0) {
    ncomb=0;
    for (int i=0;i<nbin-1;i++){
      if (xx[i][1]>0){
	fe=xe[i][1]/xx[i][1]; //fractional error in background
      }else{
	fe=1;
      }
      if (fe>femax){

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
  TH1F *h1X2[20];
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
