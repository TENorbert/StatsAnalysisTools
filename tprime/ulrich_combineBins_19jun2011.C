#include <iostream>
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"

void combineBins(int mass){ //mass = mass of tprime quark

  //define some parameters
  char fname[100]={"statres_305ipb_v0_2d.root"}; //input file name
  char oname[256];={}; //output file name
  char sname[100]; //name of signal histogram
  sprintf(oname,"tprime_%i_mujets_statres_305ipb_v0_2d_merged_29jun2011.root",mass);
  sprintf(sname,"TPrime%i_HtvsMfit",mass);
  char bname[10][100]={ //array of data and background histograms
    "Data_HtvsMfit", //data histogram must be first in this list
    "TTjets_HtvsMfit",
    "Wjets_HtvsMfit",
    "QCD_HtvsMfit",
    "WW_HtvsMfit",
    "WZ_HtvsMfit",
    "ZZ_HtvsMfit",
    "Zjets_HtvsMfit",
    "SingleToptW_HtvsMfit",
    "SingleTopT_HtvsMfit"};
  
  int nb=10; //number of histograms in list
  float femax=0.25; //max fractional error in each bin of background histogram

  TFile *f = TFile::Open(fname);
  if (f==NULL) {
    printf("Cannot open file '%s'\n",fname);
    return;
  }

  TH2F* hs; f->GetObject(sname,hs); 
  if (hs==NULL) {
    printf("Cannot find histogram '%s' in '%s'\n",sname,fname);
    return;
  }

  TH2F *hb = (TH2F*)hs->Clone(); 
  hb->SetName("hb");
  hb->Reset();
  TH2F *hX[10];
  for (int i=0;i<nb;i++){
    f->GetObject(bname[i],hX[i]); 
    if (hX[i]==NULL) {
      printf("Cannot find histogram '%s' in '%s'\n",bname[i],fname);
      return;
    }
    //sum all background histograms into hb; do not add the data histogram
    if (i>0) hb->Add(hX[i]); 
  }

  //figure out the binning
  int nx = hs->GetNbinsX();
  int ny = hs->GetNbinsY();
  int nbin=nx*ny;
  std::cout << "number of bins: x="<<nx<<", y="<<ny<<std::endl;
  
  //book some 1d histograms with the same number of bins for diagnostics
  TH1F *h1sb = new TH1F("h1sb","h1sb",nbin,0,nbin);
  TH1F *h1s = new TH1F("h1s","h1s",nbin,0,nbin);
  TH1F *h1b = new TH1F("h1b","h1b",nbin,0,nbin);

  float xs,xb,eb;
  //xsb holds the s/b values for each bin
  //xx are the histogram contents 
  //(0=signal, 1=total background, 2=data, 3...nb-1=individual backgrounds)
  float xsb[30000],xx[30000][5],xe[30000][5];
  int ibin; 
  for (int i=0;i<nx;i++){
    for (int j=0;j<ny;j++){
      ibin=hs->GetBin(i,j);
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
    }
  }

  //sort all histogram bins in decreasing s/b
  int nswap=1;
  float xtmp;
  while (nswap>0) {
    nswap=0;
    for (int i=0;i<nbin-1;i++) {
      if (xsb[i]<xsb[i+1]){
	xtmp=xsb[i];
	xsb[i]=xsb[i+1];
	xsb[i+1]=xtmp;

	for (int j=0;j<5;j++){
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
    h1sb1->SetBinContent(i,xsb[i]);
    if (xx[i][1]>0) h1fe1->SetBinContent(i,xe[i][1]/xx[i][1]);
    h1s1->SetBinContent(i,xx[i][0]);
    h1s1->SetBinError(i,xe[i][0]);
    h1b1->SetBinContent(i,xx[i][1]);
    h1b1->SetBinError(i,xe[i][1]);
  }

  //combine bins starting with the highest s/b until the fractional error in
  //the total backround in every bin is smaller than femax
  int ncomb=1;
  float xtmp;
  float fe=0;
  while (ncomb>0) {
    ncomb=0;
    for (int i=0;i<nbin;i++){
      if (xx[i][1]>0){
	fe=xe[i][1]/xx[i][1]; //fractional error in background
      }else{
	fe=1;
      }
      if (fe>femax){
	for (int k=0;k<5;k++){ //add the next bin
	  xx[i][k]=xx[i][k]+xx[i+1][k];
	  xe[i][k]=sqrt(xe[i][k]**2+xe[i+1][k]**2);
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


  //these are the output histograms
  TFile *f2 = TFile::Open(oname,"recreate");
  TH1F *h1feb2 = new TH1F("h1fe2","h1fe2",nbin,0,nbin);
  TH1F *h1s2 = new TH1F(sname,sname,nbin,0,nbin);
  TH1F *h1b2 = new TH1F("h1b2","h1b2",nbin,0,nbin);
  TH1F *h1X2[10];
  for (int i=0;i<nb;i++){
    h1X2[i] = new TH1F(bname[i],bname[i],nbin,0,nbin);
  }
  for (int i=0;i<nbin;i++){
    h1feb2->SetBinContent(i,xe[i][1]/xx[i][1]);
    h1s2->SetBinContent(i,xx[i][0]);
    h1s2->SetBinError(i,xe[i][0]);
    h1b2->SetBinContent(i,xx[i][1]);
    h1b2->SetBinError(i,xe[i][1]);
    for (int j=0;j<nb;j++){
      h1X2[j]->SetBinContent(i,xx[i][j+2]);
      h1X2[j]->SetBinError(i,xe[i][j+2]);
    }
  }
  
  h1s2->Write();
  for (int j=0;j<nb;j++){
    h1X2[j]->Write();
  }
  
  f2->Close();
  f->Close();
}
