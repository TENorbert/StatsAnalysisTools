/*Macro rebins histograms from ifilename file and saves it to */
//Last edited 12-05-11 Michael Luk

{
gROOT->Reset();
gStyle->SetCanvasColor(kWhite);
gStyle->SetPadBorderMode(0);

 bool _save = false; //dir is plots/
 
 



 TString tpmass = "350";
 const int rebval = 5;
 double datascale = 1.;
 
 TString ifilename = "andrew_ivanov.root" ;
 TFile *data = new TFile(ifilename); 
 
 TFile *_reb = new TFile(ifilename+"_reb.root","RECREATE");
 TH1F *tp_mass = (TH1F*)data->Get("TPrime"+tpmass+"_fitMass_4j_0t");
 TH1F *wj_mass = (TH1F*)data->Get("Wjets_fitMass_4j_0t");
 TH1F *ttbar_mass = (TH1F*)data->Get("Top_fitMass_4j_0t");
 TH1F *qcd_mass = (TH1F*)data->Get("DataQCD_fitMass_4j_0t");
 TH1F *sitop_mass = (TH1F*)data->Get("SingleTopT_fitMass_4j_0t");
 TH1F *data_mass = (TH1F*)data->Get("Data_fitMass_4j_0t");



   TCanvas *myCanvas = new TCanvas("myCanvas","My Canvas", 0,0,600,400);
   myCanvas->Clear();
   myCanvas->SetFillColor(0);
   myCanvas->SetLeftMargin(0.07);
   myCanvas->SetRightMargin(0.04);
   myCanvas->SetTopMargin(0.09);
   myCanvas->SetBottomMargin(0.13);
   myCanvas->SetFrameBorderMode(0);
  



/*----------- Now Scale Plots------------*/
 qcd_mass->SetFillColor(kMagenta);
 wj_mass->SetFillColor(kRed);
 ttbar_mass->SetFillColor(kBlue);
 tp_mass->SetFillColor(kGreen);
 sitop_mass->SetFillColor(kYellow);

 qcd_mass->Rebin(rebval) ;
 wj_mass->Rebin(rebval);
 ttbar_mass->Rebin(rebval);
 tp_mass->Rebin(rebval);
 data_mass->Rebin(rebval);
 sitop_mass->Rebin(rebval);
 data_mass->Scale(datascale);

THStack *massstack = new THStack("massstack","Mass Reconstruction");
massstack->Add(qcd_mass);
massstack->Add(wj_mass);
massstack->Add(ttbar_mass);
massstack->Add(tp_mass);
 massstack->Add(sitop_mass);

massstack->SetTitle("");
data_mass->SetMarkerStyle(2);

 massstack->Draw();
 data_mass->Draw("AE && same");
 massstack->SetMaximum(25);

 if((massstack->GetMaximum()) > (data_mass->GetMaximum())){
   massstack->SetMaximum(massstack->GetMaximum()*1.05);
 }
 else{massstack->SetMaximum((data_mass->GetMaximum())*1.05);}
 
 data_mass->SetTitle("");
 data_mass->SetStats(0);

 massstack->GetXaxis()->SetTitle("Mass_{RECO} [GeV]");
 massstack->GetXaxis()->CenterTitle();
 //myCanvas->SetLogy();
 
 TLegend *leg = new TLegend(0.6,0.66,0.88,0.90);
 leg->AddEntry(qcd_mass,"QCD","f");
 leg->AddEntry(wj_mass,"W+Jets.","f");
 leg->AddEntry(ttbar_mass,"TTbar","f");
 leg->AddEntry(sitop_mass,"Single Top","f");
 leg->AddEntry(tp_mass,"T' (GeV)","f");
 
 leg->AddEntry(data_mass,"Data","p");
 leg->SetHeader("   Mass Reconstruction for Min. #Chi^{2}");
 leg->SetTextSize(0.035);
 leg->SetFillColor(0);
 leg->SetLineColor(0);
 leg->Draw("same");
  
 
TCanvas *c2 = new TCanvas("c2","My Canvas 4", 0,0,600,400);
 c2->Clear();
 c2->Divide(3,2);

  c2->cd(1);
  tp_mass->Draw();
  c2->cd(2);
  wj_mass->Draw();
  c2->cd(3);
  qcd_mass->Draw();
  c2->cd(4);
  ttbar_mass->Draw();
  c2->cd(6);
  data_mass->Draw();
  c2->cd(5);
  sitop_mass->Draw();
  
  if(_save){
    c1->SaveAs("plots/"+tpmass+"gev-stacked.eps");
    c2->SaveAs("plots/"+tpmass+"gev-all.eps");
 }
  //delete myCanvas;
  //delete c2;

  tp_mass->Write();
  wj_mass->Write();
  qcd_mass->Write();
  data_mass->Write();
  ttbar_mass->Write();
  sitop_mass->Write();
  massstack->Write();

  _reb->Write();
  _reb->Close();
}

