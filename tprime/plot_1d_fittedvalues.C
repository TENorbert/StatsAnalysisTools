// /*Macro rebins histograms from ifilename file and saves it to */
//Last edited 12-05-11 Michael Luk

TString convertDouble(double dbl){
  stringstream ss;//create a stringstream                                                         
  ss << dbl;//add number to the stream                                                         
  return ss.str();//return a string with the contents of the stream          
}

void plots(TString filename, TString massorht="mass", TString tpmass="550",TString emu="#mu+jets") {

 bool _save = true; 
 /* 
550
 Floating Parameter  InitialValue    FinalValue +/-  Error     GblCorr.
  --------------------  ------------  --------------------------  --------
                  Lumi    1.0000e+00    9.9774e-01 +/-  4.41e-02  0.492916
       alpha_btag_syst    0.0000e+00    6.2215e-02 +/-  7.23e-01  0.895540
             alpha_jes    0.0000e+00   -4.6553e-01 +/-  2.16e-01  0.815945
    alpha_lepton_eff_e    0.0000e+00    3.8358e-02 +/-  8.47e-01  0.842114
   alpha_norm_ewk_syst    0.0000e+00   -6.1907e-01 +/-  2.30e-01  0.850154
   alpha_norm_top_syst    0.0000e+00    8.7062e-02 +/-  1.01e-01  0.634419
                  xsec    0.0000e+00    1.0014e-06 +/-  3.37e-02  0.101128

 */
 double alumi = 0.99774;                                                                                                                                                                               
 double abtag = 1.+0.062215*0.05;                                                                                                                                                                
 double ajes  = 1.-0.46553*0.05;                                                                                                                                                                      
 double alep  = 1.+0.038358*0.03;                                                                                                                                                                        
 double anomewk=1.-0.061907*0.50;                                                                                                                                                                         
 double anomtop=1.+0.087062*0.114;                                                                                                                                                                       
 double xsec  = 0.0000010014*20;                                                                                                                                                                             
 TString lumi = "4.680";
 TString massorht="";

 // double tpscaling = 0.3296*5.;
 TString ifilename = filename ;
 TString jes = ""; //JESup_ or JESdown_ or nothing for nom
 TFile *data = new TFile(ifilename); 

 TH1D *data_1d = (TH1D*)data->Get("DATA");
 TH1D *tp_1d = (TH1D*)data->Get("tprime"+tpmass);//+"_ht35:fitMass_"+jes+massorht);
 TH1D *ewk_1d = (TH1D*)data->Get("ewk");//"_ht35:fitMass_"+jes+massorht);
 // TH1D *sintop_1d = (TH1D*)data->Get("sintop");//"_ht35:fitMass_"+jes+massorht);
 //TH1D *zjets_1d = (TH1D*)data->Get("Zjets");//"_ht35:fitMass_"+jes+massorht);
 TH1D *top_1d = (TH1D*)data->Get("top");//"_ht35:fitMass_"+jes+massorht);
 // TH1D *qcd_1d = (TH1D*)data->Get("QCD");//"_ht35:fitMass_"+jes+massorht);

 tp_1d->Scale(alumi*abtag*ajes*alep*xsec);
 ewk_1d->Scale(alumi*abtag*ajes*alep*anomewk);
 top_1d->Scale(alumi*abtag*ajes*alep*anomtop);



 TCanvas *myCanvas = new TCanvas("myCanvas","My Canvas", 0,0,600,400);
 myCanvas->Clear();
 myCanvas->SetFillColor(0);
 myCanvas->SetLeftMargin(0.1);//07);
 myCanvas->SetRightMargin(0.04);
 myCanvas->SetTopMargin(0.09);
 myCanvas->SetBottomMargin(0.13);
 myCanvas->SetFrameBorderMode(0);
 gStyle->SetOptStat(0);



 // /*----------- Now Scale Plots------------*/
 

 tp_1d->SetLineStyle(2);
 // tp_1d->Scale(tpscaling);
 ewk_1d->SetFillColor(kBlue);
 //zjets_1d->SetFillColor(kYellow);
 //sintop_1d->SetFillColor(7);
 top_1d->SetFillColor(8);
 //ewk_1d->SetFillColor(kOrange);

 if(mht == "ht"){
  tp_1d->Rebin(2);
 top_1d->Rebin(2);
 ewk_1d->Rebin(2);
 data_1d->Rebin(2);
 }

 double  errors[100]={0.0};
 tp_1d->SetError(errors);
 top_1d->SetError(errors);
 //zjets_1d->SetError(errors);
 ewk_1d->SetError(errors);
 //sintop_1d->SetError(errors);

 THStack *stack = new THStack("stack","Stack "+jes+massorht+" Reconstruction");
 //stack->Add(qcd_1d);
 //stack->Add(zjets_1d);
 //stack->Add(sintop_1d);



 stack->Add(ewk_1d);
 stack->Add(top_1d);
 //stack->Add(tp_1d);
 
 stack->SetTitle("");
 data_1d->SetMarkerStyle(2);


 top_1d->SetMarkerStyle(0);
 //qcd_1d->SetMarkerStyle(0);
 ewk_1d->SetMarkerStyle(0);
 //zjets_1d->SetMarkerStyle(0);
 //sintop_1d->SetMarkerStyle(0);
 tp_1d->SetMarkerStyle(0);

 stack->Draw("H");
 stack->GetXaxis()->SetLabelFont(62);
 stack->GetYaxis()->SetLabelFont(62);
 if(mht == "mass"){
 stack->GetYaxis()->SetTitle("Events /20GeV");
 stack->GetXaxis()->SetTitle("Mass_{RECO} [GeV]");
 }

 else if( mht == "ht"){
 stack->GetYaxis()->SetTitle("Events /40GeV");
 stack->GetXaxis()->SetTitle("H_{T} [GeV] Fitted Objects");
 }
 stack->GetYaxis()->SetTitleOffset(1.);
 stack->GetYaxis()->SetTitleFont(62);
 stack->GetYaxis()->SetTitleSize(0.05);
 stack->GetXaxis()->SetTitleSize(0.05);

 stack->GetXaxis()->SetLabelSize(0.045);
 stack->GetYaxis()->SetLabelSize(0.045);

 tp_1d->Draw("same H");
 data_1d->Draw("same P");
 data_1d->SetMarkerSize(0.5);
 data_1d->SetMarkerStyle(20);
 if((stack->GetMaximum()) > (data_1d->GetMaximum())){
   stack->SetMaximum(stack->GetMaximum()*1.05);
 }
 else{stack->SetMaximum((data_1d->GetMaximum())*1.05);}
 
 data_1d->SetTitle("");
 data_1d->SetStats(0);

 // stack->GetXaxis()->SetTitle(massorht+"_{RECO}");
 //stack->GetXaxis()->CenterTitle();


 //integral part
    /*                                                                                                           
 double counter[4] = {0.};
 for(int i = 0 ; i < data_1d->GetNbinsX(); ++i){
   counter[0] += data_1d->GetBinContent(i);
   counter[1] += ewk_1d->GetBinContent(i)+top_1d->GetBinContent(i);
   counter[2] += top_1d->GetBinContent(i);
   counter[3] += ewk_1d->GetBinContent(i);
   }*/
 TLegend *leg = new TLegend(0.63,0.6,0.85,0.85);
 leg->AddEntry(data_1d,"Data "+lumi+" fb^{-1}","pl");//+convertDouble(counter[0])+"Ev"+" (bg tot = "+convertDouble(counter[1])+"Ev)","p");
 leg->AddEntry(ewk_1d,"Ewk ","f");//+convertDouble(counter[3])+"Ev","f");
 leg->AddEntry(top_1d,"t#bar{t} ","f");//+convertDouble(counter[2])+"Ev","f");
 leg->AddEntry(tp_1d,"x20 t'#bar{t}' "+tpmass+"GeV","l");
 leg->SetHeader(emu);

 TText *text = new TLatex(0.6,0.935,"CMS preliminary #int L = "+lumi+" fb^{-1}");
 text->SetTextSize(0.04);
 text->SetNDC();
 text->Draw();
 /* TText *text = new TLatex(0.8,0.85,"e+jets");
 text->SetTextSize(0.04);
 text->SetNDC();
 text->Draw();
 */
 
 // leg->SetHeader("   "+massorht+" Reconstruction");
 leg->SetTextSize(0.035);
 leg->SetFillColor(0);
 leg->SetLineColor(0);
 leg->Draw("same");
  
 if(_save){
   myCanvas->SaveAs(emu+"_"+massorht+"_"+lumi+"ifb_fitted_hist.pdf");
 }

}

