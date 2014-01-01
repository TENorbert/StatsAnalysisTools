//// run new Hybrid on CMS model 
// run frequentist on CMS model

// standard demo for hypothesis test 


#include "TFile.h"
#include "RooWorkspace.h"
#include "RooAbsPdf.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooStats/ModelConfig.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TCanvas.h"
#include "TLine.h"

#include "RooStats/HybridCalculator.h"
#include "RooStats/FrequentistCalculator.h"
#include "RooStats/ToyMCSampler.h"
#include "RooStats/HypoTestPlot.h"

#include "RooStats/NumEventsTestStat.h"
#include "RooStats/ProfileLikelihoodTestStat.h"
#include "RooStats/SimpleLikelihoodRatioTestStat.h"
#include "RooStats/RatioOfProfiledLikelihoodsTestStat.h"
#include "RooStats/MaxLikelihoodEstimateTestStat.h"

#include "RooStats/HypoTestInverter.h"
#include "RooStats/HypoTestInverterResult.h"
#include "RooStats/HypoTestInverterPlot.h"

using namespace RooFit;
using namespace RooStats;


bool plotHypoTestResult = false; 

bool readResult = false;

void HypoTestInvDemo(const char * fileName ="GausModel_b.root",
                     const char * wsName = "w",
                     const char * modelSBName = "model_sb",
                     const char * modelBName = "model_b",
                     const char * dataName = "data_obs",                  
                     int type = 0,  // calculator type 
                     int testStatType = 0, // test stat type
                     int npoints = 10,   
                     int ntoys=1000,
                     bool useCls = true )
{ 
   /*
    type = 0 Freq calculator 
    type = 1 Hybrid 

    testStatType = 0 LEP
                 = 1 Tevatron 
                 = 2 PL


   */

   if (fileName==0) { 
      std::cout << "give input filename " << std::endl;
      return; 
   }
   TFile * file = new TFile(fileName); 

   RooWorkspace * w = dynamic_cast<RooWorkspace*>( file->Get(wsName) );
   if (!w) {      
      return; 
   }
   w->Print();


   RooAbsData * data = w->data(dataName); 
   if (!data) { 
      Error("HypoTestDemo","Not existing data %s",dataName);
   }

   
   // get models from WS
  // get the modelConfig out of the file
  ModelConfig* bModel = (ModelConfig*) w->obj(modelBName);
  ModelConfig* sbModel = (ModelConfig*) w->obj(modelSBName);


   SimpleLikelihoodRatioTestStat slrts(*bModel->GetPdf(),*sbModel->GetPdf());
   slrts.SetNullParameters(*bModel->GetSnapshot());
   slrts.SetAltParameters(*sbModel->GetSnapshot());

   RatioOfProfiledLikelihoodsTestStat 
   ropl(*bModel->GetPdf(), *sbModel->GetPdf(), sbModel->GetSnapshot());
   ropl.SetSubtractMLE(false);
   
   ProfileLikelihoodTestStat profll(*sbModel->GetPdf());
   profll.SetOneSided(0);

   TestStatistic * testStat = &slrts;
   if (testStatType == 1) testStat = &ropl;
   if (testStatType == 2) testStat = &profll;
  
   
   HypoTestCalculatorGeneric *  hc = 0;
   if (type == 0) hc = new FrequentistCalculator(*data, *sbModel, *bModel);
   else new HybridCalculator(*data, *sbModel, *bModel);

   ToyMCSampler *toymcs = (ToyMCSampler*)hc->GetTestStatSampler();
   //toymcs->SetNEventsPerToy(1);
   toymcs->SetTestStatistic(testStat);


    if (type == 1) { 
      HybridCalculator *hhc = (HybridCalculator*) hc;
      hhc->SetToys(ntoys,ntoys); 
      // hhc->ForcePriorNuisanceAlt(*pdfNuis);
      // hhc->ForcePriorNuisanceNull(*pdfNuis);
   } 
   else 
      ((FrequentistCalculator*) hc)->SetToys(ntoys,ntoys); 

  // Get the result
   RooMsgService::instance().getStream(1).removeTopic(RooFit::NumIntegration);


   TStopwatch tw; tw.Start(); 
   const RooArgSet * poiSet = sbModel->GetParametersOfInterest();
   RooRealVar *poi = (RooRealVar*)poiSet->first();

   // fit the data first
   sbModel->GetPdf()->fitTo(*data);
   double poihat  = poi->getVal();
   //poi->setVal(30);
   //poi->setError(10);


   HypoTestInverter calc(*hc);
   // GENA: for two-sided interval
   //calc.SetConfidenceLevel(0.95);
   // GENA: for 95% upper limit
   calc.SetConfidenceLevel(0.90);

   calc.UseCLs(useCls);
   calc.SetVerbose(true);

   // can spped up using proof
   ProofConfig pc(*w, 2, "workers=2", kFALSE);
   //ProofConfig pc(*w, 30, "localhost", kFALSE);
   //ToyMCSampler * toymcs = dynamic_cast<ToyMCSampler *> (calc.GetHypoTestCalculator()->GetTestStatSampler() );
   // GENA: disable proof for now
   //toymcs->SetProofConfig(&pc);    // enable proof

   
   if (npoints > 0) {
     // GENA
     double poimin = TMath::Max(poihat -   4 * poi->getError(), 0.0);
     //poimin = poihat;
     double poimax = poihat +  4 * poi->getError();
     poimin = 0; 
     poimax = 20;
     //double poimin = poi->getMin();
     //double poimax = poi->getMax();
     std::cout << "Doing a fixed scan  in interval : " << poimin << " , " << poimax << std::endl;
     calc.SetFixedScan(npoints,poimin,poimax);
   }

   HypoTestInverterResult * r = calc.GetInterval();

   // write to a file the results
   TString resultFileName = (useCls) ? "CLs_" : "Cls+b_";
   resultFileName += fileName;

   // GENA
   //TFile * file = new TFile(resultFileName,"RECREATE");
   file = new TFile(resultFileName,"RECREATE");
   r->Write();
   file->Close();                                                                     

  double ulError = r->UpperLimitEstimatedError();
  double upperLimit = r->UpperLimit();
  std::cout << "The computed upper limit is: " << upperLimit << std::endl;
  std::cout << "an estimated error on this upper limit is: " << ulError << std::endl;

  // check using interpolation
  // double interpLimit = r->FindInterpolatedLimit(1.-r->ConfidenceLevel() );
  // cout << "The computer interpolated limits is " << interpLimit << endl;

  const int nEntries = r->ArraySize();

  std::vector<Double_t> xArray(nEntries);
  std::vector<Double_t> yArray(nEntries);
  std::vector<Double_t> yErrArray(nEntries);
  for (int i=0; i<nEntries; i++) {
    xArray[i] = r->GetXValue(i);
    yArray[i] = r->GetYValue(i);
    yErrArray[i] = r->GetYError(i);
    std::cout << xArray[i] << " , " << yArray[i] << " err = " << yErrArray[i] << std::endl;
  }
 

   // see expected result (bands)
   TGraph * g0 = new TGraph(nEntries);
   TGraphAsymmErrors * g1 = new TGraphAsymmErrors(nEntries);
   TGraphAsymmErrors * g2l = new TGraphAsymmErrors(nEntries);
   TGraphAsymmErrors * g2u = new TGraphAsymmErrors(nEntries);
   double p[7]; 
   double q[7];
   p[0] = ROOT::Math::normal_cdf(-2);
   p[1] = ROOT::Math::normal_cdf(-1.5);
   p[2] = ROOT::Math::normal_cdf(-1);
   p[3] = 0.5;
   p[4] = ROOT::Math::normal_cdf(1);
   p[5] = ROOT::Math::normal_cdf(1.5);
   p[6] = ROOT::Math::normal_cdf(2);
   for (int i=0; i<nEntries; i++) {
      SamplingDistribution * s = r->GetExpectedDistribution(i);
      // GENA
      //const std::vector<double> & values = s->GetSamplingDistribution();
      const std::vector<Double_t> & cValues = s->GetSamplingDistribution();
      std::vector<Double_t> values;
      for (std::vector<Double_t>::const_iterator val = cValues.begin();
	   val != cValues.end();
	   ++val) values.push_back(*val);
      TMath::Quantiles(values.size(), 7, &values[0],q,p,false);
      double p0 = q[3];
      double p2l =  q[1];
      double p2u =  q[5];
      g0->SetPoint(i, r->GetXValue(i), p0 ) ;
      g1->SetPoint(i, r->GetXValue(i),  p0);
      g2l->SetPoint(i, r->GetXValue(i), p2l);
      g2u->SetPoint(i, r->GetXValue(i), p2u);
      //g2->SetPoint(i, r->GetXValue(i), s->InverseCDF(0.50));
      g1->SetPointEYlow(i, q[3] - q[2]); // -1 sigma errorr   
      g1->SetPointEYhigh(i, q[4] - q[3]);//+1 sigma error

      g2l->SetPointEYlow(i, q[1]-q[0]);   // -2 -- -1 sigma error
      g2l->SetPointEYhigh(i, q[2]-q[1]);

      g2u->SetPointEYlow(i, q[5]-q[4]);
      g2u->SetPointEYhigh(i, q[6]-q[5]);


      if (plotHypoTestResult) { 
         HypoTestResult * hr = new HypoTestResult();
         hr->SetNullDistribution( r->GetBackgroundDistribution() );
         hr->SetAltDistribution( r->GetSignalAndBackgroundDistribution(i) );
         new TCanvas();
         HypoTestPlot * pl = new HypoTestPlot(*hr);
         pl->Draw();
      }
  }

   HypoTestInverterPlot *plot = new HypoTestInverterPlot("result","",r);
   TGraphErrors * g = plot->MakePlot();

   g->Draw("APL");   
   g2l->SetFillColor(kYellow);
   g2l->Draw("3");
   g2u->SetFillColor(kYellow);
   g2u->Draw("3");
   g1->SetFillColor(kGreen);
   g1->Draw("3");
   g0->SetLineColor(kBlue);
   g0->SetLineStyle(2);
   g0->SetLineWidth(2);
   g0->Draw("L");

   //g1->Draw("P");
   //g2->Draw("P");
   g->SetLineWidth(2);
   g->Draw("PL");   

   // GENA: two-sided interval
   //double alpha = 1.-r->ConfidenceLevel();
   // GENA: upper limit
   double alpha = (1.-r->ConfidenceLevel())/2.0;
   double x1 = g->GetXaxis()->GetXmin();
   double x2 = g->GetXaxis()->GetXmax();
   TLine * line = new TLine(x1, alpha, x2,alpha);
   line->SetLineColor(kRed);
   line->Draw();

   // see the expected limit and -1 +1 sigma bands
   // SamplingDistribution * limits = r->GetUpperLimitDistribution();

   // std::cout << " expected limit (median) " << limits->InverseCDF(0.50) << std::endl;
   // std::cout << " expected limit (-1 sig) " << limits->InverseCDF((ROOT::Math::normal_cdf(-1))) << std::endl;
   // std::cout << " expected limit (+1 sig) " << limits->InverseCDF((ROOT::Math::normal_cdf(+1))) << std::endl;
   
   tw.Print();

}

int main() {
   HypoTestInvDemo("tprimeModel.root");
}
