//// run new Hybrid on CMS model 
// run frequentist on CMS model
//
// Author: Lorenzo Moneta
//
// June 2011: first version, Lorenzo Moneta
// June 2011: interface changes, Gena Kukartsev
//
//
//

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
bool useProof = false;
int nworkers = 2;
bool readResult = false;


// internal routine to run the inverter
HypoTestInverterResult * RunInverter(RooWorkspace * w, const char * modelSBName, const char * modelBName, const char * dataName,
                                     int type,  int testStatType, int npoints, double poimin, double poimax, int ntoys, bool useCls );






std::vector<Double_t> HypoTestInvDemo(RooWorkspace * pWs,
                     const char * modelSBName = "model_sb",
                     const char * modelBName = "model_b",
                     const char * dataName = "data_obs",                  
                     int calculatorType = 0,  // calculator type 
                     int testStatType = 3, // test stat type
                     int npoints = 10,   
                     double poimin = 1,  // use default is poimin >= poimax
                     double poimax = 0,
                     int ntoys=1000,
                     bool useCls = true,
		     std::string suffix = "test") { 
  //
  // Return a vector of numbers (terrible design, I know) ordered as
  //  - observed limit
  //  - observed limit error
  //  - interpolated limit (same as observed?)
  //  - expected limit -2 sigma
  //  - expected limit -1 sigma
  //  - expected limit median
  //  - expected limit +1 sigma
  //  - expected limit +2 sigma
  //

  std::vector<Double_t> result;

  // check that workspace is present
  if (!pWs){
    std::cout << "No workspace found, null pointer" << std::endl;
    return result;
  }

   HypoTestInverterResult * r = 0; 
   r = RunInverter(pWs, modelSBName, modelBName, dataName, calculatorType, testStatType, npoints, poimin, poimax,  ntoys, useCls );    
   if (!r) { 
     std::cerr << "Error running the HypoTestInverter - Exit " << std::endl;
     return result;
   }
      		

   double ulError = r->UpperLimitEstimatedError();
   double upperLimit = r->UpperLimit();
   std::cout << "The computed upper limit is: " << upperLimit << std::endl;
   std::cout << "an estimated error on this upper limit is: " << ulError << std::endl;
   result.push_back(upperLimit);
   result.push_back(ulError);

   // check using interpolation
   double interpLimit = r->FindInterpolatedLimit(1.-r->ConfidenceLevel() );
   cout << "The computer interpolated limits is " << interpLimit << endl;
   result.push_back(interpLimit);

   const int nEntries = r->ArraySize();


   cout << "List of points " << endl;
   std::vector<Double_t> xArray(nEntries);
   std::vector<Double_t> yArray(nEntries);
   std::vector<Double_t> yErrArray(nEntries);
   for (int i=0; i<nEntries; i++) {
      xArray[i] = r->GetXValue(i);
      yArray[i] = r->GetYValue(i);
      yErrArray[i] = r->GetYError(i);
      std::cout << xArray[i] << " , " << yArray[i] << " err = " << yErrArray[i] << std::endl;
   }
 
   const char *  limitType = (useCls) ? "CLs" : "Cls+b";
   const char * scanType = (npoints < 0) ? "auto" : "grid";

 
   
   HypoTestInverterPlot *plot = new HypoTestInverterPlot("Result","tprime",r);
   TCanvas c1;
   //plot->Draw("CLb 2CL");  // plot all and Clb
   plot->Draw("2CL");  // plot all
   TString resultFileName = TString::Format("%s_%s_ts%d_scan_",limitType,scanType,testStatType);
   resultFileName += suffix;
   resultFileName += ".pdf";
   c1.SaveAs(resultFileName);

  if (plotHypoTestResult) { 
     TCanvas * c2 = new TCanvas();
      c2->Divide( 2, TMath::Ceil(nEntries/2));
      for (int i=0; i<nEntries; i++) {
         c2->cd(i+1);
         SamplingDistPlot * pl = 0;//plot->MakeTestStatPlot(i);
         pl->SetLogYaxis(true);
         pl->Draw();
      }
    }


   // see the expected limit and -1 +1 sigma bands
   SamplingDistribution * limits = r->GetUpperLimitDistribution();
   double * x = const_cast<double *>(&limits->GetSamplingDistribution()[0] );
   double q[5], p[5];
   p[0] = ROOT::Math::normal_cdf(-2,1);
   p[1] = ROOT::Math::normal_cdf(-1,1);
   p[2] = 0.5;
   p[3] = ROOT::Math::normal_cdf(1,1);
   p[4] = ROOT::Math::normal_cdf(2,1);
   TMath::Quantiles(limits->GetSize(), 5 , x, q, p, false);

   std::cout << " expected limit (median) " << q[2] << std::endl;
   std::cout << " expected limit (-1 sig) " << q[1] << std::endl;
   std::cout << " expected limit (+1 sig) " << q[3] << std::endl;
   std::cout << " expected limit (-2 sig) " << q[0] << std::endl;
   std::cout << " expected limit (+2 sig) " << q[4] << std::endl;
   result.push_back(q[0]);
   result.push_back(q[1]);
   result.push_back(q[2]);
   result.push_back(q[3]);
   result.push_back(q[4]);



   if (pWs != NULL) {
     
     // write to a file the results
     //const char *  limitType = (useCls) ? "CLs" : "Cls+b";
     //const char * scanType = (npoints < 0) ? "auto" : "grid";
     TString resultFileName = TString::Format("%s_%s_ts%d_",limitType,scanType,testStatType);      
     //resultFileName += fileName;
     resultFileName += suffix;
     resultFileName += ".root";
     
     TFile * fileOut = new TFile(resultFileName,"RECREATE");
     r->Write();
     fileOut->Close();                                                                     
   }   

   return result;
}


// internal routine to run the inverter
HypoTestInverterResult *  RunInverter(RooWorkspace * w, const char * modelSBName, const char * modelBName, 
                                      const char * dataName, int type,  int testStatType, 
                                      int npoints, double poimin, double poimax, 
                                      int ntoys, bool useCls ) 
{

   std::cout << "Running HypoTestInverter on the workspace " << w->GetName() << std::endl;

   w->Print();


   RooAbsData * data = w->data(dataName); 
   if (!data) { 
      Error("HypoTestDemo","Not existing data %s",dataName);
      return 0;
   }
   else 
      std::cout << "Using data set " << dataName << std::endl;

   
   // get models from WS
   // get the modelConfig out of the file
   ModelConfig* bModel = (ModelConfig*) w->obj(modelBName);
   ModelConfig* sbModel = (ModelConfig*) w->obj(modelSBName);

   if (!bModel) {
      Error("HypoTestDemo","Not existing ModelConfig %s",modelBName);
      return 0;
   }
   if (!sbModel) {
      Error("HypoTestDemo","Not existing ModelConfig %s",modelSBName);
      return 0;
   }
   // check the model 
   if (!sbModel->GetPdf()) { 
      Error("HypoTestDemo","Model %s has no pdf ",modelSBName);
      return 0;
   }
   if (!sbModel->GetParametersOfInterest()) {
      Error("HypoTestDemo","Model %s has no poi ",modelSBName);
      return 0;
   }
   // if (!sbModel->GetSnapshot() ) { 
   //    Error("HypoTestDemo","Model %s has no snapshot ",modelSBName);
   //    return 0;
   // }


   SimpleLikelihoodRatioTestStat slrts(*sbModel->GetPdf(),*bModel->GetPdf());
   if (sbModel->GetSnapshot()) slrts.SetNullParameters(*sbModel->GetSnapshot());
   if (bModel->GetSnapshot()) slrts.SetAltParameters(*bModel->GetSnapshot());

   RatioOfProfiledLikelihoodsTestStat 
      ropl(*sbModel->GetPdf(), *bModel->GetPdf(), sbModel->GetSnapshot());
   ropl.SetSubtractMLE(false);
   
   ProfileLikelihoodTestStat profll(*sbModel->GetPdf());
   if (testStatType == 3) profll.SetOneSided(1);

   TestStatistic * testStat = &slrts;
   if (testStatType == 1) testStat = &ropl;
   if (testStatType == 2 || testStatType == 3) testStat = &profll;
  
   
   HypoTestCalculatorGeneric *  hc = 0;
   if (type == 0) hc = new FrequentistCalculator(*data, *bModel, *sbModel);
   else hc = new HybridCalculator(*data, *bModel, *sbModel);

   ToyMCSampler *toymcs = (ToyMCSampler*)hc->GetTestStatSampler();
   toymcs->SetNEventsPerToy(1);
   toymcs->SetTestStatistic(testStat);


   if (type == 1) { 
      HybridCalculator *hhc = (HybridCalculator*) hc;
      hhc->SetToys(ntoys,ntoys); 
      if (bModel->GetPriorPdf()) hhc->ForcePriorNuisanceAlt(*bModel->GetPriorPdf());
      if (sbModel->GetPriorPdf()) hhc->ForcePriorNuisanceNull(*sbModel->GetPriorPdf());
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
   calc.SetConfidenceLevel(0.95);

   calc.UseCLs(useCls);
   calc.SetVerbose(true);

   // can spped up using proof
   if (useProof && nworkers > 1) { 
      ProofConfig pc(*w, nworkers, "workers=2", kFALSE);
      //ProofConfig pc(*w, 30, "localhost", kFALSE);
      //ToyMCSampler * toymcs = dynamic_cast<ToyMCSampler *> (calc.GetHypoTestCalculator()->GetTestStatSampler() );
      toymcs->SetProofConfig(&pc);    // enable proof
   }

   
   if (npoints > 0) {
      if (poimin >= poimax) { 
      //poimin = TMath::Max(poihat -   4 * poi->getError(), 0.0);
         poimin = int(poihat);
         poimax = int(poihat +  4 * poi->getError());
      }
      std::cout << "Doing a fixed scan  in interval : " << poimin << " , " << poimax << std::endl;
      calc.SetFixedScan(npoints,poimin,poimax);
   }
   else { 
      poi->setMax(10*int( (poihat+ 8 *poi->getError() )/10 ) );
   }

   HypoTestInverterResult * r = calc.GetInterval();


   return r; 
}



void HypoTestInvDemo(const char * fileName ="GausModel_b.root",
                     const char * wsName = "w",
                     const char * modelSBName = "model_sb",
                     const char * modelBName = "model_b",
                     const char * dataName = "data_obs",                  
                     int calculatorType = 0,  // calculator type 
                     int testStatType = 3, // test stat type
                     int npoints = 10,   
                     double poimin = 1,  // use default is poimin >= poimax
                     double poimax = 0,
                     int ntoys=1000,
                     bool useCls = true ) {
   /*
    type = 0 Freq calculator 
    type = 1 Hybrid 

    testStatType = 0 LEP
                 = 1 Tevatron 
                 = 2 PL
                 = 3 PL one sided

   */

   if (fileName==0) { 
      std::cout << "give input filename " << std::endl;
      return; 
   }
   TFile * file = new TFile(fileName); 

   RooWorkspace * w = dynamic_cast<RooWorkspace*>( file->Get(wsName) );
  
   HypoTestInvDemo(w,
		   modelSBName,
		   modelBName,
		   dataName,
		   calculatorType,
		   testStatType,
		   npoints,
		   poimin,
		   poimax,
		   ntoys,
		   useCls);

   return;  
}



int main() {
   HypoTestInvDemo("GausModel_b2.root");
}
