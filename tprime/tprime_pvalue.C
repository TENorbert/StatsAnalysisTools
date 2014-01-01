//////////////////////////////////////////////////////////////////////////
//
// tprime_pvalue.C
//
// Routine for getting p-value out of a RooStats model
//
//
// Gena Kukartsev
//
// December 2011: first version
// 

#include <iostream>
#include "TFile.h"
#include "TStopwatch.h"
#include "RooWorkspace.h"
#include "RooAbsData.h"
#include "RooRealVar.h"
#include "RooArgSet.h"
#include "RooRandom.h"
#include "RooStats/ModelConfig.h"
#include "RooStats/ProfileLikelihoodTestStat.h"
#include "RooStats/ToyMCSampler.h"
#include "RooStats/FrequentistCalculator.h"



double GetPValue(std::string inFileName,
		 std::string workspaceName = "ejets",
		 std::string modelConfigName = "ModelConfig",
		 std::string dataName = "obsData"){

  // Check if input file exists
  TFile * pFile = TFile::Open( inFileName.c_str() );

  // if input file was specified but not found, quit
  if( !pFile ){
    std::cout << "file not found" << std::endl;
    return -1.0;
  }

  // get the workspace out of the file
  RooWorkspace * pWs = (RooWorkspace *) pFile->Get( workspaceName.c_str() );
  if(!pWs){
    std::cout << "workspace not found" << std::endl;
    return -1.0;
  }

  // get the modelConfig out of the file
  RooStats::ModelConfig * pSbModel = (RooStats::ModelConfig *) pWs->obj( modelConfigName.c_str() );

  // get the modelConfig out of the file
  RooAbsData * pData = pWs->data( dataName.c_str() );

  // make sure ingredients are found
  if(!pData || !pSbModel){
    pWs->Print();
    std::cout << "data or ModelConfig was not found" << std::endl;
    return -1.0;
  }

  // make b-only model from SbModel by copying and setting POI to zero
  RooStats::ModelConfig * pBModel = (RooStats::ModelConfig *)pSbModel->Clone();
  pBModel->SetName("BModel");      
  RooRealVar * var = dynamic_cast<RooRealVar *>(pBModel->GetParametersOfInterest()->first());
  if (!var) return -1.0;
  double oldval = var->getVal();
  var->setVal(0);
  pBModel->SetSnapshot( RooArgSet(*var)  );
  var->setVal(oldval);

  // POI
  const RooArgSet * poiSet = pSbModel->GetParametersOfInterest();
  RooRealVar *poi = (RooRealVar*)poiSet->first();
  
  std::cout << "POI initial value:   " << poi->GetName() << " = " << poi->getVal()   << std::endl;  
  
  // test statistic
  RooStats::ProfileLikelihoodTestStat profll(*pSbModel->GetPdf());
  //profll.SetOneSided(1);
  profll.SetOneSided(0);
  std::string minimizerType = ROOT::Math::MinimizerOptions::DefaultMinimizerType();
  profll.SetMinimizer( minimizerType.c_str() );
  //profll.SetPrintLevel(mPrintLevel);
  profll.SetReuseNLL(true);
  profll.SetStrategy(0);
  RooStats::TestStatistic * testStat = &profll;

  // time
  TStopwatch tw;

  // create calculator
  //RooStats::HypoTestCalculatorGeneric *  pHtc = new RooStats::FrequentistCalculator(*pData, *pBModel, *pSbModel);
  RooStats::HypoTestCalculatorGeneric *  pHtc = new RooStats::FrequentistCalculator(*pData, *pSbModel, *pBModel);

  // get MC sampler
  RooStats::ToyMCSampler * toymcs = (RooStats::ToyMCSampler*)pHtc->GetTestStatSampler();
  if (toymcs) { 

    //toymcs->SetNEventsPerToy(1);

    toymcs->SetTestStatistic(testStat);
    
    toymcs->SetGenerateBinned(true);
    
    toymcs->SetUseMultiGen(true);
    
    RooRandom::randomGenerator()->SetSeed(0); 
    
  }

  //pWs->Print();

  ((RooStats::FrequentistCalculator*) pHtc)->SetToys(100,0);

  double pvalue = pHtc->GetHypoTest()->NullPValue();
  //double pvalue_err = pHtc->GetHypoTest()->NullPValueError();
  double pvalue_err = 0;

  std::cout << "p-value: " << pvalue 
	    << " +/- " << pvalue_err << std::endl;

  //std::cout << "all cool" << std::endl;

  pFile->Close();

  return pvalue;
}
