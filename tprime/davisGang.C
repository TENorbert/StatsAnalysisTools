#include "TFile.h"
#include "TObject.h"
#include "TKey.h"
#include "iostream"
#include "fstream"
#include "sstream"
#include "iomanip"
#include "TH1.h"
#include "TH2.h"
#include "vector"
#include "algorithm"
#include "TString.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "THStack.h"

std::vector<std::string> getHistogramNames( std::string, std::string, std::string);
TH2D * getHistogram2D(TFile *file, std::string label);

TString convertDouble(double number)
{
  stringstream ss;//create a stringstream                                                                                                             
  ss << number;//add number to the stream                                                                                                             
  return ss.str();//return a string with the contents of the stream                                                                                   
}

std::string convertInt(int number)
{
  stringstream ss;//create a stringstream                                                                                                             
  ss << number;//add number to the stream                                                                                                            
  return ss.str();//return a string with the contents of the stream                                                                                   
}


void
rebin( int mass,
       double scale,
       double maxErr,
       TString inputRootFilename,
       TString outputRootFilename ){

  //
  // Davis-style bin ganging                                                                                                                         
  //                                                                                                                                                   
  std::cout << "Doing 2D rebinning to 1D, Davis-style\n"
            << "mass = " << mass << std::endl
            << "scale factor (optional) = " << scale << std::endl
            << "max uncertainty = " << maxErr << std::endl
            << "input file with 2D templates: " << inputRootFilename << std::endl
            << "output file: " << outputRootFilename << std::endl;

  //void rebin(  std::string inputRootFilename, std::string outputRootFilename)
  //{

  // First open up map given and get the number of columns
  TString error = convertDouble(maxErr);
  TString mapFilename = "map."+error+".txt";
  std::ifstream mapFile(mapFilename );

  int numberOfcolumns = 0;

  if (mapFile.is_open() )
    {
      std::string firstLine;
      getline( mapFile, firstLine);
      std::stringstream firstLineStream( firstLine );

      std::istream_iterator<std::string> begin(firstLineStream);
      std::istream_iterator<std::string> end;
      std::vector< std::string> vstrings(begin, end);
      
      std::cout<< " Number of columns: " << vstrings.size() << std::endl;
      numberOfcolumns= vstrings.size();
    }

  mapFile.close();

  //----------------------------------------------
  // Read all values from the map, in this case each line is a y vector
  //-----------------------------------------------------------------
  
  // then put all the values in a map of ints (make sure you know which one is X and which one is Y !
  std::vector< std::vector< int > > finalMap;
  std::vector< int > yVector;
  mapFile.open(mapFilename);
  int number;
  int numbers = 0;

  if ( mapFile.is_open() )
    {
      while ( mapFile >> number)
	{
	  yVector.push_back( number );
	  //	  std::cout << number << " " ;
	  numbers++;

	  if ( numbers % numberOfcolumns == 0 )
	    {
	      finalMap.push_back( yVector );
	      //	      std::cout<< std::endl;
	      yVector.clear();
	    }

	}
    }

  // Now let's see what the map is like and maybe count the number of bins and prints the map again (to compare with the one in the ascii file)
  std::vector< std::vector< int > >::iterator  xVector = finalMap.begin();
  int totalNumberOfBins = 1;
 
  for ( ;
	xVector != finalMap.end();
	xVector++)
    {
      std::vector<int>::iterator yVector=xVector->begin();
      for ( ;
	    yVector != xVector->end();
	    yVector++)
	{
	  if ( (*yVector) > totalNumberOfBins )
	    totalNumberOfBins = (*yVector);
	}
    }
  std::cout << " original numbers of bins: " << numbers << " total number of bins after ganging: " << totalNumberOfBins <<  std::endl;

  // Open file with the 2D histograms
  TFile *rootFile = TFile::Open ( inputRootFilename );

  // get the histogram names
  std::vector<std::string> listOfHistogramNames = getHistogramNames( "/", "",convertInt(mass));

  // Create the new output file
  TFile * newFile = new TFile( outputRootFilename, "RECREATE" );
  std::vector<std::string>::iterator nh = listOfHistogramNames.begin();

  double hlower = 0.5;
  double hupper = (double )  totalNumberOfBins + 0.5;

  std::vector< TH1D *> outputHistograms;

  for (;
       nh != listOfHistogramNames.end();
       ++nh)
    {
          
      TString newHistoName( (*nh).c_str() );
      // Adding Ganged to the output histogram to distinguish it from the original (it can probably be removed)
      //      newHistoName += "_Ganged";

      TH1D *gangedHisto = new TH1D( newHistoName.Data(), newHistoName.Data() , totalNumberOfBins, hlower, hupper);

      TH2D *currentHisto = getHistogram2D( rootFile, (*nh) );
      int xNbins = currentHisto->GetXaxis()->GetNbins();
      int yNbins = currentHisto->GetYaxis()->GetNbins();

      for( int xN = 1; xN <= xNbins; ++xN)
	{
	  for (int yN = 1; yN <= yNbins; ++yN )
	    {
	      double content = currentHisto->GetBinContent( xN, yN );
	      double errorSquared   = currentHisto->GetBinError( xN, yN )*currentHisto->GetBinError( xN, yN);

	      content += gangedHisto->GetBinContent( finalMap.at(xN-1).at(yN-1) );
	      double errorInBinSquared = gangedHisto->GetBinError( finalMap.at(xN -1).at( yN-1) )*gangedHisto->GetBinError( finalMap.at(xN-1).at(yN-1) ) ;
	      
	      gangedHisto->SetBinContent( finalMap.at(xN-1).at(yN-1), content );
	      gangedHisto->SetBinError( finalMap.at(xN-1).at(yN-1), sqrt( errorSquared + errorInBinSquared) );
	    }
	}
      
      newFile->cd();
      gangedHisto->Write();
      outputHistograms.push_back( gangedHisto );
    }
         
  rootFile->Close();
  newFile->Close();
  
}

//---------------------------------------------------------------------------------------------
// --------------         Helper Functions ---------------------------------------------------
//---------------------------------------------------------------------------------------------

std::vector<std::string> getHistogramNames( std::string directory, std::string matched , std::string mass)
// GETS ALL THE HISTOGRAM NAMES FROM THE DIRECTORY GIVEN WITH THE STRING MATCHED GIVEN
{
  std::vector<std::string> listOfHistograms;

  gDirectory->cd(directory.c_str());
  size_t found;
  size_t foundtp;

  TIter next(gDirectory->GetListOfKeys());
  TKey *key;

  while ((key= (TKey*) next())){

    std::string nameHisto = key->GetName();
    if ( matched.size() == 0 ) {
      foundtp = nameHisto.find(mass);
      if ((foundtp != string::npos ||  nameHisto.find("prime")==string::npos )){
	listOfHistograms.push_back( nameHisto);
	std::cout<<nameHisto<<std::endl;
      }
    }
    else {
      found = nameHisto.find(matched);
      if (found != string::npos ){//  (foundtp != string::npos ||  nameHisto.find("prime")==string::npos )){
	listOfHistograms.push_back( nameHisto);
      }
    }
  }

  gDirectory->cd();
  return listOfHistograms;
}

//---------------------------------------------------------------------------------------------

TH2D * getHistogram2D(TFile *file, std::string label)
// Gets the histogram from the file with a certain label
{
  
  file->cd();

  TH2D * histogram = (TH2D *) file->Get( label.c_str() );
  
  if ( histogram == NULL )
    cout << " Histogram: " << label.c_str() << " not found. " << std::endl;

  return histogram;

}
