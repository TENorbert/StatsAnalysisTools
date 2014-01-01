//v1 mmhl 260811 - template sanity check - 1D at the moment
#include "TTree.h"
#include "TH1F.h"
#include "TFile.h"
#include "TString.h"
#include <math.h>
using namespace std;
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <stdarg.h>

int main(int argc, char **argv){
  if(argc != 2){ cout<<" only one argument! "<<endl; return 1;}
  TString xmlfilename = argv[argc-1];  
  ifstream xmlfile;
  xmlfile.open(xmlfilename);
  if(!xmlfile){
    cout<<"--file "<<xmlfilename<<" not found!"<<endl; return 1;
  }
  
  int inputcounter = 0;   
  string templine ; int filenamelength; 
  char buffer[128]; size_t length; 
  
  ifstream emuxmlfile;
  TString emuxmlfilename;
  TFile *rootfile;
  
  while(!xmlfile.eof()){
    getline(xmlfile,templine);
    if(inputcounter==0 && templine.find("<Input>") != string::npos &&  templine.find("</Input>") != string::npos){
      inputcounter++;

      filenamelength =  templine.find("</Input>") -  templine.find("<Input>") - 7; //6 is length of "</Input>"
      length=templine.copy(buffer,filenamelength,templine.find("<Input>")+7);
      buffer[length]='\0';
      cout<<" opening 1st  file    listed: " << string(buffer)<<endl;
      emuxmlfilename = string(buffer);
      emuxmlfile.open(emuxmlfilename);
      if(!emuxmlfile){cout<<"--file "<<string(buffer)<<" not found!"<<endl; return 1;}
    }
    
  }
  
  //now get into channel xml file and store histograms
  if(inputcounter==0) {cout<<"no input file found"<<endl; return 1;}
  inputcounter=0;
  int incounter[4]={0};
  bool skipper = false;
  int histcount=0;
  const int numhists = 100;
  TString histnames[numhists] ;
  while(!emuxmlfile.eof()){
    getline(emuxmlfile,templine);
    //this avoids the commented region
    if(templine.find("-->")!=string::npos) skipper = false;
    if(templine.find("<!--")!=string::npos){
      skipper = true;
    }
    if(skipper)continue;

    if(inputcounter==0 && templine.find("InputFile") != string::npos && templine.find(".root")!=string::npos){
      inputcounter++;
      filenamelength =  templine.find(".root") -  templine.find("InputFile") - 6; 
      length=templine.copy(buffer,filenamelength,templine.find("InputFile")+11);
      buffer[length]='\0';
      cout<<" opening 1st rootfile listed: " << string(buffer)<<endl;
      emuxmlfilename = string(buffer);
      rootfile = new TFile(emuxmlfilename);

      if(rootfile->IsZombie()){
	cout<<" Error opening file "<<endl;
	delete rootfile;
	return 1;
      }
    }

    
    if(incounter[0]==0  && templine.find("\"Data")!=string::npos && templine.find("HistoPath")!=string::npos){
      incounter[0]++;
      filenamelength = templine.find("HistoPath=") - templine.find("\"Data") -3;
      length=templine.copy(buffer,filenamelength,templine.find("\"Data")+1);
      buffer[length]='\0';
      cout<<" data template name         : "<<string(buffer)<<endl;
      histnames[histcount]=string(buffer);
      histcount+=2;
    }
    
    if(inputcounter!=0 && incounter[1]==0 && templine.find("signal")
       && (templine.find("\"Tprime")!=string::npos || templine.find("\"TPrime")!=string::npos)){
      incounter[1]++;
      if(templine.find("\"Tprime")!=string::npos){
	filenamelength = templine.find("\">") - templine.find("\"Tprime") -1;
	length=templine.copy(buffer,filenamelength,templine.find("\"Tprime")+1);
      }
      else if(templine.find("\"TPrime")!=string::npos){
	filenamelength = templine.find("\">") - templine.find("\"TPrime") -1;
	length=templine.copy(buffer,filenamelength,templine.find("\"TPrime")+1);
      }
      buffer[length]='\0';
      cout<<" signal template name       : "<<string(buffer)<<endl;
      histnames[histcount]=string(buffer);
      
      getline(emuxmlfile,templine);
      if(templine.find("Sample")!=string::npos) continue;
      if(templine.find("jes")   !=string::npos && templine.find("High=") !=string::npos && templine.find("Low=") !=string::npos){
	filenamelength = templine.find("HistoNameLow") - templine.find("High=")-8;
	length= templine.copy(buffer,filenamelength,templine.find("High=")+6);
	buffer[length]='\0';
	cout<<" --Jes Up template name     : "<<string(buffer)<<endl;
	histnames[histcount+1]=string(buffer);

	filenamelength = templine.find("\"/>") - templine.find("Low=")-5;
	length= templine.copy(buffer,filenamelength,templine.find("Low=")+5);
	buffer[length]='\0';
	cout<<" --Jes Down template name   : "<<string(buffer)<<endl;
	histnames[histcount-1]=string(buffer);
      }
      histcount+=3;
    }

    if(inputcounter!=0 && incounter[2]==0 && templine.find("top")
       &&  templine.find("\"Top")!=string::npos){
      incounter[2]++;
      filenamelength = templine.find("\">") - templine.find("\"Top") -1;
      length=templine.copy(buffer,filenamelength,templine.find("\"Top")+1);
      buffer[length]='\0';
      cout<<" top template name          : "<<string(buffer)<<endl;
      histnames[histcount]=string(buffer);

      getline(emuxmlfile,templine);
      if(templine.find("Sample")!=string::npos) continue;
      if(templine.find("jes")   !=string::npos && templine.find("High=") !=string::npos && templine.find("Low=") !=string::npos){
	filenamelength = templine.find("HistoNameLow") - templine.find("High=")-8;
	length= templine.copy(buffer,filenamelength,templine.find("\"Top")+1);
	buffer[length]='\0';
	cout<<" --Jes Up template name     : "<<string(buffer)<<endl;
	histnames[histcount+1]=string(buffer);

	filenamelength = templine.find("\"/>") - templine.find("Low=")-5;
	length= templine.copy(buffer,filenamelength,templine.find("Low=")+5);
	buffer[length]='\0';
	cout<<" --Jes Down template name   : "<<string(buffer)<<endl;
	histnames[histcount-1]=string(buffer);
      }
      histcount+=3;
    }

    if(inputcounter!=0 && incounter[3]==0 && templine.find("ewk")
       &&  templine.find("\"Ewk")!=string::npos){
      incounter[3]++;
      filenamelength = templine.find("\">") - templine.find("\"Ewk") -1;
      length=templine.copy(buffer,filenamelength,templine.find("\"Ewk")+1);
      buffer[length]='\0';
      cout<<" SM  template name          : "<<string(buffer)<<endl;

      histnames[histcount]=string(buffer);
      getline(emuxmlfile,templine);
      if(templine.find("Sample")!=string::npos) continue;
      if(templine.find("jes")   !=string::npos && templine.find("High=") !=string::npos && templine.find("Low=") !=string::npos){
	filenamelength = templine.find("HistoNameLow") - templine.find("High=")-8;
	length= templine.copy(buffer,filenamelength,templine.find("\"Ewk")+1);
	buffer[length]='\0';
	cout<<" --Jes Up template name     : "<<string(buffer)<<endl;
	histnames[histcount+1]=string(buffer);

	filenamelength = templine.find("\"/>") - templine.find("Low=")-5;
	length= templine.copy(buffer,filenamelength,templine.find("Low=")+5);
	buffer[length]='\0';
	cout<<" --Jes Down template name   : "<<string(buffer)<<endl;
	histnames[histcount-1]=string(buffer);
      }
      // histcount+=3;
    }
  }

  
  
  int masspt = 1;
  TH1F *hist[numhists];
  double entrycounter[numhists];
  double overflowcounter[numhists];
  double bgcounter[3]; //3 for jes95,100,105
  double bgoverflowcounter[3];
  cout.setf(ios::fixed,ios::floatfield);
  cout.precision(5);
  
  for(int i = 0 ; i < 3 ; i++){
    bgcounter[i]=0.; //resets bg counter
    bgoverflowcounter[i]=0.;
  }
  cout<<endl;
  //cout<<"filename  *** "<<filename<<" *** "<<endl;
  cout<<"template name         "<<'\t'<< "|| "<<'\t'<<"JES095 nev"<<'\t'<<"over+underflow"<<'\t'<<"||"<<'\t'<<'\t'
      << "|| "<<'\t'<<"NOM nev   "<<'\t'<< "over+underflow"<<'\t'<<"||"<<'\t'<<'\t'
      << "|| "<<'\t'<<"JES105 nev"<<'\t'<<"over+underflow"<<'\t'<<"||"<<endl;
   
  cout<<"----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------"<<endl;
  if(masspt > 6) return 1;
  
  // cout
  
  for(int l = 0 ; l < histcount+2 ; l++){
    if(l==1 || l ==4 ) cout<<"----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------"<<endl;
    if((l-1)%3==0 || l==0){
      if(l==0)cout<<histnames[l]<<endl;
      else{cout<<histnames[l+1]<<endl;}
      
      cout<<"====>"<<'\t'<<'\t'<<'\t';
    }
    hist[l] = (TH1F*)rootfile->Get(histnames[l]);
    entrycounter[l]=0.;
    overflowcounter[l]=0.;
    for(int ib = 0 ; ib < hist[l]->GetSize(); ib++){
      entrycounter[l] += hist[l]->GetBinContent(ib);
      //cout<<ib<<" "<< hist[l]->GetBinContent(ib)<<" "<< hist[l]->GetSize()<<endl;
      if(ib==0 || ib==hist[l]->GetSize()-1) overflowcounter[l]+= hist[l]->GetBinContent(ib);
    }
    
    //saving the bg together
    if(l>=4){
      bgcounter[l%3]         +=entrycounter[l];
      bgoverflowcounter[l%3] += overflowcounter[l];
    }
    if(l!=0 )cout<<"|| "<<'\t'<<entrycounter[l]<< '\t'<<overflowcounter[l]<<'\t'<<'\t'<<"|| "<<'\t'<<'\t';
    //else if(l > 6) cout<<"|| "<<entrycounter[l]<<" ("<<entrycounter[l]+bgcounter[l%3] <<")"<< " "<<overflowcounter[l]<<'\t'<<'\t'<<"|| "<<'\t';
 
    else{
      cout<<"|| "<<'\t'<<entrycounter[l]<< '\t'<<overflowcounter[l]<<'\t'<<'\t'<<"|| "<<'\t'<<'\t';
      cout<<"|| "<<'\t'<<entrycounter[l]<< '\t'<<overflowcounter[l]<<'\t'<<'\t'<<"|| "<<'\t'<<'\t';
      cout<<"|| "<<'\t'<<entrycounter[l]<< '\t'<<overflowcounter[l]<<'\t'<<'\t'<<"|| "<<'\t'<<'\t';
    }
    if(l%3==0)cout<<endl;
    
    //outputs total bg
    if(l==histcount+1){
      cout<<"----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------"<<endl;
      cout<<"BG_fitMass  "<<endl;
      cout<<"====>" <<'\t'<<'\t'<<'\t'<<"|| "<<'\t'<<bgcounter[1]<< '\t'<<bgoverflowcounter[1]<<'\t'<<'\t'<<"|| "<<'\t'<<'\t';
      cout<<"|| "<<'\t'<<bgcounter[2]<< '\t'<<bgoverflowcounter[2]<<'\t'<<'\t'<<"|| "<<'\t'<<'\t';
      cout<<"|| "<<'\t'<<bgcounter[0]<< '\t'<<bgoverflowcounter[0]<<'\t'<<'\t'<<"|| "<<'\t'<<endl;
    }
    
  }

  cout<<endl;
  delete rootfile;
}
