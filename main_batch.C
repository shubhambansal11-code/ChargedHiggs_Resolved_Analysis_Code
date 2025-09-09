#include "main/EventLoop.h"
#include "TApplication.h"
#include <TSystemDirectory.h>
#include <TList.h>
#include <TCollection.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1F.h>
#include <TTree.h>
#include <TString.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <ctime>

//int main(int argc, char *argv[]){
int main(int argc, char** argv){
  clock_t begin = clock();
  /*TApplication theApp("evsel", &argc, argv);
  TString path          =  theApp.Argv(1);
  TString SampleName    =  theApp.Argv(2);
  TString WP            =  theApp.Argv(3);
  TString OUTPUTDIR     =  theApp.Argv(4);*/

  TString path          =  TString(argv[1]);
  TString SampleName    =  TString(argv[2]);
  TString WP            =  TString(argv[3]);
  TString OUTPUTDIR     =  TString(argv[4]);
  //bool    batchMode     =  atoi(theApp.Argv(5));
  
  TString MCDataPeriode = "";
  std::vector<std::string> TreeNames;
  TreeNames.push_back("Nominal");
  if(path.Contains("MC16a"))MCDataPeriode = "MC16a";
  if(path.Contains("MC16d"))MCDataPeriode = "MC16d";
  if(path.Contains("MC16e"))MCDataPeriode = "MC16e";
  TString OutFile = SampleName+"_"+WP+"_"+MCDataPeriode+".root";
  
  //OutFile.ReplaceAll(".root", "_"+WP+"_"+MCDataPeriode+".root");
  //TString OutFileName = OUTPUTDIR +"/Plotfiles/"+OutFile;
  TString OutFileName = OUTPUTDIR + "/"+ OutFile;
  //std::cout<<"Outfile name"<<OutFileName<<std::endl; 
  TFile* outfile = TFile::Open(OutFileName,"RECREATE");
  for(unsigned int i=0; i < TreeNames.size(); i++){
     TString TreeName  = TString(TreeNames.at(i));
      std::cout<<"Tree name"<<TreeName<<std::endl;
     TChain* mych_data = new TChain(TreeName);
     mych_data->Add(OUTPUTDIR+"/"+path+"/"+SampleName);
     std::cout<<mych_data->GetEntries()<<"  "<<path+"/"+SampleName<<std::endl;
     EventLoop *eventLoop = new EventLoop(mych_data,SampleName,MCDataPeriode,TreeName,WP);
     eventLoop->Loop();
     TDirectory* subdir = outfile->mkdir(TreeNames.at(i).c_str());
     std::cout<<"Sub directory"<<subdir<<std::endl;
     std::cout<<"Tree name.at(i)"<<""<<TreeNames.at(i)<<std::endl;
     eventLoop->Write(subdir,TreeNames.at(i));
     delete mych_data;
     delete eventLoop;
  }

  std::cout<<"Finished looping over the events"<<std::endl;  
  outfile->Close();
  return 0;
}