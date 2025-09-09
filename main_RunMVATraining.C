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
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TGraph.h>

#include "TMVA/Factory.h"
#include "TMVA/Reader.h"
#include "TMVA/DataLoader.h"


void TrainMVA(TChain *fChain);
template <typename DataType> void RunCutOptimisation(TChain *fChain, TString observable);

int main(int argc, char** argv){
  clock_t begin = clock();
  TApplication theApp("evsel", &argc, argv);
  TString TreeName   =  theApp.Argv(1);
  TString SampleName =  theApp.Argv(2);
  ///TString observable =  theApp.Argv(3);  
  ///TString type       =  theApp.Argv(4);

  TChain* mych_data  = new TChain ("Nominal/"+TreeName);
  //TChain* mych_data  = new TChain (TreeName);
  mych_data->Add("inputs/MVATraining/"+SampleName);
  TrainMVA(mych_data);
  /*if(type=="int"){
     RunCutOptimisation<int>(mych_data,observable);
  }else if(type=="double"){
     RunCutOptimisation<double>(mych_data,observable);
  }*/
  return 0;
}

template <typename DataType> void RunCutOptimisation(TChain *fChain, TString observable){
    DataType m_observable;
    int      is_Signal;
    int      step_size = 20;
    double   x_range   = 20.;
    fChain->SetBranchAddress(observable, &m_observable);
    fChain->SetBranchAddress("is_Signal", &is_Signal);
    std::vector<int> selected_events_sig(step_size, 0);
    std::vector<int> selected_events_bkg(step_size, 0);
    for (Long64_t jentry=0; jentry<fChain->GetEntriesFast();jentry++) {
        fChain->GetEntry(jentry);   
        for(unsigned int i=0; i<step_size; i++){
            DataType m_cut = (DataType) ((i)*(x_range/step_size));
            if(is_Signal==1){
                if(m_observable >= m_cut && m_observable < m_cut + 1)selected_events_sig[i]++;
            }else if(is_Signal==0){
                if(m_observable >= m_cut && m_observable < m_cut + 1)selected_events_bkg[i]++;
            }
        }
    }
    /// find maximum 
    for(unsigned int i=0; i<step_size; i++){
        DataType m_cut      = (DataType) ((i)*(x_range/step_size));
        double significance = selected_events_sig[i]/sqrt(selected_events_bkg[i]);
        std::cout<< significance << "for >=" << m_cut <<std::endl;
    }
}

void TrainMVA(TChain *fChain){
   TFile *outputFile            = TFile::Open("output/TMVAResults_qqbb_incEtaHW.root", "RECREATE" );
   TMVA::Factory *factory	= new TMVA::Factory("TMVAClassificationCategory", outputFile, "!V:!Silent:Transformations=I;P;G");
   TMVA::DataLoader *dataloader = new TMVA::DataLoader("dataset");

   /*
   dataloader->AddVariable("nJet", "nJet","",'D');
   dataloader->AddVariable("nBJet", "nBJet","",'D');
   dataloader->AddVariable("HT_bjets", "HT_bjets","",'D');
   dataloader->AddVariable("maxEta_bjets", "maxEta_bjets","",'D');
   dataloader->AddVariable("maxPT_bjets", "maxPT_bjets","",'D');
   dataloader->AddVariable("Wleptonic_pT", "Wleptonic_pT","",'D');
   dataloader->AddVariable("Wleptonic_Eta", "Wleptonic_Eta","",'D');
   dataloader->AddVariable("Lep_PT", "Lep_PT","",'D');
   dataloader->AddVariable("MET", "MET","",'D');
   dataloader->AddVariable("Lepton_Charge", "Lepton_Charge","",'D');
   */

   dataloader->AddVariable("btagjH1", "btagjH1","",'D');
   dataloader->AddVariable("btagjH2", "btagjH2","",'D');
   dataloader->AddVariable("btagjW1", "btagjW1","",'D');
   dataloader->AddVariable("btagjW2", "btagjW2","",'D'); 
   dataloader->AddVariable("H_mass", "H_mass","",'D');
   dataloader->AddVariable("Wp_mass", "Wp_mass","",'D');
   dataloader->AddVariable("Phi_HW", "Phi_HW","",'D');
   dataloader->AddVariable("Eta_HW", "Eta_HW","",'D');
   dataloader->AddVariable("H_pT/mass_VH",   "H_pT/mass_VH","",'D');
   dataloader->AddVariable("Wp_pT/mass_VH",   "Wp_pT/mass_VH","",'D');
   //dataloader->AddVariable("pTjH1", "pTjH1","",'D');
   //dataloader->AddVariable("pTjH2", "pTjH2","",'D');
   //dataloader->AddVariable("pTjWp1", "pTjWp1","",'D');
   //dataloader->AddVariable("pTjWp2", "pTjWp2","",'D');
   //dataloader->AddVariable("dRjjH", "dRjjH","",'D');
   // dataloader->AddVariable("dRjjWp", "dRjjWp","",'D');

   /*m_myTree->Branch("H_pT",           &m_H_pT);
   m_myTree->Branch("Wp_pT",            &m_Wp_pT);
   m_myTree->Branch("dRjjH",            &m_dRjjH);
   m_myTree->Branch("dRjjWp",           &m_dRjjWp);
   m_myTree->Branch("mass_VH",          &m_mass_VH);
   m_myTree->Branch("is_Signal",        &m_is_Signal);
   m_myTree->Branch("EventWeight",	&EventWeight); */

   std::cout<<"Going to read input tree"<<std::endl;
   TTree *InTree = fChain->CloneTree(0);

   std::cout<<"Going to start Training !!!"<<std::endl;
   float signalWeight     = 1.0;
   float backgroundWeight = 1.0;

   dataloader->AddSignalTree    (fChain,     signalWeight);
   dataloader->AddBackgroundTree(fChain,     backgroundWeight);
   ///TCut mySigDef = "is_Signal == 1 &&  nJet > 4 && nBJet > 1";
   ///TCut myBkgDef = "is_Signal == 0 &&  nJet > 4 && nBJet > 1";

   TCut mySigDef = "is_Signal == 1 && H_mass < 250. && Wp_mass < 160";
   TCut myBkgDef = "is_Signal == 0 && H_mass < 250. && Wp_mass < 160";

   dataloader->PrepareTrainingAndTestTree(mySigDef, myBkgDef, "SplitMode=Random:NormMode=NumEvents:!V:nTrain_Signal=0:nTest_Signal=0:nTrain_Background=0:nTest_Background=0");
   factory->BookMethod(dataloader,TMVA::Types::kBDT, "WpH_Tagger_qqbb_etaHWinc","!H:!V:NTrees=600:MinNodeSize=2.0%:nCuts=40:BoostType=Grad:Shrinkage=0.1:UseBaggedBoost:GradBaggingFraction=0.25:SeparationType=GiniIndex:DoBoostMonitor:MaxDepth=5:NegWeightTreatment=IgnoreNegWeightsInTraining");
   factory->TrainAllMethods();
   factory->TestAllMethods();
   factory->EvaluateAllMethods();
   std::cout<<"Finished Training !!!"<<std::endl;
   outputFile->Close();
}
