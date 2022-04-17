#include "includes/ana3.h"
#include "includes/ana3.C"

#include "includes/HistManager.h"
#include "includes/AnaPars.h"

#define inFile_f "../../data/step3/%s_%s_s3.root"
#define outFile_f "../../data/step4/%s_s4.root"
#define inTreeName "CRT"
#define outTreeName "CRT"

using namespace std;
HistManager HM;

void CRT_step4(TString runName) {

  HM.AddHistBox("th1f", 1, "cosTheta", "costTheta", "cos", "", 500, 0, 1);
  HM.AddHistBox("th1f", 1, "cos", "costTheta", "cos", "", 500, 0, 1);

  TString inFileNameT = Form(inFile_f, runName.Data(), "T");
  TString inFileNameB = Form(inFile_f, runName.Data(), "B");
  TString outFileName = Form(outFile_f, runName.Data());

  TFile *inFileT = new TFile(inFileNameT), *inFileB = new TFile(inFileNameB);
  TTree *inTreeT = (TTree*)inFileT->Get(inTreeName), *inTreeB = (TTree*)inFileB->Get(inTreeName);
  Long64_t maxEvT = inTreeT->GetEntriesFast(), maxEvB = inTreeB->GetEntriesFast();

  cout<<endl<<"----> runName : "<<runName<<endl;
  cout<<"----> Input files : "<<inFileNameT<<" , "<<inFileNameB<<endl;
  cout<<"----> Output file : "<<outFileName<<endl;
  cout<<"----> Entries : "<<maxEvB<<" + "<<maxEvT<<endl;

  TFile *outFile = new TFile(outFileName, "RECREATE");
  outFile->cd();
  TTree* CRTs4 = new TTree(outTreeName, outTreeName);          
  CRTs4->SetAutoSave(1000);
  Long64_t iTrig_out;
  int iSc_out[2];
  double Z_out[2], Q_out[4], T_out[4], X_out[4];
  CRTs4->Branch("iTrig", &iTrig_out, "iTrig/L");
  CRTs4->Branch("iScint", &iSc_out, "iScint/I");
  CRTs4->Branch("Z", &Z_out, "Z/D");
  CRTs4->Branch("Q", &Q_out, "Q[4]/D");
  CRTs4->Branch("T", &T_out, "T[4]/D");
  CRTs4->Branch("X", &X_out, "X[4]/D");

  ana3 anaT(inTreeT), anaB(inTreeB);

  Long64_t jT{0}, jB{0};
  while(1) {

    if (anaT.LoadTree(jT) < 0 || anaB.LoadTree(jB) < 0) {break;}
    anaT.fChain->GetEntry(jT); anaB.fChain->GetEntry(jB); 
    Long64_t iEve = anaT.iTrig;   
    if (iEve > anaB.iTrig) { jB++; continue; }
    if (iEve < anaB.iTrig) { jT++; continue; }
    jB++; jT++;

    iTrig_out = iEve;
    Z_out[0] = anaT.Z; Z_out[1] = anaB.Z;
    iSc_out[0] = anaT.iSc; iSc_out[1] = anaB.iSc;
    for (int i = 0; i < 2; i++) {
      Q_out[i] = anaT.Q[i];
      Q_out[i+2] = anaB.Q[i];
      T_out[i] = anaT.T[i];
      T_out[i+2] = anaB.T[i];
    }

    //aggiungere tagli in tempo top/bottom e distribuzioni tempi top bottom-top
    //AGGIUNGERE HISTO DI DIAGNOSTICA (hit vs iscint e modulo, cariche, zeta)
    //METTERE CORREZIONI Z_offset da CSV
    //DOMANI LO FINISCO E MARTEDì SI PROVA

    CRTs4->Fill();

  }
  
  CRTs4->Write();
  outFile->Close();

}





