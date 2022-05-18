#include "includes/ana3.h"
#include "includes/ana3.C"

#include "includes/HistManager.h"
#include "includes/AnaPars.h"
#include "includes/CsvHandler.h"


#define inFile_f "../../data/step3/%s_%s_s3.root"
#define outFile_f "../../data/step4/%s_s4.root"
#define inTreeName "CRT"
#define outTreeName "CRT"

#define crtModulesDistance 10.0
#define crtOffsetBottomX -0.0

using namespace std;
HistManager HM;
CsvHandler CSV;

void CRT_step4(TString runName) {

  if (runName == "") {cout<<"runName is empty!!"<<endl; return;}

  double **zetaCalT = CSV.InitMatrix(1, scintNum); 
  double **zetaCalB = CSV.InitMatrix(1, scintNum); 
  CSV.Read(CSV.GetFirstFile(lutPrefix3p + runName + "_T" + lutChEqName + "*"),      ',', zetaCalT, 1, scintNum);
  CSV.Read(CSV.GetFirstFile(lutPrefix3p + runName + "_B" + lutChEqName + "*"),      ',', zetaCalB, 1, scintNum);

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
  double Z_out[2], Q_out[4], T_out[4], X_out[2];
  CRTs4->Branch("iTrig", &iTrig_out, "iTrig/L");
  CRTs4->Branch("iScint", &iSc_out, "iScint/I");
  CRTs4->Branch("Z", &Z_out, "Z[2]/D");
  CRTs4->Branch("X", &X_out, "X[2]/D");
  CRTs4->Branch("Q", &Q_out, "Q[4]/D");
  CRTs4->Branch("T", &T_out, "T[4]/D");

  HM.SetOutFile(outFile);
  HM.SetUseFolders(0);
  HM.AddHistBox("th1f", 1, "thetaR", "thetaR", "", "", 500, -pi, pi);
  HM.AddHistBox("th1f", 1, "thetaZYvert", "thetaZY", "", "", 500, -1.5, 1.5); 
  HM.AddHistBox("th1f", 1, "thetaZY", "thetaZY", "", "", 500, -1.5, 1.5);
  HM.AddHistBox("th1f", 1, "thetaXY", "thetaXY", "", "", 500, -pi/2, pi/2);
  HM.AddHistBox("th2f", 1, "heatMap", "heatMap", "Z", "cm", "iScint + iMod*(1 + scintNo)", "", 320, -scintL, scintL, 2*scintNum + 2, - 1, 2*scintNum + 1 );
  HM.AddHistBox("th1f", 1, "mtDiff", "mtDiff", "T", "ns", 1200, -60, 60);
  HM.AddHistBox("th1f", 1, "mtTop", "mtT", "T", "ns", 1200, -60, 60);
  HM.AddHistBox("th1f", 1, "mtBtm", "mtB", "T", "ns", 1200, -60, 60);
  HM.AddHistBox("th1f", 1, "mintDiffCorr", "mintDiffCorr", "T", "ns", 1200, -60, 60);
  HM.AddHistBox("th2f", 1, "zetaTop_Btm", "Z", "Zbtm", "cm", "Ztop", "cm", 640, -scintL, scintL, 640, -scintL, scintL);
  HM.AddHistBox("th2f", 1, "xTop_Btm", "x", "xbtm", "cm", "xtop", "cm", scintNum, 0, scintNum, scintNum, 0, scintNum);
  HM.SetUseFolders(1);
  HM.AddHistBox("th1f", 4, "zetaTopPerThetaXY", "Z", "Z", "cm", 640, -scintL, scintL);
  HM.AddHistBox("th1f", 2*scintNum, "zetaPerScint", "Z", "Z", "cm", 640, -scintL, scintL);

  ana3 anaT(inTreeT), anaB(inTreeB);

  Long64_t jT{0}, jB{0};
  while(1) {

    //senza dover abilitare riga 44 in step3
    if (anaT.LoadTree(jT) < 0 || anaB.LoadTree(jB) < 0) {break;}
    anaT.fChain->GetEntry(jT); anaB.fChain->GetEntry(jB); 
    iTrig_out = anaT.iTrig;   
    if (iTrig_out > anaB.iTrig) { jB++; continue; }
    if (iTrig_out < anaB.iTrig) { jT++; continue; }
    jB++; jT++;

    // //con riga 44 abilitata in step3
    // if (anaT.LoadTree(jT) < 0) {break;}
    // anaT.fChain->GetEntry(jT);
    // iTrig_out = anaT.iTrig; 
    // if (anaB.LoadTree(anaB.fChain->GetEntryNumberWithIndex(iTrig_out)) < 0) {continue;}
    // anaB.fChain->GetEntryWithIndex(iTrig_out);
    // jT++;

    double xT = 2.5*((double)anaT.iSc - 0.5*((double)scintNum - 1)); //mettere in pars --> scintW = 2.5 se non c'è già
    double xB = 2.5*((double)anaB.iSc - 0.5*((double)scintNum - 1)) + crtOffsetBottomX;
    //double zT = anaT.pZ, zB = anaB.pZ; // sono sistemate le pseudoZeta????
    double zT = anaT.Z, zB = anaB.Z; ////pZ
    double mtT{0}, mtB{0};

    iSc_out[0] = anaT.iSc; iSc_out[1] = anaB.iSc;
    Z_out[0] = zT; Z_out[1] = zB;
    X_out[0] = xT; X_out[1] = xB;
    for (int i = 0; i < 2; i++) {
      Q_out[i] = anaT.Q[i];
      Q_out[i+2] = anaB.Q[i];
      T_out[i] = anaT.T[i];
      T_out[i+2] = anaB.T[i];
      mtT += 0.5 * anaT.T[i];
      mtB += 0.5 * anaB.T[i];
    }

    HM.Fill1d("mtDiff", 0, mtB - mtT);
    HM.Fill1d("mtTop", 0, mtT);
    HM.Fill1d("mtBtm", 0, mtB);

    double maxchi2 = max( max(anaB.X2[1], anaB.X2[0]), max(anaT.X2[1], anaT.X2[0]) );

    if (TMath::Abs(mtT-mtB) > 7) {continue;}
    if (maxchi2 > 20) {continue;}

    double thetaR = TMath::ATan( (xT - xB) / (zT - zB) );
    double thetaXY = TMath::ATan( (xT - xB) / crtModulesDistance ); 
    double thetaZY = TMath::ATan( (zT - zB) / crtModulesDistance );

    if ( iSc_out[0] == iSc_out[1] )  {HM.Fill1d("thetaZYvert", 0, thetaZY);}

    HM.Fill1d("thetaZY", 0, thetaZY);
    HM.Fill1d("thetaXY", 0, thetaXY);
    HM.Fill1d("thetaR", 0, thetaR);

    HM.Fill2d("heatMap", 0, zT, anaT.iSc);
    HM.Fill2d("heatMap", 0, zB, anaB.iSc + 1 + scintNum);

    HM.Fill1d("zetaPerScint", (double)iSc_out[0], zT);
    HM.Fill1d("zetaPerScint", (double)iSc_out[1] + scintNum, zB);

    HM.Fill2d("zetaTop_Btm", 0, zT, zB);
    HM.Fill2d("xTop_Btm", 0, (double)iSc_out[0], (double)iSc_out[1]);

    double tt = TMath::Abs(thetaXY);
    int idx = (tt > 0.03) + (tt > 0.07) + (tt > 0.12); 

    HM.Fill1d("zetaTopPerThetaXY", idx, zT);

    





    






    CRTs4->Fill();

  }
  
  HM.ProcessBoxes();
  CRTs4->Write();
  outFile->Close();

  cout<<endl<<endl<<"                         root "<<outFileName<<endl;
}





