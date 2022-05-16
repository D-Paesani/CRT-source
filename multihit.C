// ad esempio
// CRT->Draw("Z:iSc", "Z!=0 && iTrig < 33e3 && nHitBar == 5", "L")
// per vedere il tempo servirebbe di convertire i th3f a graph2d, o farlo da questo codice qui

#include <fstream>
#include <chrono>
#include <TLine.h>
#include <iostream>
#include <list>

#include "TApplication.h"
#include "TSpectrum.h"
#include "TGraphErrors.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TH1F.h"
#include "TF1.h"
#include <TF1Convolution.h>

#include "includes/Analysis.h"
#include "includes/langaus.h"
#include "includes/AnaPars.h"
#include "includes/CsvHandler.h"
#include "includes/MipSelection.h"
#include "includes/templ2charge.h"
#include "includes/HistManager.h"
#include "includes/NumberingHelper.h"

#include <TSpline.h>

using namespace std;

CsvHandler CSV;
MipSelection Selection;
HistManager HM;

Long64_t nentries, nbytes, nb, ientry, jentry;

double **timeOffset, **zetaOffset;

double_t Z_out[8] = {0}, T_out[8] = {0};
int iSc_out[8] = {0, 1, 2, 3, 4, 5, 6, 7};
Long64_t jTrig_out, nhb_out;

TTree *CRTs3;
#define tree_out_name "CRT"
#define tree_index_out "iTrig"
#define reindex_tree 0


void Analysis::LoopOverEntries() {

  nentries = fChain->GetEntriesFast();
  Long64_t etp = min(maxEvToProcess3, nentries);
  cout << "Number of events to process: " << etp << endl << endl;
  nbytes = 0, nb = 0;

  jTrig_out = 0;

  for (jentry = 0; jentry < etp; jentry++) {

    ientry = LoadTree(jentry);

    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry); 
    nbytes += nb;

    if (!(jentry%10000)) {cout << Form( "     processing evt %lld / %lld  ( %.0f%% )", jentry, etp, (float)(100*jentry/etp) ) << endl;}

    double teT[16] = {0};

    for(int hit = 0; hit < nCry; hit++){

      if ( (Vmax[hit] > 1800) || (Qval[hit] < 80) || templChi2[hit] > 10 ) continue;

      if ( (modulSel=="B" && iMod[hit]==0) || (modulSel=="T" && iMod[hit]==1) ) {continue;} // TOP = 0, BTM = 1

      int hitSide=iSide[hit], hitScint = iScint[hit], hitN = hitSide*scintNum + hitScint;

      teT[hitN] = templTime[hit] - timeOffset[hitSide][hitScint];
    }

    nhb_out = 0;
    TGraphErrors g;
    TCanvas c("c", "c");
    for(int iSc=0; iSc<8; iSc++){
      if ( (teT[iSc] != 0)  && (teT[iSc+8] != 0) ){
        Z_out[iSc] = teT[iSc] - teT[iSc+8] - zetaOffset[0][iSc];
        g.AddPoint(iSc*2.5, Z_out[iSc]);
        g.SetPointError(g.GetN()-1, 0.8, 1.6);
        T_out[iSc] = (teT[iSc] + teT[iSc+8])/2;
        nhb_out++;
      }
      else Z_out[iSc] = 0;
    }

    if (nhb_out > 2){
      c.cd();
      g.Draw();
      c.SaveAs(Form("../../data/step3/multihit/%lli_%lli.png", nhb_out, jTrig_out));
      CRTs3->Fill();
      jTrig_out++;
    }

  }

  cout<<endl;
}


void Analysis::Loop(){

  cout<<endl<<endl<<"::::::::::::::::::::: CRT analysis step 3 :::::::::::::::::::::"<<endl<<endl;

  if (fChain == 0) return;

  cout<<"Creating histograms :"<<endl;

  cout<<"...done"<<endl<<endl;

  cout<<"Loading cutGs :";
  Selection.loadCutG();
  cout<<"...done"<<endl<<endl;

  cout<<"Retrieving calibration data from [" + lutPrefix3p + "] :"<<endl;
  timeOffset  = CSV.InitMatrix(2, scintNum);
  zetaOffset  = CSV.InitMatrix(1, scintNum);

  CSV.Read(CSV.GetFirstFile(lutPrefix3p + calName + lutTimeOffsName + "*"),  ',', timeOffset,  2, scintNum);
  CSV.Read(CSV.GetFirstFile(lutPrefix3p + calName + lutZetaOffName + "*"),   ',', zetaOffset,  1, scintNum);

  cout<<"...done"<<endl<<endl;

  outFile->cd();
  CRTs3 = new TTree(tree_out_name, tree_out_name);
  CRTs3->SetAutoSave(1000);
  CRTs3->Branch("iTrig",   &jTrig_out,  "iTrig/L");
  CRTs3->Branch("Z",       &Z_out,      "Z[8]/D");
  CRTs3->Branch("T",       &T_out,      "T[8]/D");
  CRTs3->Branch("iSc",       &iSc_out,      "iSc[8]/I");
  CRTs3->Branch("nHitBar",   &nhb_out,    "nHitBar/I");

  Analysis::LoopOverEntries();

  cout<<"...done"<<endl;

  outFile->cd();
  CRTs3->Write();

  cout<<endl<<endl<<"::::::::::::::::::::: analysis done :::::::::::::::::::::"<<endl<<endl;
}



#define inFile_f "../../data/step2/%s_s2.root"
#define outFile_f "../../data/step3/%s_multihit.root"

void multihit(TString run_name, TString calib_name, TString mod_select = "") {

  gErrorIgnoreLevel = kFatal; //kPrint, kInfo, kWarning, kError, kBreak, kSysError, kFatal

  if (mod_select != "T" && mod_select != "B" && mod_select != "") {cout<<"modSel must be empty, 'T' or 'B' !!"<<endl; return;}
  //opzione "" per legacy

  TString inFileName = Form(inFile_f, run_name.Data());
  TString runName = mod_select == "" ? run_name : run_name + "_" + mod_select;
  TString calibName = mod_select == "" ? calib_name : calib_name + "_" + mod_select;

  TString outFileName = Form(outFile_f, runName.Data());

  TFile *fileOut = new TFile(outFileName, "RECREATE");

  cout<<endl<<"------------> Launching step3:"<<endl;
  cout<<"----> runName : "<<runName<<endl;
  cout<<"----> Input file: "<<inFileName<<endl;
  cout<<"----> Output file: "<<outFileName<<endl;
  cout<<"----> Calibration files: "<<calibName<<endl;
  cout<<"----> Module selector: "<<mod_select<<endl<<endl;


  Analysis *a = new Analysis(inFileName, fileOut, runName, calibName, mod_select);

  a->Loop();

}
