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

#include "TFile.h"
#include "TTree.h"

#include "../../CRT_root/HistManager.h"

using namespace std;

void HistManager::CreateHistDict(vector <Double_t> pars){

  hist_dict = {};
  for(int i=0; i<20; i++){
    auto pair = CreatePair(Form("QmipZslices_%i", i),  1, 8, 2, "Qmip per Z slices", "Side %i - Scint %i", "Q" ,"pC", 500, 50, 1050);
    hist_dict.insert(pair);
  }

}

void find_valleys(){
  TFile *f = new TFile("../data/step3/run199_s3.root");
  TTree *tree = (TTree*)f->Get("CRT");

  TFile *outf = new TFile("../data/step3/run199_s3_valleys.root", "recreate");
  HistManager *HM = new HistManager(outf, -1);

  HM->CreateHistDict({});


  for(int iSl = 0; iSl < 20; iSl++){
    for(int iSc = 0; iSc < 8; iSc ++){
      for(int iSd = 0; iSd < 2; iSd++){
        TString histName = Form("QmipZslices_%i_%i_%i", iSl, iSd, iSc);
        TH1F *temp = HM->GetHist(Form("QmipZslices_%i", iSl), iSc, iSd);
        tree->Draw(Form("crt_Q[%i]>>%s", iSd, histName.Data()), Form("crt_Z > -80+%i*8 && crt_Z< -80+%i*8+8 && crt_iSc == %i", iSl, iSl, iSc));
      }
    }
  }

  auto hist_dict = HM->hist_dict;
  for (auto& hist_matrix: hist_dict) {
    hist_matrix.second->draw_all();
  }
  HM->CloseOutFile();

}