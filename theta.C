#include <iostream>
#include <list>
#include "includes/step3out.h"
#include <TH1F.h>
#include <TMath.h>

using namespace std;

const Long64_t maxEvToProcess = 1e6;
const int debug = 0;

void theta(){
  TString filename[2] = {
    "../data/step3/run233_T_s3.root", //master
    "../data/step3/run233_B_s3.root"  //slave
  };

  step3out *inst[2];

  auto *h = new TH1F("sdiff", "sdiff", 15, -7.5, 7.5); //sarebbe tangent theta * DY / 1.5

  for(int i=0; i<2; i++) inst[i] = new step3out(filename[i]);

  inst[1]->fChain->BuildIndex("iTrig");

  Long64_t nentries = inst[0]->fChain->GetEntriesFast();
  Long64_t etp = min(maxEvToProcess, nentries);
  cout << "Number of events to process: " << etp << endl << endl;


  Long64_t nbytes = 0, nb = 0;

  for (Long64_t jentry = 0; jentry < etp; jentry++) {

    Long64_t ientry = inst[0]->LoadTree(jentry);

    if (ientry < 0) break;
    nb = inst[0]->fChain->GetEntry(jentry);
    nbytes += nb;

    if (!(jentry%10000)) {cout << Form( "     processing evt %lld / %lld  ( %.0f%% )", jentry, etp, (float)(100*jentry/etp) ) << endl;}

    if(debug) cout << inst[0]->iTrig << "\t";
    Long64_t slave_entry = inst[1]->fChain->GetEntryNumberWithIndex(inst[0]->iTrig);

    if (slave_entry < 0) {
     if(debug) cout << endl;
     continue;
    }
    if(inst[1]->LoadTree(slave_entry) < 0) break;
    inst[1]->fChain->GetEntry(slave_entry);

    if(debug) cout << inst[1]->iTrig << endl;

    h->Fill((inst[0]->iSc - inst[1]->iSc));
  }
  h->Draw("e");
}
