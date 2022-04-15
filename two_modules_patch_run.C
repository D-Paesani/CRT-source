#define two_modules_patch_cxx
#include "two_modules_patch.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void two_modules_patch::Loop()
{
//   In a ROOT session, you can do:
//      root> .L two_modules_patch.C
//      root> two_modules_patch t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch

   TFile *f = new TFile("run233_new_s2.root", "recreate");

   TTree *out = new TTree("CRT", "CRT");
   out->SetAutoSave(100000);

   out->Branch("ntrig", &ntrig_out, "ntrig/I");
   out->Branch("evnum", &evnum_out, "evnum/I");
   out->Branch("nsample", &nsample_out, "nsample/I");     // function of nCry?
   out->Branch("time", &time_out, "time[nsample]/D"); // function of nCry?
   out->Branch("nCry", &nCry_out, "nCry/I");
   out->Branch("iDAQ", &iDAQ_out, "iDAQ[nCry]/I");
   out->Branch("iScint", &iScint_out, "iScint[nCry]/I");
   out->Branch("iSide", &iSide_out, "iSide[nCry]/I");
   out->Branch("iMod", &iMod_out, "iMod[nCry]/I");
   out->Branch("iMax", &iMax_out, "iMax[nCry]/I");
   out->Branch("Vmax", &Vmax_out, "Vmax[nCry]/D");
   out->Branch("Qval", &Qval_out, "Qval[nCry]/D");
   out->Branch("Tval", &Tval_out, "Tval[nCry]/D");
   out->Branch("pedL", &pedL_out, "pedL[nCry]/D");
   out->Branch("pedH", &pedH_out, "pedH[nCry]/D");
   out->Branch("wave", &wave_out, "wave[nCry][1024]/D");
   out->Branch("bline", &bline_out, "bline[nCry]/D");
   out->Branch("templTime", &templTime_out, "templTime[nCry]/D");

   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      if (jentry % 10000 == 0) cout << jentry << endl;
      ntrig_out = ntrig;

      evnum_out = evnum;
      nsample_out = nsample;
      std::copy(time, time+nsample, time_out); 
      nCry_out = nCry;

      for(int hit = 0; hit < nCry; hit++){
         int ch = iDAQ[hit];

         iDAQ_out[hit] = iDAQ[hit];
         iScint_out[hit] = (ch%4)*2;
         iSide_out[hit] = (ch/4)%2;
         iMod_out[hit] = ch/8;
         iMax_out[hit] = iMax[hit];
         Vmax_out[hit] = Vmax[hit];
         Qval_out[hit] = Qval[hit];
         Tval_out[hit] = Tval[hit];
         pedL_out[hit] = pedL[hit];
         pedH_out[hit] = pedH[hit];
         std::copy(wave[hit], wave[hit]+nsample, wave_out[hit]);
         bline_out[hit] = bline[hit];
         templTime_out[hit] = templTime[hit];
      }

      out->Fill();
   }
   out->Write();
}

void two_modules_patch_run(){
   two_modules_patch *a = new two_modules_patch();
   a->Loop();

}
