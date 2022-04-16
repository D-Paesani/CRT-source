//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Apr 14 17:56:33 2022 by ROOT version 6.24/02
// from TTree CRT/CRT
// found on file: run233_s2.root
//////////////////////////////////////////////////////////

#ifndef two_modules_patch_h
#define two_modules_patch_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

#define maxNcry 16
#define maxNsample 1024
// Header file for the classes stored in the TTree if any.

class two_modules_patch {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           ntrig;
   Int_t           evnum;
   Int_t           nsample;
   Double_t        time[1024];   //[nsample]
   Int_t           nCry;
   Int_t           iDAQ[16];   //[nCry]
   Int_t           iScint[16];   //[nCry]
   Int_t           iSide[16];   //[nCry]
   Int_t           iMax[16];   //[nCry]
   Double_t        Vmax[16];   //[nCry]
   Double_t        Qval[16];   //[nCry]
   Double_t        Tval[16];   //[nCry]
   Double_t        pedL[16];   //[nCry]
   Double_t        pedH[16];   //[nCry]
   Double_t        wave[16][1024];   //[nCry]
   Double_t        bline[16];   //[nCry]
   Double_t        templTime[16];   //[nCry]

   // List of branches
   TBranch        *b_ntrig;   //!
   TBranch        *b_evnum;   //!
   TBranch        *b_nsample;   //!
   TBranch        *b_time;   //!
   TBranch        *b_nCry;   //!
   TBranch        *b_iDAQ;   //!
   TBranch        *b_iScint;   //!
   TBranch        *b_iSide;   //!
   TBranch        *b_iMax;   //!
   TBranch        *b_Vmax;   //!
   TBranch        *b_Qval;   //!
   TBranch        *b_Tval;   //!
   TBranch        *b_pedL;   //!
   TBranch        *b_pedH;   //!
   TBranch        *b_wave;   //!
   TBranch        *b_bline;   //!
   TBranch        *b_lognTime;   //!
   TBranch        *b_lognChi2;   //!
   TBranch        *b_lognFit;   //!
   TBranch        *b_lognErr;   //!
   TBranch        *b_templTime;   //!
   TBranch        *b_templChi2;   //!
   TBranch        *b_templFit;   //!
   TBranch        *b_templErr;   //!


   Int_t           ntrig_out;
   Int_t           evnum_out;
   Int_t           nsample_out;
   Double_t        time_out[maxNsample];   //[num]
   Int_t           nCry_out;
   Int_t           iDAQ_out[maxNcry];
   Int_t           iScint_out[maxNcry];
   Int_t           iSide_out[maxNcry];
   Int_t           iMod_out[maxNcry];
   Int_t           iMax_out[maxNcry];
   Double_t        Vmax_out[maxNcry];
   Double_t        Qval_out[maxNcry];
   Double_t        Tval_out[maxNcry];
   Double_t        pedL_out[maxNcry];
   Double_t        pedH_out[maxNcry];
   Double_t        wave_out[maxNcry][maxNsample];
   Double_t        bline_out[maxNcry];
   Double_t        templTime_out[maxNcry];

   two_modules_patch(TTree *tree=0);
   virtual ~two_modules_patch();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef two_modules_patch_cxx
two_modules_patch::two_modules_patch(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("../data/step2/run233_orig_s2.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("../data/step2/run233_orig_s2.root");
      }
      f->GetObject("CRT",tree);

   }
   Init(tree);
}

two_modules_patch::~two_modules_patch()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t two_modules_patch::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t two_modules_patch::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void two_modules_patch::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("ntrig", &ntrig, &b_ntrig);
   fChain->SetBranchAddress("evnum", &evnum, &b_evnum);
   fChain->SetBranchAddress("nsample", &nsample, &b_nsample);
   fChain->SetBranchAddress("time", time, &b_time);
   fChain->SetBranchAddress("nCry", &nCry, &b_nCry);
   fChain->SetBranchAddress("iDAQ", iDAQ, &b_iDAQ);
   fChain->SetBranchAddress("iScint", iScint, &b_iScint);
   fChain->SetBranchAddress("iSide", iSide, &b_iSide);
   fChain->SetBranchAddress("iMax", iMax, &b_iMax);
   fChain->SetBranchAddress("Vmax", Vmax, &b_Vmax);
   fChain->SetBranchAddress("Qval", Qval, &b_Qval);
   fChain->SetBranchAddress("Tval", Tval, &b_Tval);
   fChain->SetBranchAddress("pedL", pedL, &b_pedL);
   fChain->SetBranchAddress("pedH", pedH, &b_pedH);
   fChain->SetBranchAddress("wave", wave, &b_wave);
   fChain->SetBranchAddress("bline", bline, &b_bline);
   fChain->SetBranchAddress("templTime", templTime, &b_templTime);
   Notify();
}

Bool_t two_modules_patch::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void two_modules_patch::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t two_modules_patch::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef two_modules_patch_cxx
