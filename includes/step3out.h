//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sat Apr 16 21:15:58 2022 by ROOT version 6.24/02
// from TTree CRT/CRT
// found on file: ../data/step3/run233_T_s3.root
// AND THEN CONVENIENTLY MODIFIED
//////////////////////////////////////////////////////////

#ifndef step3out_h
#define step3out_h

#include <TString.h>
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class step3out {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Long64_t        iTrig;
   Int_t           iSc;
   Double_t        Z;
   Double_t        pZ;
   Double_t        Q[2];
   Double_t        T[2];
   Double_t        X2[2];

   // List of branches
   TBranch        *b_iTrig;   //!
   TBranch        *b_iSc;   //!
   TBranch        *b_Z;   //!
   TBranch        *b_pZ;   //!
   TBranch        *b_Q;   //!
   TBranch        *b_T;   //!
   TBranch        *b_X2;   //!

   step3out(TString filename);
   virtual ~step3out();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

//#ifdef step3out_cxx
step3out::step3out(TString filename) : fChain(0)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   TFile *f = new TFile(filename);
   TTree *tree = (TTree*) f->Get("CRT");

   Init(tree);
}

step3out::~step3out()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t step3out::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}

Long64_t step3out::LoadTree(Long64_t entry)
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

void step3out::Init(TTree *tree)
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

   fChain->SetBranchAddress("iTrig", &iTrig, &b_iTrig);
   fChain->SetBranchAddress("iSc", &iSc, &b_iSc);
   fChain->SetBranchAddress("Z", &Z, &b_Z);
   fChain->SetBranchAddress("pZ", &pZ, &b_pZ);
   fChain->SetBranchAddress("Q", Q, &b_Q);
   fChain->SetBranchAddress("T", T, &b_T);
   fChain->SetBranchAddress("X2", X2, &b_X2);
   Notify();
}

Bool_t step3out::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void step3out::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t step3out::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
//#endif // #ifdef step3out_cxx
