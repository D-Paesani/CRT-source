//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Nov 30 10:37:50 2020 by ROOT version 5.34/14
// from TTree Wavefull/Wavefull
// found on file: ../roottople/provamod0.root
//////////////////////////////////////////////////////////

#ifndef analysis_CRT_h
#define analysis_CRT_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TString.h>
#include <TGraphErrors.h>
#include <iostream>

#include <string>
#include <fstream>
#include <iomanip>
#include <stdlib.h>
#include <sstream>

#define NCHAN 16

const int maxNcry    =   51;
const int maxNsample = 1024;
float           Qcut = 2;

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class analysis_CRT {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain
   // Declaration of leaf types
   Int_t           ntrig;
   Int_t           evnum;
   Int_t           num;
   Double_t        time[1024];   //[num]
   Double_t        onda[NCHAN][1024];   //[num]
   Double_t        ped[NCHAN];
   Double_t        pedh[NCHAN];
   Double_t        charge[NCHAN];
   Double_t        chargepeak[NCHAN];
   Int_t           imax[NCHAN];
   Double_t        tave[NCHAN];

   // List of branches
   TBranch        *b_ntrig;   //!
   TBranch        *b_evnum;   //!
   TBranch        *b_num;   //!
   TBranch        *b_time;   //!
   TBranch        *b_onda[NCHAN];   //!
   TBranch        *b_ped[NCHAN];   //!
   TBranch        *b_pedh[NCHAN];   //!
   TBranch        *b_charge[NCHAN];   //!
   TBranch        *b_chargepeak[NCHAN];   //!
   TBranch        *b_imax[NCHAN];   //!
   TBranch        *b_tave[NCHAN];   //!

   //analysis_CRT(TTree *tree=0);
   analysis_CRT(TString fName="",TTree *tree=0);
   virtual ~analysis_CRT();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   //virtual void     Loop();
   virtual void     Loop(TString OutputFile="", int evflag=0);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
   // Additional user functions
   virtual void  BookOutput(TTree *CRT=0);
   virtual void  GetValues (int Ichan);

   //
   // Output nutple
   //

   TString         OutputFile;

   Int_t           ntrig_out;
   Int_t           evnum_out;
   Int_t           nsample;
   Double_t        time_out[maxNsample];   //[num]
   Int_t           nCry;
   Int_t           iDAQ[maxNcry];
   Int_t           iScint[maxNcry];
   Int_t           iSide[maxNcry];
   Int_t           iMod[maxNcry];
   Int_t           iMax[maxNcry];
   Double_t        Vmax[maxNcry];
   Double_t        Qval[maxNcry];
   Double_t        Tval[maxNcry];
   Double_t        pedL[maxNcry];
   Double_t        pedH[maxNcry];
   Double_t        wave[maxNcry][maxNsample];
   Double_t        bline[maxNcry];
};

#endif

#ifdef analysis_CRT_cxx
//analysis_CRT::analysis_CRT(TTree *tree) : fChain(0) 
analysis_CRT::analysis_CRT(TString fName,TTree *tree)
{

// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
/*   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("../roottople/provamod0.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("../roottople/provamod0.root");
      }
      f->GetObject("Wavefull",tree);

   }
*/

   if (!fName.IsNull()) {
    TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject(fName.Data());
    if (!f || !f->IsOpen()) {
      f = new TFile(fName.Data());
    }
    f->GetObject("Wavefull",tree);
  }
  Init(tree);

}

analysis_CRT::~analysis_CRT()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t analysis_CRT::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}

Long64_t analysis_CRT::LoadTree(Long64_t entry)
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

void analysis_CRT::Init(TTree *tree)
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
   fChain->SetBranchAddress("num", &num, &b_num);
   fChain->SetBranchAddress("time", time, &b_time);
   for(int iCh=0; iCh<NCHAN; iCh++){
      int nboard = iCh/8;
      int board_ch = iCh%8;
      fChain->SetBranchAddress(Form("b%i_onda%i", nboard, board_ch),       &(onda[iCh]),       &(b_onda[iCh]));
      fChain->SetBranchAddress(Form("b%i_ped%i", nboard, board_ch),        &(ped[iCh]),        &(b_ped[iCh]));
      fChain->SetBranchAddress(Form("b%i_pedh%i", nboard, board_ch),       &(pedh[iCh]),       &(b_pedh[iCh]));
      fChain->SetBranchAddress(Form("b%i_charge%i", nboard, board_ch),     &(charge[iCh]),     &(b_charge[iCh]));
      fChain->SetBranchAddress(Form("b%i_chargepeak%i", nboard, board_ch), &(chargepeak[iCh]), &(b_chargepeak[iCh]));
      fChain->SetBranchAddress(Form("b%i_imax%i", nboard, board_ch),       &(imax[iCh]),       &(b_imax[iCh]));
      fChain->SetBranchAddress(Form("b%i_tave%i", nboard, board_ch),       &(tave[iCh]),       &(b_tave[iCh]));
   }
   Notify();

}

Bool_t analysis_CRT::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void analysis_CRT::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t analysis_CRT::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef analysis_CRT_cxx
