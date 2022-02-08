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
   Double_t        b0_onda0[1024];   //[num]
   Double_t        b0_ped0;
   Double_t        b0_pedh0;
   Double_t        b0_charge0;
   Double_t        b0_chargepeak0;
   Int_t           b0_imax0;
   Double_t        b0_tave0;
   Double_t        b0_onda1[1024];   //[num]
   Double_t        b0_ped1;
   Double_t        b0_pedh1;
   Double_t        b0_charge1;
   Double_t        b0_chargepeak1;
   Int_t           b0_imax1;
   Double_t        b0_tave1;
   Double_t        b0_onda2[1024];   //[num]
   Double_t        b0_ped2;
   Double_t        b0_pedh2;
   Double_t        b0_charge2;
   Double_t        b0_chargepeak2;
   Int_t           b0_imax2;
   Double_t        b0_tave2;
   Double_t        b0_onda3[1024];   //[num]
   Double_t        b0_ped3;
   Double_t        b0_pedh3;
   Double_t        b0_charge3;
   Double_t        b0_chargepeak3;
   Int_t           b0_imax3;
   Double_t        b0_tave3;
   Double_t        b0_onda4[1024];   //[num]
   Double_t        b0_ped4;
   Double_t        b0_pedh4;
   Double_t        b0_charge4;
   Double_t        b0_chargepeak4;
   Int_t           b0_imax4;
   Double_t        b0_tave4;
   Double_t        b0_onda5[1024];   //[num]
   Double_t        b0_ped5;
   Double_t        b0_pedh5;
   Double_t        b0_charge5;
   Double_t        b0_chargepeak5;
   Int_t           b0_imax5;
   Double_t        b0_tave5;
   Double_t        b0_onda6[1024];   //[num]
   Double_t        b0_ped6;
   Double_t        b0_pedh6;
   Double_t        b0_charge6;
   Double_t        b0_chargepeak6;
   Int_t           b0_imax6;
   Double_t        b0_tave6;
   Double_t        b0_onda7[1024];   //[num]
   Double_t        b0_ped7;
   Double_t        b0_pedh7;
   Double_t        b0_charge7;
   Double_t        b0_chargepeak7;
   Int_t           b0_imax7;
   Double_t        b0_tave7;
   Double_t        b1_onda0[1024];   //[num]
   Double_t        b1_ped0;
   Double_t        b1_pedh0;
   Double_t        b1_charge0;
   Double_t        b1_chargepeak0;
   Int_t           b1_imax0;
   Double_t        b1_tave0;
   Double_t        b1_onda1[1024];   //[num]
   Double_t        b1_ped1;
   Double_t        b1_pedh1;
   Double_t        b1_charge1;
   Double_t        b1_chargepeak1;
   Int_t           b1_imax1;
   Double_t        b1_tave1;
   Double_t        b1_onda2[1024];   //[num]
   Double_t        b1_ped2;
   Double_t        b1_pedh2;
   Double_t        b1_charge2;
   Double_t        b1_chargepeak2;
   Int_t           b1_imax2;
   Double_t        b1_tave2;
   Double_t        b1_onda3[1024];   //[num]
   Double_t        b1_ped3;
   Double_t        b1_pedh3;
   Double_t        b1_charge3;
   Double_t        b1_chargepeak3;
   Int_t           b1_imax3;
   Double_t        b1_tave3;
   Double_t        b1_onda4[1024];   //[num]
   Double_t        b1_ped4;
   Double_t        b1_pedh4;
   Double_t        b1_charge4;
   Double_t        b1_chargepeak4;
   Int_t           b1_imax4;
   Double_t        b1_tave4;
   Double_t        b1_onda5[1024];   //[num]
   Double_t        b1_ped5;
   Double_t        b1_pedh5;
   Double_t        b1_charge5;
   Double_t        b1_chargepeak5;
   Int_t           b1_imax5;
   Double_t        b1_tave5;
   Double_t        b1_onda6[1024];   //[num]
   Double_t        b1_ped6;
   Double_t        b1_pedh6;
   Double_t        b1_charge6;
   Double_t        b1_chargepeak6;
   Int_t           b1_imax6;
   Double_t        b1_tave6;
   Double_t        b1_onda7[1024];   //[num]
   Double_t        b1_ped7;
   Double_t        b1_pedh7;
   Double_t        b1_charge7;
   Double_t        b1_chargepeak7;
   Int_t           b1_imax7;
   Double_t        b1_tave7;

   // List of branches
   TBranch        *b_ntrig;   //!
   TBranch        *b_evnum;   //!
   TBranch        *b_num;   //!
   TBranch        *b_time;   //!
   TBranch        *b_b0_onda0;   //!
   TBranch        *b_b0_ped0;   //!
   TBranch        *b_b0_pedhigh0;   //!
   TBranch        *b_b0_charge0;   //!
   TBranch        *b_b0_chargepeak0;   //!
   TBranch        *b_b0_imax0;   //!
   TBranch        *b_b0_tave0;   //!
   TBranch        *b_b0_onda1;   //!
   TBranch        *b_b0_ped1;   //!
   TBranch        *b_b0_pedhigh1;   //!
   TBranch        *b_b0_charge1;   //!
   TBranch        *b_b0_chargepeak1;   //!
   TBranch        *b_b0_imax1;   //!
   TBranch        *b_b0_tave1;   //!
   TBranch        *b_b0_onda2;   //!
   TBranch        *b_b0_ped2;   //!
   TBranch        *b_b0_pedhigh2;   //!
   TBranch        *b_b0_charge2;   //!
   TBranch        *b_b0_chargepeak2;   //!
   TBranch        *b_b0_imax2;   //!
   TBranch        *b_b0_tave2;   //!
   TBranch        *b_b0_onda3;   //!
   TBranch        *b_b0_ped3;   //!
   TBranch        *b_b0_pedhigh3;   //!
   TBranch        *b_b0_charge3;   //!
   TBranch        *b_b0_chargepeak3;   //!
   TBranch        *b_b0_imax3;   //!
   TBranch        *b_b0_tave3;   //!
   TBranch        *b_b0_onda4;   //!
   TBranch        *b_b0_ped4;   //!
   TBranch        *b_b0_pedhigh4;   //!
   TBranch        *b_b0_charge4;   //!
   TBranch        *b_b0_chargepeak4;   //!
   TBranch        *b_b0_imax4;   //!
   TBranch        *b_b0_tave4;   //!
   TBranch        *b_b0_onda5;   //!
   TBranch        *b_b0_ped5;   //!
   TBranch        *b_b0_pedhigh5;   //!
   TBranch        *b_b0_charge5;   //!
   TBranch        *b_b0_chargepeak5;   //!
   TBranch        *b_b0_imax5;   //!
   TBranch        *b_b0_tave5;   //!
   TBranch        *b_b0_onda6;   //!
   TBranch        *b_b0_ped6;   //!
   TBranch        *b_b0_pedhigh6;   //!
   TBranch        *b_b0_charge6;   //!
   TBranch        *b_b0_chargepeak6;   //!
   TBranch        *b_b0_imax6;   //!
   TBranch        *b_b0_tave6;   //!
   TBranch        *b_b0_onda7;   //!
   TBranch        *b_b0_ped7;   //!
   TBranch        *b_b0_pedhigh7;   //!
   TBranch        *b_b0_charge7;   //!
   TBranch        *b_b0_chargepeak7;   //!
   TBranch        *b_b0_imax7;   //!
   TBranch        *b_b0_tave7;   //!
   TBranch        *b_b1_onda0;   //!
   TBranch        *b_b1_ped0;   //!
   TBranch        *b_b1_pedhigh0;   //!
   TBranch        *b_b1_charge0;   //!
   TBranch        *b_b1_chargepeak0;   //!
   TBranch        *b_b1_imax0;   //!
   TBranch        *b_b1_tave0;   //!
   TBranch        *b_b1_onda1;   //!
   TBranch        *b_b1_ped1;   //!
   TBranch        *b_b1_pedhigh1;   //!
   TBranch        *b_b1_charge1;   //!
   TBranch        *b_b1_chargepeak1;   //!
   TBranch        *b_b1_imax1;   //!
   TBranch        *b_b1_tave1;   //!
   TBranch        *b_b1_onda2;   //!
   TBranch        *b_b1_ped2;   //!
   TBranch        *b_b1_pedhigh2;   //!
   TBranch        *b_b1_charge2;   //!
   TBranch        *b_b1_chargepeak2;   //!
   TBranch        *b_b1_imax2;   //!
   TBranch        *b_b1_tave2;   //!
   TBranch        *b_b1_onda3;   //!
   TBranch        *b_b1_ped3;   //!
   TBranch        *b_b1_pedhigh3;   //!
   TBranch        *b_b1_charge3;   //!
   TBranch        *b_b1_chargepeak3;   //!
   TBranch        *b_b1_imax3;   //!
   TBranch        *b_b1_tave3;   //!
   TBranch        *b_b1_onda4;   //!
   TBranch        *b_b1_ped4;   //!
   TBranch        *b_b1_pedhigh4;   //!
   TBranch        *b_b1_charge4;   //!
   TBranch        *b_b1_chargepeak4;   //!
   TBranch        *b_b1_imax4;   //!
   TBranch        *b_b1_tave4;   //!
   TBranch        *b_b1_onda5;   //!
   TBranch        *b_b1_ped5;   //!
   TBranch        *b_b1_pedhigh5;   //!
   TBranch        *b_b1_charge5;   //!
   TBranch        *b_b1_chargepeak5;   //!
   TBranch        *b_b1_imax5;   //!
   TBranch        *b_b1_tave5;   //!
   TBranch        *b_b1_onda6;   //!
   TBranch        *b_b1_ped6;   //!
   TBranch        *b_b1_pedhigh6;   //!
   TBranch        *b_b1_charge6;   //!
   TBranch        *b_b1_chargepeak6;   //!
   TBranch        *b_b1_imax6;   //!
   TBranch        *b_b1_tave6;   //!
   TBranch        *b_b1_onda7;   //!
   TBranch        *b_b1_ped7;   //!
   TBranch        *b_b1_pedhigh7;   //!
   TBranch        *b_b1_charge7;   //!
   TBranch        *b_b1_chargepeak7;   //!
   TBranch        *b_b1_imax7;   //!
   TBranch        *b_b1_tave7;   //!

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
   virtual void  GetValues (int Iboard, int Ichan, int evflag=0);
   virtual void  getTlogn  ();
   virtual void  getTemplateFit(TGraphErrors* gt);
   virtual void  getTemplate2Fit();

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
   Int_t           iMax[maxNcry];
   Double_t        Vmax[maxNcry];
   Double_t        Qval[maxNcry];
   Double_t        Tval[maxNcry];
   Double_t        pedL[maxNcry];
   Double_t        pedH[maxNcry];
   Double_t        wave[maxNcry][maxNsample];
   Double_t        bline[maxNcry];
   Double_t        lognTime[maxNcry];
   Double_t        lognChi2[maxNcry];
   Double_t        lognFit[maxNcry][4];
   Double_t        lognErr[maxNcry][4];
   Double_t        templTime[maxNcry];
   Double_t        templChi2[maxNcry];
   Double_t        templFit[maxNcry][3];
   Double_t        templErr[maxNcry][3];
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
   fChain->SetBranchAddress("b0_onda0", b0_onda0, &b_b0_onda0);
   fChain->SetBranchAddress("b0_ped0", &b0_ped0, &b_b0_ped0);
   fChain->SetBranchAddress("b0_pedh0", &b0_pedh0, &b_b0_pedhigh0);
   fChain->SetBranchAddress("b0_charge0", &b0_charge0, &b_b0_charge0);
   fChain->SetBranchAddress("b0_chargepeak0", &b0_chargepeak0, &b_b0_chargepeak0);
   fChain->SetBranchAddress("b0_imax0", &b0_imax0, &b_b0_imax0);
   fChain->SetBranchAddress("b0_tave0", &b0_tave0, &b_b0_tave0);
   fChain->SetBranchAddress("b0_onda1", b0_onda1, &b_b0_onda1);
   fChain->SetBranchAddress("b0_ped1", &b0_ped1, &b_b0_ped1);
   fChain->SetBranchAddress("b0_pedh1", &b0_pedh1, &b_b0_pedhigh1);
   fChain->SetBranchAddress("b0_charge1", &b0_charge1, &b_b0_charge1);
   fChain->SetBranchAddress("b0_chargepeak1", &b0_chargepeak1, &b_b0_chargepeak1);
   fChain->SetBranchAddress("b0_imax1", &b0_imax1, &b_b0_imax1);
   fChain->SetBranchAddress("b0_tave1", &b0_tave1, &b_b0_tave1);
   fChain->SetBranchAddress("b0_onda2", b0_onda2, &b_b0_onda2);
   fChain->SetBranchAddress("b0_ped2", &b0_ped2, &b_b0_ped2);
   fChain->SetBranchAddress("b0_pedh2", &b0_pedh2, &b_b0_pedhigh2);
   fChain->SetBranchAddress("b0_charge2", &b0_charge2, &b_b0_charge2);
   fChain->SetBranchAddress("b0_chargepeak2", &b0_chargepeak2, &b_b0_chargepeak2);
   fChain->SetBranchAddress("b0_imax2", &b0_imax2, &b_b0_imax2);
   fChain->SetBranchAddress("b0_tave2", &b0_tave2, &b_b0_tave2);
   fChain->SetBranchAddress("b0_onda3", b0_onda3, &b_b0_onda3);
   fChain->SetBranchAddress("b0_ped3", &b0_ped3, &b_b0_ped3);
   fChain->SetBranchAddress("b0_pedh3", &b0_pedh3, &b_b0_pedhigh3);
   fChain->SetBranchAddress("b0_charge3", &b0_charge3, &b_b0_charge3);
   fChain->SetBranchAddress("b0_chargepeak3", &b0_chargepeak3, &b_b0_chargepeak3);
   fChain->SetBranchAddress("b0_imax3", &b0_imax3, &b_b0_imax3);
   fChain->SetBranchAddress("b0_tave3", &b0_tave3, &b_b0_tave3);
   fChain->SetBranchAddress("b0_onda4", b0_onda4, &b_b0_onda4);
   fChain->SetBranchAddress("b0_ped4", &b0_ped4, &b_b0_ped4);
   fChain->SetBranchAddress("b0_pedh4", &b0_pedh4, &b_b0_pedhigh4);
   fChain->SetBranchAddress("b0_charge4", &b0_charge4, &b_b0_charge4);
   fChain->SetBranchAddress("b0_chargepeak4", &b0_chargepeak4, &b_b0_chargepeak4);
   fChain->SetBranchAddress("b0_imax4", &b0_imax4, &b_b0_imax4);
   fChain->SetBranchAddress("b0_tave4", &b0_tave4, &b_b0_tave4);
   fChain->SetBranchAddress("b0_onda5", b0_onda5, &b_b0_onda5);
   fChain->SetBranchAddress("b0_ped5", &b0_ped5, &b_b0_ped5);
   fChain->SetBranchAddress("b0_pedh5", &b0_pedh5, &b_b0_pedhigh5);
   fChain->SetBranchAddress("b0_charge5", &b0_charge5, &b_b0_charge5);
   fChain->SetBranchAddress("b0_chargepeak5", &b0_chargepeak5, &b_b0_chargepeak5);
   fChain->SetBranchAddress("b0_imax5", &b0_imax5, &b_b0_imax5);
   fChain->SetBranchAddress("b0_tave5", &b0_tave5, &b_b0_tave5);
   fChain->SetBranchAddress("b0_onda6", b0_onda6, &b_b0_onda6);
   fChain->SetBranchAddress("b0_ped6", &b0_ped6, &b_b0_ped6);
   fChain->SetBranchAddress("b0_pedh6", &b0_pedh6, &b_b0_pedhigh6);
   fChain->SetBranchAddress("b0_charge6", &b0_charge6, &b_b0_charge6);
   fChain->SetBranchAddress("b0_chargepeak6", &b0_chargepeak6, &b_b0_chargepeak6);
   fChain->SetBranchAddress("b0_imax6", &b0_imax6, &b_b0_imax6);
   fChain->SetBranchAddress("b0_tave6", &b0_tave6, &b_b0_tave6);
   fChain->SetBranchAddress("b0_onda7", b0_onda7, &b_b0_onda7);
   fChain->SetBranchAddress("b0_ped7", &b0_ped7, &b_b0_ped7);
   fChain->SetBranchAddress("b0_pedh7", &b0_pedh7, &b_b0_pedhigh7);
   fChain->SetBranchAddress("b0_charge7", &b0_charge7, &b_b0_charge7);
   fChain->SetBranchAddress("b0_chargepeak7", &b0_chargepeak7, &b_b0_chargepeak7);
   fChain->SetBranchAddress("b0_imax7", &b0_imax7, &b_b0_imax7);
   fChain->SetBranchAddress("b0_tave7", &b0_tave7, &b_b0_tave7);
   fChain->SetBranchAddress("b1_onda0", b1_onda0, &b_b1_onda0);
   fChain->SetBranchAddress("b1_ped0", &b1_ped0, &b_b1_ped0);
   fChain->SetBranchAddress("b1_pedh0", &b1_pedh0, &b_b1_pedhigh0);
   fChain->SetBranchAddress("b1_charge0", &b1_charge0, &b_b1_charge0);
   fChain->SetBranchAddress("b1_chargepeak0", &b1_chargepeak0, &b_b1_chargepeak0);
   fChain->SetBranchAddress("b1_imax0", &b1_imax0, &b_b1_imax0);
   fChain->SetBranchAddress("b1_tave0", &b1_tave0, &b_b1_tave0);
   fChain->SetBranchAddress("b1_onda1", b1_onda1, &b_b1_onda1);
   fChain->SetBranchAddress("b1_ped1", &b1_ped1, &b_b1_ped1);
   fChain->SetBranchAddress("b1_pedh1", &b1_pedh1, &b_b1_pedhigh1);
   fChain->SetBranchAddress("b1_charge1", &b1_charge1, &b_b1_charge1);
   fChain->SetBranchAddress("b1_chargepeak1", &b1_chargepeak1, &b_b1_chargepeak1);
   fChain->SetBranchAddress("b1_imax1", &b1_imax1, &b_b1_imax1);
   fChain->SetBranchAddress("b1_tave1", &b1_tave1, &b_b1_tave1);
   fChain->SetBranchAddress("b1_onda2", b1_onda2, &b_b1_onda2);
   fChain->SetBranchAddress("b1_ped2", &b1_ped2, &b_b1_ped2);
   fChain->SetBranchAddress("b1_pedh2", &b1_pedh2, &b_b1_pedhigh2);
   fChain->SetBranchAddress("b1_charge2", &b1_charge2, &b_b1_charge2);
   fChain->SetBranchAddress("b1_chargepeak2", &b1_chargepeak2, &b_b1_chargepeak2);
   fChain->SetBranchAddress("b1_imax2", &b1_imax2, &b_b1_imax2);
   fChain->SetBranchAddress("b1_tave2", &b1_tave2, &b_b1_tave2);
   fChain->SetBranchAddress("b1_onda3", b1_onda3, &b_b1_onda3);
   fChain->SetBranchAddress("b1_ped3", &b1_ped3, &b_b1_ped3);
   fChain->SetBranchAddress("b1_pedh3", &b1_pedh3, &b_b1_pedhigh3);
   fChain->SetBranchAddress("b1_charge3", &b1_charge3, &b_b1_charge3);
   fChain->SetBranchAddress("b1_chargepeak3", &b1_chargepeak3, &b_b1_chargepeak3);
   fChain->SetBranchAddress("b1_imax3", &b1_imax3, &b_b1_imax3);
   fChain->SetBranchAddress("b1_tave3", &b1_tave3, &b_b1_tave3);
   fChain->SetBranchAddress("b1_onda4", b1_onda4, &b_b1_onda4);
   fChain->SetBranchAddress("b1_ped4", &b1_ped4, &b_b1_ped4);
   fChain->SetBranchAddress("b1_pedh4", &b1_pedh4, &b_b1_pedhigh4);
   fChain->SetBranchAddress("b1_charge4", &b1_charge4, &b_b1_charge4);
   fChain->SetBranchAddress("b1_chargepeak4", &b1_chargepeak4, &b_b1_chargepeak4);
   fChain->SetBranchAddress("b1_imax4", &b1_imax4, &b_b1_imax4);
   fChain->SetBranchAddress("b1_tave4", &b1_tave4, &b_b1_tave4);
   fChain->SetBranchAddress("b1_onda5", b1_onda5, &b_b1_onda5);
   fChain->SetBranchAddress("b1_ped5", &b1_ped5, &b_b1_ped5);
   fChain->SetBranchAddress("b1_pedh5", &b1_pedh5, &b_b1_pedhigh5);
   fChain->SetBranchAddress("b1_charge5", &b1_charge5, &b_b1_charge5);
   fChain->SetBranchAddress("b1_chargepeak5", &b1_chargepeak5, &b_b1_chargepeak5);
   fChain->SetBranchAddress("b1_imax5", &b1_imax5, &b_b1_imax5);
   fChain->SetBranchAddress("b1_tave5", &b1_tave5, &b_b1_tave5);
   fChain->SetBranchAddress("b1_onda6", b1_onda6, &b_b1_onda6);
   fChain->SetBranchAddress("b1_ped6", &b1_ped6, &b_b1_ped6);
   fChain->SetBranchAddress("b1_pedh6", &b1_pedh6, &b_b1_pedhigh6);
   fChain->SetBranchAddress("b1_charge6", &b1_charge6, &b_b1_charge6);
   fChain->SetBranchAddress("b1_chargepeak6", &b1_chargepeak6, &b_b1_chargepeak6);
   fChain->SetBranchAddress("b1_imax6", &b1_imax6, &b_b1_imax6);
   fChain->SetBranchAddress("b1_tave6", &b1_tave6, &b_b1_tave6);
   fChain->SetBranchAddress("b1_onda7", b1_onda7, &b_b1_onda7);
   fChain->SetBranchAddress("b1_ped7", &b1_ped7, &b_b1_ped7);
   fChain->SetBranchAddress("b1_pedh7", &b1_pedh7, &b_b1_pedhigh7);
   fChain->SetBranchAddress("b1_charge7", &b1_charge7, &b_b1_charge7);
   fChain->SetBranchAddress("b1_chargepeak7", &b1_chargepeak7, &b_b1_chargepeak7);
   fChain->SetBranchAddress("b1_imax7", &b1_imax7, &b_b1_imax7);
   fChain->SetBranchAddress("b1_tave7", &b1_tave7, &b_b1_tave7);
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
