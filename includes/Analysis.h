#pragma once

#ifndef Analysis_h
#define Analysis_h

#define ROOT_CLASS ana
#include "ana.h"
#include "ana.C"

#define TREE_NAME "CRT"



class Analysis: public ROOT_CLASS
{
public:
  Analysis(TString, TFile*, TString, TString, TString);
  virtual void Loop() override; //if not needed, fill with ROOT_CLASS::Loop();
  void LoopOverEntries();
  void ProcessHistos();
  void ProcessPlots();
  TFile *outFile;
  TString inFileName;
  TString runName;
  TString calName;
  TString modulSel;
  TTree *inTree;
  Int_t iMod[16]; //////sarà poi da 32
  TBranch *b_iMod;

  TTree * GetTree(TString infileName) {
  TTree *tree;
  TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject(infileName.Data());
      if (!f || !f->IsOpen()) {
         f = new TFile(infileName.Data());
      }
      f->GetObject(TREE_NAME, tree);
  inTree = tree;
  return tree;
  }


};

#endif

Analysis::Analysis(TString infileName, TFile *fileout, TString runN, TString calN, TString modulSelect) : ROOT_CLASS(GetTree(infileName)) { 

  outFile = fileout;
  inFileName = infileName;
  runName = runN;
  calName = calN;
  modulSel = modulSelect;

  if (modulSel != "") { inTree->SetBranchAddress("iMod", iMod, &b_iMod); } //questo può anche stare nella tree class e darà unknown branch se non lo trova

};



