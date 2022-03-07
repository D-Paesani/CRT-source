#include <TTree.h>
#include <TFile.h>
#include <TList.h>

void mergeTrees(){
  TFile *files[12];
  for(int i=0; i<12; i++) files[i] = new TFile(Form("run_97_%i.root", i));

  TFile *outf = new TFile("run97.root", "CREATE");

  TList *lst = new TList;

  for(auto &f: files){
    TTree *tree = (TTree*)f->Get("tree");
    lst->Add(tree);
  }

  TTree *newtree = TTree::MergeTrees(lst);
  newtree->SetName("tree");
  outf->cd();
  newtree->Write();
}
