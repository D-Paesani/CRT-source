#include <fstream>
#include <chrono>
#include <TLine.h>
#include <iostream>
#include <list>

#include "TApplication.h"
#include "TSpectrum.h"
#include "TGraphErrors.h"
#include "TGraphSmooth.h"
#include "TSpline.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TH1F.h"
#include "TF1.h"

#include "HistManager.h"
#include "NumberingHelper.h"


using namespace std;

HistManager HM;

void namerHistBox(int n, TString& hTag, TString& hTitleTag) {

  hTag = Form("hhhhhhhh_%d_", n); 
  hTitleTag = Form("[%d] titolo: ", n);
}


void usage(){

  TFile *ff = new TFile("aaa.root", "RECREATE");
  HM.SetOutFile(ff);
  
  HM.AddHistBox("th1f", 5, "histo1", "histo1title", "XXX", "s", 10, 0, 10);

  HM.SetNamerFun(&namerHistBox);
  HM.SetProcessFun(HM.GetProcDef());
  
  HM.AddHistBox("th1f",  5, "histo2", "histo2title", "XXX", "s", 10, 0, 10);
  HM.AddHistBox("th1f",  5, "histo3", "ThisdefClas", "XXX", "s", 10, 0, 10, HM.GetProcDef());
  HM.AddHistBox("th2f",  5, "histo4", "ThisdefUser", "XXX", "s", "YYY", "s", 10, 0, 10, 10, 0, 10, HM.GetProcUser(), HM.GetNamerUser());
  HM.AddHistBox("th2f", 16, "MATRIX", "matrix     ", "XXX", "s", "YYY", "s", 10, 0, 10, 10, 0, 10, HM.GetProcUser(), &NamerMatrix);
  HM.AddHistBox("th2f", 16, "ARRAYS", "arrayyyyyyy", "XXX", "s", "YYY", "s", 10, 0, 10, 10, 0, 10, HM.GetProcUser(), &NamerArray);


  HM.Fill1d("histo1", 2, 2);
  HM.Fill1d("histo1", 2, 2);
  HM.Fill1d("histo1", 2, 3);
  HM.Fill2d("histo4", 3, 2, 2);
  HM.Fill2d("histo4", 3, 4, 4);

  HM.GetHistos("histo1")[2]->Fill(9);
  HM.GetHist("histo4", 3)->Fill(5, 5);

  HM.ProcessBoxes(); 

  HM.AddHistBox("th2f", 16, "last", "last", "XXX", "s", "YYY", "s", 10, 0, 10, 10, 0, 10, HM.GetProcUser(), &NamerMatrix);
  HM.ProcessBox("last");


  HM.CloseOutFile();

}

