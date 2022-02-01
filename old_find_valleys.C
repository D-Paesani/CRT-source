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

#include "includes/Analysis.h"
#include "includes/langaus.h"
#include "includes/AnaPars.h"
#include "includes/CsvHandler.h"
#include "includes/MipSelection.h"
#include "includes/templ2charge.h"
#include "includes/HistManager.h"

using namespace std;

CsvHandler CSV;
MipSelection Selection; 
HistManager HM;

Long64_t nentries, nbytes, nb, ientry, jentry, jTrig_out;
double **timeOffset, **timeOffset_out, **zetaOffset, **zetaOffset_out, **chargeEqual, **chargeEqualErr, **chargeEqual_out, **chargeEqualErr_out;
double tDiff = 0, zeta=0;
double *intQ, *pkV, *teQ, *teT, *teA, *teB, *teX2, *ped;
list<double ** > arrayList = {&intQ, &pkV, &teQ, &teT, &teA, &teB, &teX2, &ped};
void InitVectors() { for(double** &arr: arrayList) { *arr = new double[2*scintNum](); } }
int iSc_out; double_t Z_out; double_t Q_out[2], X2_out[2], T_out[2];
TTree *CRTs3;

double qvalleys[20];

void qmip_zsl_proc(TH1* histObj, int histN, TString& histTag, int& histSkipFlag) {

  int iSd = (histN/8)%2, iSc = histN%8, zsl = histN/16;
  histTag = Form("%d_%d_slice%d",  iSd, iSc, zsl);

  TSpectrum s(2);
  s.Search(histObj, 2, "");
  double *pks = s.GetPositionX(); 
  double qpeak2 = *std::max_element(pks, pks + 2);
  double qpeak1 = *std::min_element(pks, pks + 2);
/*  TH1F *temp = (TH1F*)histObj->Clone();
  temp->GetXaxis()->SetLimits(qpeak1, qpeak2);
  temp->Scale(-1);
  TSpectrum ss(1);
  ss.Search(temp);
  qvalleys[zsl] = *(ss.GetPositionX());*/
}

void createHistBoxes() {
  HM.HistBoxes = {
    HM.AddHistBox("qmip_zsl",    1,    2*8*20,   "Mip charges",      "charge", "pC",    qBins, qFrom, qTo,      1, 0, 0,   &qmip_zsl_proc),
  };
}


void fill_mip(int iScHit) {

  double z_from_left = zeta + 80;
  if (z_from_left < 0 || z_from_left>160) return;
  int zsl = z_from_left/8;

  HM.Fill1d("qmip_zsl", iScHit+16*zsl, intQ[iScHit]);
  HM.Fill1d("qmip_zsl", iScHit+scintNum+16*zsl, intQ[iScHit+scintNum]);

}



void Analysis::ProcessPlots() {}


void Analysis::LoopOverEntries() {

  nentries = fChain->GetEntriesFast(); 
  Long64_t etp = min(maxEvToProcess3, nentries);
  cout << "Number of events to process: " << etp << endl << endl;
  nbytes = 0, nb = 0;

  for (jentry = 0; jentry < etp; jentry++) {

    ientry = LoadTree(jentry);

    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   
    nbytes += nb;

    if (!(jentry%10000)) {cout << Form( "     processing evt %lld / %lld  ( %.0f%% )", jentry, etp, (float)(100*jentry/etp) ) << endl;}

    InitVectors();
    int skipFlag = 0, m = 0, iScHit = -1;

    for(int hit = 0; hit < nCry; hit++){

      int hitSide=iSide[hit], hitScint = iScint[hit], hitN = hitSide*scintNum + hitScint; 
      double chCal = chEqReference/chargeEqual[hitSide][hitScint]; // si può anche preparare fuori dal for
      //chCal = (chCal>0.8 && chCal<1.2)?chCal:1; // da togliere
      chCal = 1; //// no offline eq

      intQ[hitN] = Qval[hit]*chCal;
      pkV[hitN] = Vmax[hit];
      ped[hitN] = pedL[hit];
      teA[hitN] = templFit[hit][0];
      teQ[hitN] = templ2charge(teA[hitN])*chCal;
      teB[hitN] = templFit[hit][2];
      teT[hitN] = templTime[hit] - timeOffset[hitSide][hitScint];
      teX2[hitN] = templChi2[hit];

      if (Selection.isSaturated(pkV[hitN])) {skipFlag = 1; continue;}
    }

    if (skipFlag) {continue;}

    for(int isc = 0; isc < scintNum; isc++) {

      if ( Selection.hitPrecheck(isc, iScint, nCry) && Selection.isChargeGood(intQ, isc) ) {

        if ( Selection.isX2Good(teX2, isc) && !Selection.isShared(intQ, isc) ) { iScHit = isc; m++; }
      }
    }

    if ( m != 1 ) {continue;}

    if ( !Selection.isTimeGood(teT[iScHit])) {continue;}

    tDiff = teT[iScHit] - teT[scintNum+iScHit];
    zeta = tDiff*scintVp/2-zetaOffset[0][iScHit];

    if ( !Selection.isZetaGood(zeta) ) {continue;}

    fill_mip(iScHit);
  }

  cout<<endl;
}


void Analysis::Loop(){

  cout<<endl<<endl<<"::::::::::::::::::::: CRT analysis step 3 :::::::::::::::::::::"<<endl<<endl;

  if (fChain == 0) return;

  cout<<"Creating histograms:"<<endl;
  HM.SetOutFile(outFile);
  createHistBoxes();
  cout<<"...done"<<endl<<endl;

  cout<<"Retrieving calibration data from [" + lutPrefix3p + "] ..."<<endl;
  timeOffset  = CSV.InitMatrix(2, scintNum); timeOffset_out = CSV.InitMatrix(2, scintNum);
  chargeEqual = CSV.InitMatrix(2, scintNum); chargeEqual_out = CSV.InitMatrix(2, scintNum);
  chargeEqualErr = CSV.InitMatrix(2, scintNum); chargeEqualErr_out = CSV.InitMatrix(2, scintNum);
  zetaOffset  = CSV.InitMatrix(1, scintNum); zetaOffset_out  = CSV.InitMatrix(1, scintNum);
  CSV.Read(CSV.GetFirstFile(lutPrefix3p + calName + lutChEqName + "*"),      ',', chargeEqual, 2, scintNum);
  CSV.Read(CSV.GetFirstFile(lutPrefix3p + calName + lutTimeOffsName + "*"),  ',', timeOffset,  2, scintNum);
  CSV.Read(CSV.GetFirstFile(lutPrefix3p + calName + lutZetaOffName + "*"),   ',', zetaOffset,  1, scintNum);
  cout<<"...done"<<endl<<endl;

  outFile->cd();

  Analysis::LoopOverEntries();
  HM.ProcessBoxes();

  cout<<"...done"<<endl;

  outFile->Close();

  cout<<endl<<endl<<"::::::::::::::::::::: analysis done :::::::::::::::::::::"<<endl<<endl;
}


int main(int argc, char*argv[]) { 

  if (argc != 5) {
    printf("Usage: %s [infile_name] [outfile_name] [run_name] [calib_name]\n", argv[0]);
    exit(-1);
  }

  Analysis::Run(argv[1], argv[2], argv[3], argv[4], -1);
}
