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

#include "includes/Analysis.h"
#include "includes/langaus.h"
#include "includes/AnaPars.h"
#include "includes/CsvHandler.h"
#include "includes/MipSelection.h"
#include "includes/templ2charge.h"
#include "includes/HistManager.h"
#include "includes/NumberingHelper.h"

#define n_slice 20



// NULLA DI DEFINITIVO 

using namespace std;

CsvHandler CSV;
MipSelection Selection;
HistManager HM;



double pars[4]; //expo+expo(2)
double par_errors[4];

double landau_mpv[2*scintNum][n_slice];
double landau_mpv_err[2*scintNum][n_slice];

void charge_pre_draw(TH1* hist, int histN, int& histSkipFlag)
{
  TSpectrum s(2);
  s.Search(hist);
  double *pks = s.GetPositionX();
  double qpeak = *std::max_element(pks, pks + 2);

  double qmax = qpeak + 500, qmin = qpeak;

  TF1 l = TF1("l", "landau", qmin, qmax);
  l.SetParameters(hist->Integral()/2, 400, 100);
  hist->Fit(&l, "R");
  float pk = l.GetMaximumX(), sigma = l.GetParameter(2);
  TF1 l2 = TF1("l", "landau", pk-60, pk+3*sigma);
  l2.SetParameters(l.GetParameter(0), l.GetParameter(1), sigma);
  hist->Fit(&l2, "R");

  pk = l2.GetMaximumX(); sigma = l2.GetParameter(2);
  TF1 l3 = TF1("l", "landau", pk-1.2*sigma, pk+4*sigma);
  l3.SetParameters(l2.GetParameter(0), l2.GetParameter(1), sigma);
  hist->Fit(&l3, "R");

  TLine *ln = new TLine(130, 0, 130, 4000);
  TString name = TString(hist->GetName());

  int slice_tmp = atoi(((TObjString*)name.Tokenize("_")->At(2))->GetString());

  landau_mpv[histN][slice_tmp] = l.GetParameter(1);
  landau_mpv_err[histN][slice_tmp] = l.GetParError(1);

  hist->Draw("same");
  ln->Draw("same");

}


void createHistBoxes(){

  for(int iSl=0; iSl<n_slice; iSl++){
    HM.AddHistBox(
      "th1f", 2*scintNum, Form("Qmip_zslice_%i", iSl), 
      Form("Charge (MIP) for %i cm < Z < %i cm ", -80 + (int)(iSl*160/n_slice), -80 + (int)((iSl+1)*160/n_slice)), 
      "Q", "pC", 125, 50, 2050,
      &charge_pre_draw
    );
  }
}

void FillVectors(unordered_map<double*, double> dict, int n){
  for(auto &pair: dict) pair.first[n] = pair.second;
}


double *intQ, *teT, *teX2, **chargeEqual, **zetaOffset, **timeOffset, zeta;
int hitSide, hitScint;

void fill_mip(int iScHit) {

  int slice_tmp = (zeta+80)/(160/n_slice);
  if (slice_tmp >= 0 && slice_tmp < n_slice){
    for(int isd = 0; isd < sideNum; isd++){
      HM.Fill1d(Form("Qmip_zslice_%i", slice_tmp), iScHit, intQ[iScHit]);
      HM.Fill1d(Form("Qmip_zslice_%i", slice_tmp), iScHit+scintNum, intQ[iScHit+scintNum]);
    }
  }

}

void Analysis::LoopOverEntries(){


  Long64_t nentries = fChain->GetEntriesFast();
  Long64_t etp = min(maxEvToProcess3, nentries);
  cout << "Number of events to process: " << etp << endl << endl;

  Long64_t nbytes = 0, nb = 0;
   
  // LOOP OVER ENTRIES
  for (Long64_t jentry=0; jentry<etp;jentry++) { //la parte interna al loop andrebbe messa una una funzione così come le parti prima e dopo, così la parte delicata sta in Loop nel .h
    Long64_t ientry = LoadTree(jentry);
    if (jentry%5000 == 0) 
      cout << Form("Processing event n.%lld of %lld: %i%%", jentry, nentries, (int)((double)jentry/nentries * 100)) << endl;
    if (ientry < 0) break;

    nb = fChain->GetEntry(jentry);   
    nbytes += nb;
    
    list<double **> arr_list = {&intQ, &teT, &teX2};

    for(double** &arr: arr_list) {
      *arr = new double[2*scintNum];
    }    
  
    int skipFlag = 0, m = 0, iScHit = -1;

    // LOOP OVER HITS
    for(int hit=0; hit<nCry; hit++){
      hitSide =iSide[hit];
      hitScint = iScint[hit];

      int hitN = hitSide*scintNum + hitScint; 

      double chCal = enableOfflineEq ? chEqReference/chargeEqual[hitSide][hitScint] : 1;
  
      if (Selection.isSaturated(Vmax[hit])) {skipFlag = 1; continue;}

      FillVectors({
        {intQ, Qval[hit] * chCal}, 
        {teX2, templChi2[hit]}, {teT, templTime[hit] - timeOffset[hitSide][hitScint]},
      }, hitN);

    }
  
    for(int isc = 0; isc <scintNum; isc++){

      if (skipFlag) {continue;}

      if ( Selection.hitPrecheck(isc, iScint, nCry) && Selection.isChargeGood(intQ, isc) ) {
        if ( Selection.isX2Good(teX2, isc) && !Selection.isShared(intQ, isc) ) { iScHit = isc; m++; }
      }

    }


    if ( m != 1 ) {continue;}

    if ( !Selection.isTimeGood(teT[iScHit])) {continue;}

    zeta = ( teT[iScHit] - teT[scintNum+iScHit])*scintVp/2- zetaOffset[0][iScHit];

    if ( !Selection.isZetaGood(zeta) ) {continue;}
    //if ( !Selection.mipCutG(intQ, zeta, iScHit) ) {continue;}

    fill_mip(iScHit);
    
  }

}

void Analysis::ProcessPlots(){

  double z_points[n_slice], z_points_err[n_slice];
  for(int iSl=0; iSl<n_slice; iSl++){
    z_points[iSl] = -80 + (160/n_slice)*iSl;
    z_points_err[iSl] = (double)(160/n_slice)/3;
  }

  outFile->cd();
  outFile->mkdir("zeta_q");
  for(int iSc=0; iSc<scintNum; iSc++){
    TCanvas *qz_c = new TCanvas(Form("qz_c_%i", iSc), "qz_c");
    qz_c->cd();
    for(int iSd=0; iSd<2; iSd++){
      TGraphErrors *qz_side0 = new TGraphErrors(n_slice, z_points, landau_mpv[iSc], z_points_err, landau_mpv_err[iSc]);
      TGraphErrors *qz_side1 = new TGraphErrors(n_slice, z_points, landau_mpv[iSc+scintNum], z_points_err, landau_mpv_err[iSc+scintNum]);
 

      TF1 l = TF1("l", "expo(0)+expo(2)", -80, 80);
      l.SetParLimits(3, 0.01, 0.09);
      l.SetParLimits(1, 0.0001, 0.03);
      
      qz_side1->Fit(&l, "R", "", -65, 75);
      TLatex *lt1 = new TLatex(-40, 1300 - 300, Form("BAL (side %i): %.1f +/- %.1f cm ", 1, 1/l.GetParameter(1), 1/(l.GetParameter(1)*l.GetParameter(1)) *l.GetParError(1) ));
      TLatex *lt2 = new TLatex(-40, 1300 - 350, Form("TAL (side %i): %.1f +/- %.1f cm ", 1, 1/l.GetParameter(3), 1/(l.GetParameter(3)*l.GetParameter(3)) *l.GetParError(3) ));

      l.SetParLimits(3, -0.09, -0.01);
      l.SetParLimits(1, -0.03, -0.0001);
      qz_side0->Fit(&l, "R", "", -75, 66);

      TLatex *lt3 = new TLatex(-40, 1300 - 400, Form("BAL (side %i): %.1f +/- %.1f cm ", 0, -1/l.GetParameter(1), 1/(l.GetParameter(1)*l.GetParameter(1)) *l.GetParError(1) ));
      TLatex *lt4 = new TLatex(-40, 1300 - 450 , Form("TAL (side %i): %.1f +/- %.1f cm", 0, -1/l.GetParameter(3), 1/(l.GetParameter(3)*l.GetParameter(3)) *l.GetParError(3) ));

      lt1->SetTextSize(0.03);
      lt2->SetTextSize(0.03);
      lt3->SetTextSize(0.03);
      lt4->SetTextSize(0.03);

      qz_side0->SetName("qz_side0");
      qz_side1->SetName("qz_side1");

      qz_side0->GetXaxis()->SetRange(-75, 75);
      qz_side1->GetXaxis()->SetRange(-75, 75);


      qz_side0->GetYaxis()->SetRange(-50, 1300);
      qz_side1->GetYaxis()->SetRange(-75, 1300);

      qz_side0->Draw("APE");
      qz_side1->Draw("p e same");
      
      lt1->Draw();
      lt2->Draw();

      lt3->Draw();
      lt4->Draw();

    }
    outFile->cd("zeta_q");
    qz_c->Write();
  }



}

void Analysis::Loop(){
  if (fChain == 0) return;

  HM.SetOutFile(outFile);
  HM.SetNamerFun(&NamerMatrix);
  createHistBoxes();

  cout<<"Retrieving calibration data from [.../" + lutPrefix3p + "] :"<<endl;
  timeOffset  = CSV.InitMatrix(2, scintNum);
  chargeEqual = CSV.InitMatrix(2, scintNum);
  zetaOffset  = CSV.InitMatrix(1, scintNum);
  CSV.Read(CSV.GetFirstFile("../" + lutPrefix3p + calName + lutChEqName + "*"),      ',', chargeEqual, 2, scintNum);
  CSV.Read(CSV.GetFirstFile("../" + lutPrefix3p + calName + lutTimeOffsName + "*"),  ',', timeOffset,  2, scintNum);
  CSV.Read(CSV.GetFirstFile("../" + lutPrefix3p + calName + lutZetaOffName + "*"),   ',', zetaOffset,  1, scintNum);
  cout<<"...done"<<endl<<endl;


  Analysis::LoopOverEntries();
  HM.ProcessBoxes(); 

  Analysis::ProcessPlots();

  outFile->Close();

}

void attLength(){
    Analysis::Run("../data/step2/run205_s2.root", "../data/step3/att205.root", "run205", "run205", -1);
}
