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


using namespace std;

CsvHandler CSV;
MipSelection Selection;
HistManager HM;


//Pars
  TString run_name  = "run182";
  TString in_path   = "./data/step2/";
  TString out_path  = "./data/template/";
  TString cal_name  = "run182"; //not used. 
  TString out_pre   = "genTemp_";

  TString argv1 = in_path + run_name + "_s2.root";
  TString argv2 = out_path + out_pre + run_name + ".root";
  TString argv3 = run_name;
  TString argv4 = cal_name;

  int     start_time = 0;
  int     stop_time = 300;
  double  templ_offs = 100;

  int     ti_bins = 1600;
  double  ti_from = 0;
  double  ti_to = 800;
  int     amp_bins = 1200;
  double  amp_from = -0.1;
  double  amp_to = 1.2;

  double  charge_min = 300;
  double  charge_max = 800;
  double  chi2_max = 10;

  double  time_cut_low = -100;
  double  time_cut_high = 100;

  const    TString preCut = Form("Qval > %f && Qval < %f && templChi2 > 0 && templChi2 < %f", charge_min, charge_max, chi2_max);

  Long64_t max_evts = 100000000;
//Pars


TDirectory* spline_dir;
TDirectory* templDraw_dir;
TDirectory* splineDraw_dir;
TDirectory* preProcessing_dir;

Long64_t nentries, nbytes, nb, ientry, jentry, jTrig_out;
double **timeOffset, **zetaOffset, **zetaOffset_out, **chargeEqual, **chargeEqualErr, **chargeEqual_out, **chargeEqualErr_out;
double teTOffset[2*scintNum], pkTOffset[2*scintNum];
double tDiff = 0, zeta=0;
double *intQ, *pkV, *teQ, *teT, *teA, *teB, *teX2, *ped, *tPk;
list<double ** > arrayList = {&intQ, &pkV, &teQ, &teT, &teA, &teB, &teX2, &ped, &tPk};
void InitVectors() { for(double** &arr: arrayList) { *arr = new double[2*scintNum](); } }
int iSc_out; double_t Z_out; double_t Q_out[2], X2_out[2], T_out[2];


void fuzzyTemp_proc(TH1* histObj, int histN, int& histSkipFlag) {

  TString histName = histObj->GetName();

  TCanvas *templDraw_can = new TCanvas(histName);
  templDraw_dir->cd();  
  histObj->SetTitle(histName);
  histObj->SetDrawOption("zcol");
  histObj->Draw("zcol");
  templDraw_can->SetLogz();
  templDraw_can->Write();

  TCanvas *spline_can = new TCanvas(histName + "_spline"); 
  spline_can->cd();

  TProfile *teProf = ((TH2*)histObj)->ProfileX();
  TSpline5 *teSpline = new TSpline5(teProf);
  TGraphErrors *teSplGr = (TGraphErrors*)(((TH2*)histObj)->ProfileX());

  teProf->SetName(histName + "_profile");  
  teSpline->SetName(histName + "_spline");
  teSpline->SetLineColor(kOrange);
  teProf->Draw();
  teSpline->Draw("L same");
  
  splineDraw_dir->cd();  
  spline_can->Write();

  teSplGr->SetMarkerStyle(8);
  teSplGr->SetMarkerSize(.5);
  teSplGr->SetMarkerColor(kBlue);
  teSplGr->SetLineColor(kOrange);

  spline_dir->cd();
  teSplGr->Write(histName + "_graph");
}


void teTimes_proc(TH1* histObj, int histN, int& histSkipFlag) {

  TString histName = histObj->GetName();

  TCanvas cc(histName + "_cut_teTime", histName + "_cut_teTime"); cc.cd();
  histObj->SetTitle(histName);
  histObj->Draw();
  
  TLine l  = TLine(time_cut_low,  histObj->GetMaximum(), time_cut_low, 0); 
  TLine ll = TLine(time_cut_high, histObj->GetBinContent(histObj->GetMaximumBin()), time_cut_high, 0);   
  l.SetLineColor(kRed); ll.SetLineColor(kRed); l.Draw("same"); ll.Draw("same");

  cc.Write();

}

void createHistBoxes() {
    
  HM.AddHistBox("th2f", 2*scintNum, "fuzzyTempl","Fuzzy template", "Time", "ns", "Normalised pulse", "", ti_bins, ti_from, ti_to, amp_bins, amp_from, amp_to, &fuzzyTemp_proc);
  HM.AddHistBox("th1f", 2*scintNum, "teTimes", "Reco times",     "Time", "ns",  600, -150, 150, &teTimes_proc);
  HM.AddHistBox("th1f", 2*scintNum, "tiDiff", "Tpeak - Ttempl", "Time", "ns",  600, -150, 150);

}


void Analysis::ProcessPlots() {}


void Analysis::LoopOverEntries() {

  preProcessing_dir->cd();

  for (int k = 0; k < 2*scintNum; k++) {

    int iSd = GetSide(k), iSc = GetScint(k); 
    TString histTag = Form("_%d_%d",  iSd, iSc);

    TH1F teT_temp = TH1F("teT_temp", "teT_temp", 100, 10, 500);
    TH1F pkT_temp = TH1F("pkT_temp", "pkT_temp", 100, 10, 500);

    fChain->Draw("templTime>>teT_temp",  preCut + Form(" && iSide == %i && iScint == %i", iSd, iSc), "goff");
    fChain->Draw("Tval>>pkT_temp",       preCut + Form(" && iSide == %i && iScint == %i", iSd, iSc), "goff");

    gStyle->SetOptFit(1);

    double tpeak = teT_temp.GetBinCenter(teT_temp.GetMaximumBin());
    double tmax = tpeak + 30, tmin = tpeak - 30;
    TF1 timeFit = TF1("g", "gaus", tmin, tmax); timeFit.SetParameter(1, tpeak); timeFit.SetParameter(2, 2);
    teT_temp.Fit(&timeFit, "R");
    teT_temp.Fit(&timeFit, "R");
    teTOffset[k] = timeFit.GetParameter(1); 
    teT_temp.Write( "teTime" +  histTag);

    tpeak = pkT_temp.GetBinCenter(pkT_temp.GetMaximumBin());
    tmax = tpeak + 30, tmin = tpeak - 30;
    timeFit = TF1("g", "gaus", tmin, tmax); timeFit.SetParameter(1, tpeak); timeFit.SetParameter(2, 2);
    pkT_temp.Fit(&timeFit, "R");
    pkT_temp.Fit(&timeFit, "R");
    pkTOffset[k] = timeFit.GetParameter(1); 
    pkT_temp.Write("pkTime" + histTag);

  }

  cout<<endl<<endl<<endl;

  nentries = fChain->GetEntriesFast(); 
  Long64_t etp = min(max_evts, nentries);
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

      int hitSide = iSide[hit], hitScint = iScint[hit], hitN = GetChan(hitSide, hitScint);
      //double chCal = chEqReference/chargeEqual[hitSide][hitScint];
      //chCal = (chCal>0.8 && chCal<1.2)?chCal:1;

      double _intQ = Qval[hit];//*chCal;
      double _pkV = Vmax[hit];
      double _pkT = Tval[hit] - pkTOffset[hitN];
      double _ped = pedL[hit];
      double _teA= templFit[hit][0];
      double _teQ = templ2charge(_teA);//*chCal;
      double _teB = templFit[hit][2];
      double _teT = templTime[hit] - teTOffset[hitN];
      double _teX2 = templChi2[hit];

      HM.Fill1d("tiDiff", hitN, _pkT - _teT);
      HM.Fill1d("teTimes", hitN, _teT);

      if ( //add rms
          !Selection.isSaturated(_pkV)   &&
          _intQ > charge_min             &&
          _intQ < charge_max             &&
          _teX2 > 0                      &&
          _teX2 < chi2_max               &&
          _teT < time_cut_high           &&
          _teT > time_cut_low            &&
          _pkT - _teT > -5               &&
          _pkT - _teT < 5 
        ) 
      { } else { continue; }  
      
      for(int tt=start_time; tt<stop_time; tt++){ 

        double ttt = (double)ana::time[tt] - teTOffset[hitN] - _teT + templ_offs;
        //double ampl = ana::wave[hit][tt]/_teA;  // ma fa cagare, bisogna mettere dei tagli
        double ampl = ana::wave[hit][tt] / _pkV; //bisogna migliorare la normalizzazione

        HM.Fill2d( "fuzzyTempl", hitN, ttt, ampl);
        
      }
    } 
  }
  cout<<endl;
}


void Analysis::Loop(){

  cout<<endl<<endl<<"::::::::::::::::::::: CRT gen template :::::::::::::::::::::"<<endl<<endl;

  if (fChain == 0) return;

  cout<<"Creating histograms:"<<endl;
  HM.SetOutFile(outFile);
  HM.SetNamerFun(&NamerMatrix);
  createHistBoxes();
  cout<<"...done"<<endl<<endl;

  // cout<<"Retrieving calibration data from [" + lutPrefix3p + "] ..."<<endl;
  // timeOffset  = CSV.InitMatrix(2, scintNum);
  // chargeEqual = CSV.InitMatrix(2, scintNum);
  // chargeEqualErr = CSV.InitMatrix(2, scintNum);
  // chargeEqual_out = CSV.InitMatrix(2, scintNum);
  // chargeEqualErr_out = CSV.InitMatrix(2, scintNum);
  // zetaOffset  = CSV.InitMatrix(1, scintNum);
  // zetaOffset_out  = CSV.InitMatrix(1, scintNum);
  // CSV.Read(CSV.GetFirstFile(lutPrefix3p + calName + lutChEqName + "*"),      ',', chargeEqual, 2, scintNum);
  // CSV.Read(CSV.GetFirstFile(lutPrefix3p + calName + lutTimeOffsName + "*"),  ',', timeOffset,  2, scintNum);
  // //CSV.Read(CSV.GetFirstFile(lutPrefix3p + calName + lutZetaOffName + "*"),   ',', zetaOffset,  1, scintNum);
  // cout<<"...done"<<endl<<endl;

  spline_dir = outFile->mkdir("splines");
  splineDraw_dir = outFile->mkdir("profiles");
  templDraw_dir = outFile->mkdir("fuzzyTemplDraw");
  preProcessing_dir = outFile->mkdir("preProcessing");

  Analysis::LoopOverEntries();
  HM.ProcessBoxes(); 
  Analysis::ProcessPlots();

  outFile->Close();

  cout<<endl<<endl<<"::::::::::::::::::::: template done :::::::::::::::::::::"<<endl<<endl;
}


int main(int argc, char*argv[]) { 

  if (argc != 5) {
    printf("Usage: %s [infile_name] [outfile_name] [run_name] [calib_name]\n", argv[0]);
    exit(-1);
  }

  Analysis::Run(argv1, argv2, argv3, argv4, -1);
}



