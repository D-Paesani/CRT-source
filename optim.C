#include <fstream>
#include <chrono>
#include <TLine.h>
#include <iostream>
#include <list>
#include <algorithm>
#include <unistd.h>
#include <stdio.h>

#include "TApplication.h"
#include "TSpectrum.h"
#include "TGraphErrors.h"
#include "TGraphSmooth.h"
#include "TSpline.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TH1F.h"
#include "TF1.h"
#include "TButton.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TApplication.h"
#include "TAxis.h"
#include "TMarker.h"
#include "TRandom.h"

#include "includes/Analysis.h"
#include "includes/langaus.h"
#include "includes/AnaPars.h"
#include "includes/CsvHandler.h"
#include "includes/MipSelection.h"
#include "includes/templ2charge.h"
#include "includes/HistManager.h"
#include "includes/NumberingHelper.h"

using namespace std;

//Pars

  Long64_t max_evts = 50000000;

  #define doFitTemplate 0 // <------
  #define isRun182 1

  TString run_name  = "run183";
  TString in_path   = "data/step2/";
  TString out_path  = "data/template/";
  TString out_pre   = "provaSplines_";

  TString splines_path   = "data/template/provaSplines_run182.root";
  TString splines_format   = "splines/fuzzyResamp_%d_%d_spline";

  TString argv1 = in_path + run_name + "_s2.root";
  TString argv2 = out_path + out_pre + run_name + ".root";
  TString argv3 = run_name;
  TString argv4 = "";

  const int samNo = 1024;
  double ex[samNo] = {0.05}; 
  const double dig_res = 0.5;
  double ey[samNo] = {5};
  double  templ_offs = 300;

  double resolution = 0, resolution_error = 0;
  double constant_fraction_current;

  const int digi_time = 4;
  int     ti_bins = 1600;
  double  ti_from = 0;
  double  ti_to = 800;
  int     resam_fact = (int)ti_bins/(ti_to-ti_from);

  int     amp_bins = 1200;
  double  amp_from = -0.1;
  double  amp_to = 1.2;

  double  charge_min = 300;
  double  charge_max = 800;


//Pars

double teTOffset[2*scintNum] = {0}, pkTOffset[2*scintNum] = {0};

TDirectory *spline_dir, *splineGr_dir, *samples_dir, *templDraw_dir, *splineDraw_dir, *templResampDraw_dir, *preProcessing_dir;

CsvHandler CSV;
MipSelection Selection;
HistManager HM;


void fuzzyTemp_proc(TH1* histObj, int histN, int& histSkipFlag) { //questa genera il tempalte tal TH2

  TString histName = histObj->GetName(); //histObj è il ptr allo i-esimo TH2F 

  TCanvas *templDraw_can = new TCanvas(histName);
  templDraw_dir->cd();  
  histObj->SetTitle(histName);
  histObj->SetDrawOption("zcol");
  histObj->Draw("zcol");
  templDraw_can->SetLogz();
  templDraw_can->Write();

  TCanvas *spline_can = new TCanvas(histName + "_spline"); 
  spline_can->cd();

  TProfile *teProf = ((TH2*)histObj)->ProfileX(); //faccio il profio
  TSpline5 *teSpline = new TSpline5(teProf); //spline from profile
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

  splineGr_dir->cd();
  teSplGr->Write(histName + "_graph");

  spline_dir->cd();
  teSpline->Write();
}

void times_proc(TH1* histObj, int histN, int& histSkipFlag) { 

   gStyle->SetOptFit(1);
   double tpeak = histObj->GetBinCenter(histObj->GetMaximumBin());
   double tmax = tpeak + 2, tmin = tpeak - 2;
   TF1 timeFit = TF1("g", "gaus", tmin, tmax); timeFit.SetParameter(1, tpeak); timeFit.SetParameter(2, 2);
   histObj->Fit(&timeFit, "R");
   double mean = timeFit.GetParameter(1), sigma = timeFit.GetParameter(2);
   timeFit = TF1("g", "gaus", mean - 3*sigma, mean + 3*sigma); timeFit.SetParameter(1, mean); timeFit.SetParameter(2, sigma);
   histObj->Fit(&timeFit, "R");
   TString name(histObj->GetName());
   if (name.Contains("tdiff") && (!name.Contains("fit")) && histN == 0) {
     resolution = timeFit.GetParameter(2);
     resolution_error = timeFit.GetParError(2);
     cout << "CF CORRENTE" << constant_fraction_current << endl;
     cout << "RESO CORRENTE" << timeFit.GetParameter(2) << endl;
   }
   cout<<"----->"<<histObj->GetName()<<endl<<endl;
}


void createHistBoxes() {

  HM.AddHistBox("th1f", 2*scintNum, "chargeRaw", "Q", "Q", "pC", qBins, qFrom, qTo);
  HM.AddHistBox("th1f", 2*scintNum, "chargeMip", "Q", "Q", "pC", qBins, qFrom, qTo);
  HM.AddHistBox("th2f", 2*scintNum, "fuzzy", "Fuzzy template", "Time", "ns", "Normalised pulse", "", ti_bins, ti_from, ti_to, amp_bins, amp_from, amp_to, &fuzzyTemp_proc);
  HM.AddHistBox("th2f", 2*scintNum, "fuzzyResamp", "Fuzzy template with resampling", "Time", "ns", "Normalised pulse", "", ti_bins, ti_from, ti_to, amp_bins, amp_from, amp_to, &fuzzyTemp_proc);
  HM.AddHistBox("th1f", 2*scintNum, "recoTime", "pseudotime", "Time", "ns",  1000, -100, 900);
  HM.AddHistBox("th1f", 2*scintNum, "recoTimeF", "time Fit", "Time", "ns",  1000, -100, 900);
  HM.AddHistBox("th1f", 2*scintNum, "timeModBin", "Flatness over bin", "Time", "ns",  100, 0, 1.25);
  HM.AddHistBox("th1f", 2*scintNum, "timeModBinF", "Flatness over bin", "Time", "ns",  100, 0, 1.25);
  HM.AddHistBox("th2f", 2*scintNum, "slewing", "slewing", "Q", "pC", "T", "s", qBins, qFrom, qTo, 800, 200, 400);
  HM.AddHistBox("th2f", 2*scintNum, "slewingFit", "slewing", "Q", "pC", "T", "s", qBins, qFrom, qTo, 800, 200, 400);
  HM.AddHistBox("th1f", scintNum, "tsum", "td", "Time", "ns",  1000, 0, 1000, &times_proc, &NamerArray);
  HM.AddHistBox("th1f", scintNum, "tdiff", "td", "Time", "ns",  1600, -20, 20, &times_proc, &NamerArray);
  HM.AddHistBox("th1f", scintNum, "tdifffit", "tdfit", "Time", "ns",  1600, -20, 20, &times_proc, &NamerArray);
  HM.AddHistBox("th1f", 2*scintNum, "bLineRms", "Base line rms", "", "mV", 200, 0.0, 10);
  HM.AddHistBox("th1f", 2*scintNum, "bLine", "Base line", "", "mV", 200, -5, 5); 
  HM.AddHistBox("th1f", 2*scintNum, "chi2", "chi2", "", "mV", 200, 0, 10); 



}



void Analysis::LoopOverEntries() {

  TSpline5 *spls[2*scintNum];
  TFile *spline_file;
  if(doFitTemplate) {
    spline_file = new TFile(splines_path);
    cout<<"Loading splines..."<<endl;
    for ( int k = 0; k < 2*scintNum; k++) {
      spls[k] = (TSpline5*)spline_file->Get(Form(splines_format, GetSide(k), GetScint(k))); 
    };
    cout<<"...done"<<endl<<endl;
  }

  preProcessing_dir->cd();
  for (int k = 0; k < 2*scintNum; k++) {

    int iSd = GetSide(k), iSc = GetScint(k); 

    int hitN = GetChan(iSd, iSc);

    if(isRun182) {
        if (iSd == 0 && iSc == 0) {}
        else if (iSd == 0 && iSc == 2) {hitN = GetChan(1,0);}
        else {continue;}
    }

    TString histTag = Form("_%d_%d",  iSd, iSc);

    TH1F pkT_temp = TH1F("pkT_temp", "pkT_temp", 100, 10, 500);

    fChain->Draw("Tval>>pkT_temp",  Form("Qval > 200 && iSide == %i && iScint == %i", iSd, iSc), "goff");

    gStyle->SetOptFit(1);

    double tpeak = pkT_temp.GetBinCenter(pkT_temp.GetMaximumBin());
    double tmax = tpeak + 30, tmin = tpeak - 30;
    TF1 timeFit = TF1("g", "gaus", tmin, tmax); timeFit.SetParameter(1, tpeak); timeFit.SetParameter(2, 2);
    pkT_temp.Fit(&timeFit, "R");
    pkT_temp.Fit(&timeFit, "R");
    pkTOffset[hitN] = timeFit.GetParameter(1);
    pkT_temp.Write("pkTime" + histTag);
  }

  Long64_t nentries, nbytes, nb, ientry, jentry;
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

    double IntQ[2*scintNum] = {0}, PkV[2*scintNum] = {0}, PkT[2*scintNum] = {0}, Ped[2*scintNum] = {0}, RcT[2*scintNum] = {0}, RcTf[2*scintNum] = {0};
    int skipFlag = 0;

    for(int hit = 0; hit < nCry; hit++) { //ogni evt è impacchttato con un tag relativo alla hit

      int hitSide = iSide[hit], hitScint = iScint[hit], hitN = GetChan(hitSide, hitScint);

      if(isRun182) { //da eliminare
        if (hitSide == 0 && hitScint == 0) {}
        else if (hitSide == 1 && hitScint == 0) { continue; } 
        else if (hitSide == 0 && hitScint == 2) { hitSide = 1; hitScint = 0; hitN = GetChan(1,0);} 
        else {continue;}
      }

      IntQ[hitN] = Qval[hit]; 
      PkV[hitN] = Vmax[hit];
      PkT[hitN] = Tval[hit] - pkTOffset[hitN];
      Ped[hitN] = pedL[hit];
      double intQ = Qval[hit]; 
      double pkV= Vmax[hit];
      double pkT= Tval[hit] - pkTOffset[hitN]; 

      if (Selection.isSaturated(pkV)) {skipFlag = 1; continue;} //sto scartando i saturati < 1800 mV
      //if (intQ<150) {continue;}
      // if (hitScint !=0) {continue;}

      double tmin = pkT - 80;
      double tmax = pkT + 25;
      double norm = pkV; //norm = intQ/1.55; 

      TGraphErrors wgr =  TGraphErrors(400, ana::time, ana::wave[hit], ex, ey); //1) faccio un tgraph con l'onda //PARTIRE DA QUA

      //addon baseline pezzotta
        double blTmp{0}, brmsTmp{0};
        const int baseSam = 25;
        int pkbin = (int)pkT/digi_time; 
        for (int k = 0; k < baseSam; k++) {
          double v = wgr.GetPointY(pkbin - (int)100/digi_time - k);
          blTmp += v;
          brmsTmp += v*v;
        } 
        blTmp = blTmp/baseSam;
        brmsTmp = TMath::Sqrt(TMath::Abs(brmsTmp/baseSam - blTmp*blTmp)); 
        double EY = TMath::Sqrt(0.48*0.48/12 + brmsTmp*brmsTmp);
        HM.Fill1d("bLineRms", hitN, brmsTmp);
        HM.Fill1d("bLine", hitN, blTmp);  
        for (int i = 0; i < 400; i++) {
            wgr.SetPointX(i, wgr.GetPointX(i) - pkTOffset[hitN]);
            wgr.SetPointError(i, 0, EY);
        }
      //addon
      TSpline5 wsp = TSpline5("wsp", &wgr);

      auto spf = [&wsp](double *x, double *){ return wsp.Eval(x[0]); }; 
      TF1 fitf = TF1("fitf", spf, tmin, tmax, 0);
      norm = fitf.GetMaximum(tmin, tmax);
      double th = norm*constant_fraction_current;
      double rcT = fitf.GetX(th);
      TMarker tp = TMarker(rcT, wsp.Eval(rcT), 2);
      TLine t1 = TLine(tmin, 0, tmin, pkV);
      TLine t2 = TLine(tmax, 0, tmax, pkV);
      TLine t3 = TLine(tmin, norm, tmax, norm);

      if ( gRandom->Uniform(0, etp/50) < 1 ) { //plot per diagnostica
        samples_dir->cd();
        wgr.SetTitle(Form("pkV=%f CF=%f rt=%f", pkV, th, rcT));
        TCanvas cc(Form("tim_%lld", jentry)); cc.cd();
        wgr.SetLineWidth(1); wgr.SetMarkerStyle(20); wgr.SetMarkerSize(.4); wgr.SetMarkerColor(kBlue); wgr.Draw("AP"); 
        tp.SetMarkerSize(3); tp.SetMarkerColor(kRed); tp.Draw("same"); 
        t1.SetLineColor(kRed); t2.SetLineColor(kRed); t1.Draw("same"); t2.Draw("same");
        t3.SetLineColor(kOrange); t3.Draw("same");
        //wsp.SetLineColor(kOrange); wsp.Draw("same"); 
        fitf.SetLineColor(kSpring); fitf.Draw("same");
        cc.Write(); 
      }

      if ( !doFitTemplate && intQ > charge_min && intQ < charge_max ) { //genero il template
        for(int itime = (int)(ti_from/digi_time); itime < (int)(ti_to/digi_time); itime++) {  //faccio un for sull'onda campionata
          double wtime = (double)ana::time[itime] - rcT + 200 - pkTOffset[hitN]; //il tempo in ns viene allineato con il t_reco
          double wampl = ana::wave[hit][itime]/norm;  //ampiezza
          HM.Fill2d( "fuzzy", hitN, wtime, wampl); //fllo un th2 con  l'onda allineata in tempo  (opzione 1)
          int subn = digi_time * resam_fact;
          for(int k = 0; k < subn; k++){   //faccio il resampling (opzione 2)
            double subt = (double)ana::time[itime] - digi_time + (double)digi_time*((double)k/(double)subn) - 200;
            HM.Fill2d( "fuzzyResamp", hitN, subt - rcT - pkTOffset[hitN] + 500, wsp.Eval(subt)/norm );
          } 
        } 
      }

      double rcTf = -9999, chi2;
      if (doFitTemplate) { //fitto il template
        auto spfn = [&](Double_t *x, Double_t *par){ return par[0]*(spls[hitN]->Eval(x[0]-par[1]))+par[2]; };
        TF1 fitfn("fitfn", spfn, pkT - 120, pkT - 15, 3);  //70 18

        fitfn.SetParameter(0, pkV);
        //fitfn.SetParLimits(0, pkV*0.9, pkV*1.1);
        fitfn.SetParameter(1, pkT - 320); //320
        //fitfn.SetParLimits(1, pkT - 300, pkT - 280);
        fitfn.SetParameter(2,  0.);
        fitfn.SetParLimits(2, -5, 5);

        //TGraphErrors wgrn =  TGraphErrors(200, ana::time, ana::wave[hit], ex, ey);
        //TGraph wgrn =  TGraphErrors(200, ana::time, ana::wave[hit]);
        TGraphErrors wgrn = wgr;
        gStyle->SetOptFit(111);
        int fitr = wgrn.Fit( "fitfn", "REMQ" );        
        fitr = wgrn.Fit( "fitfn", "REMQ" );        
        if (fitr==-1) {cout<<"--------------------------------------------------------fitfailed:"<<fitr<<endl;}
        else {
          rcTf = fitfn.GetParameter(1) + templ_offs; 
          //rcTf = fitfn.GetX(th) + templ_offs; 
          chi2 = fitfn.GetChisquare()/fitfn.GetNDF();
          }

        // fitPar[0] = fitf.GetParameter(0);
        // fitPar[1] = fitf.GetParameter(1);
        // fitPar[2] = fitf.GetParameter(2);
        // fitErr[0] = fitf.GetParError(0);
        // fitErr[1] = fitf.GetParError(1);
        // fitErr[2] = fitf.GetParError(2);

        if ( gRandom->Uniform(0, etp/50) < 1 ) { //plot per diagnostica
        //if ( fitr == -1 ) {
          samples_dir->cd();
          wgrn.SetTitle("fit");
          TCanvas cc(Form("fit_%lld", jentry)); cc.cd();
          wgrn.SetLineWidth(1); wgrn.SetMarkerStyle(20); wgrn.SetMarkerSize(.4); wgrn.SetMarkerColor(kBlue); wgrn.Draw("P"); 
          fitfn.SetLineColor(kRed); fitfn.Draw("same"); 
          cc.Write();
        }
      }

      RcTf[hitN] = rcTf;
      RcT[hitN] = rcT;

      HM.Fill1d("recoTime", hitN, rcT);
      HM.Fill1d("recoTimeF", hitN, rcTf);
      HM.Fill1d("timeModBin", hitN, rcT/digi_time - ((int)rcT/digi_time));
      HM.Fill1d("timeModBinF", hitN, rcTf/digi_time - ((int)rcTf/digi_time));
      HM.Fill2d("slewing", hitN, intQ, rcT);
      HM.Fill2d("slewingFit", hitN, intQ, rcTf);
      HM.Fill1d("chargeRaw", hitN, intQ);
      HM.Fill1d("chi2", hitN, chi2);


    } //for hit
    
    if (skipFlag) {continue;}  

    for(int isc = 0; isc < scintNum; isc++) { 

      double tdiff = RcT[isc] - RcT[isc + scintNum];

      if (Selection.isChargeGood(IntQ, isc) &&
          Selection.isChargeGood(IntQ, isc + scintNum) &&
          PkT[isc] < 50 &&
          PkT[isc+scintNum] < 50 &&
          PkT[isc] > -50 &&
          PkT[isc+scintNum] > -50 &&
          tdiff < 50 &&
          tdiff > -50 &&
          IntQ[isc] > 200 &&
          IntQ[isc + scintNum] > 200
        )
      {
        HM.Fill1d("tsum", isc, RcT[isc] + RcT[isc + scintNum] );
        HM.Fill1d("tdiff", isc, tdiff );
        HM.Fill1d("tdifffit", isc, RcTf[isc] - RcTf[isc + scintNum]); 
        HM.Fill1d("chargeMip", isc, IntQ[isc]); HM.Fill1d("chargeMip", isc+scintNum, IntQ[isc+scintNum]);
      }   

    }

  } //for entries
  cout<<endl;
}

void Analysis::ProcessPlots() {}

void Analysis::Loop(){

  gErrorIgnoreLevel = 6001;
  cout<<endl<<endl<<"::::::::::::::::::::: CRT gen template :::::::::::::::::::::"<<endl<<endl;

  if (fChain == 0) return;

  cout<<"Creating histograms:"<<endl;
  HM.SetOutFile(outFile);
  HM.SetNamerFun(&NamerMatrix);
  createHistBoxes();
  cout<<"...done"<<endl<<endl;

  spline_dir = outFile->mkdir("splines");
  splineGr_dir = outFile->mkdir("splinesGr");
  splineDraw_dir = outFile->mkdir("profiles");
  templDraw_dir = outFile->mkdir("fuzzyTemplDraw");
  samples_dir = outFile->mkdir("randomSpecimens");
  preProcessing_dir = outFile->mkdir("preProcessing");

  Analysis::LoopOverEntries();
  HM.ProcessBoxes(); 
  Analysis::ProcessPlots();

  outFile->Close();

  cout<<endl<<endl<<"::::::::::::::::::::: done :::::::::::::::::::::"<<endl<<endl;
}

int main(int argc, char*argv[]) { 

  if (argc != 6) {
    printf("Usage: %s [infile_name] [outfile_name] [run_name] [calib_name] [CF]\n", argv[0]);
    exit(-1);
  }

  ofstream of("/home/ruben/Documents/frascati/CRT-analysis/data/reso.txt", std::ios_base::app);
  constant_fraction_current = atof(argv[5]);
  Analysis::Run(argv1, argv2, argv3, argv4, -1);
  of << constant_fraction_current << " " << resolution << " " << resolution_error << endl;
  of.close();
}



