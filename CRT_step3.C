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
#include <TF1Convolution.h>

#include "includes/Analysis.h"
#include "includes/langaus.h"
#include "includes/AnaPars.h"
#include "includes/CsvHandler.h"
#include "includes/MipSelection.h"
#include "includes/templ2charge.h"
#include "includes/HistManager.h"
#include "includes/NumberingHelper.h"

#include <TSpline.h>

using namespace std;

CsvHandler CSV;
MipSelection Selection; 
HistManager HM;

Long64_t nentries, nbytes, nb, ientry, jentry, jTrig_out;
double **peakTimeOffset, **timeOffset, **timeOffset_out, **zetaOffset, **zetaOffset_out, **chargeEqual;
double **chargeEqualErr, **chargeEqual_out, **chargeEqualErr_out, **barLen_out, **ped_out, **pedErr_out, **timeDiff_out, **timeDiffErr_out;
double tDiff = 0, zeta=0, pseudotDiff = 0, pseudoZeta = 0;
double *intQ, *pkV, *teQ, *teT, *teA, *teB, *teX2, *ped, *rcT;
list<double ** > arrayList = {&intQ, &pkV, &teQ, &teT, &teA, &teB, &teX2, &ped, &rcT};
void InitVectors() { for(double** &arr: arrayList) { *arr = new double[2*scintNum](); } }
int iSc_out; double_t Z_out; double_t Q_out[2], X2_out[2], T_out[2], pZ_out;
TTree *CRTs3;
#define tree_out_name "CRT"
#define tree_index_out "iTrig"
#define reindex_tree 0


 double flat(const double *x, const double *par){
   double ampl = par[0];
   double len = par[1];
   if (x[0] < -len || x[0] > len) return 0;
   else return ampl;
 } 



void chargeMip_proc(TH1* histObj, int histN, int& histSkipFlag) {   

  gStyle->SetOptFit(1); 

  TSpectrum s(2);
  s.Search(histObj, 1, "nodraw"); // a volte ne servono 2 in autotrigger
  double *pks = s.GetPositionX();
  double qpeak = *std::max_element(pks, pks + 2);

  //double min_tmp = histObj->GetXaxis()->GetXmin();
  //double max_tmp = histObj->GetXaxis()->GetXmax();
  //histObj->GetXaxis()->SetRangeUser(390,700);
  //double qpeak = histObj->GetBinCenter(histObj->GetMaximumBin());
  histObj->GetXaxis()->SetRangeUser(qFrom, qTo);
  double qmax = qpeak + 500, qmin = qpeak-80; float pk,sigma;

  TF1 l1 = TF1("l", "landau", qmin, qmax);                 l1.SetParameters(histObj->Integral()/2, 400, 100);                histObj->Fit(&l1, "RQ"); 
  pk = l1.GetMaximumX(); sigma = l1.GetParameter(2);
  TF1 l2 = TF1("l", "landau", pk-40, pk+3*sigma);          l2.SetParameters(l1.GetParameter(0), l1.GetParameter(1), sigma);  histObj->Fit(&l2, "RQ");
  pk = l2.GetParameter(1); sigma = l2.GetParameter(2); 
  TF1 l3 = TF1("l", "landau", pk-1*sigma, pk+4*sigma);     l3.SetParameters(l2.GetParameter(0), l2.GetParameter(1), sigma);  histObj->Fit(&l3, "RQ");
  pk = l3.GetParameter(1); sigma = l3.GetParameter(2);

  if (!centerMode) {qmin = pk-1.5*sigma; qmax = pk+4*sigma;} // in autotrigger 0.8 sigma a sx
  else {qmin = pk-2*sigma; qmax = pk+6*sigma;}

  TF1 f("f", langaufun, qmin, qmax, 4);
   //Fit parameters:
   //par[0]=Width (scale) parameter of Landau density
   //par[1]=Most Probable (MP, location) parameter of Landau density
   //par[2]=Total area (integral -inf to inf, normalization constant)
   //par[3]=Width (sigma) of convoluted Gaussian function

  f.SetParameters(sigma, pk, histObj->Integral(), 30);
  f.SetParName(0, "Landau Scale");
  f.SetParName(1, "Landau MPV");
  f.SetParName(2, "Integral");
  f.SetParName(3, "Gaussian Sigma");
  histObj->Fit("f", "RQ")   ;

  int iSd = (int)((histN+1)>scintNum), iSc = histN - (iSd==1)*scintNum; 

  chargeEqual_out[iSd][iSc] = f.GetParameter(1);
  chargeEqualErr_out[iSd][iSc] = f.GetParError(1);

}

double beta(double *xx, double *par){
  double x = xx[0], N = par[0], scale = par[1];
  if (x > scale * 2.28) return 0;
  return N * TMath::Sqrt( x*x + scale*2*0.511*x ) * (2.28*scale - x) * (2.28*scale-x) * (x + scale*0.511);
}

void srChargeMip_proc(TH1* histObj, int histN, int& histSkipFlag) {   

  gStyle->SetOptFit(1);

  histObj->GetXaxis()->SetRangeUser(50, 500);
  double peak = histObj->GetBinCenter(histObj->GetMaximumBin());
  TF1 *betafit = new TF1("betafit", beta, peak, 500, 2);
  betafit->SetParameter(1, 200);
  histObj->Fit("betafit", "QR");

  double guess_endpoint = betafit->GetParameter(1) * 2.28;
  betafit = new TF1("f", beta, guess_endpoint * 0.6 , guess_endpoint * 0.9, 2);
  betafit->SetParameter(1, guess_endpoint);
  histObj->Fit("f", "RQ")   ;

  guess_endpoint = betafit->GetParameter(1) * 2.28;
  betafit = new TF1("f", beta, guess_endpoint * 0.6 , guess_endpoint * 0.9, 2);
  betafit->SetParameter(1, guess_endpoint);
  histObj->Fit("f", "RQ");

  betafit->SetParName(0, "Normalization");
  betafit->SetParName(1, "Charge-Energy conversion");
  histObj->Fit("f", "RQ");

  int iSd = (int)((histN+1)>scintNum), iSc = histN - (iSd==1)*scintNum; 

  chargeEqual_out[iSd][iSc] = betafit->GetParameter(1);
  chargeEqualErr_out[iSd][iSc] = betafit->GetParError(1);

}


void timeMip_proc(TH1* histObj, int histN, int& histSkipFlag) {   

  gStyle->SetOptFit(1);

  double tpeak = histObj->GetBinCenter(histObj->GetMaximumBin());
  double tmax = tpeak + 40, tmin = tpeak - 40;
  //TF1 timeFit = TF1("logn", "[0]*ROOT::Math::lognormal_pdf(x,log([1]),log([2]))", tmin, tmax);  
  TF1 timeFit = TF1("g", "gaus", tmin, tmax); timeFit.SetParameter(1, tpeak); timeFit.SetParameter(2, 2);

  histObj->Fit(&timeFit, "RQ")   ;

  double mean = timeFit.GetParameter(1), sigma = timeFit.GetParameter(2);

  timeFit = TF1("g", "gaus", mean - sigma, mean + sigma); timeFit.SetParameter(1, mean); timeFit.SetParameter(2, sigma);

  histObj->Fit(&timeFit, "QR");

  int iSd = (int)((histN+1)>scintNum), iSc = histN - (iSd==1)*scintNum;

  if(TString(histObj->GetName()).Contains("timeDiffMip")){
    timeDiff_out[0][iSc] = timeFit.GetParameter(1);
    timeDiffErr_out[0][iSc] = timeFit.GetParError(1);
  }
  else timeOffset_out[iSd][iSc] = timeFit.GetParameter(1);

}

void pedMip_proc(TH1* histObj, int histN, int& histSkipFlag) {   

  gStyle->SetOptFit(1);

  double qpeakPed = histObj->GetBinCenter(histObj->GetMaximumBin());
  TF1 l = TF1("l", "gaus", qpeakPed-3, qpeakPed+3);
  histObj->Fit(&l, "RQ")   ;

  int iSd = (int)((histN+1)>scintNum), iSc = histN - (iSd==1)*scintNum; 

  ped_out[iSd][iSc] = l.GetParameter(2);
  pedErr_out[iSd][iSc] = l.GetParError(2);

}

void zetaMip_proc(TH1* histObj, int histN, int& histSkipFlag) {   

  gStyle->SetOptFit(1);

  TF1 zFit;
  if (centerMode) {
    zFit = TF1("l", "gaus", -20, 20);
    histObj->Fit(&zFit, "RQ")   ;
    histObj->Fit(&zFit, "RQ")   ;
  } else {
    zFit = TF1("l", "( TMath::TanH( (x-[0])/[1] ) - TMath::TanH( (x-[2]) /[1] ) ) * ([3])", -90, 90); // non fitta mai
    zFit.SetParLimits(0, -200, -20);
    zFit.SetParLimits(1, 0, 20);
    zFit.SetParLimits(2, -10, 200);
    zFit.SetParLimits(3, 1, 1000);
    histObj->Fit(&zFit, "QR0");
    histObj->Fit(&zFit, "QR");

    //TF1 *flat_f = new TF1("step_f", "[0]*(x<[1])*(x>-[1])", -150, 150); ///////////anche così
   TF1 *flat_f = new TF1("flat_f", flat, -100, 100, 2);
   TF1 *gauss_f = new TF1("gauss_f", "gaus", -100, 100);
   TF1Convolution *f_conv = new TF1Convolution(flat_f, gauss_f, -150, 150, true);
   f_conv->SetNofPointsFFT(1000);
   //cout << "ZOFF: " << (zFit.GetParameter(2) + zFit.GetParameter(0))/2 << endl;
   //TF1 *f = new TF1("f", *f_conv, (zFit.GetParameter(2) + zFit.GetParameter(0))/2 - 90, (zFit.GetParameter(2) + zFit.GetParameter(0))/2 + 90, f_conv->GetNpar());
   TF1 *f = new TF1("f", *f_conv, -90, 90, f_conv->GetNpar());
   f->SetParameters(200., 76, 6, -66, 6);
   f->SetParName(0, "Ampl_flat");
   f->SetParName(1, "Len");
   f->SetParName(2, "Ampl_gauss");
   f->SetParName(3, "Mean_reso");
   f->SetParName(4, "Reso");

//   TF1 *flat = new TF1("f", "pol0", -80, 80);

   histObj->Fit("f", "RQ");

*/
   barLen_out[0][histN] = 0;//f->GetParameter(1);
  }

  zetaOffset[0][histN] = centerMode ? zFit.GetParameter(1) : (zFit.GetParameter(2) + zFit.GetParameter(0))/2 ;

}

void qSharing_proc(TH1* histObj, int histN, int& histSkipFlag) {   
  if (histN == 0 || histN == scintNum - 1) {histSkipFlag=1;}
}



void createHistBoxes() {

    HM.AddHistBox("th1f", 2*scintNum, "chargeRaw",      "Raw charges",      "charge", "pC",    qBins, 20, qTo);
    HM.AddHistBox("th1f", 2*scintNum, "chargeMip",      "MIP charges",      "charge", "pC",    qBins, qFrom, qTo, strontiumMode ? &srChargeMip_proc : &chargeMip_proc);

    HM.AddHistBox("th1f", scintNum,   "chargeRawPerScint",  "MIP charges",      "charge", "pC",    qBins, qFrom, qTo, &chargeMip_proc);

    HM.AddHistBox("th1f", 2*scintNum, "chargeTeMip",    "MIP template q",   "charge", "pC",    qBins, qFrom, qTo);
    HM.AddHistBox("th1f", 2*scintNum, "pedMip",         "Pedestal",         "charge", "pC",    100, -10, 10, &pedMip_proc);
    HM.AddHistBox("th1f", 2*scintNum, "voltPeak",       "Wave peak",        "ampl", "V",       100, 0, 2000);
    HM.AddHistBox("th1f", 2*scintNum, "timeMip",        "MIP times",        "time", "ns",      100, -30, 30, &timeMip_proc);
    HM.AddHistBox("th1f", scintNum,   "timeDiffMip",    "MIP time difference",        "time_difference", "ns",      1000, -50, 50, &timeMip_proc);
    HM.AddHistBox("th1f", scintNum,   "zetaMip",        "MIP zetas",        "zeta", "cm",      320, -scintL, scintL, &zetaMip_proc, &NamerArray);
    HM.AddHistBox("th1f", scintNum,   "pseudoZeta",        "MIP zetas",        "zeta", "cm",      320, -scintL, scintL, &zetaMip_proc, &NamerArray);
    HM.AddHistBox("th2f", 2*scintNum, "q_chi2",         "MIP q vs chi2",    "charge", "pC", "chi2", "",           qBins/2, qFrom, qTo, 100, 0, 40);
    HM.AddHistBox("th2f", 2*scintNum, "zeta_q",         "MIP q vs Z",       "zeta", "cm", "charge", "pC",         160, -scintL, scintL, qBins/2, qFrom, qTo);
    HM.AddHistBox("th2f", scintNum,   "qSharing",       "Sharing",          "Q_i", "pC", "Q_neighbours", "pC",    50, qFrom, qTo, 50, qFrom, qTo, &qSharing_proc);

}


void fill_raw(int hitN) {

  HM.Fill1d("chargeRaw", hitN, intQ[hitN]); 
  HM.Fill1d("chargeRaw", hitN, intQ[hitN]);
  HM.Fill1d("voltPeak", hitN, pkV[hitN]); 
  HM.Fill1d("voltPeak", hitN, pkV[hitN]);

}

void fill_mip(int iScHit) {

  HM.Fill1d("chargeMip", iScHit, intQ[iScHit]);
  HM.Fill1d("chargeMip", iScHit+scintNum, intQ[iScHit+scintNum]);
  HM.Fill1d("chargeTeMip", iScHit, teQ[iScHit]);
  HM.Fill1d("chargeTeMip", iScHit+scintNum, teQ[iScHit+scintNum]);
  HM.Fill1d("timeMip", iScHit, teT[iScHit]);
  HM.Fill1d("timeMip", iScHit+scintNum, teT[iScHit+scintNum]);
  HM.Fill1d("timeDiffMip", iScHit, teT[iScHit] - teT[iScHit+scintNum]);
  HM.Fill1d("pedMip", iScHit, ped[iScHit]);
  HM.Fill1d("pedMip", iScHit+scintNum, ped[iScHit+scintNum]);
  HM.Fill1d("zetaMip", iScHit, zeta);

  HM.Fill2d("zeta_q", iScHit, zeta, intQ[iScHit]);
  HM.Fill2d("zeta_q", iScHit+scintNum, zeta, intQ[iScHit+scintNum]);

  HM.Fill1d("pseudoZeta", iScHit, pseudoZeta);

  jTrig_out = jentry; iSc_out = iScHit;
  Z_out = zeta;
  pZ_out = pseudoZeta;
  Q_out[0] = intQ[iScHit]; Q_out[1] = intQ[iScHit+scintNum];
  T_out[0] = teT[iScHit]; T_out[1] = teT[iScHit+scintNum];
  X2_out[0] = teX2[iScHit]; X2_out[1] = teX2[iScHit+scintNum];
  CRTs3->Fill();
}



void Analysis::ProcessPlots() {

  TDirectory* calib_dir = outFile->mkdir("calibration");
  calib_dir->cd();

  //ChEq
    TGraphErrors *chEqGraph = new TGraphErrors(2*scintNum); chEqGraph->SetTitle("Charge equal");
    double qMean = 0;
    for (int k = 0; k < scintNum; k++) { 
      chEqGraph->SetPoint(k, (float)k+1-0.07, chargeEqual_out[0][k]);
      chEqGraph->SetPoint(k+scintNum,  (float)k+1+0.07, chargeEqual_out[1][k]);
      chEqGraph->SetPointError(k, 0, chargeEqualErr_out[0][k]);
      chEqGraph->SetPointError(k+scintNum, 0, chargeEqualErr_out[0][k]);
      qMean += chargeEqual_out[0][k] + chargeEqual_out[1][k];
    }
    qMean = qMean/(2*scintNum);
    TCanvas *eq_c = new TCanvas("chargeEqual", "chargeEqual"); eq_c->cd(); 
    chEqGraph->SetLineColor(kBlue); chEqGraph->SetMarkerColor(kBlue); chEqGraph->SetMarkerSize(1.4); chEqGraph->SetMarkerStyle(25); 
    chEqGraph->GetXaxis()->SetRangeUser(0, scintNum+1); chEqGraph->Draw("AP"); 
    TLine *line = new TLine(0.5, qMean, scintNum+0.5, qMean); line->SetLineColor(kRed); line->Draw("same");
    line = new TLine(0.5, 1.1*qMean, scintNum+0.5, 1.1*qMean); line->Draw("same");
    line = new TLine(0.5, 0.9*qMean, scintNum+0.5, 0.9*qMean); line->Draw("same");
    eq_c->Write("chargeEqual");
  //ChEq

  //ChEqOld
    chEqGraph = new TGraphErrors(2*scintNum); chEqGraph->SetTitle("chargeEqual_s3p1 (loaded from calib)");
    qMean = 0;
    for (int k = 0; k < scintNum; k++) { 
      chEqGraph->SetPoint(k, (float)k+1-0.07, chargeEqual[0][k]);
      chEqGraph->SetPoint(k+scintNum,  (float)k+1+0.07, chargeEqual[1][k]);
      qMean += chargeEqual[0][k] + chargeEqual[1][k];
    }
    qMean = qMean/(2*scintNum);
    eq_c = new TCanvas("chargeEqual_s3p1", "chargeEqual_s3p1"); eq_c->cd(); 
    chEqGraph->SetLineColor(kBlue); chEqGraph->SetMarkerColor(kBlue); chEqGraph->SetMarkerSize(1.4); chEqGraph->SetMarkerStyle(25); 
    chEqGraph->GetXaxis()->SetRangeUser(0, scintNum+1); chEqGraph->Draw("AP");  
    line = new TLine(0.5, qMean, scintNum+0.5, qMean); line->SetLineColor(kRed); line->Draw("same");
    line = new TLine(0.5, 1.1*qMean, scintNum+0.5, 1.1*qMean); line->Draw("same");
    line = new TLine(0.5, 0.9*qMean, scintNum+0.5, 0.9*qMean); line->Draw("same");
    eq_c->Write("chargeEqual_s3p1");
  //ChEqOld

  //TimeOff
    TGraphErrors *tim = new TGraphErrors(2*scintNum);  tim->SetTitle("timeOffset");
    for (int k = 0; k < scintNum; k++) { 
      tim->SetPoint(k, (float)k+1-0.07, timeOffset_out[0][k]);
      tim->SetPoint(k+scintNum,  (float)k+1+0.07, timeOffset_out[1][k]); 
    }
    TCanvas *tim_c = new TCanvas("timeOff", "timeOff"); tim_c->cd(); 
    tim->SetLineColor(kBlue); tim->SetMarkerColor(kBlue); tim->SetMarkerSize(1.4); tim->SetMarkerStyle(25); 
    tim->GetXaxis()->SetRangeUser(0, scintNum+1);
    tim->Draw("AP"); tim_c->Write("timeOff");
  //TimeOff

  //TimeOffOld
    tim = new TGraphErrors(2*scintNum);  tim->SetTitle("timeOffset_s3p1 (loaded from calib)");
    for (int k = 0; k < scintNum; k++) { 
      tim->SetPoint(k, (float)k+1-0.07, timeOffset[0][k]);
      tim->SetPoint(k+scintNum,  (float)k+1+0.07, timeOffset[1][k]); 
    }
    tim_c = new TCanvas("timeOff", "timeOff"); tim_c->cd(); 
    tim->SetLineColor(kBlue); tim->SetMarkerColor(kBlue); tim->SetMarkerSize(1.4); tim->SetMarkerStyle(25); 
    tim->GetXaxis()->SetRangeUser(0, scintNum+1);
    tim->Draw("AP"); tim_c->Write("timeOff_s3p1");
  //TimeOffOld

  //zetaOff
    TGraphErrors *zOffGraph = new TGraphErrors(scintNum); zOffGraph->SetTitle("Zeta offset");
    for (int k = 0; k < scintNum; k++) { 
      zOffGraph->SetPoint(k, k+1, zetaOffset_out[0][k]);
    }
    line = new TLine(0.5, 0, scintNum+0.5, 0); 
    TCanvas *zet_c = new TCanvas("zetaOffset", "zetaOffset"); zet_c->cd(); 
    zOffGraph->SetLineColor(kBlue); zOffGraph->SetMarkerColor(kBlue); zOffGraph->SetMarkerSize(1.4); zOffGraph->SetMarkerStyle(25); 
    zOffGraph->GetXaxis()->SetRangeUser(0, scintNum+1);
    line->SetLineColor(kRed);
    zOffGraph->Draw("AP"); line->Draw("same"); zet_c->Write("zetaOffset");
  //zetaOff

  //zetaOffOld
    zOffGraph = new TGraphErrors(scintNum); zOffGraph->SetTitle("zetaOffset_s3p2 (loaded from calib)");
    for (int k = 0; k < scintNum; k++) {
      zOffGraph->SetPoint(k, k+1, zetaOffset[0][k]);
    }
    line = new TLine(0.5, 0, scintNum+0.5, 0); 
    zet_c = new TCanvas("zetaOffset", "zetaOffset"); zet_c->cd(); 
    zOffGraph->SetLineColor(kBlue); zOffGraph->SetMarkerColor(kBlue); zOffGraph->SetMarkerSize(1.4); zOffGraph->SetMarkerStyle(25); 
    zOffGraph->GetXaxis()->SetRangeUser(0, scintNum+1);
    line->SetLineColor(kRed);
    zOffGraph->Draw("AP"); line->Draw("same"); zet_c->Write("zetaOffset_s3p2");
  //zetaOffOld

  //pedGraph
    TGraphErrors *pedGraph = new TGraphErrors(2*scintNum); pedGraph->SetTitle("Pedestal");
    for (int k = 0; k < scintNum; k++) {
      pedGraph->SetPoint(k, (float)k+1-0.07, ped_out[0][k]);
      pedGraph->SetPointError(k, 0, pedErr_out[0][k]);
      pedGraph->SetPoint(k+scintNum,  (float)k+1+0.07, ped_out[1][k]);
      pedGraph->SetPointError(k+scintNum, 0, pedErr_out[1][k]);
    }
    pedGraph->Draw("AP");
    pedGraph->Write();
  //pedGraph

  //aggiungere plot complessivi zeta con parametri del fit fatto con convoluzione
}


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

      if ( (modulSel=="B" && iMod[hit]==0) || (modulSel=="T" && iMod[hit]==1) ) {continue;} // TOP = 0, BTM = 1

      int hitSide=iSide[hit], hitScint = iScint[hit], hitN = hitSide*scintNum + hitScint;


      if(isRun182) { //da eliminare
        if (hitSide == 0 && hitScint == 0) {}
        else if (hitSide == 1 && hitScint == 0) { continue; } 
        else if (hitSide == 0 && hitScint == 2) { hitSide = 1; hitScint = 0; hitN = GetChan(1,0);} 
        else {continue;}
      }

      double chCal = enableOfflineEq ? chEqReference/chargeEqual[hitSide][hitScint] : 1;
      //chCal = 1.15;

      intQ[hitN] = Qval[hit]*chCal;
      pkV[hitN] = Vmax[hit];
      ped[hitN] = pedL[hit];
      teA[hitN] = templFit[hit][0];
      teQ[hitN] = templ2charge(teA[hitN])*chCal;
      teB[hitN] = templFit[hit][2];
      teT[hitN] = templTime[hit] - timeOffset[hitSide][hitScint];
      teX2[hitN] = templChi2[hit]; 

      if (Selection.isSaturated(pkV[hitN])) {skipFlag = 1; continue;}

      TGraph wgr =  TGraph(400, time, wave[hit]);

      //addon
      TSpline5 wsp = TSpline5("wsp", &wgr);

      auto spf = [&wsp](double *x, double *){ return wsp.Eval(x[0]); };
      double offset = peakTimeOffset[hitSide][hitScint];
      TF1 fitf = TF1("fitf", spf, offset - 100, offset + 100, 0);
      double th = fitf.GetMaximum(offset - 100, offset + 100)*0.15;
      rcT[hitN] = fitf.GetX(th) - offset;

      fill_raw(hitN);
    }

    if (skipFlag == 1) {continue;}

    for(int isc = 0; isc < scintNum; isc++) {

      HM.Fill1d("chargeRawPerScint", isc, intQ[isc] + intQ[isc+scintNum]);

      if (Selection.hitPrecheck(isc, iScint, nCry) && Selection.isChargeGood(intQ, isc) ) {

          if( (isc+1)%scintNum > 1 ) { HM.Fill2d("qSharing", isc, intQ[isc] + intQ[isc+scintNum],  intQ[isc-1] + intQ[isc+scintNum-1] + intQ[isc+1] + intQ[isc+scintNum+1]); }
        //Fill
        if ( !Selection.isShared(intQ, isc) ) { iScHit = isc; m++; }

      }
    }

    if ( m != 1 ) {continue;}

    if ( !Selection.isTimeGood(teT[iScHit])) {continue;}

    tDiff = teT[iScHit] - teT[scintNum+iScHit]; 

    zeta = tDiff*scintVp/2-zetaOffset[0][iScHit];

    pseudoZeta = (rcT[iScHit] - rcT[iScHit + scintNum]) * scintVp / 2;

    //if ( !Selection.isZetaGood(zeta) ) {continue;}
    //if ( 1 || !Selection.mipCutG(intQ, zeta, iScHit) ) {continue;} // ritorna già 1 se cutG non è enabled, non serve 1 || 

    fill_mip(iScHit);    
  }

  cout<<endl;
}


void Analysis::Loop(){

  cout<<endl<<endl<<"::::::::::::::::::::: CRT analysis step 3 :::::::::::::::::::::"<<endl<<endl;

  if (fChain == 0) return;

  cout<<"Creating histograms :"<<endl;
  HM.SetOutFile(outFile);
  HM.SetNamerFun(&NamerMatrix);
  createHistBoxes();
  cout<<"...done"<<endl<<endl;

  cout<<"Loading cutGs :";
  Selection.loadCutG();
  cout<<"...done"<<endl<<endl;

  cout<<"Retrieving calibration data from [" + lutPrefix3p + "] :"<<endl;
  timeOffset  = CSV.InitMatrix(2, scintNum); timeOffset_out = CSV.InitMatrix(2, scintNum);
  peakTimeOffset  = CSV.InitMatrix(2, scintNum);
  chargeEqual = CSV.InitMatrix(2, scintNum); chargeEqual_out = CSV.InitMatrix(2, scintNum);
  chargeEqualErr = CSV.InitMatrix(2, scintNum); chargeEqualErr_out = CSV.InitMatrix(2, scintNum);
  zetaOffset  = CSV.InitMatrix(1, scintNum); zetaOffset_out  = CSV.InitMatrix(1, scintNum);
  timeDiff_out  = CSV.InitMatrix(1, scintNum);
  timeDiffErr_out  = CSV.InitMatrix(1, scintNum);
  barLen_out = CSV.InitMatrix(1, scintNum); 

  ped_out = CSV.InitMatrix(2, scintNum);
  pedErr_out = CSV.InitMatrix(2, scintNum);

  CSV.Read(CSV.GetFirstFile(lutPrefix3p + calName + lutChEqName + "*"),      ',', chargeEqual, 2, scintNum);
  CSV.Read(CSV.GetFirstFile(lutPrefix3p + calName + lutTimeOffsName + "*"),  ',', timeOffset,  2, scintNum);
  CSV.Read(CSV.GetFirstFile(lutPrefix3p + calName + lutZetaOffName + "*"),   ',', zetaOffset,  1, scintNum); 
  CSV.Read(CSV.GetFirstFile(lutPrefix3p + calName + lutPeakTimeOffsName + "*"),   ',', peakTimeOffset,  2, scintNum);

  cout<<"...done"<<endl<<endl;

  outFile->cd();
  CRTs3 = new TTree(tree_out_name, tree_out_name);          
  CRTs3->SetAutoSave(1000);
  CRTs3->Branch("iTrig",   &jTrig_out,  "iTrig/L");
  CRTs3->Branch("iSc",     &iSc_out,    "iSc/I");
  CRTs3->Branch("Z",       &Z_out,      "Z/D");
  CRTs3->Branch("pZ",      &pZ_out,     "pZ/D");
  CRTs3->Branch("Q",       &Q_out,      "Q[2]/D");
  CRTs3->Branch("T",       &T_out,      "T[2]/D");
  CRTs3->Branch("X2",      &X2_out,     "X2[2]/D");

  Analysis::LoopOverEntries();
  HM.ProcessBoxes();
  Analysis::ProcessPlots();

  cout<<endl<<"Writing data to [" + lutPrefix3 + "] :"<<endl;
  CSV.Write(lutPrefix3 + runName + lutChEqName + ".csv", ',', chargeEqual_out, 2, scintNum, 4);
  CSV.Write(lutPrefix3 + runName + lutChEqErrName + ".csv", ',', chargeEqualErr_out, 2, scintNum, 4);
  CSV.Write(lutPrefix3 + runName + lutBarLenName + ".csv", ',', barLen_out, 1, scintNum, 4);
  CSV.Write(lutPrefix3 + runName + lutTimeDiffName + ".csv", ',', timeDiff_out, 1, scintNum, 4);
  CSV.Write(lutPrefix3 + runName + lutTimeDiffErrName + ".csv", ',', timeDiffErr_out, 1, scintNum, 4);

  cout<<"...done"<<endl;

  outFile->cd();
  if (reindex_tree == 1) { CRTs3->BuildIndex(tree_index_out); }
  CRTs3->Write();

  HM.CloseOutFile();
  cout<<endl<<endl<<"::::::::::::::::::::: analysis done :::::::::::::::::::::"<<endl<<endl;
}



#define inFile_f "../../data/step2/%s_s2.root"
#define outFile_f "../../data/step3/%s_s3.root"

void CRT_step3(TString run_name, TString calib_name, TString mod_select = "") {

  gErrorIgnoreLevel = kFatal; //kPrint, kInfo, kWarning, kError, kBreak, kSysError, kFatal

  if (mod_select != "T" && mod_select != "B" && mod_select != "") {cout<<"modSel must be empty, 'T' or 'B' !!"<<endl; return;}
  //opzione "" per legacy

  TString inFileName = Form(inFile_f, run_name.Data());
  TString runName = mod_select == "" ? run_name : run_name + "_" + mod_select;
  TString calibName = mod_select == "" ? calib_name : calib_name + "_" + mod_select;

  TString outFileName = Form(outFile_f, runName.Data());


  TFile *fileOut = new TFile(outFileName, "RECREATE");

  cout<<endl<<"------------> Launching step3:"<<endl;
  cout<<"----> runName : "<<runName<<endl;
  cout<<"----> Input file: "<<inFileName<<endl;
  cout<<"----> Output file: "<<outFileName<<endl;
  cout<<"----> Calibration files: "<<calibName<<endl;
  cout<<"----> Module selector: "<<mod_select<<endl<<endl;


  Analysis *a = new Analysis(inFileName, fileOut, runName, calibName, mod_select);

  a->Loop();

}
