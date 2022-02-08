// ancora non programmato



#define analysis_CRT_cxx
#include "analysis_CRT.h"
//#include <TH2.h>
//#include <TStyle.h>
//#include <TCanvas.h>
//#include <string>
//#include <stdlib.h>
#include <iostream> // Already in analysis_mod.h
#include <TButton.h>
#include <TMath.h>
#include "TF1.h"
#include "logn.h"
#include <TSpline.h>
#include <algorithm>
#include <TCanvas.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TSystem.h>
#include <unistd.h>
#include <stdio.h>
#include <TApplication.h>
#include <TAxis.h>

TSpline3 *sp3;
TSpline5 ***sp5_arr;

using namespace std;

int plot = 1;
int templ = 2; // 1 è quello vero - 2 è il mio accrocco
int go_ahead_flag = 0;

void go_ahead(){
  go_ahead_flag = 1;
}

TCanvas *c;
TPad *plots_pad;

//=============================================================================
// Fit function for template
//=============================================================================

Double_t spline_fit(Double_t *x, Double_t *par)
{
  double f1;
  f1=par[0]*sp3->Eval(x[0]-par[1])+par[2];
  return f1;  
}

int sideTmp, scintTmp;
Double_t spline5_fit(Double_t *x, Double_t *par)
{
  //if(plot) cout << Form("Fitting with template2 for side %i - scint %i", sideTmp, scintTmp) << endl;
  double f1;
  f1=par[0]*sp5_arr[sideTmp][scintTmp]->Eval(x[0]-par[1])+par[2];
  return f1;  
}

//=============================================================================
// User variables
//=============================================================================

int   const Nboard = 2;
int   const Nchan  = 8;

int         Itmp;
//int         ItmpLaser;
float       Qtmp, Ttmp, Btmp, Brmstmp, Vtmp;
float       ped1tmp, ped2tmp, lognTimetmp=0, lognChi2tmp=0;
Double_t    wavetmp[maxNsample];
int         BoaTmp, ChaTmp, DaqTmp, ScintTmp, SideTmp;

int evtmp;

float       Time, Chi2;
float       fitPar[3]={0.}, fitErr[3]={0.}, fitTmea=0., fitChi2=0.;
float       lognfitPar[4]={0.}, lognfitErr[4]={0.};

ofstream fit_log("fitlog");

//=============================================================================
// Loop function
//=============================================================================

//void analysis_CRT::Loop()
void analysis_CRT::Loop(TString OutputFile, int evflag)
{
//   In a ROOT session, you can do:
//      Root > .L analysis_CRT.C
//      Root > analysis_CRT t
//      Root > t.GetEntry(12); // Fill t data members with entry number 12
//      Root > t.Show();       // Show values of entry 12
//      Root > t.Show(16);     // Read and show values of entry 16
//      Root > t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   int skipfit = 0;
   //  float Qcut=0;
   if( evflag == 0) { // study for cosmics;
     Qcut   = 2;// prima era 10;
     skipfit =0;
   }else if( evflag==1){
     Qcut = -100.;
     skipfit=1;
   }else{
     cout << "Fatal error" << endl;
     return;
   }

   // **************************************************************************
   // Input file: template to fit timing
   // **************************************************************************
   TGraphErrors* gt;
   if (templ==1){
     //TFile* f0 = new TFile("splines1_CRT.root"); //di eleonora - su pcdiociaiuti
     TFile* f0 = new TFile("splines_pisani_182.root");
     gt = (TGraphErrors*)f0->Get("gtempl");

     sp3 = new TSpline3("sp3", gt, 0, 0., 800.);
   }
   else{
     TFile* f1 = new TFile("waves_splines_run205_coninterp.root");

     sp5_arr = new TSpline5**[2];

     for(int isd=0; isd<2; isd++){
        sp5_arr[isd] = new TSpline5*[8];

        for(int isc=0; isc<8; isc++){
          sp5_arr[isd][isc] = (TSpline5*)f1->Get(Form("spline_interp_wave_%i_%i", isd, 3)); // per il 183, che era scint3 con lato 0 in 0, 0, e lato 1 in 0, 2
        }
     }
   }
   // **************************************************************************
   // Output file
   // **************************************************************************

   cout << "Selected output file: " << OutputFile << endl;
   TFile *outFile = new TFile(OutputFile,"recreate");
   TTree *CRT = new TTree("CRT","CRT");
   BookOutput(CRT); // Book histos and ntuple

   // **************************************************************************
   // Loop on entries
   // **************************************************************************
   Long64_t nentries = 100000; //fChain->GetEntriesFast();

   if(evflag==1) nentries = 5000;

   cout<<"nentries: "<<nentries<<endl;

   Long64_t nbytes = 0, nb = 0;

   if(plot){
    c = new TCanvas("c", "c");
    plots_pad = new TPad("pad1", "pad1",0.0, 0.1, 1.0, 1.0);
    TPad *pad1 = new TPad("pad2", "pad2",0.0, 0, 1.0, 0.1);
    plots_pad->Draw();
    pad1->Draw();
    pad1->cd();
    TButton *but1 = new TButton("Next", "go_ahead()", 0.45, 0.15, 0.55, 0.95);
    but1->Draw();
    plots_pad->cd();
   }

   for (Long64_t jentry=0; jentry<nentries;jentry++) {

      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;

      if( jentry%1000==0 ) cout << "Number of processed events: " << jentry << endl;
      
      ntrig_out = ntrig;
      evnum_out = evnum;
      nsample   = num;
      evtmp = evnum;
      std::copy(time, time+nsample, time_out);

      int cryTot = 0;

      for( int Iboard=0; Iboard<Nboard; Iboard++ )	{
        for( int Ichan=0; Ichan<Nchan; Ichan++ )	   {


          GetValues(Iboard,Ichan,evflag);
          // Get values for a given board/channel

          if( Qtmp>5 ){

            iDAQ[cryTot]   = DaqTmp;
            iScint[cryTot] = ScintTmp;
            iSide[cryTot]  = SideTmp;
            Qval[cryTot]   = Qtmp;
            Tval[cryTot]   = Ttmp;
            pedL[cryTot]   = ped1tmp;
            pedH[cryTot]   = ped2tmp;
            iMax[cryTot]   = Itmp;
            Vmax[cryTot]   = Vtmp;
            std::copy(wavetmp,wavetmp+nsample,wave[cryTot]);
            bline[cryTot]  = Btmp;

            // Timing from template fit
            if( skipfit==0){
              if(plot){
                if (templ == 1) getTemplateFit(gt);
                else getTemplate2Fit();
              }
              else {
                if (templ == 1) getTemplateFit(gt);
                else getTemplate2Fit(); 
              }
              //getTlogn();
            }

            templTime[cryTot] = fitTmea;
            templChi2[cryTot] = fitChi2;
            std::copy(fitPar,fitPar+3,templFit[cryTot]);
            std::copy(fitErr,fitErr+3,templErr[cryTot]);
      
            lognTime[cryTot] = lognTimetmp;
            lognChi2[cryTot] = lognChi2tmp;
            std::copy(lognfitPar,lognfitPar+4,lognFit[cryTot]);
            std::copy(lognfitErr,lognfitErr+4,lognErr[cryTot]);

            cryTot++;

          }
          
        } // Loop on Iboard
      } // Loop on Ichan
        
      nCry = cryTot;
      
      CRT->Fill();
          
    } // end loop on events

    cout << "Closing output file" << endl;

    outFile->cd();
    CRT->Write();
    outFile->Close();
}

//=============================================================================
// Ntuple booking
//=============================================================================
void analysis_CRT::BookOutput(TTree *CRT)
{
  // Create ROOT tree output

  CRT->SetAutoSave(1000);

  CRT->Branch("ntrig"     ,&ntrig_out, "ntrig/I");
  CRT->Branch("evnum"     ,&evnum_out, "evnum/I");
  CRT->Branch("nsample"   ,&nsample,   "nsample/I");       // function of nCry?
  CRT->Branch("time"      ,&time_out,  "time[nsample]/D"); // function of nCry?
  CRT->Branch("nCry"      ,&nCry,      "nCry/I");
  CRT->Branch("iDAQ"      ,&iDAQ,      "iDAQ[nCry]/I");
  CRT->Branch("iScint"    ,&iScint,    "iScint[nCry]/I");
  CRT->Branch("iSide"     ,&iSide,     "iSide[nCry]/I");
  CRT->Branch("iMax"      ,&iMax,      "iMax[nCry]/I");
  CRT->Branch("Vmax"      ,&Vmax,      "Vmax[nCry]/D");
  CRT->Branch("Qval"      ,&Qval,      "Qval[nCry]/D");
  CRT->Branch("Tval"      ,&Tval,      "Tval[nCry]/D");
  CRT->Branch("pedL"      ,&pedL,      "pedL[nCry]/D");
  CRT->Branch("pedH"      ,&pedH,      "pedH[nCry]/D");
  CRT->Branch("wave"      ,&wave,      "wave[nCry][1024]/D");
  CRT->Branch("bline"     ,&bline,     "bline[nCry]/D");
  CRT->Branch("lognTime"  ,&lognTime,  "lognTime[nCry]/D");
  CRT->Branch("lognChi2"  ,&lognChi2,  "lognChi2[nCry]/D");
  CRT->Branch("lognFit"   ,&lognFit,   "lognFit[nCry][4]/D");
  CRT->Branch("lognErr"   ,&lognErr,   "lognErr[nCry][4]/D");
  CRT->Branch("templTime" ,&templTime, "templTime[nCry]/D");
  CRT->Branch("templChi2" ,&templChi2, "templChi2[nCry]/D");
  CRT->Branch("templFit"  ,&templFit,  "templFit[nCry][3]/D");
  CRT->Branch("templErr"  ,&templErr,  "templErr[nCry][3]/D");

}

//=============================================================================
// Get values for a given board/channel
//=============================================================================
void analysis_CRT::GetValues(int Iboard, int Ichan, int evflag)
{

  int   nVal;
  float TAmax, TWmin, TWmax;
  float Time, sum, sumsq;
  int   QbinMin, QbinMax;
  float const camp = 0.25;
  float dt = 1/camp;
  float const Qpeak_min      =  60.; //  -60 ns from Tave
  float const Qpeak_max      = 190.; // +190 ns from Tave

  // Read map [do it once, not for all channels!!!]

  FILE *mapFile;
  mapFile = fopen("CRT_map.dat", "r");
  for ( int Iloop=0; Iloop<Nboard*Nchan; Iloop++ ){
    fscanf(mapFile, "%d %d %d %d ", &BoaTmp, &ChaTmp, &ScintTmp, &SideTmp );
    //cout << Iboard << " " << Ichan << " " << BoaTmp << " " << ChaTmp << endl;
    if( BoaTmp==Iboard && ChaTmp==Ichan ) break;
  }
  fclose(mapFile);
  
  DaqTmp = 8*BoaTmp + ChaTmp;

  // cryTmp still to be defined

  if( Iboard==0 ){
    if( Ichan==0 ){
      Ttmp    = b0_tave0;
      ped1tmp = b0_ped0;
      ped2tmp = b0_pedh0;
      Itmp    = b0_imax0;
      std::copy(b0_onda0,b0_onda0+nsample,wavetmp);
    }
    else if( Ichan==1 ){
      Ttmp    = b0_tave1;
      ped1tmp = b0_ped1;
      ped2tmp = b0_pedh1;
      Itmp    = b0_imax1;
      std::copy(b0_onda1,b0_onda1+nsample,wavetmp);
    }
    else if( Ichan==2 ){
      Ttmp    = b0_tave2;
      ped1tmp = b0_ped2;
      ped2tmp = b0_pedh2;
      Itmp    = b0_imax2;
      std::copy(b0_onda2,b0_onda2+nsample,wavetmp);
    }
    else if( Ichan==3 ){
      Ttmp    = b0_tave3;
      ped1tmp = b0_ped3;
      ped2tmp = b0_pedh3;
      Itmp    = b0_imax3;
      std::copy(b0_onda3,b0_onda3+nsample,wavetmp);
    }
    else if( Ichan==4 ){
      Ttmp    = b0_tave4;
      ped1tmp = b0_ped4;
      ped2tmp = b0_pedh4;
      Itmp    = b0_imax4;
      std::copy(b0_onda4,b0_onda4+nsample,wavetmp);
    }
    else if( Ichan==5 ){
      Ttmp    = b0_tave5;
      ped1tmp = b0_ped5;
      ped2tmp = b0_pedh5;
      Itmp    = b0_imax5;
      std::copy(b0_onda5,b0_onda5+nsample,wavetmp);
    }
    else if( Ichan==6 ){
      Ttmp    = b0_tave6;
      ped1tmp = b0_ped6;
      ped2tmp = b0_pedh6;
      Itmp    = b0_imax6;
      std::copy(b0_onda6,b0_onda6+nsample,wavetmp);
    }
    else if( Ichan==7 ){
      Ttmp    = b0_tave7;
      ped1tmp = b0_ped7;
      ped2tmp = b0_pedh7;
      Itmp    = b0_imax7;
      std::copy(b0_onda7,b0_onda7+nsample,wavetmp);
    }
  }
  else if( Iboard==1 ){
    if( Ichan==0 ){
      Ttmp    = b1_tave0;
      ped1tmp = b1_ped0;
      ped2tmp = b1_pedh0;
      Itmp    = b1_imax0;
      std::copy(b1_onda0,b1_onda0+nsample,wavetmp);
    }
    else if( Ichan==1 ){
      Ttmp    = b1_tave1;
      ped1tmp = b1_ped1;
      ped2tmp = b1_pedh1;
      Itmp    = b1_imax1;
      std::copy(b1_onda1,b1_onda1+nsample,wavetmp);
    }
    else if( Ichan==2 ){
      Ttmp    = b1_tave2;
      ped1tmp = b1_ped2;
      ped2tmp = b1_pedh2;
      Itmp    = b1_imax2;
      std::copy(b1_onda2,b1_onda2+nsample,wavetmp);
    }
    else if( Ichan==3 ){
      Ttmp    = b1_tave3;
      ped1tmp = b1_ped3;
      ped2tmp = b1_pedh3;
      Itmp    = b1_imax3;
      std::copy(b1_onda3,b1_onda3+nsample,wavetmp);
    }
    else if( Ichan==4 ){
      Ttmp    = b1_tave4;
      ped1tmp = b1_ped4;
      ped2tmp = b1_pedh4;
      Itmp    = b1_imax4;
      std::copy(b1_onda4,b1_onda4+nsample,wavetmp);
    }
    else if( Ichan==5 ){
      Ttmp    = b1_tave5;
      ped1tmp = b1_ped5;
      ped2tmp = b1_pedh5;
      Itmp    = b1_imax5;
      std::copy(b1_onda5,b1_onda5+nsample,wavetmp);
    }
    else if( Ichan==6 ){
      Ttmp    = b1_tave6;
      ped1tmp = b1_ped6;
      ped2tmp = b1_pedh6;
      Itmp    = b1_imax6;
      std::copy(b1_onda6,b1_onda6+nsample,wavetmp);
    }
    else if( Ichan==7 ){
      Ttmp    = b1_tave7;
      ped1tmp = b1_ped7;
      ped2tmp = b1_pedh7;
      Itmp    = b1_imax7;
      std::copy(b1_onda7,b1_onda7+nsample,wavetmp);
    }
  }
  else{
    cout << "ERROR: Board/Channel not found. Resetting arrays! " 
      << Iboard << " " << Ichan << " " << endl;
    Ttmp    = -999.;
    ped1tmp = -999.;
    ped2tmp = -999.;
    Itmp    = -999 ;
    std::fill_n(wavetmp,nsample,-999.);
  }

  if( Ttmp!=-999. ){          // Channel exists
    //
    // Baseline evaluation
    //
  
    if(evflag==0) {
      TAmax = time_out[Itmp];   // Time (in ns) of the max amplitude   
      TWmin = TAmax-110.;       // Time window for baseline evaluation:
      TWmax = TAmax- 60.;       // 50 ns before the start of the wave
    }else if( evflag == 1){
      Itmp = 63;
      TAmax = time_out[Itmp];   // Time (in ns) of the max amplitude   
      TWmin = TAmax-110.;       // Time window for baseline evaluation:
      TWmax = TAmax- 60.;       // 50 ns before the start of the wave
    }
    
    nVal = 0;
    sum  = 0.;
    sumsq = 0.;
    for(int Ibin=0; Ibin<Itmp; Ibin++ ){
      Time = time_out[Ibin];
      if ( Time>TWmin && Time<TWmax ){
	      nVal++;
        sumsq += wavetmp[Ibin]*wavetmp[Ibin];
      	sum = sum + wavetmp[Ibin];
      }
    }
    if( nVal!=0 ){
      Btmp = sum / nVal;
      Brmstmp = TMath::Sqrt(sumsq/nVal - Btmp*Btmp);
    }

    ped1tmp = 0;
    // Subtract baseline to the waveform and evaluate pedestal
    for(int Ibin=0; Ibin<nsample; Ibin++ ){
      wavetmp[Ibin] = wavetmp[Ibin] - Btmp;
      if(Ibin < 50) ped1tmp += wavetmp[Ibin]; //fino a 200 ns
    }
    ped1tmp *= dt/50;
    //
    // Charge with baseline subtracted event-by-event
    //
    sum  = 0.;
    QbinMin = Itmp-Qpeak_min*camp;
    QbinMax = Itmp+Qpeak_max*camp;
   
    if( QbinMin<0 ) QbinMin = 0;
    if( QbinMax>maxNsample ) QbinMax = maxNsample;
    for(int Ibin=QbinMin; Ibin<QbinMax; Ibin++ ){
      sum = sum + wavetmp[Ibin]-Btmp;
    }
    Qtmp = sum*dt/50.;
    Vtmp=wavetmp[Itmp];
  }
  else{
    Btmp    = -999.;
    Qtmp    = -999.;
    Vtmp    = -999.;
  }
 
}

//=============================================================================
// Get time from template fit to waveform
//=============================================================================

void analysis_CRT::getTemplateFit(TGraphErrors *gt)
{

  TGraphErrors *gwf;
  TF1 *fitf;

  double ex[maxNsample] = {0.05}; //X-error not useful .. just for display graphs
  double dig_res = 0.5;
  double ey[maxNsample]; // = {TMath::Sqrt(Brmstmp*Brmstmp + dig_res*dig_res)};

  //cout << TMath::Sqrt(Brmstmp*Brmstmp + dig_res*dig_res) << endl;
  std::fill_n(ey, maxNsample, 1.25); //TMath::Sqrt(Brmstmp*Brmstmp + dig_res*dig_res));
  int t_peak = 0;
  double v_max = 0.;
  double CF = 0.2;

  fitTmea = 0.;
  fitChi2 = 0.;
  std::fill_n(fitPar, 3, 0.);
  std::fill_n(fitErr, 3, 0.);

  float Vpk = wavetmp[Itmp];

  if (Qtmp > 10)
  {

    TGraphErrors gwf(120, time_out, wavetmp, ex, ey);

    v_max = wavetmp[Itmp];
    t_peak = time_out[Itmp];
    sideTmp = SideTmp;
    scintTmp = ScintTmp;

    TF1 fitf("fitf", spline_fit, t_peak-80, t_peak -20, 3); //-40 -15 ok

    // Three parameter fit ( 0: norm, 1:t0, 2:baseline )
    //Set start parameters and limits

      fitf.SetParameter(0,v_max);
      fitf.SetParLimits(0, v_max*0.9, v_max*1.1);
      fitf.SetParameter(1, t_peak - 50); //per quello di daniele fa -120
      fitf.SetParLimits(1, t_peak-80, t_peak-30); //per quello di daniele fa -150 -80
      fitf.SetParameter(2,  0.);
      fitf.SetParLimits(2,-10,10);


    int res = gwf.Fit("fitf", "REMQ");
    fit_log << res << endl;

      fitTmea = fitf.GetParameter(1); //fitf.GetX(v_max * CF);
      fitPar[0] = fitf.GetParameter(0);
      fitPar[1] = fitf.GetParameter(1);
      fitPar[2] = fitf.GetParameter(2);
      fitErr[0] = fitf.GetParError(0);
      fitErr[1] = fitf.GetParError(1);
      fitErr[2] = fitf.GetParError(2);


      if (fitf.GetNDF() != 0)
        fitChi2 = fitf.GetChisquare() / fitf.GetNDF();

    if(plot){
      cout << "Tpeak: " << t_peak << endl;

      cout << "Ev n. " << evtmp << " - Side n. " << SideTmp << " - Scint n. " << ScintTmp << endl;
      cout << "Fit Time parameter: " << fitPar[1] << endl;
      cout << "Fit Chi2: " << fitChi2 << endl; 

      gwf.Draw("Al");
      gwf.SetTitle(Form("Template Fit 2020 - ev n.%i - side %i - scint %i", evtmp, SideTmp, ScintTmp));
      gwf.SetMarkerStyle(22);
      gwf.SetMarkerSize(0.5);
      gwf.Draw("pe");
      //fitf.SetRange(t_peak-40, t_peak +15);
      fitf.Draw("same");
      gwf.GetXaxis()->SetRangeUser(time[Itmp]-100, time[Itmp]+200);
      gStyle->SetOptFit(1);
      while (go_ahead_flag == 0)
      {
        gPad->Modified();
        gPad->Update();
        gSystem->ProcessEvents();
        usleep(10000);
      }
      go_ahead_flag = 0;
    }

  }
}

void analysis_CRT::getTemplate2Fit()
{

  TGraphErrors *gwf;
  TF1 *fitf;

  double ex[maxNsample] = {0.05}; //X-error not useful .. just for display graphs
  double dig_res = 0.5;
  double ey[maxNsample]; // = {TMath::Sqrt(Brmstmp*Brmstmp + dig_res*dig_res)};

  //cout << TMath::Sqrt(Brmstmp*Brmstmp + dig_res*dig_res) << endl;
  double sigmav = TMath::Sqrt(Brmstmp*Brmstmp + dig_res*dig_res);
  std::fill_n(ey, maxNsample, sigmav);
  int t_peak = 0;
  double v_max = 0.;
  double CF = 0.4;

  fitTmea = 0.;
  fitChi2 = 0.;
  std::fill_n(fitPar, 3, 0.);
  std::fill_n(fitErr, 3, 0.);

  float Vpk = wavetmp[Itmp];

  if (Vpk > 100 && time[Itmp] > 200 && time[Itmp] < 450 && Vpk < 1800)
  {

    TGraphErrors gwf(120, time_out, wavetmp, ex, ey);

    t_peak = time_out[Itmp];
    sideTmp = SideTmp;
    scintTmp = ScintTmp;

    TF1 fitf("fitf", spline5_fit, t_peak - 50, t_peak -15, 3); //-40 -15 ok

    fitf.SetParameter(0, Vpk);
    fitf.SetParLimits(0, Vpk*0.9, Vpk*1.1);
    fitf.SetParameter(1, t_peak-230); //partendo dal templ di Daniele è -120, sennò partendo da RbnFit è -230
    fitf.SetParLimits(1, t_peak-300, t_peak-150);
    fitf.SetParameter(2, 0.);
    fitf.SetParLimits(2, -5*sigmav, 5*sigmav);

    int res = gwf.Fit("fitf", "REMQ");

      fitTmea = fitf.GetParameter(1) ; //fitf.GetX(Vpk * CF);
      fitPar[0] = fitf.GetParameter(0);
      fitPar[1] = fitf.GetParameter(1);
      fitPar[2] = fitf.GetParameter(2);
      fitErr[0] = fitf.GetParError(0);
      fitErr[1] = fitf.GetParError(1);
      fitErr[2] = fitf.GetParError(2);


      if (fitf.GetNDF() != 0)
        fitChi2 = fitf.GetChisquare() / fitf.GetNDF();

    if(plot){
      cout << "Tpeak: " << t_peak << endl;

      cout << "Ev n. " << evtmp << " - Side n. " << SideTmp << " - Scint n. " << ScintTmp << endl;
      cout << "Fit Time parameter: " << fitPar[1] << endl;
      cout << "Fit Chi2: " << fitChi2 << endl; 

      gwf.Draw("Al");
      gwf.SetTitle(Form("NEW Template Fit - ev n.%i - side %i - scint %i", evtmp, SideTmp, ScintTmp));
      gwf.SetMarkerStyle(22);
      gwf.SetMarkerSize(0.5);
      gwf.Draw("pe");
      //fitf.SetRange(t_peak-40, t_peak +15);
      fitf.Draw("same");
      gwf.GetXaxis()->SetRangeUser(time[Itmp]-100, time[Itmp]+200);
      gStyle->SetOptFit(1);
      while (go_ahead_flag == 0)
      {
        gPad->Modified();
        gPad->Update();
        gSystem->ProcessEvents();
        usleep(10000);
      }
      go_ahead_flag = 0;
    }

  }
}


//=============================================================================
// Get time from logn fit to waveform
//=============================================================================
void analysis_CRT::getTlogn() // a and b are nothing, just to avoid changing the .h
{


  int Imin = 50;
  float Vmax = wavetmp[Itmp];
  float Vmin = 0.;

  float minVal = 0.03*Vmax;
  float Tpeak = time_out[Itmp];
  float t_min = 200.;
  float t_max = 450.;

  // Itmp: Max of waveform
  for (int iLoop = Itmp; iLoop > Itmp - 20; iLoop--)
  {
    if (wavetmp[iLoop] < minVal && Imin == 50)
    {
      Imin = iLoop;
      //cout<<"min "<<iLoop<<" "<<Itmp<<endl;
      Vmin = wavetmp[iLoop];
    }
  }

  if (Vmax > 50 && time[Itmp] > t_min && time[Itmp] < t_max && Vmax < 1400)
  {
    Double_t errors[1024] = {0.};
    Double_t empty[1024] = {0.};

    TGraphErrors g_crystal(100, time, wavetmp, empty, errors);
    // TF1: Last variable is number of fit parameters
    TF1 f_logn("f_logn", logn, time[Imin], time[Itmp + 10], 4);

    f_logn.SetParNames("#eta", "#sigma", "t_{0}", "Norm");
    //f_logn->SetParLimits(0,-5.,-0.001); //prove con giani
    f_logn.SetParLimits(0, -2., -0.001);
    f_logn.SetParLimits(1, 10., 300);
    //f_logn->SetParLimits(2,200,t_peak+50.);
    //if( Iboard==0 ) f_logn->SetParLimits(2,300,350);
    //if( Iboard==1 ) f_logn->SetParLimits(2,250,300);
    //f_logn->SetParLimits(3,1.,1000000);
    f_logn.SetParLimits(3, 100. * Vmax, 200. * Vmax);
    f_logn.SetParameters(-1., 50., Tpeak - 10., Vmax * 140);

    //f_pol->SetParLimits(2, 1800, 1900);
    //g_crystal->Fit("f_pol", "RBOQ","",time[i_min], time[i_max-1]);
    //f_logn.FixParameter(2, Tpeak);

    int result;
    /*
    result = g_crystal.Fit("f_logn", "RBOQ", "", time[Imin], time[Itmp - 1]);
    result = g_crystal.Fit("f_logn", "RBOQ");
    f_logn.SetParLimits(2, 200, Tpeak + 50.); //prove
    double eta_min = f_logn.GetParameter(0) * 1.25;
    double eta_max = f_logn.GetParameter(0) * 0.75;
    f_logn.SetParLimits(0, eta_min, eta_max); //prove
    result = g_crystal.Fit("f_logn", "RBOQ");
    */
    g_crystal.Fit("f_logn", "RBOQ");
    double const_frac = 0.3; //1;
    //tmean[Iboard][Ichan] = f_pol->GetX((v_max + v_min)/2.*const_frac);
    lognTimetmp = f_logn.GetX(Vmax * const_frac); //-timemu[Iboard][Ichan];
    lognChi2tmp = f_logn.GetChisquare() / f_logn.GetNDF();
    lognfitPar[0] = f_logn.GetParameter(0);
    lognfitPar[1] = f_logn.GetParameter(1);
    lognfitPar[2] = f_logn.GetParameter(2);
    lognfitPar[3] = f_logn.GetParameter(3);

    lognfitErr[0] = f_logn.GetParError(0);
    lognfitErr[1] = f_logn.GetParError(1);
    lognfitErr[2] = f_logn.GetParError(2);
    lognfitErr[3] = f_logn.GetParError(3);
    //tedge= f_pol->GetX(v_max*const_frac)-4*int(f_pol->GetX(v_max*const_frac)/4);

    if(plot){
      cout << "Ev n. " << evtmp << " - Side n. " << SideTmp << " - Scint n. " << ScintTmp << endl;
      g_crystal.Draw("Al");
      g_crystal.SetTitle(Form("Logn Fit - ev n.%i - side %i - scint %i", evtmp, SideTmp, ScintTmp));
      g_crystal.SetMarkerStyle(22);
      g_crystal.SetMarkerSize(0.5);
      g_crystal.Draw("pe");
      f_logn.SetRange(200, 600);
      f_logn.Draw("same");
      g_crystal.GetXaxis()->SetRangeUser(time[Itmp]-100, time[Itmp]+200);
      gStyle->SetOptFit(1);
      while (go_ahead_flag == 0)
      {
        gPad->Modified();
        gPad->Update();
        gSystem->ProcessEvents();
        usleep(10000);
      }
      go_ahead_flag = 0;
    }

  }
}


int main(int argc, char *argv[]){
TApplication *myapp = new TApplication("myapp", 0, 0);

analysis_CRT *a = new analysis_CRT( Form("run%s_CRTNew.root", argv[1]));

a->Loop(Form("run%s_ana.root", argv[1]), 1);
myapp->Run(true);
}
