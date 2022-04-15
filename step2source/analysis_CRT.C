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

using namespace std;

//=============================================================================
// Fit function for template
//=============================================================================

int sideTmp, scintTmp, modTmp;

//=============================================================================
// User variables
//=============================================================================

int const Nchan = 16;
const double CF = 0.15;

int Itmp;
// int         ItmpLaser;
float Qtmp, Ttmp, Btmp, Brmstmp, Vtmp;
float ped1tmp, ped2tmp, lognTimetmp = 0, lognChi2tmp = 0;
Double_t wavetmp[maxNsample];
int BoaTmp, ChaTmp, DaqTmp, ScintTmp, SideTmp, ModTmp, PseudoTtmp;

int scint[NCHAN], side[NCHAN], mod[NCHAN];

int evtmp;

float Time, Chi2, pseudoT;

//=============================================================================
// Loop function
//=============================================================================

// void analysis_CRT::Loop()
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
  // by  b_branchname->GetEntry(ientry); //read only this branch
  if (fChain == 0)
    return;

  // **************************************************************************
  // Output file
  // **************************************************************************

  cout << "Selected output file: " << OutputFile << endl;
  TFile *outFile = new TFile(OutputFile, "recreate");
  TTree *CRT = new TTree("CRT", "CRT");
  BookOutput(CRT); // Book histos and ntuple

  // **************************************************************************
  // Loop on entries
  // **************************************************************************
  Long64_t nentries = fChain->GetEntriesFast();

  cout << "nentries: " << nentries << endl;

  Long64_t nbytes = 0, nb = 0;

  FILE *mapFile;
  mapFile = fopen("CRT_map.dat", "r");
  for (int Iloop = 0; Iloop < Nchan; Iloop++)
  {
    fscanf(mapFile, "%d %d %d %d", &ChaTmp, &ScintTmp, &SideTmp, &ModTmp);
    scint[ChaTmp] = ScintTmp;
    side[ChaTmp] = SideTmp;
    mod[ChaTmp] = ModTmp;
  }
  fclose(mapFile);

  for (Long64_t jentry = 0; jentry < nentries; jentry++)
  {

    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0)
      break;
    nb = fChain->GetEntry(jentry);
    nbytes += nb;

    if (jentry % 1000 == 0)
      cout << "Number of processed events: " << jentry << endl;

    ntrig_out = ntrig;
    evnum_out = evnum;
    nsample = num;
    evtmp = evnum;
    std::copy(time, time + nsample, time_out);

    int cryTot = 0;

    for (int Ichan = 0; Ichan < Nchan; Ichan++)
    {

      GetValues(Ichan);
      // Get values for a given board/channel

      if (Qtmp > 50)
      {

        iDAQ[cryTot] = DaqTmp;
        iScint[cryTot] = ScintTmp;
        iSide[cryTot] = SideTmp;
        Qval[cryTot] = Qtmp;
        Tval[cryTot] = Ttmp;
        pedL[cryTot] = ped1tmp;
        pedH[cryTot] = ped2tmp;
        iMax[cryTot] = Itmp;
        Vmax[cryTot] = Vtmp;
        std::copy(wavetmp, wavetmp + nsample, wave[cryTot]);
        bline[cryTot] = Btmp;
        pseudot[cryTot] = PseudoTtmp;
        cryTot++;
      }

  }   // Loop on Ichan

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

  CRT->Branch("ntrig", &ntrig_out, "ntrig/I");
  CRT->Branch("evnum", &evnum_out, "evnum/I");
  CRT->Branch("nsample", &nsample, "nsample/I");     // function of nCry?
  CRT->Branch("time", &time_out, "time[nsample]/D"); // function of nCry?
  CRT->Branch("nCry", &nCry, "nCry/I");
  CRT->Branch("iDAQ", &iDAQ, "iDAQ[nCry]/I");
  CRT->Branch("iScint", &iScint, "iScint[nCry]/I");
  CRT->Branch("iSide", &iSide, "iSide[nCry]/I");
  CRT->Branch("iMod", &iMod, "iMod[nCry]/I");
  CRT->Branch("iMax", &iMax, "iMax[nCry]/I");
  CRT->Branch("Vmax", &Vmax, "Vmax[nCry]/D");
  CRT->Branch("Qval", &Qval, "Qval[nCry]/D");
  CRT->Branch("Tval", &Tval, "Tval[nCry]/D");
  CRT->Branch("pedL", &pedL, "pedL[nCry]/D");
  CRT->Branch("pedH", &pedH, "pedH[nCry]/D");
  CRT->Branch("wave", &wave, "wave[nCry][1024]/D");
  CRT->Branch("bline", &bline, "bline[nCry]/D");
  CRT->Branch("templTime", &pseudot, "templTime[nCry]/D"); // per compatibilità step 3
}

//=============================================================================
// Get values for a given board/channel
//=============================================================================
void analysis_CRT::GetValues(int Ichan)
{

  int nVal;
  float TAmax, TWmin, TWmax;
  float Time, sum, sumsq;
  int QbinMin, QbinMax;
  float const camp = 0.25;
  float dt = 1 / camp;
  float const Qpeak_min = 60.;  //  -60 ns from Tave
  float const Qpeak_max = 190.; // +190 ns from Tave

  // Read map [do it once, not for all channels!!!]

  DaqTmp = Ichan;
  ScintTmp = scint[Ichan];
  SideTmp = side[Ichan];
  ModTmp = mod[Ichan];

  Ttmp = tave[Ichan];
  ped1tmp = ped[Ichan];
  ped2tmp = pedh[Ichan];
  Itmp = imax[Ichan];
  std::copy(onda[Ichan], onda[Ichan] + nsample, wavetmp);

  if (Ttmp != -999.)
  { // Channel exists
    //
    // Baseline evaluation
    //

    TAmax = time[Itmp];  // Time (in ns) of the max amplitude
    TWmin = TAmax - 110.; // Time window for baseline evaluation:
    TWmax = TAmax - 60.;  // 50 ns before the start of the wave

    nVal = 0;
    sum = 0.;
    sumsq = 0.;

    for (int Ibin = 0; Ibin < Itmp; Ibin++)
    {
      Time = time[Ibin];
      if (Time > TWmin && Time < TWmax)
      {
        nVal++;
        sumsq += wavetmp[Ibin] * wavetmp[Ibin];
        sum = sum + wavetmp[Ibin];
      }
    }
    if (nVal != 0)
    {
      Btmp = sum / nVal;
      Brmstmp = TMath::Sqrt(sumsq / nVal - Btmp * Btmp);
    }

    ped1tmp = 0;
    // Subtract baseline to the waveform and evaluate pedestal
    for (int Ibin = 0; Ibin < nsample; Ibin++)
    {
      wavetmp[Ibin] = wavetmp[Ibin] - Btmp;
      if (Ibin < 50)
        ped1tmp += wavetmp[Ibin]; // fino a 200 ns
    }
    ped1tmp *= dt / 50;
    //
    // Charge with baseline subtracted event-by-event
    //
    sum = 0.;
    QbinMin = Itmp - Qpeak_min * camp;
    QbinMax = Itmp + Qpeak_max * camp;

    if (QbinMin < 0)
      QbinMin = 0;
    if (QbinMax > maxNsample)
      QbinMax = maxNsample;
    for (int Ibin = QbinMin; Ibin < QbinMax; Ibin++)
    {
      sum = sum + wavetmp[Ibin] - Btmp;
    }
    Qtmp = sum * dt / 50.;
    Vtmp = wavetmp[Itmp];

    if (Qtmp>50){
      TGraph wgr(400, time, wavetmp); //1) faccio un tgraph con l'onda //PARTIRE DA QUA

      TSpline5 wsp = TSpline5("wsp", &wgr);

      auto spf = [&wsp](double *x, double *){ return wsp.Eval(x[0]); };
      TF1 fitf = TF1("fitf", spf, TAmax - 50, TAmax + 50, 0);
      double norm = fitf.GetMaximum(TAmax - 50, TAmax + 50);
      double th = norm*CF;
      PseudoTtmp = fitf.GetX(th);
    }
    else PseudoTtmp = -999;
  }
  else
  {
    Btmp = -999.;
    Qtmp = -999.;
    Vtmp = -999.;
  }
}

int main(int argc, char *argv[])
{
  TApplication *myapp = new TApplication("myapp", 0, 0);

  analysis_CRT *a = new analysis_CRT(Form("run%s_CRTNew.root", argv[1]));

  a->Loop(Form("run%s_ana.root", argv[1]), 1);
  myapp->Run(true);
}
