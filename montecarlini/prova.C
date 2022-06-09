#include <iostream>
#include <fstream>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <chrono>
#include "TSystem.h"
#include <TH2.h>
#include <TStyle.h>
#include "TROOT.h"
#include "TCanvas.h"
#include "TRootCanvas.h"
#include <TProfile.h>
#include <TPaveText.h>
#include <TPaveStats.h>
#include <TF1.h>
#include <TTree.h>
#include <TFile.h>
#include "TMarker.h"
#include <string>
#include <TGraphErrors.h>
#include <TMultiGraph.h>
#include <vector>
#include <TString.h>
#include "TTimer.h"
#include "TRandom.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TLatex.h"
#include "TGClient.h"
#include <TDatime.h>
#include <time.h>
#include <TLine.h>

using namespace std;

int nEvents = 5e6;

float hGen = 1000;
float lGen = 1000;
float lBar0 = 160;
float lBar1 = lBar0;
float hBar1 = 42;
float xBar1 = 0;

TF1 *muonAngle = new TF1("th_pdf", "2/ TMath::Pi() * cos(x) * cos(x)", -TMath::Pi()/2, TMath::Pi()/2);

int nBins = nEvents/5e3;
TH1F anggen("anggen", "anggen", nBins, -2, 2);
TH1F coseni("coseni", "coseni", nBins, 0, 1);
TH1F angulus("ang", "ang", nBins, -2, 2);
TH1F angulusFromgen("angulusFromgen", "angulusFromgen", nBins, -2, 2);
TH1F zeta1("zeta1", "zeta1", nBins, -1.5*lBar1/2, 1.5*lBar1/2);
TH1F zeta0("zeta0", "zeta0", nBins, -1.5*lBar0/2, 1.5*lBar0/2);
TH1F angulussm("angulussm", "angulussm", nBins, -2, 2);
TH1F zeta1sm("zeta1sm", "zeta1sm", nBins, -1.5*lBar1/2, 1.5*lBar1/2);
TH1F zeta0sm("zeta0sm", "zeta0sm", nBins, -1.5*lBar0/2, 1.5*lBar0/2);

double smearZ(double z) {
    
    double resol = 2;

    return gRandom->Gaus(z, resol); 
}

void prova() {

    long hitCount{0};

    gRandom->SetSeed();

    for (int jj = 0; jj<nEvents; jj++){

        int isHit = 0;

        double xMuon = gRandom->Uniform(-lGen,lGen); 
        double aMuon = muonAngle->GetRandom();
        double tanMuon = TMath::Tan(aMuon);
        double cosMuon = TMath::Cos(aMuon);

        anggen.Fill(aMuon);
        coseni.Fill(cosMuon);

        double x1 = xMuon - (hGen - hBar1)*tanMuon;
        double x0 = xMuon - (hGen)*tanMuon;

        if (2*TMath::Abs(x0) > lBar0) continue;
        if (x1 - xBar1 > lBar1/2) continue;
        if (x1 - xBar1 < -lBar1/2) continue;

        hitCount++;

        zeta0.Fill(x0);
        zeta1.Fill(x1);
        angulusFromgen.Fill(aMuon);
        angulus.Fill(TMath::ATan((x1 - x0)/hBar1));

        x1 = smearZ(x1-xBar1) + xBar1;
        x0 = smearZ(x0);

        zeta0sm.Fill(x0);
        zeta1sm.Fill(x1);

        double asmear = TMath::ATan((x1 - x0)/hBar1);

        angulussm.Fill(asmear);
    }

    TFile outfile("prova.root", "recreate");
    outfile.cd();

    coseni.Write();
    anggen.Write();
    zeta0sm.Write();
    zeta1sm.Write();
    angulussm.Write();
    zeta0.Write();
    zeta1.Write();
    angulus.Write();
    angulusFromgen.Write();

    cout<<"generated "<<nEvents<<" cosmics, "<<hitCount<<" on target"<<"\n";

}
























