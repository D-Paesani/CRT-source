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



double pi = TMath::Pi();

int nEvents = 2e7;
int doVeto = 1;
int doSave = 1;
int doGame = 1;
double enCut = 1;
#define nShow 50

float hGen = 3.5;
float lGen = 100;
float lBar = 2.5;
float hBar = 1.5;
float lVeto = 2.5;
float hVeto = 1.8;

double xIn, xOut, yIn, yOut, dx, xMuon, aMuon, yL, yR, xM, en, xVeto, ang, an;
int isHit, isVeto, hitCount;

TF1 *muonAngle = new TF1("th_pdf", "2/ TMath::Pi() * cos(x) * cos(x)", -pi/2, pi/2);

TH1F *xHist = new TH1F("xHist", "xHist", 1000, -lBar, lBar);
TH1F *dxHist = new TH1F("dxHist", "dxHist", 1000, 0, 1.1*sqrt(lBar*lBar+hBar*hBar));
TH2F *hitMap = new TH2F("map", "map", 500, -lBar, lBar, 100, -hBar, 2*hBar);
TH1F *enHist = new TH1F("energy", "energy", 1000, 0, 10);

auto barArea = new TBox(-lBar/2, 0, lBar/2, hBar);
auto h1 = new TH2F("hit", "hit", 1, -3, 3, 1, -1, 5);
auto h2 = new TH2F("hits", "hits", 1, -2, 2, 3, -1, 2);
auto genLine = new TLine(-3, hGen, 3, hGen);
auto vetoLine = new TLine(-lVeto/2, hVeto, lVeto/2, hVeto);

TCanvas *c1 = new TCanvas("hits", "hits",200,10,900,900);
auto c1t = new TPad("hit","hits", 0.05, 0.5, 0.95, 0.95);
auto c1b = new TPad("hits","hits",  0.05, 0.05, 0.95, 0.45);



void dxGet() {

    isHit = 0;
    isVeto = 0;
    dx = -1;

    ang = 1/TMath::Tan(aMuon);
    an = 1/TMath::Tan(abs(aMuon));
    xM = xMuon;
    if (aMuon < 0) {xM = -xM;} 

    yR = hGen - an * (xM - lBar/2);
    yL = hGen - an * (xM + lBar/2);
    xVeto = xIn = xM - (hGen - hVeto)/an;

    if (yR < 0 || yL > hBar) {return;}
    if ((abs(xVeto) < lVeto/2) && (doVeto)) {isVeto=1;}

    if (yR > hBar) {
        isHit = 20;
        xIn = xM - (hGen - hBar)/an;
        yIn = hBar;
    } else {    
        isHit = 10;
        xIn = lBar/2;
        yIn=yR;
    }

    if (yL < 0) {
        isHit = isHit + 2;
        xOut = xM - (hGen)/an;
        yOut = 0;
    } else {    
        isHit = isHit + 1;
        xOut = -lBar/2;
        yOut = yL;
    }

    if (!isHit) {return;}

    dx = sqrt((xIn-xOut)*(xIn-xOut)+(yIn-yOut)*(yIn-yOut));

    if (aMuon < 0) { 
        xIn = - xIn;
        xOut = - xOut;
    } 
}



void Game() {

    if (!doGame || hitCount>50) {return;}

    c1t->cd();
    h1->Draw();
    barArea->Draw("same");
    genLine->Draw("same");
    if (doVeto) {vetoLine->Draw();}

    cout<<"hitNo: "<<hitCount<<"\n";     
    cout << "muonOrigin: " << xMuon << " , " << hGen << endl;
    cout << "muonAngle: " << aMuon << endl;
    cout << "inPoint: " << xIn << " , " << yIn << endl;
    cout << "outPoint: " << xOut << " , " << yOut << endl;

    TF1 l("l", "[0] - (x - [1]) * [2]", -100, 100);
    l.SetParameter(0, hGen);
    l.SetParameter(1, xMuon);
    l.SetParameter(2, -ang);
    l.Draw("same");

    TMarker genPoint(xMuon,hGen,20);
    genPoint.SetMarkerColor(kRed);
    genPoint.Draw();

    if (!isVeto) {

        TLine *ll = new TLine(xIn, yIn, xOut, yOut);
        ll->SetLineColor(kGreen);
        ll->SetLineWidth(5);
        ll->SetLineStyle(1);
        ll->Draw();

        TMarker inPoint(xIn,yIn,20);
        inPoint.SetMarkerColor(kGreen);
        inPoint.Draw();
        TMarker outPoint(xOut,yOut,20);
        outPoint.SetMarkerColor(kGreen);
        outPoint.Draw();

        c1b->cd();
        TLine *lll = new TLine(xIn, yIn, xOut, yOut);
        lll->SetLineWidth(1);
        lll->SetLineColor(kBlue);
        lll->Draw("same");
        c1t->cd();
    
    }
    
    gPad->Modified();
    gPad->Update();
    gSystem->ProcessEvents();

    char ttt;
    //cin >> ttt;

    std::this_thread::sleep_for(std::chrono::milliseconds(50 + 1500*(hitCount==1))); 
}



void BarraCosmicQ() {

    hitCount = 0;
    
    c1b->Draw();
    c1t->Draw();
    c1b->cd();
    h2->Draw();
    //barArea->Draw("same");

    for (int jj = 0; jj<nEvents; jj++){

        xMuon = gRandom->Uniform(-lGen,lGen); 
        aMuon = muonAngle->GetRandom();
        dxGet();

        if (isHit > 0) { 

            hitCount++;
            
            en = gRandom->Landau(2*dx, 0.3/2*0.5);

            xHist->Fill(xOut);
            dxHist->Fill(dx);
            hitMap->Fill(xIn,yIn);
            hitMap->Fill(xOut,yOut);  
            enHist->Fill(en);  

            Game();
        }
    } 

    auto ccc = new TCanvas("HitDist","HitDist");
    ccc->cd();
    xHist->Draw();
    xHist->GetXaxis()->SetTitle("x [m]");
    xHist->GetYaxis()->SetTitle("Hits");

    auto cc = new TCanvas("dxHist","dxDist");
    cc->cd();
    dxHist->Draw();
    dxHist->GetXaxis()->SetTitle("dx [m]");
    dxHist->GetYaxis()->SetTitle("Entries");

    auto cccc = new TCanvas("map","map");
    cccc->cd();
    hitMap->Draw("colz");

    auto ccccc = new TCanvas("e","e");
    ccccc->cd();
    enHist->Draw("colz");

    if (doSave) {
    TFile *f = new TFile("BarraCosmicQ.root", "RECREATE");
    dxHist->Write();
    xHist->Write();
    enHist->Write();
    c1->Write();
    f->Close();
    }

    cout<<"generated "<<nEvents<<" cosmics, "<<hitCount<<" on target"<<"\n";

}
























