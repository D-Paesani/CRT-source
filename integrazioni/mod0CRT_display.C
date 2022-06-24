// ANCORA NON FINITOOO
#include "module0_map.h"

#include <stdio.h>
#include <stdlib.h>
#include <TH1F.h>
#include <TBox.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TH2F.h>
#include <TButton.h>
#include <TList.h>
#include <TFile.h>
#include <TSystem.h>
#include <TF1.h>
#include <TRandom.h>
#include <TMath.h>
#include <TGraph.h>
#include <TLine.h>
#include <TROOT.h>
#include <TApplication.h>

using namespace std;

int plot = 1;
double MeVperCm = 6.167; // 21 MeV / 3.4 cm
double landauScalePerCm = MeVperCm * 0.08; // a caso, va fatto meglio

int go_ahead_flag = 0;

void go_ahead(){
  go_ahead_flag = 1;
}

class Scintillator{
    public:
        double yDown, yUp, slope, dx, en, xLeft, xRight;
        TString name;
        int noHitInfo, isHit;

        Scintillator(TString name_, double xCenter, double xLength, double yCenter, double yLength, double slope_):
            name(name_),
            xLeft(xCenter - xLength/2), xRight(xCenter + xLength/2), yDown(yCenter - yLength/2), yUp(yCenter + yLength/2), 
            slope(slope_)
        {};
        void ProcessEvent(int checkOnly, double xGen, double yGen, double m, TGraph &g){

            double yLeft = yGen - m * (xLeft - xGen);
            double yRight = yGen - m * (xRight - xGen);

            noHitInfo = checkOnly;
            double yLRmin = min(yLeft, yRight);
            double yLRmax = max(yLeft, yRight);

            //yLRmin: -inf - 0 - yDown - 1 - yUp - 2 - +inf
            //yLRmax: -inf - 0 - yDown - 1 - yUp - 2 - +inf
            //X: no Cross - !: imposs
            //L, R, T, B: left, right, top and bottom sides
            //          yLRmin
            //         0      1      2
            //       0 X      !      !
            //yLRmax 1 LB/RB  LR     !
            //       2 TB    TL/TR   X

            double xCross1, xCross2, yCross1, yCross2;

            if(yLRmin > yUp || yLRmax < yDown){
                isHit = 0;
            } //yLRmin not 2 - yLRmax not 0  

            else{
                isHit = 1;
                if (yLRmax < yUp) { //yLRmax a 1 (yDown..yUp)
                    // LB/RB - LR
                    if (yLeft > yRight){
                        xCross1 = xLeft;
                        yCross1 = yLeft;
                    }
                    else {
                        xCross1 = xRight;
                        yCross1 = yRight;
                    }
                }
                else { //yLRmax a 2 (yUp..+inf)
                    // TB - TR/TL
                    double xUp = (yGen - yUp)/m + xGen;
                    xCross1 = xUp;
                    yCross1 = yUp;
                }

                if (yLRmin < yDown){ //yLRmin a 0 (-inf..yDown)
                    // LB/RB - TB
                    double xDown = (yGen - yDown)/m + xGen;
                    xCross2 = xDown;
                    yCross2 = yDown;
                }
                else { //yLRmin a 1 (yDown..yUp)
                    // LR - TL/TR
                    if (yLeft < yRight){
                        xCross2 = xLeft;
                        yCross2 = yLeft;
                    }
                    else {
                        xCross2 = xRight;
                        yCross2 = yRight;
                    }
                }

                dx = sqrt( (xCross1-xCross2)*(xCross1-xCross2) + (yCross1-yCross2)*(yCross1-yCross2) );
                en = gRandom->Landau(MeVperCm*dx, landauScalePerCm*dx);
                en *= 1 + slope * (xCross1 + xCross2)/2; // ma va integrato o basta così?
            
                if(plot){
                    g.AddPoint(xCross1, yCross1);
                    g.AddPoint(xCross2, yCross2);
                }
            
            }

        }

        TBox * GetBox(){
            TBox *b = new TBox(xLeft, yDown, xRight, yUp);
            return b;
        }
};


void Qsimul_CRTandScintillatorvsMod0(){

    int nev = 1e5;

    double yGen=90, yOffsetCRT = 60, yOffsetPadel = 0, semiTotLengthCrt = 10, semiLengthPadel = 12.5;
    
    double min_abs_m = (yOffsetCRT - yOffsetPadel)/(semiTotLengthCrt + semiLengthPadel);

    double max_abs_th = TMath::ACos(1/TMath::Sqrt(1+pow(min_abs_m, -2))) + 0.01; //10 mrad di sicurezza
    double max_abs_xgen = yGen / min_abs_m - semiLengthPadel + 2; // 10 cm di sicurezza

    double m, xGen, th, en;

    double xGenMin=-100, xGenMax=100;

    int RowTmp, ColTmp, ChaTmp, tmp, row[20], col[20];
    FILE *mapFile;
    mapFile = fopen("mod0_dirac_map.dat", "r");

    for ( int iCh=0; iCh<20; iCh++ ){
        fscanf(mapFile, "%d %d %d %d %d", &tmp, &ChaTmp, &RowTmp, &ColTmp, &tmp );
        row[ChaTmp] = RowTmp;
        col[ChaTmp] = ColTmp;
    }

    TF1 *th_pdf = new TF1("th_pdf", "2/ TMath::Pi() * cos(x) * cos(x)", -max_abs_th, max_abs_th);
    
    TH2F *h;
    TLine *line;
    TCanvas *c;
    TPad *plots_pad;
    TButton *but1;

    TFile *f;
    TTree *tree;
    Int_t           evnum_out = 0;
    Int_t           nHits_out;
    Int_t           iDAQ_out[20];
    Double_t        Qval_out[20];
    Int_t           crtBar_out;

    if(plot){
        c = new TCanvas("c", "c");
        plots_pad = new TPad("pad1", "pad1",0.0, 0.1, 1.0, 1.0);
        TPad *pad1 = new TPad("pad2", "pad2",0.0, 0, 1.0, 0.1);
        plots_pad->Draw();
        pad1->Draw();
        pad1->cd();
        but1 = new TButton("Next", "go_ahead()", 0.45, 0.15, 0.55, 0.95);
        but1->Draw();
        plots_pad->cd();
        line = new TLine(xGenMin, yGen, xGenMax, yGen);
        h = new TH2F("tmp", "tmp", 1000, -50, 50, 1000, -10, 100);
        h->Draw();
        line->Draw();
    }
    else{
        f = new TFile("sim.root", "RECREATE");
        tree = new TTree("mod0","mod0");   
        tree->SetAutoSave(1000);  
        tree->Branch("evnum"     ,&evnum_out,     "evnum/I");
        tree->Branch("nHits"     ,&nHits_out,     "nHits/I");
        tree->Branch("iDAQ"      ,&iDAQ_out,      "iDAQ[nHits]/I");
        tree->Branch("Qval"      ,&Qval_out,      "Qval[nHits]/D");
        tree->Branch("crtBar"    ,&crtBar_out,    "crtBar/I");
    }

    Scintillator *crtBar[8];
    // inizializzazione barre crt
    for(int i=0; i<8; i++){
        crtBar[i] = new Scintillator(Form("crtBar_%i", i), - semiTotLengthCrt + 1.25 + 2.5*i, 2.5, yOffsetCRT, 1.5, 0);
        if(plot){
            auto box = crtBar[i]->GetBox();
            box->SetFillColor(kRed);
            box->SetLineColor(kBlack);
            box->SetLineWidth(1.);
            box->Draw("l same");
        }
    }

    // inizializzazione cristalli mod0
    Scintillator *crystal[20];
    double slope[20] = {0}; // come facciamo???
    for(int i=0; i<20; i++){

        RowTmp = row[i];
        ColTmp = col[i];

        crystal[i] = new Scintillator(Form("mod0cry_%i", i), 0, 20, 20 + ycry[RowTmp][ColTmp]/10, 3.4, slope[i]);
        if(plot){
            auto box = crystal[i]->GetBox();
            box->SetFillColor(kBlue);
            box->SetLineColor(kBlack);
            box->SetLineWidth(1.);
            box->Draw("l same");
        }
    }

    Scintillator *padel = new Scintillator("paletta" ,0, 2*semiLengthPadel, yOffsetPadel, 5, 0);
    if(plot){
            auto box = padel->GetBox();
            box->SetFillColor(kGreen);
            box->SetLineColor(kBlack);
            box->SetLineWidth(1.);
            box->Draw("l same");
        }


    for(int iEv=0; iEv<nev; iEv++){
        xGen = gRandom->Uniform(-max_abs_xgen, max_abs_xgen);

        th = th_pdf->GetRandom();

        m = 1/TMath::Tan(th);

        TGraph toDraw;
        if (plot){
            toDraw.SetMarkerStyle(8);
            toDraw.SetMarkerSize(0.5);
        }

        int nHitBar = 0, hitBar = 99;
        for(int i=0; i<8; i++){
            // si può mettere soglia sull'energia depositata
            auto scint = crtBar[i];
            scint->ProcessEvent(1, xGen, yGen, m, toDraw);
            if (scint->isHit){
                hitBar = i;
                nHitBar++;
            }
        }
        if (nHitBar != 1) continue;

        int iHitMod0 = 0;
        for(int i=0; i<20; i++){
            auto scint = crystal[i];
            scint->ProcessEvent(0, xGen, yGen, m, toDraw);
            if (scint->isHit){
                en = scint->en;
                iDAQ_out[iHitMod0] = i;
                Qval_out[iHitMod0] = en * 173; // pC/MeV - misurato da MPV di Qval run98 iDAQ==0 sapendo che il deposito è 21 MeV
                iHitMod0++;
            }
        }

        if (iHitMod0==0) continue;

        padel->ProcessEvent(1, xGen, yGen, m, toDraw);
        if (! padel->isHit) continue;                           

        if(plot){

            TF1 l("l", "[0] - (x - [1]) * [2]", -100, 100);
            l.SetParameter(0, yGen);
            l.SetParameter(1, xGen);
            l.SetParameter(2, m);
            l.Draw("same");

            toDraw.Draw("p");

            while (go_ahead_flag == 0)
            {
                gPad->Modified();
                gPad->Update();
                gSystem->ProcessEvents();
                usleep(10000);
            }
            go_ahead_flag = 0;
        }
        else{
            nHits_out = iHitMod0;
            crtBar_out = hitBar;
            evnum_out++;
            tree->Fill();
        }
    }

    tree->Write();
    f->Close();

    cout << "Generation efficiency: " << evnum_out/(double)nev << endl;
    //exit(0); uncomment to measure execution time with 'time root -x ...'
}


int main(){
    TApplication *myapp = new TApplication("myapp", 0, 0);

    Qsimul_CRTandScintillatorvsMod0();

    myapp->Run(true);
}
