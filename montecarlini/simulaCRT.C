#include <fstream>
#include "TTree.h"
#include "TF1.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TH1.h"
#include "TFile.h"
#include <stdlib.h>
#include <iostream>
#include "TH2D.h"
#include "TCanvas.h"
#include "TMath.h"
#include <string>
#include "TMath.h"
#include "TH1.h"
#include "TSpectrum.h"
#include "TVirtualFitter.h"
#include "TTime.h"
#include "TLegend.h"
#include "TLine.h"
#include "TStyle.h"
#include <TDatime.h>
#include "TRandom.h"
#include "TRandom3.h"
#include "TRandom2.h"
#include "TBenchmark.h"


using namespace std;



void simulaCRT(){

  float NPEA = 30000; // INVENTATO
  float resmip = 0.1; // INVENTATO
  int nevents = 1000000;
  cout << " Npe A " << NPEA << endl;
  cout << " Mip Core resolution @ center " << resmip << endl;

   auto form2 = new TF1("form2","[0]*TMath::Landau(x,[1],[2])",0,100000);
   float Norm = 100;
   float MaxProb = NPEA;
   float Tail = sqrt(resmip*MaxProb);
   form2->SetParameters(Norm,MaxProb,Tail);

   auto attFunc = new TF1("attFunc", "[0]*(exp(-x/[1]) + [3]*exp(-x/[2]))", 0, 160);
   float lBAL = 380;
   float lTAL = 30;
   float Ir_over_Id = 0.8;
   float p0 = 1/(1+Ir_over_Id);

   attFunc->SetParameters(p0, lBAL, lTAL, Ir_over_Id);

   TRandom2 *r = new TRandom2();
   float convFact = 0.1; // INVENTATO
   float cut = 1150;// INVENTATO
   // now convolute with a Gaussian "Spread" of NPEA 
   auto PEnoCut = new TH1F("PEnoCut","Photoelectrons on A after light transport (no cut)",500,0,4000);
   auto PE = new TH1F("PE","Photoelectrons on A after light transport (> 1150)",500,0,4000);

   auto Z = new TH1F("Z", "Z", 500, -100, 100);

   for (int jj =0; jj<nevents; jj++){
     float NPval = form2->GetRandom(); //n photons from scintillation
     float zval = r->Uniform(-80,80); // uniform in cm
     float NPval_att = NPval*attFunc->Eval(zval+80);
     float NPE = NPval_att * convFact;
     float smearval = sqrt(NPE); // std of Poissonian
     float gauss_NPE = gRandom->Gaus(0,smearval); // Poissonian
     NPE += gauss_NPE;
     PEnoCut->Fill(NPE);
     if(NPE > cut){
       PE->Fill(NPE);
       Z->Fill(zval);
     }
   }


   TCanvas *c1 = new TCanvas("c1", "simulaCRT",200,10,900,900);
    
   auto padPE = new TPad("pad1","The pad with PE histogram", 0.05, 0.5, 0.45, 0.95);
   auto padPEnoCut = new TPad("pad2","The pad with cutted PE histogram", 0.5, 0.5, 0.95, 0.95);  
   auto padZ = new TPad("pad3","The pad with Z",  0.05, 0.05, 0.95, 0.45);
   padPE->Draw();
   padPEnoCut->Draw();
   padZ->Draw();
    
   padPE->cd();
   PE->Draw("e");
   padPEnoCut->cd();
   PEnoCut->Draw("e");
   padZ->cd();
   Z->Draw("e");

    //
   // in a ROOT file and save the formula, function and histogram
   //
   TFile myfile("simulaCRT.root","RECREATE");

   //form2->Write();
   PE->Write();
   PEnoCut->Write();
   Z->Write();
   //h2f->Write();
   gBenchmark->Show("simulaCRT");
}
