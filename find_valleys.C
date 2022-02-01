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

#include "TFile.h"
#include "TTree.h"

#include "includes/HistManager.h"
#include "includes/NumberingHelper.h"
#include "includes/AnaPars.h"

using namespace std;

int nslices = 15;
double zmin = -80, zmax=80;
double slice_width = (zmax-zmin)/nslices;
int slicestart = 2, sliceend = 13;
int ngoodslices = sliceend - slicestart;

HistManager HM;

void createHistBoxes() {
    
  for(int iSl=slicestart; iSl<sliceend; iSl++){

    double curr_z_low_edge = zmin + iSl*slice_width;
    double curr_z_high_edge = curr_z_low_edge + slice_width;
    HM.AddHistBox(
      "th1f", 2*scintNum,
      Form("QmipZslices_%i", iSl), Form("Qmip per %.0f cm < Z < %.0f cm", curr_z_low_edge, curr_z_high_edge),
      "Q" ,"pC", qBins/2, qFrom, qTo
    );
  }
}

void find_valleys(){
  TFile *f = new TFile("../data/step3/run205_s3.root");

  TFile *outf = new TFile("../data/step3/run205_valleys.root", "recreate");

  HM.SetOutFile(outf);
  HM.SetNamerFun(&NamerMatrix);
  createHistBoxes();

  double curr_z_low_edge, curr_z_high_edge;

  for(int iSc = 0; iSc < 8; iSc ++){
    for(int iSd = 0; iSd < 2; iSd++){

      double valleys[ngoodslices], valleys_err[ngoodslices], zetas[ngoodslices], zetas_err[ngoodslices];
      std::fill_n(zetas_err, ngoodslices, 0);

      TH2F *curr_qz_hist = (TH2F*)f->Get(Form("zeta_q/zeta_q_%i_%i", iSd, iSc));

      for(int iSl = slicestart; iSl < sliceend; iSl++){

        TH1F *curr_q_hist = (TH1F*)HM.GetHist(Form("QmipZslices_%i", iSl), GetChan(iSd, iSc));

        curr_z_low_edge = zmin + iSl*slice_width;
        curr_z_high_edge = curr_z_low_edge + slice_width;
        curr_qz_hist->GetXaxis()->SetRangeUser(curr_z_low_edge, curr_z_high_edge);
        curr_q_hist->Add(curr_qz_hist->ProjectionY());
        curr_q_hist->Smooth(6);
        TSpectrum smax(2);
        smax.Search(curr_q_hist, 5, "noMarkov");
        double *pks = smax.GetPositionX();
        double qmin = *std::min_element(pks, pks + 2) + 60;
        double qmax = *std::max_element(pks, pks + 2) - 60;

        TF1 g = TF1("g", "gaus + [3]", qmin, qmax);
        g.SetParameters(-30, (qmin+qmax)/2, (qmax-qmin)/4, curr_q_hist->GetMaximum());
        g.SetParLimits(0, -10, -1000);
        g.SetParLimits(2, 10, 300);
        cout << curr_q_hist->GetMinimum() << endl;
        curr_q_hist->Fit(&g, "R");
        valleys[iSl-slicestart] = g.GetParameter(1);
        valleys_err[iSl-slicestart] = g.GetParError(1);
        zetas[iSl-slicestart] = (curr_z_low_edge + curr_z_high_edge)/2;
        curr_q_hist->Write();
      }

      TCanvas *c = new TCanvas(Form("c%i%i", iSd, iSc));
      curr_qz_hist->GetXaxis()->SetRangeUser(-160, 160);

      curr_qz_hist->Draw("zcol");
      TGraphErrors *ge = new TGraphErrors(ngoodslices, zetas, valleys, zetas_err, valleys_err);
      ge->SetMarkerColor(kRed);
      ge->SetLineColor(kRed);
      ge->SetMarkerStyle(22);
      ge->SetMarkerSize(1);
      ge->SetName(Form("ge%i%i", iSd, iSc));
      ge->Draw("p same");
      ge->Write();
      c->Write();
    }
  }

  HM.CloseOutFile();

}
