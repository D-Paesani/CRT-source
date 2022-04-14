#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TF1.h"
#include "TCanvas.h"

void asymm_vs_q(){
  TString filename(
    "/home/ruben/Documents/frascati/CRT-analysis/data/plot/"
    "asymm_q_stronzio_centro.root" //oppure _lato0.root
  );
  auto *f = new TFile(filename);
  auto *c = (TCanvas*)f->Get("c1");
  c->Draw();
  auto *h = (TH2F*)c->GetListOfPrimitives()->At(1)->Clone();
  h->Draw("zcol");

  int n_slices = 20;
  double qmin=150, qmax=350;
  double slice_w = (qmax - qmin)/n_slices;
  auto **proj = new TH1F*[n_slices];
  auto *out_f = new TFile("sasymm_vs_q.root", "recreate");
  out_f->cd();
  auto *g = new TGraphErrors();
  for(int i=0; i<n_slices; i++){
    h->SetAxisRange(qmin + slice_w*i, qmin + slice_w*(i+1));
    proj[i] = (TH1F*)h->ProjectionY()->Clone();
    double mean = proj[i]->GetMean(), sigma = proj[i]->GetRMS();
    TF1 f("f", "gaus", mean-1.5*sigma, mean+1.5*sigma);
    proj[i]->Fit(&f, "R");
    mean = f.GetParameter(1), sigma = f.GetParameter(2);
    f = TF1("f", "gaus", mean-1.5*sigma, mean+1.5*sigma);
    proj[i]->Fit(&f, "R");
    proj[i]->Write(Form("asymm_q%.0f", qmin + slice_w*(i+0.5)));
    if(f.GetProb() > 0.02){
      g->AddPoint(qmin + slice_w*(i+0.5), f.GetParameter(2));
      g->SetPointError(g->GetN()-1, slice_w/3, f.GetParError(2));
    }
  }
  TF1 func("f", "TMath::Sqrt([0]*[0] + 1/([1]*x))/TMath::Sqrt(2)", 50, 600);
  func.SetParName(0, "Constant Term");
  func.SetParName(1, "NPE_per_pC");
  func.SetParameters(1, 1);
  gStyle->SetOptFit();
  auto *c1 = new TCanvas("asymm_q_canvas", "c");
  c1->cd();
  for(int i=0; i<3; i++) g->Fit(&func, "R");
  g->Draw("AP");
  g->Write("asymm_q");
  c1->Write();
}

int main(){
  asymm_vs_q();
}
