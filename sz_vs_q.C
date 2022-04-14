#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TF1.h"
#include "TCanvas.h"

void sz_vs_q(){
  TString filename(
    "/home/ruben/Documents/frascati/CRT-analysis/data/plot/"
    "zq_stronzio_centro.root" //oppure _lato0.root
  );
  auto *f = new TFile(filename);
  auto *c = (TCanvas*)f->Get("Canvas_1");
  //auto *c = (TCanvas*)f->Get("c1");
  c->Draw();
  auto *h = (TH2F*)c->GetListOfPrimitives()->At(1)->Clone();
  h->Draw("zcol");

  int n_slices = 21;
  double qmin=100, qmax=200;
  double slice_w = (qmax - qmin)/n_slices;
  auto **proj = new TH1F*[n_slices];
  auto *out_f = new TFile("sz_vs_q.root", "recreate");
  out_f->cd();
  auto *g = new TGraphErrors();
  for(int i=0; i<n_slices; i++){
    h->SetAxisRange(qmin + slice_w*i, qmin + slice_w*(i+1));
    proj[i] = (TH1F*)h->ProjectionY()->Clone();
    double mean = proj[i]->GetMean(), sigma = proj[i]->GetRMS();
    TF1 f("f", "gaus", mean-2*sigma, mean+2*sigma);
    proj[i]->Fit(&f, "RQ");
    proj[i]->Write(Form("z_q%.0f", qmin + slice_w*(i+0.5)));
    if(f.GetProb() > 0.05 || 1){
      g->AddPoint((1.28)*(qmin + slice_w*(i+0.5)), 2*f.GetParameter(2)/12.2);
      g->SetPointError(g->GetN()-1, slice_w/3 * 1.28, 2*f.GetParError(2)/12.2);
    }
  }
  g->AddPoint(500, 0.220);
  g->SetPointError(g->GetN()-1, 15, 0.002);
  TF1 func("f", "3.5*TMath::Sqrt(2)/TMath::Sqrt([0]*x)", 50, 600);
  func.SetParName(0, "NPE_per_pC");
  func.SetParameters(1, 1, 1);
  gStyle->SetOptFit();
  auto *c1 = new TCanvas("sz_q_canvas", "c");
  c1->cd();
  for(int i=0; i<3; i++) g->Fit(&func, "R");
  g->Draw("AP");
  g->Write("sz_q");
  c1->Write();
  cout << func.Eval(480) << endl;

}

int main(){
  sz_vs_q();
}
