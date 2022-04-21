
#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TF1.h"
#include "TCanvas.h"

void add_from_file(TGraphErrors *g, TString filename, TString canvasname, double qmin, double qmax, double slice_w, int checkfit){
  auto *f = new TFile(filename);
  auto *c = (TCanvas*)f->Get(canvasname);
  c->Draw();
  auto *h = (TH2F*)c->GetListOfPrimitives()->At(1)->Clone();
  h->Draw("zcol");

  int n_slices = -(int)((qmin-qmax)/slice_w);

  auto **proj = new TH1F*[n_slices];
  for(int i=0; i<n_slices; i++){
    h->SetAxisRange(qmin + slice_w*i, qmin + slice_w*(i+1));
    proj[i] = (TH1F*)h->ProjectionY()->Clone();
    auto *tempx = (TH1F*)h->ProjectionX();
    double mean = proj[i]->GetMean(), sigma = proj[i]->GetRMS();
    TF1 f("f", "gaus", mean-2*sigma, mean+2*sigma);
    proj[i]->Fit(&f, "RQ");
    if((f.GetProb() > 0.05) || (!checkfit) ){
      g->AddPoint(tempx->GetMean(), 2*f.GetParameter(2)/12.5 / TMath::Sqrt(2) );
      g->SetPointError(g->GetN()-1, tempx->GetRMS(), 2*f.GetParError(2)/12.5 / TMath::Sqrt(2));
    }
    else{
      cout << "Fit failed from :" << qmin + slice_w*i << "to " << qmin + slice_w*(i+1) << " pC" << endl;
    }
  }
  f->Close();
}

void sz_vs_q(){

  auto *g = new TGraphErrors();

  add_from_file(
    g, "/home/ruben/Documents/frascati/CRT-analysis/data/plot/"
    "zq_stronzio_centro.root",
    "c1",
    150, 300, 3, 0
  );

/*  add_from_file(
    g, "/home/ruben/Documents/frascati/CRT-analysis/data/plot/"
    "zq_183.root",
    "c1",
    400, 700, 150, 0
  );*/

  TF1 func("f", "TMath::Sqrt( [0]*[0] + [1]*[1]/x + [2]*[2]/(x*x))", 50, 1000);
  func.SetParameters(0.1, 5, 50);
  func.SetParLimits(0, 0, 1);
  func.SetParLimits(1, 0, 10);
  func.SetParLimits(2, 0, 100);
  func.SetParName(0, "a");
  func.SetParName(1, "b");
  func.SetParName(2, "c");

  gStyle->SetOptFit();
  auto *c1 = new TCanvas("sz_q_canvas", "c");
  c1->cd();
  for(int i=0; i<3; i++) g->Fit(&func, "R");
  g->Draw("AP");
  cout << func.Eval(480) << endl;

}

int main(){
  sz_vs_q();
}
