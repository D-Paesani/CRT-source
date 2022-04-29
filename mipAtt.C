
#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TF1.h"
#include "TCanvas.h"

void add_from_file(TGraphErrors *g, TString filename, TString canvasname, double qmin, double qmax, double slice_w, int checkfit){
  auto *inf = new TFile(filename);
  auto *outf = new TFile("mioattout.root", "recreate");
  auto *c = (TCanvas*)inf->Get(canvasname);
  c->Draw();
  auto *h = (TH2F*)c->GetListOfPrimitives()->At(1)->Clone();
  h->Draw("zcol");

  int n_slices = -(int)((qmin-qmax)/slice_w);

  auto **proj = new TH1F*[n_slices];
  for(int i=1; i<n_slices-3; i++){
    h->SetAxisRange(qmin + slice_w*i, qmin + slice_w*(i+1));
    proj[i] = (TH1F*)h->ProjectionY()->Clone();
    auto *tempx = (TH1F*)h->ProjectionX();
    double peak = proj[i]->GetMean(), sigma = proj[i]->GetRMS();
    TF1 f("f", "landau", peak-3*sigma, peak+3*sigma);
    proj[i]->Fit(&f, "RQ");
    outf->cd();
    peak = f.GetParameter(1); sigma = f.GetParameter(2);
    f = TF1("f", "landau", peak-sigma, peak+2*sigma);
    proj[i]->Fit(&f, "RQ");

    peak = f.GetParameter(1); sigma = f.GetParameter(2);
    f = TF1("f", "landau", peak-1.5*sigma, peak+2*sigma);

    proj[i]->Fit(&f, "RQ");
    proj[i]->Write();

    if( ( (f.GetProb() > 0.05) || (!checkfit) ) && f.GetParError(1) < 100 ){
      g->AddPoint(tempx->GetMean(), f.GetParameter(1));
      g->SetPointError(g->GetN()-1, slice_w/3, f.GetParError(1));
    }
    else{
      cout << "Fit failed from :" << qmin + slice_w*i << "to " << qmin + slice_w*(i+1) << " pC" << endl;
    }
  }
  inf->Close();
  outf->Close();
}

void mipAtt(){

  auto *g = new TGraphErrors();

  add_from_file(
    g, "/home/ruben/Documents/frascati/CRT-analysis/data/plot/"
    "zqmip1.root",
    "c1",
    -100, 100, 10, 1
  );

/*  add_from_file(
    g, "/home/ruben/Documents/frascati/CRT-analysis/data/plot/"
    "zq_183.root",
    "c1",
    400, 700, 150, 0
  );*/

  auto *c1 = new TCanvas("sz_q_canvas", "c");
  c1->cd();
  g->Draw("AP");

  TF1 func("f", "[0] * ( TMath::Exp((x-80)/380) + [1]*TMath::Exp((x-80)/[2]) )", -80, 70);
  func.SetParameters(100, 3, 20);
  func.SetParLimits(1, 0.1, 10);
  func.SetParLimits(2, 1, 50);
  func.SetParName(0, "Normalization");
  func.SetParName(1, "L_{I}/L_{D}");
  func.SetParName(2, "TAL");
  for(int i=0; i<3; i++) g->Fit(&func, "R");
  gStyle->SetOptFit(1);

}

int main(){
  mipAtt();
}
