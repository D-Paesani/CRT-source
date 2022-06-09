#include "../../CRT-analysis/CRT-source/includes/HistManager.h"
#include "module0_map.h"
#include <TSpectrum.h>
#include <TLine.h>

HistManager HM;
double landau_mu_ch_crt[20][8], landau_mu_err_ch_crt[20][8];
double landau_mu_ch[20]={0}, landau_mu_err_ch[20]={0};

void NamerMatrix(int hN, TString& hTag, TString& hTitleTag) {
    hTag = Form("_%d_%d", hN%20, hN/20);
    hTitleTag = Form("[sd%d][sc%d] ", hN%20, hN/20);
}

void q_proc(TH1* hist, int histN) {
  gErrorIgnoreLevel = kError;
  TSpectrum s(1);
  s.Search(hist, 2, "nodraw");
  gErrorIgnoreLevel = kWarning;
  double qpeak = *(s.GetPositionX());

  double qmax = qpeak + 1200, qmin = qpeak - 600;

  TF1 l1 = TF1("l", "landau", qmin, qmax);
  l1.SetParameters(hist->Integral()/2, qpeak, 1000);
  hist->Fit(&l1, "Rq");
  float pk = l1.GetMaximumX(), sigma = l1.GetParameter(2);

  pk = l1.GetParameter(1); sigma = l1.GetParameter(2);
  TF1 l3 = TF1("l", "landau", pk-2*sigma, pk+6*sigma);
  l3.SetParameters(l1.GetParameter(0), l1.GetParameter(1), sigma);
  hist->Fit(&l3, "Rq");

  pk = l3.GetParameter(1); sigma = l3.GetParameter(2);
  TF1 l4 = TF1("l", "landau", pk-1.5*sigma, pk+4*sigma);
  l3.SetParameters(l3.GetParameter(0), l3.GetParameter(1), sigma);
  hist->Fit(&l4, "Rq");

  int iDaq = histN%20, crtBar = histN/20;
  if (TString(hist->GetName()).Contains("_crt")){
    landau_mu_ch_crt[iDaq][crtBar] = l4.GetParameter(1);
    landau_mu_err_ch_crt[iDaq][crtBar] = l4.GetParError(1);
  }
}

void createHistBoxes() {
  HM.AddHistBox("th1f", 20, "q_ch", "Charge", "Q", "pC", 400, 400, 20400);
  HM.AddHistBox("th1f", 20*8, "q_ch_crt","Charge", "Q", "pC", 200, 400, 20400);
  HM.AddHistBox("th2f", 1, "q_norm_crt", "Normalized peak charge vs. CRT bar", "crtBar", "", "Charge / Mean Charge", "", 8, 0, 8, 20, 0.95, 1.05);
}

void lru(){

  int row[20], col[20], BoaTmp, ChaTmp, RowTmp, ColTmp, SipmTmp;
  double minxcry, maxxcry, minycry, maxycry;
  FILE *mapFile;
  mapFile = fopen("mod0_dirac_map.dat", "r");
  for ( int Iloop=0; Iloop<20; Iloop++ ){
    fscanf(mapFile, "%d %d %d %d %d", &BoaTmp, &ChaTmp, &RowTmp, &ColTmp, &SipmTmp );
    row[ChaTmp] = RowTmp;
    col[ChaTmp] = ColTmp;

    if(Iloop==0){
      minxcry = xcry[RowTmp][ColTmp];
      maxxcry = xcry[RowTmp][ColTmp];
      minycry = ycry[RowTmp][ColTmp];
      maxycry = ycry[RowTmp][ColTmp];
    }

    if (xcry[RowTmp][ColTmp] < minxcry) minxcry = xcry[RowTmp][ColTmp];
    if (xcry[RowTmp][ColTmp] > maxxcry) maxxcry = xcry[RowTmp][ColTmp];
    if (ycry[RowTmp][ColTmp] < minycry) minycry = ycry[RowTmp][ColTmp];
    if (ycry[RowTmp][ColTmp] > maxycry) maxycry = ycry[RowTmp][ColTmp];
  }
  fclose(mapFile);
  double x_canvas_coord_norm = maxxcry + wcry - minxcry;
  double y_canvas_coord_norm = maxycry + wcry - minycry;

  //TString outFileName = "out_lru.root";
  TString outFileName = "out_lru_sim.root";

  TFile *outFile = new TFile(outFileName, "RECREATE");
  cout << "Writing to file: " << outFileName.Data() << endl;

  HM.SetOutFile(outFile);
  HM.SetNamerFun(&NamerMatrix);

  createHistBoxes();

  //TString inFileName = "run_98_all_out.root";
  TString inFileName = "sim.root";
  TFile *inFile = new TFile(inFileName);
  cout << "Reading from file: " << inFileName.Data() << endl;

  TTree *tree = (TTree*)inFile->Get("mod0");

  int histN;
  double x2[20*8], y2[20*8], y2_err[20*8], w;
  TH1* tmp;

  for(int iCh=0; iCh<20; iCh++){
    tmp = HM.GetHist("q_ch", iCh);
    tree->Draw(Form("Qval>>h_%i(400, 400, 20400)", iCh), Form("iDAQ==%i", iCh), "goff");
    tmp->Add((TH1F*)gDirectory->Get(Form("h_%i", iCh)));
    q_proc(tmp, iCh);

    RowTmp = row[iCh];
    ColTmp = col[iCh];
    for(int iBar=0; iBar<8; iBar++){
      histN = iCh+20*iBar;

      tmp = HM.GetHist("q_ch_crt", histN);
      tree->Draw(Form("Qval>>h_%i_%i(200, 400, 20400)", iCh, iBar), Form("iDAQ==%i && crtBar==%i", iCh, iBar), "goff");
      tmp->Add((TH1F*)gDirectory->Get(Form("h_%i_%i", iCh, iBar)));
      q_proc(tmp, histN);

      w = 1/(landau_mu_err_ch_crt[iCh][iBar]*landau_mu_err_ch_crt[iCh][iBar]);
      landau_mu_ch[iCh] += landau_mu_ch_crt[iCh][iBar] * w;
      landau_mu_err_ch[iCh] += w;


      x2[histN] = (iCh==14) ? 1000 : iCh + 0.125*iBar;
      y2[histN] = landau_mu_ch_crt[iCh][iBar];

      y2_err[histN] = landau_mu_err_ch_crt[iCh][iBar];

    }
    landau_mu_ch[iCh] /= landau_mu_err_ch[iCh];
    landau_mu_err_ch[iCh] = pow(landau_mu_err_ch[iCh], -0.5);
  }

  double x1[8] = {0, 1, 2, 3, 4, 5, 6, 7}, y1[20][8], y1_err[20][8], q_norm, q_norm_err;

  for(int iCh=0; iCh<20; iCh++){
    if(iCh==14 || iCh==13 || iCh==7) continue;

    for(int iBar=0; iBar<8; iBar++){

      q_norm = landau_mu_ch_crt[iCh][iBar]/landau_mu_ch[iCh];
      q_norm_err = q_norm * TMath::Sqrt(
        pow(landau_mu_err_ch_crt[iCh][iBar]/landau_mu_ch_crt[iCh][iBar], 2)
          +
        pow(landau_mu_err_ch[iCh]/landau_mu_ch[iCh], 2)
      );

      if(iCh!=14){
        ((TH2F*)HM.GetHist("q_norm_crt", 0))->Fill(iBar, q_norm, 1/(q_norm_err*q_norm_err));
      }
      
      y1[iCh][iBar] = q_norm;
      y1_err[iCh][iBar] = q_norm_err;
    }
  }

  cout << "Processing plots" << endl;
  HM.ProcessBoxes();

  TProfile *q_norm_crt_prof = ((TH2F*)HM.GetHist("q_norm_crt", 0))->ProfileX();
  q_norm_crt_prof->SetName("q_norm_crt_prof");
  q_norm_crt_prof->SetTitle("Normalized charge vs. crtBar");
  outFile->cd();
  q_norm_crt_prof->Write();

  TCanvas *c1 = new TCanvas("q_graph1", "Q vs crtBar vs CH");

  gStyle->SetTitleFontSize(0.08);
  for(int iCh=0; iCh<20; iCh++){
    if(iCh==14 || iCh==13 || iCh==7) continue;

    RowTmp = row[iCh];
    ColTmp = col[iCh];
    double xlow = (xcry[RowTmp][ColTmp] - minxcry)/x_canvas_coord_norm;
    double ylow = (ycry[RowTmp][ColTmp] - minycry)/y_canvas_coord_norm;
    double xhigh = (xcry[RowTmp][ColTmp] + wcry - minxcry)/x_canvas_coord_norm;
    double yhigh = (ycry[RowTmp][ColTmp] + wcry - minycry)/y_canvas_coord_norm - 0.01;

    c1->cd();

    TPad *p = new TPad(Form("p%i", iCh), "", xlow, ylow, xhigh, yhigh);
    p->Draw();
    p->cd();

    TGraphErrors *q_graph = new TGraphErrors(8, x1, y1[iCh], 0, y1_err[iCh]);
    q_graph->SetMarkerStyle(22);
    q_graph->SetTitle(Form("CH%i Norm. Charge vs. CRT bar", iCh));
    auto xaxis = q_graph->GetXaxis();
    xaxis->SetRangeUser(-1, 8);
    xaxis->SetLabelSize(0.1);
    xaxis->SetLabelOffset(0.02);

    auto yaxis = q_graph->GetYaxis();
    yaxis->SetTitle("Normalized Charge");
    yaxis->SetTitleSize(0.08);
    yaxis->SetLabelSize(0.1);    
    q_graph->Draw("AP");
  }

  outFile->cd();
  c1->Write();


  TGraphErrors *q_graph = new TGraphErrors(20*8, x2, y2, 0, y2_err);

  TCanvas *c2 = new TCanvas("q_graph2", "Q vs crtBar vs CH");
  c2->cd();
  q_graph->SetName("q_graph2");
  q_graph->SetTitle("Q vs crtBar vs CH");
  q_graph->GetXaxis()->SetTitle("DAQ channel + (CRT bar / 8)");
  q_graph->GetYaxis()->SetTitle("Charge landau par. [pC]");
  q_graph->SetMarkerStyle(22);
  q_graph->GetXaxis()->SetRangeUser(0, 19);

  q_graph->Draw("AP");
  for (int iBar=0; iBar<8; iBar++){
    for(int iCh=0; iCh<20; iCh++){
      histN = iCh+20*iBar;
      TLatex*t = new TLatex(x2[histN], y2[histN] + y2_err[histN]+100, Form("%i", iBar));
      t->SetTextAlign(22);
      t->SetTextColor(iCh%2 ? kRed : kBlue);
      t->SetTextSize(0.015);
      c2->cd();
      t->Draw();
    }
  }
  outFile->cd();
  c2->Write();

  outFile->Close();
}

int main(){
  lru();
}
