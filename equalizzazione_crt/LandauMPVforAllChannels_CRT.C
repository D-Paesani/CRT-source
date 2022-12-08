double MPV(TH1* histObj) {

  gStyle->SetOptFit(1);

  TSpectrum s(3);
  s.Search(histObj, 1); // a volte ne servono 2 in autotrigger
  double *pks = s.GetPositionX();
  double qpeak = *std::min_element(pks, pks + 2);

  //double min_tmp = histObj->GetXaxis()->GetXmin();
  //double max_tmp = histObj->GetXaxis()->GetXmax();
  //histObj->GetXaxis()->SetRangeUser(390,700);
  //double qpeak = histObj->GetBinCenter(histObj->GetMaximumBin());

  double qmax = qpeak + 200, qmin = qpeak-70;
  float pk,sigma;

  TF1 l1 = TF1("l", "landau", qmin, qmax);
  l1.SetParameters(histObj->Integral()/2, qpeak, 100);
  histObj->Fit(&l1, "RQ");

  pk = l1.GetMaximumX(); sigma = l1.GetParameter(2);
  TF1 l2 = TF1("l", "landau", pk-80, pk+2*sigma);
  l2.SetParameters(l1.GetParameter(0), l1.GetParameter(1), sigma);
  histObj->Fit(&l2, "RQ");

  pk = l2.GetParameter(1); sigma = l2.GetParameter(2);
  TF1 l3 = TF1("l", "landau", pk-1.8*sigma, pk+2.5*sigma);
  l3.SetParameters(l2.GetParameter(0), l2.GetParameter(1), sigma);

  histObj->Fit(&l3, "RQ");

  return l3.GetParameter(1);

}

void LandauMPVforAllChannels_CRT(){
  auto *f = new TFile("../../data/step3/run358_T_s3.root");

  TH1F *charge[8];
  for(int i=0; i<8; i++) charge[i] = new TH1F(Form("q%i", i), Form("q%i", i), 200, 150, 950);

  auto *c = new TCanvas("c", "c");
  c->Divide(2, 4);

  auto *tree = (TTree*)f->Get("CRT");

  ofstream of("landau_crt_T_358.csv");

  for(int i=0; i<8; i++){
    //tree->Draw(Form("Q[%i]>>q%i", i%2, i), Form("iSc==%i && Z > -20 && Z < 20", 2*(i/2)), "goff");
    tree->Draw(Form("Q[%i]>>q%i", i%2, i), Form("iSc==%i", 2*(i/2)), "goff");
    c->cd(i+1);
    charge[i]->Draw();
    of << i%2 << " " << 2*(i/2) << " " << MPV(charge[i]) << endl;

  }

  of.close();

}
