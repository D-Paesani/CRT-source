double MPV(TH1* histObj) {   

  gStyle->SetOptFit(1); 

  TSpectrum s(2);
  s.Search(histObj, 1); // a volte ne servono 2 in autotrigger
  double *pks = s.GetPositionX();
  double qpeak = *std::max_element(pks, pks + 2);

  //double min_tmp = histObj->GetXaxis()->GetXmin();
  //double max_tmp = histObj->GetXaxis()->GetXmax();
  //histObj->GetXaxis()->SetRangeUser(390,700);
  //double qpeak = histObj->GetBinCenter(histObj->GetMaximumBin());

  double qmax = qpeak + 300, qmin = qpeak-100;
  float pk,sigma;

  TF1 l1 = TF1("l", "landau", qmin, qmax);                 
  l1.SetParameters(histObj->Integral()/2, qpeak, 200);                
  histObj->Fit(&l1, "RQ"); 

  pk = l1.GetMaximumX(); sigma = l1.GetParameter(2);
  TF1 l2 = TF1("l", "landau", pk-80, pk+3*sigma);          
  l2.SetParameters(l1.GetParameter(0), l1.GetParameter(1), sigma);  
  histObj->Fit(&l2, "RQ");

  pk = l2.GetParameter(1); sigma = l2.GetParameter(2); 
  TF1 l3 = TF1("l", "landau", pk-0.8*sigma, pk+4*sigma);     
  l3.SetParameters(l2.GetParameter(0), l2.GetParameter(1), sigma);  
  
  
  histObj->Fit(&l3, "RQ");
 
  return l3.GetParameter(1);
}

void LandauMPVforAllChannels(){
  auto *f = new TFile("../../../data/dirac/step2/run_193_all.root");

  TH1F *charge[20];
  for(int i=0; i<20; i++) charge[i] = new TH1F(Form("q%i", i), Form("q%i", i), 400, 100, 900);

  auto *c = new TCanvas("c", "c");
  c->Divide(5, 4);

  auto *tree = (TTree*)f->Get("mod0");

  ofstream of("landau_193.csv");

  for(int i=0; i<20; i++){
    tree->Draw(Form("Qval>>q%i", i), Form("iDAQ==%i", i), "goff");
    c->cd(i+1);
    charge[i]->Draw();
    of << i << " " << MPV(charge[i]) << endl;

  }

  of.close();

}
