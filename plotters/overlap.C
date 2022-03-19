void overlap(){
  TFile *f = new TFile("../../data/step3/run204_s3.root");
  int colors[8] = {kPink+1, kBlack, kRed, kBlue, kGreen, kViolet, kCyan, kOrange};
  TLegend *lg = new TLegend();
  gStyle->SetOptStat(0);

  for(int i=0; i<8; i++){
    auto *h = (TH1F*)f->Get(Form("chargeMip/chargeMip_1_%i", i));
    h->GetFunction("l")->SetBit(TF1::kNotDraw);
    h->Scale(10000/h->Integral());
    h->Draw("HISTO SAME");
    h->SetLineColor(colors[i]);
    lg->AddEntry(h, Form("Scint. %i", i));
  }
  lg->Draw();
}
