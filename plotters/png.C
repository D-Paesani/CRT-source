{
  TFile *f = new TFile("run205_s3_Qgreaterthan200pC.root");
  gStyle->SetOptFit(1);
  int i = 0;
  TH1F *tmp;

    TCanvas *c = new TCanvas("c", "c");
    c->cd();
  for(int i=0; i<8; i++){
    cout << i << endl;
    auto tmp = (TH1F*)f->Get(Form("zetaMip/zetaMip_%i", i));
    tmp->Draw();
    gPad->Update();
    tmp->GetXaxis()->SetRangeUser(-110, 110);
    auto ps = (TPaveStats *)tmp->GetListOfFunctions()->FindObject("stats");
    ps->SetX1NDC(0.3);
    ps->SetX2NDC(0.7);
    ps->SetY1NDC(0.55);
    ps->SetY2NDC(0.15);
    gPad->Modified();
    c->SaveAs(Form("z_run205_qgt200pC_scint%i.png", i));
    //c->Delete();
  }
}
