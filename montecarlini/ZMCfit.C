double flat(const double *x, const double *par){
  double ampl = par[0];
  double len = par[1];
  if (x[0] < -len || x[0] > len) return 0;
  else return ampl;
}


void ZMCfit(){

   /*
   // Construction of histogram to fit.
   TH1F *h_Data = new TH1F("h_Data", "Unif convoluted by Gaussian", 100, -100, 100);
   for (int i = 0; i < 1e6; i++) {
      // Gives a alpha of -0.3 in the Unif.
      double x = gRandom->Uniform(-80, 80);
      x += gRandom->Gaus(0., 3.);
      // Probability density function of the addition of two variables is the
      // convolution of two density functions.
      h_Data->Fill(x);
   }
   */

   TFile *inFile = new TFile("../../data/step3/run205_s3_Qgreaterthan200pC.root");

   TH1F *h_Data = (TH1F*) inFile->Get("zetaMip/zetaMip_0");
   TF1 *flat_f = new TF1("flat_f", flat, -100, 100, 2);
   TF1 *gauss_f = new TF1("gauss_f", "gaus", -100, 100);
   TF1Convolution *f_conv = new TF1Convolution(flat_f, gauss_f, -80, 80, true);
   f_conv->SetNofPointsFFT(1000);
   TF1 *f = new TF1("f", *f_conv, -80, 80, f_conv->GetNpar());
   f->SetParameters(6., 74, 6, 0, 3);

   // Fit.
   new TCanvas("c", "c", 800, 1000);
   h_Data->Fit("f", "R");
   h_Data->Draw();
}
