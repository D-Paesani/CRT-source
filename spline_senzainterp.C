#include <fstream>
#include <TApplication.h>
#include "Analysis.h"
#include <TSpline.h>
#include <TButton.h>
#include <TSystem.h>

using namespace std;

double landau_mu[sideNum][scintNum], landau_mu_err[sideNum][scintNum];

// ------------------------------------------------------------ TO DO
// modificare is_mip (prima o poi) per raffinare il taglio
// sistemare il codice per essere pronti a equalizzare
// plot sigma/MPV
// asymm vds. Q

// Histograms to create
void Analysis::CreateHistDict(vector <Double_t> pars){
  hist_dict = {
    CreatePair("interp_wave", 2, 8, 2, "Interpolated waveform histogram", "%i%i",   "Time", "ns", 1600, 0, 800, "Voltage", "mV", 500, 0, 1.2),
  };
}

int go_ahead_flag = 0;
int plot = 0;

TCanvas *c;
TPad *plots_pad;

void go_ahead(){
  go_ahead_flag = 1;
}

void prof_draw_single(TH1 *obj, int ndim, TFile *outfile)
{
  obj->SetFillStyle(0);

  obj->Draw("zcol");

  TProfile *prof = ((TH2*)obj)->ProfileX();
  prof->SetName(Form("%s_profile", obj->GetName()));

  TSpline5 *sp5 = new TSpline5(prof);
  sp5->SetName(Form("spline_%s", obj->GetName()));
  sp5->SetLineColor(kRed);
  sp5->Draw("l same");

  outfile->cd();

  prof->Write();
  sp5->Write();

  obj->Write();
}

void Analysis::Loop(){
  if (fChain == 0) return;

  Long64_t nentries = fChain->GetEntriesFast();
  CreateHistDict({});

  Long64_t nbytes = 0, nb = 0;
  int sideTmp, scintTmp, canv=0, ncanvas=3;
  double Qtmp, Ttmp, Chi2tmp, Vtmp;
  double interp_time_tmp;
  double vp = 13; //cm/ns
  
  double xerr[200]={0}, yerr[200]={0.75};

  double interp_wave[2400], interp_time[2400];

  // LOOP OVER ENTRIES

  if(plot){
    c = new TCanvas("c", "c");
    plots_pad = new TPad("pad1", "pad1",0.0, 0.1, 1.0, 1.0);
    TPad *pad1 = new TPad("pad2", "pad2",0.0, 0, 1.0, 0.1);
    plots_pad->Draw();
    pad1->Draw();
    pad1->cd();
    TButton *but1 = new TButton("Next", "go_ahead()", 0.45, 0.15, 0.55, 0.95);
    but1->Draw();
    plots_pad->cd();
  }

  for (Long64_t jentry=0; jentry<20000;jentry++) { //la parte interna al loop andrebbe messa una una funzione così come le parti prima e dopo, così la parte delicata sta in Loop nel .h
    Long64_t ientry = LoadTree(jentry);

    if (ientry < 0) break;
    if(jentry%1000 == 0)
    cout << Form("Processing event n.%lld of %lld: %i%%", jentry, nentries, (int)((double)jentry/nentries * 100)) << endl;

    nb = fChain->GetEntry(jentry);   
    nbytes += nb;

    // LOOP OVER HIT
    for(int hit=0; hit<nCry; hit++){
      sideTmp=iSide[hit];
      scintTmp = iScint[hit];
      Qtmp = Qval[hit];
      double DTtmp = templTime[hit] - Tval[hit];

      if(Qtmp > 300 && Qtmp < 800 && DTtmp > -130 && DTtmp < -120 && templChi2[hit] < 10 && templTime[hit] < 200 && templTime[hit] > 100){

        TGraphErrors meas(500, time, wave[hit], xerr, yerr); //300 ns from 200 to 500
        TSpline5 s("grs", &meas);
        s.SetLineColor(kRed);
        s.Draw("same");

        for(int i=0; i<300; i++){
            interp_time_tmp = 0.5*i;
            interp_wave[i] = s.Eval(interp_time_tmp)/Vmax[hit]; //va usato picco spline
            interp_time[i] = interp_time_tmp;
            if(!plot && Qtmp > 300 && Qtmp < 800) GetHist(
                "interp_wave", scintTmp, sideTmp
              )->Fill(4*i - templTime[hit] + 100 - 4, wave[hit][i]/Vmax[hit]);
        }

        if (plot) {
          TGraph g(2400, interp_time, interp_wave);

          meas.SetMarkerStyle(22);
          meas.SetMarkerSize(1.5);
          meas.SetTitle(Form("Spline waveform fit - side %i - scint %i", sideTmp, scintTmp));
          meas.GetXaxis()->SetTitle("Time (ns)");
          meas.GetYaxis()->SetTitle("Voltage (mV)");
          meas.Draw("ape");

          g.Draw("p same");
          g.SetMarkerColor(kGreen);
          g.SetMarkerStyle(21);
          g.SetMarkerSize(0.8);


          while (go_ahead_flag == 0)
          {
            gPad->Modified();
            gPad->Update();
            gSystem->ProcessEvents();
            usleep(10000);
          }
          go_ahead_flag = 0;

        }
        
      }

    }
  
  }


  hist_dict["interp_wave"]->draw_single = &prof_draw_single;

  for (auto& hist_matrix: hist_dict) {
    hist_matrix.second->draw_all();
  }

  CloseOutFile();

}

int main(int argc, char*argv[]){ 

  if (argc != 3)
  {
    printf("Usage: %s [infile_name] [outfile_name]\n", argv[0]);
    exit(-1);
  }

  Analysis::Run(argv[1], argv[2], 1);
}

void spline_senzainterp(){
  Analysis::Run("run182_ana.root", "waves_splines_run182_fromRubenFit_new2.root", -1);
}
