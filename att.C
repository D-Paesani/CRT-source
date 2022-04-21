#include <TMath.h>

void att(){
    double q0[7], q0_err[7], q1[7], q1_err[7], pos[7], ex[7];
    int i=0;
    ifstream a0("att0.csv");
    while(a0 >> q0[i] >> q0_err[i] >> pos[i]) {ex[i] = 0.5; i++;}
    auto *g0 = new TGraphErrors(7, pos, q0, ex, q0_err);
    i=0;
    ifstream a1("att1.csv");
    while(a1 >> q1[i] >> q1_err[i] >> pos[i]) i++;
    cout << "Sbilanciamento: " << q0[0]/q1[0] << endl;
/*    for(int i=0; i<7; i++){
      q1[i] *= q0[0]/q1[0];
    }*/
    auto *g1 = new TGraphErrors(7, pos, q1, ex, q1_err);
    g1->Draw("AP");
    g0->Draw("P");

    TF1 l = TF1("l", "[0] * (TMath::Exp(-x/[1]) + [2]*TMath::Exp(-x/[3]) )", 0, 160);
    l.SetParameters(1, 300, 0.5, 20);
    //l.FixParameter(1, 0.0026);

    g0->Fit(&l, "R", "", 0, 160);
    TLatex *lt[4];
    lt[0] = new TLatex(80, 200, Form("BAL (side %i): %.1f +/- %.1f cm ", 1, 1/l.GetParameter(1), 1/(l.GetParameter(1)*l.GetParameter(1)) *l.GetParError(1) ));
    lt[1] = new TLatex(80, 60, Form("TAL (side %i): %.1f +/- %.1f cm ", 1, 1/l.GetParameter(3), 1/(l.GetParameter(3)*l.GetParameter(3)) *l.GetParError(3) ));

    l = TF1("l", "[0] * (TMath::Exp(-(160-x)/[1]) + [2]*TMath::Exp(-(160-x)/[3]) )", 0, 160);
    l.SetParameters(10, 300, 0.5, 20);

    //l.FixParameter(1, -0.0026);
    g1->Fit(&l, "R", "", 0, 160);

    lt[2] = new TLatex(80, 100, Form("BAL (side %i): %.1f +/- %.1f cm ", 0, -1/l.GetParameter(1), 1/(l.GetParameter(1)*l.GetParameter(1)) *l.GetParError(1) ));
    lt[3] = new TLatex(80, 50 , Form("TAL (side %i): %.1f +/- %.1f cm", 0, -1/l.GetParameter(3), 1/(l.GetParameter(3)*l.GetParameter(3)) *l.GetParError(1) ));

/*    for(int i=0; i<4; i++){
      lt[i]->Draw();
    }*/
}
