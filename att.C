#include <TMath.h>

void att(){
    double q0[7], q0_err[7], q1[7], q1_err[7], pos[7], ex[7];
    int i=0;
    ifstream a0("../data/att0_new.csv");
    while(a0 >> q0[i] >> q0_err[i] >> pos[i]) {ex[i] = 0.5; i++;}
    i=0;
    ifstream a1("../data/att1_new.csv");
    while(a1 >> q1[i] >> q1_err[i] >> pos[i]) i++;
    cout << "Sbilanciamento: " << q0[0]/q1[0] << endl;
    for(int i=0; i<7; i++){
      q0_err[i] = q0[i] * 0.025;
      q1_err[i] = q1[i] * 0.025;
    }

    auto *g0 = new TGraphErrors(7, pos, q0, ex, q0_err);
    auto *g1 = new TGraphErrors(7, pos, q1, ex, q1_err);

    gStyle->SetOptFit(1);

    g1->Draw("AP");
    g0->Draw("P");

    TF1 l = TF1("l", "[0] * (TMath::Exp(-x/380) + [1]*TMath::Exp(-x/[2]) )", 0, 160);
    l.SetParameters(100, 3, 30);
    l.SetParLimits(1, 0.2, 5);
    l.SetParName(0, "Normalization");
    l.SetParName(1, "L_{I}/L_{D}");
    l.SetParName(2, "TAL");

    g0->Fit(&l, "R", "", 0, 160);
    g0->Fit(&l, "R", "", 0, 160);
    g0->Fit(&l, "R", "", 0, 160);

    l = TF1("l", "[0] * (TMath::Exp(-(160-x)/380) + [1]*TMath::Exp(-(160-x)/[2]) )", 0, 160);
    l.SetParameters(100, 3, 30);
    l.SetParLimits(1, 0.2, 5);
    l.SetParName(0, "Normalization");
    l.SetParName(1, "L_{I}/L_{D}");
    l.SetParName(2, "TAL");

    g1->Fit(&l, "R", "", 0, 160);
    g1->Fit(&l, "R", "", 0, 160);
    g1->Fit(&l, "R", "", 0, 160);

}
