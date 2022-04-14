{
    double q0[15], q0_err[15], q1[15], q1_err[15], pos[15], ex[15];
    int i=0;
    ifstream a0("att0.csv");
    while(a0 >> q0[i] >> q0_err[i] >> pos[i]) {i++; ex[i] = 0.5;}
    auto *g0 = new TGraphErrors(15, pos, q0, ex, q0_err);
    i=0;
    ifstream a1("att1.csv");
    while(a1 >> q1[i] >> q1_err[i] >> pos[i]) i++;
    for(int i=0; i<15; i++){
      q1[i] *= q0[8]/q1[11];
    }
    auto *g1 = new TGraphErrors(15, pos, q1, ex, q1_err);
    g0->Draw("AP");
    g1->Draw("P");

    TF1 l = TF1("l", "expo(0)+expo(2)", 0, 160);
    l.SetParLimits(3, 0.01, 0.09);
    l.SetParLimits(1, 0.0001, 0.03);

    g1->Fit(&l, "R", "", 15, 160);
    TLatex *lt[4];
    lt[0] = new TLatex(80, 200, Form("BAL (side %i): %.1f +/- %.1f cm ", 1, 1/l.GetParameter(1), 1/(l.GetParameter(1)*l.GetParameter(1)) *l.GetParError(1) ));
    lt[1] = new TLatex(80, 150, Form("TAL (side %i): %.1f +/- %.1f cm ", 1, 1/l.GetParameter(3), 1/(l.GetParameter(3)*l.GetParameter(3)) *l.GetParError(3) ));

    l.SetParLimits(3, -0.09, -0.01);
    l.SetParLimits(1, -0.03, -0.0001);
    g0->Fit(&l, "R", "", 0, 145);

    lt[2] = new TLatex(80, 100, Form("BAL (side %i): %.1f +/- %.1f cm ", 0, -1/l.GetParameter(1), 1/(l.GetParameter(1)*l.GetParameter(1)) *l.GetParError(1) ));
    lt[3] = new TLatex(80, 50 , Form("TAL (side %i): %.1f +/- %.1f cm", 0, -1/l.GetParameter(3), 1/(l.GetParameter(3)*l.GetParameter(3)) *l.GetParError(1) ));

    for(int i=0; i<4; i++){
      lt[i]->Draw();
    }
}
