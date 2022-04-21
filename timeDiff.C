{
    double tdiff[7], tdiff_err[7], pos_err[7], pos[7];
    int i=0;
    ifstream a("time_diff_new2.csv");
    while(a >> tdiff[i] >> tdiff_err[i] >> pos[i]){ pos_err[i]=1; i++;};
    auto *g = new TGraphErrors(6, pos, tdiff, pos_err, tdiff_err);
/*    TF1 *f = new TF1("f", "2*(x - [0]/2)/[1]", 0, 160, 2);
    f->SetParameters(160, 14);
    gStyle->SetOptFit(1);
    g->Fit(f, "R");*/
    g->Draw();
}
