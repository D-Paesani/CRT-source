{
    double tdiff[15], tdiff_err[15], pos_err[15], pos[15];
    int i=0;
    ifstream a("time_diff_new.csv");
    while(a >> tdiff[i] >> tdiff_err[i] >> pos[i]){ i++; pos_err[i]=0.5;};
    auto *g = new TGraphErrors(15, pos, tdiff, pos_err, tdiff_err);
/*    TF1 *f = new TF1("f", "2*(x - [0]/2)/[1]", 0, 160, 2);
    f->SetParameters(160, 14);
    gStyle->SetOptFit(1);
    g->Fit(f, "R");*/
    g->Draw();
}
