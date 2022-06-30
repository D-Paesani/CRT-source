
const int Nboard = 1;
const int Nchan = 20;

int BoaTmp, ChaTmp, DaqTmp, HitTmp, RowTmp, ColTmp, SipmTmp, CryTmp;
int rowID[Nboard][Nchan], colID[Nboard][Nchan], sipmID[Nboard][Nchan];

int get_iCry(int jRow, int jCol)
{
  int nCol[7] = {6, 7, 8, 9, 8, 7, 6};

  int iVal = jCol;
  for (int iRow = 0; iRow < jRow; iRow++)
  {
    iVal = iVal + nCol[iRow];
  }

  return iVal;
}

void LandauMPVforAllChannels2()
{
  auto *f = new TFile("../../data/mon_run178.root");

  TH1F *charge[20];
  for (int i = 0; i < 20; i++)
    charge[i] = new TH1F(Form("q%i", i), Form("q%i", i), 100, 0, 15000);

  auto *c = new TCanvas("c", "c");
  c->Divide(5, 4);

  ofstream of("landau2.csv");

  FILE *mapFile;
  mapFile = fopen("mod0_dirac_map.dat", "r");
  for (int Iloop = 0; Iloop < Nboard * Nchan; Iloop++)
  {
    fscanf(mapFile, "%d %d %d %d %d", &BoaTmp, &ChaTmp, &RowTmp, &ColTmp, &SipmTmp);
    rowID[BoaTmp][ChaTmp] = RowTmp;
    colID[BoaTmp][ChaTmp] = ColTmp;
    sipmID[BoaTmp][ChaTmp] = SipmTmp;
  }
  fclose(mapFile);

  for (int i = 0; i < 20; i++)
  {

    RowTmp = rowID[0][i];
    ColTmp = colID[0][i];
    SipmTmp = sipmID[BoaTmp][ChaTmp];

    CryTmp = get_iCry(RowTmp, ColTmp);

    TString h_name = Form("eneEqPreliminary/eneEqPreliminary_%i%c_r%i_c%i_s%i", CryTmp, SipmTmp ? 'R': 'L', RowTmp, ColTmp, SipmTmp);

    charge[i] = (TH1F *)f->Get(h_name);

    c->cd(i + 1);

    charge[i]->Draw();
    double mpv = charge[i]->GetFunction("land")->GetParameter(1);

    of << i << " " << mpv << endl;
  }

  of.close();
}
