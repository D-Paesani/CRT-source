import ROOT
import pandas as pd
from array import array

def timeDiff(nrun, chan0, chan1):
  file = ROOT.TFile(f"../data/stronzio_10.22/run{nrun}_ana.root")
  tree = file.Get("CRT")
  h = ROOT.TH1F(f"tdiff{nrun}", "tdiff", 400, -20, 20);
  tree.Draw(f"templTime[0] - templTime[1]>>tdiff{nrun}", f"iDAQ[0]=={chan0} && iDAQ[1]=={chan1}", "goff")
  if h.GetEntries() == 0: tree.Draw(f"templTime[0] - templTime[1]>>tdiff{nrun}", f"iDAQ[0]=={chan0} && iDAQ[1]=={chan1}", "goff")
  gaus = ROOT.TF1(f"gaus{nrun}", "gaus", h.GetMean() - 0.5*h.GetRMS(), h.GetMean() + 0.5*h.GetRMS())
  h.Fit(f"gaus{nrun}", "RQ0")
  pars = gaus.GetParameters()
  mean, sigma = pars[1], pars[2]
  par_errors = gaus.GetParErrors()
  mean_err, sigma_err = par_errors[1], par_errors[2]
  return (mean, sigma, mean_err, sigma_err)


df = pd.read_csv("../data/stronzio_10.22/run_stronzio_ottobre2022.csv", sep=',')

tree = ROOT.TTree("timeDiff", "timeDiff")

mean = array('f', [0.] )
sigma = array('f', [0.] )
mean_err = array('f', [0.] )
sigma_err = array('f', [0.] )
scint = array('i', [0])
mod = array('i', [0])
dist = array('i', [0])

tree.Branch('mean', mean, 'mean/F')
tree.Branch('sigma', sigma, 'sigma/F')
tree.Branch('mean_err', mean_err, 'mean_err/F')
tree.Branch('sigma', sigma_err, 'sigma_err/F')
tree.Branch("scint", scint, "scint/I")
tree.Branch("mod", mod, "mod/I")
tree.Branch("dist", dist, "dist/I")

for i in range(len(df)):
  row = df.loc[i]
  mean[0], sigma[0], mean_err[0], sigma_err[0] = timeDiff(row.run, row.caen, row.caen+4)
  scint[0] = row.scint
  mod[0] = row.modulo
  dist[0] = row.distfromside0
  tree.Fill()

tree.SaveAs("str.root")
