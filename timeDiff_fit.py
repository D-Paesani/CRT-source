import ROOT
import pandas as pd

file = ROOT.TFile("timeDiff_stronzio_10.22.root")
tree = file.Get("timeDiff")

vp = {0: {}, 1: {}}
vp_err = {0: {}, 1: {}}

for mod in [0, 1]:
  for scint in [3, 5, 7]:
    n = tree.Draw(f"mean:dist:mean_err:mean_err*0.2/mean_err", f"mod=={mod} && scint=={scint} && dist<125 && dist > 15", "goff")
    mean = tree.GetV1();
    dist = tree.GetV2();
    mean_err = tree.GetV3()
    dist_err = tree.GetV4()
    g = ROOT.TGraphErrors(n, dist, mean, dist_err, mean_err)
    c = ROOT.TCanvas(f"timeDiff_{mod}_{scint}")
    g.Draw()
    input("")
    g.Fit("pol1")
    f = g.GetFunction("pol1")
    print(f.GetProb())
    vp[mod][scint] = 2/f.GetParameter(1)
    vp_err[mod][scint] = 2* f.GetParError(1)/f.GetParameter(1)/f.GetParameter(1)

pd.DataFrame(vp).to_csv("vp.csv")
pd.DataFrame(vp_err).to_csv("vp_err.csv")
