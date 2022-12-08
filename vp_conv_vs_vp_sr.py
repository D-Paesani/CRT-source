import pandas as pd
import matplotlib.pyplot as plt
import ROOT
from array import array

vp_sr = pd.read_csv("vp_stronzio_10.22.csv")
vp_err_sr = pd.read_csv("vp_err_stronzio_10.22.csv")

conv_B_vp = pd.read_csv("../data/calibration/luts_s3/run284_B_barLen.csv", sep=",", header=None).T
conv_T_vp = pd.read_csv("../data/calibration/luts_s3/run284_T_barLen.csv", sep=",", header=None).T

conv_B_vp_err = pd.read_csv("../data/calibration/luts_s3/run275_B_barLenErr.csv", sep=",", header=None).T
conv_T_vp_err = pd.read_csv("../data/calibration/luts_s3/run275_T_barLenErr.csv", sep=",", header=None).T

conv_B_vp.columns = ["len_nominalvp"]
conv_T_vp.columns = ["len_nominalvp"]

conv_B_vp_err.columns = ["len_nominalvp"]
conv_T_vp_err.columns = ["len_nominalvp"]

len_old_B = conv_B_vp["len_nominalvp"].loc[[2, 4, 6]]
len_old_T = conv_T_vp["len_nominalvp"].loc[[2, 4, 6]]

vp_cosmics_B = 14.3 / len_old_B * 160
vp_err_cosmics_B = conv_B_vp_err["len_nominalvp"].loc[[2, 4, 6]]/len_old_B * vp_cosmics_B

vp_cosmics_T = 14.3 / len_old_T * 160
vp_err_cosmics_T = conv_T_vp_err["len_nominalvp"].loc[[2, 4, 6]]/len_old_T * vp_cosmics_T

c_B = ROOT.TCanvas("c_B", "Bottom Module Vp comparison")
c_T = ROOT.TCanvas("c_T", "Top module Vp comparison")

c_B.cd()
g = ROOT.TGraphErrors(3, array("f", vp_cosmics_B), array("f", vp_sr["0"]), array("f", vp_err_cosmics_B), array("f", vp_err_sr["0"]))
g.Draw("AP")
g.GetXaxis().SetTitle("V_{p} from cosmics [cm/ns]")
g.GetYaxis().SetTitle("V_{p} from Sr/Y source [cm/ns]")
g.SetTitle("Bottom module")
g.Fit("pol1")
c_B.SaveAs("vp_B_comparison.root")

c_T.cd()
g = ROOT.TGraphErrors(3, array("f", vp_cosmics_T), array("f", vp_sr["1"]), array("f", vp_err_cosmics_T), array("f", vp_err_sr["1"]))
g.Draw("AP")
g.GetXaxis().SetTitle("V_{p} from cosmics [cm/ns]")
g.GetYaxis().SetTitle("V_{p} from Sr/Y source [cm/ns]")
g.SetTitle("Tottom module")
g.Fit("pol1")
c_T.SaveAs("vp_T_comparison.root")

input()

'''
plt.plot(14.3 / conv_B.len_nominalvp * 160 , sr["0"])
plt.xlabel("Vp from cosmics [cm/ns]")
plt.ylabel("Vp from Sr source [cm/ns]")
plt.show()


len_old = dt_real * 14.3

160 = dt_real * vp

-> ->

160/len_old = vp/14.3
'''
