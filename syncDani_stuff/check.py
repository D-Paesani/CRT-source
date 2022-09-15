import ROOT
import matplotlib.pyplot as plt
import numpy as np

dirac_f = ROOT.TFile("../data/roottople/run165sd.root")
dirac_tree = dirac_f.Get("mod0")

crt_f = ROOT.TFile("../data/crt/crt_run165sd.root")
crt_tree = crt_f.Get("CRT")

crt_bar = []
dirac_bar = []

crt_evnum = []
dirac_evnum = []

for entry in crt_tree:
    crt_bar.append(entry.iScint[1])
    crt_evnum.append(entry.iTrig-706)

for entry in dirac_tree:
   if entry.crtBarBot < 10:
     dirac_bar.append(entry.crtBarBot)
     dirac_evnum.append(entry.evnum)

# plt.plot(crt_evnum, crt_bar)
# plt.plot(dirac_evnum, dirac_bar)
# plt.show()

