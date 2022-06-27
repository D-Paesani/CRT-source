import ROOT
import matplotlib.pyplot as plt
import numpy as np

dirac_f = ROOT.TFile("run_165_all.root")
dirac_tree = dirac_f.Get("mod0")

crt_f = ROOT.TFile("../../step2/run165_ana.root")
crt_tree = crt_f.Get("CRT")

crt_bar = []
dirac_bar = []

crt_evnum = []
dirac_evnum = []

for entry in crt_tree:
  for hit in range(entry.nCry):
     if entry.iMod[hit] == 1 and entry.iSide[hit]==0:
       crt_bar.append(entry.iScint[hit])
       crt_evnum.append(entry.evnum-706)

for entry in dirac_tree:
   if entry.crtBarBot < 10:
     dirac_bar.append(entry.crtBarBot)
     dirac_evnum.append(entry.evnum)

plt.plot(crt_evnum, crt_bar)
plt.plot(dirac_evnum, dirac_bar)
plt.show()
