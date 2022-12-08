import ROOT
import matplotlib.pyplot as plt
import numpy as np

dirac_f = ROOT.TFile("../../data/dirac/step2/run_165all_step2_sync.root")
dirac_tree = dirac_f.Get("mod0")

crt_f = ROOT.TFile("../../data/step2/run284_s2.root")
crt_tree = crt_f.Get("CRT")

crt_bar = []
dirac_bar = []

crt_evnum = []
dirac_evnum = []
crt_entry = []


for n, entry in enumerate(crt_tree):
  if n%1000 == 0: print(n)
  if entry.nCry > 4: continue
  for hit in range(entry.nCry):
     if entry.iMod[hit] == 1 and entry.iSide[hit]==0:
       crt_bar.append(entry.iScint[hit])
       crt_evnum.append(entry.evnum)
       crt_entry.append(n)

print(f"Fraction of good crt events: {len(crt_evnum)/crt_tree.GetEntries()}")

for n, entry in enumerate(dirac_tree):
   if n%1000 == 0: print(n)
   if entry.crtBarBot < 10:
     dirac_bar.append(entry.crtBarBot)
     dirac_evnum.append(entry.iTrigGlobal-1)


eff = []

crt_evnum = np.asarray(crt_evnum)
good = 0
for i, n in enumerate(dirac_evnum):
  inds = np.where(crt_evnum == n)[0]
  if len(inds) != 1: continue
  else: ind = inds[0]
  if crt_bar[ind] == dirac_bar[i]:
    good += 1

  if i%1000 == 0:
    if good < 1000: eff.append(good/1000)
    good = 0
    print(i)

'''
plt.plot(crt_evnum,   crt_bar, color="blue")
plt.plot(dirac_evnum, dirac_bar, color="red")

plt.scatter(crt_evnum,   crt_bar, color="blue")
plt.scatter(dirac_evnum, dirac_bar, color="red")

plt.xlabel("eventi CRT / MOD0")
plt.ylabel("barra CRT secondo CRT / MOD0")

plt.show()
'''

plt.plot(eff)
plt.xlabel("Events / 1000")
plt.ylabel("Efficiency")
plt.show()
