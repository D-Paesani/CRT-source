import ROOT
import matplotlib.pyplot as plt
import numpy as np
from array import array
import sys

n_boards = 2

dirac_files = [ ROOT.TFile(f"../../data/dirac/step2/run_23all_step2_dirac{i}.root") for i in range(n_boards) ]

dirac_trees = [file.Get("mod0") for file in dirac_files]

outf = ROOT.TFile("out.root", "recreate")

new_trees = []
to_keep = [array("i", [0])]*len(dirac_trees)

for n, tree in enumerate(dirac_trees): # 2 volte
  tree.BuildIndex("iTrigGlobal")
  new_tree = tree.CloneTree(0)
  new_tree.Branch("keep", to_keep[n], "keep/I")
  new_trees.append(new_tree)

n_trigs = int(min([tree.GetMaximum("iTrigGlobal") for tree in dirac_trees]))

for i in range(n_trigs+1):
  entries = []
  all_OK = 1
  for tree in dirac_trees:
    entry = tree.GetEntryNumberWithIndex(i)
    if entry <= 0:
      all_OK = 0
    else:
      entries.append(entry)
  if all_OK == 1:
    for n, tree in enumerate(new_trees):
      dirac_trees[n].GetEntry(entries[n])
      to_keep[n][0] = 1
      tree.Fill()

clean_trees = []

outf.cd()
for n, tree in enumerate(new_trees):
  clean_tree = tree.CopyTree("keep==1")
  clean_tree.SetName(f"b{n}")
  clean_tree.SetTitle(f"b{n}")
  clean_trees.append(clean_tree)
  clean_tree.Write()


master = clean_trees[0].CloneTree()
master.SetName("mod0_sync")
master.SetTitle("mod0_sync")

for tree in clean_trees:
  master.AddFriend(tree)

master.Write()
outf.Write()
outf.Close()


outf = ROOT.TFile("out.root", "update")
outf.cd()
ROOT.gDirectory.Delete("mod0;*");
outf.Write()
outf.Close()

# mod0_sync.Draw("b0.Qval:b1.Qval")
