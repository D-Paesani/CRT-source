import ROOT
from array import array

def get_matched_trees(crt_tree, dirac_tree):

  evnum = array('i', [0])

  crt_newtree = crt_tree.CloneTree(0)
  dirac_newtree = dirac_tree.CloneTree(0)

  crt_newtree.Branch("evnum", evnum, "evnum/I")
  crt_newtree.Branch("evnum", evnum, "evnum/I")


  events = newtree.GetEntries()

  crt_lost = 0
  dirac_lost = 0

  for i in range(events):

    crt_tree.GetEntry(i)
    dirac_tree.GetEntry(i)


# matcho (se non è tutto ok continue)

#riempio la branch nuova per newtree crt e newtree dirac con un indice da contatore

  newtree.Fill() # per crt e dirac

newtree.Write()




# gli alberi per frammento vanno in una lista, che viene mandata a MergeTree (per crt e dirac)




