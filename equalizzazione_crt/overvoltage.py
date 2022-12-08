import pandas as pd
import numpy as np
import sys

run_old = sys.argv[1]
run_new = sys.argv[2]

p0, p1 = 1, 0.176

mpv = pd.read_csv(f"mpv{run_old}.csv", header=None, sep=" ")
print(mpv)

old_volt = pd.read_csv(f"volt{run_old}.csv", header=None)
print(old_volt)

vop = pd.read_csv("vop.csv", header=None)
print(vop)

b = (440/mpv)
c = (p0 + p1*(old_volt - vop))
a = b*c

new_volt = vop + (a - p0)/p1

new_volt = new_volt.round(2).iloc[:, 0]

print(f"NEW VOLTAGES:\n{new_volt}")
new_volt.to_csv(f"volt{run_new}.csv", header=None, index=False)
