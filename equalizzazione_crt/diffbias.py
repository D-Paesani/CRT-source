import pandas as pd
import numpy as np
import sys

run_old = sys.argv[1]
run_new = sys.argv[2]

p0, p1 = 1, 0.176

mpv = pd.read_csv(f"mpv{run_old}.csv", header=None, sep=" ")
print(mpv)

target_charge = 4200

a = target_charge/mpv * (p0)

biasdiff = (a - p0)/p1

biasdiff = biasdiff.round(2).iloc[:, 1]

print(f"BIAS DIFFERENCES WRT VOP:\n{biasdiff}")
biasdiff.to_csv(f"biasdiff{run_new}.csv", header=None, index=False)
