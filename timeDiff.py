import os
import pandas as pd

df = pd.read_csv("../qnew.org", header=None, sep=" ")
f = open("../time_diff_new2.csv", "w")
for i in range(len(df)):
  nrun = df.iloc[i, 0]
  pos = df.iloc[i, 1]
  #os.system(f"./CRT_step3.sh run{nrun} run236");
  tdiff = pd.read_csv(f"../../data/calibration/luts_s3/run{nrun}_timeDiff.csv", header=None, sep=",").iloc[0, 7]
  tdiff_err = pd.read_csv(f"../../data/calibration/luts_s3/run{nrun}_timeDiffErr.csv", header=None, sep=",").iloc[0, 7]
  f.write(f"{tdiff} {tdiff_err} {pos}\n")
  print(f"{tdiff} {tdiff_err} {pos}")
f.close()
