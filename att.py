import os
import pandas as pd

df = pd.read_csv("q.org", header=None, sep=" ")
f0 = open("att0.csv", "w")
f1 = open("att1.csv", "w")

for i in range(len(df)):
  nrun = df.iloc[i, 0]
  pos = df.iloc[i, 1]
  #os.system(f"./CRT_step3.sh run{nrun} run217");
  q0 = pd.read_csv(f"../../data/calibration/luts_s3/run{nrun}_chargEq.csv", header=None, sep=",").iloc[0, 0]
  q0_err = pd.read_csv(f"../../data/calibration/luts_s3/run{nrun}_chargEqErr.csv", header=None, sep=",").iloc[0, 0]
  q1 = pd.read_csv(f"../../data/calibration/luts_s3/run{nrun}_chargEq.csv", header=None, sep=",").iloc[1, 0]
  q1_err = pd.read_csv(f"../../data/calibration/luts_s3/run{nrun}_chargEqErr.csv", header=None, sep=",").iloc[1, 0]
  f0.write(f"{q0} {q0_err} {pos}\n")
  f1.write(f"{q1} {q1_err} {pos}\n")

f0.close()
f1.close()
