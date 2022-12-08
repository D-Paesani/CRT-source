
import pandas as pd
import matplotlib.pyplot as plt
import sys

run = sys.argv[1]

ref = float(sys.argv[2])/100  #percentuale - es 10 -> plotta rette a q_ref +/- 10%

q_ref = 350

df = pd.read_csv(f"mpv{run}.csv", header=None, sep=" ")

df.columns = ["q"]

plt.scatter(list(range(20)), df.q)

plt.ylim(0, 500)

plt.xticks(list(range(20)))

plt.plot(list(range(20)), [q_ref]*20)

plt.plot(list(range(20)), [q_ref*(1+ref)]*20)

plt.plot(list(range(20)), [q_ref*(1-ref)]*20)

plt.text(7, q_ref*1.20, f"{q_ref} pC +/- {100*ref:.0f}%")


plt.xlabel("Dirac Channel");

plt.ylabel("Cosmic-rays Landau MPV [pC]")

plt.show()
