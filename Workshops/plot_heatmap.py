import sys
import numpy as np
import matplotlib.pyplot as plt

# usage: python plot_heatmap.py [csv_file]
fname = sys.argv[1] if len(sys.argv) > 1 else \
    "Entanglement_entropy_renyi_1_n_24_J2_0.0_Deltaz_1.0_OBC_stot_.csv"

data = np.loadtxt(fname, delimiter=",")
t = data[:, 0]
S = data[:, 1:]                      # rows: time, columns: bonds 1..N-1
bonds = np.arange(1, S.shape[1] + 1)

plt.figure(figsize=(6, 5))
plt.imshow(S, origin="lower", aspect="auto", cmap="viridis",
           extent=[bonds[0] - 0.5, bonds[-1] + 0.5, t[0], t[-1]])
plt.colorbar(label=r"$S_{vN}$")
plt.xlabel("bond")
plt.ylabel("t")
plt.title("Entanglement entropy after local quench")
plt.tight_layout()
plt.savefig(fname.replace(".csv", "_heatmap.png"), dpi=200)
plt.show()
